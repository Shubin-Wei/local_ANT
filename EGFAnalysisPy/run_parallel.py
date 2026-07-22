"""
Parallel driver for EGFAnalysisPy.

Same per-file processing as run.py, but distributed across worker processes.
Design notes (TensorFlow 1.x specific):
  - TF sessions are not picklable, so each worker builds its OWN Pick() once,
    in the pool initializer, and reuses it for all files it handles.
  - The main process never imports/initializes TF (it only lists files and
    manages the pool); workers are started with the 'spawn' method to avoid
    the TF-fork deadlock.
  - Native thread counts are pinned to 1 per process so that parallelism comes
    from processes, not from oversubscribed BLAS/FFT thread pools.

Usage:
    NWORKERS=16 python run_parallel.py        # or just: python run_parallel.py
    python run_parallel.py --config config_line2.py -w 8
"""

import os

# Pin native threads BEFORE importing numpy/scipy/tf. With 'spawn' each child
# re-imports this module, so these are set in every process.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
# Silence TensorFlow C++ INFO/WARNING banner (AVX2/XLA/etc.) which is written
# straight to fd 2 and would otherwise clutter the terminal at worker startup.
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

import sys
import argparse
import importlib.util
import random
import traceback
import multiprocessing as mp

import numpy as np
from tqdm import tqdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DISPER_PICKER_DIR = os.path.join(SCRIPT_DIR, 'DisperPicker')
for _p in (SCRIPT_DIR, DISPER_PICKER_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Per-worker globals (populated by _init_worker in each process).
_pick = None
_config = None
_T = None


def _resolve_config_path(config_path):
    if config_path is None:
        return os.path.join(SCRIPT_DIR, 'config.py')
    if os.path.isabs(config_path):
        return config_path
    cwd_path = os.path.abspath(config_path)
    if os.path.exists(cwd_path):
        return cwd_path
    script_path = os.path.join(SCRIPT_DIR, config_path)
    if os.path.exists(script_path):
        return script_path
    return cwd_path


def _load_config(config_path):
    """Load a config file as module name 'config' for DisperPicker imports."""
    resolved = _resolve_config_path(config_path)
    if not os.path.exists(resolved):
        raise FileNotFoundError(f'Config file not found: {resolved}')

    config_dir = os.path.dirname(resolved)
    if config_dir not in sys.path:
        sys.path.insert(0, config_dir)

    spec = importlib.util.spec_from_file_location('config', resolved)
    module = importlib.util.module_from_spec(spec)
    sys.modules['config'] = module
    spec.loader.exec_module(module)
    return module.Config(), resolved


def _init_worker(config_path):
    """Build one Pick() (TF model) per worker process."""
    global _pick, _config, _T
    import logging
    _config, _ = _load_config(config_path)
    # Route this worker's stdout/stderr (pick.py prints, library logging, numpy
    # warnings) to its own log file so N workers don't interleave on the
    # terminal. Genuine task errors are still surfaced on the real terminal via
    # sys.__stderr__ in process_one; the main process prints the progress line.
    logdir = os.path.join(_config.result_path, 'logs')
    os.makedirs(logdir, exist_ok=True)
    logf = open(os.path.join(logdir, f'worker_{os.getpid()}.log'), 'a', buffering=1)
    sys.stdout = logf
    sys.stderr = logf
    for h in logging.getLogger('EGFAnalysisTimeFreq').handlers:
        if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler):
            h.stream = logf
    from pick import Pick  # heavy import (tensorflow) kept out of the parent
    _pick = Pick()
    _T = np.linspace(_config.range_T[0], _config.range_T[1], _config.range_T[2])


def _load_label(path_dir, key, n):
    p = path_dir + '/' + key[0] + '.' + key[1] + '.txt'
    if os.path.exists(p):
        return np.loadtxt(p)[:n]
    return np.zeros((1, n))


def process_one(task):
    """Process a single file. Returns (file_name, group_v, phase_v) or None."""
    index, file_num, file_name_full = task
    config, pick, T = _config, _pick, _T
    import EGFAnalysisTimeFreq
    import qc
    from plot.plot_test2 import plot_test

    file_name = os.path.splitext(file_name_full)[0]
    print(f'[{index + 1}/{file_num}] {file_name}', flush=True)
    try:
        batch_group_image = np.zeros((1, config.input_size[0], config.input_size[1]))
        batch_phase_image = np.zeros((1, config.input_size[0], config.input_size[1]))

        batch_group_v = batch_phase_v = None
        if config.test:
            key = file_name.split('_')
            batch_group_v = _load_label(config.test_data_path + '/group_velocity', key, config.input_size[1])
            batch_phase_v = _load_label(config.test_data_path + '/phase_velocity', key, config.input_size[1])
            batch_group_v = np.array([batch_group_v[:, 1]])
            batch_phase_v = np.array([batch_phase_v[:, 1]])

        file_path = os.path.join(config.waveform_path, file_name_full)

        # 1. No time-variable filtering
        gfcn = EGFAnalysisTimeFreq.gfcn_analysis(
            DataFileName=file_path, isEGF=config.isEGF, StartT=config.StartT, EndT=config.EndT,
            DeltaT=config.DeltaT, StartV=config.StartV, EndV=config.EndV, DeltaV=config.DeltaV,
            GreenFcnObjects=config.GreenFcnObjects, WinAlpha=config.WinAlpha,
            NoiseTime=config.NoiseTime, MinSNR=config.MinSNR,
            MinWavelengthRatio=getattr(config, 'MinWavelengthRatio',
                                       getattr(config, 'min_wavelength_ratio', 2.0)))

        sta_info = [gfcn.StaDist, gfcn.Longitude_A, gfcn.Latitude_A, gfcn.Longitude_B, gfcn.Latitude_B]
        if gfcn.StaDist < config.MinStaDist:
            print(f'  Skipped {file_name}: station distance {gfcn.StaDist:.3f} km < {config.MinStaDist} km', flush=True)
            return None

        group_image = gfcn.GroupVelocityImgCalculate()
        group_image = group_image[:config.input_size[0], :config.input_size[1]]
        if True not in np.isnan(group_image):
            batch_group_image[0, config.input_size[0] - group_image.shape[0]:, :group_image.shape[1]] = group_image
        snr = gfcn.SNR_T
        phase_image = gfcn.PhaseVelocityImgCalculate(
            TimeVariableFilter=gfcn.TimeVariableFilterType.no, WinPeriodNum=config.WinPeriodNum,
            WinMinTime=config.WinMinTime, FilterKaiserPara=config.FilterKaiserPara,
            MaxFilterLengthLog=config.MaxFilterLengthLog)
        phase_image = phase_image[:config.input_size[0], :config.input_size[1]]
        if True not in np.isnan(phase_image):
            batch_phase_image[0, config.input_size[0] - phase_image.shape[0]:, :group_image.shape[1]] = phase_image

        pick.mean_confidence_G = 0
        group_velocity, phase_velocity, _, _ = pick.pick(
            group_image=batch_group_image, phase_image=batch_phase_image, sta_info=[sta_info],
            snr=[snr], file_list=[file_name], ct=0.01, save_result=False,
            group_velocity_label=batch_group_v, phase_velocity_label=batch_phase_v)

        # 2. Time-variable filtering
        if np.count_nonzero(group_velocity[0]) > 0:
            gfcn.GroupDisperCurve = group_velocity[0]
            phase_image = gfcn.PhaseVelocityImgCalculate(
                TimeVariableFilter=gfcn.TimeVariableFilterType.obs, WinPeriodNum=config.WinPeriodNum,
                WinMinTime=config.WinMinTime, FilterKaiserPara=config.FilterKaiserPara,
                MaxFilterLengthLog=config.MaxFilterLengthLog)
            phase_image = phase_image[:config.input_size[0], :config.input_size[1]]
            if True not in np.isnan(phase_image):
                batch_phase_image[0, config.input_size[0] - phase_image.shape[0]:, :group_image.shape[1]] = phase_image
        pick.mean_confidence_G = 0.2
        min_wavelength_ratio = getattr(config, 'MinWavelengthRatio',
                                       getattr(config, 'min_wavelength_ratio', 2.0))
        group_velocity, phase_velocity, prob_G, prob_C = pick.pick(
            group_image=batch_group_image, phase_image=batch_phase_image, sta_info=[sta_info],
            snr=[snr], file_list=[file_name], ct=min_wavelength_ratio, save_result=True,
            group_velocity_label=batch_group_v, phase_velocity_label=batch_phase_v)

        # 3. Quality control
        if np.count_nonzero(snr) == 0:
            snr = np.ones(config.input_size[1])
        group_velocity = qc.qc(group_velocity[0], snr,
                               v_range=[config.StartV, config.EndV], diff_range=[-0.07, 0.08],
                               upward=0.4, each_stage_upward=0.3, min_len=config.min_len, skip=False)
        phase_velocity = qc.qc(phase_velocity[0], snr,
                               v_range=[config.StartV, config.EndV], diff_range=[-0.1, 0.1],
                               upward=0, each_stage_upward=0.2, min_len=config.min_len, skip=True)
        format_G = np.array([T, group_velocity, snr, prob_G]).T
        format_C = np.array([T, phase_velocity, snr, prob_C]).T

        # sparse plotting (same gate as run.py)
        if random.random() <= config.random_plot_ratio:
            fig_name = '{}/qc_plot/{}'.format(config.result_path, file_name)
            plot_test(batch_group_image[0], group_velocity, batch_phase_image[0], phase_velocity, fig_name)

        disp_name = '{}/qc_result/GDisp.{}.txt'.format(config.result_path, file_name)
        with open(disp_name, 'w') as f:
            f.write(f'{sta_info[1]:.8f}    {sta_info[2]:.8f}\n')
            f.write(f'{sta_info[3]:.8f}    {sta_info[4]:.8f}\n')
            np.savetxt(f, format_G, fmt="%1.2f  %1.3f  %.3f  %.3f")
        disp_name = '{}/qc_result/CDisp.{}.txt'.format(config.result_path, file_name)
        with open(disp_name, 'w') as f:
            f.write(f'{sta_info[1]:.8f}    {sta_info[2]:.8f}\n')
            f.write(f'{sta_info[3]:.8f}    {sta_info[4]:.8f}\n')
            np.savetxt(f, format_C, fmt="%1.2f  %1.3f  %.3f  %.3f")

        return (file_name, group_velocity, phase_velocity)
    except Exception as e:
        # Short message to the real terminal; full traceback to the worker log.
        print(f'  ERROR on {file_name}: {e}', file=sys.__stderr__, flush=True)
        traceback.print_exc(file=sys.stderr)
        return None


def parse_args():
    parser = argparse.ArgumentParser(
        description='Parallel driver for EGFAnalysisPy dispersion picking.')
    parser.add_argument(
        '-w', '--workers', type=int,
        default=int(os.environ.get('NWORKERS', min(16, mp.cpu_count()))),
        help='Number of worker processes (default: $NWORKERS or min(16, CPU count)).')
    parser.add_argument(
        '-n', '--max-files', type=int,
        default=int(os.environ.get('MAXFILES', 0)),
        help='Process only the first N files; 0 means all '
             '(default: $MAXFILES or 0). Useful for quick tests.')
    parser.add_argument(
        '-c', '--config',
        default=os.environ.get('CONFIG', os.path.join(SCRIPT_DIR, 'config.py')),
        help='Path to config.py-compatible file '
             '(default: $CONFIG or config.py next to run_parallel.py).')
    return parser.parse_args()


def main():
    args = parse_args()
    config, config_path = _load_config(args.config)
    file_list = sorted(os.listdir(config.waveform_path))
    if args.max_files > 0:
        file_list = file_list[:args.max_files]
    file_num = len(file_list)
    nworkers = args.workers
    print(f'\nConfig: {config_path}', flush=True)
    print(f'File number: {file_num}   Workers: {nworkers}\n', flush=True)

    # Pre-create output dirs once, in the parent, to avoid worker race conditions.
    for d in ('qc_plot', 'qc_result', 'pick_result', 'plot', 'logs'):
        os.makedirs(f'{config.result_path}/{d}', exist_ok=True)
    print(f'(per-worker output: {config.result_path}/logs/worker_*.log)\n', flush=True)

    tasks = [(i, file_num, file_list[i]) for i in range(file_num)]
    ctx = mp.get_context('spawn')
    results = []
    skipped = 0
    with ctx.Pool(processes=nworkers, initializer=_init_worker,
                  initargs=(config_path,)) as pool:
        # imap_unordered streams results as they finish -> smooth progress.
        # Order does not matter for the statistics (mean/overlay are order-invariant).
        pbar = tqdm(pool.imap_unordered(process_one, tasks, chunksize=1),
                    total=file_num, desc='Picking', unit='file', dynamic_ncols=True)
        for r in pbar:
            if r is None:
                skipped += 1
            else:
                results.append(r)
            pbar.set_postfix(ok=len(results), skipped=skipped)

    # 4. Statistics (drop skipped/failed)
    import qc
    all_G = np.array([r[1] for r in results])
    all_C = np.array([r[2] for r in results])
    print(f'\nProcessed {len(results)}/{file_num} files ({skipped} skipped/failed)', flush=True)

    # qc.plot_disp / qc.mean dump large arrays to stdout; route that chatter to a
    # log file so the terminal summary stays clean. The PNGs are still written.
    stats_log = os.path.join(config.result_path, 'logs', 'statistics.log')
    saved_stdout = sys.stdout
    with open(stats_log, 'w') as lf:
        sys.stdout = lf
        try:
            qc.plot_disp(all_C, config.range_T, config.result_path + '/all_C.png')
            qc.mean(all_C, config.range_T, config.result_path + '/mean_C.png')
            qc.plot_disp(all_G, config.range_T, config.result_path + '/all_G.png')
            qc.mean(all_G, config.range_T, config.result_path + '/mean_G.png')
        finally:
            sys.stdout = saved_stdout

    print(f'Statistics figures written to {config.result_path}/ '
          '(all_C.png, mean_C.png, all_G.png, mean_G.png)', flush=True)


if __name__ == '__main__':
    main()
