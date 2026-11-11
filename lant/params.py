"""What every configuration key means and what its value rests on (one table for the step reports, `lant plan`
and the configuration written by `lant init`).

DOC[key] = (kind, text).  kind says what the default (or the rule behind "auto") rests on:
    required    the user's knowledge; there is no default
    literature  a value or a method taken from the cited paper or program
    theory      follows from a stated relation (sampling theorem, wavelength, filter width, ...)
    rule        "auto": derived from the data by the stated rule; the rule's coefficients are empirical unless a
                source is cited
    empirical   chosen by experience on test data; no literature value was found (change it with care and say so)
    choice      a choice of the user without a right value (what to draw, what to leave out, where to write)
    technical   speed, memory or file handling only; the results do not depend on it
The derivation of an "auto" value with the numbers of the present data is written by the rule itself
(.lant/auto_parameters.json) and shown in the step reports.
"""

DOC = {
    # ------------------------------------------------------------------ project
    "project.name": ("choice", "project name (figures and records only)"),
    "project.array_type": (
        "rule",
        "auto: from the station geometry (how many times the station cloud is longer than "
        "wide; below 3 areal, above 10 line, in between shown for confirmation); or "
        "surface / line",
    ),
    "project.output_dir": ("required", "every result, figure and record is written here"),
    # ------------------------------------------------------------------ data
    "data.raw_dir": ("required", "raw-data directory, searched recursively; any layout and file names"),
    "data.format": (
        "technical",
        "auto: ObsPy detects the format of each file (MSEED, SAC, GSE2, SEG2, ...); or one ObsPy format name",
    ),
    "data.file_patterns": ("technical", "file-name patterns that are tried"),
    "data.ignore_patterns": ("technical", "file-name patterns that are not read"),
    "data.station_file": ("required", "station table in any common form (CSV, text, Excel, StationXML)"),
    "data.component": ("choice", "Z: vertical-component correlations (Rayleigh waves); the only one supported"),
    "data.channel": (
        "rule",
        "auto: the vertical channel of each station from the channel name or the SAC header; or a "
        'pattern such as "*3" when the names do not tell',
    ),
    "data.network_default": ("choice", "network code for record headers without one"),
    "data.rename": ("choice", "{header name: table name} where the automatic matching fails"),
    "data.day_start_utc_hour": (
        "rule",
        "processing days start at this UTC hour; auto: the hour at which >= 60 % of the "
        "records longer than 12 h start, else 0",
    ),
    "data.day_label": ("choice", 'how a processing day is labelled: utc_start, or a fixed offset such as "UTC+8"'),
    "data.velocity_range_kms": (
        "choice",
        "approximate Rayleigh-wave phase velocities of the area [slow, fast] in km/s; "
        "sets the lag window, the band and the F-J scan.  The default covers soft "
        "soil to hard rock; a narrower range saves computing time (the reference "
        "step reports the measured velocities)",
    ),
    "data.response": (
        "choice",
        "instrument response file or directory (StationXML, dataless SEED, RESP): the "
        "responses are removed in the clean step (Bensen et al. 2007).  Without it every "
        "station must carry the same instrument model",
    ),
    "data.line": (
        "choice",
        "when the station table has a line column: the line this project processes (one project per line)",
    ),
    "data.min_trace_length_s": ("empirical", "pieces of record shorter than this are not used"),
    "data.exclude_stations": ("choice", "stations left out entirely (an instrument or coupling problem you know of)"),
    "data.start_date": ("choice", "first processing day used (YYYY-MM-DD)"),
    "data.end_date": ("choice", "last processing day used (YYYY-MM-DD)"),
    # ------------------------------------------------------------------ clean
    "clean.write_archive": (
        "technical",
        "write the cleaned day files (about the size of the raw data); false: later "
        "steps read the raw files, and the 1-s periodic signal is not tested "
        "(response removal and moved nodes need the archive)",
    ),
    "clean.nproc": ("technical", "processes"),
    "clean.encoding": ("technical", "miniSEED encoding of the archive"),
    "clean.min_free_disk_factor": (
        "technical",
        "refuse to write the archive when the free disk is below this factor x its estimated size",
    ),
    # ------------------------------------------------------------------ correlation
    "correlation.fs_max": (
        "choice",
        "highest processing rate (Hz); the top of the band is 0.4 fs / (1 + sqrt(ln 2 / "
        "alpha)): raise it for small, high-frequency arrays",
    ),
    "correlation.fs": (
        "rule",
        "rate after decimation; auto: the highest whole-number rate <= fs_max that every raw rate "
        "reaches by an integer factor",
    ),
    "correlation.seg_len_s": (
        "rule",
        "segment length; auto: max(2 x maxlag (circular correlation), 10 longest periods)",
    ),
    "correlation.seg_step_s": ("literature", "segment step; auto: half the segment (50 % overlap, Seats et al. 2012)"),
    "correlation.maxlag_s": (
        "rule",
        "lag window; auto: longest distance / (0.6 x slowest velocity) + 4 longest periods "
        "(the slowest arrival, then a stretch of noise for the SNR)",
    ),
    "correlation.fmin_hz": (
        "rule",
        "bottom of the correlation band; auto: where the longest pair is 2 wavelengths at "
        "the slowest velocity, widened to the half-gain point of the Gaussian filter there",
    ),
    "correlation.fmax_hz": (
        "theory",
        "top of the correlation band; auto: dispersion.fmax_hz widened to the half-gain "
        "point of its Gaussian filter, not above 0.4 fs (anti-alias filter)",
    ),
    "correlation.smooth_bins": (
        "literature",
        "smoothing of the coherence denominator in frequency bins; 1 = none (cross-coherence as in Xu et al. 2021)",
    ),
    "correlation.max_distance_km": (
        "choice",
        "only pairs up to this distance are correlated (none: all pairs); for "
        "very large arrays, where all pairs do not fit into memory",
    ),
    "correlation.max_over_std": (
        "literature",
        "a segment is dropped when its peak exceeds this many standard deviations "
        "of the station-day (NoisePy max_over_std, Jiang & Denolle 2020)",
    ),
    "correlation.remove_1s_comb": (
        "literature",
        "subtract the waveform folded modulo 1 s, hour by hour (coherent logger "
        "noise at 1-s intervals, Takagi et al. 2015); auto: for the stations "
        "where the clean step detects it (F test of the folded records); "
        "true: all stations; false: none",
    ),
    "correlation.nproc": (
        "technical",
        "days processed at the same time (memory: one day of every station per process)",
    ),
    "correlation.threads_per_process": ("technical", "FFT threads per process"),
    # ------------------------------------------------------------------ stack
    "stack.block_hours": (
        "rule",
        "sub-stack length of the jackknife; auto: the longest of 24, 12, 8, 6, 4, 3, 2, 1 h "
        "that gives the median station at least 8 sub-stacks",
    ),
    # (see also correlation.max_distance_km)
    "stack.station_qc.enabled": (
        "choice",
        "station diagnosis: sub-stacks of a station that differ from the rest in "
        "correlation or time shift are listed (never removed without you)",
    ),
    "stack.station_qc.partners": (
        "empirical",
        "each station is compared through its pairs with this many nearest stations",
    ),
    "stack.station_qc.min_partners": (
        "empirical",
        "pairs with data needed to judge a station sub-stack (fewer: not judged, and reported so)",
    ),
    "stack.exclude_station_days": ("choice", 'station-days you remove: [station, YYYYDDD], [station, "YYYYDDD-"], ...'),
    "stack.exclude_station_hours": ("choice", "station-hours you remove: [station, YYYYDDD, hour]"),
    "stack.exclude_candidates": (
        "choice",
        "also remove the candidates of the diagnosis (set by `lant confirm stack --exclude-candidates`)",
    ),
    # ------------------------------------------------------------------ reference
    "reference.source": (
        "literature",
        "fj: frequency-Bessel spectrum of this array (Wang et al. 2019); file: a curve "
        "you give (as AFTAN and Yao's EGFAnalysisTimeFreq take a reference curve)",
    ),
    "reference.file": ("choice", 'the borrowed curve: npz written by lant, or text "freq_Hz c_km/s"'),
    "reference.zone_km": (
        "choice",
        "none: one reference curve for the whole array; a size in km: one F-J curve per "
        "square zone of that size from the pairs whose path midpoint lies in it, for "
        "arrays across which the structure changes strongly; a zone "
        "without a clear curve uses that of the whole array",
    ),
    "reference.fj_mode": (
        "empirical",
        "max_peak: the strongest ridge (areal arrays); track_fundamental: the ridge "
        "followed from an unambiguous start (lines)",
    ),
    "reference.normal_dispersion": (
        "choice",
        "track_fundamental: the curve may not rise with frequency; false where a "
        "velocity reversal makes the ridge really rise",
    ),
    "reference.fj_rmax_lambda": (
        "empirical",
        "pairs up to this many wavelengths enter the F-J spectrum of a frequency",
    ),
    "reference.fj_c_nominal": (
        "rule",
        "velocity that turns fj_rmax_lambda into a distance; auto: the curve picked "
        "in a first pass with every pair (two passes)",
    ),
    "reference.fj_cmin": ("rule", "slowest scanned velocity; auto: 0.5 x the slowest of data.velocity_range_kms"),
    "reference.fj_cmax": ("rule", "fastest scanned velocity; auto: 1.5 x the fastest of data.velocity_range_kms"),
    "reference.fj_dc": ("rule", "velocity step of the scan; auto: the 1-2-2.5-5 step not above fj_cmin / 200"),
    "reference.fj_nfft": ("technical", "FFT length of the spectra; auto: power of 2 >= 4 x the lag window, >= 8192"),
    "reference.fj_fit_fmax": (
        "rule",
        "highest frequency of the fitted picks; auto: the top of the longest run of clear picks",
    ),
    "reference.track_start_band_hz": (
        "rule",
        "where the ridge tracking starts; auto: of starts every 1/3 octave, the one followed over the most "
        "frequencies (on a tie the slower)",
    ),
    # ------------------------------------------------------------------ dispersion
    "dispersion.fmin_hz": (
        "rule",
        "lowest measured frequency; auto: where the longest pair is qc.lam_c wavelengths of the reference curve",
    ),
    "dispersion.fmax_hz": (
        "theory",
        "highest measured frequency; auto: 0.4 fs / (1 + sqrt(ln 2 / alpha)), the highest "
        "whose Gaussian filter stays inside the anti-alias band; the QC decides how far "
        "the data reach",
    ),
    "dispersion.per_octave": ("literature", "measurement frequencies per octave, log-spaced (Xu et al. 2021)"),
    "dispersion.vmin": ("rule", "bottom of the velocity axis; auto: 0.6 x the slowest reference velocity in the band"),
    "dispersion.vmax": ("rule", "top of the velocity axis; auto: 1.3 x the fastest"),
    "dispersion.dv": ("rule", "step of the velocity axis; auto: the 1-2-2.5-5 step not above vmin / 100"),
    "dispersion.alpha": (
        "literature",
        "relative width of the Gaussian filter exp(-alpha ((f - fc) / fc)^2), as in FTAN "
        "(Bensen et al. 2007); larger = narrower",
    ),
    "dispersion.phase_constant": (
        "theory",
        "phase of the far-field 2-D Green's function, pi/4 (Lin et al. 2008 eq. 7)",
    ),
    "dispersion.min_dist_km": ("choice", "pairs shorter than this are not measured (0: lam_min alone rules)"),
    "dispersion.lam_min": (
        "literature",
        "paths shorter than this many reference wavelengths are not measured (Luo et "
        "al. 2015: one wavelength can be enough)",
    ),
    "dispersion.anchor_max_lambda": (
        "empirical",
        "the cycle count is fixed where the path is lam_min to this many "
        "wavelengths long (few cycles: the reference cannot be a cycle off)",
    ),
    "dispersion.njack": (
        "empirical",
        "delete-a-group jackknife groups of the sub-stacks; auto: the sub-stacks holding >= 25 % of the fullest one",
    ),
    "dispersion.min_jackknife_groups": (
        "empirical",
        "fewer groups with data: no error estimate, the sample is rejected",
    ),
    "dispersion.nproc": ("technical", "processes"),
    "dispersion.chunk": ("technical", "pairs per task"),
    "dispersion.side_off_bands": ("choice", "bands 'lo-hi,...' where the side-lobe rule is not applied"),
    "dispersion.qc.snr_c": (
        "literature",
        "SNR at the picked arrival (peak / rms of the trailing noise); the literature "
        "uses 4-10 (Bensen et al. 2007: 10 for global data)",
    ),
    "dispersion.qc.jk_c_rel": ("empirical", "jackknife standard error allowed, relative to the velocity"),
    "dispersion.qc.ucons_c": (
        "empirical",
        "disagreement allowed between the measured group velocity and the one implied by the phase curve (0 = off)",
    ),
    "dispersion.qc.side_c": ("empirical", "height of a competing envelope peak allowed, relative to the pick"),
    "dispersion.qc.dev_pre": ("empirical", "deviation of the group pick allowed from the reference (0 = off)"),
    "dispersion.qc.lam_c": (
        "literature",
        "accepted samples need a path of this many of their own wavelengths (Bensen "
        "et al. 2007 use 3; Luo et al. 2015 show 1 can be enough)",
    ),
    "dispersion.qc.band_edge": (
        "choice",
        "a frequency is not measured where the amplitude spectrum of the "
        "coherent signal (all pairs together; signal window minus the later lags, so "
        "that the noise level does not matter) is below this fraction of its mean "
        "over the filter band (0.5 = -6 dB, the usual definition of a band edge): "
        "beyond the band the noise excites, the filter passes only energy of other "
        "frequencies (0 = not applied)",
    ),
    "dispersion.qc.min_run": ("empirical", "accepted samples must form runs of at least this many frequencies"),
    "dispersion.qc.snr_g": ("empirical", "first gate: the largest SNR of a frequency must reach this"),
    # ------------------------------------------------------------------ export
    "export.tmin_s": ("rule", "shortest period passed on; auto: 1 / dispersion.fmax_hz"),
    "export.tmax_s": ("rule", "longest period passed on; auto: 1 / dispersion.fmin_hz"),
    "export.one_segment": (
        "literature",
        "one continuous segment per pair, the longest (as curves are picked in Yao's EGFAnalysisTimeFreq)",
    ),
    "export.curve_figures": (
        "choice",
        "single-pair figures drawn: this many spread over distance plus the conspicuous ones (0 none, -1 all)",
    ),
    # ------------------------------------------------------------------ inversion
    "inversion.wave": ("choice", "rayleigh_phase (the only one supported)"),
    "inversion.median_min_pairs": ("empirical", "pairs a period needs to enter the average curve of the 1-D model"),
    "inversion.depth_sequence_km": (
        "rule",
        "depth nodes; auto: steps of max(shortest wavelength / 4, depth / 4) down "
        "to half the longest wavelength (rule of thumb, Foti et al. 2018)",
    ),
    "inversion.fine_step_km": ("choice", "extra depth nodes every this many km down to fine_zmax_km (0 = none)"),
    "inversion.fine_zmax_km": ("choice", "see fine_step_km"),
    "inversion.grid_spacing_deg": (
        "rule",
        "grid spacing [lat, lon]; auto: about half the median nearest-station distance after refinement",
    ),
    "inversion.grid_refine": ("empirical", "every grid interval is divided by this once the extent is fixed"),
    "inversion.line_strip_cells": ("empirical", "lines: cells kept on each side of the stations across the line"),
    "inversion.start_model.source": ("choice", "invert: 1-D model from this array's average curve; file: depth Vs"),
    "inversion.start_model.file": ("choice", "the starting model of your own (two columns: depth km, Vs km/s)"),
    "inversion.start_model.smooth": ("empirical", "1-D inversion: weight of the model curvature"),
    "inversion.start_model.damp": ("empirical", "1-D inversion: damping"),
    "inversion.start_model.iterations": ("empirical", "1-D inversion: iterations"),
    "inversion.start_model.minthk": ("technical", "1-D forward: sublayers per node interval"),
    "inversion.weight": (
        "empirical",
        "DSurfTomo smoothing weight (Fang et al. 2015); set by hand in the usual range "
        "(below 10, mostly below 5), not optimised; lines 1.5",
    ),
    "inversion.damp": (
        "empirical",
        "DSurfTomo damping (Fang et al. 2015); set by hand in the usual range 0.5-1.5, not optimised; lines 1.5",
    ),
    "inversion.sublayers": ("literature", "DSurfTomo: sublayers between depth nodes for the kernels (its manual: 2-5)"),
    "inversion.vmin": ("rule", "lower Vs bound; auto: 0.6 x the slowest Vs of the starting model"),
    "inversion.vmax": ("rule", "upper Vs bound; auto: 1.4 x the fastest"),
    "inversion.iterations": ("empirical", "DSurfTomo iterations"),
    "inversion.sparsity": (
        "rule",
        "DSurfTomo: fraction of its matrix reserved for nonzero entries; auto: twice the "
        "filling estimated from the paths and the grid, at least 0.05 (too little makes "
        "DSurfTomo stop with a segmentation fault; more only costs memory)",
    ),
    "inversion.threshold": ("empirical", "DSurfTomo: residuals beyond this many quartiles get zero weight"),
    "inversion.noise": (
        "rule",
        "Gaussian noise added to the synthetic travel times of the checkerboards, a fraction of "
        "each time; auto: the relative rms of the residuals of the final main model that "
        "DSurfTomo's rejection rule keeps, so that the synthetic data are as noisy as the real "
        "data fit (Rawlinson & Spakman 2016); the checkerboards then run after the main project",
    ),
    "inversion.solver": ("choice", "dsurftomo (<= 60 periods) or pdsurftomo (parallel, <= 80 periods)"),
    "inversion.threads": ("technical", "threads per project"),
    "inversion.njobs": ("technical", "projects run at the same time"),
    # ------------------------------------------------------------------ checks
    "checks.checkerboards": (
        "rule",
        "horizontal checkerboards [{name, cell_km}]; auto: cells of 0.12, 0.18 and 0.25 x "
        "the aperture (lines: 0.165, 0.24, 0.34 x the line length, along the line)",
    ),
    "checks.vertical": (
        "rule",
        "depth-varying checkerboards; auto: cell 0.18 x aperture, sign reversals at 0.21, 0.46, "
        "0.71 x the deepest updated node",
    ),
    "checks.amplitude": ("empirical", "checkerboard perturbation, a fraction of the 1-D Vs of each depth"),
    "checks.controls": (
        "empirical",
        "control inversions: reg_w<w>_d<d> other regularisation, start_grad<f> other starting model",
    ),
    "checks.precheck.min_kept": (
        "empirical",
        "a project is invalid when DSurfTomo's rejection rule would keep less "
        "than this fraction of the data in the first iteration",
    ),
    "checks.precheck.max_median_iqr": (
        "empirical",
        "... or when the median residual is beyond this many interquartile ranges",
    ),
    "checks.iteration_residuals": ("choice", "also evaluate the model of every iteration (one forward run each)"),
    # ------------------------------------------------------------------ plots
    "plots.formats": ("choice", "figure formats"),
    "plots.dpi": ("choice", "resolution of the png figures"),
    "plots.topography": (
        "choice",
        "relief under the maps; auto: SRTM 1-arc-second tiles (Farr et al. 2007) from the "
        "GMT data server, kept in ~/.cache/lant/srtm1 (only the tile names are requested; "
        "maps up to 3 degrees wide); or a grid file / GMT remote grid; null: no relief",
    ),
    "plots.vs_range": ("rule", "colour range; auto: 2nd-98th percentile of the model under the stations, to 0.1 km/s"),
    "plots.vs_cpt": ("choice", "GMT colour table"),
    "plots.slice_depths_km": ("rule", "depths of the slices; auto: up to 8 updated depth nodes, evenly chosen"),
    "plots.sections": (
        "rule",
        "sections [{name, start, end}]; auto: two through the array centre along its principal "
        "axes (areal) or along the line",
    ),
    "plots.section_zmax_km": ("choice", "depth of the sections (default: the deepest updated node)"),
    "plots.contour_step": ("rule", "contour interval; auto: about a tenth of the colour range"),
    "plots.font_scale": ("choice", "font size factor"),
    "plots.pair_example": (
        "choice",
        "pair(s) of the single-pair figures (stack time, dispersion): a name or a list; "
        "default: three pairs at the 20 / 50 / 80 % quantiles of the distances of the "
        "pairs with long curves",
    ),
}

KINDS = ("required", "literature", "theory", "rule", "empirical", "choice", "technical")


def flat(d, prefix="", free=("data.rename",)):
    """[(dotted key, value)] of a configuration mapping, in its order (free mappings are values)."""
    out = []
    for k, v in d.items():
        key = prefix + k
        if isinstance(v, dict) and v and key not in free and not key.startswith("_"):
            out += flat(v, key + ".", free)
        elif not key.startswith("_"):
            out.append((key, v))
    return out


def describe(key):
    """(kind, text) of a key (also for a key below a documented one)."""
    k = key
    while k:
        if k in DOC:
            return DOC[k]
        k = k.rpartition(".")[0]
    return ("choice", "")


def _full_template(
    raw_dir="/path/to/raw_data", station_file="/path/to/stations.csv", output_dir="/path/to/output", line=None
):
    """Text of a full configuration file for `lant init`: the three required paths, then every other key with its
    default, commented out, with what it means and what its value rests on."""
    import yaml

    from . import config as C

    cfg = C.defaults("surface")
    cfg["project"]["array_type"] = "auto"
    req = {".".join(k) for k in C.REQUIRED}
    dump = lambda v: yaml.safe_dump(v, default_flow_style=True, width=10**6).strip().replace("\n...", "")
    fill = {"project.output_dir": output_dir, "data.raw_dir": raw_dir, "data.station_file": station_file}
    L = [
        "# local_ANT configuration, written by `lant init`.",
        "#",
        "# 1. Fill in the three REQUIRED paths below.",
        "# 2. `lant plan -c THIS_FILE` shows what every step will do with which values.",
        "# 3. `lant run all -c THIS_FILE` runs the chain; it stops where a judgement is yours and says what to do.",
        "#",
        '# Every other key is listed further down with its default, commented out.  "auto" means: derived from',
        "# your data by the rule in the comment (the step reports show the derivation with your numbers); a value",
        "# you write here always wins.  In square brackets: what the default rests on -",
        "#   literature (cited source) | theory | rule (derived from the data) | empirical (experience, no literature",
        "#   value) | choice (yours, no right value) | technical (speed and files only)",
        "",
        "# ---------------------------------------------------------------- REQUIRED",
    ]
    by_sec = {}
    for key, v in flat(cfg):
        by_sec.setdefault(key.split(".")[0], []).append((key, v))
    for sec in ("project", "data"):
        L.append("%s:" % sec)
        if sec == "data" and line:
            L.append("  line: %s   # %s" % (dump(str(line)), DOC["data.line"][1]))
        for key, v in by_sec[sec]:
            if key in req:
                L.append("  %s: %s%s# %s" % (key.split(".", 1)[1], dump(fill[key]), " " * 3, DOC[key][1]))
    L += [
        "",
        "# ---------------------------------------------------------------- everything else (defaults)",
        "# To change a value: uncomment its line AND the section line above it (keep the indentation); keys of",
        "# `project` and `data` go into the sections above.",
        "# Line arrays (project.array_type line, or auto on an elongated array) use other defaults for:",
        "#   " + ", ".join("%s: %s" % (k, dump(v)) for k, v in flat(C.LINE_OVERRIDES)),
    ]
    for sec, items in by_sec.items():
        L.append("")
        L.append("# %s:" % sec)
        open_sub = None
        for key, v in items:
            if key in req:
                continue
            parts = key.split(".")[1:]
            if len(parts) > 1:
                if open_sub != parts[0]:
                    L.append("#   %s:" % parts[0])
                    open_sub = parts[0]
                ind = "    "
            else:
                open_sub, ind = None, "  "
            kind_, text = describe(key)
            L.append("# %s# [%s] %s" % (ind, kind_, text))
            L.append("# %s%s: %s" % (ind, parts[-1], dump(v)))
    return "\n".join(L) + "\n"


def template(
    raw_dir="/path/to/raw_data",
    station_file="/path/to/stations.csv",
    output_dir="/path/to/output",
    line=None,
    full=False,
):
    """A compact starting point; the complete template uses the same defaults and descriptions."""
    if full:
        return _full_template(raw_dir, station_file, output_dir, line)
    import yaml

    from . import config as C

    cfg = C.defaults()
    dump = lambda v: yaml.safe_dump(v, default_flow_style=True).strip().replace("\n...", "")
    lines = [
        "# lant configuration: three required paths; uncomment optional settings as needed.",
        "# Full template: lant init --full -o full.yaml. Effective values: lant plan -c THIS_FILE.",
        "project:",
        "  output_dir: " + dump(output_dir),
        "  # array_type: auto  # derive surface/line from station geometry",
        "data:",
        "  raw_dir: " + dump(raw_dir),
        "  station_file: " + dump(station_file),
    ]
    if line:
        lines.append("  line: " + dump(str(line)))
    for key in ("velocity_range_kms", "day_label", "response", "exclude_stations", "start_date", "end_date"):
        lines.append("  # %s: %s  # %s" % (key, dump(cfg["data"][key]), describe("data." + key)[1]))
    lines += ["# plots:", "#   formats: " + dump(cfg["plots"]["formats"]) + "  # add pdf for vector output"]
    return "\n".join(lines) + "\n"
