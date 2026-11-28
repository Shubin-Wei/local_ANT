# local_ANT

**Ambient-noise surface-wave tomography for dense arrays: from continuous records to a 3-D shear-velocity model.**

[![tests](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml/badge.svg)](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

[中文](README.zh.md) | English

lant (`local_ANT`) is a command-line tool for high-frequency ambient-noise imaging with dense (nodal) arrays, areal or linear. It takes continuous vertical-component records and a station table, and gives Rayleigh-wave phase velocities between station pairs and a 3-D Vs model. The steps are run one at a time; each writes the settings it used, the values it derived (with their rules), a log and figures to check before the next step.

## Steps

| Command | What it does | Main references | Results |
|---|---|---|---|
| `lant clean` | cuts the records into UTC days per station; removes the median and, if given, the instrument response; low-pass and decimation | Bensen et al. (2007) | `1_clean` |
| `lant cc` | segments (those with gaps or large amplitudes left out), cross-coherence of every station pair, linear stack | Bensen et al. (2007); Seats et al. (2012); Xu et al. (2021) | `2_cc` |
| `lant disp` | reference curve (array F-J spectrum, a 1-D model or a curve file); two-station phase velocities with jackknife errors and quality control | Wang et al. (2019); Levshin & Ritzwoller (2001); Lin et al. (2008) | `3_disp` |
| `lant inv` | 1-D starting model; direct 3-D inversion of the travel times with ray tracing | Fang et al. (2015) (DSurfTomo, in Python) | `4_inv` |
| `lant checker` | checkerboard resolution test | Rawlinson & Spakman (2016) | `5_checker_*` |

## Install

Linux, with conda (or mamba) and git:

```bash
git clone https://github.com/Shubin-Wei/local_ANT.git
cd local_ANT
conda env create -f environment.yml
conda activate lant
lant check
```

GMT 6 (installed by the environment) draws the maps; without it everything else works.

## Quick start

You need:

- continuous records in any format ObsPy reads (miniSEED, SAC, ...), in any folder layout;
- a station table with a header line, `network station longitude latitude elevation` (the network column may be left out), or StationXML; the names must be those in the record headers;
- the slowest and fastest Rayleigh phase velocity you expect in the area (km/s).

```bash
lant init -o my.yaml --raw /data/raw --stations stations.txt --out /data/out --velocity 0.3 3.0
# read my.yaml (lant init wrote suggested values and their rules into it), then:
lant clean -c my.yaml     # check 1_clean/figures
lant cc    -c my.yaml     # check 2_cc/figures
lant disp  -c my.yaml     # check 3_disp/figures/fj_reference.png above all
lant inv   -c my.yaml
lant checker -c my.yaml --size 2
```

What to look at after each step: [tutorial](docs/tutorial.md).

## Results

```
OUT/1_clean/       days/*.h5 (cleaned station-days), availability.csv, stations_used.csv, record_index.csv
OUT/2_cc/          stack.h5, days/*.h5 (sub-stacks)
OUT/3_disp/        dispersion.csv, measurements.h5, reference.csv, median_vs_reference.csv, fj.h5
OUT/4_inv/         vs_model.csv, start_model.csv, residuals.csv, dsurftomo/
OUT/5_checker_*/   recovery.csv
every step also:   params.yaml, log.txt, summary.txt, figures/
```

Large arrays are HDF5, everything else CSV; `lant export-sac` writes the stacked correlations as SAC files. `4_inv/dsurftomo/` is a project in DSurfTomo's format.

## Documentation

- [Tutorial](docs/tutorial.md): how to run each step and what to check
- [Methods](docs/methods.md): principles, formulas, sources and limits
- [Parameters](docs/parameters.md): every setting, its default and its basis
- [History](docs/history.md)

## Citation

Please cite the methods you use (listed in [methods](docs/methods.md)), in particular Fang et al. (2015) for the 3-D inversion, and this software ([CITATION.cff](CITATION.cff)).

## License

MIT. The inversion (`lant/inversion/dsurf_*.py`) is a Python port of DSurfTomo (Fang et al. 2015; MIT licence).
