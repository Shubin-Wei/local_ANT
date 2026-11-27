# local_ANT

**Ambient-noise surface-wave tomography for dense arrays: from continuous records to a 3-D shear-velocity model.**

[![tests](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml/badge.svg)](https://github.com/Shubin-Wei/local_ANT/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-blue)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

[中文](README.zh.md) | English

lant (`local_ANT`) is a command-line tool for high-frequency ambient-noise imaging with dense (nodal) arrays, areal or linear. It does the basic processing chain of the literature, nothing more:

- **General**: any record format ObsPy reads, any folder layout, any array; no function made for one data set.
- **Easy**: one configuration file with a comment on every setting; `lant init` suggests the numbers for your array.
- **Transparent**: every step writes the settings it used, every value it derived with the rule that gave it, its log and its figures. Intermediate results are plain HDF5 and CSV files; the inversion project is in DSurfTomo's format.
- **Simple**: three steps, run and checked one at a time.
- **Correct**: every method has a source ([methods](docs/methods.md)); unit tests check the formulas against exact answers.

## Steps

| Command | What it does | Basis |
|---|---|---|
| `lant cc` | reads the records, correlates every station pair (cross-coherence of overlapping segments), stacks | Bensen et al. (2007); Seats et al. (2012); Xu et al. (2021) |
| `lant disp` | reference curve (array F-J spectrum, a 1-D model, or a curve); two-station phase velocities with phase-matched filtering, jackknife errors and quality control | Wang et al. (2019); Levshin & Ritzwoller (2001); Lin et al. (2008); Luo et al. (2015) |
| `lant inv` | 1-D starting model; direct 3-D inversion of the travel times with ray tracing | Fang et al. (2015) (DSurfTomo's method, in Python) |
| `lant checker` | checkerboard resolution test | Lévêque et al. (1993); Rawlinson & Spakman (2016) |

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

```bash
lant init -o my.yaml --raw /data/raw --stations stations.txt --out /data/out --velocity 0.3 3.0
# read my.yaml, then:
lant cc   -c my.yaml      # check 1_cc/figures, then
lant disp -c my.yaml      # check 2_disp/figures/fj_reference.png above all, then
lant inv  -c my.yaml
lant checker -c my.yaml --size 2
```

The station table needs a header line: `network station longitude latitude elevation` (the network column may be left out). `--velocity` is the slowest and fastest Rayleigh phase velocity you expect in the area.

What to look at after each step: [tutorial](docs/tutorial.md).

## Output

```
OUT/1_cc/       stack.h5, days/*.h5, stations_used.csv, record_index.csv
OUT/2_disp/     dispersion.csv, measurements.h5, reference.csv, median_vs_reference.csv, fj.h5
OUT/3_inv/      vs_model.csv, start_model.csv, residuals.csv, dsurftomo/
OUT/4_checker_*/recovery.csv
each step also: params.yaml, log.txt, summary.txt, figures/
```

## Documentation

- [Tutorial](docs/tutorial.md): how to run each step and what to check
- [Methods](docs/methods.md): principles, formulas, sources and limits
- [Parameters](docs/parameters.md): every setting, its default and its basis
- [History](docs/history.md)

## Citation

Please cite the methods you use (listed in [methods](docs/methods.md)), in particular Fang et al. (2015) for the 3-D inversion, and this software ([CITATION.cff](CITATION.cff)).

## License

MIT. The inversion (`lant/inversion/dsurf_*.py`) is a Python port of DSurfTomo (Fang et al. 2015; MIT licence).
