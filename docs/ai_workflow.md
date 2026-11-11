# Running with an AI assistant

[中文](ai_workflow.zh.md) | English

This file has two readers: users who want an AI assistant (Claude Code, Codex and the like) to run the chain for them, and the assistant asked to do so.

## 1. How the user starts

Open the assistant in the repository directory and tell it three things:

```
Process my data with this project. Raw data in /data/raw, station table /data/stations.csv, results to /data/result.
```

If you know them, add: the rough range of surface-wave velocities, whether there are response files, whether all instruments are one model, which stations had problems in the field. The assistant does the rest by this file and stops to ask where a judgement is yours.

## 2. Rules for the assistant

1. **Use only the commands of the project** (`lant …`). Do not write processing code of your own and do not edit `lant/` to suit the data; when you find a defect in the program, stop and tell the user.
2. **Do not change the raw data; do not write results into the repository.** All output goes to the user's output directory.
3. **Judgements belong to the user.** When `lant run all` stops for a judgement (exit code 3), give the user the items under "FOR YOU TO JUDGE" in `records/steps/<step>.md` and the figures to look at, and wait. **Never run `lant confirm` on your own**, unless the user explicitly said "confirm this step if you find nothing wrong".
4. **Do not guess what only the user knows**: instrument models, what happened to a station in the field, the time zone, whose reference curve to borrow.
5. **Do not tune for a nicer result**: no looser quality control for more points, no regularisation chosen by checkerboard or residual, no removal of stations that "look bad".
6. Every parameter that differs from its default needs a basis (data statistics, a figure, literature), recorded in `records/decisions.md`.
7. **Report numbers and figures only**, without words such as "stable", "reliable", "passed"; say what was not checked and what the data limit.

## 3. Who decides what

| | Items |
|---|---|
| The user gives | the three paths; if known: velocity range, response files, time zone, stations to exclude |
| The program derives, the report shows the derivation | sampling rate, band, lag window, segment length, sub-stack length, F-J scan, measurement band and velocity axis, depth nodes, grid, velocity bounds, checkerboards, sections, colour scale |
| The program lists, the user judges | clean: ambiguous names, coordinate mismatches, unknown orientation, unclear array type, responses not removed and models unknown; stack: remove the candidates or not; reference: is the picked ridge the fundamental mode; inversion: starting model, deepest node, velocity bounds |
| The assistant may decide (and records in decisions.md) | processes and threads (from cores and the memory estimate); no clean archive when the disk is short |
| The assistant must ask first | all judgements above; `data.rename`; excluding a whole station; borrowing a reference curve or starting model; any change of a quality-control threshold or of the regularisation; deleting existing results |

## 4. The procedure

```bash
lant check                                                   # 0. installation; install what it lists
lant run all --raw DIR --stations FILE --out OUT             # 1. writes OUT/config.yaml and runs
lant status -c OUT/config.yaml --json                        # 2. after every stop: state of each step, reason, next command
```

`lant status --json` is the only interface you need to parse: per step the state (stands, stale, waiting for confirmation, not run, failed), the reason, and the `next` command. The loop:

1. Run the command given as `next`.
2. After the stop read `OUT/records/steps/<step>.md` and the files under "LOOK AT" (look at the figures if you can read images; otherwise ask the user to).
3. Check by the table below; report numbers, figure locations and the items to judge.
4. When the user has decided, carry it out (`lant confirm …` or a change of the configuration) and add a line to `decisions.md`.
5. Back to 1, until `figures` stands.

| Step | Look at | Common cases |
|---|---|---|
| clean | stations matched, stations without data, `coverage.png`, whether responses were removed, `pulse_1s.csv` | names do not match: ask for the naming rule, use `data.rename`; channel not recognised: `data.channel` |
| cc | pairs, share of dropped segments, memory estimate | not enough memory: `correlation.max_distance_km` or a smaller `nproc` |
| stack | `candidates.txt`, `station_qc/station_qc.png` | more candidates than expected by chance (1): the chain stops; show them, the user decides. One candidate: listed in the report only |
| reference | `fj.png`: is the ridge continuous, is it the slowest branch; measured velocities against `velocity_range_kms` | possibly a higher mode: give your reasons, propose `track_fundamental`, a narrower velocity range or a borrowed curve; the user decides |
| disp | accepted samples, rejections by reason, `3_dispersion_curves.png`, `median_vs_reference.txt` | few points: report which criterion rejects them; do not loosen it yourself |
| inversion | `4_start_model.png`, `4_kernels.png`, number of periods, depth nodes | more periods than the solver takes: change `per_octave` as the message says |
| run, evaluate | share kept by the precheck, residuals, checkerboard recovery, control differences | main project invalid in the precheck: stop and report |
| figures | slices, sections | show the checkerboards along with them |

## 5. Records

One line per decision in `OUT/records/decisions.md`:

```
| date | step | key or item | decision | basis | who |
| 2026-11-09 | stack | 3 candidate station sub-stacks | kept | correlation just below the limit, shifts inside the limit (03_stack/candidates.txt) | user |
```

Confirmations are written by `lant confirm` into the `confirmed:` section at the end of the configuration; the other files in `records/` are written by the program and are not edited by hand.

## 6. When something goes wrong

| What you see | What to do |
|---|---|
| Configuration error (also for removed old keys) | the message names the key and the reason; follow it |
| A step fails | read the last entry of `run_log.jsonl` and the traceback with `--debug`; a data problem: ask the user; a defect of the program: stop and report |
| A step is stale | the configuration or an earlier result changed; run again from the step `next` names; confirmations have to be given again |
| The output directory is locked | another run is active; do not delete the lock file, find that process first |
| The result looks unreasonable | do not "repair" it with parameters; report what you see, the numbers and possible causes, and ask the user |
