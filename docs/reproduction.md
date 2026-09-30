# Reproduction guide

How to recompute the derived data in `data/products/` from public raw data with `code/run.py`, and how to compare a re-run with
the released data.

All 33 `run.py --all` stages were re-run end to end from the original raw downloads on 2026-09-30 (macOS workstation, NumPy
2.5.3). Every result value equalled the released data; the only differences were run records (times, code hashes, input paths,
download counts). The code has changed since then only in comments, docstrings, console messages and I/O metadata (see "Changes
since the full re-run" below).

## Ways to use the repository

1. **Use the data** (no installation): `data/events/events_merged.csv` (one row per rain event) and the per-stage products in
   `data/products/`; column dictionaries in `data/README.md`.
2. **Re-run the mooring stages from the released pixel cache**: P1/P1b download the mooring records (~0.55 GB); P2 can then reuse
   `data/products/p2/cmorph_pixels.csv` and `p2_glodap.json` as its cache (`p2_rim_test.py --resume-from data/products/p2`, or
   point `[upstream] p2_dir` there for the downstream stages), avoiding ~41 GB of CMORPH downloads. This shortcut has not been
   tested from the released copies.
3. **Full re-run from raw data** (this guide): about 12 h of wall clock with the raw archive in place (measured), 1-1.5 days when
   everything has to be downloaded (estimated, bandwidth-bound); ~100 GB of disk for the global stage and ~160 GB for all raw
   archives.

## Requirements

- Python 3.12 and the packages in `environment/requirements.txt` (P1, P1b and the download stages need only the standard library).
- A NASA Earthdata account and a token file for the SMAP, SPURS and PISTON downloads (stages `p5-sss`, `p7a`, `raw-spurs2`,
  `raw-p7e`); all other sources are anonymous.
- The two files of the Witte et al. (2026a) Code Ocean capsule for `p4c` (and for the self-test of `p3-global`, which fetches them
  itself): `python3 code/tools/fetch_witte_capsule.py`.
- Memory: peak memory of most stages is ≤ 3 GiB. `p3-global` (11 GiB requested; 2 workers of up to 4.6 GB each) and `dw33`
  (10 GiB requested; up to 5.0 GB per worker) must run one after the other; with 12 GiB available, only stages of ≤ 1 GiB fit
  beside `p3-global` and ≤ 2 GiB beside `dw33`.
- Disk: the P3 input cache is ~57 GB plus ~42.6 GB of ERA5 originals (the downloader checks for ≥ 80 GB free); the raw archive of all
  stages grows to ~160 GB; stage outputs are < 5 GB.
- Network: downloads are rate-limited per stage (`[options] raw_rate_mbps`, default 2 MB/s); the reference runs used ≤ 4 MB/s in
  total.

## Setup

```
python3.12 -m venv .venv && .venv/bin/python -m pip install --only-binary=:all: -r environment/requirements.txt
source .venv/bin/activate
cd code
python3 -m unittest tests/test_public.py       # static self-tests (20, seconds, no data)
cp config.example.toml config.local.toml       # edit [paths]; token_file = path only; out_root = new directory per run
python3 run.py --config config.local.toml --list
python3 run.py --config config.local.toml --all --dry-run
```

## Running

```
python3 run.py --config config.local.toml --stage p1              # one stage (default mode in stages.py)
python3 run.py --config config.local.toml --all                   # all 33 in_all stages in dependency order; completed stages skipped
python3 tools/fill_upstream_sha.py <out_root>                     # optional, only with check_upstream_sha = true (below)
python3 tools/diff_expected.py <out_root> --report diff.tsv       # compare the re-run with the released data
```

`run.py` never deletes files and refuses a non-empty output directory unless `--resume` is given. Each stage writes
`<out>/_repro_run.json` (command, config SHA-256, script SHA-256, staged inputs, options, times, exit code) and `_repro_stdout.txt`.
The child process gets `REPRO_CONFIG` (absolute path of the configuration) and `REPRO_OUTPUT_DIR` (stage output directory); an
optional `REPRO_ATTEMPT_ID` is written into the provenance tag of archived raw files, and `REPRO_RECORD_HOST=1` records the host
name in the run record. Scripts that need an Earthdata token read only the configured file path and put its content only into
requests to Earthdata hosts.

## What the I/O layer does

The stage scripts contain the analysis; `code/pipeline/repro_paths.py`, `repro_io.py` and `raw_sync.py` handle input and output:

- **Paths**: every data location comes from the configuration (`[paths]`, `[upstream]`); without a configuration, paths resolve to
  `/__REPRO_UNCONFIGURED__/<key>`, so a misconfigured run fails immediately instead of writing to a default disk.
- **Raw downloads are kept**: CMORPH hourly files and the GLODAP archive are written to the raw archive (`<raw_root>`) and read from
  there instead of being deleted after extraction; ERA5 u/v originals are kept after deriving wind speed; every archived file gets a
  line in `<raw_root>/manifest.jsonl` (URL, bytes, SHA-256, time). `raw-sync` copies the originals that the P1/P1b/P3 downloaders keep
  in their own caches into the same archive.
- **The raw archive is append-only**: `raw_archive.py` and `raw-sync` stop instead of replacing an archived file of different size;
  `p5_sss_sat.py` saves a second OPeNDAP subset of the same name as `<name>.dup<time>`.
- **Input-SHA gates** (below) and the **anemometer-height switch** (below) are read from the configuration.
- **Module checks**: `dw13a` wraps `p2_rim_test.py` and `dw33` reuses the engine of `p3_global.py`; both check at start that the
  imported module is the file in `code/pipeline/` (`repro_paths.module_sha`). `p3-global` records the SHA-256 of `p3_global.py` in
  its static fields, and `dw33` refuses static fields written by a different version of that file.

Analysis logic, thresholds, random seeds and decision rules are not affected by any of these.

## Dependency graph

Arrow = "runs before"; `raw-sync` may run any time after the P1/P1b/P3 downloads.

```
p1 -> p1b -> p2 -> p2-sens
             |-> p4-mech -> p4b -> p4e
             |      '-------'--> p4f -> p7c (p7c also reads p4-mech, p4b)
             |      '-> cal-dw13-meta (also reads the p1b station list)
             |-> p3c <-- p3-global <-- p3-fetch-extra <-- p3-fetch
             |    |-> p4-explore (also reads p1, p1b, p2, p3-global)
             |    '-> p5-beta (also reads p2, p3-global, p4-mech)
p1, p1b -> p5-sss -> p7a, p7d          p3-global -> p4d        p4c (no upstream)
p5-sss, p7a, p7c -> p8a
raw-spurs2 -> p6, p8b                  raw-p7e -> p7e, piston-track
p5-sss, p7a, p6, p8a, p7e, p8b -> cal-surface-diag
p3-global -> dw33 -> dw33-post (also reads p4-mech, p3c)
p2, cal-dw13-meta -> dw13a
```

## Stages

Stage IDs are the names of the output directories (`data/products/<stage>/`). Six stages are not part of `--all`: `env-check`
(environment self-test), `p6-plan`, `p7e-plan`, `p8b-plan` (metadata and availability surveys that do not read salinity values),
`dw33-bench` (one-day benchmark of the wind-binned engine) and `raw-backfill` (backfill of originals deleted by earlier script
versions).

| Stage | Content |
|---|---|
| p1, p1b | rain events at 8 tropical TAO/RAMA moorings (p1) and 5 subtropical/midlatitude moorings plus new periods (p1b): 646 events, event and control-window salinity changes, gates D1/D2/D4/D5 |
| p2, p2-sens | mooring test: CMORPH-driven RIM-3 and S20 predictions at the observed depths, ratio R_RIM(1 m), discrimination gate D3; four sensitivity variants |
| p3-fetch, p3-fetch-extra, p3-global, p3c | global recalculation for the year 2000: inputs, reproduction of the Witte et al. global numbers (P3a), substitution curves U(s) (P3b), decision table (P3c) |
| p4-explore | exploratory analyses (few-cluster inference, 0.5 m channel, latitude concentration, descriptives) |
| p4-mech, p4b, p4e | mechanism diagnostics at the moorings (mixing, stability, rain mismatch, deep freshwater budget, budget closure) |
| p4c, p4d, p4f | RIM-3 freshwater conservation: synthetic check, global uniform scaling, wind trend |
| p5-beta | chemical slope r_beta = beta_obs / beta_model (seven candidates) |
| p5-sss, p7a, p7d, p8a | satellite surface salinity collocation (SMAP RSS L2C V6, JPL L2B V5): ratios, time since rain, wind strata |
| p6, p7e, p8b | in situ near-surface profiles: SPURS-2 (salinity snake, Lady Amber layer ratio), PISTON, SPURS-1 |
| p7c | wind dependence of the mooring ratio (H1-H6) |
| raw-spurs2, raw-p7e, raw-sync | raw downloads and archive |
| piston-track | PISTON ship positions (descriptive; not used in any estimate) |
| cal-dw13-meta, dw13a | sensor depth and anemometer height from public deployment metadata; anemometer-height sensitivity of the mooring test |
| cal-surface-diag | C1: tropical satellite / in situ surface ratio; C2: denominator shares of each test |
| dw33, dw33-post | global RIM-3 fluxes accumulated by pixel wind speed; substitution of the mooring ratios by wind class |

## Resources and run times

Reference run = the run that produced the released data (with downloads). Full re-run = the end-to-end re-run of 2026-09-30 with
the raw archive already in place, so almost nothing was downloaded; it is the stage time of a complete run, not a benchmark. The
same numbers are in the `ref_run` field of `code/stages.py`.

| Stage | Main inputs | Resources (requested) | Reference run | Full re-run |
|---|---|---|---|---|
| p1 | PMEL ERDDAP MAPCO2, NDBC OceanSITES GTMBA (anonymous) | 1 GiB | 1,176 s; 95.4 MB downloaded | 22 s |
| p1b | as p1, new moorings | 1 GiB | 1,919 s; 455.9 MB | 95 s |
| p2 | P1/P1b caches, CMORPH 25,400 hours, GLODAP | 2 GiB | 20,868 s (CPU 1,951 s; 41.1 GB, bandwidth-bound) | 985 s |
| p2-sens | P2 pixel cache, 20,411 more CMORPH hours | 2 GiB | 15,083 s (33.2 GB) | 890 s |
| p3-fetch | CMORPH 2000, ERA5 u/v (NSF NCAR mirror), HYCOM, OISST, Watson (RECCAP2), GLODAP | 1 GiB, ≥ 80 GB free | 12,559 s; 57.2 GB | 601 s (conversion only) |
| p3-fetch-extra | ERA5 2001-01, WOA09 basin mask | 1 GiB | 781 s; 0.98 GB | 50 s |
| p3-global | P3 cache, capsule (self-test) | 2 workers, 11 GiB | 25,040 s (4.39 GB per worker) | 25,373 s |
| p3c, p4-explore, p4c, p4d, p4f, dw33-post, cal-surface-diag | small products | ≤ 1.5 GiB | seconds | 0.3-23 s |
| p4-mech, p4b, p4e, p5-beta, p7c | caches, P2 pixel cache | 1.25-1.5 GiB | 88-209 s each | 98-242 s |
| p5-sss | CMR, PO.DAAC OPeNDAP RSS SMAP L2C V6 (token), CMORPH | 1.5 GiB | 8,699 s (13.2 GB CMORPH) | 338 s |
| p7a | JPL SMAP L2B CAP V5 (token), CMORPH | 1.5 GiB | 1,944 s | 122 s |
| raw-spurs2, raw-p7e | PO.DAAC SPURS-2, NASA ASDC PISTON, PO.DAAC SPURS-1 (token) | 0.5 GiB | 1,274 s (2.53 GB); PISTON/SPURS-1 similar | 21 s, 42 s (verification only) |
| p6, p7e, p8b, piston-track | raw archive (read-only) | ≤ 3 GiB | 4-35 s | 0.5-22 s |
| p7d, p8a | P5/P7a/P7c products | ≤ 2 GiB | 1-197 s | 1 s, 90 s |
| cal-dw13-meta | NCEI OCADS and OceanSITES metadata (anonymous) | 0.5 GiB | minutes | 306 s |
| dw13a | P2 directory, D13M | 2 GiB | 214 s | 215 s |
| dw33 | P3 cache and static fields | 2 workers, 10 GiB | 4.05 h | 14,670 s |

Totals of the full re-run: 45,077 stage-seconds (12.5 h), of which `p3-global` 7.0 h and `dw33` 4.1 h; 11.7 h of wall clock in
batches: P3 download (0.2 h) in parallel with the mooring chain P1 … P8b (0.5 h) → `p3-global` (7.0 h; `raw-sync`, 0.5 GiB, beside
it) → `dw33` (4.1 h) in parallel with the satellite and post-processing stages (0.5 h) → the last gated stages (minutes), plus the
manual input-SHA gate steps if `check_upstream_sha = true`. `p3-global` and `dw33` are the critical path and cannot overlap (memory). Without the raw archive, add
the download times of the reference run (P3 download 3.5 h, P2 5.8 h, P2-sens 4.2 h, P5 2.4 h, P1/P1b 0.9 h; estimated 1-1.5 days in
total, bandwidth-bound).

## Decisions before a full re-run

1. **P1/P1b data**: a fresh download may return revised real-time GTMBA records, which can change the event set; P1b, P2 and P8a
   stop (exit 3) if the fixed counts (204, 646, D4 = 34) change. Alternative: seed `<out>/p1/cache` and `<out>/p1b/cache` from an
   archive of the original downloads and run with `--resume` (the scripts skip cached files).
2. **P2 pixel cache**: re-extract from CMORPH originals (closest to "from scratch") or point `[upstream]` at the released
   `data/products/p2/cmorph_pixels.csv` (not from scratch, much faster).
3. **P3 cache**: `p3-fetch` re-downloads ~57 GB; ERA5 originals are kept (another ~42.6 GB).
4. **P5/P7a caches**: new directories re-query CMR/OPeNDAP; PO.DAAC reprocessing could change results.
5. **Anemometer height** (`[options] wind_height`): `uniform_4m` reproduces the released data; `per_station`/`per_source` give
   corrected variants (use a separate out_root; see below).
6. **Deployment metadata** (`cal-dw13-meta`): re-fetching reads the metadata as of the run date (NCEI may have added files).

## Input-SHA gates

Five stages can check the SHA-256 of their upstream inputs before they compute anything: `p7d` reads P5; `p8a` reads
P5/P7a/P7c; `cal-surface-diag` reads P5/P6/P8a/P7e/P8b; `dw33-post` reads P4/P3B/P3C; `dw13a` reads D13M. The expected values
are the SHA-256 values of the upstream products of the reference run (`[upstream_sha]` overrides them).

**Default: off.** With `[options] check_upstream_sha = false` (the default, also when the key is absent) the stages do not compare:
they record the SHA-256 of the files they read in their outputs and continue, so `run.py --all` runs through. Compare the results
with the released data afterwards (`tools/diff_expected.py`).

**Optional: on.** With `check_upstream_sha = true` the stages only run on the exact upstream products they were written for. A
re-run of an upstream stage writes new files (at least the run time inside them changes), so a full re-run then stops at each of
these stages (exit 3). To continue:

1. check the new upstream products (for example with `tools/diff_expected.py <out_root> --present-only`);
2. run `python3 tools/fill_upstream_sha.py <out_root>` (read-only; `--only p8a,cal12` restricts the output); it prints, for every
   gate, either a commented line (unchanged) or a TOML line with the new SHA-256 and the upstream stage it depends on;
3. paste the printed lines into the `[upstream_sha]` table of your configuration and continue with `run.py --resume`.

Do not fill the table automatically: in this mode the stop is the point at which a person confirms that the new upstream products
are the intended ones.

## Anemometer height

The pipeline converts mooring wind to 10 m assuming a 4 m anemometer everywhere. Recorded heights at group B moorings are
3.26-3.84 m (WHOTS), 3.29-3.55 m (SOFS), 3.29-3.35 m (Stratus), 3.8-4.2 m (KEO, Papa). `wind_height = "uniform_4m"` (default) returns
the original conversion factor unchanged (bit-identical results). `"per_station"` uses per-station medians (WHOTS 3.355, SOFS 3.3125,
Stratus 3.32 m; `[wind_height_m]` overrides them). `"per_source"` uses the height of the file/variable that supplied each hourly value
and is implemented only by stage `dw13a`; `run.py` refuses other wind-affected stages under `per_source` and refuses `dw13a` under
`per_station`. `run.py --list` marks the stages affected by the switch. With measured heights R_RIM(1 m) changes from 0.29616 to
0.29742 and no decision changes (`data/products/dw13a/`).

## Comparing a re-run with the released data

`tools/diff_expected.py OUT_ROOT [--present-only] [--report FILE]` compares every product listed in `code/expected/products.tsv`
with the reference-run SHA-256 and with the released copy in `data/products/`:

- a file with the reference-run SHA-256 is `same`; `expect_identical = yes` marks the 31 files that were byte-identical in the full
  re-run and are expected to stay so (five more JSON files were byte-identical then, but the current scripts write fewer or
  reworded text fields in them; see "Changes since the full re-run");
- JSON products (`compare = json`) with a different SHA-256 are compared leaf by leaf with the released copy (numbers within
  1e-12 + 1e-9·|old|; the Chinese keys and strings written by the scripts are matched to the released English ones through
  `docs/translation-table.tsv`). The leaf paths in the `ignore` column are skipped: run records (run times, code hashes, input paths
  and input hashes, download counts) and text fields whose key or wording in the reference-run outputs differs from what the
  current scripts write (some ignore entries name keys that exist only in those outputs). No numeric field is skipped. Any other
  difference is a failure (`differs`); a match is `same(json-leaves)`;
- `cmorph_pixels.csv` (P2, P2-sens) is compared as a multiset of rows (`compare = rows`), because its row order follows the order in
  which download threads finish;
- files that are not released (logs, plan outputs, self-test records) are `differs(expected)` when their SHA-256 differs;
- rows of stages outside `--all` (`dw33-bench`) are `not-run` when the stage was not run; `--present-only` skips stages that have not
  produced a directory yet (for comparisons during a re-run).

Exit code 1 on a missing file or an unexpected difference. On the released products themselves (`python3 tools/diff_expected.py
../data/products`) it reports no failures (`tests/test_public.py` checks this). On the complete outputs of the full re-run of
2026-09-30 it exits 0 (36 same, 32 same(json-leaves), 2 same(rows-multiset), 5 differs(expected) for unreleased logs and plan
files, 1 not-run). Counts, event numbers and decision categories must be identical; any change of a gate or interpretation class
should be reported.

## Changes since the full re-run

After the re-run of 2026-09-30, the code was changed only in comments, docstrings and console messages, in the names of the
environment variables (`REPRO_*`), in the User-Agent strings of the downloaders, in the title attribute of the derived ERA5 wind
files, and in how `dw13a` and `dw33` verify the imported module (`repro_paths.module_sha`, which reads the current file). Because
the scripts' own SHA-256 values changed, the `code_sha256` fields in new outputs differ from the released ones (they are run records
and are skipped by the comparison), and a `p3-global` static-field directory written by an earlier version of `p3_global.py` is
rebuilt rather than reused.

Later changes (2026-10-01), also without effect on any computed value (checked by a syntax-tree comparison of the scripts):

- a few descriptive fields are no longer written (references to planning documents in most summaries, a follow-up note in
  `p4-mech`, a wording note in `p5-beta`, a label in `p3c`, a condition note in the P2 exit record), some descriptive strings are worded
  without references to documents outside this repository, and five keys were renamed (`params_fixed`,
  `first_round_verdict_unchanged`, `GR_band_abs_s_minus_1_P30`, `s6_8p2`, `smooth_judgment`). The released copies contain exactly
  what the current scripts write (translated); `tools/diff_expected.py` still matches the reference-run outputs through the
  historical keys in `docs/translation-table.tsv` and the `ignore` column of `code/expected/products.tsv`;
- the input-SHA gates are off by default (`[options] check_upstream_sha`, see "Input-SHA gates");
- `dw33-post` no longer writes the Markdown key table `dw33_keys.md` (its values are all in `dw33_post.json`).

## Known issues

- **Hard-coded upstream constants** (not updated by a re-run): `p4d` (M_RIM, M1), `p4f` (wind cuts), `dw33` (class ratios and cuts),
  `dw33-post` (coverage), `p8a` (cuts). If P2/P4/P4b values change, these stages keep the reference-run constants (dw33-post
  detects its mismatch and stops). Check them after re-running P2/P4.
- **Fixed counts as gates**: P1b (204), P2 (646, D4 = 34), P7c, P7d, P8a, cal-surface-diag and dw33-post stop if upstream values
  change (intended).
- **Script revisions of the reference run**: P1 was run with an earlier revision of `p1_events.py` that lacked one optional argument
  used from P1b on; the P3 download ran an earlier revision of `p3_fetch.py` whose ERA5 conversion needed a manual directory creation
  (fixed in the current version); P3c ran a copy of `p3c_decide.py` that took its input directory from the command line. None of
  these affects any value.
- **Numerical environment**: some small stages of the reference run used NumPy 2.4.2 and three stages ran on Windows; last-digit
  differences are possible. In the full re-run (macOS, NumPy 2.5.3) all result values were identical, including those of the Windows
  stages.
- **Text written by the scripts**: code comments, run-time messages and many descriptive strings in the products are in Chinese; the
  released copies are translated (`docs/translation-table.tsv`).
- **raw-sync**: the step that copies the P5 cache into the archive (`--mode after-p5`) is not part of `--all`.
- **Raw archive is append-only**: `raw-sync` (and `raw_archive.py`) stop with an error instead of replacing an archived file whose
  size differs from the new copy. Files such as the P3 download manifests and logs change with every download run, so a second
  `raw-sync` into the same `raw_root` after a new `p3-fetch` stops (exit 3); point `raw-sync` at a new `raw_root` or move the old files
  aside after checking them.
- **Derived ERA5 wind files carry the code version**: `p3-fetch` writes its version string into the global attribute `created_by`
  of each `era5_ws10/*.nc`, so their SHA-256 changes with the code version although the data are unchanged. Compare these files by
  content, not by SHA-256.
- **NumPy deprecation**: `p3_fetch_extra.py` sets the shape of an array by assignment, which NumPy 2.5 reports as a
  DeprecationWarning; it works with the pinned versions but will fail once NumPy removes it.
- **Peak memory per batch** was recorded only inside `p3-global` and `dw33` (per month); the other stages' peaks come from the
  reference-run records.
