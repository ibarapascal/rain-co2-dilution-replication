# Rain freshening at ocean moorings and the RIM-3 rain dilution scheme: derived data and code

This repository contains derived data on rain-induced freshening of the upper ocean and the scripts that compute them from public
raw data. Version 0.9.1.

## What the data are

- **Rain events at moorings.** 646 rain events (≥ 10 mm in 24 h after at least 24 h without rain) at 13 moorings of the NOAA PMEL
  carbon mooring network (MAPCO2) with co-located rain gauges, 2004-2025: 8 tropical TAO/RAMA moorings (group A) and 5 subtropical
  and midlatitude moorings (group B). For each event: rain, pre-rain salinity, observed salinity change at the MAPCO2 channel
  (nominal 0.5 m), at ~1 m and at 5 m, pCO₂ change, and the same quantities for up to five rain-free control windows.
- **Comparison with RIM-3.** The rain dilution scheme RIM-3 and the surface relation S20 of Witte et al. (2026), driven by CMORPH
  rain at the mooring pixel, evaluated at the observed depths; ratios of observed to modeled freshening with cluster-bootstrap
  intervals; sensitivity variants (thresholds, windows, gauge forcing, anemometer height).
- **Mechanism diagnostics.** Stratification, mixing (law-of-the-wall reference), rain mismatch between gauge and satellite, deep
  freshwater budgets to 10-25 m, freshwater conservation of RIM-3, wind dependence, chemical slope of pCO₂ versus salinity.
- **Surface and cross-platform tests.** Satellite surface salinity (SMAP RSS L2C V6 and JPL L2B V5) collocated with post-rain
  overpasses at the moorings; near-surface profiles from the SPURS-2 and PISTON field campaigns and the SPURS-1 mooring.
- **Global recalculation.** A reimplementation of the global calculation of the rain effect on the air-sea CO₂ flux for the year 2000
  (CMORPH, ERA5, HYCOM, OISST, Watson et al. fCO₂, GLODAP): reproduction of the original global numbers, flux as a function of a
  uniform scaling s of the RIM-3 freshening, band/basin/month tables, and fluxes accumulated by pixel wind speed.

`data/README.md` documents every file and column.

## Repository layout

| Path | Content |
|---|---|
| `data/events/events_merged.csv` | one row per rain event (646 rows); per-event quantities of all mooring stages |
| `data/products/<stage>/` | outputs of each pipeline stage (per-event and per-overpass tables, band/basin/month tables, summaries); `MANIFEST.tsv` lists every file with its SHA-256 and release treatment |
| `data/raw_index/raw_manifest.jsonl.gz` | index of the 65,882 third-party raw files the pipeline downloads (URL, bytes, SHA-256) |
| `data/SHA256SUMS` | checksums of all files in `data/` |
| `code/run.py`, `code/stages.py` | entry point and stage table (39 stages; 33 run by `--all`) |
| `code/pipeline/` | stage scripts and the I/O layer (central paths, raw archive, input-SHA gates) |
| `code/tools/` | comparison with the released data, input-SHA gate helper, merged event table, third-party capsule download |
| `code/expected/products.tsv` | reference SHA-256 and comparison rules used by `code/tools/diff_expected.py` |
| `code/tests/test_public.py` | static self-tests |
| `docs/reproduction.md` | how to re-run and compare; stage run times and memory; known issues |
| `docs/third-party-data.md` | third-party sources, identifiers, terms and required acknowledgments |
| `docs/translation-table.tsv` | the Chinese strings and keys written by the scripts and their English form in the released data |
| `environment/requirements.txt` | pinned Python packages |

## Reproducing the data

Environment: Python 3.12 and `environment/requirements.txt`; a NASA Earthdata account (token file) for the satellite and
field-campaign downloads; ~160 GB of disk for the raw archive and ~100 GB for the global stage; 12 GiB of memory for the global
stages (run one after the other).

```
python3.12 -m venv .venv && .venv/bin/python -m pip install --only-binary=:all: -r environment/requirements.txt
source .venv/bin/activate
cd code
python3 -m unittest tests/test_public.py              # static self-tests, seconds, no data
cp config.example.toml config.local.toml              # set raw_root, fast_root, out_root, token_file
python3 run.py --config config.local.toml --list      # stages, dependencies, network/token needs
python3 run.py --config config.local.toml --all       # all stages in dependency order
python3 tools/diff_expected.py <out_root>             # compare the re-run with data/products/
```

Run time: about 12 h of wall clock with the raw downloads in place (measured; `p3-global` 7.0 h with 2 workers of up to 4.6 GB,
`dw33` 4.1 h with 2 workers of up to 5.0 GB, all other stages ≤ 3 GiB and seconds to minutes each), 1-1.5 days including the
downloads (estimated, bandwidth-bound). Per-stage times and memory are in `docs/reproduction.md` and in the `ref_run` field of
`code/stages.py`.

Five stages can check the SHA-256 of their upstream inputs against the upstream products of the reference run. This check is off
by default (`[options] check_upstream_sha = false`): the stages record the SHA-256 of what they read and continue. With
`check_upstream_sha = true` they stop (exit 3) on a mismatch until the new values are entered in the configuration
(`tools/fill_upstream_sha.py`, described in `docs/reproduction.md`).

`tools/diff_expected.py` compares every product with the released copy: files byte for byte where the pipeline output is
deterministic, JSON summaries leaf by leaf (numbers within 1e-12 + 1e-9·|x|, skipping run records such as times, code hashes and
input paths), and the pixel cache as a multiset of rows. All stages were re-run from the original raw downloads on 2026-09-30 and
matched the released data value by value.

## Third-party data and licenses

All raw data are public: MAPCO2 moorings (NOAA PMEL), GTMBA and OceanSITES moorings (NOAA PMEL, NDBC, WHOI, University of Hawaii,
IMOS), CMORPH CDR, ERA5, HYCOM, OISST, the Watson et al. fCO₂ product, GLODAPv2.2016b, WOA09, SMAP RSS and JPL salinity, SPURS-1/2 and
PISTON field data. They are downloaded by the pipeline and not redistributed here, except a small CMORPH pixel extract, GLODAP values at
the moorings and PISTON ship positions. Identifiers, terms and required acknowledgments are in
[`docs/third-party-data.md`](docs/third-party-data.md). The original code of Witte et al. (2026) is available from Code Ocean
(DOI 10.24433/CO.9378898.v1); `code/tools/fetch_witte_capsule.py` downloads the two files the pipeline uses.

Code (`code/`): MIT License (`LICENSE`). Data and documentation (`data/`, `docs/`): CC BY 4.0 (`LICENSE-DATA`). The terms of the
third-party sources apply in addition.

## How to cite

Archived on Zenodo: https://doi.org/10.5281/zenodo.23067878 (all versions; each release also has its own DOI, listed on the
Zenodo record). Citation metadata are in `CITATION.cff`.
