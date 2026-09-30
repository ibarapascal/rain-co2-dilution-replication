# Changelog

Changes to the released data and code, newest first. Dates are JST.

## [0.9.1] — metadata fix

- `CITATION.cff`: `license` is now a single identifier (`CC-BY-4.0`), so that Zenodo can read the file; the code in `code/` stays
  under the MIT License (`LICENSE`), as stated in the README. Version 0.9.0 was not archived on Zenodo for this reason.
- `docs/third-party-data.md`: acknowledgment wording for GTMBA, IMOS and ERA5 brought in line with the current text requested by
  each provider.
- Data and code are unchanged.

## [0.9.0] — first public version

### Contents
- `data/events/events_merged.csv`: 646 rain events at 13 MAPCO2 moorings (2004-2025), one row per event with the quantities of all
  mooring stages (event definition, observed and RIM-3/S20 freshening, pCO₂ change, mechanism diagnostics, chemical slope, wind
  dependence, anemometer-height sensitivity, freshwater budget).
- `data/products/`: 73 files from 28 pipeline stages: per-event and per-overpass tables, summaries, the CMORPH pixel cache of the
  mooring test, and the global recalculation for 2000 (U(s) curves, band, basin and month tables, fluxes by pixel wind speed).
  `MANIFEST.tsv` gives the SHA-256 of each file as written by the pipeline and as released, and the release treatment (Chinese text
  translated, local paths removed, descriptive text as written by the current scripts; numbers unchanged).
- `data/raw_index/raw_manifest.jsonl.gz`: index of the 65,882 third-party raw files (URL, bytes, SHA-256).
- `code/`: pipeline of 39 stages (33 in `run.py --all`) from raw downloads to all products, with one entry point (`run.py`), a stage
  table with measured run times and memory, a central configuration, an append-only raw archive, optional input-SHA gates against
  upstream drift (off by default), an anemometer-height switch, a comparison tool for re-runs (`tools/diff_expected.py`) and static self-tests.
- `docs/`: reproduction guide, third-party data and required acknowledgments, translation table.

### Verified
- All 33 `run.py --all` stages were run end to end from the original raw downloads on 2026-09-30 (macOS, NumPy 2.5.3). Every result
  value equals the released data: 36 files are byte-identical, 32 JSON products agree leaf by leaf, the two CMORPH pixel caches agree
  as row multisets; the only differences are run records (run times, code hashes, input paths, download counts).
  `tools/diff_expected.py` exits 0 on these outputs and on the released products themselves.
- After that re-run the code changed only in comments, docstrings, console messages, environment-variable names, User-Agent
  strings, one netCDF title attribute and the check that `dw13a` and `dw33` import the pipeline's own copy of the module they wrap
  (`repro_paths.module_sha`). Later the scripts stopped writing a few descriptive fields, reworded others and renamed five keys
  (the released copies match the current scripts), the input-SHA gates became optional (off by default) and `dw33-post` no longer
  writes `dw33_keys.md` (`docs/reproduction.md`, "Changes since the full re-run"). A syntax-tree comparison confirms that the
  computations are unchanged.
- `code/tests/test_public.py`: 23 tests pass.
