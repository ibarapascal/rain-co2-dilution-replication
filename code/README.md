# code/

| Path | What it is |
|---|---|
| `run.py` | single entry point: `--list`, `--stage ID [--mode M]`, `--all`, `--dry-run`, `--resume`, `--force-deps` |
| `stages.py` | stage table (39 stages, 33 in `--all`): commands, modes, dependencies, completion markers, resources, measured run times (`ref_run`), anemometer-height effect |
| `config.example.toml` | central configuration template (paths, options, upstream overrides, input-SHA gates, off by default) |
| `pipeline/` | the stage scripts (analysis) and the I/O layer: `repro_paths.py` (configuration, upstream locations, input-SHA gates, anemometer-height switch, `module_sha`), `repro_io.py` (downloads into the raw archive; the archive is never overwritten), `raw_sync.py` (copy caches into the archive) |
| `pipeline/raw_archive_lists/` | granule lists (URL, bytes, checksums) for the SPURS-1/2 and PISTON downloads |
| `expected/products.tsv` | reference-run SHA-256 of every product and the comparison rule per file (used by `tools/diff_expected.py`) |
| `tools/` | `diff_expected.py` (compare a re-run with the released data), `fill_upstream_sha.py` (optional: expected values for the input-SHA gates when `check_upstream_sha = true`), `build_event_table.py` (646-row merged event table), `fetch_witte_capsule.py` (third-party capsule files) |
| `tests/test_public.py` | static self-tests (standard library, seconds) |

The stage scripts keep their Chinese comments and some Chinese run-time messages; the implementation choices listed in each script
header use letter codes (I, J, K, L, M, ...) that the comments refer to. See `../docs/reproduction.md` for how to run the pipeline
and compare the results.
