"""stages.py — stage table of the pipeline (single source of truth for commands, dependencies and resources);
read by run.py and tests/test_public.py.

Template placeholders (expanded by run.py):
  {out}           output directory of this stage, <out_root>/<stage id>
  {py}            interpreter ([options] python; default = the interpreter running run.py)
  {tag}           provenance tag written to <raw_root>/manifest.jsonl (repro-<stage>:<REPRO_ATTEMPT_ID or local>)
  {path:K}        repro_paths.path(K)       {up:K}  repro_paths.upstream(K)
  {raw:N}         <raw_root>/N              {opt:K} [options] K (e.g. raw_rate_mbps)
Fields
  script   script in pipeline/; modes = argv per mode (without interpreter and script); default = default mode
  deps     stages that must be complete first (run.py --all runs in topological order; a single stage only checks done files)
  inputs   small upstream files that run.py copies into the stage directory ("out" or "_inputs") for scripts with fixed layouts
  done     completion marker files
  net/token  needs network / needs a NASA Earthdata token (token file location from the local config; content only in request headers)
  machine  kind of machine the reference run used (descriptive); res = (RAM GiB, CPU millicores, estimated s, limit s), a
           resource request with headroom over the measured peak
  ref_run  wall clock and peak memory measured in the reference runs (descriptive)
  in_all   part of `run.py --all` (metadata surveys, backfills and environment tests are not)
  wind     effect of the anemometer-height switch ([options] wind_height): direct = the stage converts mooring wind itself;
           inherit = reads affected upstream products; const = reads an affected quantity hard-coded as a constant (does not follow
           the switch; docs/reproduction.md, known issues); source = the stage implements the per-source heights itself (dw13a)
  modes    usually argv lists; multi-step stages use {"steps": [argv1, argv2, ...]} (run in order, stop at the first non-zero exit)
Change Log
  2026-09-27 first version (integration of all stage scripts).
  2026-09-29 p8a, p8b, p8b-plan, piston-track, cal-dw13-meta, cal-surface-diag, dw33-bench, dw33, dw33-post; wind field.
  2026-09-29c dw13a (multi-step mode).
  2026-10-01 public dataset release: neutral titles; run notes as ref_run.
"""

P1ARGS = ["--p1-events", "{up:p1_events}", "--p1b-dir", "{up:p1b_dir}"]
P8A_ARGS = ["--p5-dir", "{up:p5_dir}", "--p7a-dir", "{up:p7a_dir}", "--p7c-events", "{up:p7c_events}"] + P1ARGS + ["--out", "{out}"]
DW33_ARGS = ["--cache-root", "{path:p3_cache}", "--static-dir", "{path:p3_work}/static", "--work-dir", "{path:dw33_work}",
             "--out", "{out}"]
LA_ROOT = "{raw:spurs2}/SPURS2_LADYAMBER"
DW13A_META = ["--meta", "{up:dw13_meta}"]
DW13A_STEPS = {
    "selftest": ["--selftest"],
    "nominal": ["--stage", "nominal", "--resume-from", "{up:p2_dir}"] + DW13A_META + ["--out", "{out}/nominal"],
    "gate": ["--gate", "--orig", "{up:p2_dir}", "--nominal", "{out}/nominal", "--out", "{out}"],
    "measured": ["--stage", "measured", "--resume-from", "{up:p2_dir}"] + DW13A_META + ["--out", "{out}/measured"],
    "compare": ["--compare", "--orig", "{up:p2_dir}", "--nominal", "{out}/nominal", "--measured", "{out}/measured", "--out", "{out}"],
}

STAGES = [
    {"id": "env-check", "script": "venv_selftest.py", "title": "Environment self-test (interpreter and package versions, write permission on the data disk)",
     "modes": {"run": ["--rw-dir", "{path:fast_root}"]}, "default": "run", "deps": [], "done": ["selftest.json"],
     "net": False, "token": False, "machine": "workstation", "res": (0.5, 500, 60, 600),
     "ref_run": "seconds", "in_all": False},

    {"id": "p1", "script": "p1_events.py", "title": "P1 tropical TAO/RAMA moorings (8): mooring download -> rain events -> gates D1/D2/D5",
     "modes": {"full": [], "smoke": ["--smoke"]}, "default": "full", "deps": [],
     "done": ["p1_summary.json", "p1_events.csv"], "net": True, "token": False, "machine": "workstation",
     "res": (1.0, 1000, 1800, 7200), "ref_run": "1,176 s with downloads (95 MB); 22 s with the download cache in place", "in_all": True},

    {"id": "p1b", "script": "p1b_extend.py", "title": "P1b one-time sample expansion (new periods + 5 group B moorings) -> combined set of 646 events, D4 = 34",
     "modes": {"full": ["--r-rule", "per-var", "--p1-events", "{up:p1_events}"],
               "plan": ["--plan", "--r-rule", "per-var", "--p1-events", "{up:p1_events}"], "selftest": ["--selftest"]},
     "default": "full", "deps": ["p1"], "done": ["p1b_summary.json", "p1b_events_new.csv", "p1b_stations.json"],
     "net": True, "token": False, "machine": "workstation", "res": (1.0, 1000, 7200, 21600),
     "ref_run": "1,919 s with downloads (456 MB), peak ~0.3 GB; 95 s with the cache in place", "in_all": True},

    {"id": "p2", "script": "p2_rim_test.py", "title": "P2 main set: CMORPH-driven RIM-3/S20, gate D3, R_RIM(1 m), exit decision",
     "modes": {"full": ["--part", "primary"] + P1ARGS, "plan": ["--plan"] + P1ARGS, "selftest": ["--selftest"],
               "smoke": ["--smoke"] + P1ARGS},
     "default": "full", "deps": ["p1", "p1b"], "done": ["p2_summary.json", "p2_d3.json", "p2_events.csv", "cmorph_pixels.csv"],
     "net": True, "token": False, "machine": "workstation", "res": (2.0, 2000, 14400, 57600),
     "ref_run": "20,868 s with downloads (CPU 1,951 s; 41.1 GB, bandwidth-bound); 985 s with the raw archive in place", "wind": "direct", "in_all": True},

    {"id": "p2-sens", "script": "p2_rim_test.py", "title": "P2 sensitivity part (K21, four items; --resume-from the main-set pixel cache)",
     "modes": {"full": ["--part", "sensitivity", "--resume-from", "{up:p2_dir}"] + P1ARGS}, "default": "full",
     "deps": ["p2"], "done": ["p2_sensitivity.json"], "net": True, "token": False, "machine": "workstation",
     "res": (2.0, 2000, 10800, 43200), "ref_run": "15,083 s with downloads (CPU 1,079 s; 33.2 GB), peak ~0.57 GiB; 890 s with the raw archive in place", "wind": "direct", "in_all": True},

    {"id": "p3-fetch", "script": "p3_fetch.py", "title": "P3 global inputs download (CMORPH 2000, ERA5 -> ws10, HYCOM, OISST, Watson, GLODAP)",
     "modes": {"full": ["--cache", "{path:p3_cache}"], "verify": ["--cache", "{path:p3_cache}", "--verify"],
               "selftest": ["--selftest"]},
     "default": "full", "deps": [], "ensure_dirs": ["{path:p3_cache}"], "done": ["fetch_summary.json"],
     "net": True, "token": False, "machine": "workstation", "res": (1.0, 1000, 43200, 108000),
     "ref_run": "12,559 s, 57.2 GB (~4.7 MB/s); ERA5 originals (~42.6 GB) are also kept in the cache; 601 s with the raw archive in place (conversion only)", "in_all": True},

    {"id": "p3-fetch-extra", "script": "p3_fetch_extra.py", "title": "P3 supplementary download (first 24 ERA5 steps of 2001-01, WOA09 basin mask)",
     "modes": {"full": ["--cache", "{path:p3_cache}"]}, "default": "full", "deps": ["p3-fetch"],
     "done": ["fetch_extra_summary.json"], "net": True, "token": False, "machine": "workstation", "res": (1.0, 1000, 3600, 14400),
     "ref_run": "781 s, peak 0.98 GB; 50 s with the raw archive in place", "in_all": True},

    {"id": "p3-global", "script": "p3_global.py", "title": "P3 global recalculation: P3a reproduction gate + P3b substitution curves (monthly, parallel)",
     "modes": {"full": ["--mode", "full", "--workers", "2", "--cache-root", "{path:p3_cache}", "--work-dir", "{path:p3_work}"],
               "selftest": ["--mode", "selftest", "--cache-root", "{path:p3_cache}"],
               "bench": ["--mode", "bench", "--cache-root", "{path:p3_cache}"]},
     "default": "full", "deps": ["p3-fetch", "p3-fetch-extra"], "ensure_dirs": ["{path:p3_work}"],
     "done": ["selftest.json", "p3a_repro.json", "p3b_curves.json", "p3_run.json", "p3_bands.csv", "p3_regions.csv", "p3_monthly.csv"],
     "net": False, "token": False, "machine": "workstation", "res": (11.0, 2000, 27000, 72000),
     "ref_run": "25,040-25,373 s with 2 workers, peak 4.4-4.6 GB per worker (11 GiB requested)", "in_all": True},

    {"id": "p3c", "script": "p3c_decide.py", "title": "P3c decision table (M1, S1-S7, E1, gates G-X/G-R/G-N/G-0)",
     "modes": {"run": ["{out}"]}, "default": "run", "deps": ["p2", "p3-global"], "stage_dir": "out",
     "inputs": {"p3a_repro.json": "{up:p3_dir}/p3a_repro.json", "p3b_curves.json": "{up:p3b_curves}",
                "p2/p2_summary.json": "{up:p2_dir}/p2_summary.json", "p2/p2_events.csv": "{up:p2_dir}/p2_events.csv"},
     "done": ["p3c_decision.json"], "net": False, "token": False, "machine": "laptop",
     "res": (0.5, 1000, 60, 600), "ref_run": "under 1 s",
     "wind": "inherit", "in_all": True},

    {"id": "raw-sync", "script": "raw_sync.py", "title": "Copy raw originals into the raw archive (P1/P1b caches, P3 cache incl. ERA5 originals)",
     "modes": {"run": ["--stages", "A1,A2,A3,A4,A6", "--raw-root", "{path:raw_root}", "--out", "{out}", "--tag", "{tag}"],
               "after-p5": ["--stages", "A5", "--raw-root", "{path:raw_root}", "--out", "{out}", "--tag", "{tag}"]},
     "default": "run", "deps": ["p1", "p1b", "p3-fetch", "p3-fetch-extra"], "done": ["raw_sync_status.json"],
     "net": False, "token": False, "machine": "workstation", "res": (0.5, 1000, 7200, 28800),
     "ref_run": "copy only; minutes", "in_all": True},

    {"id": "p4-explore", "script": "p4_explore.py", "title": "Exploratory analyses PX (few-cluster inference, 0.5 m channel, latitude concentration, descriptives)",
     "modes": {"run": ["{out}/_inputs", "--out", "{out}/p4_explore.json"]}, "default": "run",
     "deps": ["p1", "p1b", "p2", "p3-global", "p3c"], "stage_dir": "_inputs",
     "inputs": {"p2/p2_summary.json": "{up:p2_dir}/p2_summary.json", "p2/p2_events.csv": "{up:p2_dir}/p2_events.csv",
                "extra/p1/p1_events.csv": "{up:p1_events}", "extra/p1b/p1b_events_new.csv": "{up:p1b_dir}/p1b_events_new.csv",
                "p3b_curves.json": "{up:p3b_curves}", "p3c_decision.json": "{up:p3c}"},
     "done": ["p4_explore.json"], "net": False, "token": False, "machine": "laptop",
     "res": (0.5, 1000, 60, 600), "ref_run": "under 10 s", "wind": "inherit", "in_all": True},

    {"id": "p4-mech", "script": "p4_mech.py", "title": "P4 first-round mechanism diagnostics (D-U/D-Z/D-T/D-B/D-P/D-S; criteria fixed before running)",
     "modes": {"full": P1ARGS + ["--p2-dir", "{up:p2_dir}"], "plan": ["--plan"] + P1ARGS + ["--p2-dir", "{up:p2_dir}"],
               "selftest": ["--selftest"]},
     "default": "full", "deps": ["p1", "p1b", "p2"], "done": ["p4_mech_summary.json", "p4_mech_events.csv"],
     "net": True, "token": False, "machine": "workstation", "res": (1.5, 1000, 300, 1800),
     "ref_run": "88 s, peak 0.69 GB", "wind": "direct", "in_all": True},

    {"id": "p4b", "script": "p4b_posthoc.py", "title": "P4b second-round mechanism diagnostics (E1 null, E2 stability correction, E3 deep budget, RIM fine-grid freshwater)",
     "modes": {"full": P1ARGS + ["--p2-dir", "{up:p2_dir}", "--p4-dir", "{up:p4_dir}"], "selftest": ["--selftest"]},
     "default": "full", "deps": ["p4-mech"], "done": ["p4b_summary.json", "p4b_events.csv"],
     "net": True, "token": False, "machine": "workstation", "res": (1.5, 1000, 600, 1800),
     "ref_run": "130 s, peak 0.77 GB", "wind": "direct", "in_all": True},

    {"id": "p4c", "script": "p4c_rim_conservation.py", "title": "P4c RIM-3 freshwater conservation check (synthetic input; needs xarray, pandas)",
     "modes": {"run": ["--out", "{out}/p4c.json"]}, "default": "run", "deps": [], "done": ["p4c.json"],
     "net": False, "token": False, "machine": "laptop", "res": (1.0, 1000, 120, 900),
     "ref_run": "about 20 s", "in_all": True},

    {"id": "p4d", "script": "p4d_global_conservation.py", "title": "P4d uniform approximation of conservation scaling in the global calculation (reads P3B)",
     "modes": {"run": ["{up:p3b_curves}", "{out}/p4d.json"]}, "default": "run", "deps": ["p3-global"], "done": ["p4d.json"],
     "net": False, "token": False, "machine": "laptop", "res": (0.5, 1000, 60, 600),
     "ref_run": "under 1 s; M_RIM and M1 are constants in the script (docs/reproduction.md, known issues)", "wind": "const", "in_all": True},

    {"id": "p4e", "script": "p4e_budget.py", "title": "P4e freshwater budget closure and R1/R2 discrimination",
     "modes": {"full": P1ARGS + ["--p2-dir", "{up:p2_dir}", "--p4-dir", "{up:p4_dir}", "--p4b-dir", "{up:p4b_dir}"],
               "selftest": ["--selftest"]},
     "default": "full", "deps": ["p4b"], "done": ["p4e_summary.json", "p4e_events.csv"],
     "net": True, "token": False, "machine": "workstation", "res": (1.25, 1000, 300, 1200),
     "ref_run": "94 s, peak 0.80 GB", "wind": "direct", "in_all": True},

    {"id": "p4f", "script": "p4f_wind_trend.py", "title": "P4f conservation factor versus wind trend (reads P4E, P4BE; imports p4c)",
     "modes": {"run": ["{up:p4_events}", "{up:p4b_dir}/p4b_events.csv", "--p4sum", "{up:p4_dir}/p4_mech_summary.json",
                       "--out", "{out}/p4f.json"]},
     "default": "run", "deps": ["p4-mech", "p4b"], "done": ["p4f.json"], "net": False, "token": False,
     "machine": "laptop", "res": (0.5, 1000, 60, 600), "ref_run": "about 2 s",
     "wind": "const", "in_all": True},

    {"id": "p5-beta", "script": "p5_beta_mech.py", "title": "P5 r_beta mechanism test (seven candidates; PyCO2SYS carbonate baseline)",
     "modes": {"full": P1ARGS + ["--p2-dir", "{up:p2_dir}", "--p4-events", "{up:p4_events}", "--p3b-curves", "{up:p3b_curves}",
                                 "--p3c", "{up:p3c}"],
               "plan": ["--plan"] + P1ARGS + ["--p2-dir", "{up:p2_dir}", "--p4-events", "{up:p4_events}",
                                              "--p3b-curves", "{up:p3b_curves}", "--p3c", "{up:p3c}"],
               "selftest": ["--selftest"]},
     "default": "full", "deps": ["p2", "p3-global", "p3c", "p4-mech"], "done": ["p5_beta_summary.json", "p5_beta_events.csv"],
     "net": True, "token": False, "machine": "workstation", "res": (1.25, 1000, 360, 1200),
     "ref_run": "209 s, peak 0.56 GB", "wind": "inherit", "in_all": True},

    {"id": "p5-sss", "script": "p5_sss_sat.py", "title": "P5 SMAP RSS L2C V6 satellite SSS collocation test (CMR -> OPeNDAP -> CMORPH -> RIM-3 footprint)",
     "modes": {"full": ["--full", "--cache", "{path:p5_cache}", "--raw-root", "{path:raw_root}", "--rate-mbps", "{opt:raw_rate_mbps}",
                        "--raw-tag", "{tag}", "--token-file", "{path:token_file}"] + P1ARGS,
               "plan": ["--plan", "--cache", "{path:p5_cache}", "--token-file", "{path:token_file}"] + P1ARGS,
               "selftest": ["--selftest"]},
     "default": "full", "deps": ["p1", "p1b"], "ensure_dirs": ["{path:p5_cache}"],
     "done": ["p5_power.json", "p5_summary.json", "p5_events.csv"], "net": True, "token": True, "machine": "workstation",
     "res": (1.5, 1000, 9000, 28800), "ref_run": "8,699 s with downloads (CMORPH 8,186 hours, 13.2 GB at ~1.5 MB/s), peak 1.14 GB; 338 s with the raw archive in place", "wind": "direct", "in_all": True},

    {"id": "raw-spurs2", "script": "raw_fetch_list.py", "title": "Download SPURS-2 originals (8 collections, 48 files; PO.DAAC; bytes and MD5 checked)",
     "modes": {"run": ["--run", "--list", "raw_archive_lists/spurs2.json", "--raw-root", "{path:raw_root}", "--rate", "{opt:raw_rate_mbps}",
                       "--wait-max-s", "21600", "--tag", "{tag}", "--out", "{out}"]},
     "default": "run", "deps": [], "done": ["raw_fetch_status.json"], "net": True, "token": True, "machine": "workstation",
     "res": (0.5, 500, 8100, 32400), "ref_run": "1,274 s, 2.53 GB; 21 s when the files are already archived (verification only)", "in_all": True},

    {"id": "p6-plan", "script": "p6_plan_meta.py", "title": "P6 metadata and availability survey (salinity values not read)",
     "modes": {"run": ["--root", "{raw:spurs2}", "--out", "{out}"]}, "default": "run", "deps": ["raw-spurs2"],
     "done": ["p6_plan.json"], "net": False, "token": False, "machine": "workstation", "res": (2.0, 1000, 300, 1800),
     "ref_run": "seconds", "in_all": False},

    {"id": "p6", "script": "p6_spurs2_profile.py", "title": "P6 SPURS-2 near-surface profile test (salinity snake 0.01-0.02 m versus RIM-3 surface)",
     "modes": {"full": ["--full", "--root", "{raw:spurs2}", "--cache", "{path:p6_cache}", "--out", "{out}"], "selftest": ["--selftest"]},
     "default": "full", "deps": ["raw-spurs2"], "ensure_dirs": ["{path:p6_cache}"],
     "done": ["p6_power.json", "p6_summary.json", "p6_events.csv"], "net": False, "token": False, "machine": "workstation",
     "res": (2.0, 1000, 300, 1800), "ref_run": "4.4 s, peak 0.37 GB",
     "in_all": True},

    {"id": "p7a", "script": "p7a_jpl_smap.py", "title": "P7a JPL SMAP L2B CAP V5 independent-chain collocation test",
     "modes": {"full": ["--full", "--cache", "{path:p7a_cache}", "--p5-cache", "{path:p5_cache}", "--p5-results", "{up:p5_dir}",
                        "--raw-root", "{path:raw_root}", "--rate-mbps", "{opt:raw_rate_mbps}", "--raw-tag", "{tag}",
                        "--token-file", "{path:token_file}"] + P1ARGS,
               "plan": ["--plan", "--cache", "{path:p7a_cache}", "--p5-cache", "{path:p5_cache}", "--p5-results", "{up:p5_dir}",
                        "--raw-root", "{path:raw_root}", "--rate-mbps", "{opt:raw_rate_mbps}", "--raw-tag", "{tag}",
                        "--token-file", "{path:token_file}"] + P1ARGS,
               "selftest": ["--selftest"]},
     "default": "full", "deps": ["p5-sss"], "ensure_dirs": ["{path:p7a_cache}"],
     "done": ["p7a_power.json", "p7a_summary.json", "p7a_events_jpl.csv", "p7a_pairs.csv"], "net": True, "token": True,
     "machine": "workstation", "res": (1.5, 1000, 2400, 7200), "ref_run": "1,944 s with downloads, peak 1.12 GB; 122 s with the raw archive in place", "wind": "direct", "in_all": True},

    {"id": "p7c", "script": "p7c_wind_dependence.py", "title": "P7c wind-dependence analysis (H1-H6; 9 main tests, Bonferroni)",
     "modes": {"full": P1ARGS + ["--p2-dir", "{up:p2_dir}", "--p4-dir", "{up:p4_dir}", "--p4b-dir", "{up:p4b_dir}"],
               "formula": ["--formula"], "selftest": ["--selftest"]},
     "default": "full", "deps": ["p4-mech", "p4b", "p4f"], "done": ["p7c_summary.json", "p7c_events.csv"],
     "net": True, "token": False, "machine": "workstation", "res": (1.25, 1000, 300, 1200),
     "ref_run": "97 s, peak 0.72 GB",
     "wind": "direct", "in_all": True},

    {"id": "p7d", "script": "p7d_p5_timebin_ci.py", "title": "P7d intervals of the P5 ratio by time since rain at overpass (three interval types)",
     "modes": {"full": ["--p5-dir", "{up:p5_dir}"], "selftest": ["--selftest-synthetic"]},
     "default": "full", "deps": ["p5-sss"], "done": ["p7d_summary.json", "p7d_bins.csv"], "net": False, "token": False,
     "machine": "workstation", "res": (1.0, 1000, 60, 600), "ref_run": "1.1 s, peak 117 MB", "wind": "inherit", "in_all": True},

    {"id": "raw-p7e", "script": "p7e_fetch.py", "title": "Download P7e originals (PISTON 40 files, SPURS-1 24 files; NASA ASDC / PO.DAAC)",
     "modes": {"run": ["--run", "--lists", "raw_archive_lists/spurs1.json,raw_archive_lists/piston.json", "--raw-root", "{path:raw_root}",
                       "--rate", "{opt:raw_rate_mbps}", "--wait-max-s", "21600", "--tag", "{tag}", "--out", "{out}"],
               "selftest": ["--selftest"]},
     "default": "run", "deps": [], "done": ["p7e_fetch_rc.json", "spurs1/raw_fetch_status.json", "piston/raw_fetch_status.json"],
     "net": True, "token": True, "machine": "workstation",
     "res": (0.5, 500, 2400, 25200), "ref_run": "exit code 4 in the reference run: 12 sidecar .md5 files disagree with their own CMR checksums, the data files are complete; 42 s when already archived",
     "in_all": True},

    {"id": "p7e-plan", "script": "p7e_plan_meta.py", "title": "P7e metadata survey",
     "modes": {"run": ["--roots", "{raw:piston},{raw:spurs1}", "--out", "{out}"]}, "default": "run", "deps": ["raw-p7e"],
     "done": [], "net": False, "token": False, "machine": "Windows", "res": (2.0, 1000, 300, 1800),
     "ref_run": "11 s", "in_all": False},

    {"id": "p7e", "script": "p7e_piston_spurs1.py", "title": "P7e PISTON 0.35 m and SPURS-1 0.865 m tests",
     "modes": {"full": ["--full", "--piston-root", "{raw:piston}", "--spurs1-root", "{raw:spurs1}", "--out", "{out}"],
               "load-check": ["--load-check", "--piston-root", "{raw:piston}", "--spurs1-root", "{raw:spurs1}", "--out", "{out}"],
               "selftest": ["--selftest", "--out", "{out}"]},
     "default": "full", "deps": ["raw-p7e"],
     "done": ["p7e_power.json", "p7e_summary.json", "p7e_piston_events.csv", "p7e_spurs1_events.csv"],
     "net": False, "token": False, "machine": "Windows", "res": (2.0, 1000, 120, 1800),
     "ref_run": "5.4 s, peak 0.10 GB", "in_all": True},

    {"id": "p8a", "script": "p8a_sat_wind_strata.py",
     "title": "P8a satellite collocation ratios stratified by wind at overpass (TL/TM/TH; reads P5/P7a/P7c; mooring caches read-only, network blocked)",
     "modes": {"full": P8A_ARGS, "check": ["--check"] + P8A_ARGS, "selftest": ["--selftest-synthetic"]},
     "default": "full", "deps": ["p1", "p1b", "p5-sss", "p7a", "p7c"], "done": ["p8a_summary.json", "p8a_bins.csv", "p8a_events_wind.csv"],
     "net": False, "token": False, "machine": "Windows", "res": (2.0, 1000, 600, 1800),
     "ref_run": "90-197 s (rebuilding the mooring series takes ~3 min), peak 0.32 GB", "wind": "direct", "in_all": True},

    {"id": "p8b-plan", "script": "p8b_plan_meta.py", "title": "P8b metadata and availability survey (Lady Amber, Revelle, central mooring; salinity values not read)",
     "modes": {"run": ["--la-root", LA_ROOT, "--spurs2-root", "{raw:spurs2}", "--out", "{out}"]}, "default": "run",
     "deps": ["raw-spurs2"], "done": [], "net": False, "token": False, "machine": "Windows", "res": (3.0, 1000, 300, 1800),
     "ref_run": "35 s", "in_all": False},

    {"id": "p8b", "script": "p8b_ladyamber_layer_ratio.py", "title": "P8b SPURS-2 Lady Amber same-platform surface (~1 cm) to 1 m layer ratio Q_L (R1/R2 discrimination)",
     "modes": {"full": ["--full", "--root", LA_ROOT, "--out", "{out}"],
               "load-check": ["--load-check", "--root", LA_ROOT, "--out", "{out}"], "selftest": ["--selftest", "--out", "{out}"]},
     "default": "full", "deps": ["raw-spurs2"], "done": ["p8b_power.json", "p8b_summary.json", "p8b_events.csv"],
     "net": False, "token": False, "machine": "Windows", "res": (3.0, 1000, 180, 1800),
     "ref_run": "24.4 s, peak 0.63 GB",
     "in_all": True},

    {"id": "piston-track", "script": "piston_track_extent.py", "title": "PISTON hourly ship positions and extent (for the site map; not used in any estimate)",
     "modes": {"run": ["--root", "{raw:piston}", "--out", "{out}"]}, "default": "run", "deps": ["raw-p7e"],
     "done": ["piston_track_extent.json", "piston_track_hourly.csv"], "net": False, "token": False, "machine": "workstation",
     "res": (1.0, 1000, 60, 600), "ref_run": "seconds", "in_all": True},

    {"id": "cal-dw13-meta", "script": "cal_dw13_meta.py",
     "title": "Deployment metadata check (SSS sensor depth in NCEI OCADS deployment XML; group B anemometer heights in OceanSITES; public metadata only)",
     "modes": {"run": ["--p1b-stations", "{up:p1b_dir}/p1b_stations.json", "--p4e", "{up:p4_events}", "--out", "{out}"]},
     "default": "run", "deps": ["p1b", "p4-mech"], "done": ["dw13_meta.json"], "net": True, "token": False,
     "machine": "laptop", "res": (0.5, 500, 900, 3600),
     "ref_run": "about 5 min (network-bound); the 227 metadata originals are kept in <out>/raw_meta/", "in_all": True},

    {"id": "dw13a", "script": "dw13a_windheight.py",
     "title": "Anemometer-height sensitivity: P2 main set with group B heights per source file/variable; all-4 m correctness gate -> measured heights -> comparison",
     "modes": dict(DW13A_STEPS, full={"steps": [DW13A_STEPS[k] for k in ("selftest", "nominal", "gate", "measured", "compare")]}),
     "default": "full", "deps": ["p2", "cal-dw13-meta"],
     "done": ["dw13a_gate.json", "dw13a_compare.json", "measured/p2_summary.json", "measured/dw13a_heights.json"],
     "net": False, "token": False, "machine": "workstation", "res": (2.0, 1000, 2400, 10800),
     "ref_run": "214 s", "wind": "source", "in_all": True},

    {"id": "cal-surface-diag", "script": "cal_surface_diag.py",
     "title": "C1 (tropical-subset satellite / SPURS-2 in situ surface ratio) + C2 (denominator shares of each test)",
     "modes": {"run": ["--sp", "{out}/_inputs", "--out", "{out}"]}, "default": "run", "stage_dir": "_inputs",
     "inputs": {"p5sss-out/p5_events.csv": "{up:p5_dir}/p5_events.csv",
                "p7a-out/p7a_events_jpl.csv": "{up:p7a_dir}/p7a_events_jpl.csv",
                "p6-out/p6_events.csv": "{up:p6_dir}/p6_events.csv",
                "p8a-out/p8a_events_wind.csv": "{up:p8a_dir}/p8a_events_wind.csv",
                "p7e-out/p7e_spurs1_events.csv": "{up:p7e_dir}/p7e_spurs1_events.csv",
                "p8b-out/p8b_events.csv": "{up:p8b_dir}/p8b_events.csv"},
     "deps": ["p5-sss", "p7a", "p6", "p8a", "p7e", "p8b"], "done": ["cal_c1_c2.json"], "net": False, "token": False,
     "machine": "laptop", "res": (1.0, 1000, 120, 900), "ref_run": "under 1 s", "wind": "inherit", "in_all": True},

    {"id": "dw33-bench", "script": "dw33_wind_bins.py", "title": "Wind-binned global flux, bench: one day (2000-01-01), new engine versus P3 engine term by term + annual extrapolation",
     "modes": {"run": ["--mode", "bench"] + DW33_ARGS}, "default": "run", "deps": ["p3-global"], "ensure_dirs": ["{path:dw33_work}"],
     "done": ["dw33_bench.json"], "net": False, "token": False, "machine": "workstation", "res": (5.5, 1000, 900, 3600),
     "ref_run": "one day of the global engine: 1.43 s per time step, peak RSS 3.98 GiB; 10 terms agree with the P3 engine within 1.2e-15", "in_all": False},

    {"id": "dw33", "script": "dw33_wind_bins.py", "title": "Wind-binned global flux: RIM-3 interfacial/dilution flux accumulated by pixel U10 bins (s = 0, 1 and per-pixel check cases; monthly, parallel)",
     "modes": {"full": ["--mode", "full", "--workers", "2"] + DW33_ARGS}, "default": "full", "deps": ["p3-global"],
     "ensure_dirs": ["{path:dw33_work}"], "done": ["dw33_bins.json", "dw33_run.json"], "net": False, "token": False,
     "machine": "workstation", "res": (10.0, 2000, 16000, 47000),
     "ref_run": "14,670 s with 2 workers, up to 5.0 GB per worker", "wind": "const", "in_all": True},

    {"id": "dw33-post", "script": "dw33_post.py", "title": "Wind-binned global flux, post-processing: mooring wind-class bootstrap substitution -> global U, shares, coverage",
     "modes": {"run": ["--bins", "{up:dw33_dir}/dw33_bins.json", "--p4e", "{up:p4_events}", "--p4s", "{up:p4_dir}/p4_mech_summary.json",
                       "--p3b", "{up:p3b_curves}", "--p3c", "{up:p3c}", "--out", "{out}"]},
     "default": "run", "deps": ["dw33", "p4-mech", "p3-global", "p3c"], "done": ["dw33_post.json"],
     "net": False, "token": False, "machine": "laptop", "res": (1.5, 1000, 300, 1200),
     "ref_run": "seconds", "wind": "inherit",
     "in_all": True},

    {"id": "raw-backfill", "script": "raw_archive.py", "title": "Historical backfill of deleted originals (B1 P5 CMORPH, B2 SMAP full granules, B3 P2 CMORPH, B4 ERA5)",
     "modes": {"run": ["--run", "--stages", "B1,B2,B3,B4", "--raw-root", "{path:raw_root}", "--rate-full", "4", "--rate-shared", "2",
                       "--tag", "{tag}", "--out", "{out}"],
               "plan": ["--plan", "--raw-root", "{path:raw_root}", "--out", "{out}"]},
     "default": "plan", "deps": [], "done": [], "net": True, "token": True, "machine": "workstation", "res": (0.5, 1000, 57000, 172800),
     "ref_run": "backfill of originals deleted by earlier script versions; not needed for a re-run", "in_all": False},
]

BY_ID = {s["id"]: s for s in STAGES}


def topo(ids=None):
    """按 deps 的拓扑序返回阶段 id；ids 给定时只含这些（及其顺序），不自动补上游。"""
    want = [s["id"] for s in STAGES] if ids is None else list(ids)
    order, seen, temp = [], set(), set()

    def visit(i):
        if i in seen:
            return
        if i in temp:
            raise ValueError(f"阶段依赖成环：{i}")
        temp.add(i)
        for d in BY_ID[i]["deps"]:
            visit(d)
        temp.discard(i)
        seen.add(i)
        order.append(i)
    for i in want:
        visit(i)
    return [i for i in order if i in want]
