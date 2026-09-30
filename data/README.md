# Derived data

Derived data of this repository, licensed under CC BY 4.0 (see `../LICENSE-DATA`; third-party terms in `../docs/third-party-data.md`).
No third-party raw data are redistributed, with one exception listed below (`products/p2/cmorph_pixels.csv`, a small extract
of NOAA CMORPH CDR values at the mooring pixels). Every file is listed in `SHA256SUMS`.

| Directory | Content |
|---|---|
| `events/` | `events_merged.csv`: one row per rain event (646), all per-event quantities of the mooring test and the mechanism diagnostics |
| `products/<stage>/` | released copies of the reference-run outputs of each pipeline stage (per-event and per-overpass tables, band/basin/month tables, summaries) |
| `products/MANIFEST.tsv` | for every product file: stage, product code, original SHA-256 (as written by the pipeline), released SHA-256, size, treatment |
| `raw_index/raw_manifest.jsonl.gz` | index of the 65,882 original third-party files downloaded by the pipeline (URL, bytes, SHA-256, time) |

## Conventions

- Salinity in psu (practical salinity), salinity changes ΔS in psu (negative = freshening); pCO₂ in µatm; rain in mm (accumulations)
  or mm h⁻¹ (rates); wind U10 in m s⁻¹ at 10 m; depths in m (positive down); times in UTC, ISO 8601 (`2014-11-07T13:00Z`).
- `onset_utc` is the first rain hour t₀ of an event (UTC hour bin [HH:00, HH+1:00)). Events are identified by station + onset_utc;
  `events_merged.csv` also gives a sequential `event_id` (E001-E646, in the order of `p2_events.csv`).
- ΔS(z) of an event = mean over [t₀, t₀+6 h) minus median over [t₀−6 h, t₀). "0.5 m" or `05` is the MAPCO2
  salinity channel (nominal 0.5 m); "1 m" or `1` is the 1 m record (group B: shallowest level ≤ 1.5 m, depth in `s1_depth_m`).
- Model quantities: `rim` = RIM-3 (Witte et al., 2026) evaluated at the given depth, `s20` = S20 surface relation; by default driven by
  CMORPH at the nearest pixel; `_gauge`/`G` = driven by the mooring rain gauge instead.
- Groups: `A_tropical_8` (group A, 8 tropical TAO/RAMA moorings, first round) and `B_subtropical_midlatitude_5` (group B, 5 subtropical and
  midlatitude moorings; code label "O" in the P1b table). Seasons: DJF/MAM/JJA/SON. Station names as in the MAPCO2/OceanSITES files
  (`MOSEAN/WHOTS` = WHOTS; `TAO165E` = 0°, 165°E, etc.).
- Empty cells and `nan` = not available (missing data or not applicable to that event).
- Numbers are as written by the pipeline scripts (no rounding was applied for release).

## Treatment for release

The released products are the outputs of the reference run of each stage, with three kinds of text-only treatment (listed per
file in `products/MANIFEST.tsv`, column `treatment`; numeric values were never changed):

- Chinese text written by the scripts (category labels, notes, interpretation strings, a few dictionary keys) was translated into
  English using one reviewed table (`../docs/translation-table.tsv`, original → released).
- Fields that held absolute paths of local machines were removed from summaries.
- Descriptive text fields and keys are as written by the current scripts, which differ from the reference-run outputs in a few
  descriptive fields (no longer written, or worded differently) and five keys (listed in `../docs/translation-table.tsv`).

Files that needed none of these are byte-identical to the reference-run outputs (`original_sha256` = `released_sha256`). Because
the pipeline writes Chinese labels, a re-run produces files that differ from the released copies in those text fields only;
`code/tools/diff_expected.py` maps the translated strings and keys (see
`../docs/reproduction.md`).

## `events/events_merged.csv` (646 rows × 233 columns)

Built by `code/tools/build_event_table.py` from the released tables below (joined on station + onset_utc; every table matches all 646
events exactly). Columns: `event_id`, `station`, `group`, `onset_utc`, `season`, then prefixed blocks:

| Prefix | Source file | Content (see the file's dictionary below) |
|---|---|---|
| `def_` | `products/p1/p1_events.csv` (group A, `def_round` = 1) or `products/p1b/p1b_events_new.csv` (group B, `def_round` = 2) | event definition quantities: regime, rain, pre-rain medians, observed ΔS, control windows |
| `p2_` | `products/p2/p2_events.csv` | observed and modeled freshening of the main test, pCO₂ change |
| `p4_` | `products/p4-mech/p4_mech_events.csv` | first-round mechanism diagnostics |
| `p4b_` | `products/p4b/p4b_events.csv` | second-round mechanism diagnostics |
| `p5b_` | `products/p5-beta/p5_beta_events.csv` | chemical-slope candidates |
| `p7c_` | `products/p7c/p7c_events.csv` | wind-dependence candidates |
| `dw13a_` | `products/dw13a/measured/dw13a_events_full.csv` | anemometer-height sensitivity (measured heights) |
| `p4e_<window>_` | `products/p4e/p4e_events.csv` (group A only; windows `0_6h`, `6_12h`, `12_24h`, `24_48h`) | freshwater budget per window |

Surface and cross-platform tests are per overpass or per cruise event and are therefore separate tables (`products/p5-sss`, `p6`, `p7a`,
`p7e`, `p8a`, `p8b`).

## Product files (`products/<stage>/`)

Stage names follow `code/stages.py`; product codes (P2S, P4E, ...) are short identifiers used in `products/MANIFEST.tsv`,
`code/expected/products.tsv` and the code comments. Generating script = `code/pipeline/<script>`. Implementation choices referred to
by letter codes (I1-I17, K1-K30, M1-M17, ...) are listed in the header of each script (in Chinese).

### p1 — event count, first round (`p1_events.py`)

`p1_events.csv` (204 events, group A). `station`; `regime` (rain regime defined in advance: western Pacific warm pool / Bay of Bengal
monsoon / central-eastern Pacific); `onset_utc`; `local_solar_hour` (h); `season`; `rain_first_hour_mm` (gauge rain in the onset hour, mm);
`rain_24h_mm` (gauge rain in [t₀, t₀+24 h), mm); `n_rain_hours_24h` (hours ≥ 0.4 mm h⁻¹); `pre_S1_median`, `pre_S05_median` (pre-rain
6 h median salinity at 1 m and at the MAPCO2 channel, psu); `dS1_0_6h`, `dS05_0_6h` (event ΔS, psu); `has_5m` (1 if 5 m salinity is
valid before and after); `n_pco2_in_window` (valid pCO₂ values in the window); `n_ctrl` (control windows, max 5); `ctrl_dS1_mean`
(mean ΔS(1 m) of the control windows, psu); `ctrl_onsets_utc`, `ctrl_dS1` (semicolon-separated control onsets and their ΔS(1 m)).
`p1_summary.json`: gate results D1/D2/D5 of the first round and descriptive counts.

### p1b — sample expansion (`p1b_extend.py`)

`p1b_events_new.csv` (442 events, group B): columns as P1 plus `group` (`O` = group B), `s1_depth_m` (depth of the level used as
"1 m", m) and `s1_is_substitute` (1 if that level is not at 1.0 ± 0.05 m); `regime` = `O:<station>`.
`p1b_stations.json`: screening of all MAPCO2 moorings (included/excluded with reason, levels, resolution).
`p1b_summary.json`: gates D1/D2/D4/D5 of the combined set (646 events).

### p2 — mooring test (`p2_rim_test.py`)

`p2_events.csv` (646 events): `station`, `group`, `onset_utc`, `season`, `rain_24h_mm`; `s1_substitute` (1 = substitute "1 m" level);
`n_ctrl`; `dS05_obs`, `dS1_obs`, `dS5_obs` (observed ΔS at the MAPCO2 channel, 1 m, 5 m; psu); `dS_rim_05`, `dS_rim_1`, `dS_rim_5`
(RIM-3 prediction at 0.5, "1 m" and 5 m, CMORPH-driven; psu); `dS_s20` (S20 surface prediction; psu); `dS_rim_1_gauge`, `dS_s20_gauge`
(gauge-driven; empty in the primary part); `dpCO2_Tnorm` (temperature-normalized pCO₂ change, µatm); `model_reason_cmorph` (why a
model value is missing); `wind_floor_hours` (hours at the 0.1 m s⁻¹ wind floor); `dS_rim_05_k27_explor`, `dS_rim_1_k27_explor`
(exploratory K27 variant, current-term depth factor at t = 1800 s); `model_depth_s1_layer_m` (depth at which the model was evaluated
for the "1 m" level, m).
`p2_summary.json`: main result R_RIM(1 m) with intervals, R_S20, penetration ratio P5, chemical slope β, strata, exit decision.
`p2_d3.json`: discrimination gate D3 (written before any ratio). `p2_glodap.json`: GLODAPv2.2016b surface carbonate values used at
each mooring (nearest 1° cell).
`cmorph_pixels.csv` (14.5 MB, no header; NOAA CMORPH CDR extract): one line per UTC hour and station:
`hour,station,v0,v1,m0,m1,n0,n1` where `hour` = hours since 1970-01-01T00Z, `v0`/`v1` = CMORPH rain rate at the nearest 8 km pixel for
the first/second half-hour (mm h⁻¹), `m0`/`m1` = mean over the 3×3 pixel block (mm h⁻¹), `n0`/`n1` = number of valid pixels in the
block; each hour ends with a line `hour,__DONE__` (or `hour,__MISSING__,reason`). This is the pixel cache read by P2 and all later
mooring stages (`--resume-from`); with it, the mooring stages can be re-run without downloading the ~41 GB of CMORPH files.

### p2-sens — sensitivity part

`p2_sensitivity.json`: the four sensitivity items (windows, thresholds, gauge forcing, cluster choice). Its pixel cache
(26.1 MB) is not included (size limit; see "Not included").

### p3-global — global recalculation (`p3_global.py`)

`selftest.json` (equivalence self-tests against the original capsule), `p3a_repro.json` (reproduction gate P3a), `p3b_curves.json`
(U(s) curves for s = 0-3, whole and history-only scaling, S20, 5 m variant), `p3_run.json` (run record per month).
`p3_bands.csv` (24 rows, 5° latitude bands from `lat_lo` to `lat_hi`, °N): annual totals for 2000 in PgC yr⁻¹; `W_all_PgC` wind-driven
flux; `dep_PgC` wet deposition; `turb_PgC` rain turbulence term; `skin_PgC` skin-salinity constant term (L30); `int|<case>_<s>_PgC`
interfacial term and `gi|<case>_<s>_PgC` its component linear in β (for case `whole` = all RIM-3 terms scaled by s, `hist` = only the
rain-history terms scaled, `s20` = S20, `z5` = 5 m variant; s = scale factor). Sign: positive = outgassing (dep: positive = uptake).
`p3_regions.csv`: the same quantities for s = 1 and S20 by 5° band (`lat_lo`) × WOA09 basin mask code (`basin_code`; −100 = pixel
without a basin code, inferred).
`p3_monthly.csv`: monthly totals (PgC per month) for 2000 (`month` = YYYYMM): `W_all`, `dep`, `turb`, `skin`, `int|whole_1`,
`dil|whole_1` (dilution part), `int|s20_1`, `dil|s20_1`; `n_W`, `n_D` = numbers of pixel half-hours in the wind-driven and deposition
sums. Monthly and band sums reproduce the annual wind-driven flux of −1.442 PgC yr⁻¹.

### p3c — decision (`p3c_decide.py`)

`p3c_decision.json`: substitution of s = R_RIM(1 m) into the curves, gates G-X/G-R/G-N/G-0, robustness variants S1-S7, E1, regime
ratios (`regimes.by_regime`: `warm_pool`, `SPCZ`, `bay_of_bengal`, `central_east_pacific`, `subtropics_midlatitudes`).

### p4-explore — exploratory analyses (`p4_explore.py`)

`p4_explore.json`: few-cluster inference, MAPCO2 channel inversion, stratification, latitude concentration, descriptives.

### p4-mech — first-round mechanism diagnostics (`p4_mech.py`)

`p4_mech_events.csv` (646): `obs<z>` observed ΔS at z = 05 (0.5 m), 1, 5, 10, 20, 25 m (psu; `obs05_m`, `obs1_m` = same-hour paired
versions); `rim<z>` RIM-3 at z = 0, 05, 1, 5, 10, 20, 25 m (psu); `rim1_k27` K27 variant; `rim1_m`, `rim1_m3` 3×3-mean-driven variants;
`rim1_shift1`, `rim1_shift2` model shifted by 1 or 2 half-hours (`shift*_pad_zero` = padded half-hours); `low<z>` wall-layer
(law-of-the-wall) theory at z = 0, 1, 5, 10 m, CMORPH-driven, and `lowg<z>` gauge-driven; `s0_<z>` pre-rain salinity at z (psu); `z1`
depth of the "1 m" level (m); `u10` event mean U10 (m s⁻¹); `zeta` = z1/L (Obukhov); `g6`, `c6` gauge and CMORPH rain in [t₀, t₀+6 h)
(mm); `lam` = ln((g6+0.1)/(c6+0.1)); `rho` = low1/rim1; `e5proxy` = rim1/rim05 (near-surface mixing proxy); `rain24` (mm);
`pcum_cmorph_model`, `pcum_gauge_model`, `pcum_gauge_obs` mean cumulative rain in the window at model / observation times (mm);
`model_ok` (model valid). `p4_mech_summary.json`: the six pillars and the verdict.

### p4b — second-round mechanism diagnostics (`p4b_posthoc.py`)

`p4b_events.csv` (646): `z1`, `g6`, `c6`, `lam`, `obs1`, `rim1`, `low1` as above; `G6sim` gauge rain used by the null simulation (mm);
`rimG1` gauge-driven RIM-3 at the "1 m" level; `fdn1` non-uniform-grid numerical solution at z1; `stab0/1/5/10` stability-corrected
theory (β_h = 5) at 0, z1, 5, 10 m; `stab78_1` with β_h = 7.8; `stabG0`, `stabG1` gauge-driven; `rho_s` = stab1/rim1; `hfine100`,
`hfine10` RIM-3 freshwater content integrated on a fine grid to 100 m and 10 m (m of freshwater); `h3pt10` three-point version to 10 m;
`pcum_cmorph_model` (mm). `p4b_e1_null.csv` (6,000 rows): null simulation E1 — `setting`, simulation index `sim`, ratios by λ class
(`R0`, `R1`, `R2`) and overall (`Rall`), normalized class ratios (`n0`-`n2`) and the shape statistic `D`. `p4b_summary.json`: E1-E3.

### p4c, p4d, p4f — model derivations

`p4c.json` (RIM-3 freshwater conservation from synthetic input), `p4d.json` (uniform conservation scaling of the global curve),
`p4f.json` (conservation factor versus wind trend).

### p4e — freshwater budget (`p4e_budget.py`)

`p4e_events.csv` (387 rows = group A events × windows 0-6, 6-12, 12-24, 24-48 h): `window` (time window after t₀); `D` deepest level used (10/20/25 m);
`has20` (20 m available); `h10`, `hD` observed column freshwater height to 10 m and to D (m of freshwater; closure ratio
Q = 1000·h/rain); `h05D` the same starting at 0.5 m; `rain` gauge rain accumulated since t₀ (mm); `f05`, `f1` fractional freshening −ΔS/S₀ at 0.5 m and 1 m; `fR0g`, `fR0c` RIM-3 surface fractional
freshening, gauge- and CMORPH-driven. `p4e_summary.json`: closure ratios and the R1/R2 test.

### p5-beta — chemical slope (`p5_beta_mech.py`)

`p5_beta_events.csv` (646): `dS05`, `dS1` (psu), `dp` (pCO₂ change, µatm) and window variants `_0-3`, `_3-6`, `_6-12` (h); `dp_raw`
(not temperature-normalized); `dp_c0413`, `dp_c0433` (alternative temperature coefficients); `dp_isochem`; `dT` (K); `S0_05`, `T_pre`
(°C), `pco2_pre`, `pco2_air`, `dpco2_pre` (µatm); `u10`; `g6`, `c6`; control windows `n_ctrl_valid`, `ctrl_mean_x` (ΔS), `ctrl_mean_y`
(ΔpCO₂); model slopes `bm*` (µatm psu⁻¹) for the candidates (E1 local, E2 rain DIC 0/12/50, E3 constants, E4 local baseline, E4b);
`E4_ok`; `gx_central`, `gx_upper` gas-exchange estimates; `rain_max6`; `rain_insufficient`; `s1_substitute`; `_glodap_pco2`; `in_D4`
(1 = strong-freshening subset). `i` = internal row index. `p5_beta_summary.json`: candidates (a)-(g) and their overall mapping.

### p5-sss — RSS SMAP collocation (`p5_sss_sat.py`)

`p5_events.csv` (87 events with a valid post-rain overpass): `post_gid` (SMAP granule), `post_utc`, `dt_post_h` (hours after t₀);
`n_ref` reference overpasses; `node` (am/pm); `n_looks_post`; `sat_post`, `sat_ref` satellite SSS (psu); `y` = sat_post − sat_ref
(psu); `F_post`, `F_ref` footprint-mean RIM-3 dilution factor; `S0` reference salinity; `x` = S0·(F_post − F_ref), the RIM-3 surface
prediction (psu); `x_pix0` central-pixel version; `sig` formal uncertainty of y (psu); `imerg_c_post` IMERG rain at the central cell
(mm h⁻¹); `rain_foot_cmorph_post` footprint CMORPH rain; `dropped_w_post` weight of dropped pixels; `obs1`, `obs05`, `rim1`, `rim05`
same-time mooring values (psu); `s1_depth`. `p5_power.json` (written before the ratio), `p5_summary.json` (R_sfc and sensitivity).

### p6 — SPURS-2 in situ surface (`p6_spurs2_profile.py`)

`p6_events.csv` (29 events): `cruise`, `onset`, `acc24` (24 h rain, mm), `cluster` (cruise × ISO week), `disp_km` (ship displacement);
per level `<L>_obs` (observed ΔS), `<L>_rim` (RIM-3), `<L>_s0`, `<L>_npre`, `<L>_npost` for L = `snake` (salinity snake, 0.01-0.02 m),
`ssp005`...`ssp110` (profiler at 0.05-1.10 m), `usps2`, `usps3` (2 m, 3 m), `tsg5` (thermosalinograph, 5 m), `snake_met` (snake with
the meteorological file's surface salinity); `rain_hours_0_12`.

### p7a — JPL SMAP independent chain (`p7a_jpl_smap.py`)

`p7a_events_jpl.csv` (96): as `p5_events.csv` for the JPL product. `p7a_pairs.csv` (85 same-orbit pairs): `post_rev` orbit, `y_R`,
`x_R` (RSS), `y_J`, `x_J` (JPL), `x_J60` (JPL model with a 60 km kernel), `d` = y_R − y_J, `imerg_post`, `cmorph_foot_post_J`.

### p7c — wind dependence (`p7c_wind_dependence.py`)

`p7c_events.csv` (646): `z1`, `sub` (substitute level), `u10`, `rain24`, `g6`, `c6`, `obs1`, `obs5`, `rim0/1/5`; `rim1_lin` linearized
RIM-3; `G0/G1/G5` conservative Gaussian model; `_w03`, `_w612` = windows 0-3 h and 6-12 h; `rimG1` gauge-driven; `M` RIM-3 freshwater
multiplier; `ctrl_mean1`, `n_ctrl_valid`; `pretrend1` pre-rain trend at 1 m; `s0_1`; `pcum_cmorph_model`, `pcum_gauge_obs`.
`p7c_summary.json`: H1-H6 and the overall mapping. `p7c_formula.json`: formula analysis on synthetic events.

### p7d — satellite ratio by time since rain

`p7d_bins.csv` (14 classes): `key`, `kind`, `label`, `n`, cluster counts, `sum_num`, `sum_den`, `R` (ratio), interval bounds for
station × season (`ss_`), station (`st_`), jackknife-t (`jk_`) and leave-one-out range (`loo_`), `frac_den_nonneg`, `p5_gate_met`,
`category` (interpretation class), `station_same_side`, `jackknife_same_side`.

### p7e — PISTON and SPURS-1 (`p7e_piston_spurs1.py`)

`p7e_power.json` (written before the ratios), `p7e_summary.json`; `p7e_piston_events.csv` (22) and `p7e_spurs1_events.csv` (13): as P6 with levels `so_top` (SurfOtter top, 0.35 m), `so_1m`, `tsg5`
(PISTON) and `s086` (0.865 m), `s21` (2.1 m), `tsg5` (SPURS-1); `<L>_z` level depth (m).

### p8a — satellite ratios by wind at overpass (`p8a_sat_wind_strata.py`)

`p8a_events_wind.csv` (183): `product` (RSS/JPL), `u_tr` (mooring U10 at overpass), `u10_event` (event mean), classes `bin_t`
(terciles L/M/H at 6.03/9.33 m s⁻¹), `bin_8` (≤ 8 / > 8 m s⁻¹), `bin_e` (event-level). `p8a_bins.csv`: ratios per class as in P7d plus
`slope` with bounds.

### p8b — Lady Amber layer ratio (`p8b_ladyamber_layer_ratio.py`)

`p8b_events.csv` (30): as P6 with `lat` and levels `snake`, `la1m` (Lady Amber 1 m), `tsg5`.

### Other

`piston-track/piston_track_hourly.csv` (hourly ship positions of the two PISTON cruises: `cruise`, `i`, `lat`, `lon`) and
`piston_track_extent.json`. `cal-dw13-meta/dw13_meta.json` (sensor depth and anemometer height from public deployment metadata) and
`raw_meta_SHA256SUMS.txt` (checksums of the 227 metadata files read; the files themselves are not redistributed).
`cal-surface-diag/cal_c1_c2.json` (C1 tropical satellite/in situ ratio, C2 denominator shares). `dw33/dw33_bins.json` (global fluxes
accumulated by pixel wind bins, gC) and `dw33_run.json`; `dw33-post/dw33_post.json` (wind-class substitution). `dw13a/`: comparison,
correctness gate, and the main-set outputs recomputed with measured anemometer heights (`measured/`) and with 4 m everywhere
(`nominal/p2_summary.json`; `nominal/p2_events.csv` is byte-identical to `p2/p2_events.csv` and is not repeated);
`dw13a_events_full.csv`: per-event anemometer height at onset (`z_onset_m`), kinds of height sources in the window, and the
recomputed model values.

## `raw_index/raw_manifest.jsonl.gz`

One JSON object per original third-party file (65,882): `dataset`, `url`, `path` (relative path in the raw archive), `bytes`,
`sha256`, `status` (`ok` downloaded, `copied` from an earlier cache, `md5_mismatch` only for 12 sidecar `.md5` files whose own CMR
checksum differed), `t` (UTC time), `source_job` (the pipeline step that archived the file); some rows add `md5`, `granule`,
`collection`, `short_name`, `etag`, `expect_*`. The field giving the local source path of copied files and the run identifiers in
`source_job` and in four log paths were removed. Use it to check that a re-download is the same version.

## Not included (size limit 20 MB or not derived)

| File | Size | SHA-256 | How to obtain |
|---|---|---|---|
| P2 sensitivity pixel cache `cmorph_pixels.csv` (stage p2-sens) | 26.1 MB | `5de465db22453175eefda824dc759dc5f097ebac91da4a47f478426145837b67` | re-run stage `p2-sens` (downloads ~33 GB of CMORPH) |
| Original raw-download index `manifest.jsonl` | 36.0 MB | `194d00c185ddc7a8c4a18dbd294ceb879e14bcb424f4528bc59a27209d4f07e5` | released here without the local-path field and gzipped (`raw_index/`) |
| P3 and dw33 monthly intermediates and static fields | several GB | recorded in the run records | re-run `p3-global` / `dw33` |
| P5/P6/P7a caches (satellite subsets, CMORPH windows) | GB | — | re-run the stages (`raw_index` pins the originals) |
| Third-party originals (CMORPH, ERA5, HYCOM, OISST, SMAP, SPURS, PISTON, moorings, GLODAP, WOA09, Watson, metadata XML) | ~160 GB | `raw_index/raw_manifest.jsonl.gz` | download from the providers (`../docs/third-party-data.md`) |
| Metadata-survey outputs (`p6_plan.json` 10 MB etc.) and logs | — | — | run records; not needed for any result |
