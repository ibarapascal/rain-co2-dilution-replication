"""p4_explore.py — 探索性补分析；不进任何判门。

功能：只读已有的小产物（P2S、P2E、P1E、P1bE、P3B、P3C），用 numpy 秒级算不需要原始数据的补分析；
输出 p4_explore.json。E1、E9、E10 是描述量或既有曲线上的算术。

用法：python3 p4_explore.py [IN_DIR] [--out p4_explore.json]
  IN_DIR 须含 p2/p2_summary.json、p2/p2_events.csv、extra/p1/p1_events.csv、extra/p1b/p1b_events_new.csv、
  p3b_curves.json、p3c_decision.json（run.py 摆放）。
bootstrap 与 P2 同算法（p2_rim_test.cluster_boot：簇键排序、PCG64(20260926)、integers(0,K,(B,K))、Σnum/Σden、百分位）。
每项分析的输入字段与输出键：
  E0 validate          P2E dS1_obs,dS_rim_1,dS05_obs,dpCO2_Tnorm,dS5_obs,dS_rim_5；站×季簇 → 复算 R_RIM(1 m)、β_obs、P5 的 95% CI，
                       与 P2S 比对（应逐位一致，验证本脚本的抽样实现）                         → validate.*
  E1 descriptives      P2E dS1_obs、rain_24h_mm → 646 事件 1 m 观测 ΔS 均值、24 h 雨量中位数→ descriptives.*
  E2 few_clusters      同 E0 的量，改以站为簇（12 簇）的 bootstrap；另报站×季与站两种簇的 delete-one-cluster
                       jackknife t 区间（df＝G−1）；逐量报 n 与 G                                   → few_clusters.*
  E3 loso              R_RIM(1 m) 逐站留一                                                           → loso.*
  E4 channel_05m       P2E dS05_obs,dS1_obs,dS_rim_05,dS_rim_1,dS_s20,dpCO2_Tnorm；P1E/P1bE pre_S05_median,pre_S1_median：
                       观测与模型的 Σ0.5 m／Σ1 m；「同传感器样」事件（|雨前 0.5−1 m 中位差| ≤0.003 且 |ΔS05−ΔS1| ≤0.003）
                       的比例（分群）及剔除后的比值、R_S20(0.5 m)、β_obs                           → channel_05m.*
  E5 mixing_proxy      P2E dS_rim_1,dS_rim_05,dS_rim_1_k27_explor,model_depth_s1_layer_m：模型 1 m／0.5 m 衰减比
                       （Kz∝U² 的单调代理；只取 1.0 m 层且模型 0.5 m 淡化 ≤ −0.005 psu 的事件）三分位上的 R_RIM(1 m)
                       与 E1 变体比值                                                                 → mixing_proxy.*
  E6 rain_terciles     P2E rain_24h_mm 三分位上的 R_RIM(1 m)                                           → rain_terciles.*
  E7 joint_beta        同一次站×季重抽里同时算 s＝R_RIM(1 m)、β_obs（D4 子集）、β_model（按观测 ΔS05 加权，P2S
                       beta_model_by_station pCO₂ 斜率）→ r_β 的区间与 S6 的 U 在 (s, r_β) 联合不确定度下的区间；
                       U_S6 按 P3C 同式 U＝Dep＋[I(s)＋(r−1)G(s)](1−c)（P3B whole 曲线）             → joint_beta.*
  E8 p5_power          P5 观测−模型配对差的 bootstrap SE 与 MDE（2.80×SE）                           → p5_power.*
  E9 decomposition     P3B/P3C：U(s)＝Dep＋常数项＋随 s 变化的稀释分量；s＝1 与 M1 的分量及其相对变化；
                       U、Δ 的 PgC 值及占 GCB 2024 海洋碳汇 2.9 PgC yr⁻¹ 的比例                        → decomposition.*
  E10 thresholds       G-R 带 ±2.5 pp 按实算斜率换算的 |s−1|；SESOI 0.35 按实算曲线对应的 U 与总雨效变动  → thresholds.*
  E11 lat_concentration P3B bands whole_I_int_PgC：M1 下 Δ（界面项部分）落在 |纬度| ≤10°、≤20° 的比例；s＝1 界面项同口径 → lat.*
  E12 min_effect       E0 重抽中 R_RIM(1 m) ≥0.65 的次数（最小效应检验 H0: R ≥ 0.65 的 bootstrap p）  → min_effect.*
固定 seed 20260926；全部运行 <10 s、内存 <200 MB。
Change log：2026-09-27 初版。
"""
import csv
import json
import math
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
DEF_IN = os.path.join(_rp.path("out_root"), "p4-explore", "_inputs")  # [repro] 输入目录由 run.py 摆放
args = [a for a in sys.argv[1:] if not a.startswith("--")]
IN = os.path.abspath(args[0]) if args else DEF_IN
OUT = os.path.join(_rp.path("out_root"), "p4-explore", "p4_explore.json")  # [repro] 输出写阶段目录
if "--out" in sys.argv:
    OUT = sys.argv[sys.argv.index("--out") + 1]
B, SEED, MDE_Z = 10000, 20260926, 2.80
GCB_SOCEAN = 2.9  # 文献常数：GCB 2024 的 2023 年 S_OCEAN

S = json.load(open(os.path.join(IN, "p2", "p2_summary.json")))
C = json.load(open(os.path.join(IN, "p3b_curves.json")))
P3C = json.load(open(os.path.join(IN, "p3c_decision.json")))  # [repro] 由 run.py 摆进输入目录
EV = list(csv.DictReader(open(os.path.join(IN, "p2", "p2_events.csv"))))
PRE = {}
for f in ("extra/p1/p1_events.csv", "extra/p1b/p1b_events_new.csv"):
    for r in csv.DictReader(open(os.path.join(IN, f))):
        PRE[(r["station"], r["onset_utc"])] = r


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def r5(x, nd=5):
    return None if x is None or not math.isfinite(x) else round(float(x), nd)


for r in EV:
    for k in ("dS05_obs", "dS1_obs", "dS5_obs", "dS_rim_05", "dS_rim_1", "dS_rim_5", "dS_s20", "dpCO2_Tnorm",
              "rain_24h_mm", "dS_rim_1_k27_explor", "dS_rim_05_k27_explor", "model_depth_s1_layer_m"):
        r[k] = fnum(r[k])
    p = PRE[(r["station"], r["onset_utc"])]
    r["pre05"], r["pre1"] = fnum(p["pre_S05_median"]), fnum(p["pre_S1_median"])
    r["isA"] = r["group"].startswith("A")


def rows(num, den, pred=lambda r: True):
    return [r for r in EV if r[num] is not None and r[den] is not None and pred(r)]


def draws(K, seed=SEED):
    return np.random.Generator(np.random.PCG64(seed)).integers(0, K, size=(B, K))


def cluster_ratio(rs, num, den, level="season", seed=SEED, extra=None):
    """P2 同算法的簇 bootstrap；level＝season（站×季）或 station（站）。extra＝(num2, den2) 共用抽样的配对量。"""
    keys = [(r["station"], r["season"]) if level == "season" else r["station"] for r in rs]
    uk = sorted(set(keys))
    kid = {k: j for j, k in enumerate(uk)}
    K = len(uk)
    idx = np.array([kid[k] for k in keys], int)
    sn = np.bincount(idx, weights=[r[num] for r in rs], minlength=K)
    sd = np.bincount(idx, weights=[r[den] for r in rs], minlength=K)
    d = draws(K, seed)
    with np.errstate(divide="ignore", invalid="ignore"):
        b = sn[d].sum(1) / sd[d].sum(1)
    pt = sn.sum() / sd.sum()
    bf = b[np.isfinite(b)]
    out = {"n": len(rs), "G": K, "point": r5(pt), "ci95": [r5(np.percentile(bf, 2.5)), r5(np.percentile(bf, 97.5))],
           "se_boot": r5(np.std(bf, ddof=1))}
    # delete-one-cluster jackknife，t(G−1) 区间
    th = np.array([(sn.sum() - sn[g]) / (sd.sum() - sd[g]) for g in range(K)])
    se_j = math.sqrt((K - 1) / K * np.sum((th - th.mean()) ** 2))
    tq = stats.t.ppf(0.975, K - 1)
    out["jackknife_t"] = {"se": r5(se_j), "ci95": [r5(pt - tq * se_j), r5(pt + tq * se_j)], "df": K - 1}
    if extra:
        sn2 = np.bincount(idx, weights=[r[extra[0]] for r in rs], minlength=K)
        sd2 = np.bincount(idx, weights=[r[extra[1]] for r in rs], minlength=K)
        with np.errstate(divide="ignore", invalid="ignore"):
            b2 = sn2[d].sum(1) / sd2[d].sum(1)
        diff = (b - b2)[np.isfinite(b - b2)]
        out["paired"] = {"point": r5(sn2.sum() / sd2.sum()), "ci95": [r5(np.percentile(b2, 2.5)), r5(np.percentile(b2, 97.5))],
                         "diff_ci95": [r5(np.percentile(diff, 2.5)), r5(np.percentile(diff, 97.5))],
                         "diff_se": r5(np.std(diff, ddof=1)), "diff_mde": r5(MDE_Z * np.std(diff, ddof=1))}
    out["_boot"] = b
    return out


def strip(d):
    return {k: (strip(v) if isinstance(v, dict) else v) for k, v in d.items() if not k.startswith("_")}


R = {"script": "src/p4_explore.py", "status": "探索性（不进判门）", "B": B, "seed": SEED, "inputs_dir": IN}

# ---------------- E0 validate
r1 = rows("dS1_obs", "dS_rim_1")
r05 = rows("dS05_obs", "dS_rim_05")
d4 = [r for r in EV if r["dS05_obs"] is not None and r["dS05_obs"] <= -0.2 and r["dpCO2_Tnorm"] is not None]
r5m = rows("dS5_obs", "dS1_obs", lambda r: r["dS_rim_5"] is not None and r["dS_rim_1"] is not None)
v_R = cluster_ratio(r1, "dS1_obs", "dS_rim_1")
v_b = cluster_ratio(d4, "dpCO2_Tnorm", "dS05_obs")
v_P5 = cluster_ratio(r5m, "dS5_obs", "dS1_obs", extra=("dS_rim_5", "dS_rim_1"))
pr = S["primary"]
R["validate"] = {
    "R_RIM_1m": {"this": strip(v_R)["ci95"], "p2s": pr["R_RIM_1m"]["ci95"], "n": v_R["n"], "G": v_R["G"]},
    "beta_D4": {"this": strip(v_b)["ci95"], "p2s": pr["beta_D4subset"]["ci95"], "n": v_b["n"], "G": v_b["G"]},
    "P5": {"this": strip(v_P5)["ci95"], "p2s": pr["P5"]["ci95"], "n": v_P5["n"], "G": v_P5["G"],
           "diff_this": v_P5["paired"]["diff_ci95"], "diff_p2s": pr["P5"]["P5_RIM3"]["diff_main_minus_this_ci95"]},
}

# ---------------- E1 descriptives
R["descriptives"] = {"mean_dS1_obs_646": r5(np.mean([r["dS1_obs"] for r in r1]), 6), "n": len(r1),
                     "rain24_median_mm": r5(np.median([r["rain_24h_mm"] for r in EV]), 3)}

# ---------------- E2 few clusters
fc = {}
for nm, (rs, a, b_, ex) in {"R_RIM_1m": (r1, "dS1_obs", "dS_rim_1", None),
                            "beta_D4": (d4, "dpCO2_Tnorm", "dS05_obs", None),
                            "P5_obs": (r5m, "dS5_obs", "dS1_obs", ("dS_rim_5", "dS_rim_1")),
                            "R_S20_05m": (rows("dS05_obs", "dS_s20"), "dS05_obs", "dS_s20", None)}.items():
    fc[nm] = {"station_x_season": strip(cluster_ratio(rs, a, b_, "season", extra=ex)),
              "station": strip(cluster_ratio(rs, a, b_, "station", extra=ex))}
for g, pred in (("A", lambda r: r["isA"]), ("B", lambda r: not r["isA"])):
    fc["R_RIM_1m_" + g] = {"station": strip(cluster_ratio([r for r in r1 if pred(r)], "dS1_obs", "dS_rim_1", "station"))}
    fc["beta_D4_" + g] = {"station_x_season": strip(cluster_ratio([r for r in d4 if pred(r)], "dpCO2_Tnorm", "dS05_obs")),
                          "station": strip(cluster_ratio([r for r in d4 if pred(r)], "dpCO2_Tnorm", "dS05_obs", "station"))}
R["few_clusters"] = fc

# ---------------- E3 LOSO
loso = {}
for st in sorted({r["station"] for r in r1}):
    sub = [r for r in r1 if r["station"] != st]
    loso[st] = r5(sum(r["dS1_obs"] for r in sub) / sum(r["dS_rim_1"] for r in sub))
R["loso"] = {"by_left_out_station": loso, "min": min(loso.values()), "max": max(loso.values())}

# ---------------- E4 0.5 m channel
both = [r for r in EV if None not in (r["dS05_obs"], r["dS1_obs"], r["dS_rim_05"], r["dS_rim_1"])]


def same_like(r):
    return (r["pre05"] is not None and r["pre1"] is not None and abs(r["pre05"] - r["pre1"]) <= 0.003
            and abs(r["dS05_obs"] - r["dS1_obs"]) <= 0.003)


def ratio(rs, a, b_):
    den = sum(r[b_] for r in rs)
    return r5(sum(r[a] for r in rs) / den) if rs and den else None


ch = {"n_both": len(both),
      "obs_05_over_1": ratio(both, "dS05_obs", "dS1_obs"), "model_05_over_1": ratio(both, "dS_rim_05", "dS_rim_1"),
      "obs_mean_05": r5(np.mean([r["dS05_obs"] for r in both])), "obs_mean_1": r5(np.mean([r["dS1_obs"] for r in both])),
      "obs_median_05": r5(np.median([r["dS05_obs"] for r in both])), "obs_median_1": r5(np.median([r["dS1_obs"] for r in both])),
      "by_group": {}}
for g, pred in (("A", lambda r: r["isA"]), ("B", lambda r: not r["isA"])):
    sub = [r for r in both if pred(r)]
    fl = [r for r in sub if same_like(r)]
    keep = [r for r in sub if not same_like(r)]
    ch["by_group"][g] = {"n": len(sub), "n_same_like": len(fl), "frac_same_like": r5(len(fl) / len(sub), 3),
                         "obs_05_over_1": ratio(sub, "dS05_obs", "dS1_obs"), "model_05_over_1": ratio(sub, "dS_rim_05", "dS_rim_1"),
                         "obs_05_over_1_excl_same": ratio(keep, "dS05_obs", "dS1_obs")}
fl_all = [r for r in both if same_like(r)]
keep_all = [r for r in both if not same_like(r)]
ch["n_same_like"] = len(fl_all)
ch["frac_same_like"] = r5(len(fl_all) / len(both), 3)
ch["obs_05_over_1_excl_same"] = ratio(keep_all, "dS05_obs", "dS1_obs")
ch["model_05_over_1_excl_same"] = ratio(keep_all, "dS_rim_05", "dS_rim_1")
keep_ids = {id(r) for r in EV if not (r["dS05_obs"] is not None and r["dS1_obs"] is not None and same_like(r))}
rr = cluster_ratio([r for r in r1 if id(r) in keep_ids], "dS1_obs", "dS_rim_1")
ch["R_RIM_1m_excl_same"] = {k: strip(rr)[k] for k in ("n", "G", "point", "ci95")}
rs20 = cluster_ratio([r for r in rows("dS05_obs", "dS_s20") if id(r) in keep_ids], "dS05_obs", "dS_s20")
ch["R_S20_05m_excl_same"] = {k: strip(rs20)[k] for k in ("n", "G", "point", "ci95")}
d4k = [r for r in d4 if id(r) in keep_ids]
bk = cluster_ratio(d4k, "dpCO2_Tnorm", "dS05_obs")
ch["beta_D4_excl_same"] = {k: strip(bk)[k] for k in ("n", "G", "point", "ci95")}
ch["d4_same_like_by_group"] = {"A": sum(1 for r in d4 if r["isA"] and id(r) not in keep_ids),
                               "B": sum(1 for r in d4 if (not r["isA"]) and id(r) not in keep_ids)}
R["channel_05m"] = ch

# ---------------- E5 mixing proxy（模型 1 m／0.5 m 衰减比三分位）
px = [r for r in r1 if r["model_depth_s1_layer_m"] == 1.0 and r["dS_rim_05"] is not None and r["dS_rim_05"] <= -0.005]
for r in px:
    r["q_mod"] = r["dS_rim_1"] / r["dS_rim_05"]
qs = np.array([r["q_mod"] for r in px])
cuts = np.percentile(qs, [100 / 3, 200 / 3])
mp = {"n": len(px), "n_excluded_small_model": sum(1 for r in r1 if r["model_depth_s1_layer_m"] == 1.0) - len(px),
      "tercile_cuts": [r5(c, 4) for c in cuts], "terciles": {}}
tot_den = sum(r["dS_rim_1"] for r in px)
for i, (lo, hi) in enumerate(((-np.inf, cuts[0]), (cuts[0], cuts[1]), (cuts[1], np.inf))):
    sub = [r for r in px if lo < r["q_mod"] <= hi]
    cr = strip(cluster_ratio(sub, "dS1_obs", "dS_rim_1"))
    k27 = [r for r in sub if r["dS_rim_1_k27_explor"] is not None]
    ck = strip(cluster_ratio(k27, "dS1_obs", "dS_rim_1_k27_explor"))
    mp["terciles"]["T%d" % (i + 1)] = {"q_mod_median": r5(np.median([r["q_mod"] for r in sub]), 3), "n": cr["n"], "G": cr["G"],
                                       "R_RIM_1m": cr["point"], "ci95": cr["ci95"],
                                       "R_E1_1m": ck["point"], "E1_ci95": ck["ci95"],
                                       "share_model_den": r5(sum(r["dS_rim_1"] for r in sub) / tot_den, 3),
                                       "frac_A": r5(sum(r["isA"] for r in sub) / len(sub), 3)}
R["mixing_proxy"] = mp

# ---------------- E6 rain terciles
rc = np.percentile([r["rain_24h_mm"] for r in r1], [100 / 3, 200 / 3])
rt = {"cuts_mm": [r5(c, 2) for c in rc], "terciles": {}}
for i, (lo, hi) in enumerate(((-np.inf, rc[0]), (rc[0], rc[1]), (rc[1], np.inf))):
    cr = strip(cluster_ratio([r for r in r1 if lo < r["rain_24h_mm"] <= hi], "dS1_obs", "dS_rim_1"))
    rt["terciles"]["T%d" % (i + 1)] = {k: cr[k] for k in ("n", "G", "point", "ci95")}
R["rain_terciles"] = rt

# ---------------- E7 joint (s, r_beta)
bm_st = {k: v["pCO2_slope"] for k, v in S["beta_model_by_station"].items()}
allr = [r for r in EV]
keys = sorted({(r["station"], r["season"]) for r in allr})
kid = {k: j for j, k in enumerate(keys)}
K = len(keys)


def binc(rs, f):
    return np.bincount([kid[(r["station"], r["season"])] for r in rs], weights=[f(r) for r in rs], minlength=K)


sn1, sd1 = binc(r1, lambda r: r["dS1_obs"]), binc(r1, lambda r: r["dS_rim_1"])
bn, bd = binc(d4, lambda r: r["dpCO2_Tnorm"]), binc(d4, lambda r: r["dS05_obs"])
bmn = binc(d4, lambda r: bm_st[r["station"]] * r["dS05_obs"])
d = draws(K)
with np.errstate(divide="ignore", invalid="ignore"):
    sb = sn1[d].sum(1) / sd1[d].sum(1)
    bo = bn[d].sum(1) / bd[d].sum(1)
    bmod = bmn[d].sum(1) / bd[d].sum(1)
    rb = bo / bmod
sg = np.array(C["global"]["whole"]["s"])
Iw = np.array(C["global"]["whole"]["I_int_pct"])
Gw = np.array(C["global"]["whole"]["G_int_pct"])
c, dep = C["c"], C["deposition_pct"]


def U_S6(s, r):
    return dep + (np.interp(s, sg, Iw) + (r - 1) * np.interp(s, sg, Gw)) * (1 - c)


ok = np.isfinite(sb) & np.isfinite(rb)
u6 = U_S6(np.clip(sb[ok], 0, 3), rb[ok])
rb_pt = (bn.sum() / bd.sum()) / (bmn.sum() / bd.sum())
s_pt = sn1.sum() / sd1.sum()
R["joint_beta"] = {"K_union_clusters": K, "n_s": len(r1), "n_beta": len(d4), "G_beta": int(np.count_nonzero(bd)),
                   "r_beta_point": r5(rb_pt), "r_beta_ci95": [r5(np.percentile(rb[ok], 2.5)), r5(np.percentile(rb[ok], 97.5))],
                   "beta_model_pooled_point": r5(bmn.sum() / bd.sum()),
                   "S6_U_point": r5(U_S6(s_pt, rb_pt)),
                   "S6_U_joint_ci95": [r5(np.percentile(u6, 2.5)), r5(np.percentile(u6, 97.5))],
                   "corr_s_rbeta": r5(np.corrcoef(sb[ok], rb[ok])[0, 1], 3), "n_draws_used": int(ok.sum())}

# ---------------- E8 P5 power
R["p5_power"] = {"n": v_P5["n"], "G": v_P5["G"], "diff_point": r5(v_P5["point"] - v_P5["paired"]["point"]),
                 "diff_se": v_P5["paired"]["diff_se"], "diff_mde": v_P5["paired"]["diff_mde"],
                 "diff_mde_rel_to_P5_RIM": r5(v_P5["paired"]["diff_mde"] / v_P5["paired"]["point"], 3)}

# ---------------- E9 decomposition（P3B 曲线上的算术）
U = np.array(C["global"]["whole"]["U_pp"])
U0, U1 = float(np.interp(0, sg, U)), float(np.interp(1, sg, U))
M1 = P3C["M1"]
Um, Uci = M1["U_curve"], M1["U_curve_ci"]
W = abs(C["wind_PgC"])
dec = {"Dep": r5(dep, 4), "U0": r5(U0, 4), "const_term_U0_minus_Dep": r5(U0 - dep, 4),
       "dil_component_s1": r5(U1 - U0, 4), "dil_component_M1": r5(Um - U0, 4),
       "dil_component_M1_ci": [r5(Uci[0] - U0, 4), r5(Uci[1] - U0, 4)],
       "dil_component_rel_change": r5((Um - U0) / (U1 - U0) - 1, 4),
       "dil_component_rel_change_ci": [r5((Uci[0] - U0) / (U1 - U0) - 1, 4), r5((Uci[1] - U0) / (U1 - U0) - 1, 4)],
       "share_U1": {"Dep": r5(dep / U1, 3), "const": r5((U0 - dep) / U1, 3), "dil": r5((U1 - U0) / U1, 3)},
       "share_UM1": {"Dep": r5(dep / Um, 3), "const": r5((U0 - dep) / Um, 3), "dil": r5((Um - U0) / Um, 3)},
       "PgC": {"U1": r5(U1 / 100 * W, 4), "U_M1": r5(Um / 100 * W, 4), "Delta": r5(M1["Delta"] / 100 * W, 4),
               "Delta_ci": [r5(M1["Delta_ci"][0] / 100 * W, 4), r5(M1["Delta_ci"][1] / 100 * W, 4)]},
       "frac_of_GCB_SOCEAN": {"U1": r5(U1 / 100 * W / GCB_SOCEAN, 4), "U_M1": r5(Um / 100 * W / GCB_SOCEAN, 4),
                              "Delta": r5(abs(M1["Delta"]) / 100 * W / GCB_SOCEAN, 4)},
       "note": "U(s) 在 [0,1] 上近似线性；随 s 变化的分量＝U(s)−U(0)；常数项＝U(0)−Dep（湍流＋皮层常数项，均乘 1−c）"}
R["decomposition"] = dec

# ---------------- E10 thresholds
slope = U1 - U0
Iint0, Iint1 = float(np.interp(0, sg, Iw)), float(np.interp(1, sg, Iw))
R["thresholds"] = {"curve_slope_pp_per_s": r5(slope, 4), "GR_band_abs_s_minus_1_actual": r5(2.5 / slope, 3),
                   "GR_band_abs_s_minus_1_P30": r5(2.5 / 6.671, 3),
                   "SESOI_035_in_U_pp": r5(0.35 * slope, 3), "total_rain_slope_pp": r5(Iint1 - Iint0, 3),
                   "SESOI_035_in_total_rain_pp": r5(0.35 * (Iint1 - Iint0), 3)}

# ---------------- E11 latitudinal concentration
BD = C["bands"]
lo, hi = np.array(BD["lat_lo"]), np.array(BD["lat_hi"])
Ib = np.array(BD["whole_I_int_PgC"])  # [band][s]，正＝出海（吸收增量为负）
dI = np.array([np.interp(M1["s"], sg, Ib[j]) - np.interp(1.0, sg, Ib[j]) for j in range(len(lo))])
I1 = np.array([np.interp(1.0, sg, Ib[j]) for j in range(len(lo))])
lat = {}
for w in (10, 20):
    m = (lo >= -w) & (hi <= w)
    lat["le%d" % w] = {"frac_Delta_interfacial": r5(dI[m].sum() / dI.sum(), 3), "frac_I_s1": r5(I1[m].sum() / I1.sum(), 3)}
R["lat"] = lat

# ---------------- E12 min effect
bR = v_R["_boot"][np.isfinite(v_R["_boot"])]
R["min_effect"] = {"H0": "R_RIM(1 m) ≥ 0.65", "n_boot_ge_0.65": int(np.sum(bR >= 0.65)), "B": int(len(bR)),
                   "p_upper_bound": "<1e-4" if np.sum(bR >= 0.65) == 0 else r5(np.mean(bR >= 0.65))}

os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
json.dump(R, open(OUT, "w"), ensure_ascii=False, indent=1)
print("wrote", os.path.abspath(OUT))
