#!/usr/bin/env python3
"""p4f_wind_trend.py — 探索性分析：RIM-3 淡水不守恒倍数 M 能否解释 R_RIM(1 m) 随风速的下降。

背景：第一轮 D-U（p4-mech 阶段）已见 R_RIM(1 m) 按 U10 三分位 0.37→0.27→0.16；
P4c（p4c 阶段）给出 RIM-3 淡水倍数 M＝c1·√π·w/d0，d0 随风速增大，所以 M 随风速
下降。若「R<1 全由不守恒造成」，守恒缩放后的模型应与观测一致，逐事件预测 R_pred,i＝1/M_i，随风速**上升**——
与观测方向相反。本脚本把这一冲突定量化。

输入（只读已有产物，不读原始数据、不联网）：
  - p4-out/p4_mech_events.csv（p4-mech 阶段，sha256 c5feedaa…）：u10（0–6 h 模型 10 m 风均值）、obs1、rim1、c6、
    station、season、onset_utc。
  - p4b-out/p4b_events.csv（p4b 阶段，sha256 b02aa1e6…）：hfine100（CMORPH 驱动 RIM-3 在 0–100 m 细网格上
    −ΔF 的积分，m）、pcum_cmorph_model（模型时刻 CMORPH 累计雨量窗内均值，mm）。二者按每个事件的实际风速、
    d0 表与 24 h 雨史算出（P4b N12），是本脚本逐事件 M_i 的来源：M_i＝hfine100·1000/pcum_cmorph_model。
  - p4c_rim_conservation.py：只 import d0_table 与 C1，作解析近似的交叉核对。

映射：
  F1 逐事件 M_i 如上；580/646 事件可算（其余 hfine100 缺测或累计雨量为 0，与 P4b M_RIM 同集）。
  F2 「若 R<1 全由不守恒造成」＝守恒缩放后的 RIM-3（c1→c1/M_i，线性区等价于整个剖面除以 M_i）在 1 m 与观测相等：
     obs1_i ≈ rim1_i/M_i。于是三分位 k 的合并预测 R_pred,k＝Σ(rim1_i/M_i)/Σrim1_i（按模型 1 m 淡化加权的 1/M 均值），
     观测 R_obs,k＝Σobs1/Σrim1；残差比 ρ_k＝R_obs,k/R_pred,k＝Σobs1/Σ(rim1/M)。ρ≡1 即守恒解释在该档完全成立。
  F3 1 m 比值与整层积分并不完全对应，本映射依赖三条近似（结果切片须写明）：
     (a) 整层等比缩放＝「幅度整体放大、形状不错」（R1 读法）；若是 R2（海面对、剖面太深），1 m 的缩放不等于 1/M；
     (b) 线性区：低风大雨时 d0/(d0+x) 饱和，缩放 c1 与缩放 ΔS 不等价（P4c：最多差约 15%）；
     (c) M_i 是 0–6 h 窗内对时间平均的整层量，含当前项（c2，淡水不足）与历史项（c1，淡水过剩）的混合；1 m 的 rim1
         也是同窗均值，但两者在时间上的加权不同。
  F4 风速三分位：沿用第一轮切点（646 事件上的 1/3、2/3 分位 6.031179／9.325534 m s⁻¹），不在 580 子集上重算。
  F5 推断：站×季整簇 bootstrap（与 p2.cluster_boot 同算法：簇键排序、PCG64(20260926)、B=10000、Σnum/Σden、百分位），
     三档与对比量在同一次抽样中算；另报以站为簇的区间作对照。
  判读规则：
     方向：sign(R_pred,2−R_pred,0) 与 sign(R_obs,2−R_obs,0) 相反 → 「方向冲突」，相同 → 「方向一致」。
     幅度：ln(ρ_2/ρ_0) 的站×季 95% 区间不含 0 → 「守恒解释不能解释风速趋势」；含 0 → 「趋势差异不显著」。
     水平：三档 ρ_k 的区间都含 1 → 「守恒可解释各档水平」；否则逐档报告哪档在 1 之外。
  另报（描述）：580 事件不分档的合并 R_obs、R_pred、ρ；Spearman(u10, M_i)；各档合并 M_k＝Σhfine100·1000/Σpcum；解析近似 M_a＝C1·√π·0.95/d0(u10, c6/6)
  （d0 取 p4c 的胶囊表；w 取 0.95、雨强取 0–6 h CMORPH 均值，只作量级与方向核对）。
  V0：在 646 全集上按同算法复算第一轮 D-U 三档点值与区间，须与 p4_mech_summary.json 一致（|Δ|≤1e-5）。

用法（laptop，约 2 s）：python3 p4f_wind_trend.py <p4_mech_events.csv> <p4b_events.csv> [--p4sum p4_mech_summary.json] [--out p4f.json]
依赖：numpy、scipy（后者仅 p4c 的 d0 插值）。
Change Log：2026-09-27 初版（p4f-2026-09-27a）。
"""
import argparse
import csv
import hashlib
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

VERSION = "p4f-2026-09-27a"
B, SEED = 10000, 20260926
CUTS = (6.031179, 9.325534)            # 第一轮 D-U 切点（p4_mech_summary.json D_U.cuts）
W_TYP = 0.95                           # P4c：现场有效 w


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return math.nan
    return v


def load(p4_csv, p4b_csv):
    a = {(r["station"], r["onset_utc"]): r for r in csv.DictReader(open(p4_csv, encoding="utf-8"))}
    b = {(r["station"], r["onset_utc"]): r for r in csv.DictReader(open(p4b_csv, encoding="utf-8"))}
    assert set(a) == set(b), "两张事件表的事件集不一致"
    rows = []
    for k, ra in a.items():
        rb = b[k]
        o1, m1 = num(ra["obs1"]), num(ra["rim1"])
        assert abs(o1 - num(rb["obs1"])) < 1e-9 and abs(m1 - num(rb["rim1"])) < 1e-9, f"obs1/rim1 不一致 {k}"
        h, p = num(rb["hfine100"]), num(rb["pcum_cmorph_model"])
        M = h * 1000.0 / p if (np.isfinite(h) and np.isfinite(p) and p > 0) else math.nan
        rows.append({"station": ra["station"], "season": ra["season"], "u10": num(ra["u10"]), "obs1": o1,
                     "rim1": m1, "c6": num(ra["c6"]), "h": h, "p": p, "M": M})
    return rows


def cluster_draws(keys):
    uk = sorted(set(keys))
    kid = {k: j for j, k in enumerate(uk)}
    idx = np.array([kid[k] for k in keys], int)
    rng = np.random.Generator(np.random.PCG64(SEED))
    return idx, rng.integers(0, len(uk), size=(B, len(uk))), len(uk)


def boot_ratio(idx, draw, K, numv, denv):
    sn = np.bincount(idx, weights=np.asarray(numv, float), minlength=K)
    sd = np.bincount(idx, weights=np.asarray(denv, float), minlength=K)
    with np.errstate(divide="ignore", invalid="ignore"):
        return sn[draw].sum(axis=1) / sd[draw].sum(axis=1)


def ci(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return [round(float(np.percentile(a, 2.5)), 5), round(float(np.percentile(a, 97.5)), 5)]


def tercile_lab(u):
    return np.where(u <= CUTS[0], 0, np.where(u <= CUTS[1], 1, 2))


def block(rows, keyname):
    """三档 R_obs、R_pred、ρ、M_k 与对比量（同一次抽样）。keyname＝'ss'（站×季）或 'st'（站）。"""
    keys = [(r["station"], r["season"]) if keyname == "ss" else r["station"] for r in rows]
    idx, draw, K = cluster_draws(keys)
    lab = tercile_lab(np.array([r["u10"] for r in rows]))
    o = np.array([r["obs1"] for r in rows])
    m = np.array([r["rim1"] for r in rows])
    Mi = np.array([r["M"] for r in rows])
    h = np.array([r["h"] for r in rows]) * 1000.0
    p = np.array([r["p"] for r in rows])
    cells, bs = [], {}
    for k in range(3):
        s = lab == k
        z = lambda v: np.where(s, v, 0.0)  # noqa: E731
        bs[f"Robs{k}"] = boot_ratio(idx, draw, K, z(o), z(m))
        bs[f"Rpred{k}"] = boot_ratio(idx, draw, K, z(m / Mi), z(m))
        bs[f"rho{k}"] = boot_ratio(idx, draw, K, z(o), z(m / Mi))
        bs[f"M{k}"] = boot_ratio(idx, draw, K, z(h), z(p))
        sub = [r for j, r in enumerate(rows) if s[j]]
        cells.append({
            "k": k, "n": int(s.sum()), "G_station_season": len({(r["station"], r["season"]) for r in sub}),
            "G_station": len({r["station"] for r in sub}),
            "u10_median": round(float(np.median([r["u10"] for r in sub])), 3),
            "R_obs": round(float(o[s].sum() / m[s].sum()), 5), "R_obs_ci95": ci(bs[f"Robs{k}"]),
            "R_pred": round(float((m[s] / Mi[s]).sum() / m[s].sum()), 5), "R_pred_ci95": ci(bs[f"Rpred{k}"]),
            "rho": round(float(o[s].sum() / (m[s] / Mi[s]).sum()), 5), "rho_ci95": ci(bs[f"rho{k}"]),
            "M_pooled": round(float(h[s].sum() / p[s].sum()), 4), "M_pooled_ci95": ci(bs[f"M{k}"]),
            "M_median": round(float(np.median(Mi[s])), 4),
        })
    with np.errstate(divide="ignore", invalid="ignore"):
        lnrr = np.log(bs["rho2"] / bs["rho0"])
        dobs = bs["Robs2"] - bs["Robs0"]
        dpred = bs["Rpred2"] - bs["Rpred0"]
    contrast = {
        "R_obs_top_minus_bottom": round(cells[2]["R_obs"] - cells[0]["R_obs"], 5), "R_obs_diff_ci95": ci(dobs),
        "R_pred_top_minus_bottom": round(cells[2]["R_pred"] - cells[0]["R_pred"], 5), "R_pred_diff_ci95": ci(dpred),
        "ln_rho2_over_rho0": round(math.log(cells[2]["rho"] / cells[0]["rho"]), 5), "ln_rho2_over_rho0_ci95": ci(lnrr),
        "rho2_over_rho0": round(cells[2]["rho"] / cells[0]["rho"], 5),
        "frac_draws_dobs_lt0": round(float(np.mean(dobs < 0)), 4), "frac_draws_dpred_gt0": round(float(np.mean(dpred > 0)), 4),
    }
    return {"clusters": K, "cells": cells, "contrast": contrast}


def spearman(x, y):
    from scipy import stats
    r = stats.spearmanr(x, y)
    return round(float(r.statistic if hasattr(r, "statistic") else r.correlation), 4)


def judge(ss):
    c, t = ss["cells"], ss["contrast"]
    dir_obs = np.sign(c[2]["R_obs"] - c[0]["R_obs"])
    dir_pred = np.sign(c[2]["R_pred"] - c[0]["R_pred"])
    direction = "方向冲突" if dir_obs * dir_pred < 0 else "方向一致"
    lo, hi = t["ln_rho2_over_rho0_ci95"]
    magnitude = "守恒解释不能解释风速趋势" if (lo > 0 or hi < 0) else "趋势差异不显著"
    out1 = [cc["k"] for cc in c if not (cc["rho_ci95"][0] <= 1.0 <= cc["rho_ci95"][1])]
    level = "守恒可解释各档水平" if not out1 else f"ρ 的区间不含 1 的档：{out1}"
    return {"direction": direction, "magnitude": magnitude, "level": level}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("p4_events")
    ap.add_argument("p4b_events")
    ap.add_argument("--p4sum")
    ap.add_argument("--out")
    a = ap.parse_args()
    from p4c_rim_conservation import C1, d0_table
    D0 = d0_table()

    rows_all = load(a.p4_events, a.p4b_events)
    out = {"script": os.path.basename(__file__), "version": VERSION, "seed": SEED, "B": B, "cuts_u10": list(CUTS),
           "inputs": {"p4_mech_events.csv": sha(a.p4_events), "p4b_events.csv": sha(a.p4b_events),
                      "p4c_rim_conservation.py": sha(os.path.join(HERE, "p4c_rim_conservation.py")),
                      "self": sha(os.path.abspath(__file__))},
           "n_all": len(rows_all)}

    # V0：646 全集复算第一轮 D-U（只有 R_obs）
    idx, draw, K = cluster_draws([(r["station"], r["season"]) for r in rows_all])
    lab = tercile_lab(np.array([r["u10"] for r in rows_all]))
    o = np.array([r["obs1"] for r in rows_all])
    m = np.array([r["rim1"] for r in rows_all])
    v0 = []
    for k in range(3):
        s = lab == k
        b = boot_ratio(idx, draw, K, np.where(s, o, 0.0), np.where(s, m, 0.0))
        v0.append({"k": k, "n": int(s.sum()), "point": round(float(o[s].sum() / m[s].sum()), 5), "ci95": ci(b)})
    out["V0_full646"] = v0
    if a.p4sum:
        ref = json.load(open(a.p4sum))["D_U"]["cells"]
        dmax = max(max(abs(v0[k]["point"] - ref[k]["point"]), abs(v0[k]["ci95"][0] - ref[k]["ci95"][0]),
                       abs(v0[k]["ci95"][1] - ref[k]["ci95"][1])) for k in range(3))
        out["V0_max_abs_diff_vs_p4"] = dmax
        out["V0_pass"] = bool(dmax <= 1e-5)

    rows = [r for r in rows_all if np.isfinite(r["M"]) and np.isfinite(r["u10"])]
    out["n_M"] = len(rows)
    Mv = np.array([r["M"] for r in rows])
    out["M_pooled_all"] = round(float(sum(r["h"] for r in rows) * 1000 / sum(r["p"] for r in rows)), 5)
    out["M_event_quantiles"] = {q: round(float(np.quantile(Mv, q)), 4) for q in (0.05, 0.25, 0.5, 0.75, 0.95)}
    out["spearman_u10_M"] = spearman([r["u10"] for r in rows], Mv)
    # 全部 580 事件合并（不分档）：守恒解释能否说明 R 的平均水平
    idx, draw, K = cluster_draws([(r["station"], r["season"]) for r in rows])
    o = np.array([r["obs1"] for r in rows])
    m = np.array([r["rim1"] for r in rows])
    out["pooled_580"] = {
        "R_obs": round(float(o.sum() / m.sum()), 5), "R_obs_ci95": ci(boot_ratio(idx, draw, K, o, m)),
        "R_pred": round(float((m / Mv).sum() / m.sum()), 5), "R_pred_ci95": ci(boot_ratio(idx, draw, K, m / Mv, m)),
        "rho": round(float(o.sum() / (m / Mv).sum()), 5), "rho_ci95": ci(boot_ratio(idx, draw, K, o, m / Mv)),
        "clusters": K}
    ss = block(rows, "ss")
    st = block(rows, "st")
    out["station_season"] = ss
    out["station"] = st
    out["judgment"] = judge(ss)

    # 解析近似交叉核对（描述）
    ana = []
    for k in range(3):
        sub = [r for r in rows if tercile_lab(np.array([r["u10"]]))[0] == k]
        Ma = np.array([C1 * math.sqrt(math.pi) * W_TYP / D0(r["u10"], max(r["c6"] / 6.0, 0.0)) for r in sub])
        wts = np.array([r["rim1"] for r in sub])
        ana.append({"k": k, "M_analytic_median": round(float(np.nanmedian(Ma)), 3),
                    "R_pred_analytic": round(float(np.nansum(wts / Ma) / np.nansum(wts)), 4),
                    "d0_median": round(float(np.nanmedian([D0(r["u10"], max(r["c6"] / 6.0, 0.0)) for r in sub])), 3)})
    out["analytic_crosscheck"] = ana

    txt = json.dumps(out, ensure_ascii=False, indent=1)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            f.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
