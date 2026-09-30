#!/usr/bin/env python3
"""p4e_budget.py — P4e 探索性分析：E3′（每事件最深可用层的淡水收支闭合比）与 E5（R1「RIM-3 幅度整体放大」
vs R2「海面幅度对、剖面画得太深」的淡水守恒判别）。第一轮 negative_L 与 P4b「不可评」都不改；全部结果属探索性分析。

方法与判读规则在运行前写定；判读在 judge_e3p()／judge_e5()。

用法：
  python p4e_budget.py --selftest          # 合成数据自测（无网络、秒级）
  python p4e_budget.py                     # 正式（先跑自测，不过即退出 4）
  --p1-events / --p1b-dir / --p2-dir        同 p4_mech（默认 P1/P1b/P2 产物）
  --p4-dir DIR                              第一轮输出目录（p4_mech_events.csv 逐事件核对 V1）
  --p4b-dir DIR                             P4b 输出目录（p4b_summary.json，V-P4B 复现 18 个闭合块）
  --out DIR                                 输出目录（默认 $REPRO_OUTPUT_DIR）
依赖：numpy、scipy＋同目录 p1_events.py、p1b_extend.py、p2_rim_test.py、p4_mech.py、p4b_posthoc.py（只 import，不改）。
数据：只读现有缓存（P1/P1b 缓存、P2 主集 cmorph_pixels.csv、第一轮 CSV、P4b summary）；无新下载。
产物（<out>/）：p4e_selftest.json、p4e_summary.json（含判读）、p4e_events.csv（逐事件×窗派生量）、p4e_log.txt。
退出码：0 跑完（无论判读）；3 异常或 V1／V-P4B 核对不过；4 自测不过。

实现选择（Q 条）：
  Q1 事件：pm.build＋pm.event_mech 重建 646 事件；pb.check_v1 与第一轮 CSV 逐事件核对，不过 → p4e_abort.json、退出 3。
  Q2 逐事件量（event_e）：S0_z＝p2.pre_median（1 m 用 event_mech 的 s0_1，同 P4b）；观测 ΔS_z(W)＝p2.obs_delta，z∈{0.5,1,5,10,20,25}；
     雨量计驱动 RIM-3＝pb.gauge_forcing_ext（[H−30,H+48)）＋p2.rim_factor，小时值两半步平均，窗减雨前 6 h 中位（pm.delta_at）；
     CMORPH 驱动只前两窗；海面 f_R0＝−ΔF(z=0)；雨量分母＝P4b 同码（箱中点累计的窗内平均）。键名 _obs/_modg/_modc/_s0s/_rainmid
     与 P4b 相同，以便直接调用 pb.budget_rows／closure_block 做 V-P4B。
  Q3 E3′ 事件集＝pb.budget_rows(D=10,'gauge') 的准入＋雨量 >0；D_e：25 m 观测与模型有值且可积 → 25（有 20 m 用 25w20 权重），
     否则 20 m 可积 → 20，否则 10；权重 pb.W_D。Q＝Σh·1000/Σr̄；Q10 与 ΔQ 同一次抽样。
  Q4 E5：f_a＝−ΔS_a/S0_a（a＝0.5 m：sss05；a＝1 m：s1）；H＝g(ℓ,a)·(f_R0−f_a)，不截断；g(ℓ,a)＝[ℓ(1−e^{−a/ℓ})−a e^{−a/ℓ}]/(1−e^{−a/ℓ})，
     LIN＝a/2。E5a 的 Q_obs 用 0.5 m 起算权重 W05；E5b 用 E3′ 的 h。成员 𝓕＝LIN、EXP{0.20,0.10,0.05,0.02,0.01}；
     Q+H 的重抽由 Q_h、Q_Δ 两个比值数组线性组合（同一抽样，精确）。ℓ* 扫描 logspace(log10 0.005, log10 5, 61)＋LIN。
  Q5 推断：pm.boot_pairs（站×季、站两种簇，B=10000，seed 20260926）；逐站剔除 pm.jackknife_station；可评 n≥40 且 ≥8 簇。
  Q6 V-P4B：pb.budget_rows＋pb.closure_block 在本脚本逐事件量上重算 P4b 的 18 个闭合块，n 与点估计须与 p4b_summary.json 一致
     （|差|≤1e−4），不过 → p4e_abort.json、退出 3。V-sub（D_e=20/25 子集对 P4b D20/D25 块）只报。
  Q7 判读：judge_e3p、judge_e5（七类，按序取第一个成立者）、smooth 层（同分类只用 LIN）、ell_star。

Change Log：
  2026-09-27 初版。
"""

import argparse
import csv
import json
import math
import os
import sys
import time
import traceback
from datetime import datetime, timezone

import p1_events as p1
import p2_rim_test as p2
import p4_mech as pm
import p4b_posthoc as pb

VERSION = "p4e-2026-09-27a"
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
P4B_DIR_DEFAULT = _rp.upstream("p4b_dir")  # [repro] 读 p4b 阶段输出目录
NAN = float("nan")
isnum, rnd = pm.isnum, pm.rnd

# ---- 事先写定的常数 ----
WINDOWS = pb.WINDOWS                                   # [0,6) [6,12) [12,24) [24,48)
MAIN_WIN = (0, 6)
T_U = 1.30                                             # 雨量计漏接上限 k＝真/计 ≤1.30
T_SET = (1.0, 1.30, 1.5)
SHAPES = (("LIN", None), ("EXP0.20", 0.20), ("EXP0.10", 0.10), ("EXP0.05", 0.05), ("EXP0.02", 0.02), ("EXP0.01", 0.01))
ELL_SCAN_N = 61
ELL_SCAN_LO, ELL_SCAN_HI = 0.005, 5.0
W05 = {10: {0.5: 0.75, 1.0: 2.25, 5.0: 4.5, 10.0: 2.5},
       20: {0.5: 0.75, 1.0: 2.25, 5.0: 4.5, 10.0: 7.5, 20.0: 5.0},
       "25w20": {0.5: 0.75, 1.0: 2.25, 5.0: 4.5, 10.0: 7.5, 20.0: 7.5, 25.0: 2.5},
       "25n20": {0.5: 0.75, 1.0: 2.25, 5.0: 4.5, 10.0: 10.0, 25.0: 7.5}}
B_W5 = {1.0: 3.0, 5.0: 2.0}
B_W05_5 = {0.5: 0.75, 1.0: 2.25, 5.0: 2.0}
ZS = (0.5, 1.0, 5.0, 10.0, 20.0, 25.0)
SKEY = {0.5: "sss05", 1.0: "s1", 5.0: "s5", 10.0: "s10", 20.0: "s20", 25.0: "s25"}
VP4B_TOL = 1e-4


def _np():
    import numpy
    return numpy


def ell_scan():
    np = _np()
    return np.logspace(math.log10(ELL_SCAN_LO), math.log10(ELL_SCAN_HI), ELL_SCAN_N)


# ======================================================================== 形状与权重（Q4）
def g_shape(ell, a):
    """0–a m 内从 f_R0 过渡到 f_a 的剖面，超出 f_a 的面积系数：H＝g·(f_R0−f_a)。ell=None → 线性（a/2）。"""
    if ell is None or math.isinf(ell):
        return a / 2.0
    x = a / ell
    one_m = -math.expm1(-x)
    return (ell * one_m - a * math.exp(-x)) / one_m


def w05_for(D, has20):
    return W05["25w20" if has20 else "25n20"] if D == 25 else W05[D]


def deepest(obs, mod, s0s):
    """Q3：返回 (D, has20, h_obs_D)；1/5/10 m 都不可积时返回 (None, False, NaN)。"""
    has20 = isnum(obs.get(20.0)) and isnum(mod.get(20.0)) and isnum(s0s.get(20.0)) and s0s[20.0] > 0
    cands = []
    if isnum(obs.get(25.0)) and isnum(mod.get(25.0)):
        cands.append((25, has20))
    if has20:
        cands.append((20, True))
    cands.append((10, False))
    for D, h20 in cands:
        w = pb.weights_for(D, h20)
        ho, hm = pb.h_col(obs, s0s, w), pb.h_col(mod, s0s, w)
        if isnum(ho) and isnum(hm):
            return D, h20, ho
    return None, False, NAN


# ======================================================================== 逐事件量（Q2）
def rainmid_all(ser, i):
    """与 pb.event_extra 同码：从 t0 起 48 h 雨量计，各窗的箱中点累计平均。"""
    np = _np()
    out = {}
    rr = np.asarray(ser["rain"][i:i + pb.EXT_HI_H], float)
    if len(rr) == pb.EXT_HI_H and np.all(np.isfinite(rr)):
        rr = np.maximum(rr, 0.0)
        out = {w: pb.rain_mid_window(list(rr), w) for w in WINDOWS}
    elif len(rr):
        rr2 = np.where(np.isfinite(rr), np.maximum(rr, 0.0), np.nan)
        for w in WINDOWS:
            seg = rr2[:w[1]]
            if len(seg) == w[1] and np.all(np.isfinite(seg)):
                out[w] = pb.rain_mid_window(list(seg), w)
    return out


def event_e(stn, e, row, pixc):
    np = _np()
    ser, i, H = stn["ser"], e["i"], e["hour"]
    z1 = row["z1"]
    r = {k: row.get(k) for k in ("station", "group", "season", "onset_utc", "z1")}
    s0s = {z: (p2.pre_median(ser[SKEY[z]], i) if SKEY[z] in ser else NAN) for z in ZS}
    s0s[1.0] = row.get("s0_1")
    r["_s0s"] = s0s
    r["_obs"] = {w: {z: (p2.obs_delta(ser[SKEY[z]], i, w) if SKEY[z] in ser else NAN) for z in ZS} for w in WINDOWS}
    zmod = {1.0: z1, 5.0: 5.0, 10.0: 10.0, 20.0: 20.0, 25.0: 25.0}
    pre = list(range(6))

    def model(P, U, wins):
        hr = {z: p2.rim_factor(P, U, zmod[z]).reshape(-1, 2).mean(axis=1) for z in pb.DEEP_Z}
        h0 = p2.rim_factor(P, U, 0.0).reshape(-1, 2).mean(axis=1)
        mod, f0 = {}, {}
        for w in wins:
            post = list(range(6 + w[0], 6 + w[1]))
            mod[w] = {z: s0s[z] * pm.delta_at(hr[z], pre, post) if isnum(s0s[z]) else NAN for z in pb.DEEP_Z}
            d0 = pm.delta_at(h0, pre, post)
            f0[w] = -d0 if isnum(d0) else NAN
        return mod, f0

    r["_modc"], r["_fR0c"] = {}, {}
    P, U, _ = p2.event_forcing(stn, H, "cmorph", pixc)
    if P is not None:
        r["_modc"], r["_fR0c"] = model(P, U, WINDOWS[:2])
    r["_modg"], r["_fR0g"] = {}, {}
    Pe, Ue = pb.gauge_forcing_ext(stn, H)
    if Pe is not None:
        r["_modg"], r["_fR0g"] = model(Pe, Ue, WINDOWS)
    r["_rainmid"] = rainmid_all(ser, i)
    return r


# ======================================================================== 事件集（Q3、Q4）
def e3p_rows(rows, w):
    """E3′ 事件集（A 群）：pb.budget_rows(D=10,'gauge') 的准入＋雨量 >0；附最深可用层。"""
    out = []
    for r in rows:
        if r["group"] != p2.GROUP_A:
            continue
        mod = r["_modg"].get(w)
        if not mod:
            continue
        obs, s0s = r["_obs"][w], r["_s0s"]
        ho10, hm10 = pb.h_col(obs, s0s, pb.W_D[10]), pb.h_col(mod, s0s, pb.W_D[10])
        if not (isnum(ho10) and isnum(hm10) and isnum(obs.get(1.0)) and isnum(mod.get(1.0))):
            continue
        rain = r["_rainmid"].get(w, NAN)
        if not (isnum(rain) and rain > 0):
            continue
        D, h20, hD = deepest(obs, mod, s0s)
        out.append({"station": r["station"], "season": r["season"], "onset_utc": r["onset_utc"], "D": D, "has20": h20,
                    "h": hD, "h10": ho10, "rain": rain, "obs": obs, "s0s": s0s,
                    "fR0g": r["_fR0g"].get(w, NAN), "fR0c": r["_fR0c"].get(w, NAN)})
    return out


def e5_rows(X, a, forcing):
    """E5 行：a=0.5（E5a，0.5 m 起算）或 a=1.0（E5b）；forcing 'gauge'/'cmorph' 取 f_R0。"""
    out = []
    for x in X:
        fR0 = x["fR0g"] if forcing == "gauge" else x["fR0c"]
        if not isnum(fR0):
            continue
        if a == 0.5:
            oa, sa = x["obs"].get(0.5), x["s0s"].get(0.5)
            if not (isnum(oa) and isnum(sa) and sa > 0):
                continue
            fa = -oa / sa
            h = pb.h_col(x["obs"], x["s0s"], w05_for(x["D"], x["has20"]))
        else:
            fa = -x["obs"][1.0] / x["s0s"][1.0]
            h = x["h"]
        if not (isnum(h) and isnum(fa)):
            continue
        out.append({"station": x["station"], "season": x["season"], "h": h, "fa": fa, "fR0": fR0, "d": fR0 - fa,
                    "rain": x["rain"], "D": x["D"]})
    return out


# ======================================================================== 统计（Q5）
def _meta(X):
    return {"n": len(X), "G_station_season": len({(x["station"], x["season"]) for x in X}),
            "G_station": len({x["station"] for x in X}), "stations": sorted({x["station"] for x in X})}


def _evaluable(meta):
    return meta["n"] >= pm.BUDGET_MIN_N and meta["G_station_season"] >= pm.BUDGET_MIN_G


def ratio_family(X, nums, den):
    """nums: {名: 每行分子函数}；den: 每行分母函数。返回 {名: 点、站×季 CI、站簇 CI、逐站剔除}＋两种簇的重抽数组。"""
    np = _np()
    res = {}
    if len(X) < 5:
        return res, None, None
    dv = [den(x) for x in X]
    pairs = {k: ([f(x) for x in X], dv) for k, f in nums.items()}
    bss, _ = pm.boot_pairs([(x["station"], x["season"]) for x in X], pairs)
    bst, _ = pm.boot_pairs([x["station"] for x in X], pairs)
    for k, f in nums.items():
        nv = pairs[k][0]
        pt = sum(nv) / sum(dv) if sum(dv) else NAN
        res[k] = {"point": rnd(pt), "ci95_ss": pm.ci95(bss[k])[0], "ci95_station": pm.ci95(bst[k])[0],
                  "jackknife_station": pm.jackknife_station(
                      [dict(x, _n=f(x), _d=den(x)) for x in X], lambda rr: pm.ratio_of(rr, "_n", "_d"))}
    return res, bss, bst


def e3p_block(X):
    np = _np()
    b = _meta(X)
    b["evaluable"] = _evaluable(b)
    b["D_counts"] = {str(D): sum(1 for x in X if x["D"] == D) for D in (25, 20, 10)}
    if len(X) < 5:
        return b
    fam, bss, bst = ratio_family(X, {"Q": lambda x: x["h"] * 1000, "Q10": lambda x: x["h10"] * 1000},
                                 lambda x: x["rain"])
    b.update(fam)
    b["dQ_deep_minus_10"] = {"point": rnd(fam["Q"]["point"] - fam["Q10"]["point"]),
                             "ci95_ss": pm.ci95(bss["Q"] - bss["Q10"])[0], "ci95_station": pm.ci95(bst["Q"] - bst["Q10"])[0]}
    # 各层段贡献（描述）：Σ w_z f_z ·1000/Σr̄ 按层
    seg = {}
    for x in X:
        w = pb.weights_for(x["D"], x["has20"])
        for z, wt in w.items():
            seg[z] = seg.get(z, 0.0) + wt * (-x["obs"][z] / x["s0s"][z])
    den = sum(x["rain"] for x in X)
    b["layer_contrib"] = {f"{z:g}m": rnd(v * 1000 / den) for z, v in sorted(seg.items())}
    b["reading"] = judge_e3p(b)
    return b


def e5_block(Y, a):
    """E5：Q_obs 与 𝓕 各成员 Q+H、ℓ* 扫描、描述量；同一次抽样。"""
    np = _np()
    b = _meta(Y)
    b["a_m"] = a
    b["evaluable"] = _evaluable(b)
    if len(Y) < 5:
        return b
    fam, bss, bst = ratio_family(Y, {"Qh": lambda x: x["h"] * 1000, "Qd": lambda x: x["d"] * 1000},
                                 lambda x: x["rain"])
    qh, qd = fam["Qh"]["point"], fam["Qd"]["point"]
    b["Q_obs"] = fam["Qh"]
    b["Q_delta_f_per_rain"] = fam["Qd"]
    sfa = sum(x["fa"] for x in Y)
    b["surface_over_a_ratio"] = rnd(sum(x["fR0"] for x in Y) / sfa) if sfa else None
    mem = {}
    for name, ell in SHAPES:
        g = g_shape(ell, a)
        ss, st = bss["Qh"] + g * bss["Qd"], bst["Qh"] + g * bst["Qd"]
        mem[name] = {"ell_m": ell, "g_m": rnd(g, 6), "H_over_rain": rnd(g * qd), "point": rnd(qh + g * qd),
                     "ci95_ss": pm.ci95(ss)[0], "ci95_station": pm.ci95(st)[0],
                     "jackknife_station": pm.jackknife_station(
                         [dict(x, _n=(x["h"] + g * x["d"]) * 1000, _d=x["rain"]) for x in Y],
                         lambda rr: pm.ratio_of(rr, "_n", "_d"))}
    b["members"] = mem
    fav = min(mem, key=lambda k: (mem[k]["ci95_ss"][0] if mem[k]["ci95_ss"][0] is not None else math.inf))
    b["most_favorable_to_R2"] = fav
    b["H_needed_to_reach_T_U_over_rain"] = rnd(T_U - qh)
    # ℓ* 扫描
    scan = []
    for ell in ell_scan():
        g = g_shape(float(ell), a)
        scan.append({"ell_m": rnd(float(ell), 5), "g_m": rnd(g, 6), "point": rnd(qh + g * qd),
                     "ci_lo_ss": pm.ci95(bss["Qh"] + g * bss["Qd"])[0][0]})
    b["scan"] = scan
    b["ell_star"] = {str(T): ell_star(scan, mem["LIN"]["ci95_ss"][0], T) for T in T_SET}
    return b


def ell_star(scan, lin_lo, T):
    """最小网格 ℓ，使其自身及所有更大的 ℓ′ 与 LIN 的站×季 CI 下端都 >T；不存在返回 None。"""
    if lin_lo is None or lin_lo <= T:
        return None
    k_star = None
    for k in range(len(scan) - 1, -1, -1):
        lo = scan[k]["ci_lo_ss"]
        if lo is not None and lo > T:
            k_star = k
        else:
            break
    return scan[k_star]["ell_m"] if k_star is not None else "≥LIN（网格上无）"


def b_block(rows, w):
    """B 群描述：0–5 m Q5、0.5 m 起算 Q5 与 Q5+H（LIN、EXP0.01），雨量计驱动 f_R0。"""
    X, Y = [], []
    for r in rows:
        if r["group"] != p2.GROUP_B:
            continue
        obs, s0s = r["_obs"][w], r["_s0s"]
        rain = r["_rainmid"].get(w, NAN)
        if not (isnum(rain) and rain > 0):
            continue
        h5 = pb.h_col(obs, s0s, B_W5)
        if isnum(h5):
            X.append({"station": r["station"], "season": r["season"], "h": h5, "rain": rain})
        h05 = pb.h_col(obs, s0s, B_W05_5)
        fR0 = r["_fR0g"].get(w, NAN)
        if isnum(h05) and isnum(fR0):
            fa = -obs[0.5] / s0s[0.5]
            Y.append({"station": r["station"], "season": r["season"], "h": h05, "d": fR0 - fa, "rain": rain})
    out = {"Q5": _meta(X), "Q5_from05": _meta(Y)}
    if len(X) >= 5:
        fam, _, _ = ratio_family(X, {"Q5": lambda x: x["h"] * 1000}, lambda x: x["rain"])
        out["Q5"].update(fam["Q5"])
    if len(Y) >= 5:
        fam, bss, _ = ratio_family(Y, {"Qh": lambda x: x["h"] * 1000, "Qd": lambda x: x["d"] * 1000}, lambda x: x["rain"])
        out["Q5_from05"]["Q_obs"] = fam["Qh"]
        for name, ell in (("LIN", None), ("EXP0.01", 0.01)):
            g = g_shape(ell, 0.5)
            out["Q5_from05"][f"QH_{name}"] = {"point": rnd(fam["Qh"]["point"] + g * fam["Qd"]["point"]),
                                              "ci95_ss": pm.ci95(bss["Qh"] + g * bss["Qd"])[0]}
    return out


# ======================================================================== 判读（Q7）
def _cls(lo, hi, T=T_U):
    if lo is None or hi is None:
        return None
    if lo <= 1 <= hi:
        return "闭合"
    if hi < 1:
        return "缺失"
    return "多出且超过漏接上限" if lo > T else "多出"


def judge_e3p(b):
    if not b.get("evaluable"):
        return {"status": "not_evaluable"}
    q = b.get("Q", {})
    c_ss = _cls(*q.get("ci95_ss", [None, None]))
    c_st = _cls(*q.get("ci95_station", [None, None]))
    return {"status": c_ss, "robust": "稳健" if c_ss == c_st else "依赖簇定义", "station_cluster_class": c_st}


def judge_e5(b, only=None, T=T_U):
    """E5 七类（按序取第一个成立者）；only='LIN' 时为平滑层（成员只用 LIN）。"""
    if not b or not b.get("evaluable"):
        return {"class": 1, "status": "不可评"}
    q_lo, q_hi = b["Q_obs"]["ci95_ss"]
    mem = b["members"] if only is None else {only: b["members"][only]}
    lows = {k: v["ci95_ss"][0] for k, v in mem.items()}
    fav = min(lows, key=lambda k: lows[k] if lows[k] is not None else math.inf)
    L_min, P_fav = lows[fav], mem[fav]["point"]
    base = {"member": fav, "L_min": L_min, "P_fav": P_fav, "Q_obs_ci95_ss": [q_lo, q_hi], "T": T}
    if q_lo is None or q_hi is None or L_min is None or P_fav is None:
        return dict(base, **{"class": 1, "status": "不可评"})
    if q_lo > T:
        return dict(base, **{"class": 2, "status": "不可区分（观测自身超出）"})
    if L_min > T:
        return dict(base, **{"class": 3, "status": "R2 与淡水守恒不相容，支持 R1"})
    if q_lo <= 1 <= q_hi and P_fav <= T:
        return dict(base, **{"class": 4, "status": "不能区分"})
    if P_fav > T:
        return dict(base, **{"class": 5, "status": "倾向 R1（不显著）"})
    if q_hi < 1:
        return dict(base, **{"class": 6, "status": "不能区分（观测缺淡水）"})
    return dict(base, **{"class": 7, "status": "不能区分（观测略多，在漏接内）"})


def e5_wording(main, smooth, ell_star_T):
    if main.get("class") == 3:
        return "事后诊断：在最有利于 R2 的剖面假设下 R2 仍违背淡水守恒，支持 R1（幅度整体放大）"
    if smooth.get("class") == 3:
        return (f"R2 只有在海面淡水集中于 e-folding <ℓ*={ell_star_T} m 的薄皮层时才与淡水守恒相容；"
                "若 0–0.5 m 剖面平滑（线性），R2 与守恒不相容（条件句，事后诊断）")
    return "淡水收支不能区分 R1 与 R2"


# ======================================================================== 分析
def analyze(rows):
    S = {"E3p": {}, "E5": {}, "B": {}}
    e5 = {"E5a_gauge": {}, "E5b_gauge": {}, "E5a_cmorph": {}, "E5b_cmorph": {}}
    for w in WINDOWS:
        tag = f"{w[0]}-{w[1]}h"
        X = e3p_rows(rows, w)
        S["E3p"][tag] = e3p_block(X)
        e5["E5a_gauge"][tag] = e5_block(e5_rows(X, 0.5, "gauge"), 0.5)
        e5["E5b_gauge"][tag] = e5_block(e5_rows(X, 1.0, "gauge"), 1.0)
        if w in WINDOWS[:2]:
            e5["E5a_cmorph"][tag] = e5_block(e5_rows(X, 0.5, "cmorph"), 0.5)
            e5["E5b_cmorph"][tag] = e5_block(e5_rows(X, 1.0, "cmorph"), 1.0)
        S["B"][tag] = b_block(rows, w)
    S["E5"]["blocks"] = e5
    # 主判
    mt = f"{MAIN_WIN[0]}-{MAIN_WIN[1]}h"
    main_ver = "E5a_gauge" if e5["E5a_gauge"][mt].get("evaluable") else "E5b_gauge"
    mb = e5[main_ver][mt]
    jm, js = judge_e5(mb), judge_e5(mb, only="LIN")
    es = mb.get("ell_star", {}).get(str(T_U))
    S["E5"]["main"] = {"window": mt, "version": main_ver, "judgment": jm, "smooth_judgment": js, "ell_star_T_U": es,
                       "wording": e5_wording(jm, js, es)}
    sens = {}
    for ver, blocks in e5.items():
        for tag, blk in blocks.items():
            for T in T_SET:
                sens[f"{ver}|{tag}|T={T}"] = {"main_class": judge_e5(blk, T=T).get("status"),
                                               "smooth_class": judge_e5(blk, only="LIN", T=T).get("status")}
    S["E5"]["sensitivity_classes"] = sens
    S["E5"]["sensitivity_disagreeing_with_main"] = sorted(k for k, v in sens.items() if k.endswith(f"T={T_U}")
                                                          and v["main_class"] != jm.get("status"))
    return S


def check_p4b(rows, p4b_summary):
    """Q6 V-P4B：重算 P4b 的 18 个闭合块。"""
    ref = p4b_summary["E3"]["blocks"]
    res, ok = {}, True
    for forcing, wins in (("cmorph", WINDOWS[:2]), ("gauge", WINDOWS)):
        for w in wins:
            for D in (10, 20, 25):
                key = f"{forcing}_{w[0]}-{w[1]}h_D{D}"
                cb = pb.closure_block(pb.budget_rows(rows, D, forcing, w))
                rc = ref.get(key, {}).get("closure_obs_over_gauge", {})
                good = cb.get("n") == rc.get("n") and isnum(cb.get("point")) and isnum(rc.get("point")) \
                    and abs(cb["point"] - rc["point"]) <= VP4B_TOL
                res[key] = {"this": [cb.get("n"), cb.get("point")], "p4b": [rc.get("n"), rc.get("point")], "match": bool(good)}
                ok = ok and good
    return {"pass": ok, "blocks": res}


def v_sub(rows, p4b_summary):
    ref = p4b_summary["E3"]["blocks"]
    out = {}
    for w in WINDOWS:
        X = e3p_rows(rows, w)
        for D in (20, 25):
            sub = [x for x in X if x["D"] == D]
            pt = rnd(sum(x["h"] for x in sub) * 1000 / sum(x["rain"] for x in sub)) if sub else None
            rc = ref.get(f"gauge_{w[0]}-{w[1]}h_D{D}", {}).get("closure_obs_over_gauge", {})
            out[f"{w[0]}-{w[1]}h_D{D}"] = {"this": [len(sub), pt], "p4b": [rc.get("n"), rc.get("point")]}
    return out


def run(args, out_dir, log):
    t0 = time.monotonic()
    stations, events, val, prim, plan = pm.build(args, out_dir, log)
    data, cov = pm.load_pixels(args.p2_dir, prim)
    sbn = {s["name"]: s for s in stations}
    pixc, pixm = pm.PixView(data, "c"), pm.PixView(data, "m")
    rows = [pm.event_mech(sbn[e["_stn"]], e, pixc, pixm)[0] for e in events]
    v1 = pb.check_v1(rows, os.path.join(args.p4_dir, "p4_mech_events.csv"))
    log.log(f"V1 {json.dumps(v1, ensure_ascii=False)}", echo=True)
    if not v1["pass"]:
        p2.jdump({"V1": v1}, os.path.join(out_dir, "p4e_abort.json"))
        raise RuntimeError("V1 不过")
    ext = [event_e(sbn[e["_stn"]], e, rows[k], pixc) for k, e in enumerate(events)]
    log.log(f"event_e 完成 {len(ext)} 事件 {time.monotonic() - t0:.0f}s", echo=True)
    s4b = json.load(open(os.path.join(args.p4b_dir, "p4b_summary.json"), encoding="utf-8"))
    vp = check_p4b(ext, s4b)
    log.log(f"V-P4B pass={vp['pass']}", echo=True)
    if not vp["pass"]:
        p2.jdump({"V_P4B": vp}, os.path.join(out_dir, "p4e_abort.json"))
        raise RuntimeError("V-P4B 不过")
    S = analyze(ext)
    out = {"script": "p4e_budget.py", "version": VERSION,
           "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "n_events": len(events),
           "rebuild_validation": {k: val[k] for k in ("n_ref", "n_rebuilt", "d4_count")}, "cmorph_pixel_cache": cov,
           "downloads": {"cache_http_requests": plan.get("cache_http_requests"), "new_data_bytes_expected": 0},
           "constants": {"T_U": T_U, "T_SET": T_SET, "shapes": SHAPES, "main_window": MAIN_WIN,
                         "ell_scan": [ELL_SCAN_LO, ELL_SCAN_HI, ELL_SCAN_N], "B": pm.B, "seed": pm.SEED},
           "V1_round1_csv": v1, "V_P4B": vp, "V_sub_info": v_sub(ext, s4b)}
    out.update(S)
    out["runtime_s"] = round(time.monotonic() - t0, 1)
    p2.jdump(out, os.path.join(out_dir, "p4e_summary.json"))
    cols = ["station", "group", "season", "onset_utc", "window", "D", "has20", "h10", "hD", "h05D", "rain", "f05", "f1",
            "fR0g", "fR0c"]
    with open(os.path.join(out_dir, "p4e_events.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=cols)
        wr.writeheader()
        for w in WINDOWS:
            for x in e3p_rows(ext, w):
                o5, s5 = x["obs"].get(0.5), x["s0s"].get(0.5)
                f05 = -o5 / s5 if isnum(o5) and isnum(s5) and s5 > 0 else NAN
                h05 = pb.h_col(x["obs"], x["s0s"], w05_for(x["D"], x["has20"]))
                vals = {"station": x["station"], "group": p2.GROUP_A, "season": x["season"], "onset_utc": x["onset_utc"],
                        "window": f"{w[0]}-{w[1]}h", "D": x["D"], "has20": int(x["has20"]), "h10": x["h10"], "hD": x["h"],
                        "h05D": h05, "rain": x["rain"], "f05": f05, "f1": -x["obs"][1.0] / x["s0s"][1.0],
                        "fR0g": x["fR0g"], "fR0c": x["fR0c"]}
                wr.writerow({k: (rnd(v, 7) if isinstance(v, float) else v) for k, v in vals.items()})
    print(json.dumps({"E3p": {k: v.get("reading") for k, v in S["E3p"].items()}, "E5_main": S["E5"]["main"]},
                     ensure_ascii=False), flush=True)
    return 0


# ======================================================================== 自测（合成数据，无网络）
def selftest(out_dir=None):
    np = _np()
    from scipy import integrate
    res, fails = [], []
    t_start = time.monotonic()

    def check(cond, msg):
        res.append(("PASS" if cond else "FAIL") + " " + msg)
        if not cond:
            fails.append(msg)

    # 1 g(ℓ,a) 对数值积分与极限
    errs = []
    for a in (0.5, 1.0):
        for ell in (0.01, 0.05, 0.2, 1.0, 5.0):
            prof = lambda z: (math.exp(-z / ell) - math.exp(-a / ell)) / (1 - math.exp(-a / ell))
            num = integrate.quad(prof, 0, a, epsabs=1e-13, epsrel=1e-11, points=[min(a, 5 * ell)])[0]
            errs.append(abs(g_shape(ell, a) - num) / num)
    check(max(errs) < 1e-8, f"g(ℓ,a) 与数值积分一致（最大相对误差 {max(errs):.1e}）")
    check(abs(g_shape(1e4, 0.5) - 0.25) < 1e-5 and g_shape(None, 0.5) == 0.25, "g：ℓ→∞ → a/2（LIN）")
    check(abs(g_shape(1e-3, 0.5) / 1e-3 - 1) < 1e-9, "g：ℓ→0 → ℓ")
    gs = [g_shape(e, 0.5) for e in (0.01, 0.02, 0.05, 0.1, 0.2, 1.0, None)]
    check(all(gs[k] < gs[k + 1] for k in range(len(gs) - 1)), "g 随 ℓ 单调增，LIN 最大")
    # 2 权重和
    check(all(abs(sum(w.values()) - d) < 1e-12 for w, d in ((W05[10], 10), (W05[20], 20), (W05["25w20"], 25),
                                                           (W05["25n20"], 25), (B_W5, 5), (B_W05_5, 5))),
          "E5a／B 群权重和＝积分深度")
    # 3 均匀剖面：Q＝1（两种顶部口径、三种 D）
    s0s = {z: 35.0 for z in ZS}
    ok = True
    for D, h20, wts in ((10, False, pb.W_D[10]), (20, True, pb.W_D[20]), (25, True, pb.W_D["25w20"]),
                        (25, False, pb.W_D["25n20"])):
        c = 0.002
        obs = {z: -c * 35.0 for z in ZS}
        rain = c * D * 1000
        ok = ok and abs(pb.h_col(obs, s0s, wts) * 1000 / rain - 1) < 1e-12
        ok = ok and abs(pb.h_col(obs, s0s, w05_for(D, h20)) * 1000 / rain - 1) < 1e-12
    check(ok, "均匀淡水分数剖面 → Q＝1（0–1 m 与 0–0.5 m 起算，D=10/20/25）")
    # 4 已知 R2 剖面：0–0.5 m 指数（ℓ=0.1）、以下均匀 → Q_obs(W05)＋H(EXP0.10)＝真积分
    f0, fa, D = 0.02, 0.004, 10
    true_int = fa * D + (f0 - fa) * g_shape(0.1, 0.5)
    obs = {z: -fa * 35.0 for z in ZS}
    hq = pb.h_col(obs, s0s, W05[10])
    check(abs(hq + g_shape(0.1, 0.5) * (f0 - fa) - true_int) < 1e-14, "已知剖面：Q_obs＋H(真形状)＝真积分")
    # 数值核：直接积分合成剖面
    num = integrate.quad(lambda z: fa + (f0 - fa) * (math.exp(-z / 0.1) - math.exp(-5)) / (1 - math.exp(-5)), 0, 0.5)[0] + fa * 9.5
    check(abs(num - true_int) / true_int < 1e-9, "已知剖面：解析 H 与剖面数值积分一致")
    # 5 最深层选择
    mod = {z: -0.1 for z in pb.DEEP_Z}
    o_all = {z: -0.05 for z in ZS}
    check(deepest(o_all, mod, s0s)[:2] == (25, True), "最深层：25 与 20 都有 → D=25（25w20）")
    o_no20 = {**o_all, 20.0: NAN}
    check(deepest(o_no20, mod, s0s)[:2] == (25, False), "最深层：只有 25 → D=25（25n20）")
    o_no25 = {**o_all, 25.0: NAN}
    check(deepest(o_no25, mod, s0s)[:2] == (20, True), "最深层：只有 20 → D=20")
    o_10 = {**o_all, 25.0: NAN, 20.0: NAN}
    check(deepest(o_10, mod, s0s)[:2] == (10, False), "最深层：无深层 → D=10")
    check(deepest({**o_10, 5.0: NAN}, mod, s0s)[0] is None, "最深层：5 m 缺 → 不可积")
    # 6 E5 块：簇完全相同 → CI 退化为点；H 在 f_R0＝f_a 时为 0；ℓ* 与解析一致
    Y = []
    for st in range(5):
        for se in range(3):
            for _ in range(4):
                Y.append({"station": f"S{st}", "season": se, "h": 0.012, "fa": 0.003, "fR0": 0.003, "d": 0.0,
                          "rain": 10.0, "D": 10})
    b0 = e5_block(Y, 0.5)
    check(b0["evaluable"] and all(abs(m["point"] - 1.2) < 1e-9 and abs(m["H_over_rain"]) < 1e-12
                                  for m in b0["members"].values()), "E5：f_R0＝f_a → H＝0，Q+H＝Q_obs＝1.2")
    check(all(m["ci95_ss"] == [1.2, 1.2] for m in b0["members"].values()), "E5：簇完全相同 → CI 退化为点")
    Y2 = [dict(y, fR0=0.043, d=0.04) for y in Y]          # Qd＝0.04·1000/10＝4 → Q+H＝1.2＋4g
    b2 = e5_block(Y2, 0.5)
    g_need = (T_U - 1.2) / 4.0                             # 1.2＋4g＞1.3 ⇔ g＞0.025
    grid = [s["ell_m"] for s in b2["scan"]]
    exp_star = next(e for e in grid if g_shape(e, 0.5) > g_need + 1e-12)
    check(b2["ell_star"][str(T_U)] == exp_star, f"ℓ* 扫描＝解析（g>{g_need:.4f} 的最小网格 ℓ＝{exp_star}，得 {b2['ell_star'][str(T_U)]}）")
    check(b2["most_favorable_to_R2"] == "EXP0.01", "Σ(f_R0−f_a)>0 → 最有利成员＝EXP0.01")
    b3 = e5_block([dict(y, fR0=0.004, d=0.001) for y in Y], 0.5)   # Qd＝0.1 → LIN 只加 0.025
    check(b3["ell_star"][str(T_U)] is None and b3["members"]["LIN"]["point"] < T_U, "ℓ*：连 LIN 都不超 T → None")
    # 7 判读七类与平滑层
    def fake(qlo, qhi, members):
        return {"evaluable": True, "Q_obs": {"ci95_ss": [qlo, qhi]},
                "members": {k: {"ci95_ss": [lo, lo + 0.5], "point": pt} for k, (lo, pt) in members.items()}}
    lin_hi = ("LIN", (1.45, 1.8))
    cases = [
        (fake(1.4, 2.0, dict([lin_hi, ("EXP0.01", (1.4, 1.7))])), 2),
        (fake(1.1, 1.8, dict([lin_hi, ("EXP0.01", (1.35, 1.6))])), 3),
        (fake(0.7, 1.7, dict([lin_hi, ("EXP0.01", (0.71, 1.2))])), 4),
        (fake(0.9, 1.9, dict([lin_hi, ("EXP0.01", (0.95, 1.35))])), 5),
        (fake(0.3, 0.8, dict([("LIN", (0.5, 1.0)), ("EXP0.01", (0.31, 0.55))])), 6),
        (fake(1.05, 1.6, dict([("LIN", (1.1, 1.3)), ("EXP0.01", (1.06, 1.25))])), 7),
    ]
    got = [judge_e5(b)["class"] for b, _ in cases]
    check(got == [c for _, c in cases], f"E5 判读七类分支（得 {got}）")
    check(judge_e5({"evaluable": False})["class"] == 1, "E5 判读：不可评")
    sm = judge_e5(cases[2][0], only="LIN")
    check(sm["class"] == 3 and "薄皮层" in e5_wording(judge_e5(cases[2][0]), sm, 0.12), "平滑层：LIN 下 3 类 → 条件句措辞")
    check(e5_wording({"class": 4}, {"class": 4}, None) == "淡水收支不能区分 R1 与 R2", "两层都非 3 → 不能区分")
    check(judge_e3p({"evaluable": True, "Q": {"ci95_ss": [0.7, 1.6], "ci95_station": [0.8, 1.5]}})["status"] == "闭合"
          and judge_e3p({"evaluable": True, "Q": {"ci95_ss": [1.35, 2.0], "ci95_station": [0.9, 2.1]}})["robust"] == "依赖簇定义"
          and judge_e3p({"evaluable": True, "Q": {"ci95_ss": [1.35, 2.0], "ci95_station": [1.4, 2.1]}})["status"] == "多出且超过漏接上限"
          and judge_e3p({"evaluable": True, "Q": {"ci95_ss": [0.1, 0.6], "ci95_station": [0.1, 0.6]}})["status"] == "缺失",
          "E3′ 读法分支（闭合／依赖簇定义／多出且超过漏接上限／缺失）")
    # 8 端到端：p2 两站夹具（观测＝0.5×RIM 各深度）→ event_mech → event_e → V-P4B 自洽 → analyze
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        fx = p2._e2e_fixture(td)
        pixd = {h: {nm: (v[0], v[1], v[0], v[1], 9, 9) for nm, v in d.items()} for h, d in fx["cm"].items()}
        pixc, pixm = pm.PixView(pixd, "c"), pm.PixView(pixd, "m")
        for stn in fx["stations"]:
            n_ = stn["n"]
            Pfull = np.zeros(2 * n_)
            for hh in range(n_):
                Pfull[2 * hh:2 * hh + 2] = pixd[stn["h0"] + hh][stn["name"]][:2]
            Uh = np.repeat(stn["ser"]["wind"] * p2.WIND_FACTOR, 2)
            for key, zt in (("s10", 10.0), ("s20", 20.0), ("s25", 25.0)):
                if stn["group"] == p2.GROUP_A:
                    F = np.ones(2 * n_)
                    F[48:] = p2.rim_factor(Pfull, Uh, zt)
                    stn["ser"][key] = 35.0 * (1 - 0.5 * (1 - F.reshape(-1, 2).mean(axis=1)))
                else:
                    stn["ser"][key] = np.full(n_, np.nan)
        sbn = {s["name"]: s for s in fx["stations"]}
        evs = fx["ev"][p2.PRIMARY_THRESHOLD]
        rows_e = [pm.event_mech(sbn[e["_stn"]], e, pixc, pixm)[0] for e in evs]
        ext = [event_e(sbn[e["_stn"]], e, rows_e[k], pixc) for k, e in enumerate(evs)]
        X = e3p_rows(ext, MAIN_WIN)
        check(len(X) > 10 and all(x["D"] == 25 and x["has20"] for x in X), f"端到端 E3′：A 站全部取 D=25（n={len(X)}）")
        cb10 = pb.closure_block(pb.budget_rows(ext, 10, "gauge", MAIN_WIN))
        q10 = sum(x["h10"] for x in X) * 1000 / sum(x["rain"] for x in X)
        check(cb10["n"] == len(X) and abs(cb10["point"] - rnd(q10)) < 1e-9, "端到端 V-P4B 自洽：pb.closure_block(D10)＝本脚本 Q10")
        cb25 = pb.closure_block(pb.budget_rows(ext, 25, "gauge", MAIN_WIN))
        q25 = sum(x["h"] for x in X) * 1000 / sum(x["rain"] for x in X)
        check(abs(cb25["point"] - rnd(q25)) < 1e-9, "端到端：D=25 子集＝pb D25 块")
        # 观测＝0.5×模型（同权重） → Q_obs＝0.5×Σh_mod/Σr̄
        hm = sum(pb.h_col(r["_modg"][MAIN_WIN], r["_s0s"], pb.W_D["25w20"]) for r in ext
                 if r["group"] == p2.GROUP_A and r["_modg"].get(MAIN_WIN) and any(x["onset_utc"] == r["onset_utc"] for x in X))
        check(abs(q25 / (0.5 * hm * 1000 / sum(x["rain"] for x in X)) - 1) < 0.02, "端到端：Q_obs ≈ 0.5×模型淡水/雨量（夹具真比 0.5）")
        e0 = next(r for r in ext if r["group"] == p2.GROUP_A and r["_modg"])
        stn0 = sbn[e0["station"]]
        ev0 = next(e for e in evs if e["onset_utc"] == e0["onset_utc"] and e["_stn"] == e0["station"])
        Pe, Ue = pb.gauge_forcing_ext(stn0, ev0["hour"])
        d0 = pm.delta_at(p2.rim_factor(Pe, Ue, 0.0).reshape(-1, 2).mean(axis=1), list(range(6)), list(range(6, 12)))
        check(abs(e0["_fR0g"][MAIN_WIN] + d0) < 1e-15 and e0["_fR0g"][MAIN_WIN] > 0, "端到端：f_R0＝−ΔF(0)（手算）且 >0")
        check(all(abs(r["_fR0g"][w] - r["_fR0c"][w]) < 1e-12 for r in ext if r["_modg"] and r["_modc"] for w in WINDOWS[:2]),
              "端到端：夹具雨量计＝CMORPH → 两种驱动 f_R0 相同")
        Y = e5_rows(X, 0.5, "gauge")
        y0 = Y[0]
        x0 = X[0]
        fa_hand = -x0["obs"][0.5] / x0["s0s"][0.5]
        check(len(Y) > 5 and abs(y0["fa"] - fa_hand) < 1e-15 and abs(y0["d"] - (x0["fR0g"] - fa_hand)) < 1e-15,
              f"端到端 E5a 行：f_a、Δf 手算一致（n={len(Y)}）")
        try:
            Se = analyze(ext)
            err = None
        except Exception as ex:
            Se, err = None, repr(ex) + traceback.format_exc()[-600:]
        check(err is None, f"端到端 analyze 无异常（{err}）")
        if Se:
            mb = Se["E5"]["blocks"]["E5a_gauge"]["0-6h"]
            lin = mb["members"]["LIN"]
            hand = (sum(y["h"] for y in Y) + 0.25 * sum(y["d"] for y in Y)) * 1000 / sum(y["rain"] for y in Y)
            check(abs(lin["point"] - rnd(hand)) < 1e-9, "端到端：LIN 成员 Q+H＝手算")
            check("judgment" in Se["E5"]["main"] and "smooth_judgment" in Se["E5"]["main"] and Se["B"]["0-6h"]["Q5"]["n"] > 0,
                  "端到端：主判、平滑层与 B 群描述字段齐")
    out = {"script": "p4e_budget.py", "version": VERSION, "n_checks": len(res), "n_fail": len(fails), "checks": res,
           "pass": not fails, "runtime_s": round(time.monotonic() - t_start, 1)}
    if out_dir:
        p2.jdump(out, os.path.join(out_dir, "p4e_selftest.json"))
    for line in res:
        print(line)
    print(f"自测 {len(res) - len(fails)}/{len(res)} 通过（{out['runtime_s']} s）", flush=True)
    return 0 if not fails else 4


def main(argv=None):
    ap = argparse.ArgumentParser(description="P4e 探索性分析 E3′／E5")
    ap.add_argument("--out")
    ap.add_argument("--p1-events", default=p2.P1_EVENTS_DEFAULT)
    ap.add_argument("--p1b-dir", default=p2.P1B_DIR_DEFAULT)
    ap.add_argument("--p2-dir", default=pm.P2_DIR_DEFAULT)
    ap.add_argument("--p4-dir", default=pb.P4_DIR_DEFAULT)
    ap.add_argument("--p4b-dir", default=P4B_DIR_DEFAULT)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    out_dir = args.out or os.environ.get("REPRO_OUTPUT_DIR")
    if args.selftest:
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        return selftest(out_dir)
    if not out_dir:
        print("需要 --out 或环境变量 REPRO_OUTPUT_DIR", file=sys.stderr)
        return 3
    os.makedirs(out_dir, exist_ok=True)
    rc = selftest(out_dir)
    if rc:
        return rc
    log = p1.Log(os.path.join(out_dir, "p4e_log.txt"))
    log.log(f"=== start {VERSION} p2_dir={args.p2_dir} p4_dir={args.p4_dir} p4b_dir={args.p4b_dir} out={out_dir}", echo=True)
    try:
        rc = run(args, out_dir, log)
    except Exception:
        log.log("FATAL：\n" + traceback.format_exc(), echo=True)
        return 3
    log.log("=== done")
    log.close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
