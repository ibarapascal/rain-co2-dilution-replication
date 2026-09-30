#!/usr/bin/env python3
"""p5_beta_mech.py — P5：化学斜率比 r_β＝β_obs/β_model≈1.73 的机制检验（探索性；判读规则、候选、阈值在运行前写定）。

全部结果属探索性分析，不进 P2/P3 判门，不改 X10 的口径。

用法：
  python p5_beta_mech.py --selftest       # 合成数据自测（无网络；秒级；无 PyCO2SYS 时跳过碳酸盐项并记 skipped）
  python p5_beta_mech.py --plan           # 先自测，再重建事件、做 V0／V-S6／V-P4 校验，只报输入可用性计数（不算任何新比值）
  python p5_beta_mech.py                  # 正式运行：先自测，再重建 646 事件并做 V0 校验，再算各候选
  --p1-events PATH / --p1b-dir DIR        同 p2_rim_test（P1/P1b 缓存以符号链接只读引用）
  --p2-dir DIR                            P2 主集输出（读 p2_summary.json、p2_events.csv、p2_glodap.json）
  --p4-events PATH                        P4 第一轮 p4_mech_events.csv（取 u10、g6、c6；sha256 c5feedaa…）
  --p3b-curves PATH / --p3c PATH          S6′ 换算用（P3B whole 曲线、P3C M1 s；V-S6 须复现 8.16555）
  --out DIR                               输出目录（默认 $REPRO_OUTPUT_DIR）
依赖：numpy、scipy、PyCO2SYS 1.8.3.4＋同目录 p1_events.py、p1b_extend.py、p2_rim_test.py（只 import，不改）。
数据：只读现有缓存，无新下载；MAPCO2 ERDDAP 年块 CSV 里本就有 pCO2_air 列（P1 查询即含）。
产物（<out>/）：p5_beta_selftest.json、--plan：p5_beta_plan.json；正式：p5_beta_summary.json（含逐候选 verdict 与总映射）、p5_beta_events.csv、p5_beta_log.txt。
退出码：0 跑完（无论判读）；3 异常（含 V0／V-S6 校验不一致、事件对不上）；4 自测不过。

检验项：
  V0   复算 β_obs(D4)＝23.33897 [16.85003, 29.85017]（p2.ratio_stats 同抽样）、逐事件 ΔpCO₂ 与 ΔS(0.5 m) 对 P2 CSV、
       站 β_model 重算对 P2（1e-3）；V-S6 复算 S6 U＝8.16555（1e-3）；V-P4 P4 CSV 646 事件全对上
  M0   主 r_β＝ΣΔpCO₂/Σβ_st·ΔS(0.5 m)（D4）：站×季簇 bootstrap、站簇 bootstrap、站级 jackknife-t（ln r）、逐站剔除
  (a)  气体交换：A1 雨前海气 ΔpCO₂ 符号分组；A2 逐事件放气/吸气量订正（中心 h=1 m、k×1；上界 h=0.5 m、k×3）；
       A3 时间窗 0–3/3–6/6–12 h；A4 U10 中位数分组
  (b)  平流/锋面：B1 |ΔT|≥0.5 °C（另 0.3）筛除；B2 雨量不足筛除；B3 并集；B4 对照窗背景 ΔpCO₂–ΔS 斜率
  (c)  温度归一：0.0413、0.0433、站点 PyCO2SYS 等化学系数；不归一（信息）
  (d)  选择/构成：D1 A/B 群；D2 逐站剔除、剔 BOBOA；D3 ΔS 强度分层、阈值扫描、按 1 m 选样；D4 选择偏向蒙特卡洛
  (e)  碳酸盐体系：E1 局部导数（0.5 psu）；E2 雨水 DIC 0/12/50；E3 常数 Lueker2000/Millero2010；
       E4 逐事件观测基线（观测雨前 pCO₂、T、S，TA 按 GLODAP TA/S 缩放）局部导数；E4b 同基线 9 psu 拟合
  (f)  分母层位：F1 用 ΔS(1 m) 作分母（A 群真 1 m 另报）
  (g)  背景漂移/日变化：G1 事件减本事件对照窗均值后的 r_β；对照窗 ΔpCO₂ 均值
  Z    组合（描述）：G1＋B3＋E4
  S6′  各订正 r′ 代入 P3C 同式得 U（描述）

实现选择（T 条）：
  T1 事件集＝P2 主集 646 事件（p2.build_all 重建、K1 逐条校验）；D4＝ΔS(0.5 m, 0–6 h)≤−0.2 且 ΔpCO₂ 有效（34 事件）。
  T2 ΔS、ΔT 同 I8（窗内有效值均值 − 雨前 6 h 有效值中位数）；ΔpCO₂ 同 K18（系数 c 可变；c=None＝不归一，不要求 T）。
  T3 r_β 写成单一比值 ΣΔpCO₂／Σ(β_st·ΔS)（与 β_obs／β_model 合并式恒等）；β_st＝站 GLODAP 表层 pCO₂ 斜率（K19 重算）。
  T4 推断：站×季整簇 bootstrap，簇键排序、PCG64(20260926)、B=10000，同一父集内所有变体共用同一抽样矩阵；子集用掩码
     （抽中簇在子集内无事件时该抽样为 NaN，丢弃并报有效数）；百分位 95% 区间。另报站簇 bootstrap（同 seed）、
     逐站剔除点值、站级 jackknife-t（ln r，df＝站数−1）。
  T5 海气 ΔpCO₂_pre＝雨前 6 h pCO₂_sw（原值）中位数 − [t0−24 h, t0+6 h) pCO₂_air 中位数；无大气值则 A1/A2 缺。
  T6 气体交换量：δ＝−RF·(pCO₂_pre/DIC)·K0·ΔpCO₂_pre·(m·k·t̄/h)；k＝0.251·U10²·(Sc/660)^−0.5 cm/h（Wanninkhof 2014，
     Sc 用其海水 CO₂ 多项式、T＝T_pre），t̄＝窗内平均经过时间（0–6 h 取 3 h），RF、K0、DIC 取 E4 逐事件基线
     （E4 不可算时用站 GLODAP 值）；U10 取 P4 的 u10（K9 换算、[t0,t0+6 h)）。订正后分子＝ΔpCO₂−δ。
  T7 雨量不足：0.5 m×(−ΔS/S0)＞1.30×max(G6, C6)/1000（G6 雨量计、C6 CMORPH 中心像元，P4 同窗）；两者皆缺＝不判为可疑。
  T8 对照窗：本事件的 ctrl_hours（同 UTC 时刻、前 24 h 无雨，P1 I9），按同式求 ΔS(0.5 m)、ΔpCO₂(K18)；
     B4 按（站, 对照起点）去重，簇键（站, 对照起点季节）；斜率过原点 Σxy/Σ(β_st·x²)。
  T9 D4 蒙特卡洛：真值 x＝观测 ΔS(0.5 m)（全 646 事件中 ΔpCO₂、ΔS 有效者），y＝m·β_st·x；噪声对 (e_x, e_y) 从本事件
     有效对照窗等概率抽取（本事件无则用同站池）；耦合＝同一对照窗，去耦＝两次独立抽取；按 x+e_x≤−0.2 选样，
     r̂＝Σy*/Σβ_st·x*；B=10000、PCG64(20260926＋偏移)。噪声叠在已含噪的观测上，偏向量级为保守上界。
  T10 E4：TA_e＝TA_GLODAP·S0/S_GLODAP，DIC_e 由 (TA_e, pCO₂_pre) 在 (S0, T_pre) 解出（PyCO2SYS par 1/4，常数 17）；
     稀释同 Witte（TA·SR、DIC·SR＋25·(1−SR)），局部导数＝S0 至 S0−0.5 步长 0.1 的线性拟合斜率，T 固定 T_pre。
  T11 C 的站点等化学系数＝ln(pCO₂(T+0.5)/pCO₂(T−0.5))，TA/DIC/S 取站 GLODAP。
  T12 判读函数 verdict()；不做多重比较校正（逐候选报告）。

Change Log：
  2026-09-27 初版。
"""

import argparse
import csv
import glob
import gzip
import json
import math
import os
import sys
import time
import traceback
from datetime import datetime, timezone

import p1_events as p1
import p2_rim_test as p2

VERSION = "p5beta-2026-09-27a"
HERE = os.path.dirname(os.path.abspath(__file__))
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
P2_DIR_DEFAULT = _rp.upstream("p2_dir")  # [repro] 以下四个默认值读上游阶段产物
P4_EVENTS_DEFAULT = _rp.upstream("p4_events")
P3B_DEFAULT = _rp.upstream("p3b_curves")
P3C_DEFAULT = _rp.upstream("p3c")
NAN = float("nan")

# ---- 统计与窗口 ----
B = p2.BOOT_B
SEED = p2.BOOT_SEED
WIN = p2.PRIMARY_WIN
WINDOWS = {"0-3": (0, 3), "3-6": (3, 6), "6-12": (6, 12)}
TBAR_H = {"0-6": 3.0, "0-3": 1.5, "3-6": 4.5, "6-12": 9.0}
D4_THR = p1.D4_DS_THRESH
T_COEF_MAIN = p2.T_COEF
T_COEF_ALT = {"c0413": 0.0413, "c0433": 0.0433}

# ---- 判读阈值 ----
R_MAJOR = 1.25
PHI_PARTIAL = 0.25
MIN_N, MIN_G = 8, 4
CLEAN_MIN_N, CLEAN_MIN_G = 15, 5
P_MC_MAJOR = 0.05

# ---- 候选参数 ----
DT_SCREEN, DT_SCREEN_ALT = 0.5, 0.3
UNDERCATCH = 1.30
LENS_MIN_M = 0.5
GX_CASES = {"central": {"h_m": 1.0, "k_mult": 1.0}, "upper": {"h_m": 0.5, "k_mult": 3.0}}
AIR_WIN = (-24, 6)
DIC_RAIN_ALT = (0.0, 12.0, 50.0)
K_CARB_ALT = {"Lueker2000": 10, "Millero2010": 14}
LOCAL_DS, LOCAL_STEP = 0.5, 0.1
THRESH_SWEEP = (-0.1, -0.15, -0.2, -0.3, -0.4)
STRATA = (("weak_-0.2_to_-0.1", -0.2, -0.1), ("mid_-0.4_to_-0.2", -0.4, -0.2), ("strong_le_-0.4", -99.0, -0.4))

# ---- 校验容差 ----
V0_TOL = 1e-4
BM_TOL = 1e-3
DP_TOL = 6e-4
DS_TOL = 6e-5
S6_TOL = 1e-3


def _np():
    import numpy
    return numpy


def isnum(x):
    return x is not None and isinstance(x, (int, float)) and math.isfinite(x)


def rnd(x, nd=5):
    return round(float(x), nd) if isnum(x) else None


# ======================================================================== 观测量（T2）
def obs_delta(arr, i, win):
    return p2.obs_delta(arr, i, win)


def pre_median(arr, i):
    return p2.pre_median(arr, i)


def dp_coef(stn, i, win=WIN, coef=T_COEF_MAIN):
    """ΔpCO₂（K18 同式，系数 coef）；coef=None 为不归一（不要求 T 有效）。coef=0.0423 时与 p2.dpco2_norm 逐值相同。"""
    np = _np()
    pc = stn["ser"]["pco2"]
    if coef is None:
        pre = pc[max(0, i - p1.PRE_H):i]
        post = pc[i + win[0]:i + win[1]]
        pre, post = pre[np.isfinite(pre)], post[np.isfinite(post)]
        if not len(pre) or not len(post):
            return NAN
        return float(post.mean() - np.median(pre))
    t = stn["ser"]["sst"]
    tref = pre_median(t, i)
    if not isnum(tref):
        return NAN

    def norm(sl):
        a, b = pc[sl], t[sl]
        ok = np.isfinite(a) & np.isfinite(b)
        return a[ok] * np.exp(coef * (tref - b[ok]))

    pre = norm(slice(max(0, i - p1.PRE_H), i))
    post = norm(slice(i + win[0], i + win[1]))
    if not len(pre) or not len(post):
        return NAN
    return float(post.mean() - np.median(pre))


def seg_median(arr, a, b):
    np = _np()
    x = arr[max(0, a):max(0, b)]
    x = x[np.isfinite(x)]
    return float(np.median(x)) if len(x) else NAN


# ======================================================================== 气体交换（T6）
def schmidt_co2(t):
    """Wanninkhof 2014 海水 CO₂ Schmidt 数多项式（−2–40 °C）。"""
    return 2116.8 - 136.25 * t + 4.7353 * t ** 2 - 0.092307 * t ** 3 + 0.0007555 * t ** 4


def k_w14_ms(u10, t):
    """k（m/s）＝0.251·U10²·(Sc/660)^−0.5 cm/h。"""
    return 0.251 * u10 ** 2 * (schmidt_co2(t) / 660.0) ** -0.5 / 100.0 / 3600.0


def delta_gx(rf, pco2, dic, k0, dpco2_pre, u10, t, h_m, k_mult, tbar_h):
    """窗内平均的气体交换 pCO₂ 变化（µatm）；ΔpCO₂_pre>0（过饱和）→ 放气 → δ<0。"""
    if not all(isnum(v) for v in (rf, pco2, dic, k0, dpco2_pre, u10, t)):
        return NAN
    k = k_mult * k_w14_ms(u10, t)
    d_dic = -k * tbar_h * 3600.0 * k0 * dpco2_pre / h_m            # µmol/kg（K0 以 µmol kg⁻¹ µatm⁻¹ 计）
    return rf * pco2 / dic * d_dic


# ======================================================================== 碳酸盐（T10、T11、E）
def have_pyco2():
    try:
        import PyCO2SYS  # noqa: F401
        return True
    except Exception:
        return False


def pyco2_dilution_slope(ta, dic, sss, t, sals, dic_rain=p2.DIC_RAIN, kc=p2.K_CARBONIC):
    import PyCO2SYS as pyco2
    from scipy import stats
    np = _np()
    sals = np.asarray(sals, float)
    sr = sals / sss
    res = pyco2.sys(par1=ta * sr, par2=dic * sr + dic_rain * (1 - sr), par1_type=1, par2_type=2,
                    salinity=sals, temperature=t, opt_k_carbonic=kc)
    return float(stats.linregress(x=sals, y=res["pCO2"]).slope)


def witte_sals(sss):
    np = _np()
    return np.array([sss - k for k in range(10)]) if sss > 9 else np.arange(sss, 0, -1)


def local_sals(sss):
    np = _np()
    n = int(round(LOCAL_DS / LOCAL_STEP))
    return np.array([sss - LOCAL_STEP * k for k in range(n + 1)])


def station_betas(glodap):
    """站 β 变体（pCO₂ 口径）＋等化学温度系数。glodap：{站: {TAlk, TCO2, salinity, temperature}}。"""
    import PyCO2SYS as pyco2
    np = _np()
    out = {}
    for st, g in glodap.items():
        if not all(isnum(g.get(k)) for k in ("TAlk", "TCO2", "salinity", "temperature")):
            out[st] = None
            continue
        ta, dic, s, t = g["TAlk"], g["TCO2"], g["salinity"], g["temperature"]
        d = {"main": pyco2_dilution_slope(ta, dic, s, t, witte_sals(s)),
             "E1_local": pyco2_dilution_slope(ta, dic, s, t, local_sals(s))}
        for dr in DIC_RAIN_ALT:
            d[f"E2_dicrain{dr:g}"] = pyco2_dilution_slope(ta, dic, s, t, witte_sals(s), dic_rain=dr)
        for nm, kc in K_CARB_ALT.items():
            d[f"E3_{nm}"] = pyco2_dilution_slope(ta, dic, s, t, witte_sals(s), kc=kc)
        r = pyco2.sys(par1=np.array([ta, ta]), par2=np.array([dic, dic]), par1_type=1, par2_type=2,
                      salinity=s, temperature=np.array([t - 0.5, t + 0.5]), opt_k_carbonic=p2.K_CARBONIC)
        d["T_coef_isochem"] = float(math.log(r["pCO2"][1] / r["pCO2"][0]))
        r0 = pyco2.sys(par1=ta, par2=dic, par1_type=1, par2_type=2, salinity=s, temperature=t,
                       opt_k_carbonic=p2.K_CARBONIC)
        d["glodap_pCO2"] = float(r0["pCO2"])
        d["glodap_RF"] = float(r0["revelle_factor"])
        d["glodap_K0"] = float(r0["k_CO2"])
        d["glodap_DIC"] = dic
        out[st] = d
    return out


def event_baseline(ta_g, s_g, s0, t_pre, pco2_pre):
    """E4：逐事件观测基线 → (β_local, β_9psu, DIC_e, RF, K0)；输入不全返回 None。"""
    import PyCO2SYS as pyco2
    if not all(isnum(v) for v in (ta_g, s_g, s0, t_pre, pco2_pre)):
        return None
    ta = ta_g * s0 / s_g
    r = pyco2.sys(par1=ta, par2=pco2_pre, par1_type=1, par2_type=4, salinity=s0, temperature=t_pre,
                  opt_k_carbonic=p2.K_CARBONIC)
    dic = float(r["dic"])
    if not isnum(dic):
        return None
    return {"beta_local": pyco2_dilution_slope(ta, dic, s0, t_pre, local_sals(s0)),
            "beta_9psu": pyco2_dilution_slope(ta, dic, s0, t_pre, witte_sals(s0)),
            "DIC": dic, "TA": ta, "RF": float(r["revelle_factor"]), "K0": float(r["k_CO2"])}


# ======================================================================== MAPCO2 大气 pCO₂（T5）
def read_mapco2_col(cache, erddap_id, col, lo, hi):
    acc = {}
    for path in sorted(glob.glob(os.path.join(cache, "erddap", f"{erddap_id}_*.csv.gz"))):
        if path.endswith("_info.csv.gz"):
            continue
        with gzip.open(path, "rt", encoding="utf-8") as f:
            rd = csv.reader(f)
            hdr = next(rd)
            next(rd)
            if col not in hdr:
                continue
            it, ic = hdr.index("time"), hdr.index(col)
            for r in rd:
                v = p1._num(r[ic])
                if v is None or not (lo < v < hi):
                    continue
                h = p1.iso_to_hour(r[it])
                s = acc.setdefault(h, [0.0, 0])
                s[0] += v
                s[1] += 1
    return {h: s / c for h, (s, c) in acc.items()}


def attach_air(stations, cache):
    np = _np()
    cov = {}
    for stn in stations:
        d = read_mapco2_col(cache, stn["erddap"], "pCO2_air", 50.0, 2000.0)
        a = np.full(stn["n"], np.nan)
        for h, v in d.items():
            i = h - stn["h0"]
            if 0 <= i < stn["n"]:
                a[i] = v
        stn["ser"]["pco2_air"] = a
        cov[stn["name"]] = int(np.isfinite(a).sum())
    return cov


# ======================================================================== 逐事件量
def event_vars(stn, e, p4row, bst, glo):
    """一个事件的全部 P5 派生量（观测部分只在正式运行里调用）。bst：站 β 变体；glo：站 GLODAP 值。"""
    np = _np()
    ser = stn["ser"]
    i = e["i"]
    r = {"station": stn["name"], "group": stn["group"], "season": e["season"], "onset_utc": e["onset_utc"],
         "local_solar_hour": e.get("local_solar_hour", NAN), "i": i}
    r["dS05"] = obs_delta(ser["sss05"], i, WIN)
    r["dS1"] = obs_delta(ser["s1"], i, WIN)
    r["dp"] = dp_coef(stn, i, WIN, T_COEF_MAIN)
    for nm, w in WINDOWS.items():
        r[f"dS05_{nm}"] = obs_delta(ser["sss05"], i, w)
        r[f"dp_{nm}"] = dp_coef(stn, i, w, T_COEF_MAIN)
    r["dp_raw"] = dp_coef(stn, i, WIN, None)
    for nm, c in T_COEF_ALT.items():
        r[f"dp_{nm}"] = dp_coef(stn, i, WIN, c)
    r["dp_isochem"] = dp_coef(stn, i, WIN, bst["T_coef_isochem"]) if bst else NAN
    r["dT"] = obs_delta(ser["sst"], i, WIN)
    r["S0_05"] = pre_median(ser["sss05"], i)
    r["T_pre"] = pre_median(ser["sst"], i)
    r["pco2_pre"] = pre_median(ser["pco2"], i)
    r["pco2_air"] = seg_median(ser["pco2_air"], i + AIR_WIN[0], i + AIR_WIN[1]) if "pco2_air" in ser else NAN
    r["dpco2_pre"] = r["pco2_pre"] - r["pco2_air"] if isnum(r["pco2_pre"]) and isnum(r["pco2_air"]) else NAN
    for k in ("u10", "g6", "c6"):
        v = p4row.get(k) if p4row else None
        try:
            r[k] = float(v)
        except (TypeError, ValueError):
            r[k] = NAN
    # 对照窗（T8）
    pairs = []
    for hc in e.get("ctrl_hours", []):
        c = hc - stn["h0"]
        x = obs_delta(ser["sss05"], c, WIN)
        y = dp_coef(stn, c, WIN, T_COEF_MAIN)
        if isnum(x) and isnum(y):
            pairs.append((x, y, int(c), p1.season_of(hc)))
    r["ctrl_pairs"] = pairs
    r["n_ctrl_valid"] = len(pairs)
    r["ctrl_mean_x"] = float(np.mean([p[0] for p in pairs])) if pairs else NAN
    r["ctrl_mean_y"] = float(np.mean([p[1] for p in pairs])) if pairs else NAN
    # β（T3、E）
    r["bm"] = bst["main"] if bst else NAN
    for k in ("E1_local", "E3_Lueker2000", "E3_Millero2010") + tuple(f"E2_dicrain{d:g}" for d in DIC_RAIN_ALT):
        r["bm_" + k] = bst[k] if bst else NAN
    eb = event_baseline(glo.get("TAlk"), glo.get("salinity"), r["S0_05"], r["T_pre"], r["pco2_pre"]) \
        if (glo and have_pyco2()) else None
    r["bm_E4_local"] = eb["beta_local"] if eb else NAN
    r["bm_E4b_9psu"] = eb["beta_9psu"] if eb else NAN
    r["E4_ok"] = int(eb is not None)
    rf = eb["RF"] if eb else (bst["glodap_RF"] if bst else NAN)
    k0 = eb["K0"] if eb else (bst["glodap_K0"] if bst else NAN)
    dic = eb["DIC"] if eb else (bst["glodap_DIC"] if bst else NAN)
    pco2b = r["pco2_pre"] if isnum(r["pco2_pre"]) else (bst["glodap_pCO2"] if bst else NAN)
    for nm, cs in GX_CASES.items():
        r[f"gx_{nm}"] = delta_gx(rf, pco2b, dic, k0, r["dpco2_pre"], r["u10"], r["T_pre"], cs["h_m"], cs["k_mult"],
                                 TBAR_H["0-6"])
    # 雨量不足（T7）
    rain = max([v for v in (r["g6"], r["c6"]) if isnum(v)], default=NAN)
    r["rain_max6"] = rain
    if isnum(r["dS05"]) and isnum(r["S0_05"]) and isnum(rain):
        r["rain_insufficient"] = int(LENS_MIN_M * (-r["dS05"] / r["S0_05"]) > UNDERCATCH * rain / 1000.0)
    else:
        r["rain_insufficient"] = 0
    return r


# ======================================================================== 推断（T4）
class Boot:
    """站×季（或任意键）整簇 bootstrap；抽样矩阵与 p2.cluster_boot 同式（同键集＝同抽样）。"""

    def __init__(self, keys, seed=SEED, b=B):
        np = _np()
        uk = sorted(set(keys))
        self.kid = {k: j for j, k in enumerate(uk)}
        self.idx = np.array([self.kid[k] for k in keys], int)
        self.K = len(uk)
        rng = np.random.Generator(np.random.PCG64(seed))
        self.draw = rng.integers(0, self.K, size=(b, self.K))

    def ratio(self, num, den, mask=None):
        np = _np()
        num, den = np.asarray(num, float), np.asarray(den, float)
        m = np.isfinite(num) & np.isfinite(den)
        if mask is not None:
            m &= np.asarray(mask, bool)
        n = int(m.sum())
        if n == 0:
            return {"n": 0, "G": 0, "point": None, "draws": None}
        sn = np.bincount(self.idx[m], weights=num[m], minlength=self.K)
        sd = np.bincount(self.idx[m], weights=den[m], minlength=self.K)
        with np.errstate(divide="ignore", invalid="ignore"):
            bs = sn[self.draw].sum(axis=1) / sd[self.draw].sum(axis=1)
        d = float(den[m].sum())
        pt = float(num[m].sum()) / d if d != 0 else NAN
        return {"n": n, "G": int(len(set(self.idx[m].tolist()))), "point": pt, "draws": bs}


def summ(res, label=""):
    np = _np()
    out = {"label": label, "n": res["n"], "G": res["G"], "point": rnd(res["point"])}
    bs = res.get("draws")
    if bs is not None:
        f = bs[np.isfinite(bs)]
        out["n_boot_finite"] = int(len(f))
        if len(f) > 100:
            out["ci95"] = [rnd(np.percentile(f, 2.5)), rnd(np.percentile(f, 97.5))]
            out["se_boot"] = rnd(float(np.std(f, ddof=1)))
    return out


def ln_ratio_ci(ra, rb):
    """同一抽样矩阵下 ln(ra/rb) 的点值与区间。"""
    np = _np()
    if ra["draws"] is None or rb["draws"] is None or not isnum(ra["point"]) or not isnum(rb["point"]) \
            or ra["point"] <= 0 or rb["point"] <= 0:
        return {"point": None}
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.log(ra["draws"]) - np.log(rb["draws"])
    d = d[np.isfinite(d)]
    out = {"point": rnd(math.log(ra["point"] / rb["point"])), "n_boot_finite": int(len(d))}
    if len(d) > 100:
        out["ci95"] = [rnd(np.percentile(d, 2.5)), rnd(np.percentile(d, 97.5))]
    return out


def rows_arr(rows, f):
    np = _np()
    out = []
    for r in rows:
        try:
            v = f(r)
        except (KeyError, TypeError):
            v = NAN
        out.append(v if isnum(v) else NAN)
    return np.array(out, float)


def station_robust(rows, num, den, mask=None):
    """C-9 少簇替代：站簇 bootstrap、逐站剔除、站级 jackknife-t（ln r）。"""
    np = _np()
    from scipy import stats
    m = np.isfinite(num) & np.isfinite(den)
    if mask is not None:
        m &= np.asarray(mask, bool)
    sts = [r["station"] for r in rows]
    bt = Boot(sts)
    res = summ(bt.ratio(num, den, m), "站簇 bootstrap")
    used = sorted({s for s, k in zip(sts, m) if k})
    loso = {}
    for s in used:
        mm = m & np.array([x != s for x in sts])
        d = den[mm].sum()
        loso[s] = rnd(num[mm].sum() / d) if mm.any() and d != 0 else None
    vals = [v for v in loso.values() if isnum(v) and v > 0]
    out = {"station_boot": res, "loso": loso,
           "loso_min": rnd(min(vals)) if vals else None, "loso_max": rnd(max(vals)) if vals else None}
    G = len(used)
    pt = num[m].sum() / den[m].sum() if m.any() and den[m].sum() != 0 else NAN
    if G >= 3 and len(vals) == G and isnum(pt) and pt > 0:
        th = np.log(np.array(vals))
        se = math.sqrt((G - 1) / G * float(((th - th.mean()) ** 2).sum()))
        tq = float(stats.t.ppf(0.975, G - 1))
        out["jackknife_t_lnr"] = {"G": G, "df": G - 1, "se": rnd(se),
                                  "ci95_r": [rnd(math.exp(math.log(pt) - tq * se)), rnd(math.exp(math.log(pt) + tq * se))]}
    return out


# ======================================================================== 判读（T12）
def phi_of(r_corr, r_base):
    if not (isnum(r_corr) and isnum(r_base)) or r_base <= 1.0 or r_corr <= 0:
        return None
    return 1.0 - math.log(r_corr) / math.log(r_base)


def verdict(corr, base_point, sig_ok, min_n=MIN_N, min_g=MIN_G):
    """corr：summ() 结果（订正后 r′）；base_point：同一事件集上的未订正 r；sig_ok：特征检验 True/False/None。"""
    if corr is None or corr.get("n", 0) < min_n or corr.get("G", 0) < min_g or not isnum(corr.get("point")) \
            or "ci95" not in corr:
        return {"label": "不可评", "phi": None}
    phi = phi_of(corr["point"], base_point)
    if sig_ok is False:
        return {"label": "不支持（特征不符）", "phi": rnd(phi)}
    if corr["point"] <= R_MAJOR and corr["ci95"][0] <= 1.0:
        return {"label": "主要来源", "phi": rnd(phi)}
    if phi is not None and phi >= PHI_PARTIAL:
        return {"label": "部分", "phi": rnd(phi)}
    return {"label": "不支持", "phi": rnd(phi)}


def verdict_mc(mc, r_obs):
    if not mc or not isnum(mc.get("mean_r")):
        return {"label": "不可评", "phi": None}
    phi = math.log(mc["mean_r"]) / math.log(r_obs) if mc["mean_r"] > 0 and r_obs > 1 else None
    if mc["p_ge_obs"] >= P_MC_MAJOR:
        return {"label": "主要来源", "phi": rnd(phi)}
    if phi is not None and phi >= PHI_PARTIAL:
        return {"label": "部分", "phi": rnd(phi)}
    return {"label": "不支持", "phi": rnd(phi)}


RANK = {"主要来源": 3, "部分": 2, "不支持": 1, "不支持（特征不符）": 1, "不可评": 0}


def abstract_mapping(V, robust):
    """S6（8.2%）的总映射。"""
    lab = {k: v["label"] for k, v in V.items()}
    if lab.get("a") in ("主要来源", "部分"):
        return {"branch": "①"}
    art = [k for k in ("b", "d", "f", "g") if lab.get(k) == "主要来源"]
    if art:
        return {"branch": "②"}
    if lab.get("e") == "主要来源":
        return {"branch": "③"}
    part = [k for k in "bcdefg" if lab.get(k) == "部分"]
    if part:
        return {"branch": "④"}
    if robust:
        return {"branch": "⑤"}
    return {"branch": "⑥"}


# ======================================================================== 蒙特卡洛（T9）
def mc_selection(x, bm, pools, thr, r_obs, beta_mult=1.0, coupled=True, seed=SEED, b=B, chunk=500):
    np = _np()
    x, bm = np.asarray(x, float), np.asarray(bm, float)
    N = len(x)
    lens = np.array([len(p) for p in pools], int)
    offs = np.concatenate([[0], np.cumsum(lens)[:-1]])
    fx = np.concatenate([np.array([q[0] for q in p], float) for p in pools]) if lens.sum() else np.zeros(0)
    fy = np.concatenate([np.array([q[1] for q in p], float) for p in pools]) if lens.sum() else np.zeros(0)
    rng = np.random.Generator(np.random.PCG64(seed))
    rs, nsel = [], []
    done = 0
    while done < b:
        nb = min(chunk, b - done)
        u = rng.random((nb, N))
        j = offs + np.minimum((u * lens).astype(int), lens - 1)
        if coupled:
            j2 = j
        else:
            u2 = rng.random((nb, N))
            j2 = offs + np.minimum((u2 * lens).astype(int), lens - 1)
        xs = x + fx[j]
        ys = beta_mult * bm * x + fy[j2]
        sel = xs <= thr
        num = (ys * sel).sum(axis=1)
        den = (bm * xs * sel).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            rs.append(num / den)
        nsel.append(sel.sum(axis=1))
        done += nb
    r = np.concatenate(rs)
    ns = np.concatenate(nsel)
    f = r[np.isfinite(r)]
    return {"N_events": int(N), "beta_mult": rnd(beta_mult), "coupled": bool(coupled), "B": int(b),
            "mean_r": rnd(float(f.mean())), "median_r": rnd(float(np.median(f))),
            "pct2.5_97.5": [rnd(np.percentile(f, 2.5)), rnd(np.percentile(f, 97.5))],
            "p_ge_obs": rnd(float((f >= r_obs).mean())), "mean_n_selected": rnd(float(ns.mean()), 2)}


# ======================================================================== S6′
def s6_fn(curves_path, p3c_path):
    np = _np()
    C = json.load(open(curves_path, encoding="utf-8"))
    P3C = json.load(open(p3c_path, encoding="utf-8"))
    sg = np.array(C["global"]["whole"]["s"])
    Iw = np.array(C["global"]["whole"]["I_int_pct"])
    Gw = np.array(C["global"]["whole"]["G_int_pct"])
    c, dep = C["c"], C["deposition_pct"]
    s = P3C["M1"]["s"]

    def U(r):
        return float(dep + (np.interp(s, sg, Iw) + (r - 1) * np.interp(s, sg, Gw)) * (1 - c))
    return U, s, P3C["info"]["S6"]


# ======================================================================== 分析
def analyze(rows, betas_ok=True):
    """rows：全部事件的 event_vars；返回 (S, V)。纯函数（自测可用合成 rows 调用）。"""
    np = _np()
    S = {}
    ev = [r for r in rows if isnum(r["dS05"]) and isnum(r["dp"]) and isnum(r["bm"])]
    d4 = [r for r in ev if r["dS05"] <= D4_THR]
    keys = [(r["station"], r["season"]) for r in d4]
    bt = Boot(keys)
    A = lambda f: rows_arr(d4, f)  # noqa: E731
    dp, ds, bm = A(lambda r: r["dp"]), A(lambda r: r["dS05"]), A(lambda r: r["bm"])
    den = bm * ds
    base = bt.ratio(dp, den)
    r0 = base["point"]
    S["M0"] = summ(base, "r_β＝ΣΔpCO₂/Σβ_st·ΔS(0.5 m)，D4")
    S["M0"]["beta_obs"] = summ(bt.ratio(dp, ds), "β_obs")
    S["M0"]["beta_model"] = summ(bt.ratio(den, ds), "β_model pooled")
    S["M0"]["station_robust"] = station_robust(d4, dp, den)
    S["M0"]["n_stations"] = len({r["station"] for r in d4})

    def same_subset(mask):
        return bt.ratio(dp, den, mask)

    V = {}
    # ---------------- (a) 气体交换
    a = {}
    dpp = A(lambda r: r["dpco2_pre"])
    sup, und = np.isfinite(dpp) & (dpp > 0), np.isfinite(dpp) & (dpp < 0)
    rs, ru = bt.ratio(dp, den, sup), bt.ratio(dp, den, und)
    a["A1_super"], a["A1_under"] = summ(rs, "过饱和"), summ(ru, "欠饱和")
    a["A1_n_air_missing"] = int((~np.isfinite(dpp)).sum())
    a["A1_dpco2_pre_median_by_group"] = {g: rnd(float(np.nanmedian(dpp[np.array([r["group"] == g for r in d4])])))
                                         if np.isfinite(dpp[np.array([r["group"] == g for r in d4])]).any() else None
                                         for g in sorted({r["group"] for r in d4})}
    D = ln_ratio_ci(rs, ru)
    a["A1_ln_super_over_under"] = D
    ok_a1 = a["A1_super"]["n"] >= MIN_N and a["A1_under"]["n"] >= MIN_N
    sig_a = (D["point"] > 0) if (ok_a1 and isnum(D.get("point"))) else None
    a["A1_signature"] = {"evaluable": ok_a1, "direction_ok": sig_a,
                         "significant": bool(ok_a1 and "ci95" in D and D["ci95"][0] > 0)}
    for nm in GX_CASES:
        gx = A(lambda r, nm=nm: r[f"gx_{nm}"])
        m = np.isfinite(gx)
        corr = bt.ratio(dp - gx, den, m)
        b_same = same_subset(m)
        a[f"A2_{nm}"] = {"corrected": summ(corr, f"订正后 r′（{nm}）"), "base_same_subset": summ(b_same),
                         "mean_delta_gx_uatm": rnd(float(np.nanmean(gx))) if m.any() else None,
                         "mean_dp_uatm": rnd(float(dp[m].mean())) if m.any() else None,
                         "verdict": verdict(summ(corr), b_same["point"], sig_a)}
    for nm in WINDOWS:
        dpw, dsw = A(lambda r, nm=nm: r[f"dp_{nm}"]), A(lambda r, nm=nm: r[f"dS05_{nm}"])
        a[f"A3_{nm}"] = summ(bt.ratio(dpw, bm * dsw), f"窗 {nm}")
        a[f"A3_{nm}_super"] = summ(bt.ratio(dpw, bm * dsw, sup), f"窗 {nm} 过饱和")
    a["A3_ln_6-12_over_0-3"] = ln_ratio_ci(bt.ratio(A(lambda r: r["dp_6-12"]), bm * A(lambda r: r["dS05_6-12"])),
                                           bt.ratio(A(lambda r: r["dp_0-3"]), bm * A(lambda r: r["dS05_0-3"])))
    a["A3_ln_6-12_over_0-3_super"] = ln_ratio_ci(
        bt.ratio(A(lambda r: r["dp_6-12"]), bm * A(lambda r: r["dS05_6-12"]), sup),
        bt.ratio(A(lambda r: r["dp_0-3"]), bm * A(lambda r: r["dS05_0-3"]), sup))
    u = A(lambda r: r["u10"])
    if np.isfinite(u).sum() >= 2:
        umed = float(np.nanmedian(u))
        rh, rl = bt.ratio(dp, den, np.isfinite(u) & (u > umed)), bt.ratio(dp, den, np.isfinite(u) & (u <= umed))
        a["A4_u10_median"] = rnd(umed)
        a["A4_high"], a["A4_low"] = summ(rh, "U10>中位数"), summ(rl, "U10≤中位数")
        a["A4_ln_high_over_low"] = ln_ratio_ci(rh, rl)
    S["a_gas_exchange"] = a
    V["a"] = dict(a["A2_central"]["verdict"])
    V["a"]["upper_bound_verdict"] = a["A2_upper"]["verdict"]["label"]
    V["a"]["correction"] = "A2_central"
    V["a"]["r_corr"] = a["A2_central"]["corrected"]
    # ---------------- (b) 平流/锋面
    b_ = {}
    dT = A(lambda r: r["dT"])
    ins = A(lambda r: r["rain_insufficient"]) > 0.5
    b_["dT_summary"] = {"n": int(np.isfinite(dT).sum()), "median": rnd(float(np.nanmedian(dT))) if np.isfinite(dT).any() else None,
                        "n_abs_ge_0.5": int((np.abs(dT) >= DT_SCREEN).sum()), "n_abs_ge_0.3": int((np.abs(dT) >= DT_SCREEN_ALT).sum()),
                        "n_warming_gt_0": int((dT > 0).sum())}
    b_["n_rain_insufficient"] = int(ins.sum())
    sus_t = np.abs(np.nan_to_num(dT, nan=0.0)) >= DT_SCREEN
    sus_t3 = np.abs(np.nan_to_num(dT, nan=0.0)) >= DT_SCREEN_ALT
    for nm, sus in (("B1_dT0.5", sus_t), ("B1_dT0.3", sus_t3), ("B2_rain", ins), ("B3_union", sus_t | ins)):
        b_[nm] = {"clean": summ(bt.ratio(dp, den, ~sus), "剔除后"), "suspect": summ(bt.ratio(dp, den, sus), "可疑"),
                  "n_suspect": int(sus.sum())}
    # B4 背景斜率
    b_["B4_D4_controls"] = background_slope(d4)
    b_["B4_all_controls"] = background_slope(ev)
    bg = b_["B4_D4_controls"]["rho_bg"]
    c_ok = b_["B3_union"]["clean"]
    s_ok = b_["B3_union"]["suspect"]
    sig_b_parts = []
    if "ci95" in bg:
        sig_b_parts.append(bg["ci95"][0] > 1.0)
    if s_ok["n"] >= MIN_N and isnum(s_ok["point"]) and isnum(c_ok["point"]):
        sig_b_parts.append(s_ok["point"] > c_ok["point"])
    sig_b = (any(sig_b_parts) if sig_b_parts else None)
    b_["signature"] = {"parts": sig_b_parts, "ok": sig_b}
    S["b_advection"] = b_
    V["b"] = verdict(c_ok, r0, sig_b, CLEAN_MIN_N, CLEAN_MIN_G)
    V["b"]["correction"] = "B3_union clean"
    V["b"]["r_corr"] = c_ok
    # ---------------- (c) 温度归一
    c_ = {"dT_mean_D4": rnd(float(np.nanmean(dT))) if np.isfinite(dT).any() else None}
    best = None
    for nm in ("c0413", "c0433", "isochem"):
        v = A(lambda r, nm=nm: r[f"dp_{nm}"])
        m = np.isfinite(v)
        corr = summ(bt.ratio(v, den, m), nm)
        bs_ = same_subset(m)
        vd = verdict(corr, bs_["point"], None)
        c_[nm] = {"corrected": corr, "base_same_subset": summ(bs_), "verdict": vd}
        key = (RANK[vd["label"]], vd["phi"] if isnum(vd["phi"]) else -9)
        if best is None or key > best[0]:
            best = (key, nm)
    raw = A(lambda r: r["dp_raw"])
    c_["raw_info"] = {"corrected": summ(bt.ratio(raw, den), "不归一（信息）")}
    S["c_temperature"] = c_
    V["c"] = dict(c_[best[1]]["verdict"])
    V["c"]["correction"] = best[1]
    V["c"]["r_corr"] = c_[best[1]]["corrected"]
    # ---------------- (d) 选择/构成
    d_ = {}
    grp = np.array([r["group"] for r in d4])
    ra, rb = bt.ratio(dp, den, grp == p2.GROUP_A), bt.ratio(dp, den, grp == p2.GROUP_B)
    d_["D1_A"], d_["D1_B"], d_["D1_ln_A_over_B"] = summ(ra, "A 群"), summ(rb, "B 群"), ln_ratio_ci(ra, rb)
    d_["D1_A_station_robust"] = station_robust(d4, dp, den, grp == p2.GROUP_A)
    d_["D1_B_station_robust"] = station_robust(d4, dp, den, grp == p2.GROUP_B)
    stn = np.array([r["station"] for r in d4])
    d_["D2_excl_BOBOA"] = summ(bt.ratio(dp, den, stn != "BOBOA"), "剔除 BOBOA")
    d_["D2_n_by_station"] = {s: int((stn == s).sum()) for s in sorted(set(stn.tolist()))}
    # D3（全事件父集）
    kall = [(r["station"], r["season"]) for r in ev]
    ba = Boot(kall)
    E = lambda f: rows_arr(ev, f)  # noqa: E731
    dpa, dsa, bma = E(lambda r: r["dp"]), E(lambda r: r["dS05"]), E(lambda r: r["bm"])
    dena = bma * dsa
    d_["D3_all_events_info"] = summ(ba.ratio(dpa, dena), "全事件")
    for nm, lo, hi in STRATA:
        d_[f"D3_{nm}"] = summ(ba.ratio(dpa, dena, (dsa > lo) & (dsa <= hi)), nm)
    d_["D3_threshold_sweep"] = {f"{t:g}": summ(ba.ratio(dpa, dena, dsa <= t), f"ΔS≤{t:g}") for t in THRESH_SWEEP}
    ds1a = E(lambda r: r["dS1"])
    d_["D3_select_on_1m"] = summ(ba.ratio(dpa, dena, np.isfinite(ds1a) & (ds1a <= D4_THR)),
                                 "按 ΔS(1 m)≤−0.2 选样、分母仍为 0.5 m")
    # D4 蒙特卡洛
    pools, xs_, bms_ = [], [], []
    st_pool = {}
    for r in ev:
        st_pool.setdefault(r["station"], []).extend([(p[0], p[1]) for p in r["ctrl_pairs"]])
    n_own = n_station = n_none = 0
    for r in ev:
        p = [(q[0], q[1]) for q in r["ctrl_pairs"]]
        if p:
            n_own += 1
        elif st_pool.get(r["station"]):
            p = st_pool[r["station"]]
            n_station += 1
        else:
            n_none += 1
            continue
        pools.append(p)
        xs_.append(r["dS05"])
        bms_.append(r["bm"])
    d_["D4_pool_source"] = {"own_controls": n_own, "station_pool": n_station, "excluded_no_pool": n_none}
    if pools and isnum(r0):
        d_["D4_H0_coupled"] = mc_selection(xs_, bms_, pools, D4_THR, r0, 1.0, True, SEED + 1)
        d_["D4_H0_decoupled"] = mc_selection(xs_, bms_, pools, D4_THR, r0, 1.0, False, SEED + 2)
        d_["D4_H1_coupled_recovery"] = mc_selection(xs_, bms_, pools, D4_THR, r0, r0, True, SEED + 3)
    S["d_selection"] = d_
    V["d"] = verdict_mc(d_.get("D4_H0_coupled"), r0)
    V["d"]["correction"] = "D4_H0_coupled（期望 r̂）"
    loso_vals = S["M0"]["station_robust"]["loso"]
    V["d"]["loso_any_le_R_MAJOR"] = any(isnum(v) and v <= R_MAJOR for v in loso_vals.values())
    # ---------------- (e) 碳酸盐
    e_ = {}
    best = None
    variants = ["E4_local", "E4b_9psu", "E1_local", "E2_dicrain0", "E2_dicrain12", "E2_dicrain50",
                "E3_Lueker2000", "E3_Millero2010"]
    for nm in variants:
        bv = A(lambda r, nm=nm: r["bm_" + nm])
        m = np.isfinite(bv)
        corr = summ(bt.ratio(dp, bv * ds, m), nm)
        bs_ = same_subset(m)
        vd = verdict(corr, bs_["point"], None)
        e_[nm] = {"corrected": corr, "base_same_subset": summ(bs_), "verdict": vd,
                  "beta_model_pooled": summ(bt.ratio(bv * ds, ds, m))}
    e_["E4_n_ok"] = int(A(lambda r: r["E4_ok"]).sum())
    e_["pco2_pre_vs_glodap"] = pco2_level_table(d4)
    S["e_carbonate"] = e_
    V["e"] = dict(e_["E4_local"]["verdict"])
    V["e"]["correction"] = "E4_local"
    V["e"]["r_corr"] = e_["E4_local"]["corrected"]
    V["e"]["sensitivity_labels"] = {nm: e_[nm]["verdict"]["label"] for nm in variants[1:]}
    # ---------------- (f) 分母层位
    f_ = {}
    ds1 = A(lambda r: r["dS1"])
    m1 = np.isfinite(ds1)
    corr = summ(bt.ratio(dp, bm * ds1, m1), "分母 ΔS(1 m)")
    bs_ = same_subset(m1)
    inv = bt.ratio(ds1, ds, m1)
    f_["F1"] = {"corrected": corr, "base_same_subset": summ(bs_), "sum_dS1_over_sum_dS05": summ(inv, "ΣΔS1/ΣΔS05")}
    f_["F1_true_1m_A"] = summ(bt.ratio(dp, bm * ds1, m1 & (grp == p2.GROUP_A)), "A 群真 1 m")
    sig_f = (inv["point"] > 1.0) if isnum(inv["point"]) else None
    f_["signature"] = {"inversion_point": rnd(inv["point"]), "ok": sig_f}
    S["f_denominator"] = f_
    V["f"] = verdict(corr, bs_["point"], sig_f)
    V["f"]["correction"] = "F1"
    V["f"]["r_corr"] = corr
    # ---------------- (g) 背景漂移/日变化
    g_ = {}
    cx, cy = A(lambda r: r["ctrl_mean_x"]), A(lambda r: r["ctrl_mean_y"])
    mg = np.isfinite(cx) & np.isfinite(cy)
    corr = summ(bt.ratio(dp - cy, bm * (ds - cx), mg), "对照扣除后 r′")
    bs_ = same_subset(mg)
    ones = np.ones(len(d4))
    drift_y = summ(bt.ratio(cy, ones, mg), "对照窗 ΔpCO₂ 均值（µatm）")
    drift_x = summ(bt.ratio(cx, ones, mg), "对照窗 ΔS(0.5 m) 均值（psu）")
    g_["G1"] = {"corrected": corr, "base_same_subset": summ(bs_), "ctrl_mean_dp": drift_y, "ctrl_mean_dS05": drift_x,
                "event_mean_dp": summ(bt.ratio(dp, ones, mg)), "event_mean_dS05": summ(bt.ratio(ds, ones, mg))}
    sig_g = (drift_y["point"] < 0) if isnum(drift_y["point"]) else None
    g_["signature"] = {"ctrl_dp_mean": drift_y["point"], "ok": sig_g,
                       "significant": bool("ci95" in drift_y and drift_y["ci95"][1] < 0)}
    lh = A(lambda r: r["local_solar_hour"])
    g_["local_solar_hour_hist_4h"] = {f"{k * 4:02d}-{k * 4 + 4:02d}": int(((lh >= 4 * k) & (lh < 4 * k + 4)).sum())
                                      for k in range(6)}
    S["g_background"] = g_
    V["g"] = verdict(corr, bs_["point"], sig_g)
    V["g"]["correction"] = "G1"
    V["g"]["r_corr"] = corr
    # ---------------- Z 组合（描述）
    bv = A(lambda r: r["bm_E4_local"])
    clean = ~(sus_t | ins)
    mz = mg & np.isfinite(bv) & clean
    S["Z_combined_descriptive"] = {"corrected": summ(bt.ratio(dp - cy, bv * (ds - cx), mz), "G1＋B3＋E4（描述）"),
                                   "base_same_subset": summ(same_subset(mz))}
    # ---------------- 各订正的少簇替代（C-9）
    rob = {}
    for k, (numv, denv, mk) in {
            "a": (dp - A(lambda r: r["gx_central"]), den, np.isfinite(A(lambda r: r["gx_central"]))),
            "b": (dp, den, clean),
            "e": (dp, A(lambda r: r["bm_E4_local"]) * ds, np.isfinite(A(lambda r: r["bm_E4_local"]))),
            "f": (dp, bm * ds1, m1),
            "g": (dp - cy, bm * (ds - cx), mg)}.items():
        rob[k] = station_robust(d4, numv, denv, mk)
    S["corrections_station_robust"] = rob
    # ---------------- 总映射
    cand = [(k, v) for k, v in V.items() if RANK.get(v["label"], 0) >= 2]
    cand.sort(key=lambda kv: (RANK[kv[1]["label"]], kv[1]["phi"] if isnum(kv[1]["phi"]) else -9), reverse=True)
    mr = S["M0"]["station_robust"]
    corr_lo = [V[k]["r_corr"]["ci95"][0] for k in ("a", "b", "c", "e", "f", "g")
               if isinstance(V[k].get("r_corr"), dict) and "ci95" in V[k]["r_corr"]]
    robust = bool(S["M0"].get("ci95") and S["M0"]["ci95"][0] > 1.0
                  and mr["station_boot"].get("ci95") and mr["station_boot"]["ci95"][0] > 1.0
                  and isnum(mr["loso_min"]) and mr["loso_min"] > 1.0
                  and "jackknife_t_lnr" in mr and mr["jackknife_t_lnr"]["ci95_r"][0] > 1.0
                  and corr_lo and all(x > 1.0 for x in corr_lo))
    S["synthesis"] = {"ranked_candidates_ge_partial": [k for k, _ in cand],
                      "most_likely": cand[0][0] if cand else None,
                      "r_beta_robust_all": robust,
                      "s6_8p2": abstract_mapping(V, robust)}
    S["verdicts"] = V
    return S


def background_slope(rows):
    """B4：对照窗（按站＋起点去重）ΔpCO₂–ΔS(0.5 m) 过原点斜率与模型斜率之比。"""
    seen = {}
    for r in rows:
        for x, y, c, sea in r["ctrl_pairs"]:
            seen[(r["station"], c)] = (x, y, r["bm"], (r["station"], sea))
    if len(seen) < 2:
        return {"n": len(seen), "rho_bg": {"n": len(seen)}}
    vals = list(seen.values())
    bt = Boot([v[3] for v in vals])
    x = rows_arr(vals, lambda v: v[0])
    y = rows_arr(vals, lambda v: v[1])
    bm = rows_arr(vals, lambda v: v[2])
    np = _np()
    ones = np.ones(len(vals))
    return {"n_windows": len(vals),
            "b_bg_uatm_per_psu": summ(bt.ratio(x * y, x * x), "Σxy/Σx²"),
            "rho_bg": summ(bt.ratio(x * y, bm * x * x), "Σxy/Σβ_st·x²"),
            "mean_x": summ(bt.ratio(x, ones)), "mean_y": summ(bt.ratio(y, ones)),
            "sd_x": rnd(float(np.std(x, ddof=1))), "sd_y": rnd(float(np.std(y, ddof=1))),
            "corr_xy": rnd(float(np.corrcoef(x, y)[0, 1])) if np.std(x) > 0 and np.std(y) > 0 else None}


def pco2_level_table(d4):
    np = _np()
    out = {}
    for s in sorted({r["station"] for r in d4}):
        sub = [r for r in d4 if r["station"] == s]
        pp = [r["pco2_pre"] for r in sub if isnum(r["pco2_pre"])]
        out[s] = {"n": len(sub), "pco2_pre_median": rnd(float(np.median(pp))) if pp else None,
                  "glodap_pco2": rnd(sub[0].get("_glodap_pco2")),
                  "bm_main": rnd(sub[0]["bm"]),
                  "bm_E4_local_median": rnd(float(np.nanmedian([r["bm_E4_local"] for r in sub])))
                  if any(isnum(r["bm_E4_local"]) for r in sub) else None}
    return out


# ======================================================================== 主流程
def load_p4(path):
    rows = {}
    with open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            rows[(r["station"], r["onset_utc"])] = r
    return rows


def sha256(path):
    return p2.sha256_file(path)


def run(args, out_dir, log, plan_only=False):
    np = _np()
    t0 = time.monotonic()
    info = {"script": "p5_beta_mech.py", "version": VERSION,
            "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "inputs": {"p2_dir": args.p2_dir, "p4_events": args.p4_events, "p4_events_sha256": sha256(args.p4_events),
                       "p3b_curves_sha256": sha256(args.p3b_curves), "p3c_sha256": sha256(args.p3c)}}
    stations, ev, val, _prim, _extra, plan = p2.build_all(args, out_dir, log)
    events = ev[p2.PRIMARY_THRESHOLD]
    info["rebuild_validation"] = {k: val[k] for k in ("n_ref", "n_rebuilt", "d4_count")}
    info["downloads"] = {"cache_http_requests": plan.get("cache_http_requests"), "cache_links": plan.get("cache_links")}
    info["pco2_air_hours_by_station"] = attach_air(stations, os.path.join(out_dir, "cache"))
    glo_all = json.load(open(os.path.join(args.p2_dir, "p2_glodap.json"), encoding="utf-8"))["stations"]
    p2s = json.load(open(os.path.join(args.p2_dir, "p2_summary.json"), encoding="utf-8"))
    betas = station_betas(glo_all)
    # V0-β：站 β 重算
    bchk = {s: (betas[s]["main"] if betas[s] else None, (p2s["beta_model_by_station"].get(s) or {}).get("pCO2_slope"))
            for s in betas}
    bad = {s: v for s, v in bchk.items() if (v[0] is None) != (v[1] is None)
           or (v[0] is not None and abs(v[0] - v[1]) > BM_TOL)}
    if bad:
        raise RuntimeError(f"V0-β 站 β_model 重算不一致：{bad}")
    p4 = load_p4(args.p4_events)
    sbn = {s["name"]: s for s in stations}
    rows = []
    miss_p4 = 0
    for e in events:
        stn = sbn[e["_stn"]]
        k = (stn["name"], e["onset_utc"])
        pr = p4.get(k)
        miss_p4 += pr is None
        r = event_vars(stn, e, pr, betas.get(stn["name"]), glo_all.get(stn["name"]))
        r["s1_substitute"] = bool(e.get("s1_substitute"))
        r["_glodap_pco2"] = betas[stn["name"]]["glodap_pCO2"] if betas.get(stn["name"]) else NAN
        rows.append(r)
    if miss_p4:
        raise RuntimeError(f"V-P4：{miss_p4} 个事件在 P4 CSV 里找不到")
    # V0：逐事件对 P2 CSV
    p2ev = {(r["station"], r["onset_utc"]): r for r in p2.read_event_csv(os.path.join(args.p2_dir, "p2_events.csv"))}
    nbad = 0
    for r in rows:
        q = p2ev.get((r["station"], r["onset_utc"]))
        if q is None:
            nbad += 1
            continue
        for mine, col, tol in ((r["dp"], "dpCO2_Tnorm", DP_TOL), (r["dS05"], "dS05_obs", DS_TOL)):
            ref = float(q[col]) if q[col] not in ("", "None") else None
            if (ref is None) != (not isnum(mine)) or (ref is not None and abs(mine - ref) > tol):
                nbad += 1
    if nbad or len(rows) != len(p2ev):
        raise RuntimeError(f"V0 逐事件与 P2 CSV 不一致：{nbad} 处（行数 {len(rows)} vs {len(p2ev)}）")
    # V0：β_obs D4 与 P2 同抽样
    d4recs = [{"station": r["station"], "season": r["season"], "dpco2": r["dp"], "obs": {(0.5, WIN): r["dS05"]}}
              for r in rows if isnum(r["dS05"]) and r["dS05"] <= D4_THR]
    v0 = p2.ratio_stats(d4recs, lambda r: r["dpco2"], lambda r: r["obs"][(0.5, WIN)], "V0")
    ref = p2s["primary"]["beta_D4subset"]
    v0_ok = (v0["n_events"] == ref["n_events"] and abs(v0["point"] - ref["point"]) <= V0_TOL
             and all(abs(a - b) <= V0_TOL for a, b in zip(v0["ci95"], ref["ci95"])))
    info["V0"] = {"this": {k: v0[k] for k in ("n_events", "point", "ci95", "n_clusters")},
                  "p2": {k: ref[k] for k in ("n_events", "point", "ci95", "n_clusters")}, "pass": bool(v0_ok),
                  "station_beta_recomputed": {s: rnd(v[0], 4) for s, v in bchk.items()}}
    if not v0_ok:
        raise RuntimeError(f"V0 β_obs 复算不一致：{info['V0']}")
    # V-S6
    U, s_m1, s6 = s6_fn(args.p3b_curves, args.p3c)
    u_chk = U(s6["r_beta"])
    info["V_S6"] = {"U_recomputed": rnd(u_chk), "U_p3c": rnd(s6["U"]), "s_M1": rnd(s_m1),
                    "pass": bool(abs(u_chk - s6["U"]) <= S6_TOL)}
    if not info["V_S6"]["pass"]:
        raise RuntimeError(f"V-S6 不一致：{info['V_S6']}")
    info["station_betas"] = {s: ({k: rnd(v, 5) for k, v in d.items()} if d else None) for s, d in betas.items()}
    if plan_only:                                  # 只报输入可用性计数，不算任何比值、分组或筛除计数
        d4m = [r for r in rows if isnum(r["dS05"]) and r["dS05"] <= D4_THR and isnum(r["dp"])]
        fin = lambda rs, k: int(sum(isnum(r.get(k)) for r in rs))  # noqa: E731
        info["availability_counts_only"] = {
            grp: {"n": len(rs), **{k: fin(rs, k) for k in ("dpco2_pre", "pco2_air", "u10", "g6", "c6", "dT", "dS1",
                                                              "dp_raw", "dp_isochem", "bm_E4_local", "gx_central",
                                                              "ctrl_mean_y", "dS05_6-12", "dp_6-12")}}
            for grp, rs in (("all", rows), ("D4", d4m))}
        info["n_ctrl_valid_hist_D4"] = {str(k): sum(1 for r in d4m if r["n_ctrl_valid"] == k) for k in range(6)}
        info["runtime_s"] = round(time.monotonic() - t0, 1)
        p2.jdump(info, os.path.join(out_dir, "p5_beta_plan.json"))
        print(json.dumps({"plan": "ok", "V0": info["V0"]["pass"], "V_S6": info["V_S6"]["pass"]}, ensure_ascii=False),
              flush=True)
        return 0
    S = dict(info)
    S.update(analyze(rows))
    # S6′（描述）
    s6p = {"M0": rnd(U(S["M0"]["point"]))}
    for k, v in S["verdicts"].items():
        rc = v.get("r_corr")
        if isinstance(rc, dict) and isnum(rc.get("point")):
            s6p[k] = {"r": rc["point"], "U": rnd(U(rc["point"])),
                      "U_at_ci": [rnd(U(x)) for x in rc["ci95"]] if "ci95" in rc else None}
    s6p["r_eq_1"] = rnd(U(1.0))
    S["S6_prime_descriptive"] = s6p
    S["runtime_s"] = round(time.monotonic() - t0, 1)
    p2.jdump(S, os.path.join(out_dir, "p5_beta_summary.json"))
    cols = [k for k in rows[0] if k not in ("ctrl_pairs",)]
    with open(os.path.join(out_dir, "p5_beta_events.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols + ["in_D4"])
        w.writeheader()
        for r in rows:
            d = {k: (rnd(r[k], 6) if isinstance(r[k], float) else r[k]) for k in cols}
            d["in_D4"] = int(isnum(r["dS05"]) and r["dS05"] <= D4_THR and isnum(r["dp"]))
            w.writerow(d)
    print(json.dumps({"verdicts": {k: v["label"] for k, v in S["verdicts"].items()}, "synthesis": S["synthesis"]},
                     ensure_ascii=False), flush=True)
    del np
    return 0


# ======================================================================== 自测（合成数据，无网络）
def _syn_rows(seed=7, n_st=6, per=8, beta_true_mult=1.0, drift=0.0, noise=0.0, bg_slope=0.0):
    """合成事件行：y＝mult·β·x＋drift＋背景；对照窗＝(e_x, bg_slope·e_x＋drift)。"""
    import numpy as np
    rng = np.random.Generator(np.random.PCG64(seed))
    rows = []
    for s in range(n_st):
        bm = 12.0 + s
        for k in range(per):
            x = -0.2 - 0.4 * rng.random() if k % 2 == 0 else -0.05 * rng.random()
            ctrl = []
            for c in range(5):
                ex = noise * rng.standard_normal()
                ctrl.append((ex, bg_slope * ex + drift, 1000 * s + 10 * k + c, "DJF"))
            y = beta_true_mult * bm * x + drift
            r = {"station": f"S{s}", "group": p2.GROUP_A if s < 3 else p2.GROUP_B, "season": "DJF" if k < 4 else "JJA",
                 "onset_utc": f"2000-01-{k + 1:02d}T00:00Z", "local_solar_hour": float((3 * k) % 24), "i": k,
                 "dS05": x, "dS1": 1.3 * x, "dp": y, "bm": bm,
                 "dS05_0-3": x, "dS05_3-6": x, "dS05_6-12": x, "dp_0-3": y, "dp_3-6": y, "dp_6-12": y,
                 "dp_raw": y, "dp_c0413": y, "dp_c0433": y, "dp_isochem": y, "dT": -0.1 if k != 2 else 0.8,
                 "S0_05": 35.0, "T_pre": 28.0, "pco2_pre": 420.0, "pco2_air": 400.0 if s < 3 else 440.0,
                 "u10": 4.0 + k, "g6": 10.0, "c6": 8.0, "ctrl_pairs": ctrl, "n_ctrl_valid": 5,
                 "ctrl_mean_x": float(np.mean([c[0] for c in ctrl])), "ctrl_mean_y": float(np.mean([c[1] for c in ctrl])),
                 "bm_E1_local": bm, "bm_E3_Lueker2000": bm, "bm_E3_Millero2010": bm, "bm_E2_dicrain0": bm,
                 "bm_E2_dicrain12": bm, "bm_E2_dicrain50": bm, "bm_E4_local": bm * beta_true_mult,
                 "bm_E4b_9psu": bm, "E4_ok": 1, "gx_central": 0.0, "gx_upper": 0.0, "rain_max6": 10.0,
                 "rain_insufficient": 0, "s1_substitute": False, "_glodap_pco2": 400.0}
            r["dpco2_pre"] = r["pco2_pre"] - r["pco2_air"]
            rows.append(r)
    return rows


def selftest(out_dir=None):
    import numpy as np
    res = []

    def check(cond, msg):
        res.append({"ok": bool(cond), "msg": msg})

    # 1 比值与抽样：y＝2βx → r＝2（点与区间退化）
    rows = _syn_rows(beta_true_mult=2.0)
    S = analyze(rows)
    check(abs(S["M0"]["point"] - 2.0) < 1e-9 and abs(S["M0"]["ci95"][0] - 2.0) < 1e-9 and abs(S["M0"]["ci95"][1] - 2.0) < 1e-9,
          "M0：y＝2βx 时 r＝2 且区间退化")
    # 2 Boot 与 p2.cluster_boot 同抽样
    keys = [(r["station"], r["season"]) for r in rows]
    num = [r["dp"] + 0.1 * j for j, r in enumerate(rows)]
    den = [r["dS05"] for r in rows]
    ref, _ = p2.cluster_boot(keys, {"m": (num, den)})
    mine = Boot(keys).ratio(np.array(num), np.array(den))["draws"]
    check(np.allclose(ref["m"], mine, equal_nan=True), "Boot 抽样与 p2.cluster_boot 逐值相同")
    # 3 dp_coef 与 p2.dpco2_norm
    n = 60
    ser = {"pco2": np.full(n, np.nan), "sst": np.full(n, np.nan), "sss05": np.full(n, np.nan)}
    rng = np.random.Generator(np.random.PCG64(3))
    for k in range(0, n, 3):
        ser["pco2"][k] = 400 + 5 * rng.standard_normal()
        ser["sst"][k] = 28 + 0.3 * rng.standard_normal()
    stn = {"ser": ser}
    check(all(abs(dp_coef(stn, i, WIN, T_COEF_MAIN) - p2.dpco2_norm(stn, i)) < 1e-12 for i in (10, 20, 31)),
          "dp_coef(0.0423) 与 p2.dpco2_norm 逐值相同")
    pre = ser["pco2"][14:20]
    post = ser["pco2"][20:26]
    check(abs(dp_coef(stn, 20, WIN, None) - (np.nanmean(post) - np.nanmedian(pre))) < 1e-12, "dp_coef(None)＝原值差")
    # 4 G1：漂移 d 同加于事件与对照 → 对照扣除后 r 回到 1，未扣除时偏离
    rows = _syn_rows(beta_true_mult=1.0, drift=-3.0)
    S = analyze(rows)
    check(abs(S["verdicts"]["g"]["r_corr"]["point"] - 1.0) < 1e-9 and S["M0"]["point"] > 1.3,
          "G1：共同漂移下未扣除 r>1.3、扣除后 r＝1")
    check(S["verdicts"]["g"]["label"] == "主要来源" and S["g_background"]["signature"]["ok"] is True,
          "G1 判读为主要来源、特征成立")
    # 5 MC：无噪声 → r̂＝1；独立噪声 → 衰减（<1）；陡背景斜率耦合 → >1
    x = np.array([-0.5, -0.3, -0.25, -0.1, -0.05, 0.02])
    bm = np.full(len(x), 13.0)
    z = [[(0.0, 0.0)] * 5 for _ in x]
    mc0 = mc_selection(x, bm, z, D4_THR, 1.7, b=200)
    check(abs(mc0["mean_r"] - 1.0) < 1e-9, "MC：无噪声 r̂＝1")
    rng = np.random.Generator(np.random.PCG64(5))
    pool_ind = [[(0.1 * rng.standard_normal(), 0.5 * rng.standard_normal()) for _ in range(40)] for _ in x]
    mc1 = mc_selection(x, bm, pool_ind, D4_THR, 1.7, coupled=False, b=2000)
    check(mc1["mean_r"] < 1.0, "MC：独立噪声下选样使 r̂ 偏向 0（<1）")
    pool_c = []
    for _ in x:
        ee = 0.1 * rng.standard_normal(40)
        pool_c.append([(float(a), float(60.0 * a)) for a in ee])
    mc2 = mc_selection(x, bm, pool_c, D4_THR, 1.7, coupled=True, b=2000)
    check(mc2["mean_r"] > 1.0, "MC：陡背景斜率（60）耦合噪声使 r̂>1")
    # 6 气体交换
    check(abs(schmidt_co2(20.0) - 668.3) < 1.0, "Sc(20 °C)≈668（W14 表）")
    check(abs(k_w14_ms(10.0, 20.0) * 3600 * 100 * (668.34 / 660.0) ** 0.5 - 25.1) < 0.05, "k660(10 m/s)＝25.1 cm/h")
    dg = delta_gx(10.0, 450.0, 2000.0, 0.03, 50.0, 7.0, 27.0, 1.0, 1.0, 3.0)
    check(dg < 0 and delta_gx(10.0, 450.0, 2000.0, 0.03, -50.0, 7.0, 27.0, 1.0, 1.0, 3.0) > 0,
          "δ_gx：过饱和为负、欠饱和为正")
    check(abs(delta_gx(10.0, 450.0, 2000.0, 0.03, 50.0, 7.0, 27.0, 0.5, 3.0, 3.0) / dg - 6.0) < 1e-9,
          "δ_gx 与 k 倍数成正比、与 h 成反比")
    # 7 A1 分组方向：过饱和站设 r 更大
    rows = _syn_rows()
    for r in rows:
        if r["dpco2_pre"] > 0:
            r["dp"] *= 2.0
            for w in WINDOWS:
                r[f"dp_{w}"] = r["dp"]
    S = analyze(rows)
    check(S["a_gas_exchange"]["A1_ln_super_over_under"]["point"] > 0.6, "A1：过饱和组 r 更大时 ln 比>0")
    # 8 平流筛：可疑事件（|ΔT|≥0.5）设为高 r → 剔除后回到 1
    rows = _syn_rows()
    for r in rows:
        if r["dT"] >= 0.5:
            r["dp"] *= 5.0
    S = analyze(rows)
    check(abs(S["b_advection"]["B1_dT0.5"]["clean"]["point"] - 1.0) < 1e-9 and S["M0"]["point"] > 1.0,
          "B1：剔除 |ΔT|≥0.5 后 r＝1")
    # 9 雨量不足判定
    r = {"dS05": -0.5, "S0_05": 35.0, "g6": 1.0, "c6": 2.0}
    lhs = LENS_MIN_M * 0.5 / 35.0
    check(lhs > UNDERCATCH * 2.0 / 1000 and not (LENS_MIN_M * 0.1 / 35.0 > UNDERCATCH * 20.0 / 1000),
          "雨量不足判定两例")
    # 10 判读函数
    check(verdict({"n": 30, "G": 10, "point": 1.1, "ci95": [0.8, 1.4]}, 1.73, None)["label"] == "主要来源", "verdict 主要来源")
    check(verdict({"n": 30, "G": 10, "point": 1.45, "ci95": [1.1, 1.9]}, 1.73, None)["label"] == "部分", "verdict 部分")
    check(verdict({"n": 30, "G": 10, "point": 1.7, "ci95": [1.2, 2.2]}, 1.73, None)["label"] == "不支持", "verdict 不支持")
    check(verdict({"n": 5, "G": 3, "point": 1.0, "ci95": [0.5, 1.5]}, 1.73, None)["label"] == "不可评", "verdict 不可评")
    check(verdict({"n": 30, "G": 10, "point": 1.0, "ci95": [0.8, 1.2]}, 1.73, False)["label"] == "不支持（特征不符）",
          "verdict 特征不符")
    check(verdict_mc({"mean_r": 1.2, "p_ge_obs": 0.2}, 1.73)["label"] == "主要来源", "verdict_mc 主要来源")
    # 11 jackknife：各站同 r → se＝0
    rows = _syn_rows(beta_true_mult=1.5)
    sr = station_robust(rows, rows_arr(rows, lambda r: r["dp"]), rows_arr(rows, lambda r: r["bm"] * r["dS05"]))
    check(sr["jackknife_t_lnr"]["se"] < 1e-9 and abs(sr["loso_min"] - 1.5) < 1e-9, "jackknife/LOSO 常数情形")
    # 12 B4：背景斜率 60、β＝13 → Σxy/Σx²＝60
    rows = _syn_rows(noise=0.1, bg_slope=60.0)
    bg = background_slope(rows)
    check(abs(bg["b_bg_uatm_per_psu"]["point"] - 60.0) < 1e-9, "B4 过原点斜率复原 60")
    # 13 映射
    check(abstract_mapping({"a": {"label": "部分"}}, False)["branch"] == "①"
          and abstract_mapping({k: {"label": "不支持"} for k in "abcdefg"}, True)["branch"] == "⑤"
          and abstract_mapping({"a": {"label": "不支持"}, "g": {"label": "主要来源"}}, False)["branch"] == "②",
          "S6 映射分支")
    # 14 碳酸盐（需 PyCO2SYS）
    if have_pyco2():
        ta, dic, s, t = 2294.2332, 2042.6534, 34.9647, 24.2313
        b0 = pyco2_dilution_slope(ta, dic, s, t, witte_sals(s))
        check(abs(b0 - 18.0582) < 1e-3, "β_model(TAO110W GLODAP)＝18.0582（P2 值）")
        eb = event_baseline(ta, s, s, t, 450.0)
        import PyCO2SYS as pyco2
        chk = pyco2.sys(par1=eb["TA"], par2=eb["DIC"], par1_type=1, par2_type=2, salinity=s, temperature=t,
                        opt_k_carbonic=p2.K_CARBONIC)
        check(abs(float(chk["pCO2"]) - 450.0) < 0.01 and eb["beta_local"] > 0, "E4 逐事件基线反解 pCO₂ 自洽")
        sb = station_betas({"X": {"TAlk": ta, "TCO2": dic, "salinity": s, "temperature": t}})["X"]
        check(0.035 < sb["T_coef_isochem"] < 0.05, "等化学温度系数在 0.035–0.05")
        pyco2_status = "run"
    else:
        pyco2_status = "skipped（无 PyCO2SYS；workstation 正式运行前的自测会跑）"
    n_ok = sum(r["ok"] for r in res)
    out = {"version": VERSION, "n_checks": len(res), "n_ok": n_ok, "pyco2_checks": pyco2_status, "checks": res}
    if out_dir:
        p2.jdump(out, os.path.join(out_dir, "p5_beta_selftest.json"))
    print(json.dumps({"n_checks": len(res), "n_ok": n_ok, "pyco2": pyco2_status,
                      "failed": [r["msg"] for r in res if not r["ok"]]}, ensure_ascii=False), flush=True)
    return 0 if n_ok == len(res) else 4


def main(argv=None):
    ap = argparse.ArgumentParser(description="P5 r_β 机制检验")
    ap.add_argument("--out")
    ap.add_argument("--p1-events", default=p2.P1_EVENTS_DEFAULT)
    ap.add_argument("--p1b-dir", default=p2.P1B_DIR_DEFAULT)
    ap.add_argument("--p2-dir", default=P2_DIR_DEFAULT)
    ap.add_argument("--p4-events", default=P4_EVENTS_DEFAULT)
    ap.add_argument("--p3b-curves", default=P3B_DEFAULT)
    ap.add_argument("--p3c", default=P3C_DEFAULT)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--plan", action="store_true")
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
    if not have_pyco2():
        print("正式运行需要 PyCO2SYS", file=sys.stderr)
        return 3
    log = p1.Log(os.path.join(out_dir, "p5_beta_log.txt"))
    log.log(f"=== start {VERSION} p2_dir={args.p2_dir} out={out_dir}", echo=True)
    try:
        rc = run(args, out_dir, log, plan_only=args.plan)
    except Exception:
        log.log("FATAL：\n" + traceback.format_exc(), echo=True)
        return 3
    log.log("=== done")
    log.close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
