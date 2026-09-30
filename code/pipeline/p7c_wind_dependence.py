#!/usr/bin/env python3
"""p7c_wind_dependence.py — P7c 探索性分析（非盲）：矛盾 C——R_RIM(1 m) 随 U10 下降（0.37→0.27→0.16），
而「RIM 淡水不守恒」解释预测它应随风速上升（P4f：ρ 1.85→1.07→0.56）。本脚本对事先写定的候选解释（H1–H6）逐条检验。

判据在运行前写定；判读在 judge_*()。第一轮 negative_L、P4b「不可评」、P4e、P4f 都不改。

用法：
  python p7c_wind_dependence.py --formula      # 第一步公式分析（合成输入，无数据，约 1 s）→ stdout JSON（或 --out 目录下 p7c_formula.json）
  python p7c_wind_dependence.py --selftest     # 合成数据自测（无网络、秒级）
  python p7c_wind_dependence.py                # 正式（先跑自测，不过即退出 4）
  --p1-events / --p1b-dir / --p2-dir           同 p4_mech（默认 P1/P1b/P2 产物）
  --p4-dir DIR                                 第一轮输出目录（p4_mech_events.csv 逐事件核对 V1；p4_mech_summary.json 切点与 V0）
  --p4b-dir DIR                                P4b 输出目录（p4b_events.csv：hfine100、pcum_cmorph_model → 逐事件 M_i，同 P4f；rimG1 核对）
  --out DIR                                    输出目录（默认 $REPRO_OUTPUT_DIR）
依赖：numpy、scipy＋同目录 p1_events、p1b_extend、p2_rim_test、p4_mech、p4b_posthoc（只 import，不改）。
数据：只读现有缓存（P1/P1b 缓存、P2 主集 cmorph_pixels.csv、第一轮与 P4b 产物 CSV/JSON）；无新下载。
产物（<out>/）：p7c_selftest.json、p7c_formula.json、p7c_summary.json（含逐候选判读）、p7c_events.csv、p7c_log.txt。
退出码：0 跑完（无论判读）；3 异常或核对（V0／V1／V-RIM／V-P4F）不过；4 自测不过。

实现选择（Y 条）：
  Y1 事件：pm.build＋pm.event_mech 重建 646 事件；pb.check_v1 与第一轮 CSV 逐事件核对（不过退出 3）。分析集＝model_ok 且
     u10、obs1、rim1 有限（与第一轮 D-U 同集，646）。
  Y2 风档：第一轮 D-U 切点（p4_mech_summary.json D_U.cuts，须与 6.031179／9.325534 在 1e-5 内一致），k=0 低、2 高；
     修订 1：分档用 646 集 u10 的精确 1/3、2/3 分位（第一轮即如此；6 位小数的切点会把恰在分位上的事件分错档）；
     全部「趋势」量＝ln(比值_高)−ln(比值_低)，比值＝Σ分子/Σ分母（同 P2）。中档只报不判。
  Y3 RIM 重新实现 rim_parts()/rim_F()：与 p2.rim_factor 同式（K2–K9、K28）；自测对 p2.rim_factor 逐值 max|Δ|≤1e-12，
     正式运行对每事件第一轮 rim1 复算 max|Δ|≤1e-7（V-RIM，不过退出 3）。线性化 rim_F_lin＝1−(Σterm/d0＋cterm/d0c)（去饱和）。
  Y4 守恒 Gaussian 族 cons_G(κ)：同 48 个历史半步＋当前半步、同风、同 K_z＝κ·2.5e-5·U²，每个雨脉冲 p（m）按
     p/√(πK t)·exp(−z²/4Kt) 分布（深度积分＝p，守恒），历史 t＝ti（1800…86400 s），当前半步 t＝900 s；不含 d0、c1、c2、tc。
     κ=1 即「RIM 的扩散形状＋淡水守恒」。以 F_G＝1−G 走与 RIM 相同的小时化与窗差（K5、K6）。
  Y5 逐事件 M_i＝hfine100·1000/pcum_cmorph_model（P4b CSV，同 P4f）；守恒预测趋势 T_pred＝ln[Σ(rim1/M)/Σrim1]_高−同_低；
     V-P4F：580 集三档 R_obs、R_pred、ρ 须与 P4f 值在 1e-4 内一致（不过退出 3）。
  Y6 候选的逐事件量：
     H1a 形状：A＝T(Σobs5/Σobs1)−T(ΣG5/ΣG1)（参照＝守恒 Gaussian κ=1；对 RIM 的同式只作描述）；
     H1b 时间：B＝T(obs1/G1, [6,12) h)−T(obs1/G1, [0,3) h)（p2.obs_delta 与同一 G 小时序列；对 RIM 的同式只作描述）；
     H2 雨量计驱动 rimG1（p2.event_forcing 'gauge'，负值与缺测置 0，同 P2 敏感性／P4b）；
     H3 站层（Mantel–Haenszel 率比，层＝站）；H4 构成层（层＝24 h 雨量三分位〔646 集切点〕×替代层旗标 |z1−1|>0.05 m）；
     H5 对照窗订正 obs1−mean(ctrl_ds1)（P1 I9 对照，≥1 个有效对照）；H6a 分母换海面 rim0；H6b 分母换线性化 RIM；
     H6c 分母换守恒 Gaussian κ=1 的 G1（对 P4f 的 1/M 近似）。
  Y7 推断：站×季整簇 bootstrap（簇键排序、PCG64(20260926)、B=10000，与 p2.cluster_boot 同一抽样）；每个候选在「原始与订正量
     都有效」的同一事件集上、同一次抽样中算原始趋势 T_ref、订正趋势 T_c、Δ＝T_c−T_ref；另报以站为簇的 95% 区间作对照。
     多重比较：主检验 9 个（H1a、H1b、H2、H3、H4、H5、H6a、H6b、H6c），族错误率 0.05 → Bonferroni 单检验 α＝0.05/9，
     判读一律用 (α/2, 1−α/2) 百分位区间（99.44%）；95% 区间只作报告。次要／描述量不进族、不做判读或只按事先写定的规则判。
  Y8 MH：RR＝Σ_s O_hi,s·D_lo,s/D_s ÷ Σ_s O_lo,s·D_hi,s/D_s（O＝Σ观测、D＝Σ模型、D_s＝D_hi,s＋D_lo,s），只用高低两档都有事件的层；
     重抽中某层 D_s=0 即该层项为 0；RR≤0 记 NaN（计数）。
  Y9 可评：高、低两档各 n≥30 且 ≥8 个站×季簇（H3 按共同支撑事件计）；否则「不可评」。重抽 NaN 比例 >1% 标注。
  Y10 H1c（次要，95%）：H1a 同集上逐档拟合 κ_k 使 ΣG5(κ)/ΣG1(κ)＝Σobs5/Σobs1（κ 网格 0.02–50 对数等距 81 点，ln κ 上线性插值；
     非单调时取最接近 κ=1 的交点，无交点取 |差| 最小的网格点并计越界），ρ′_k＝Σobs1/ΣG1(κ_k)；T_ρ′＝ln(ρ′_高/ρ′_低)，对照 κ=1 的 T_ρ,κ1。
  Y11 描述：各档有效混合深度 h1＝Σ累计雨量(m)/Σ(−ΔS1/S0)（观测：雨量计与 CMORPH 两种分母；模型：RIM、G κ=1）、每 mm 淡化的
     风速幂指数；ln(ΣC6/ΣG6) 高−低；构成表；A／B 站群各自趋势；雨前 6 h 趋势订正（obs1−pretrend）；小时合成 R(h)；
     雨量计漏接最坏情形（高档 ×1.30、低档 ×1.0）对 H2 的界；已判「支持」的 H2／H3／H5 订正的联合（描述）。

Change Log：
  2026-09-27 初版（p7c-2026-09-27a）。
  2026-09-27 b（修订 1）：首次运行在 V-P4F 核对处中止（任何候选统计量之前）：恰在 1/3 分位上的 1 个事件因
     6 位小数切点被分到中档。改为精确分位切点（Y2）；判据与其余实现不变。
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import p1_events as p1          # noqa: E402
import p2_rim_test as p2        # noqa: E402

VERSION = "p7c-2026-09-27b"
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
P4B_DIR_DEFAULT = _rp.upstream("p4b_dir")  # [repro] 读 p4b 阶段输出目录
NAN = float("nan")

# ---- 事先写定的常数 ----
CUTS_REF = (6.031179, 9.325534)
N_PRIMARY = 9
ALPHA_FAM = 0.05
ALPHA_ONE = ALPHA_FAM / N_PRIMARY
Q_BONF = (100 * ALPHA_ONE / 2, 100 * (1 - ALPHA_ONE / 2))
Q_95 = (2.5, 97.5)
B, SEED = 10000, 20260926
CELL_MIN_N, CELL_MIN_G = 30, 8
NAN_FLAG = 0.01
SHARE_MOST = 0.5            # 解释原始趋势「大部分」的份额门
SHARE_CONFLICT = 0.8        # 「消除冲突」的份额门
H1C_POINT = 0.3             # H1c：|T_ρ′| 点估计上限
KAPPA_GRID_LO, KAPPA_GRID_HI, KAPPA_N = 0.02, 50.0, 81
T_CUR_G = 900.0             # Y4 当前半步的扩散时间
UNDERCATCH_MAX = 1.30
SUB_TOL = 0.05
P4F_REF = {"R_obs": (0.316, 0.26005, 0.14868), "R_pred": (0.17103, 0.24391, 0.26702),
           "rho": (1.84766, 1.06616, 0.55683), "n": (178, 195, 207)}
P4F_TOL = 1e-4
_CUTS = list(CUTS_REF)      # 修订 1：正式运行时改为 646 集 u10 的精确三分位（与 CUTS_REF 差 ≤1e-5）


def _np():
    import numpy
    return numpy


def isnum(x):
    return x is not None and isinstance(x, (int, float)) and math.isfinite(x)


def rnd(x, nd=5):
    return round(float(x), nd) if isnum(x) else None


# ======================================================================== RIM 重新实现（Y3）与守恒 Gaussian（Y4）
def rim_parts(P, U, z, t_cur_depth=p2.T_CURRENT_S):
    """与 p2.rim_factor 同式，返回 (term[N,48], d0[N,48], cterm[N], d0c[N])。"""
    np = _np()
    from numpy.lib.stride_tricks import sliding_window_view as swv
    P = np.asarray(P, float)
    U = np.asarray(U, float)
    L = len(P)
    tc = np.array(p2.TC, float)
    ti = np.array(p2.TI, float)
    Pw = swv(P, 48)[:L - 48]
    Uw = swv(U, 48)[:L - 48]
    irr = Pw / 1000.0 / 3600.0
    kz = p2.KZ_COEF * Uw ** 2
    d0 = p2.d0_interp(Uw, Pw)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        term = p2.C1 * irr * tc / np.sqrt(kz * ti)
        if z:
            term = term * np.exp(-(z ** 2) / (4.0 * kz * ti))
        Pc, Uc = P[48:], U[48:]
        irrc = Pc / 1000.0 / 3600.0
        kzc = p2.KZ_COEF * Uc ** 2
        d0c = p2.d0_interp(Uc, Pc)
        cterm = p2.C2 * irrc * 1800.0 / np.sqrt(kzc * p2.T_CURRENT_S)
        if z:
            cterm = cterm * np.exp(-(z ** 2) / (4.0 * kzc * t_cur_depth))
    return term, d0, cterm, d0c


def rim_F(P, U, z):
    np = _np()
    term, d0, cterm, d0c = rim_parts(P, U, z)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.prod(d0 / (d0 + term), axis=-1) * (d0c / (d0c + cterm))


def rim_F_lin(P, U, z):
    """线性化（去 d0/(d0+x) 饱和）：F_lin＝1−(Σterm/d0＋cterm/d0c)。"""
    np = _np()
    term, d0, cterm, d0c = rim_parts(P, U, z)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 - (np.sum(term / d0, axis=-1) + cterm / d0c)


def cons_G(P, U, z, kappa=1.0):
    """Y4：守恒 Gaussian 的淡水分数 G（与 rim_F 同长度），返回 F_G＝1−G。"""
    np = _np()
    from numpy.lib.stride_tricks import sliding_window_view as swv
    P = np.asarray(P, float)
    U = np.asarray(U, float)
    L = len(P)
    ti = np.array(p2.TI, float)
    Pw = swv(P, 48)[:L - 48] * 0.5 / 1000.0          # 半步雨量（m）
    Uw = swv(U, 48)[:L - 48]
    k = kappa * p2.KZ_COEF * Uw ** 2
    Pc, Uc = P[48:] * 0.5 / 1000.0, U[48:]
    kc = kappa * p2.KZ_COEF * Uc ** 2
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        g = np.sum(Pw / np.sqrt(math.pi * k * ti) * np.exp(-(z ** 2) / (4.0 * k * ti)), axis=-1)
        gc = Pc / np.sqrt(math.pi * kc * T_CUR_G) * np.exp(-(z ** 2) / (4.0 * kc * T_CUR_G))
    return 1.0 - (g + gc)


def hourly(F):
    return F.reshape(-1, 2).mean(axis=1)


def win_delta_h(Fh, win):
    """小时序列（t0−6 …）的窗差：mean[6+a, 6+b) − median[0, 6)。"""
    np = _np()
    pre, post = Fh[:6], Fh[6 + win[0]:6 + win[1]]
    if not (np.all(np.isfinite(pre)) and np.all(np.isfinite(post))):
        return NAN
    return float(np.mean(post) - np.median(pre))


# ======================================================================== 第一步：公式分析（合成，无数据）
def formula_analysis():
    np = _np()
    Us = [2.0, 3.0, 4.0, 4.5, 6.0, 7.8, 8.0, 10.0, 11.6, 14.0]
    cases = {"10mm_h_x1h": [(60, 10.0), (61, 10.0)], "5mm_h_x3h": [(60 + j, 5.0) for j in range(6)],
             "20mm_h_x0.5h": [(60, 20.0)]}
    zs = (0.0, 0.5, 1.0, 5.0)
    out = {"note": "合成事件：84 个半步，t0＝第 60 个半步，常风 U10；给出 0–6 h 窗 ΔF（雨前 6 h 中位为基线）每 mm 雨的量，"
                   "RIM＝p2 同式、RIM_lin＝去饱和、G＝守恒 Gaussian（κ=1）；ratio_G_over_RIM≈1/M。", "cases": {}}
    for cname, rain in cases.items():
        P = np.zeros(84)
        for j, v in rain:
            P[j] = v
        tot_mm = sum(v * 0.5 for _, v in rain)
        rows = []
        for u in Us:
            U = np.full(84, u)
            rec = {"U10": u}
            for z in zs:
                fr = -win_delta_h(hourly(rim_F(P, U, z)), (0, 6)) / tot_mm
                fl = -win_delta_h(hourly(rim_F_lin(P, U, z)), (0, 6)) / tot_mm
                fg = -win_delta_h(hourly(cons_G(P, U, z)), (0, 6)) / tot_mm
                rec[f"z{z}"] = {"RIM_per_mm": rnd(fr, 7), "RIMlin_per_mm": rnd(fl, 7), "G_per_mm": rnd(fg, 7),
                                "ratio_G_over_RIM": rnd(fg / fr, 4) if fr else None}
            r1_03 = -win_delta_h(hourly(rim_F(P, U, 1.0)), (0, 3))
            r1_612 = -win_delta_h(hourly(rim_F(P, U, 1.0)), (6, 12))
            g1_03 = -win_delta_h(hourly(cons_G(P, U, 1.0)), (0, 3))
            g1_612 = -win_delta_h(hourly(cons_G(P, U, 1.0)), (6, 12))
            rec["RIM_1m_over_0m"] = rnd(rec["z1.0"]["RIM_per_mm"] / rec["z0.0"]["RIM_per_mm"], 4)
            rec["RIM_5m_over_1m"] = rnd(rec["z5.0"]["RIM_per_mm"] / rec["z1.0"]["RIM_per_mm"], 4)
            rec["G_5m_over_1m"] = rnd(rec["z5.0"]["G_per_mm"] / rec["z1.0"]["G_per_mm"], 4)
            rec["RIM_1m_decay_612_over_03"] = rnd(r1_612 / r1_03, 4)
            rec["G_1m_decay_612_over_03"] = rnd(g1_612 / g1_03, 4)
            rec["d0"] = rnd(float(p2.d0_interp(np.array([u]), np.array([max(v for _, v in rain)]))[0]), 3)
            rows.append(rec)
        # 每 mm 淡化对 U 的对数斜率（4.5→11.6，即 P4f 低/高档中位风速）
        def at(u, key, z):
            return next(r for r in rows if r["U10"] == u)[f"z{z}"][key]
        slopes = {}
        for z in zs:
            for key in ("RIM_per_mm", "G_per_mm"):
                a, b_ = at(4.5, key, z), at(11.6, key, z)
                slopes[f"{key}_z{z}"] = rnd(math.log(b_ / a) / math.log(11.6 / 4.5), 3) if a and b_ and a > 0 and b_ > 0 else None
        out["cases"][cname] = {"rain_mm": tot_mm, "rows": rows, "loglog_slope_U4.5_to_11.6": slopes}
    return out


# ======================================================================== 统计核心（Y7、Y8）
class Boot:
    """一个事件集上的站×季（或站）整簇抽样；counts[B,K]＝各簇在每次重抽中的次数（与 p2.cluster_boot 同一抽样）。"""

    def __init__(self, keys):
        np = _np()
        uk = sorted(set(keys))
        kid = {k: j for j, k in enumerate(uk)}
        self.idx = np.array([kid[k] for k in keys], int)
        self.K = len(uk)
        rng = np.random.Generator(np.random.PCG64(SEED))
        draw = rng.integers(0, self.K, size=(B, self.K))
        self.counts = np.zeros((B, self.K))
        for j in range(self.K):
            self.counts[:, j] = (draw == j).sum(axis=1)

    def sums(self, v):
        """v：[n] 或 [n,m] → (点和 [m], 重抽和 [B,m])。"""
        np = _np()
        v = np.asarray(v, float)
        one = v.ndim == 1
        if one:
            v = v[:, None]
        cs = np.zeros((self.K, v.shape[1]))
        np.add.at(cs, self.idx, v)
        pt = cs.sum(axis=0)
        bs = self.counts @ cs
        return (pt[0], bs[:, 0]) if one else (pt, bs)


def pct(a, q):
    np = _np()
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if len(a) < 100:
        return [None, None]
    return [rnd(np.percentile(a, q[0])), rnd(np.percentile(a, q[1]))]


def nan_frac(a):
    np = _np()
    a = np.asarray(a, float)
    return float(np.mean(~np.isfinite(a)))


def safe_ln(x):
    np = _np()
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.asarray(x, float)
        return np.where(x > 0, np.log(np.where(x > 0, x, 1.0)), np.nan)


def trend_ratio(bt, lab, num, den):
    """ln(Σnum/Σden)_高 − ln(...)_低：返回 (点, 重抽数组)。"""
    np = _np()
    cols = []
    for k in (0, 2):
        s = lab == k
        cols += [np.where(s, num, 0.0), np.where(s, den, 0.0)]
    pt, bs = bt.sums(np.stack(cols, axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        tp = safe_ln(pt[2] / pt[3]) - safe_ln(pt[0] / pt[1])
        tb = safe_ln(bs[:, 2] / bs[:, 3]) - safe_ln(bs[:, 0] / bs[:, 1])
    return float(tp), tb


def tercile_ratios(bt, lab, num, den):
    np = _np()
    out = []
    for k in range(3):
        s = lab == k
        pt, bs = bt.sums(np.stack([np.where(s, num, 0.0), np.where(s, den, 0.0)], axis=1))
        with np.errstate(divide="ignore", invalid="ignore"):
            out.append({"k": k, "n": int(s.sum()), "point": rnd(pt[0] / pt[1]), "ci95": pct(bs[:, 0] / bs[:, 1], Q_95)})
    return out


def mh_trend(bt, lab, strata, num, den):
    """Y8：层内高/低档 MH 率比的对数；只用两档都有事件的层。返回 (点, 重抽数组, 所用层, 事件掩码)。"""
    np = _np()
    strata = np.asarray(strata)
    use = [s for s in sorted(set(strata.tolist())) if ((strata == s) & (lab == 0)).any() and ((strata == s) & (lab == 2)).any()]
    cols = []
    for s in use:
        for k in (2, 0):
            m = (strata == s) & (lab == k)
            cols += [np.where(m, num, 0.0), np.where(m, den, 0.0)]
    pt, bs = bt.sums(np.stack(cols, axis=1))

    def rr(a):
        a = np.atleast_2d(a)
        n_, d_ = 0.0, 0.0
        for j in range(len(use)):
            ohi, dhi, olo, dlo = a[:, 4 * j], a[:, 4 * j + 1], a[:, 4 * j + 2], a[:, 4 * j + 3]
            ds = dhi + dlo
            with np.errstate(divide="ignore", invalid="ignore"):
                n_ = n_ + np.where(ds != 0, ohi * dlo / ds, 0.0)
                d_ = d_ + np.where(ds != 0, olo * dhi / ds, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            return safe_ln(n_ / d_)
    mask = np.isin(strata, use) & ((lab == 0) | (lab == 2))
    return float(rr(pt[None, :])[0]), rr(bs), use, mask


def cell_ok(rows, lab, clusters):
    """Y9：高、低两档的 n 与簇数。"""
    res = {}
    for k in (0, 2):
        sel = [c for c, l_ in zip(clusters, lab) if l_ == k]
        res[k] = {"n": len(sel), "clusters": len(set(sel))}
    ok = all(res[k]["n"] >= CELL_MIN_N and res[k]["clusters"] >= CELL_MIN_G for k in (0, 2))
    return ok, {"low": res[0], "high": res[2]}


# ======================================================================== 判读
def judge_artifact(delta_bonf, delta_pt, t_ref, evaluable):
    """H2–H6b：Δ＝T_c−T_ref。"""
    if not evaluable:
        return "不可评"
    lo = delta_bonf[0]
    if lo is None:
        return "不可评"
    if lo <= 0:
        return "不支持"
    share = delta_pt / abs(t_ref) if t_ref else float("inf")
    return "支持（解释原始趋势大部分）" if share >= SHARE_MOST else "支持（解释原始趋势一部分）"


def judge_conflict(delta_bonf, share_c):
    if delta_bonf[0] is None or share_c is None:
        return "不可评"
    if delta_bonf[0] <= 0:
        return "未缓解冲突"
    return "足以消除冲突" if share_c >= SHARE_CONFLICT else "部分缓解冲突"


def judge_direction(ci_bonf, sign_pred, evaluable):
    """H1a（预言 >0）、H1b（预言 <0）。"""
    if not evaluable or ci_bonf[0] is None:
        return "不可评"
    lo, hi = ci_bonf
    if sign_pred > 0:
        return "支持" if lo > 0 else ("与预言相反" if hi < 0 else "未检出")
    return "支持" if hi < 0 else ("与预言相反" if lo > 0 else "未检出")


def judge_h1(a, b_):
    vs = [a, b_]
    if "不可评" in vs and not any(v == "支持" for v in vs):
        return "不可评"
    n_sup = sum(v == "支持" for v in vs)
    if n_sup == 2:
        return "支持"
    if n_sup == 1:
        return "部分支持" if "与预言相反" not in vs else "证据相互矛盾"
    return "不支持（与预言相反）" if "与预言相反" in vs else "不支持（未检出）"


def judge_h1c(t_pt, t_ci95, d_ci95):
    if t_ci95[0] is None or d_ci95[0] is None:
        return "不可评"
    inside = (t_ci95[0] <= 0 <= t_ci95[1]) or (-H1C_POINT <= t_ci95[0] and t_ci95[1] <= H1C_POINT)
    ok = inside and abs(t_pt) <= H1C_POINT and d_ci95[0] > 0
    return "形状拟合的守恒模型能定量解释冲突" if ok else "形状拟合的守恒模型不能定量解释冲突"


def overall_map(C):
    """对矛盾 C 的总映射。"""
    full = [k for k, v in C.items() if k.startswith(("H2", "H3", "H4", "H5", "H6")) and v.get("conflict") == "足以消除冲突"]
    part = [k for k, v in C.items() if k.startswith(("H2", "H3", "H4", "H5", "H6")) and v.get("conflict") == "部分缓解冲突"]
    h1 = C.get("H1", {}).get("verdict")
    h1c = C.get("H1c", {}).get("verdict")
    phys = h1 == "支持" and h1c == "形状拟合的守恒模型能定量解释冲突"
    if full:
        return f"冲突由 {full} 消除（测量/统计构成类解释）"
    if phys:
        return "冲突可由风速相关的混合深度错配（守恒＋形状）定量解释"
    if part or h1 in ("支持", "部分支持"):
        return f"部分解释：缓解冲突的候选 {part}；H1＝{h1}，H1c＝{h1c}"
    return "解释不了：没有候选达到写死的判据"


# ======================================================================== 逐事件（Y6）
def event_p7(stn, e, row, rec, pixc):
    np = _np()
    ser, i, H = stn["ser"], e["i"], e["hour"]
    z1 = row["z1"]
    s0 = row.get("s0_1")
    x = {k: row.get(k) for k in ("station", "group", "season", "onset_utc", "z1", "rain24", "u10", "obs1", "obs5", "rim1",
                                 "rim0", "rim5", "g6", "c6", "pcum_cmorph_model", "pcum_gauge_obs", "model_ok")}
    x["s0_1"] = s0
    x["sub"] = int(abs(z1 - 1.0) > SUB_TOL)
    for w in ((0, 3), (6, 12)):
        x[f"obs1_w{w[0]}{w[1]}"] = p2.obs_delta(ser["s1"], i, w)
    cl = [c for c in e.get("ctrl_ds1", []) if isnum(c)]
    x["ctrl_mean1"] = float(np.mean(cl)) if cl else NAN
    x["n_ctrl_valid"] = len(cl)
    a = np.asarray(ser["s1"][max(0, i - 12):max(0, i - 6)], float)
    b_ = np.asarray(ser["s1"][max(0, i - 6):i], float)
    a, b_ = a[np.isfinite(a)], b_[np.isfinite(b_)]
    x["pretrend1"] = float(np.median(b_) - np.mean(a)) if len(a) and len(b_) else NAN
    obs_h = []
    pre_med = p2.pre_median(ser["s1"], i)
    for h in range(12):
        v = ser["s1"][i + h] if i + h < len(ser["s1"]) else NAN
        obs_h.append(float(v - pre_med) if isnum(pre_med) and math.isfinite(v) else NAN)
    x["_obs_h"] = obs_h
    P, U, _ = p2.event_forcing(stn, H, "cmorph", pixc)
    x["_mod_h"] = [NAN] * 12
    for k in ("rim1_chk", "rim1_w03", "rim1_w612", "rim1_lin", "G1", "G5", "G0", "G1_w03", "G1_w612"):
        x[k] = NAN
    x["_G1k"], x["_G5k"] = None, None
    if P is not None and isnum(s0):
        Fh = hourly(rim_F(P, U, z1))
        x["rim1_chk"] = s0 * win_delta_h(Fh, (0, 6))
        x["rim1_w03"] = s0 * win_delta_h(Fh, (0, 3))
        x["rim1_w612"] = s0 * win_delta_h(Fh, (6, 12))
        x["_mod_h"] = [float(s0 * (Fh[6 + h] - np.median(Fh[:6]))) for h in range(12)]
        x["rim1_lin"] = s0 * win_delta_h(hourly(rim_F_lin(P, U, z1)), (0, 6))
        Gh = hourly(cons_G(P, U, z1))
        x["G1"] = s0 * win_delta_h(Gh, (0, 6))
        x["G1_w03"] = s0 * win_delta_h(Gh, (0, 3))
        x["G1_w612"] = s0 * win_delta_h(Gh, (6, 12))
        x["G0"] = s0 * win_delta_h(hourly(cons_G(P, U, 0.0)), (0, 6))
        s05 = p2.pre_median(ser["s5"], i)
        x["G5"] = s05 * win_delta_h(hourly(cons_G(P, U, 5.0)), (0, 6)) if isnum(s05) else NAN
        if isnum(x["obs5"]) and isnum(x["obs1"]) and isnum(s05):
            kg = kappa_grid()
            x["_G1k"] = [s0 * win_delta_h(hourly(cons_G(P, U, z1, kk)), (0, 6)) for kk in kg]
            x["_G5k"] = [s05 * win_delta_h(hourly(cons_G(P, U, 5.0, kk)), (0, 6)) for kk in kg]
    Pg, Ug, _ = p2.event_forcing(stn, H, "gauge", pixc)
    x["rimG1"] = s0 * win_delta_h(hourly(rim_F(Pg, Ug, z1)), (0, 6)) if (Pg is not None and isnum(s0)) else NAN
    return x


def kappa_solve(rm, tg, lk):
    """Y10：模型比值曲线 rm(ln κ) 与目标 tg 的交点；多个交点取 |ln κ| 最小者（最接近 κ=1）；无交点取 |rm−tg| 最小的
    网格点并记越界。返回 (ln κ, 越界旗标)。"""
    np = _np()
    d = rm - tg
    sol = []
    for j in range(len(d) - 1):
        if d[j] == 0:
            sol.append(lk[j])
        elif d[j] * d[j + 1] < 0:
            sol.append(lk[j] + (lk[j + 1] - lk[j]) * d[j] / (d[j] - d[j + 1]))
    if d[-1] == 0:
        sol.append(lk[-1])
    if sol:
        return float(min(sol, key=abs)), False
    return float(lk[int(np.argmin(np.abs(d)))]), True


def kappa_grid():
    np = _np()
    return np.exp(np.linspace(math.log(KAPPA_GRID_LO), math.log(KAPPA_GRID_HI), KAPPA_N))


# ======================================================================== 分析
def fin(*vals):
    return all(isnum(v) for v in vals)


def subset(rows, keys, extra=None):
    return [r for r in rows if fin(*[r.get(k) for k in keys]) and (extra is None or extra(r))]


def arr(rows, k):
    np = _np()
    return np.array([r[k] for r in rows], float)


def labels(rows):
    np = _np()
    u = arr(rows, "u10")
    return np.where(u <= _CUTS[0], 0, np.where(u <= _CUTS[1], 1, 2))


def ss_keys(rows):
    return [(r["station"], r["season"]) for r in rows]


def artifact_block(name, rows, num_ref, den_ref, num_c, den_c, mh_strata=None, conflict_mode="R", keys_needed=()):
    """H2–H6：同一事件集、同一次抽样的 T_ref、T_c、Δ；冲突份额在 M 可得子集上算。"""
    np = _np()
    need = [num_ref, den_ref, num_c, den_c] + list(keys_needed)
    X = subset(rows, need)
    lab = labels(X)
    res = {"name": name, "n_set": len(X)}
    if mh_strata is not None:
        strata = np.array([mh_strata(r) for r in X], dtype=object).astype(str)
        _, _, use, mask = mh_trend(Boot(ss_keys(X)), lab, strata, arr(X, num_c), arr(X, den_c))
        X = [r for r, m in zip(X, mask) if m or False]
        lab = labels(X)
        strata = np.array([mh_strata(r) for r in X], dtype=object).astype(str)
        res["mh_strata_used"] = use
        res["n_common_support"] = len(X)
    ok, cells = cell_ok(X, lab, ss_keys(X))
    res["cells"] = cells
    res["evaluable"] = ok
    for tag, keys in (("ss", ss_keys(X)), ("st", [r["station"] for r in X])):
        bt = Boot(keys)
        tr, trb = trend_ratio(bt, lab, arr(X, num_ref), arr(X, den_ref))
        if mh_strata is not None:
            tc, tcb, _, _ = mh_trend(bt, lab, strata, arr(X, num_c), arr(X, den_c))
        else:
            tc, tcb = trend_ratio(bt, lab, arr(X, num_c), arr(X, den_c))
        d, db = tc - tr, tcb - trb
        if tag == "ss":
            res.update({"T_ref": rnd(tr), "T_ref_ci95": pct(trb, Q_95), "T_c": rnd(tc), "T_c_ci95": pct(tcb, Q_95),
                        "T_c_ci_bonf": pct(tcb, Q_BONF), "delta": rnd(d), "delta_ci95": pct(db, Q_95),
                        "delta_ci_bonf": pct(db, Q_BONF), "nan_frac_delta": rnd(nan_frac(db), 4),
                        "share_of_raw_trend": rnd(d / abs(tr)) if tr else None})
            res["nan_flag"] = nan_frac(db) > NAN_FLAG
        else:
            res["delta_ci95_station_cluster"] = pct(db, Q_95)
    res["terciles_ref"] = tercile_ratios(Boot(ss_keys(X)), lab, arr(X, num_ref), arr(X, den_ref))
    if mh_strata is None:
        res["terciles_c"] = tercile_ratios(Boot(ss_keys(X)), lab, arr(X, num_c), arr(X, den_c))
    # 冲突份额（Y5：M 可得子集）
    XM = [r for r in X if isnum(r.get("M"))]
    if len(XM) >= 2 * CELL_MIN_N:
        labm = labels(XM)
        bt = Boot(ss_keys(XM))
        tr, trb = trend_ratio(bt, labm, arr(XM, num_ref), arr(XM, den_ref))
        if mh_strata is not None:
            sm = np.array([mh_strata(r) for r in XM], dtype=object).astype(str)
            tc, tcb, _, _ = mh_trend(bt, labm, sm, arr(XM, num_c), arr(XM, den_c))
        else:
            tc, tcb = trend_ratio(bt, labm, arr(XM, num_c), arr(XM, den_c))
        m1 = arr(XM, "rim1")
        tp, tpb = trend_ratio(bt, labm, m1 / arr(XM, "M"), m1)
        if conflict_mode == "rho":                     # H6c：T_ref 本身就是 P4f 的 ρ 趋势，冲突＝|T_ref|
            conf, confb = abs(tr), np.abs(trb)
        else:
            conf, confb = tp - tr, tpb - trb
        d, db = tc - tr, tcb - trb
        with np.errstate(divide="ignore", invalid="ignore"):
            shb = db / confb
        res["conflict"] = {"n_M": len(XM), "T_pred": rnd(tp), "T_ref_M": rnd(tr), "T_c_M": rnd(tc),
                           "conflict_size": rnd(conf), "delta_M": rnd(d), "delta_M_ci_bonf": pct(db, Q_BONF),
                           "share_of_conflict": rnd(d / conf) if conf else None, "share_ci95": pct(shb, Q_95),
                           "T_c_minus_T_pred": rnd(tc - tp), "T_c_minus_T_pred_ci95": pct(tcb - tpb, Q_95)}
    return res, X


def rows_A(rows, g):
    return [r for r in rows if r["group"] == g]


def analyze(rows, log=None, check_p4f=True):
    np = _np()
    S = {}
    base = subset(rows, ("u10", "obs1", "rim1"), lambda r: r.get("model_ok"))
    S["n_base"] = len(base)
    _CUTS[:] = list(CUTS_REF)
    if check_p4f:                                   # 修订 1：切点用精确分位数（第一轮 D-U 即如此），JSON 的 6 位小数只作核对
        q = [float(v) for v in np.quantile(arr(base, "u10"), [1 / 3, 2 / 3])]
        if max(abs(q[0] - CUTS_REF[0]), abs(q[1] - CUTS_REF[1])) > 1e-5:
            raise RuntimeError(f"精确切点与第一轮不一致 {q}")
        _CUTS[:] = q
    S["cuts_used"] = list(_CUTS)
    labb = labels(base)
    bt = Boot(ss_keys(base))
    S["V0_DU"] = tercile_ratios(bt, labb, arr(base, "obs1"), arr(base, "rim1"))
    tb_pt, tb_b = trend_ratio(bt, labb, arr(base, "obs1"), arr(base, "rim1"))
    S["T_raw_646"] = {"point": rnd(tb_pt), "ci95": pct(tb_b, Q_95), "ci_bonf": pct(tb_b, Q_BONF)}
    # V-P4F
    XM = subset(base, ("M",))
    lab = labels(XM)
    vp = {"n": [int((lab == k).sum()) for k in range(3)]}
    for k in range(3):
        s = lab == k
        o, m, M = arr(XM, "obs1")[s], arr(XM, "rim1")[s], arr(XM, "M")[s]
        vp.setdefault("R_obs", []).append(float(o.sum() / m.sum()))
        vp.setdefault("R_pred", []).append(float((m / M).sum() / m.sum()))
        vp.setdefault("rho", []).append(float(o.sum() / (m / M).sum()))
    dmax = max(abs(vp[q][k] - P4F_REF[q][k]) for q in ("R_obs", "R_pred", "rho") for k in range(3))
    S["V_P4F"] = {"this": {q: [rnd(v) for v in vp[q]] for q in ("R_obs", "R_pred", "rho")}, "n": vp["n"],
                  "ref": P4F_REF, "max_abs_diff": dmax, "pass": bool(dmax <= P4F_TOL and tuple(vp["n"]) == P4F_REF["n"])}
    if check_p4f and not S["V_P4F"]["pass"]:
        raise RuntimeError(f"V-P4F 不过：{S['V_P4F']}")

    C = {}
    # ---- H1a 形状（参照＝守恒 Gaussian κ=1；RIM 参照只作描述）
    X = subset(base, ("obs5", "rim5", "G1", "G5"))
    labx = labels(X)
    ok, cells = cell_ok(X, labx, ss_keys(X))
    blk = {"n_set": len(X), "cells": cells, "evaluable": ok}
    for tag, keys in (("ss", ss_keys(X)), ("st", [r["station"] for r in X])):
        b1 = Boot(keys)
        to, tob = trend_ratio(b1, labx, arr(X, "obs5"), arr(X, "obs1"))
        tg, tgb = trend_ratio(b1, labx, arr(X, "G5"), arr(X, "G1"))
        tm, tmb = trend_ratio(b1, labx, arr(X, "rim5"), arr(X, "rim1"))
        if tag == "ss":
            blk.update({"T_shape_obs": rnd(to), "T_shape_obs_ci95": pct(tob, Q_95), "T_shape_G": rnd(tg),
                        "T_shape_G_ci95": pct(tgb, Q_95), "T_shape_rim": rnd(tm), "T_shape_rim_ci95": pct(tmb, Q_95),
                        "A": rnd(to - tg), "A_ci95": pct(tob - tgb, Q_95), "A_ci_bonf": pct(tob - tgb, Q_BONF),
                        "A_vs_rim_descr": rnd(to - tm), "A_vs_rim_descr_ci95": pct(tob - tmb, Q_95),
                        "nan_frac": rnd(nan_frac(tob - tgb), 4)})
        else:
            blk["A_ci95_station_cluster"] = pct(tob - tgb, Q_95)
    b1 = Boot(ss_keys(X))
    blk["terciles_obs5_over_obs1"] = tercile_ratios(b1, labx, arr(X, "obs5"), arr(X, "obs1"))
    blk["terciles_G5_over_G1"] = tercile_ratios(b1, labx, arr(X, "G5"), arr(X, "G1"))
    blk["terciles_rim5_over_rim1"] = tercile_ratios(b1, labx, arr(X, "rim5"), arr(X, "rim1"))
    blk["verdict"] = judge_direction(blk["A_ci_bonf"], +1, ok)
    C["H1a"] = blk
    X5 = X
    # ---- H1b 时间
    X = subset(base, ("obs1_w03", "obs1_w612", "rim1_w03", "rim1_w612", "G1_w03", "G1_w612"))
    labx = labels(X)
    ok, cells = cell_ok(X, labx, ss_keys(X))
    blk = {"n_set": len(X), "cells": cells, "evaluable": ok}
    for tag, keys in (("ss", ss_keys(X)), ("st", [r["station"] for r in X])):
        b1 = Boot(keys)
        t1, t1b = trend_ratio(b1, labx, arr(X, "obs1_w612"), arr(X, "G1_w612"))
        t0, t0b = trend_ratio(b1, labx, arr(X, "obs1_w03"), arr(X, "G1_w03"))
        r1, r1b = trend_ratio(b1, labx, arr(X, "obs1_w612"), arr(X, "rim1_w612"))
        r0, r0b = trend_ratio(b1, labx, arr(X, "obs1_w03"), arr(X, "rim1_w03"))
        if tag == "ss":
            blk.update({"T_612_G": rnd(t1), "T_612_G_ci95": pct(t1b, Q_95), "T_03_G": rnd(t0), "T_03_G_ci95": pct(t0b, Q_95),
                        "Bstat": rnd(t1 - t0), "B_ci95": pct(t1b - t0b, Q_95), "B_ci_bonf": pct(t1b - t0b, Q_BONF),
                        "B_vs_rim_descr": rnd(r1 - r0), "B_vs_rim_descr_ci95": pct(r1b - r0b, Q_95),
                        "nan_frac": rnd(nan_frac(t1b - t0b), 4)})
        else:
            blk["B_ci95_station_cluster"] = pct(t1b - t0b, Q_95)
    b1 = Boot(ss_keys(X))
    blk["terciles_R_03"] = tercile_ratios(b1, labx, arr(X, "obs1_w03"), arr(X, "rim1_w03"))
    blk["terciles_R_612"] = tercile_ratios(b1, labx, arr(X, "obs1_w612"), arr(X, "rim1_w612"))
    blk["terciles_rhoG_03"] = tercile_ratios(b1, labx, arr(X, "obs1_w03"), arr(X, "G1_w03"))
    blk["terciles_rhoG_612"] = tercile_ratios(b1, labx, arr(X, "obs1_w612"), arr(X, "G1_w612"))
    blk["verdict"] = judge_direction(blk["B_ci_bonf"], -1, ok)
    C["H1b"] = blk
    C["H1"] = {"verdict": judge_h1(C["H1a"]["verdict"], C["H1b"]["verdict"])}
    # ---- H1c（次要）
    C["H1c"] = h1c_block([r for r in X5 if r.get("_G1k") is not None])
    # ---- H2–H6
    h2, X2 = artifact_block("H2 雨量计驱动", base, "obs1", "rim1", "obs1", "rimG1")
    h2["verdict"] = judge_artifact(h2["delta_ci_bonf"], h2["delta"], h2["T_ref"], h2["evaluable"])
    h2["conflict_verdict"] = judge_conflict(h2["conflict"]["delta_M_ci_bonf"], h2["conflict"]["share_of_conflict"]) \
        if "conflict" in h2 else "不可评"
    uc = h2["T_c"] - math.log(UNDERCATCH_MAX) if isnum(h2["T_c"]) else None
    h2["undercatch_worst_case_T_G"] = rnd(uc)
    C["H2"] = h2
    h3, _ = artifact_block("H3 站层 MH", base, "obs1", "rim1", "obs1", "rim1", mh_strata=lambda r: r["station"])
    C["H3"] = h3
    cuts_r = [float(q) for q in np.quantile(arr(base, "rain24"), [1 / 3, 2 / 3])]

    def comp(r):
        t = 0 if r["rain24"] <= cuts_r[0] else (1 if r["rain24"] <= cuts_r[1] else 2)
        return f"rain{t}_sub{r['sub']}"
    h4, _ = artifact_block("H4 构成层 MH", base, "obs1", "rim1", "obs1", "rim1", mh_strata=comp)
    h4["rain24_cuts"] = [rnd(c, 3) for c in cuts_r]
    C["H4"] = h4
    for r in base:
        r["obs1_cc"] = r["obs1"] - r["ctrl_mean1"] if fin(r["obs1"], r.get("ctrl_mean1")) else NAN
    C["H5"], _ = artifact_block("H5 对照窗订正", base, "obs1", "rim1", "obs1_cc", "rim1")
    C["H6a"], _ = artifact_block("H6a 分母＝海面 RIM", base, "obs1", "rim1", "obs1", "rim0")
    C["H6b"], _ = artifact_block("H6b 分母＝线性化 RIM", base, "obs1", "rim1", "obs1", "rim1_lin")
    # H6c：ρ 趋势（P4f 1/M）→ 守恒 Gaussian 分母
    for r in base:
        r["rim1_over_M"] = r["rim1"] / r["M"] if fin(r["rim1"], r.get("M")) else NAN
    h6c, _ = artifact_block("H6c 守恒 Gaussian vs 1/M", base, "obs1", "rim1_over_M", "obs1", "G1", conflict_mode="rho")
    C["H6c"] = h6c
    for key in ("H3", "H4", "H5", "H6a", "H6b", "H6c"):
        b_ = C[key]
        b_["verdict"] = judge_artifact(b_["delta_ci_bonf"], b_["delta"], b_["T_ref"], b_["evaluable"])
        b_["conflict_verdict"] = judge_conflict(b_["conflict"]["delta_M_ci_bonf"], b_["conflict"]["share_of_conflict"]) \
            if "conflict" in b_ else "不可评"
    if isnum(C["H6c"].get("delta")) and C["H6c"]["evaluable"]:   # H6c 的份额以冲突（ρ 趋势）计
        c6 = C["H6c"]["conflict"]
        share = c6.get("share_of_conflict")
        lo = C["H6c"]["delta_ci_bonf"][0]
        C["H6c"]["verdict"] = "不可评" if lo is None else ("不支持" if lo <= 0 else (
            "支持（解释冲突大部分）" if (share or 0) >= SHARE_MOST else "支持（解释冲突一部分）"))
    for key in ("H2", "H3", "H4", "H5", "H6a", "H6b", "H6c"):
        C[key]["conflict"] = C[key].get("conflict", {})
        C[key]["conflict"]["verdict"] = C[key].pop("conflict_verdict", "不可评")
    S["candidates"] = C
    S["overall"] = overall_map({k: {"conflict": v.get("conflict", {}).get("verdict"), "verdict": v.get("verdict")}
                                for k, v in C.items()})
    S["descriptive"] = descriptive(base, C)
    return S


def h1c_block(X):
    """Y10：逐档拟合 κ（同一次抽样），ρ′ 与 κ=1 对照。"""
    np = _np()
    res = {"n_set": len(X)}
    if len(X) < 2 * CELL_MIN_N:
        res["verdict"] = "不可评"
        return res
    lab = labels(X)
    ok, cells = cell_ok(X, lab, ss_keys(X))
    res["cells"] = cells
    kg = kappa_grid()
    lk = np.log(kg)
    G1k = np.array([r["_G1k"] for r in X], float)
    G5k = np.array([r["_G5k"] for r in X], float)
    good = np.all(np.isfinite(G1k), axis=1) & np.all(np.isfinite(G5k), axis=1)
    X = [r for r, g in zip(X, good) if g]
    G1k, G5k, lab = G1k[good], G5k[good], lab[good]
    bt = Boot(ss_keys(X))
    o1, o5 = arr(X, "obs1"), arr(X, "obs5")
    i1 = G1k[:, int(np.argmin(np.abs(kg - 1.0)))]                    # κ 网格中点＝1（对称对数网格）
    per = {}
    for k in (0, 1, 2):
        s = lab == k
        cols = np.concatenate([np.where(s[:, None], G1k, 0.0), np.where(s[:, None], G5k, 0.0),
                               np.stack([np.where(s, o1, 0.0), np.where(s, o5, 0.0), np.where(s, i1, 0.0)], axis=1)], axis=1)
        pt, bs = bt.sums(cols)
        per[k] = (pt, bs)

    def fit(sums):
        """sums：[..., 2*KAPPA_N+3] → (κ, ρ′, ρ_κ1, 越界旗标)。"""
        s = np.atleast_2d(sums)
        g1, g5 = s[:, :KAPPA_N], s[:, KAPPA_N:2 * KAPPA_N]
        ob1, ob5, gk1 = s[:, 2 * KAPPA_N], s[:, 2 * KAPPA_N + 1], s[:, 2 * KAPPA_N + 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_m = g5 / g1
            target = ob5 / ob1
        kap = np.full(len(s), np.nan)
        rho_p = np.full(len(s), np.nan)
        edge = np.zeros(len(s), bool)
        for j in range(len(s)):
            rm, tg = ratio_m[j], target[j]
            if not (np.all(np.isfinite(rm)) and np.isfinite(tg)):
                continue
            lkj, at_edge = kappa_solve(rm, tg, lk)
            edge[j] = at_edge
            kap[j] = math.exp(lkj)
            g1j = float(np.interp(lkj, lk, g1[j]))
            rho_p[j] = ob1[j] / g1j if g1j else np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            rho1 = ob1 / gk1
        return kap, rho_p, rho1, edge
    fits = {}
    for k in (0, 1, 2):
        pt, bs = per[k]
        fp = fit(pt[None, :])
        fb = fit(bs)
        fits[k] = (fp, fb)
        res[f"tercile{k}"] = {"n": int((lab == k).sum()), "kappa": rnd(fp[0][0], 4), "kappa_ci95": pct(fb[0], Q_95),
                              "kappa_at_edge": bool(fp[3][0]), "boot_edge_frac": rnd(float(np.mean(fb[3])), 4),
                              "rho_prime": rnd(fp[1][0]), "rho_prime_ci95": pct(fb[1], Q_95),
                              "rho_kappa1": rnd(fp[2][0]), "rho_kappa1_ci95": pct(fb[2], Q_95)}
    tp = safe_ln(fits[2][0][1][0] / fits[0][0][1][0])
    tpb = safe_ln(fits[2][1][1] / fits[0][1][1])
    t1 = safe_ln(fits[2][0][2][0] / fits[0][0][2][0])
    t1b = safe_ln(fits[2][1][2] / fits[0][1][2])
    res.update({"evaluable": ok, "T_rho_prime": rnd(float(tp)), "T_rho_prime_ci95": pct(tpb, Q_95),
                "T_rho_kappa1": rnd(float(t1)), "T_rho_kappa1_ci95": pct(t1b, Q_95),
                "delta": rnd(float(tp - t1)), "delta_ci95": pct(tpb - t1b, Q_95), "nan_frac": rnd(nan_frac(tpb - t1b), 4)})
    res["verdict"] = judge_h1c(float(tp), res["T_rho_prime_ci95"], res["delta_ci95"]) if ok else "不可评"
    return res


def descriptive(base, C):
    np = _np()
    D = {}
    lab = labels(base)
    # D1 有效混合深度 h1 与每 mm 淡化幂指数
    X = subset(base, ("s0_1", "pcum_cmorph_model", "G1"))
    labx = labels(X)
    bt = Boot(ss_keys(X))
    f_obs = -arr(X, "obs1") / arr(X, "s0_1")
    f_rim = -arr(X, "rim1") / arr(X, "s0_1")
    f_g = -arr(X, "G1") / arr(X, "s0_1")
    pc = arr(X, "pcum_cmorph_model") / 1000.0
    d1 = {"n": len(X), "h1_obs_cmorph_m": tercile_ratios(bt, labx, pc, f_obs), "h1_rim_m": tercile_ratios(bt, labx, pc, f_rim),
          "h1_G_kappa1_m": tercile_ratios(bt, labx, pc, f_g)}
    Xg = subset(X, ("pcum_gauge_obs",))
    d1["h1_obs_gauge_m"] = tercile_ratios(Boot(ss_keys(Xg)), labels(Xg), arr(Xg, "pcum_gauge_obs") / 1000.0,
                                          -arr(Xg, "obs1") / arr(Xg, "s0_1"))
    umed = [float(np.median(arr(X, "u10")[labx == k])) for k in range(3)]
    d1["u10_median"] = [rnd(u, 3) for u in umed]
    lu = math.log(umed[2] / umed[0])
    for nm in ("h1_obs_cmorph_m", "h1_rim_m", "h1_G_kappa1_m", "h1_obs_gauge_m"):
        a, b_ = d1[nm][0]["point"], d1[nm][2]["point"]
        d1[f"slope_per_mm_vs_U_{nm}"] = rnd(-math.log(b_ / a) / lu, 3) if a and b_ and a > 0 and b_ > 0 else None
    D["D1_mixing_depth"] = d1
    # D2 λ
    X = subset(base, ("g6", "c6"))
    labx = labels(X)
    bt = Boot(ss_keys(X))
    tt, ttb = trend_ratio(bt, labx, arr(X, "c6"), arr(X, "g6"))
    D["D2_cmorph_over_gauge"] = {"n": len(X), "terciles": tercile_ratios(bt, labx, arr(X, "c6"), arr(X, "g6")),
                                 "ln_trend": rnd(tt), "ln_trend_ci95": pct(ttb, Q_95)}
    # D3 构成
    comp = []
    for k in range(3):
        s = [r for r, l_ in zip(base, lab) if l_ == k]
        st = {}
        for r in s:
            st[r["station"]] = st.get(r["station"], 0) + 1
        gfrac = [r["g6"] / r["rain24"] for r in s if fin(r.get("g6"), r.get("rain24")) and r["rain24"] > 0]
        comp.append({"k": k, "n": len(s), "stations": dict(sorted(st.items(), key=lambda t: -t[1])),
                     "frac_group_A": rnd(np.mean([r["group"] == p2.GROUP_A for r in s]), 3),
                     "frac_substitute": rnd(np.mean([r["sub"] for r in s]), 3),
                     "rain24_median": rnd(np.median([r["rain24"] for r in s]), 2),
                     "c6_median": rnd(np.nanmedian([r["c6"] for r in s]), 2),
                     "g6_median": rnd(np.nanmedian([r["g6"] if isnum(r["g6"]) else np.nan for r in s]), 2),
                     "g6_over_rain24_median": rnd(np.median(gfrac), 3) if gfrac else None,
                     "seasons": {q: sum(r["season"] == q for r in s) for q in ("DJF", "MAM", "JJA", "SON")}})
    D["D3_composition"] = comp
    # D4 站群
    for g, tag in ((p2.GROUP_A, "A"), (p2.GROUP_B, "B")):
        X = [r for r in base if r["group"] == g]
        labx = labels(X)
        ok, cells = cell_ok(X, labx, ss_keys(X))
        bt = Boot(ss_keys(X))
        tt, ttb = trend_ratio(bt, labx, arr(X, "obs1"), arr(X, "rim1"))
        D[f"D4_group_{tag}"] = {"n": len(X), "cells": cells, "evaluable": ok,
                                "terciles": tercile_ratios(bt, labx, arr(X, "obs1"), arr(X, "rim1")),
                                "T": rnd(tt), "T_ci95": pct(ttb, Q_95)}
    # D5 雨前趋势订正
    X = subset(base, ("pretrend1",))
    for r in X:
        r["obs1_pt"] = r["obs1"] - r["pretrend1"]
    labx = labels(X)
    bt = Boot(ss_keys(X))
    t0_, t0b = trend_ratio(bt, labx, arr(X, "obs1"), arr(X, "rim1"))
    t1_, t1b = trend_ratio(bt, labx, arr(X, "obs1_pt"), arr(X, "rim1"))
    D["D5_pretrend"] = {"n": len(X), "T_ref": rnd(t0_), "T_c": rnd(t1_), "delta": rnd(t1_ - t0_),
                        "delta_ci95": pct(t1b - t0b, Q_95),
                        "pretrend_mean_by_tercile": [rnd(np.mean(arr(X, "pretrend1")[labx == k]), 5) for k in range(3)],
                        "ctrl_mean_by_tercile": [rnd(np.nanmean([r["ctrl_mean1"] for r, l_ in zip(base, lab) if l_ == k
                                                                 and isnum(r.get("ctrl_mean1"))]), 5) for k in range(3)]}
    # D6 小时合成 R(h)
    comp_h = []
    for k in range(3):
        s = [r for r, l_ in zip(base, lab) if l_ == k]
        row = []
        for h in range(12):
            pr = [(r["_obs_h"][h], r["_mod_h"][h]) for r in s if isnum(r["_obs_h"][h]) and isnum(r["_mod_h"][h])]
            so, sm = sum(a for a, _ in pr), sum(b_ for _, b_ in pr)
            row.append({"h": h, "n": len(pr), "obs_sum": rnd(so, 3), "mod_sum": rnd(sm, 3), "R": rnd(so / sm) if sm else None})
        comp_h.append({"k": k, "hours": row,
                       "obs_peak_hour": int(np.argmin([x["obs_sum"] if x["obs_sum"] is not None else 0 for x in row])),
                       "mod_peak_hour": int(np.argmin([x["mod_sum"] if x["mod_sum"] is not None else 0 for x in row]))})
    D["D6_hourly_composite"] = comp_h
    # D7 联合（只对已判支持的 H2/H3/H5）
    sup = [k for k in ("H2", "H3", "H5") if str(C[k].get("verdict", "")).startswith("支持")]
    D["D7_joint"] = {"supported_used": sup}
    if sup:
        num = "obs1_cc" if "H5" in sup else "obs1"
        den = "rimG1" if "H2" in sup else "rim1"
        X = subset(base, (num, den, "obs1", "rim1"))
        labx = labels(X)
        bt = Boot(ss_keys(X))
        tr, trb = trend_ratio(bt, labx, arr(X, "obs1"), arr(X, "rim1"))
        if "H3" in sup:
            tc, tcb, _, _ = mh_trend(bt, labx, np.array([r["station"] for r in X]), arr(X, num), arr(X, den))
        else:
            tc, tcb = trend_ratio(bt, labx, arr(X, num), arr(X, den))
        D["D7_joint"].update({"n": len(X), "T_ref": rnd(tr), "T_c": rnd(tc), "T_c_ci95": pct(tcb, Q_95),
                              "delta": rnd(tc - tr), "delta_ci95": pct(tcb - trb, Q_95)})
    return D


# ======================================================================== 主流程
def load_m(p4b_dir):
    out = {}
    with open(os.path.join(p4b_dir, "p4b_events.csv"), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                h, p = float(r["hfine100"]), float(r["pcum_cmorph_model"])
            except ValueError:
                h, p = NAN, NAN
            try:
                g = float(r["rimG1"])
            except ValueError:
                g = NAN
            out[(r["station"], r["onset_utc"])] = (h * 1000.0 / p if (math.isfinite(h) and math.isfinite(p) and p > 0) else NAN, g)
    return out


def run(args, out_dir, log):
    import p4_mech as pm
    import p4b_posthoc as pb
    np = _np()
    t0 = time.monotonic()
    s4 = json.load(open(os.path.join(args.p4_dir, "p4_mech_summary.json"), encoding="utf-8"))
    cuts = s4["D_U"]["cuts"]
    if max(abs(cuts[0] - CUTS_REF[0]), abs(cuts[1] - CUTS_REF[1])) > 1e-5:
        raise RuntimeError(f"切点不一致 {cuts}")
    stations, events, val, prim, plan = pm.build(args, out_dir, log)
    data, cov = pm.load_pixels(args.p2_dir, prim)
    sbn = {s["name"]: s for s in stations}
    pixc, pixm = pm.PixView(data, "c"), pm.PixView(data, "m")
    pairs = [pm.event_mech(sbn[e["_stn"]], e, pixc, pixm) for e in events]
    rows = [p[0] for p in pairs]
    v1 = pb.check_v1(rows, os.path.join(args.p4_dir, "p4_mech_events.csv"))
    log.log(f"V1 {json.dumps(v1, ensure_ascii=False)}", echo=True)
    if not v1["pass"]:
        p2.jdump({"V1": v1}, os.path.join(out_dir, "p7c_abort.json"))
        raise RuntimeError("V1 不过")
    mmap = load_m(args.p4b_dir)
    X = []
    for k, e in enumerate(events):
        x = event_p7(sbn[e["_stn"]], e, rows[k], pairs[k][1], pixc)
        x["M"], x["rimG1_p4b"] = mmap.get((x["station"], x["onset_utc"]), (NAN, NAN))
        X.append(x)
    log.log(f"event_p7 完成 {len(X)} 事件 {time.monotonic() - t0:.0f}s", echo=True)
    # V-RIM：重新实现的 RIM 对第一轮 rim1
    dm, nc = 0.0, 0
    for x in X:
        if isnum(x["rim1"]) or isnum(x["rim1_chk"]):
            nc += 1
            dm = max(dm, abs(x["rim1"] - x["rim1_chk"]) if fin(x["rim1"], x["rim1_chk"]) else float("inf"))
    vr = {"n_compared": nc, "max_abs_diff": dm, "pass": bool(dm <= 1e-7)}
    dg = [abs(x["rimG1"] - x["rimG1_p4b"]) for x in X if fin(x["rimG1"], x["rimG1_p4b"])]
    vr["rimG1_vs_p4b"] = {"n": len(dg), "max_abs_diff": max(dg) if dg else None,
                          "n_nan_mismatch": sum(isnum(x["rimG1"]) != isnum(x["rimG1_p4b"]) for x in X)}
    log.log(f"V-RIM {json.dumps(vr, ensure_ascii=False)}", echo=True)
    if not vr["pass"]:
        p2.jdump({"V_RIM": vr}, os.path.join(out_dir, "p7c_abort.json"))
        raise RuntimeError("V-RIM 不过")
    S = analyze(X, log)
    ref = s4["D_U"]["cells"]
    dv0 = max(max(abs(S["V0_DU"][k]["point"] - ref[k]["point"]), abs(S["V0_DU"][k]["ci95"][0] - ref[k]["ci95"][0]),
                  abs(S["V0_DU"][k]["ci95"][1] - ref[k]["ci95"][1])) for k in range(3))
    S["V0_DU_check"] = {"max_abs_diff_vs_p4": dv0, "pass": bool(dv0 <= 1e-5)}
    if not S["V0_DU_check"]["pass"]:
        p2.jdump({"V0": S["V0_DU_check"], "this": S["V0_DU"]}, os.path.join(out_dir, "p7c_abort.json"))
        raise RuntimeError("V0 不过")
    out = {"script": "p7c_wind_dependence.py", "version": VERSION,
           "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "n_events": len(events),
           "rebuild_validation": {k: val[k] for k in ("n_ref", "n_rebuilt", "d4_count")}, "cmorph_pixel_cache": cov,
           "downloads": {"cache_http_requests": plan.get("cache_http_requests"), "new_data_bytes_expected": 0},
           "constants": {"cuts": CUTS_REF, "B": B, "seed": SEED, "n_primary": N_PRIMARY, "q_bonf": Q_BONF,
                         "cell_min": [CELL_MIN_N, CELL_MIN_G], "share_most": SHARE_MOST, "share_conflict": SHARE_CONFLICT,
                         "h1c_point": H1C_POINT, "kappa_grid": [KAPPA_GRID_LO, KAPPA_GRID_HI, KAPPA_N]},
           "V1_round1_csv": v1, "V_RIM": vr}
    out.update(S)
    if not S["V_P4F"]["pass"]:
        raise RuntimeError("V-P4F 不过")
    out["formula"] = formula_analysis()
    out["runtime_s"] = round(time.monotonic() - t0, 1)
    p2.jdump(out, os.path.join(out_dir, "p7c_summary.json"))
    cols = ["station", "group", "season", "onset_utc", "z1", "sub", "u10", "rain24", "g6", "c6", "obs1", "obs5", "rim1", "rim0",
            "rim5", "rim1_lin", "G1", "G5", "G0", "G1_w03", "G1_w612", "rimG1", "M", "obs1_w03", "obs1_w612", "rim1_w03", "rim1_w612",
            "ctrl_mean1", "n_ctrl_valid", "pretrend1", "s0_1", "pcum_cmorph_model", "pcum_gauge_obs"]
    with open(os.path.join(out_dir, "p7c_events.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=cols)
        wr.writeheader()
        for x in X:
            wr.writerow({k: (rnd(x.get(k), 7) if isinstance(x.get(k), float) else x.get(k)) for k in cols})
    print(json.dumps({"overall": S["overall"], "verdicts": {k: v.get("verdict") for k, v in S["candidates"].items()}},
                     ensure_ascii=False), flush=True)
    return 0


# ======================================================================== 自测（合成数据，无网络）
def selftest(out_dir=None):
    np = _np()
    t0 = time.monotonic()
    res = []

    def check(cond, msg):
        res.append({"ok": bool(cond), "msg": msg})

    rng = np.random.default_rng(1)
    # 1 RIM 重新实现 = p2.rim_factor
    worst = 0.0
    for trial in range(30):
        P = np.where(rng.random(84) < 0.3, rng.gamma(1.0, 8.0, 84), 0.0)
        U = rng.uniform(0.5, 16, 84)
        for z in (0.0, 0.5, 1.0, 1.37, 5.0):
            a, b_ = rim_F(P, U, z), p2.rim_factor(P, U, z)
            m = np.isfinite(a) & np.isfinite(b_)
            check(np.array_equal(np.isfinite(a), np.isfinite(b_)), f"rim_F NaN 位置一致 t{trial} z{z}") if trial < 2 else None
            worst = max(worst, float(np.max(np.abs(a[m] - b_[m]))) if m.any() else 0.0)
    check(worst <= 1e-12, f"rim_F 对 p2.rim_factor max|Δ|={worst:.2e} ≤1e-12")
    # 2 线性化：小雨时≈、恒有 1−F_lin ≥ 1−F
    P = np.zeros(84)
    P[60] = 0.2
    U = np.full(84, 6.0)
    a, b_ = 1 - rim_F(P, U, 1.0), 1 - rim_F_lin(P, U, 1.0)
    rel = float(np.max(np.abs(b_ - a) / np.maximum(a, 1e-30))[()])
    check(rel < 1e-3, f"小雨线性化相对差 {rel:.1e}<1e-3")
    P[60:64] = 60.0
    U = np.full(84, 2.0)
    check(np.all(1 - rim_F_lin(P, U, 0.0) >= 1 - rim_F(P, U, 0.0) - 1e-15), "去饱和使淡化不减")
    # 3 守恒 Gaussian：0–200 m 积分＝窗内雨量
    P = np.zeros(84)
    P[40:44] = 12.0
    P[70] = 5.0
    U = np.full(84, 5.0)
    P1 = np.zeros(84)
    P1[60:62] = 10.0
    zz = np.linspace(0, 200, 40001)
    Gz = np.array([1 - cons_G(P, U, z)[-1] for z in zz])
    integ = float(np.sum((Gz[1:] + Gz[:-1]) / 2 * np.diff(zz)))
    rain_m = float(np.sum(P[83 - 48:83]) * 0.5 / 1000 + P[83] * 0.5 / 1000)
    check(abs(integ / rain_m - 1) < 1e-3, f"G 深度积分/雨量＝{integ / rain_m:.5f}")
    # 4 κ 单调：G5/G1 随 κ 增
    kg = kappa_grid()
    rr = [win_delta_h(hourly(cons_G(P1, U, 5.0, k)), (0, 6)) / win_delta_h(hourly(cons_G(P1, U, 1.0, k)), (0, 6)) for k in kg]
    rr = np.array(rr)
    check(np.all(np.diff(rr[np.isfinite(rr)]) > -1e-12), "单脉冲（雨前干）G5/G1 随 κ 单调不减")
    lk_ = np.log(kg)
    curve = np.sin(lk_)                                   # 非单调曲线：多交点取 |ln κ| 最小
    sol, e_ = kappa_solve(curve, 0.5, lk_)
    check(abs(sol - math.asin(0.5)) < 0.01 and not e_, f"κ 多交点取最近 κ=1：{sol:.4f}")
    sol, e_ = kappa_solve(curve, 2.0, lk_)
    check(e_, "κ 无交点记越界")
    # 5 Boot 与 p2.cluster_boot 同一抽样
    keys = [("S%d" % (j % 7), "Q%d" % (j % 3)) for j in range(60)]
    num = rng.normal(-1, 0.5, 60)
    den = rng.normal(-2, 0.3, 60)
    ref, K = p2.cluster_boot(keys, {"r": (num, den)})
    bt = Boot(keys)
    pt, bs = bt.sums(np.stack([num, den], axis=1))
    check(K == bt.K and np.allclose(bs[:, 0] / bs[:, 1], ref["r"], rtol=0, atol=1e-12), "Boot 与 p2.cluster_boot 逐次一致")
    # 6 trend_ratio 对照直接计算
    lab = np.array([j % 3 for j in range(60)])
    tp, tb = trend_ratio(bt, lab, num, den)
    direct = math.log(num[lab == 2].sum() / den[lab == 2].sum()) - math.log(num[lab == 0].sum() / den[lab == 0].sum())
    check(abs(tp - direct) < 1e-12, "trend_ratio 点估计")
    # 7 MH：层间比值不同、层内高/低同比 r → RR＝r
    strata = np.array(["a", "b", "c"])[(np.arange(60) // 3) % 3]
    r_true = 0.5
    base_ = np.where(strata == "a", 0.2, np.where(strata == "b", 0.6, 1.1))
    den2 = -np.abs(rng.normal(1, 0.2, 60))
    num2 = base_ * den2 * np.where(lab == 2, r_true, 1.0)
    mp, mb, use, _ = mh_trend(bt, lab, strata, num2, den2)
    check(abs(mp - math.log(r_true)) < 1e-12 and np.nanmax(np.abs(mb - math.log(r_true))) < 1e-9, f"MH 恒等 r：{mp:.6f}")
    # 层混杂：原始比值趋势≠0 而 MH＝0
    lab3 = np.where(strata == "a", 0, np.where(strata == "c", 2, 1))
    lab3[::5] = 0
    lab3[1::7] = 2
    num3 = base_ * den2
    tr3, _ = trend_ratio(bt, lab3, num3, den2)
    m3, _, _, _ = mh_trend(bt, lab3, strata, num3, den2)
    check(abs(tr3) > 0.3 and abs(m3) < 1e-12, f"MH 去掉层混杂：原始 {tr3:.3f}→MH {m3:.2e}")
    # 8 κ 拟合回收
    P = np.zeros(84)
    P[60:62] = 10.0
    rows = []
    for j in range(90):
        u = [4.0, 7.5, 11.0][j % 3] + 0.3 * rng.standard_normal()
        U = np.full(84, u)
        kt = [0.5, 1.0, 3.0][j % 3]
        g1 = [35 * win_delta_h(hourly(cons_G(P, U, 1.0, k)), (0, 6)) for k in kg]
        g5 = [35 * win_delta_h(hourly(cons_G(P, U, 5.0, k)), (0, 6)) for k in kg]
        o1 = 35 * win_delta_h(hourly(cons_G(P, U, 1.0, kt)), (0, 6))
        o5 = 35 * win_delta_h(hourly(cons_G(P, U, 5.0, kt)), (0, 6))
        rows.append({"station": f"S{j}", "season": "DJF", "u10": u, "obs1": o1, "obs5": o5, "_G1k": g1, "_G5k": g5})
    h = h1c_block(rows)
    check(abs(h["tercile0"]["kappa"] - 0.5) / 0.5 < 0.05 and abs(h["tercile2"]["kappa"] - 3.0) / 3.0 < 0.05,
          f"κ 回收 {h['tercile0']['kappa']}／{h['tercile2']['kappa']}")
    check(abs(h["tercile0"]["rho_prime"] - 1) < 0.05 and abs(h["tercile2"]["rho_prime"] - 1) < 0.05, "ρ′≈1（真模型即守恒 κ_k）")
    check(h["verdict"] == "形状拟合的守恒模型能定量解释冲突", f"H1c 判读 {h['verdict']}")
    # 9 常数与判读函数
    check(abs(Q_BONF[0] - 100 * 0.05 / 9 / 2) < 1e-12, "Bonferroni 百分位")
    check(judge_artifact([0.1, 0.9], 0.6, -0.86, True) == "支持（解释原始趋势大部分）", "judge_artifact 大部分")
    check(judge_artifact([0.1, 0.9], 0.3, -0.86, True) == "支持（解释原始趋势一部分）", "judge_artifact 一部分")
    check(judge_artifact([-0.1, 0.9], 0.6, -0.86, True) == "不支持", "judge_artifact 不支持")
    check(judge_artifact([0.1, 0.9], 0.6, -0.86, False) == "不可评", "judge_artifact 不可评")
    check(judge_conflict([0.2, 1.5], 0.9) == "足以消除冲突" and judge_conflict([0.2, 1.5], 0.3) == "部分缓解冲突"
          and judge_conflict([-0.2, 1.5], 0.9) == "未缓解冲突", "judge_conflict")
    check(judge_direction([0.1, 0.5], +1, True) == "支持" and judge_direction([-0.5, -0.1], +1, True) == "与预言相反"
          and judge_direction([-0.1, 0.5], -1, True) == "未检出" and judge_direction([-0.5, -0.1], -1, True) == "支持",
          "judge_direction")
    check(judge_h1("支持", "未检出") == "部分支持" and judge_h1("支持", "支持") == "支持"
          and judge_h1("未检出", "与预言相反") == "不支持（与预言相反）", "judge_h1")
    # 10 公式分析可运行、方向：RIM 每 mm 随 U 降，G/RIM 随 U 升（1 m）
    fa = formula_analysis()
    rows_f = fa["cases"]["10mm_h_x1h"]["rows"]
    r45 = next(r for r in rows_f if r["U10"] == 4.5)
    r116 = next(r for r in rows_f if r["U10"] == 11.6)
    check(r116["z1.0"]["RIM_per_mm"] < r45["z1.0"]["RIM_per_mm"], "公式：RIM(1 m) 每 mm 随 U 降")
    check(r116["z1.0"]["ratio_G_over_RIM"] > r45["z1.0"]["ratio_G_over_RIM"], "公式：G/RIM(1 m) 随 U 升（守恒预测方向）")
    # 11 端到端（合成事件表）：站混杂造出的风速趋势应被 H3 识别，H2（雨量计与 CMORPH 同）应不支持
    S_ = analyze(_synthetic_rows(rng), check_p4f=False)
    Cc = S_["candidates"]
    check(S_["T_raw_646"]["point"] < -0.3, f"合成：原始趋势 {S_['T_raw_646']['point']}")
    check(Cc["H3"]["verdict"].startswith("支持") and abs(Cc["H3"]["T_c"]) < 0.1, f"合成：H3 {Cc['H3']['verdict']} T_c={Cc['H3']['T_c']}")
    check(Cc["H2"]["verdict"] == "不支持", f"合成：H2 {Cc['H2']['verdict']}")
    check(isinstance(S_["overall"], str) and "D7_joint" in S_["descriptive"], "合成：总映射与描述块可运行")
    res = [r for r in res if r is not None]
    fails = [r for r in res if not r["ok"]]
    out = {"version": VERSION, "n": len(res), "n_fail": len(fails), "fails": fails, "runtime_s": round(time.monotonic() - t0, 2)}
    if out_dir:
        p2.jdump(out, os.path.join(out_dir, "p7c_selftest.json"))
    for r in fails:
        print("FAIL", r["msg"], flush=True)
    print(f"自测 {len(res) - len(fails)}/{len(res)} 通过（{out['runtime_s']} s）", flush=True)
    return 0 if not fails else 4


def _synthetic_rows(rng):
    """自测用合成事件表：两类站（比值 0.5 与 0.15），高风事件集中在低比值站；站内无风速趋势。"""
    np = _np()
    kg = kappa_grid()
    rows = []
    for j in range(480):
        st = f"S{j % 12}"
        low_ratio = (j % 12) >= 6
        u = rng.uniform(9.5, 13) if (low_ratio and rng.random() < 0.7) else rng.uniform(2, 13)
        if not low_ratio and rng.random() < 0.6:
            u = rng.uniform(2, 6)
        ratio = 0.15 if low_ratio else 0.5
        rim1 = -abs(rng.normal(0.2, 0.05)) / (u / 6.0)
        o1 = ratio * rim1 * (1 + 0.1 * rng.standard_normal())
        g1 = [rim1 * 0.3 * (k ** -0.5) for k in kg]
        g5 = [rim1 * 0.3 * (k ** -0.5) * min(1.0, 0.2 * k) for k in kg]
        rows.append({"station": st, "group": p2.GROUP_A if j % 2 else p2.GROUP_B, "season": ["DJF", "MAM", "JJA", "SON"][j % 4],
                     "onset_utc": f"t{j}", "z1": 1.0 if j % 3 else 1.5, "sub": 0 if j % 3 else 1, "rain24": float(rng.uniform(10, 60)),
                     "u10": float(u), "obs1": o1, "obs5": 0.4 * o1, "rim1": rim1, "rim0": 1.2 * rim1, "rim5": 0.3 * rim1,
                     "rim1_lin": 1.1 * rim1, "G1": 0.3 * rim1, "G5": 0.1 * rim1, "G0": 0.4 * rim1, "rimG1": rim1 * 0.9,
                     "M": 4.0, "G1_w03": 0.3 * rim1, "G1_w612": 0.2 * rim1, "obs1_w03": o1, "obs1_w612": 0.7 * o1, "rim1_w03": rim1, "rim1_w612": 0.7 * rim1,
                     "ctrl_mean1": 0.001 * rng.standard_normal(), "n_ctrl_valid": 5, "pretrend1": 0.001 * rng.standard_normal(),
                     "s0_1": 35.0, "pcum_cmorph_model": 10.0, "pcum_gauge_obs": 9.0, "g6": 9.0, "c6": 10.0, "model_ok": True,
                     "_obs_h": [o1] * 12, "_mod_h": [rim1] * 12, "_G1k": g1, "_G5k": g5})
    return rows


def main(argv=None):
    import p4_mech as pm
    import p4b_posthoc as pb
    ap = argparse.ArgumentParser(description="P7c 矛盾 C（R 随风速下降）候选解释检验")
    ap.add_argument("--out")
    ap.add_argument("--p1-events", default=p2.P1_EVENTS_DEFAULT)
    ap.add_argument("--p1b-dir", default=p2.P1B_DIR_DEFAULT)
    ap.add_argument("--p2-dir", default=pm.P2_DIR_DEFAULT)
    ap.add_argument("--p4-dir", default=pb.P4_DIR_DEFAULT)
    ap.add_argument("--p4b-dir", default=P4B_DIR_DEFAULT)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--formula", action="store_true")
    args = ap.parse_args(argv)
    out_dir = args.out or os.environ.get("REPRO_OUTPUT_DIR")
    if args.formula:
        fa = formula_analysis()
        txt = json.dumps(fa, ensure_ascii=False, indent=1)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
            with open(os.path.join(out_dir, "p7c_formula.json"), "w", encoding="utf-8") as f:
                f.write(txt)
        print(txt)
        return 0
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
    log = p1.Log(os.path.join(out_dir, "p7c_log.txt"))
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
