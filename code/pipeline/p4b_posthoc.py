#!/usr/bin/env python3
"""p4b_posthoc.py — P4 第二轮探索性诊断 E1–E3：第一轮没过的 P4（雨量错配）、P2a（理论排序）、P3（淡水去向）
能否用站得住的方法论证。第一轮判定（negative_L）不改；本脚本全部结果属探索性分析。

方法与判读规则在运行前写定；判读逻辑在 judge_*()／total_reading()。

用法：
  python p4b_posthoc.py --selftest          # 合成数据自测（无网络、约十几秒）
  python p4b_posthoc.py                     # 正式（先跑自测，不过即退出 4）
  --p1-events / --p1b-dir / --p2-dir        同 p4_mech（默认 P1/P1b/P2 产物）
  --p4-dir DIR                              第一轮输出目录（读 p4_mech_events.csv 做逐事件一致性核对）
  --out DIR                                 输出目录（默认 $REPRO_OUTPUT_DIR）
依赖：numpy、scipy＋同目录 p1_events.py、p2_rim_test.py、p4_mech.py（只 import，不改）。
数据：只读现有缓存（P1/P1b 缓存、P2 主集 cmorph_pixels.csv、第一轮 p4_mech_events.csv）；无新下载。
产物（<out>/）：p4b_selftest.json、p4b_summary.json（含各项判读与总判读）、p4b_e1_null.csv（每次模拟的形状量）、
  p4b_events.csv（逐事件新增派生量）、p4b_log.txt。退出码：0 跑完（无论判读）；3 异常或一致性核对不过；4 自测不过。

实现选择（N 条）：
  N1 事件与逐事件第一轮派生量：p4_mech.build＋event_mech 原样重算，逐事件与第一轮 CSV 核 obs1/rim1/g6/c6/lam/u10/low1
     （|差|≤2e-6），V0 与 P2 核（1e-4）；不一致即 p4b_abort.json、退出 3。
  N2 E1(i) 标量响应曲线：f_e(a)＝s0·Δ[RIM(a·P_G; z1)]（0–6 h 窗），a∈{0}∪logspace(−3,3,241)（第 120 点置为精确 1），
     批量 RIM（rim_factor_batch，逐行等于 p2.rim_factor）；模拟时 ln a 线性插值，0<a<1e−3 在 (0,f(0)) 与 (1e−3,…) 间线性，
     a>1e3 或落到 NaN 格点→该事件本次剔除并计数。
  N3 E1 分箱：群(A/B)×该群 ln(G6+0.1) 五分位（np.quantile 0.2..0.8，≤切点归低箱），箱内真实 λ 有放回重抽；G6＝Σ P_G[60:72]·0.5。
  N4 E1 θ=0：T＝P_G；C*＝k·P_G，k＝max((G6+0.1)e^{−ℓ*}−0.1,0)/G6（G6=0 → C*=0）；λ*＝ln((G6+0.1)/(k·G6+0.1))。
     θ=0.5：ε_g～N(θ(ℓ*−m_b), θ(1−θ)v_b)，T＝e^{−ε_g}·P_G，C* 同上。噪声 η：同站对照窗 ΔS(1 m) 合并池有放回（空池用全站池）。
     R0＝0.296（第一轮主值，常数 R0_MAIN）；敏感性 R0＝(ii) 点估计。N_SIM＝2000，PCG64(20260926)，抽取顺序固定
     （箱→事件、站→事件、ε_g），三个设定各自新建同 seed 生成器。
  N5 E1 包络：n_k＝R_k/R_all（k=0,1,2）与 D＝ln(R_2/R_0) 四量；联合 95%＝Bonferroni 各 [0.625, 99.375] 分位；逐点 95% 仅信息。
  N6 E1(ii)/(iii)：雨量计驱动 RIM-3（p2.event_forcing 'gauge'，K21 负值置 0）的 R_G 与按 λ 三分位的 R_G（p4_mech.tercile_block）。
  N7 E2 数值解：非均匀网格（顶格 0.02 m，×1.03，至 ≥2000 m），隐式欧拉，每半步 30×60 s，半步中点把 p·0.5 mm 加入顶格，
     LAPACK dgttrf/dgttrs（每半步分解一次）；K_f＝κu*z_f/φ_h(z_f·κB/u*³)；B＝g·β_S·S0·P（P 该半步雨强），S0＝1 m 雨前中位数；
     观测深度按格心线性插值，海面取顶格；小时值＝两个半步末之均；ΔS 同第一轮 delta_at。变体：中性（V2 核对）、β_h=5
     （主）、β_h=7.8、雨量计驱动 β_h=5。无雨的前段（c≡0）跳过求解。
  N8 E2 统计：ρ_s＝stab1/rim1（rim1≤−0.005）；三分位用 p4_mech.tercile_block（theory=(stab1,rim1)）；Spearman(ρ_s, obs1/rim1) 的
     站×季整簇 bootstrap（与 cluster_boot 同一 PCG64 seed 与 draw 方式）；Λ_s 用 p4_mech.pooled。
  N9 E3 权重：D=10 {1:3,5:4.5,10:2.5}；D=20 {1:3,5:4.5,10:7.5,20:5}；D=25 有 20 m {1:3,5:4.5,10:7.5,20:7.5,25:2.5}、
     无 20 m {1:3,5:4.5,10:10,25:7.5}；S0＝各深度雨前中位数（须 >0）。φ 与 Q、R1 同一次抽样（cluster_boot）。
  N10 E3 窗：[0,6)、[6,12)、[12,24)、[24,48) h；CMORPH 驱动只算前两窗（缓存只到 t0+12 h）；雨量计驱动强迫延长到
     [H−30, H+48)（同 p2.event_forcing gauge 分支：K9 风插补与下限、K21），R1 同驱动同窗。
  N11 E3 MDE＝2.8016·SE_φ，SE_φ＝(P97.5−P2.5)/3.92（站×季 bootstrap）；另报 SD。
  N12 E3 闭合比：h_obs·1000/雨量计累计（箱中点，从 t0 起，窗内小时平均；负值置 0，缺测→该事件剔除）；
     M_RIM：CMORPH RIM-3 在 z 网格（0–10 m 每 0.05、10–100 m 每 0.5）上 −Δ(F) 的梯形积分 ÷ 模型时刻 CMORPH 累计（M10）。
  N14 V3（信息，不中止）：CMORPH 0–6 h 的 φ_10 与第一轮 D_B.phi（p4_mech_summary.json）核对（n 相同、|差|≤1e-4）。
  N13 判读、总映射全部在 judge_e1/judge_e2/judge_e3/total_reading。

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

VERSION = "p4b-2026-09-27a"
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
P4_DIR_DEFAULT = _rp.upstream("p4_dir")  # [repro] 读 p4-mech 阶段输出目录
NAN = float("nan")
isnum, rnd = pm.isnum, pm.rnd

# ---- 事先写定的常数 ----
R0_MAIN = 0.296
N_SIM = 2000
SIM_SEED = 20260926
A_GRID_N = 241
ENV_JOINT = (0.625, 99.375)
ENV_POINT = (2.5, 97.5)
THETAS = (0.0, 0.5)
BETA_MAIN, BETA_ALT, GAMMA_U = 5.0, 7.8, 16.0
NSUB = 30
GRID_DZ0, GRID_R, GRID_H = 0.02, 1.03, 2000.0
V2_TOL = 0.03
LEVEL_BAND = (0.5, 2.0)
MDE_Z = p2.MDE_Z                                   # 1.959964+0.841621
PHI_REF = 0.24
MDE_BIG = 0.5
WINDOWS = ((0, 6), (6, 12), (12, 24), (24, 48))
EXT_HI_H = 48
W_D = {10: {1.0: 3.0, 5.0: 4.5, 10.0: 2.5},
       20: {1.0: 3.0, 5.0: 4.5, 10.0: 7.5, 20.0: 5.0},
       "25w20": {1.0: 3.0, 5.0: 4.5, 10.0: 7.5, 20.0: 7.5, 25.0: 2.5},
       "25n20": {1.0: 3.0, 5.0: 4.5, 10.0: 10.0, 25.0: 7.5}}
ZFINE = None
V1_FIELDS = ("obs1", "rim1", "g6", "c6", "lam", "u10", "low1")
V1_TOL = 2e-6


def _np():
    import numpy
    return numpy


# ======================================================================== RIM 批量与标量响应（N2）
def rim_factor_batch(P, U, z, t_cur_depth=p2.T_CURRENT_S):
    """p2.rim_factor 的批量版：P (S,L)、U (L,)；返回 (S, L−48)。逐行与 p2.rim_factor 相同（自测核）。"""
    np = _np()
    from numpy.lib.stride_tricks import sliding_window_view as swv
    P = np.atleast_2d(np.asarray(P, float))
    U = np.asarray(U, float)
    L = P.shape[-1]
    tc = np.array(p2.TC, float)
    ti = np.array(p2.TI, float)
    Pw = swv(P, 48, axis=-1)[:, :L - 48, :]
    Uw = np.broadcast_to(swv(U, 48)[:L - 48], Pw.shape)
    irr = Pw / 1000.0 / 3600.0
    kz = p2.KZ_COEF * Uw ** 2
    d0 = p2.d0_interp(Uw, Pw)
    Pc = P[:, 48:]
    Uc = np.broadcast_to(U[48:], Pc.shape)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        term = p2.C1 * irr * tc / np.sqrt(kz * ti)
        if z:
            term = term * np.exp(-(z ** 2) / (4.0 * kz * ti))
        prior = d0 / (d0 + term)
        irrc = Pc / 1000.0 / 3600.0
        kzc = p2.KZ_COEF * Uc ** 2
        d0c = p2.d0_interp(Uc, Pc)
        cterm = p2.C2 * irrc * 1800.0 / np.sqrt(kzc * p2.T_CURRENT_S)
        if z:
            cterm = cterm * np.exp(-(z ** 2) / (4.0 * kzc * t_cur_depth))
        cur = d0c / (d0c + cterm)
    return np.prod(prior, axis=-1) * cur


def rim_multi_z(P, U, zs):
    """一条强迫、多个深度：返回 (nz, L−48)，逐深度与 p2.rim_factor(P,U,z) 相同（自测核）。"""
    np = _np()
    from numpy.lib.stride_tricks import sliding_window_view as swv
    P = np.asarray(P, float)
    U = np.asarray(U, float)
    zs = np.asarray(zs, float)
    L = len(P)
    tc = np.array(p2.TC, float)
    ti = np.array(p2.TI, float)
    Pw = swv(P, 48)[:L - 48]
    Uw = swv(U, 48)[:L - 48]
    irr = Pw / 1000.0 / 3600.0
    kz = p2.KZ_COEF * Uw ** 2
    d0 = p2.d0_interp(Uw, Pw)
    Pc, Uc = P[48:], U[48:]
    kzc = p2.KZ_COEF * Uc ** 2
    d0c = p2.d0_interp(Uc, Pc)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        term = p2.C1 * irr * tc / np.sqrt(kz * ti)
        prior = d0[None] / (d0[None] + term[None] * np.exp(-(zs[:, None, None] ** 2) / (4.0 * kz * ti)[None]))
        cterm = p2.C2 * (Pc / 1000.0 / 3600.0) * 1800.0 / np.sqrt(kzc * p2.T_CURRENT_S)
        cur = d0c[None] / (d0c[None] + cterm[None] * np.exp(-(zs[:, None] ** 2) / (4.0 * kzc * p2.T_CURRENT_S)[None]))
    return np.prod(prior, axis=-1) * cur


def trap(y, x):
    """梯形积分（不依赖 numpy.trapz／trapezoid 的版本差异）。"""
    np = _np()
    y, x = np.asarray(y, float), np.asarray(x, float)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def delta_rows(hourly, pre=slice(0, 6), post=slice(6, 12)):
    """hourly (S,H)：每行 mean(post)−median(pre)；任一非有限→NaN。"""
    np = _np()
    a, b = hourly[:, pre], hourly[:, post]
    ok = np.all(np.isfinite(a), axis=1) & np.all(np.isfinite(b), axis=1)
    out = np.where(ok, np.mean(np.where(np.isfinite(b), b, 0), axis=1) - np.median(np.where(np.isfinite(a), a, 0), axis=1),
                   np.nan)
    return out


def a_grid():
    np = _np()
    g = np.logspace(-3, 3, A_GRID_N)
    g[(A_GRID_N - 1) // 2] = 1.0
    return np.concatenate([[0.0], g])


def scale_curve(P, U, z1, s0):
    """f(a)＝s0·Δ[RIM(a·P; z1)]，a∈a_grid()。"""
    np = _np()
    A = a_grid()
    F = rim_factor_batch(A[:, None] * np.asarray(P, float)[None, :], U, z1)
    hourly = F.reshape(F.shape[0], -1, 2).mean(axis=2)
    return s0 * delta_rows(hourly)


def interp_curve(f, a):
    """在 a_grid 上的曲线 f，于任意 a≥0 取值（N2）。a 可为数组。"""
    np = _np()
    A = a_grid()
    a = np.asarray(a, float)
    out = np.full(a.shape, np.nan)
    la = np.log(A[1:])
    big = a > A[-1]
    mid = (a >= A[1]) & ~big
    small = (a > 0) & (a < A[1])
    out[a == 0] = f[0]
    out[mid] = np.interp(np.log(a[mid]), la, f[1:])
    out[small] = f[0] + (f[1] - f[0]) * a[small] / A[1]
    # 落到 NaN 格点：np.interp 会把 NaN 传播到相邻区间，保持 NaN
    return out


# ======================================================================== E1 模拟（N3–N5）
def e1_bins(groups, g6):
    np = _np()
    g6 = np.asarray(g6, float)
    lg = np.log(g6 + 0.1)
    b = np.full(len(g6), -1, int)
    for gi, g in enumerate(sorted(set(groups))):
        m = np.array([x == g for x in groups])
        cuts = np.quantile(lg[m], [0.2, 0.4, 0.6, 0.8])
        b[m] = gi * 5 + np.digitize(lg[m], cuts, right=True)
    return b


def e1_simulate(curves, g6, lam_obs, bins, stations, noise_pools, R0, theta, n_sim, seed):
    """返回 dict：每次模拟的 R_k、R_all、n_k、D，与计数。curves (n, len(a_grid))。"""
    np = _np()
    rng = np.random.Generator(np.random.PCG64(seed))
    n = len(g6)
    g6 = np.asarray(g6, float)
    lam_obs = np.asarray(lam_obs, float)
    ell = np.empty((n_sim, n))
    mb, vb = {}, {}
    for b in sorted(set(bins.tolist())):
        idx = np.nonzero(bins == b)[0]
        pool = lam_obs[idx]
        mb[b], vb[b] = float(pool.mean()), float(pool.var())
        ell[:, idx] = pool[rng.integers(0, len(pool), size=(n_sim, len(idx)))]
    eta = np.empty((n_sim, n))
    allpool = np.concatenate([np.asarray(v, float) for v in noise_pools.values() if len(v)])
    n_empty = 0
    for st in sorted(set(stations)):
        idx = np.nonzero(np.array([s == st for s in stations]))[0]
        pool = np.asarray(noise_pools.get(st, []), float)
        if not len(pool):
            pool, n_empty = allpool, n_empty + 1
        eta[:, idx] = pool[rng.integers(0, len(pool), size=(n_sim, len(idx)))]
    if theta > 0:
        m = np.array([mb[b] for b in bins])
        v = np.array([vb[b] for b in bins])
        eg = theta * (ell - m[None]) + math.sqrt(theta * (1 - theta)) * np.sqrt(v)[None] * rng.standard_normal((n_sim, n))
        aT = np.exp(-eg)
    else:
        aT = np.ones((n_sim, n))
    cstar6 = np.maximum((g6[None] + 0.1) * np.exp(-ell) - 0.1, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.where(g6[None] > 0, cstar6 / g6[None], 0.0)
    lam_s = np.log((g6[None] + 0.1) / (k * g6[None] + 0.1))
    num = np.empty((n_sim, n))
    den = np.empty((n_sim, n))
    for j in range(n):
        num[:, j] = R0 * interp_curve(curves[j], aT[:, j]) + eta[:, j]
        den[:, j] = interp_curve(curves[j], k[:, j])
    valid = np.isfinite(num) & np.isfinite(den)
    q = np.quantile(lam_s, [1 / 3, 2 / 3], axis=1)
    lab = (lam_s > q[0][:, None]).astype(int) + (lam_s > q[1][:, None]).astype(int)
    numz, denz = np.where(valid, num, 0.0), np.where(valid, den, 0.0)
    R = np.empty((n_sim, 3))
    for kk in range(3):
        m = lab == kk
        with np.errstate(divide="ignore", invalid="ignore"):
            R[:, kk] = (numz * m).sum(1) / (denz * m).sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        Rall = numz.sum(1) / denz.sum(1)
        nk = R / Rall[:, None]
        D = np.where(R[:, 0] > 0, np.log(R[:, 2] / R[:, 0]), np.nan)
    return {"R": R, "Rall": Rall, "nk": nk, "D": D, "n_invalid_event_draws": int((~valid).sum()),
            "n_clipped_cstar0": int((cstar6 == 0).sum()), "n_empty_noise_pool_stations": n_empty,
            "n_D_nan": int((~np.isfinite(D)).sum())}


def envelope(sim, obs):
    """obs：{'n0','n1','n2','D'}；返回每量的联合/逐点区间、观测是否在内、观测在 null 中的分位。"""
    np = _np()
    arrs = {"n0": sim["nk"][:, 0], "n1": sim["nk"][:, 1], "n2": sim["nk"][:, 2], "D": sim["D"]}
    out = {}
    for key, a in arrs.items():
        a = a[np.isfinite(a)]
        if len(a) < 20:
            out[key] = {"obs": rnd(obs.get(key)), "n_finite": int(len(a)), "inside_joint": False, "joint95": [None, None]}
            continue
        lo, hi = np.percentile(a, ENV_JOINT)
        plo, phi_ = np.percentile(a, ENV_POINT)
        v = obs.get(key)
        out[key] = {"obs": rnd(v), "joint95": [rnd(lo), rnd(hi)], "pointwise95": [rnd(plo), rnd(phi_)],
                    "median": rnd(float(np.median(a))), "n_finite": int(len(a)),
                    "inside_joint": bool(isnum(v) and lo <= v <= hi),
                    "obs_quantile_in_null": rnd(float(np.mean(a <= v)), 4) if isnum(v) else None}
    out["all_inside_joint"] = all(out[k]["inside_joint"] for k in ("n0", "n1", "n2", "D"))
    lv = {}
    for nm, a in (("R0", sim["R"][:, 0]), ("R1", sim["R"][:, 1]), ("R2", sim["R"][:, 2]), ("Rall", sim["Rall"])):
        a = a[np.isfinite(a)]
        lv[nm] = {"pointwise95": [rnd(np.percentile(a, 2.5)), rnd(np.percentile(a, 97.5))],
                  "median": rnd(float(np.median(a)))} if len(a) >= 20 else {"n_finite": int(len(a))}
    out["levels_info_only"] = lv
    out["counts"] = {k: sim[k] for k in ("n_invalid_event_draws", "n_clipped_cstar0", "n_empty_noise_pool_stations", "n_D_nan")}
    return out


# ======================================================================== E2 稳定修正（N7）
def fd_grid():
    np = _np()
    dz = []
    z = 0.0
    while z < GRID_H:
        d = GRID_DZ0 * GRID_R ** len(dz)
        dz.append(d)
        z += d
    dz = np.array(dz)
    zf = np.concatenate([[0.0], np.cumsum(dz)])
    return dz, zf, 0.5 * (zf[:-1] + zf[1:])


def phi_h(zeta, beta, gamma=GAMMA_U):
    np = _np()
    zeta = np.asarray(zeta, float)
    with np.errstate(invalid="ignore"):
        return np.where(zeta >= 0, 1.0 + beta * zeta, (1.0 - gamma * np.minimum(zeta, 0.0)) ** -0.5)


def fd_solve(P, U, s0, zs, beta=None, buoy_scale=1.0, K_const=None, nsub=NSUB, half_s=pm.HALF_S, grid=None):
    """半步序列强迫（P mm/h、U m/s，长 L）→ {z: c 在各半步末（长 L）}。beta=None 为中性；K_const 给定则用常 K（自测）。
    z=0 取顶格。质量守恒：Σdz·c＝累计入水（自测核）。"""
    np = _np()
    from scipy.linalg import lapack
    dz, zf, zc = grid if grid is not None else fd_grid()
    n = len(dz)
    P = np.asarray(P, float)
    U = np.asarray(U, float)
    L = len(P)
    us = pm.ustar_w(U)
    zi = zf[1:-1]
    dzc = zc[1:] - zc[:-1]
    dt = half_s / nsub
    c = np.zeros(n)
    out = {z: np.zeros(L) for z in zs}
    started = False
    for j in range(L):
        if not started and P[j] <= 0:
            continue
        started = True
        if K_const is not None:
            K = np.full(n - 1, float(K_const))
        else:
            K = pm.KAPPA * us[j] * zi
            if beta is not None and isnum(s0):
                bflux = pm.G_ACC * pm.BETA_S * s0 * max(P[j], 0.0) / 1000.0 / 3600.0 * buoy_scale
                zeta = zi * pm.KAPPA * bflux / max(us[j], 1e-12) ** 3
                K = K / phi_h(zeta, beta)
        g = dt * K / dzc
        d = dz.copy()
        d[:-1] += g
        d[1:] += g
        dl, du = -g.copy(), -g.copy()
        dl_f, d_f, du_f, du2, ipiv, info = lapack.dgttrf(dl, d, du)
        for s in range(nsub):
            if s == nsub // 2:
                c[0] += P[j] * 0.5 / 1000.0 / dz[0]
            x, info2 = lapack.dgttrs(dl_f, d_f, du_f, du2, ipiv, dz * c)
            c = x
        for z in zs:
            out[z][j] = c[0] if z == 0.0 else float(np.interp(z, zc, c))
    return out


def fd_hourly(P, U, s0, zs, **kw):
    np = _np()
    c = fd_solve(P, U, s0, zs, **kw)
    return {z: v[48:].reshape(-1, 2).mean(axis=1) for z, v in c.items()}


# ======================================================================== E3 延长强迫与收支（N9–N12）
def gauge_forcing_ext(stn, H, hi_h=EXT_HI_H):
    """p2.event_forcing 'gauge' 分支，窗口延长到 [H−30, H+hi_h)。返回 (P, U) 半步或 (None, None)。"""
    np = _np()
    i0 = H + p2.WIN_LO_H - stn["h0"]
    nh = hi_h - p2.WIN_LO_H
    pad = p2.WIND_GAP_MAX_H
    ws = stn["ser"]["wind"]
    lo, hi = i0 - pad, i0 + nh + pad
    seg = np.full(hi - lo, np.nan)
    a, b = max(lo, 0), min(hi, len(ws))
    if b > a:
        seg[a - lo:b - lo] = ws[a:b]
    w = p2.fill_wind(seg, p2.WIND_GAP_MAX_H)[pad:pad + nh] * _rp.wind_factor(stn["name"], p2.WIND_FACTOR, p2.WIND_Z0_M)  # [repro] 风高开关 [options] wind_height；缺省 uniform_4m 时 _rp.wind_factor 原样返回原系数
    if not np.all(np.isfinite(w)):
        return None, None
    w = np.maximum(w, p2.WIND_FLOOR)
    rs = stn["ser"]["rain"]
    g = np.full(nh, np.nan)
    a, b = max(i0, 0), min(i0 + nh, len(rs))
    if b > a:
        g[a - i0:b - i0] = rs[a:b]
    g = np.where(np.isfinite(g), np.maximum(g, 0.0), 0.0)
    return np.repeat(g, 2), np.repeat(w, 2)


def weights_for(D, has20):
    if D == 25:
        return W_D["25w20"] if has20 else W_D["25n20"]
    return W_D[D]


def h_col(vals, s0s, w):
    """淡水高 Σw·(−ΔS/S0)；vals/s0s: {z: 值}。任一缺→NaN。"""
    s = 0.0
    for z, wt in w.items():
        dv, s0 = vals.get(z), s0s.get(z)
        if not (isnum(dv) and isnum(s0)) or s0 <= 0:
            return NAN
        s += wt * (-dv / s0)
    return s


def rain_mid_window(r_hourly, win):
    """N12：从 t0 起的箱中点累计雨量在窗内小时上的平均（mm）；r_hourly 从 t0 起。"""
    tot, acc = 0.0, []
    for x in r_hourly:
        acc.append(tot + 0.5 * x)
        tot += x
    sel = acc[win[0]:win[1]]
    return sum(sel) / len(sel) if len(sel) == win[1] - win[0] else NAN


# ======================================================================== 逐事件派生量
DEEP_Z = (1.0, 5.0, 10.0, 20.0, 25.0)
SKEY = {1.0: "s1", 5.0: "s5", 10.0: "s10", 20.0: "s20", 25.0: "s25"}


def event_extra(stn, e, row, pixc, grid):
    """一个事件的 E1–E3 新增派生量（写回 row 的副本）。"""
    np = _np()
    r = dict(row)
    ser, i, H = stn["ser"], e["i"], e["hour"]
    z1 = r["z1"]
    s0_1 = r.get("s0_1")
    r["ctrl_ds1"] = [x for x in e.get("ctrl_ds1", []) if isnum(x)]
    r["curve"] = None
    P, U, _ = p2.event_forcing(stn, H, "cmorph", pixc)
    Pg, Ug, _ = p2.event_forcing(stn, H, "gauge", pixc)
    # ---- E1
    if Pg is not None and isnum(s0_1):
        Fg = p2.rim_factor(Pg, Ug, z1).reshape(-1, 2).mean(axis=1)
        r["rimG1"] = s0_1 * pm.delta_at(Fg, list(range(6)), list(range(6, 12)))
        r["G6sim"] = float(np.sum(Pg[pm.I_T0:pm.I_T0 + 12]) * 0.5)
        r["curve"] = scale_curve(Pg, Ug, z1, s0_1)
    else:
        r["rimG1"] = NAN
    # ---- E2
    if P is not None and isnum(s0_1):
        zs = (0.0, z1, 5.0, 10.0)

        def ds(h, z):
            return s0_1 * pm.delta_at(-h[z], list(range(6)), list(range(6, 12)))
        hn = fd_hourly(P, U, s0_1, (z1,), beta=None, grid=grid)
        r["fdn1"] = ds(hn, z1)
        h5 = fd_hourly(P, U, s0_1, zs, beta=BETA_MAIN, grid=grid)
        for zl, zz in ((0, 0.0), (1, z1), (5, 5.0), (10, 10.0)):
            s0z = s0_1 if zl in (0, 1) else r.get(f"s0_{zl}")
            r[f"stab{zl}"] = s0z * pm.delta_at(-h5[zz], list(range(6)), list(range(6, 12))) if isnum(s0z) else NAN
        h78 = fd_hourly(P, U, s0_1, (z1,), beta=BETA_ALT, grid=grid)
        r["stab78_1"] = ds(h78, z1)
        if Pg is not None:
            hg = fd_hourly(Pg, Ug, s0_1, (0.0, z1), beta=BETA_MAIN, grid=grid)
            r["stabG1"], r["stabG0"] = ds(hg, z1), ds(hg, 0.0)
        r["rho_s"] = r["stab1"] / r["rim1"] if (isnum(r.get("stab1")) and isnum(r.get("rim1"))
                                                and r["rim1"] <= pm.RHO_MIN_DEN) else NAN
        # ---- E3d：RIM 细网格质量
        Fz = rim_multi_z(P, U, ZFINE)
        hz = Fz.reshape(len(ZFINE), -1, 2).mean(axis=2)
        dF = delta_rows(hz)
        r["hfine100"] = trap(-dF, ZFINE) if np.all(np.isfinite(dF)) else NAN
        m10 = ZFINE <= 10.0 + 1e-9
        r["hfine10"] = trap(-dF[m10], ZFINE[m10]) if np.all(np.isfinite(dF[m10])) else NAN
        pick = {z: float(dF[int(np.argmin(np.abs(ZFINE - z)))]) for z in (1.0, 5.0, 10.0)}
        r["h3pt10"] = sum(w * (-pick[z]) for z, w in W_D[10].items())
    # ---- E3 收支（仅 A 群）
    if r["group"] == p2.GROUP_A:
        s0s = {z: p2.pre_median(ser[SKEY[z]], i) for z in DEEP_Z}
        s0s[1.0] = s0_1
        r["_s0s"] = s0s
        obs = {}
        for w in WINDOWS:
            obs[w] = {z: p2.obs_delta(ser[SKEY[z]], i, w) for z in DEEP_Z}
        r["_obs"] = obs
        zmod = {1.0: z1, 5.0: 5.0, 10.0: 10.0, 20.0: 20.0, 25.0: 25.0}
        modc = {}
        if P is not None:
            hr = {z: pm.rim_hourly(P, U, zmod[z]) for z in DEEP_Z}
            for w in WINDOWS[:2]:
                modc[w] = {z: s0s[z] * pm.delta_at(hr[z], list(range(6)), list(range(6 + w[0], 6 + w[1])))
                           if isnum(s0s[z]) else NAN for z in DEEP_Z}
        r["_modc"] = modc
        modg = {}
        Pe, Ue = gauge_forcing_ext(stn, H)
        if Pe is not None:
            hr = {z: p2.rim_factor(Pe, Ue, zmod[z]).reshape(-1, 2).mean(axis=1) for z in DEEP_Z}
            for w in WINDOWS:
                modg[w] = {z: s0s[z] * pm.delta_at(hr[z], list(range(6)), list(range(6 + w[0], 6 + w[1])))
                           if isnum(s0s[z]) else NAN for z in DEEP_Z}
        r["_modg"] = modg
        rr = np.asarray(ser["rain"][i:i + EXT_HI_H], float)
        if len(rr) == EXT_HI_H and np.all(np.isfinite(rr)):
            rr = np.maximum(rr, 0.0)
            r["_rainmid"] = {w: rain_mid_window(list(rr), w) for w in WINDOWS}
        else:
            r["_rainmid"] = {}
            if len(rr):
                rr2 = np.where(np.isfinite(rr), np.maximum(rr, 0.0), np.nan)
                for w in WINDOWS:
                    seg = rr2[:w[1]]
                    if len(seg) == w[1] and np.all(np.isfinite(seg)):
                        r["_rainmid"][w] = rain_mid_window(list(seg), w)
    return r


# ======================================================================== 统计
def spearman_np(x, y):
    np = _np()
    from scipy.stats import rankdata
    if len(x) < 3:
        return NAN
    rx, ry = rankdata(x), rankdata(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    d = math.sqrt(float(rx @ rx) * float(ry @ ry))
    return float(rx @ ry) / d if d > 0 else NAN


def spearman_cluster(keys, x, y):
    np = _np()
    x, y = np.asarray(x, float), np.asarray(y, float)
    uk = sorted(set(keys))
    kid = {k: j for j, k in enumerate(uk)}
    idx = np.array([kid[k] for k in keys])
    members = [np.nonzero(idx == j)[0] for j in range(len(uk))]
    rng = np.random.Generator(np.random.PCG64(pm.SEED))
    draw = rng.integers(0, len(uk), size=(pm.B, len(uk)))
    bs = np.empty(pm.B)
    for b in range(pm.B):
        sel = np.concatenate([members[k] for k in draw[b]])
        bs[b] = spearman_np(x[sel], y[sel])
    ci, nb = pm.ci95(bs)
    return {"n": len(x), "G": len(uk), "rho": rnd(spearman_np(x, y), 4), "ci95": ci, "n_boot": nb}


def phi_block(rows, hobs, hrim, o1, m1):
    """φ＝(Q−R1)/(1−R1) 同一次抽样；返回点估计、站×季/站簇 CI、逐站剔除、SE、MDE。rows 已过滤。"""
    np = _np()
    res = {"n": len(rows), "G_station_season": len({(r["station"], r["season"]) for r in rows}),
           "G_station": len({r["station"] for r in rows}), "stations": sorted({r["station"] for r in rows})}
    res["evaluable"] = res["n"] >= pm.BUDGET_MIN_N and res["G_station_season"] >= pm.BUDGET_MIN_G
    if res["n"] < 5 or res["G_station_season"] < 2:
        return res
    X = [{"station": r["station"], "season": r["season"], "ho": r[hobs], "hr": r[hrim], "o1": r[o1], "m1": r[m1]} for r in rows]
    for nm, keys in (("ss", [(x["station"], x["season"]) for x in X]), ("st", [x["station"] for x in X])):
        bb, _ = pm.boot_pairs(keys, {"Q": ([x["ho"] for x in X], [x["hr"] for x in X]),
                                     "R1": ([x["o1"] for x in X], [x["m1"] for x in X])})
        with np.errstate(divide="ignore", invalid="ignore"):
            ph = (bb["Q"] - bb["R1"]) / (1.0 - bb["R1"])
        res[f"Q_ci95_{nm}"] = pm.ci95(bb["Q"])[0]
        res[f"phi_ci95_{nm}"], res[f"phi_nboot_{nm}"] = pm.ci95(ph)
        if nm == "ss":
            ok = ph[np.isfinite(ph)]
            if len(ok) > 100:
                lo, hi = np.percentile(ok, [2.5, 97.5])
                se = (hi - lo) / 3.92
                res["se_phi_quantile"] = rnd(se)
                res["sd_phi_boot"] = rnd(float(np.std(ok, ddof=1)))
                res["MDE_80"] = rnd(MDE_Z * se)
    Q = pm.ratio_of(X, "ho", "hr")
    R1 = pm.ratio_of(X, "o1", "m1")
    res.update({"Q": rnd(Q), "R1": rnd(R1), "phi": rnd((Q - R1) / (1 - R1)) if isnum(Q) and isnum(R1) and R1 != 1 else None})

    def fphi(rr):
        q, r1 = pm.ratio_of(rr, "ho", "hr"), pm.ratio_of(rr, "o1", "m1")
        return (q - r1) / (1 - r1) if isnum(q) and isnum(r1) and r1 != 1 else NAN
    res["jackknife_phi_station"] = pm.jackknife_station(X, fphi)
    return res


def budget_rows(rows, D, forcing, win):
    """E3 的事件集与淡水高：forcing 'cmorph'（_modc）或 'gauge'（_modg）；R1 同驱动同窗。"""
    out = []
    for r in rows:
        if r["group"] != p2.GROUP_A or "_obs" not in r:
            continue
        mod = (r["_modc"] if forcing == "cmorph" else r["_modg"]).get(win)
        if not mod:
            continue
        obs = r["_obs"][win]
        has20 = isnum(obs.get(20.0)) and isnum(mod.get(20.0)) and isnum(r["_s0s"].get(20.0)) and r["_s0s"][20.0] > 0
        if D in (20,) and not has20:
            continue
        if D == 25 and not (isnum(obs.get(25.0)) and isnum(mod.get(25.0))):
            continue
        w = weights_for(D, has20)
        ho, hr = h_col(obs, r["_s0s"], w), h_col(mod, r["_s0s"], w)
        if not (isnum(ho) and isnum(hr) and isnum(obs.get(1.0)) and isnum(mod.get(1.0))):
            continue
        out.append({"station": r["station"], "season": r["season"], "ho": ho, "hr": hr, "o1": obs[1.0], "m1": mod[1.0],
                    "rainmid": r["_rainmid"].get(win, NAN)})
    return out


def closure_block(brows):
    np = _np()
    X = [x for x in brows if isnum(x["rainmid"]) and x["rainmid"] > 0]
    res = {"n": len(X), "G_station_season": len({(x["station"], x["season"]) for x in X})}
    if len(X) < 5:
        return res
    res["point"] = rnd(sum(x["ho"] for x in X) * 1000 / sum(x["rainmid"] for x in X))
    bb, _ = pm.boot_pairs([(x["station"], x["season"]) for x in X],
                          {"C": ([x["ho"] * 1000 for x in X], [x["rainmid"] for x in X])})
    res["ci95"] = pm.ci95(bb["C"])[0]
    lo, hi = res["ci95"]
    res["reading"] = None if lo is None else ("闭合（CI 含 1）" if lo <= 1 <= hi else ("缺失（上端<1）" if hi < 1 else "多出（下端>1）"))
    return res


# ======================================================================== 判读（N13）
def judge_e1(env0, env5, ii):
    if env0 is None or ii is None or not ii.get("evaluable") or any(
            env0.get(k, {}).get("obs") is None or (env0.get(k, {}).get("joint95") or [None])[0] is None
            for k in ("n0", "n1", "n2", "D")):
        return {"status": "not_evaluable"}
    inside = env0["all_inside_joint"]
    robust = bool(env5 and env5["all_inside_joint"])
    hi = (ii.get("ci95") or [None, None])[1]
    ii_pass = hi is not None and hi < pm.SESOI_LO
    if inside and ii_pass:
        st, txt = "化解", "P4 失败可由选择伪影解释，错配不能解释总体高估" + ("（稳健：θ=0.5 也在包络内）" if robust else "（θ=0.5 不在包络内）")
    elif inside:
        st, txt = "部分化解", "分档形状可由伪影解释，但雨量计驱动下总体 R 不排除 ≥0.65"
    else:
        def _side(k):
            v, (lo, hi) = env0[k]["obs"], env0[k]["joint95"]
            if v is None or lo is None:
                return "不可比"
            return "高于包络" if v > hi else ("低于包络" if v < lo else "在内")
        side = {k: _side(k) for k in ("n0", "n1", "n2", "D")}
        st, txt = "未化解", f"λ 分档形状超出纯选择伪影（各量：{side}）"
    return {"status": st, "text": txt, "i_inside_theta0": inside, "i_inside_theta05": robust, "ii_pass": ii_pass}


def judge_e2(order, level, sp, v2_ok):
    if not v2_ok:
        return {"status": "not_evaluable", "order": "not_evaluable", "level": "not_evaluable",
                "reason": "V2 中性数值解与解析解不一致"}
    res = {}
    if not order.get("evaluable") or not all(c["evaluable"] for c in (order["cells"][0], order["cells"][2])):
        res["order"] = "not_evaluable"
    else:
        pts = [c["point"] for c in order["cells"]]
        mono = all(isnum(p) for p in pts) and pts[0] < pts[1] < pts[2]
        lo = order["contrast_top_minus_bottom_ci95_ss"][0]
        los = order["contrast_top_minus_bottom_ci95_st"][0]
        if mono and lo is not None and lo > 0:
            res["order"] = "可排序" if (los is not None and los > 0) else "可排序（弱）"
        else:
            res["order"] = "不可排序"
    res["rank_corr_support"] = bool(sp and sp.get("ci95") and sp["ci95"][0] is not None and sp["ci95"][0] > 0)
    ci = level.get("ci95") if level else None
    if not ci or ci[0] is None:
        res["level"] = "not_evaluable"
    elif ci[0] >= LEVEL_BAND[0] and ci[1] <= LEVEL_BAND[1]:
        res["level"] = "量级一致（二倍以内）"
    elif ci[0] > 1:
        res["level"] = "理论仍混合过强（观测淡化多于理论）"
    elif ci[1] < 1:
        res["level"] = "理论混合过弱"
    else:
        res["level"] = "不定"
    res["status"] = res["order"]
    return res


def judge_e3(main):
    b = main.get("block") or {}
    if not b.get("evaluable") or b.get("phi_ci95_ss") is None or b["phi_ci95_ss"][0] is None:
        return {"status": "not_evaluable", "D": main.get("D")}
    lo = b["phi_ci95_ss"][0]
    jk = b.get("jackknife_phi_station") or {}
    loo_min = (jk.get("loo_range") or [None])[0] if jk.get("evaluable") else None
    mde = b.get("MDE_80")
    if lo > 0:
        st = "见到淡水"
        sub = "稳健见到" if (loo_min is not None and loo_min > 0) else "见到（依赖单站）"
    elif mde is not None and mde >= MDE_BIG:
        st, sub = "功效不足", "功效不足，不构成否定"
    elif mde is not None and mde >= PHI_REF:
        st, sub = "功效不足", "功效不足以检出壁层形状量级（φ≈0.24），不构成否定（弱）"
    else:
        st, sub = "有功效而未见", "在能检出壁层量级的功效下没有见到深层淡水（对 H_V 深层去向的事后反证）"
    return {"status": st, "text": sub, "D": main["D"], "phi": b.get("phi"), "phi_ci95_ss": b["phi_ci95_ss"],
            "MDE_80": mde, "loo_min": loo_min}


TOTAL_MAP = {
    ("化解", "见到淡水"): "机制基本闭合（事后）：错配只造分档形状、不造总体高估；淡水在 A 群更深处找得到",
    ("化解", "功效不足"): "错配排除，去向未决",
    ("化解", "有功效而未见"): "错配排除，但深层未见淡水：去向偏侧向，或 RIM-3 自身多造淡水（看 M_RIM 与观测闭合比）",
    ("部分化解", "见到淡水"): "分档形状可由伪影解释、淡水去向偏深，但总体量级对雨量驱动敏感",
    ("部分化解", "功效不足"): "错配只部分排除，去向未决",
    ("部分化解", "有功效而未见"): "错配只部分排除，深层未见",
    ("未化解", "见到淡水"): "去向偏深（H_V 部分支持），但雨量错配有伪影以外的成分，s 的全球代入需要错配订正",
    ("未化解", "功效不足"): "维持 negative_L 读法",
    ("未化解", "有功效而未见"): "negative_L 读法加强",
}


def total_reading(e1, e2, e3):
    k1, k3 = e1.get("status"), e3.get("status")
    if k1 == "not_evaluable" or k3 == "not_evaluable":
        cell = f"不可评（E1={k1}，E3={k3}）；只报可评项"
    else:
        cell = TOTAL_MAP[(k1, k3)]
    notes = []
    if e2.get("reason"):
        notes.append("E2 不可评（" + e2["reason"] + "）")
    else:
        if e2["order"] == "not_evaluable":
            notes.append("E2 排序不可评")
        if e2["order"].startswith("可排序"):
            notes.append("稳定修正恢复逐事件排序（大部分非盲）")
        if e2["level"].startswith("量级一致"):
            notes.append("零参数稳定修正理论与观测量级在二倍以内")
        if not e2["order"].startswith("可排序") and not e2["level"].startswith("量级一致"):
            notes.append("稳定修正也不能解释逐事件差异")
    return {"cell": cell, "E2_notes": notes, "label": "探索性诊断",
            "first_round_verdict_unchanged": "negative_L"}


# ======================================================================== 主分析
def analyze(rows, stations, events, ref, n_sim=N_SIM, v1=None, round1=None):
    np = _np()
    S = {}
    ok = [r for r in rows if r.get("model_ok")]
    main = pm.pooled(ok, "obs1", "rim1", "V0")
    v0_ok = True if ref is None else (abs(main["point"] - ref["point"]) <= pm.V0_TOL and main["n"] == ref["n_events"]
                                      and all(abs(a - b) <= pm.V0_TOL for a, b in zip(main["ci95"], ref["ci95"])))
    S["V0"] = {"this": {k: main[k] for k in ("n", "point", "ci95")}, "pass": bool(v0_ok)}
    S["V1_round1_csv"] = v1
    if not v0_ok or (v1 is not None and not v1["pass"]):
        raise RuntimeError(f"一致性核对不过：V0={S['V0']} V1={v1}")
    # ---------------- E1
    e1 = {}
    ii = pm.pooled([r for r in ok if isnum(r.get("rimG1"))], "obs1", "rimG1", "E1(ii) 雨量计驱动 R_RIM(1 m)")
    e1["ii"] = {k: v for k, v in ii.items() if not k.startswith("_")}
    e1["ii"]["p2_sens_reference"] = {"point": 0.249, "ci95": [0.210, 0.281], "src": "p2-sens 阶段（p2_sensitivity.json）"}
    iii = pm.tercile_block([r for r in ok if isnum(r.get("rimG1"))], lambda r: r["lam"], "E1(iii) R_G by λ", den="rimG1")
    e1["iii"] = {k: v for k, v in iii.items() if not k.startswith("_")}
    obs_dp = pm.tercile_block(ok, lambda r: r["lam"], "D-P 复算")
    Rall = pm.ratio_of([r for r in ok if isnum(r["obs1"]) and isnum(r["rim1"]) and isnum(r["lam"])], "obs1", "rim1")
    pts = [c["point"] for c in obs_dp.get("cells", [])] if obs_dp.get("evaluable") else [None] * 3
    if all(isnum(p) for p in pts) and isnum(Rall) and Rall != 0:
        obs_stats = {"n0": pts[0] / Rall, "n1": pts[1] / Rall, "n2": pts[2] / Rall,
                     "D": math.log(pts[2] / pts[0]) if pts[0] > 0 and pts[2] > 0 else NAN}
    else:
        obs_stats = {"n0": NAN, "n1": NAN, "n2": NAN, "D": NAN}
    e1["observed"] = {"R_k": pts, "R_all": rnd(Rall), **{k: rnd(v) for k, v in obs_stats.items()}}
    sim_set = [r for r in ok if r.get("curve") is not None and isnum(r.get("G6sim")) and isnum(r.get("lam"))
               and isnum(r.get("rimG1"))]
    e1["sim_n_events"] = len(sim_set)
    e1["sim_events_excluded"] = len(ok) - len(sim_set)
    env = {}
    sims = {}
    if len(sim_set) >= 30:
        curves = np.array([r["curve"] for r in sim_set])
        g6 = np.array([r["G6sim"] for r in sim_set])
        lam = np.array([r["lam"] for r in sim_set])
        bins = e1_bins([r["group"] for r in sim_set], g6)
        stn = [r["station"] for r in sim_set]
        pools = {}
        for r in sim_set:
            pools.setdefault(r["station"], []).extend(r["ctrl_ds1"])
        e1["noise_pool_sizes"] = {k: len(v) for k, v in sorted(pools.items())}
        e1["bins"] = {int(b): int((bins == b).sum()) for b in sorted(set(bins.tolist()))}
        e1["g6_vs_round1_g6_max_absdiff"] = rnd(max(abs(r["G6sim"] - r["g6"]) for r in sim_set if isnum(r["g6"])), 8)
        settings = [("theta0_R0main", 0.0, R0_MAIN), ("theta05_R0main", 0.5, R0_MAIN)]
        if isnum(ii.get("point")):
            settings.append(("theta0_R0gauge", 0.0, ii["point"]))
        for name, th, r0 in settings:
            sims[name] = e1_simulate(curves, g6, lam, bins, stn, pools, r0, th, n_sim, SIM_SEED)
            env[name] = envelope(sims[name], obs_stats)
            env[name]["theta"], env[name]["R0"] = th, rnd(r0)
    e1["envelopes"] = env
    e1["judgment"] = judge_e1(env.get("theta0_R0main"), env.get("theta05_R0main"), e1["ii"])
    S["E1"] = e1
    # ---------------- E2
    e2 = {}
    vv = [r for r in ok if isnum(r.get("fdn1")) and isnum(r.get("low1"))]
    ratio = sum(r["fdn1"] for r in vv) / sum(r["low1"] for r in vv) if vv else NAN
    v2_ok = isnum(ratio) and abs(ratio - 1) <= V2_TOL
    e2["V2_neutral_fd_over_analytic"] = {"n": len(vv), "ratio": rnd(ratio), "pass": bool(v2_ok)}
    okr = [r for r in ok if isnum(r.get("rho_s")) and isnum(r.get("stab1"))]
    order = pm.tercile_block(okr, lambda r: r["rho_s"], "E2 R_RIM(1 m) by ρ_s tercile", theory=("stab1", "rim1"))
    e2["order"] = {k: v for k, v in order.items() if not k.startswith("_")}
    spr = [r for r in okr if isnum(r["obs1"])]
    e2["spearman_rho_s_vs_r_e"] = spearman_cluster([(r["station"], r["season"]) for r in spr],
                                                   [r["rho_s"] for r in spr], [r["obs1"] / r["rim1"] for r in spr]) \
        if len(spr) >= 30 else {"n": len(spr)}
    lv = pm.pooled([r for r in ok if isnum(r.get("stab1"))], "obs1", "stab1", "Λ_s＝Σobs1/Σstab1（β_h=5）")
    e2["level"] = {k: v for k, v in lv.items() if not k.startswith("_")}
    e2["level"]["R_th_stab_pooled"] = rnd(pm.ratio_of([r for r in ok if isnum(r.get("stab1")) and isnum(r["rim1"])],
                                                      "stab1", "rim1"))
    l78 = pm.pooled([r for r in ok if isnum(r.get("stab78_1"))], "obs1", "stab78_1", "β_h=7.8")
    e2["level_beta78"] = {k: v for k, v in l78.items() if not k.startswith("_")}
    o78 = [dict(r, rho78=r["stab78_1"] / r["rim1"]) for r in ok if isnum(r.get("stab78_1")) and r["rim1"] <= pm.RHO_MIN_DEN]
    t78 = pm.tercile_block(o78, lambda r: r["rho78"], "β_h=7.8 排序", theory=("stab78_1", "rim1"))
    e2["order_beta78"] = {k: v for k, v in t78.items() if not k.startswith("_")}
    lg = pm.pooled([r for r in ok if isnum(r.get("stabG1"))], "obs1", "stabG1", "雨量计驱动 Λ_s")
    e2["level_gauge"] = {k: v for k, v in lg.items() if not k.startswith("_")}
    og = [dict(r, rhoG=r["stabG1"] / r["rimG1"]) for r in ok if isnum(r.get("stabG1")) and isnum(r.get("rimG1"))
          and r["rimG1"] <= pm.RHO_MIN_DEN]
    tg = pm.tercile_block(og, lambda r: r["rhoG"], "雨量计驱动排序（R_G by ρ_G）", den="rimG1", theory=("stabG1", "rimG1"))
    e2["order_gauge"] = {k: v for k, v in tg.items() if not k.startswith("_")}
    sv = [r for r in ok if all(isnum(r.get(k)) for k in ("stab0", "rim0", "stab1", "rim1"))]
    e2["surface"] = {"n": len(sv), "stab_surface_over_rim_surface": rnd(pm.ratio_of(sv, "stab0", "rim0")),
                     "stab_1m_over_rim_1m": rnd(pm.ratio_of(sv, "stab1", "rim1")),
                     "stab_surface_over_stab_1m": rnd(pm.ratio_of(sv, "stab0", "stab1"))}
    e2["nonblind"] = {"spearman_rho_s_vs_zeta_gauge": pm.spearman([r["rho_s"] for r in okr], [r["zeta"] for r in okr]),
                      "spearman_rho_s_vs_lam": pm.spearman([r["rho_s"] for r in okr], [r["lam"] for r in okr]),
                      "spearman_rho_s_vs_rain24": pm.spearman([r["rho_s"] for r in okr], [r["rain24"] for r in okr])}
    e2["judgment"] = judge_e2(e2["order"], e2["level"], e2["spearman_rho_s_vs_r_e"], v2_ok)
    S["E2"] = e2
    # ---------------- E3
    e3 = {"blocks": {}}
    for forcing, wins in (("cmorph", WINDOWS[:2]), ("gauge", WINDOWS)):
        for w in wins:
            for D in (10, 20, 25):
                br = budget_rows(rows, D, forcing, w)
                blk = phi_block(br, "ho", "hr", "o1", "m1")
                blk["closure_obs_over_gauge"] = closure_block(br)
                e3["blocks"][f"{forcing}_{w[0]}-{w[1]}h_D{D}"] = blk
    if round1 is not None:
        b10 = e3["blocks"]["cmorph_0-6h_D10"]
        r1b = round1.get("D_B", {})
        e3["V3_phi10_vs_round1"] = {"this": {"n": b10.get("n"), "phi": b10.get("phi"), "phi_ci95_ss": b10.get("phi_ci95_ss")},
                                    "round1": {"n": r1b.get("n"), "phi": r1b.get("phi"), "phi_ci95_ss": r1b.get("phi_ci95_ss")},
                                    "match": bool(b10.get("n") == r1b.get("n") and isnum(b10.get("phi")) and isnum(r1b.get("phi"))
                                                  and abs(b10["phi"] - r1b["phi"]) <= 1e-4)}
    b25 = e3["blocks"]["cmorph_0-6h_D25"]
    mainD = 25 if b25.get("evaluable") else 20
    e3["primary"] = {"D": mainD, "rule": "φ_25 可评（n≥40 且 ≥8 簇）则为主量，否则 φ_20"}
    e3["judgment"] = judge_e3({"D": mainD, "block": e3["blocks"][f"cmorph_0-6h_D{mainD}"]})
    late = [e3["blocks"][f"gauge_{w[0]}-{w[1]}h_D{mainD}"] for w in WINDOWS[2:]]
    e3["late_freshwater_flag"] = bool(e3["judgment"].get("status") != "见到淡水" and any(
        b.get("evaluable") and b.get("phi_ci95_ss") and b["phi_ci95_ss"][0] is not None and b["phi_ci95_ss"][0] > 0
        for b in late))
    mr = {}
    for tag, sub in (("all", ok), ("A", [r for r in ok if r["group"] == p2.GROUP_A])):
        X = [r for r in sub if all(isnum(r.get(k)) for k in ("hfine100", "hfine10", "h3pt10", "pcum_cmorph_model"))
             and r["pcum_cmorph_model"] > 0]
        den = sum(r["pcum_cmorph_model"] for r in X)
        m = {"n": len(X)}
        if X and den > 0:
            m.update({"M_RIM_fine_0_100m": rnd(sum(r["hfine100"] for r in X) * 1000 / den),
                      "M_RIM_fine_0_10m": rnd(sum(r["hfine10"] for r in X) * 1000 / den),
                      "M_RIM_3pt_0_10m": rnd(sum(r["h3pt10"] for r in X) * 1000 / den)})
            bb, _ = pm.boot_pairs([(r["station"], r["season"]) for r in X],
                                  {"M": ([r["hfine100"] * 1000 for r in X], [r["pcum_cmorph_model"] for r in X])})
            m["M_RIM_fine_0_100m_ci95"] = pm.ci95(bb["M"])[0]
            v = m["M_RIM_fine_0_100m"]
            m["reading"] = "RIM-3 大致守恒" if 0.8 <= v <= 1.25 else ("RIM-3 多造淡水" if v > 1.25 else "RIM-3 少于雨量")
        mr[tag] = m
    e3["M_RIM"] = mr
    S["E3"] = e3
    S["total"] = total_reading(e1["judgment"], e2["judgment"], e3["judgment"])
    return S, sims


def check_v1(rows, csv_path):
    """N1：与第一轮 p4_mech_events.csv 逐事件核对。"""
    if not os.path.exists(csv_path):
        return {"pass": False, "reason": f"缺 {csv_path}"}
    ref = {}
    with open(csv_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            ref[(row["station"], row["onset_utc"])] = row
    n_bad, n_cmp, worst = 0, 0, 0.0
    missing = 0
    for r in rows:
        x = ref.get((r["station"], r["onset_utc"]))
        if x is None:
            missing += 1
            continue
        for k in V1_FIELDS:
            a, b = r.get(k), x.get(k, "")
            if b == "" and not isnum(a):
                continue
            try:
                bv = float(b)
            except ValueError:
                n_bad += 1
                continue
            if not isnum(a) and not math.isfinite(bv):
                continue
            n_cmp += 1
            d = abs(float(a) - bv) if isnum(a) else float("inf")
            worst = max(worst, d)
            if d > V1_TOL:
                n_bad += 1
    return {"pass": n_bad == 0 and missing == 0 and len(ref) == len(rows), "n_rows": len(rows), "n_ref": len(ref),
            "missing": missing, "n_compared": n_cmp, "n_bad": n_bad, "max_abs_diff": worst}


def run(args, out_dir, log):
    global ZFINE
    np = _np()
    ZFINE = np.concatenate([np.arange(0.0, 10.0, 0.05), np.arange(10.0, 100.0 + 1e-9, 0.5)])
    t0 = time.monotonic()
    stations, events, val, prim, plan = pm.build(args, out_dir, log)
    data, cov = pm.load_pixels(args.p2_dir, prim)
    sbn = {s["name"]: s for s in stations}
    pixc, pixm = pm.PixView(data, "c"), pm.PixView(data, "m")
    rows = [pm.event_mech(sbn[e["_stn"]], e, pixc, pixm)[0] for e in events]
    v1 = check_v1(rows, os.path.join(args.p4_dir, "p4_mech_events.csv"))
    log.log(f"V1 {json.dumps(v1, ensure_ascii=False)}", echo=True)
    if not v1["pass"]:
        p2.jdump({"V1": v1}, os.path.join(out_dir, "p4b_abort.json"))
        raise RuntimeError("V1 不过")
    grid = fd_grid()
    t1 = time.monotonic()
    ext = []
    for k, e in enumerate(events):
        ext.append(event_extra(sbn[e["_stn"]], e, rows[k], pixc, grid))
        if k % 100 == 0:
            log.log(f"event_extra {k}/{len(events)} {time.monotonic() - t1:.0f}s", echo=True)
    ref = json.load(open(os.path.join(args.p2_dir, "p2_summary.json"), encoding="utf-8"))["primary"]["R_RIM_1m"]
    r1p = os.path.join(args.p4_dir, "p4_mech_summary.json")
    round1 = json.load(open(r1p, encoding="utf-8")) if os.path.exists(r1p) else None
    S, sims = analyze(ext, stations, events, ref, v1=v1, round1=round1)
    info = {"script": "p4b_posthoc.py", "version": VERSION,
            "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "n_events": len(events),
            "rebuild_validation": {k: val[k] for k in ("n_ref", "n_rebuilt", "d4_count")}, "cmorph_pixel_cache": cov,
            "downloads": {"cache_http_requests": plan.get("cache_http_requests"), "new_data_bytes_expected": 0},
            "fd_grid_cells": len(grid[0]), "n_sim": N_SIM}
    out = dict(info)
    out.update(S)
    out["runtime_s"] = round(time.monotonic() - t0, 1)
    p2.jdump(out, os.path.join(out_dir, "p4b_summary.json"))
    with open(os.path.join(out_dir, "p4b_e1_null.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["setting", "sim", "R0", "R1", "R2", "Rall", "n0", "n1", "n2", "D"])
        for name, sm in sims.items():
            for j in range(len(sm["D"])):
                w.writerow([name, j] + [rnd(float(x), 6) for x in (*sm["R"][j], sm["Rall"][j], *sm["nk"][j], sm["D"][j])])
    cols = ["station", "group", "season", "onset_utc", "z1", "g6", "G6sim", "c6", "lam", "obs1", "rim1", "rimG1", "fdn1",
            "low1", "stab0", "stab1", "stab5", "stab10", "stab78_1", "stabG0", "stabG1", "rho_s", "hfine100", "hfine10",
            "h3pt10", "pcum_cmorph_model"]
    with open(os.path.join(out_dir, "p4b_events.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in ext:
            w.writerow({k: (rnd(r.get(k), 6) if isinstance(r.get(k), float) else r.get(k)) for k in cols})
    print(json.dumps({"E1": S["E1"]["judgment"], "E2": S["E2"]["judgment"], "E3": S["E3"]["judgment"], "total": S["total"]},
                     ensure_ascii=False), flush=True)
    return 0


# ======================================================================== 自测（合成数据，无网络）
def selftest(out_dir=None):
    global ZFINE
    np = _np()
    ZFINE = np.concatenate([np.arange(0.0, 10.0, 0.05), np.arange(10.0, 100.0 + 1e-9, 0.5)])
    res, fails = [], []
    t_start = time.monotonic()

    def check(cond, msg):
        res.append(("PASS" if cond else "FAIL") + " " + msg)
        if not cond:
            fails.append(msg)

    rng = np.random.default_rng(5)
    # 1 批量 RIM ＝ p2.rim_factor；多深度 ＝ p2.rim_factor；标量曲线插值
    Ps = np.where(rng.random((4, 84)) < 0.25, rng.gamma(2, 4, (4, 84)), 0.0)
    Um = rng.uniform(1.5, 10, 84)
    for z in (0.0, 1.0, 5.0):
        Fb = rim_factor_batch(Ps, Um, z)
        ref = np.array([p2.rim_factor(Ps[k], Um, z) for k in range(4)])
        check(np.allclose(Fb, ref, rtol=1e-12, atol=1e-14, equal_nan=True), f"rim_factor_batch＝p2.rim_factor（z={z}）")
    Fz = rim_multi_z(Ps[0], Um, np.array([0.0, 1.0, 7.3, 25.0]))
    check(all(np.allclose(Fz[q], p2.rim_factor(Ps[0], Um, z), rtol=1e-12, atol=1e-14, equal_nan=True)
              for q, z in enumerate((0.0, 1.0, 7.3, 25.0))), "rim_multi_z＝p2.rim_factor（4 个深度）")
    Pg = np.zeros(84)
    Pg[60:66] = [3.0, 6.0, 8.0, 2.0, 1.0, 0.5]
    Pg[20:22] = [1.0, 2.0]
    f = scale_curve(Pg, Um, 1.0, 35.0)
    check(abs(f[0]) < 1e-12, "标量曲线 f(0)＝0（无雨无稀释）")
    errs = []
    for a in (0.0005, 0.013, 0.37, 1.0, 2.9, 17.0):
        F = p2.rim_factor(a * Pg, Um, 1.0).reshape(-1, 2).mean(axis=1)
        direct = 35.0 * pm.delta_at(F, list(range(6)), list(range(6, 12)))
        got = float(interp_curve(f, np.array([a]))[0])
        errs.append(abs(got - direct) / max(abs(direct), 1e-4))
    check(max(errs) < 1e-3, f"标量曲线插值相对误差 <1e-3（最大 {max(errs):.2e}）")
    check(abs(float(interp_curve(f, np.array([1.0]))[0]) - 35.0 * pm.delta_at(
        p2.rim_factor(Pg, Um, 1.0).reshape(-1, 2).mean(axis=1), list(range(6)), list(range(6, 12)))) < 1e-12,
          "a=1 为精确格点")
    # 2 E1 模拟器：已知参数下复现（用解析幂律曲线代替 RIM，快速）
    A = a_grid()
    n = 600
    scale_e = rng.gamma(2.0, 0.05, n)
    curves = np.array([-s * A ** 0.8 for s in scale_e])
    g6 = rng.gamma(1.5, 8.0, n) + 0.5
    groups = ["A" if k % 3 else "B" for k in range(n)]
    stations = [f"S{k % 12}" for k in range(n)]
    lam0 = np.zeros(n)
    bins = e1_bins(groups, g6)
    check(len(set(bins.tolist())) == 10 and bins.min() >= 0, "E1 分箱：2 群×5 箱")
    pools0 = {f"S{k}": [0.0] for k in range(12)}
    s0 = e1_simulate(curves, g6, lam0, bins, stations, pools0, 0.3, 0.0, 50, 1)
    check(np.allclose(s0["R"][:, 0], 0.3, atol=1e-9) and np.allclose(s0["Rall"], 0.3, atol=1e-9),
          "E1 零误差零噪声：R_all＝R0（λ* 全同 → 全归低档，该档 R＝R0）")
    lam_true = rng.normal(0.2, 1.2, n)
    pools = {f"S{k}": list(rng.normal(0, 0.005, 40)) for k in range(12)}
    s1 = e1_simulate(curves, g6, lam_true, bins, stations, pools, 0.3, 0.0, 400, 2)
    med = np.nanmedian(s1["nk"], axis=0)
    check(med[0] < med[1] < med[2] and np.nanmedian(s1["D"]) > 0.5, f"E1 伪影：误差大时 null 形状单调上升（中位 n_k {np.round(med, 2)}）")
    # 用模拟器自身（另一 seed）生成「观测」→ 应在联合包络内
    s_obs = e1_simulate(curves, g6, lam_true, bins, stations, pools, 0.3, 0.0, 1, 99)
    ob = {"n0": s_obs["nk"][0, 0], "n1": s_obs["nk"][0, 1], "n2": s_obs["nk"][0, 2], "D": s_obs["D"][0]}
    ev_in = envelope(s1, ob)
    check(ev_in["all_inside_joint"], "E1 同过程生成的观测落在联合包络内")
    ob_out = {"n0": ob["n0"] * 0.2, "n1": ob["n1"], "n2": ob["n2"] * 4, "D": ob["D"] + 3.0}
    check(not envelope(s1, ob_out)["all_inside_joint"], "E1 形状远比伪影陡 → 包络外")
    s5 = e1_simulate(curves, g6, lam_true, bins, stations, pools, 0.3, 0.5, 400, 2)
    check(np.nanmedian(s5["D"]) < np.nanmedian(s1["D"]), "E1 θ=0.5 的 null 形状比 θ=0 平缓")
    # 3 E2 数值解
    grid = fd_grid()
    dz, zf, zc = grid
    check(abs(zf[-1] - GRID_H) < GRID_H * 0.05 and len(dz) < 400, f"网格至 {zf[-1]:.0f} m、{len(dz)} 格")
    Kc = 2e-3
    L = 8
    P = np.zeros(L)
    P[0] = 10.0                                            # 5 mm 于半步 0 中点
    c = fd_solve(P, np.full(L, 6.0), 35.0, (1.0, 5.0), K_const=Kc, grid=grid)
    t = (L - 0.5) * pm.HALF_S
    for z in (1.0, 5.0):
        ana = 0.005 / math.sqrt(math.pi * Kc * t) * math.exp(-z * z / (4 * Kc * t))
        check(abs(c[z][-1] / ana - 1) < 0.03, f"E2 常 K 对 Gaussian 解析解 z={z}（{c[z][-1]:.4e}/{ana:.4e}）")
    def _mass(Pq, Uq, **kw):
        return float(np.sum(dz * fd_solve_field(Pq, Uq, 35.0, grid=grid, **kw)))
    check(abs(_mass(P, np.full(L, 6.0), K_const=Kc) / 0.005 - 1) < 1e-6, "E2 质量守恒（常 K）")
    # 中性数值 vs 解析壁层解（常风、时变风），小时 ΔS 级
    Pw = np.zeros(84)
    Pw[60:68] = [4.0, 8.0, 10.0, 6.0, 3.0, 1.0, 0.5, 0.2]
    Pw[30:32] = [2.0, 1.0]
    for tag, Uw in (("常风 6", np.full(84, 6.0)), ("时变风", 3.0 + 6.0 * (np.arange(84) % 7) / 6.0)):
        hn = fd_hourly(Pw, Uw, 35.0, (1.0, 5.0, 10.0), beta=None, grid=grid)
        la = pm.low_hourly(Pw, Uw, (1.0, 5.0, 10.0))
        for z in (1.0, 5.0, 10.0):
            dn = pm.delta_at(-hn[z], list(range(6)), list(range(6, 12)))
            da = pm.delta_at(-la[z], list(range(6)), list(range(6, 12)))
            check(abs(dn / da - 1) < 0.03, f"E2 中性数值 vs 解析（{tag}，z={z}）：{dn:.3e}/{da:.3e}")
    check(abs(_mass(Pw, np.full(84, 6.0), beta=BETA_MAIN) / (Pw.sum() * 0.5 / 1000) - 1) < 1e-6, "E2 质量守恒（稳定修正）")
    h0 = fd_hourly(Pw, np.full(84, 6.0), 35.0, (1.0,), beta=None, grid=grid)[1.0]
    hs = fd_hourly(Pw, np.full(84, 6.0), 35.0, (1.0,), beta=BETA_MAIN, buoy_scale=1e-9, grid=grid)[1.0]
    check(np.allclose(hs, h0, rtol=1e-6, atol=1e-15), "E2 B→0：稳定修正解＝中性数值解")
    d5 = pm.delta_at(-fd_hourly(Pw, np.full(84, 4.0), 35.0, (1.0,), beta=BETA_MAIN, grid=grid)[1.0], list(range(6)), list(range(6, 12)))
    d78 = pm.delta_at(-fd_hourly(Pw, np.full(84, 4.0), 35.0, (1.0,), beta=BETA_ALT, grid=grid)[1.0], list(range(6)), list(range(6, 12)))
    dnn = pm.delta_at(-fd_hourly(Pw, np.full(84, 4.0), 35.0, (1.0,), beta=None, grid=grid)[1.0], list(range(6)), list(range(6, 12)))
    check(dnn < 0 and abs(d78) > abs(d5) > abs(dnn), f"E2 单调：|ΔS| 中性 {dnn:.2e} < β5 {d5:.2e} < β7.8 {d78:.2e}")
    check(abs(float(phi_h(0.0, 5.0)) - 1) < 1e-15 and abs(float(phi_h(-0.5, 5.0)) - 1 / 3) < 1e-12
          and abs(float(phi_h(0.2, 5.0)) - 2.0) < 1e-12, "φ_h 手算（ζ=0、−0.5、0.2）")
    # 4 MDE 公式与功效
    sd = 0.13
    draws = rng.normal(0.1, sd, 200000)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    se = (hi - lo) / 3.92
    check(abs(se / sd - 1) < 0.02, f"SE＝分位宽/3.92 复现 SD（{se:.4f}/{sd}）")
    mde = MDE_Z * se
    est = rng.normal(mde, sd, 200000)
    power = float(np.mean(est - 1.959964 * sd > 0))
    check(abs(power - 0.80) < 0.01, f"MDE 下功效≈0.80（得 {power:.3f}）")
    # 5 权重、φ 合成、闭合
    check(all(abs(sum(w.values()) - d) < 1e-12 for w, d in ((W_D[10], 10), (W_D[20], 20), (W_D["25w20"], 25), (W_D["25n20"], 25))),
          "E3 梯形权重和＝积分深度")
    brows = []
    for st in range(5):
        for se_ in range(4):
            for _ in range(6):
                m = {1.0: -0.2, 5.0: -0.06, 10.0: -0.01, 20.0: -0.001, 25.0: -0.0005}
                o = {z: 0.3 * m[1.0] * fr + rng.normal(0, 0.001) for z, fr in ((1.0, 1), (5.0, 0.8), (10.0, 0.6), (20.0, 0.3), (25.0, 0.2))}
                brows.append({"station": f"S{st}", "season": se_, "ho": h_col(o, {z: 35.0 for z in m}, W_D["25w20"]),
                              "hr": h_col(m, {z: 35.0 for z in m}, W_D["25w20"]), "o1": o[1.0], "m1": m[1.0], "rainmid": 5.0})
    pb = phi_block(brows, "ho", "hr", "o1", "m1")
    q_true = sum(W_D["25w20"][z] * 0.06 * fr for z, fr in ((1.0, 1), (5.0, 0.8), (10.0, 0.6), (20.0, 0.3), (25.0, 0.2))) / \
        sum(W_D["25w20"][z] * -v for z, v in {1.0: -0.2, 5.0: -0.06, 10.0: -0.01, 20.0: -0.001, 25.0: -0.0005}.items())
    phi_true = (q_true - 0.3) / 0.7
    check(pb["evaluable"] and abs(pb["phi"] - phi_true) < 0.02 and pb["phi_ci95_ss"][0] > 0 and pb.get("MDE_80") is not None,
          f"E3 φ_25 合成（真值 {phi_true:.3f}，得 {pb['phi']}，MDE {pb.get('MDE_80')}）")
    cb = closure_block(brows)
    check(cb["n"] == len(brows) and cb["point"] is not None, "E3 闭合比可算")
    check(abs(rain_mid_window([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], (6, 12)) - 9.0) < 1e-12, "箱中点累计（6–12 h 窗，均匀 1 mm/h → 9 mm）")
    # 6 Spearman 簇 bootstrap
    keys = [(f"S{k % 12}", k % 4) for k in range(400)]
    x = rng.normal(size=400)
    sp1 = spearman_cluster(keys, x, x + rng.normal(0, 1, 400))
    sp0 = spearman_cluster(keys, x, rng.normal(size=400))
    check(sp1["ci95"][0] > 0 and sp0["ci95"][0] < 0 < sp0["ci95"][1], "Spearman 簇 bootstrap：正相关 CI>0、零相关含 0")
    # 7 判读分支
    envd = lambda inside: {"all_inside_joint": inside, **{k: {"obs": 1.0, "joint95": [0.5, 2.0]} for k in ("n0", "n1", "n2", "D")}}
    iiok = {"evaluable": True, "ci95": [0.2, 0.3]}
    iibad = {"evaluable": True, "ci95": [0.4, 0.7]}
    check(judge_e1(envd(True), envd(True), iiok)["status"] == "化解" and "稳健" in judge_e1(envd(True), envd(True), iiok)["text"],
          "E1 判读：包络内＋(ii) 过 → 化解（稳健）")
    check(judge_e1(envd(True), envd(False), iibad)["status"] == "部分化解", "E1 判读：包络内＋(ii) 不过 → 部分化解")
    check(judge_e1(envd(False), envd(False), iiok)["status"] == "未化解", "E1 判读：包络外 → 未化解")
    env_na = {"all_inside_joint": False, **{k: {"obs": None, "joint95": [None, None]} for k in ("n0", "n1", "n2", "D")}}
    check(judge_e1(env_na, None, iiok)["status"] == "not_evaluable", "E1 判读：观测形状不可算 → not_evaluable")

    def cell(pt, ev=True):
        return {"point": pt, "evaluable": ev}
    ordp = {"evaluable": True, "cells": [cell(0.1), cell(0.2), cell(0.4)], "contrast_top_minus_bottom_ci95_ss": [0.1, 0.5],
            "contrast_top_minus_bottom_ci95_st": [0.05, 0.5]}
    ordw = dict(ordp, contrast_top_minus_bottom_ci95_st=[-0.05, 0.5])
    ordn = dict(ordp, cells=[cell(0.2), cell(0.1), cell(0.4)])
    check(judge_e2(ordp, {"ci95": [0.8, 1.6]}, {"ci95": [0.1, 0.3]}, True) ==
          {"order": "可排序", "rank_corr_support": True, "level": "量级一致（二倍以内）", "status": "可排序"}, "E2 判读：可排序＋量级一致")
    check(judge_e2(ordw, {"ci95": [1.5, 3.0]}, None, True)["order"] == "可排序（弱）"
          and judge_e2(ordw, {"ci95": [1.5, 3.0]}, None, True)["level"].startswith("理论仍混合过强"), "E2 判读：弱排序＋理论混合过强")
    check(judge_e2(ordn, {"ci95": [0.2, 0.4]}, None, True)["order"] == "不可排序"
          and judge_e2(ordn, {"ci95": [0.2, 0.4]}, None, True)["level"] == "理论混合过弱", "E2 判读：非单调 → 不可排序；理论混合过弱")
    check(judge_e2(ordp, {"ci95": [0.8, 1.6]}, None, False)["status"] == "not_evaluable", "E2 判读：V2 不过 → not_evaluable")

    def b3(lo, mde, loo=0.1):
        return {"D": 25, "block": {"evaluable": True, "phi": 0.1, "phi_ci95_ss": [lo, 0.5], "MDE_80": mde,
                                   "jackknife_phi_station": {"evaluable": True, "loo_range": [loo, 0.3]}}}
    check(judge_e3(b3(0.05, 0.3))["text"] == "稳健见到" and judge_e3(b3(0.05, 0.3, -0.1))["text"] == "见到（依赖单站）",
          "E3 判读：见到（稳健／依赖单站）")
    check(judge_e3(b3(-0.1, 0.6))["text"].startswith("功效不足，不构成否定") and judge_e3(b3(-0.1, 0.3))["text"].endswith("（弱）")
          and judge_e3(b3(-0.1, 0.2))["status"] == "有功效而未见", "E3 判读：功效不足（两档）／有功效而未见")
    check(total_reading({"status": "化解"}, {"status": "可排序", "order": "可排序", "level": "不定"},
                        {"status": "功效不足"})["cell"] == "错配排除，去向未决", "总映射：化解＋功效不足")
    check(len(TOTAL_MAP) == 9, "总映射 3×3 齐全")
    # 8 端到端：p2 两站夹具 → event_mech → event_extra → analyze（小 n_sim）
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
        ext = [event_extra(sbn[e["_stn"]], e, rows_e[k], pixc, grid) for k, e in enumerate(evs)]
        try:
            Se, _ = analyze(ext, fx["stations"], evs, None, n_sim=60)
            err = None
        except Exception as ex:
            Se, err = None, repr(ex) + traceback.format_exc()[-600:]
        check(err is None, f"端到端 analyze 无异常（{err}）")
        if Se:
            check(Se["E2"]["V2_neutral_fd_over_analytic"]["pass"], f"端到端 V2 中性数值/解析＝{Se['E2']['V2_neutral_fd_over_analytic']['ratio']}")
            check(abs((Se["E1"]["ii"].get("point") or 0) - pm.ratio_of([r for r in ext if isnum(r.get("rimG1")) and isnum(r["obs1"])],
                                                                        "obs1", "rimG1")) < 1e-4, "端到端 E1(ii) 点估计")
            bk = Se["E3"]["blocks"]["cmorph_0-6h_D10"]
            check(bk.get("n", 0) > 0 and bk.get("phi") is not None and abs(bk["phi"]) < 0.1,
                  f"端到端 E3 φ_10≈0（各深度同比，得 {bk.get('phi')}）")
            check(Se["E3"]["M_RIM"]["all"].get("n", 0) > 0 and "total" in Se, "端到端 M_RIM 与总判读字段齐")
            check(Se["E1"]["sim_n_events"] > 0 and "theta0_R0main" in Se["E1"]["envelopes"], "端到端 E1 模拟跑通")
    out = {"script": "p4b_posthoc.py", "version": VERSION, "n_checks": len(res), "n_fail": len(fails), "checks": res,
           "pass": not fails, "runtime_s": round(time.monotonic() - t_start, 1)}
    if out_dir:
        p2.jdump(out, os.path.join(out_dir, "p4b_selftest.json"))
    for line in res:
        print(line)
    print(f"自测 {len(res) - len(fails)}/{len(res)} 通过（{out['runtime_s']} s）", flush=True)
    return 0 if not fails else 4


def fd_solve_field(P, U, s0, grid=None, **kw):
    """自测用：返回末时刻的全场 c（与 fd_solve 同一求解，只多返回场）。"""
    np = _np()
    from scipy.linalg import lapack
    dz, zf, zc = grid if grid is not None else fd_grid()
    n = len(dz)
    P = np.asarray(P, float)
    U = np.asarray(U, float)
    us = pm.ustar_w(U)
    zi = zf[1:-1]
    dzc = zc[1:] - zc[:-1]
    dt = pm.HALF_S / NSUB
    c = np.zeros(n)
    beta, Kc = kw.get("beta"), kw.get("K_const")
    for j in range(len(P)):
        if Kc is not None:
            K = np.full(n - 1, float(Kc))
        else:
            K = pm.KAPPA * us[j] * zi
            if beta is not None:
                bflux = pm.G_ACC * pm.BETA_S * s0 * max(P[j], 0.0) / 1000.0 / 3600.0
                K = K / phi_h(zi * pm.KAPPA * bflux / max(us[j], 1e-12) ** 3, beta)
        g = dt * K / dzc
        d = dz.copy()
        d[:-1] += g
        d[1:] += g
        f = lapack.dgttrf(-g.copy(), d, -g.copy())
        for s in range(NSUB):
            if s == NSUB // 2:
                c[0] += P[j] * 0.5 / 1000.0 / dz[0]
            c, _ = lapack.dgttrs(f[0], f[1], f[2], f[3], f[4], dz * c)
    return c


def main(argv=None):
    ap = argparse.ArgumentParser(description="P4b 探索性诊断 E1–E3")
    ap.add_argument("--out")
    ap.add_argument("--p1-events", default=p2.P1_EVENTS_DEFAULT)
    ap.add_argument("--p1b-dir", default=p2.P1B_DIR_DEFAULT)
    ap.add_argument("--p2-dir", default=pm.P2_DIR_DEFAULT)
    ap.add_argument("--p4-dir", default=P4_DIR_DEFAULT)
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
    log = p1.Log(os.path.join(out_dir, "p4b_log.txt"))
    log.log(f"=== start {VERSION} p2_dir={args.p2_dir} p4_dir={args.p4_dir} out={out_dir}", echo=True)
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
