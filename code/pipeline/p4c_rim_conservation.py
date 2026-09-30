#!/usr/bin/env python3
"""p4c_rim_conservation.py — 模型推导：RIM-3 淡水「不守恒约 4.5 倍」是原公式固有，还是本项目实现选择造成。

只用合成输入、固定参数、秒级，无数据读取、无网络（除 uv 首次建临时环境装 xarray）。

做什么：
  C0 源码指纹：胶囊 CO2_Rain_Flux_Toolbox.py／main.py、p2_rim_test.py、p4b_posthoc.py、p4_mech.py 的 sha256。
  C1 胶囊原函数 RIMv3（xarray，直接 import；gsw／PyCO2SYS 用空模块占位——RIMv3 不用它们）在合成事件上算海面稀释因子
     F0＝sal_diluted/sal，与 p2.rim_factor(z=0) 逐值比较。
  C2 Witte et al. 式 (4)（含深度因子 exp(−z²/4Kt)）在本文件里独立逐字转写（循环写法，不共用 p2 代码），与
     p2.rim_factor 及 p4b.rim_multi_z 在 z=0.5/1/5/10 m 比较（核 p2 深度移植）。
  C3 守恒核对：∫₀^200m (1−F) dz ÷ 已落雨量，按雨后时刻 τ 列出；与线性解析式
     M_lin＝Σ c·√π·(tc/1800)·RA_true/d0 ÷ ΣRA_true 比较（Gaussian 面积恒等式使 K_z 与深度因子里的 t 全部约掉）。
  C4 实现选择归因：当前项深度 t＝1 s（K3 主路径）／1800 s（K27）／原文字面 √(K_z·1800) 三种；z 上限 100 vs 200 m；
     p4b 口径（小时平均、0–6 h 窗减雨前中位数、÷模型时刻累计雨量均值，直接调用 p4b_posthoc 函数）。
  C5 海面等价混合深度 h_eq＝已落雨量/(1−F0)（F0 取胶囊原函数），对照 d0、守恒 Gaussian 的 √(πK_z t)。

用法（约 5 s；需 xarray，可用 uv 临时环境）：
  cd code/pipeline
  uv run --no-project --with xarray --with pandas --with scipy --with numpy python p4c_rim_conservation.py [--out X.json]
依赖：numpy、scipy、xarray、pandas；同目录 witte_capsule/、p2_rim_test.py、p4b_posthoc.py、p4_mech.py（只 import，不改）。
输出：stdout 表格；--out 给出时另写 JSON。
"""
import argparse
import hashlib
import json
import math
import os
import sys
import types

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CAP_DIR = os.path.join(HERE, "witte_capsule")
sys.path.insert(0, HERE)

WINDS = [2.0, 4.0, 6.0, 8.0, 10.0]
# 合成降雨：名称 → (雨强 mm/h, 持续半步数)
SCEN = {"pulse10_30min": (10.0, 1), "const5_3h": (5.0, 6), "const10_1h": (10.0, 2), "burst40_1h": (40.0, 2)}
N_PRE, N_POST = 48, 72                     # 48 个干半步作历史，之后 36 h
TAU_H = [0.5, 1.0, 3.0, 6.0, 12.0, 23.5, 24.5, 30.0]
ZGRID = np.concatenate([np.arange(0.0, 10.0, 0.01), np.arange(10.0, 200.0 + 1e-9, 0.1)])


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def trap(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def series(rate, n):
    P = np.zeros(N_PRE + N_POST)
    P[N_PRE:N_PRE + n] = rate
    return P


# ---------------------------------------------------------------- C1 胶囊原函数
def capsule_F0(P, U):
    for m in ("gsw", "PyCO2SYS"):
        sys.modules.setdefault(m, types.ModuleType(m))   # 占位：RIMv3 不调用
    import pandas as pd
    import xarray as xr
    sys.path.insert(0, CAP_DIR)
    import CO2_Rain_Flux_Toolbox as tb
    tt = pd.date_range("2000-01-01", periods=len(P), freq="30min")
    co = {"latitude": [0.0], "longitude": [0.0], "time": tt}
    dims = ("latitude", "longitude", "time")               # main.py L19–21：time 放最后
    pr = xr.DataArray(P[None, None, :], coords=co, dims=dims)
    wd = xr.DataArray(U[None, None, :], coords=co, dims=dims)
    sl = xr.DataArray(np.full((1, 1, len(P)), 35.0), coords=co, dims=dims)
    out = tb.RIMv3(pr, wd, sl)
    return np.asarray(out.transpose("latitude", "longitude", "time").values[0, 0], float) / 35.0


# ---------------------------------------------------------------- C2 Witte et al. 式 (4) 独立转写
C1, C2, KZC = 5.7, 1 / 11, 2.5e-5
TC = np.arange(600, 1800, 25).astype(float)               # 48 个（代码）；原文称 49 元素，当前项用 1800
TI = (-np.arange(0, 24 * 3600, 1800) + 24 * 3600).astype(float)


def d0_table():
    from scipy.interpolate import RegularGridInterpolator
    DL = np.array([[1.2, 1.2, 1.0, 1.0, 1.0, 1.1, 1.1], [1.2, 1.2, 1.0, 1.0, 1.0, 1.1, 1.1],
                   [1.6, 1.6, 1.9, 2.0, 2.0, 2.2, 2.2], [1.5, 1.5, 2.5, 2.9, 3.1, 3.5, 3.5],
                   [1.9, 1.9, 2.7, 3.6, 4.2, 4.7, 4.7], [2.4, 2.4, 3.0, 4.0, 5.0, 5.8, 5.8],
                   [2.4, 2.4, 3.0, 4.0, 5.0, 5.8, 5.8]])
    f = RegularGridInterpolator(([0, 2, 4, 6, 8, 10, 200], [0, 2, 5, 10, 20, 50, 200]), DL, bounds_error=False,
                                fill_value=np.nan)
    return lambda u, r: float(f([[u, r]])[0])


D0 = d0_table()


def eq4(P, U, z, cur_mode="code"):
    """Witte et al. 式 (4) 逐项循环。cur_mode：'code'＝c2·RA/√(Kz·1 s)，深度因子 t=1 s（K3 主路径）；
    'k27'＝量级同 code、深度因子 t=1800 s；'paper'＝c2·RA/√(Kz·1800)·exp(−z²/(4Kz·1800))（原文字面）。
    返回半步 48..L−1 的 F(z)。"""
    L = len(P)
    out = np.empty(L - 48)
    for k, t in enumerate(range(48, L)):
        f = 1.0
        for i in range(48):
            j = t - 48 + i
            kz = KZC * U[j] ** 2
            d0 = D0(U[j], P[j])
            x = C1 * (P[j] / 3.6e6 * TC[i]) / math.sqrt(kz * TI[i]) * math.exp(-z * z / (4 * kz * TI[i]))
            f *= d0 / (d0 + x)
        kz = KZC * U[t] ** 2
        d0 = D0(U[t], P[t])
        ra = P[t] / 3.6e6 * 1800.0
        if cur_mode == "paper":
            x = C2 * ra / math.sqrt(kz * 1800.0) * math.exp(-z * z / (4 * kz * 1800.0))
        else:
            tdep = 1800.0 if cur_mode == "k27" else 1.0
            x = C2 * ra / math.sqrt(kz * 1.0) * math.exp(-z * z / (4 * kz * tdep))
        out[k] = f * d0 / (d0 + x)
    return out


def eq4_grid(P, U, zs, cur_mode="code"):
    """eq4 的向量化（同式，供深度积分）：返回 (nz, L−48)。"""
    L = len(P)
    zs = np.asarray(zs, float)[:, None]
    out = np.empty((zs.shape[0], L - 48))
    d0v = np.array([D0(u, r) for u, r in zip(U, P)])
    for k, t in enumerate(range(48, L)):
        j = np.arange(t - 48, t)
        kz = KZC * U[j] ** 2
        x = C1 * (P[j] / 3.6e6 * TC) / np.sqrt(kz * TI) * np.exp(-zs ** 2 / (4 * kz * TI))
        f = np.prod(d0v[j] / (d0v[j] + x), axis=1)
        kzc, ra = KZC * U[t] ** 2, P[t] / 3.6e6 * 1800.0
        if cur_mode == "paper":
            xc = C2 * ra / np.sqrt(kzc * 1800.0) * np.exp(-zs[:, 0] ** 2 / (4 * kzc * 1800.0))
        else:
            tdep = 1800.0 if cur_mode == "k27" else 1.0
            xc = C2 * ra / np.sqrt(kzc) * np.exp(-zs[:, 0] ** 2 / (4 * kzc * tdep))
        out[:, k] = f * d0v[t] / (d0v[t] + xc)
    return out


def m_lin(P, U, t):
    """线性解析：步 t 输出时 ∫(1−F)dz ≈ Σ_hist c1√π(tc/1800)RA_true/d0 + c2√π RA_true/d0（当前项，t 无关）。"""
    s = 0.0
    for i in range(48):
        j = t - 48 + i
        if P[j] > 0:
            s += C1 * math.sqrt(math.pi) * (TC[i] / 1800.0) * P[j] * 0.5 / D0(U[j], P[j])
    if P[t] > 0:
        s += C2 * math.sqrt(math.pi) * P[t] * 0.5 / D0(U[t], P[t])
    return s                                              # mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    a = ap.parse_args()
    import p2_rim_test as p2
    import p4b_posthoc as pb
    import p4_mech as pm
    R = {"sha256": {n: sha(p) for n, p in [
        ("CO2_Rain_Flux_Toolbox.py", os.path.join(CAP_DIR, "CO2_Rain_Flux_Toolbox.py")),
        ("main.py", os.path.join(CAP_DIR, "main.py")), ("p2_rim_test.py", os.path.join(HERE, "p2_rim_test.py")),
        ("p4b_posthoc.py", os.path.join(HERE, "p4b_posthoc.py")), ("p4_mech.py", os.path.join(HERE, "p4_mech.py"))]}}
    print("C0 sha256:", {k: v[:12] for k, v in R["sha256"].items()})

    # ---- C1 / C2 逐值核对
    c1max, c2max = 0.0, 0.0
    for sc, (rate, n) in SCEN.items():
        P = series(rate, n)
        for u in WINDS:
            U = np.full(len(P), u)
            fc = capsule_F0(P, U)
            c1max = max(c1max, float(np.nanmax(np.abs(fc - p2.rim_factor(P, U, 0.0)))))
            for z in (0.5, 1.0, 5.0, 10.0):
                e = eq4(P, U, z)
                c2max = max(c2max, float(np.max(np.abs(e - p2.rim_factor(P, U, z)))),
                            float(np.max(np.abs(e - pb.rim_multi_z(P, U, [z])[0]))))
    R["C1_capsule_vs_p2_z0_maxabs"] = c1max
    R["C2_eq4_vs_p2_maxabs"] = c2max
    print(f"C1 胶囊 RIMv3 海面 F vs p2.rim_factor(z=0) max|Δ|={c1max:.2e}；C2 Witte 式(4)独立转写 vs p2/p4b 深度 max|Δ|={c2max:.2e}")

    # ---- C3 守恒核对（主路径 cur_mode=code）
    print("\nC3 M(τ)＝∫0–200m(1−F)dz ÷ 已落雨量（数值｜线性解析）；τ＝雨始后小时（步末）")
    print("scen          U  " + " ".join(f"{t:>12}" for t in TAU_H))
    R["C3"] = {}
    for sc, (rate, n) in SCEN.items():
        P = series(rate, n)
        for u in WINDS:
            U = np.full(len(P), u)
            G = eq4_grid(P, U, ZGRID)
            row, cells = {}, []
            for tau in TAU_H:
                t = N_PRE + int(round(tau * 2)) - 1          # 步 t 结束于 τ
                k = t - 48
                rain = float(np.sum(P[:t + 1]) * 0.5)
                h = trap(1 - G[:, k], ZGRID) * 1000
                ml = m_lin(P, U, t)
                row[str(tau)] = {"M": h / rain, "M_lin": ml / rain}
                cells.append(f"{h / rain:5.2f}|{ml / rain:5.2f}")
            R["C3"][f"{sc}_U{u:g}"] = row
            print(f"{sc:13s} {u:3g}  " + " ".join(f"{c:>12}" for c in cells))

    # ---- C4 实现选择归因（const5_3h 为代表；τ=1、3、6 h）
    print("\nC4 当前项口径与 z 上限对 M 的影响（const5_3h；τ=3 h 雨仍在下、τ=6 h 雨停 3 h）")
    R["C4"] = {}
    P = series(5.0, 6)
    z100 = ZGRID <= 100 + 1e-9
    for u in WINDS:
        U = np.full(len(P), u)
        rec = {}
        for mode in ("code", "k27", "paper"):
            G = eq4_grid(P, U, ZGRID, mode)
            for tau in (3.0, 6.0):
                t = N_PRE + int(tau * 2) - 1
                rain = float(np.sum(P[:t + 1]) * 0.5)
                rec[f"{mode}_tau{tau:g}"] = trap(1 - G[:, t - 48], ZGRID) * 1000 / rain
                if mode == "code":
                    rec[f"code_z100_tau{tau:g}"] = trap(1 - G[z100, t - 48], ZGRID[z100]) * 1000 / rain
        R["C4"][f"U{u:g}"] = rec
        print(f"U={u:>4g}  " + "  ".join(f"{k}={v:.3f}" for k, v in rec.items()))

    # ---- C4b p4b 口径（直接调用 p4b 函数；84 半步，t0 在下标 60）
    print("\nC4b p4b 口径 M_RIM（rim_multi_z＋delta_rows＋trap，ZFINE 0–100 m；÷ cum_model_hourly）")
    zf = np.concatenate([np.arange(0.0, 10.0, 0.05), np.arange(10.0, 100.0 + 1e-9, 0.5)])
    R["C4b_p4b_metric"] = {}
    for sc, (rate, n) in SCEN.items():
        cells = []
        for u in WINDS:
            P = np.zeros(84)
            P[60:60 + n] = rate
            U = np.full(84, u)
            hz = pb.rim_multi_z(P, U, zf).reshape(len(zf), -1, 2).mean(axis=2)
            dF = pb.delta_rows(hz)
            h = pb.trap(-dF, zf) * 1000
            den = pm.cum_model_hourly(list(P[60:72] * 0.5))
            R["C4b_p4b_metric"][f"{sc}_U{u:g}"] = h / den
            cells.append(f"U{u:g}:{h / den:5.2f}")
        print(f"{sc:13s} " + "  ".join(cells))

    # ---- C5 海面等价混合深度（胶囊原函数）
    print("\nC5 pulse10_30min：h_eq＝已落雨量/(1−F0_capsule) 与 d0、√(πK_z·age) 对照（m）")
    R["C5"] = {}
    P = series(10.0, 1)
    for u in WINDS:
        U = np.full(len(P), u)
        F0 = capsule_F0(P, U)
        kz = KZC * u * u
        rec = {"d0": D0(u, 10.0)}
        for tau in (1.0, 3.0, 6.0, 12.0):
            t = N_PRE + int(tau * 2) - 1
            i_hist = 48 - (t - N_PRE)                  # 雨步在历史窗里的下标
            ti = TI[i_hist]
            heq = 5e-3 / (1 - F0[t - 48])
            rec[f"tau{tau:g}"] = {"h_eq": heq, "sqrt_pi_Kz_ti": math.sqrt(math.pi * kz * ti),
                                  "ratio_gauss_over_heq": math.sqrt(math.pi * kz * ti) / heq, "ti_s": ti}
        R["C5"][f"U{u:g}"] = rec
        print(f"U={u:>4g} d0={rec['d0']:.2f}  " + "  ".join(
            f"τ{tau:g}h: h_eq={rec[f'tau{tau:g}']['h_eq']:.2f} Gauss={rec[f'tau{tau:g}']['sqrt_pi_Kz_ti']:.2f} "
            f"比={rec[f'tau{tau:g}']['ratio_gauss_over_heq']:.2f}" for tau in (1.0, 6.0, 12.0)))

    # ---- 守恒所需参数
    R["conserving_d0_for_c1"] = C1 * math.sqrt(math.pi)
    R["conserving_c1_for_d0"] = {f"{d:g}": d / math.sqrt(math.pi) for d in (1.0, 2.0, 3.0, 5.8)}
    print(f"\n守恒要求：d0＝c1·√π·(tc/1800)＝{C1 * math.sqrt(math.pi):.2f}·w m（表内 d0 最大 5.8）；或 c1＝d0/√π")
    if a.out:
        json.dump(R, open(a.out, "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
