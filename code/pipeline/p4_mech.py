#!/usr/bin/env python3
"""p4_mech.py — 机制诊断：R_RIM(1 m)<1 来自更强的垂向混合（s 可代入全球），还是浅层困住／侧向／雨量错配。

判据 P1a–P4、总判据与出口在运行前写定（decide()）。全部结果属探索性分析，不进 P2/P3 判门。

用法：
  python p4_mech.py --selftest            # 合成数据自测（无网络、秒级）
  python p4_mech.py --plan                # 先跑自测，再重建 646 事件、核 P2 像元缓存覆盖与深层盐度可用性（只计数，
                                          # 不算任何观测变化量或比值）；下载量应为 0
  python p4_mech.py                       # 正式诊断
  --p1-events PATH / --p1b-dir DIR        同 p2_rim_test（默认 P1/P1b 产物，缓存以符号链接只读引用）
  --p2-dir DIR                            P2 主集输出目录（读 cmorph_pixels.csv 与 p2_summary.json；默认 p2 阶段输出）
  --out DIR                               本地输出目录（默认 $REPRO_OUTPUT_DIR）
依赖：numpy、scipy＋同目录 p1_events.py、p1b_extend.py、p2_rim_test.py（只 import，不改）。
数据：只读现有缓存——P1/P1b 的 MAPCO2／GTMBA／OceanSITES 缓存（A 群 SALT 文件原本就含 1–140 m 全部层，10／20／25 m
  只是 P1 未抽取）、P2 主集 cmorph_pixels.csv（中心像元＋3×3 均值，半小时）。无新下载（B 群 10 m 需扩 p1b 发现逻辑，本版不做，M3）。
产物（<out>/）：p4_mech_selftest.json；--plan：p4_mech_plan.json；正式：p4_mech_summary.json（含 verdict）、
  p4_mech_events.csv（逐事件派生量）、p4_mech_log.txt。
退出码：0 跑完（无论判定）；3 异常（含 V0 复算主结果不一致、事件重建不一致）；4 自测不过。

诊断项：
  D-U  事件 0–6 h 平均 U10 三分位上的 R_RIM(1 m)（＋K27 变体、U×雨量交叉、ln R–ln U 斜率）           → P1a
  D-Z  理论尺度①：Obukhov 长度 L=u*³/(κB)，B=g·β_S·S0·P（雨量计 0–6 h），ζ=z1/L 三分位上的 R_RIM(1 m）          → P1b
  D-T  理论尺度②：零拟合中性壁层解 K=κu*z 的逐事件预测（解析：c=Σ p/Λ·exp(−z/Λ)，Λ=κ∫u*dt），
       R_th＝ΣΔS_LoW/ΣΔS_RIM；按逐事件预测比 ρ 三分位排序检验、R_obs/R_th 下界检验（观测 1 m 淡化不少于中性
       壁层预测的一半＝「物理上合理的垂向混合足以解释 R」的必要条件）、理论隐含的海面/1 m 比（描述）          → P2a、P2b
  D-B  A 群 1/5/10 m 逐事件淡水收支（传感器深度梯形积分，顶部按 1 m 均匀），积分比 Q10、1 m 比 R1、
       深层回收率 φ=(Q10−R1)/(1−R1)（>0＝观测淡水相对 RIM-3 更深）；形状比与壁层解参照、闭合比（观测/模型/
       壁层解 ÷ 雨量）、20/25 m 延伸（描述）                                                                → P3
  D-P  点-像元：雨量计/CMORPH 0–6 h 雨量对数比 λ 三分位上的 R_RIM(1 m)；3×3 均值驱动的 R（描述）              → P4
  D-S  0.5 m 与 1 m 按同采样小时配对的倒置比 I_m、1 m 按 MAPCO2 采样小时重构的 R、模型时间轴
       后移 30/60 min 的 R、逐站同小时「同传感器样」比例（描述，不进总判据）
  V0   复算 R_RIM(1 m) 点估计与 95% CI，须与 P2 p2_summary.json 在 1e-4 内一致，否则退出 3（验证事件、模型与抽样实现）

实现选择（M 条）：
  M1 事件集＝P2 主集 646 事件（p2.build_all 重建并逐条校验 K1）；各诊断的分析集＝该诊断所需观测与模型都有效的事件（同 K10）。
  M2 模型 RIM-3 全部沿用 p2（K2–K9、K28 实际「1 m」层深、K3 当前项 t=1 s）；新增深度 0 m（S0 取 1 m 层雨前中位数）与 10 m
     （S0 取 10 m 雨前中位数）只用于 D-T／D-B。
  M3 深层盐度：A 群 GTMBA SALT 缓存文件（与 P1 同一批文件、同 QC {1,2}、同小时化 I4）抽 10、20、25 m；B 群无 10 m（不下载）。
     D-B 只用 A 群；0–5 m 版收支（A＋B）只作描述。
  M4 U10＝浮标风 ×1.0865（K9，含 ≤6 h 插补与 0.1 m/s 下限），事件量＝[t0, t0+6 h) 12 个半步均值。
  M5 u*_w＝√(ρa·Cd/ρw)·U10，Cd＝10⁻³(2.7/U+0.142+0.0764U)（Large & Yeager 2009，U 下限 0.5 m/s 只用于 Cd），
     ρa=1.22、ρw=1025 kg m⁻³，κ=0.4。不含破波增强（Craig & Banner 1994）与层结修正——D-T 是「中性壁层」零参数基准。
  M6 壁层解：K(z,t)=κu*(t)z 的面源解 c(z,t)=Σ_j p_j/Λ_j·exp(−z/Λ_j)，Λ_j=κ∫_{t_j}^{t}u*dt'（半步中点入水，半步末取值；
     对时变 u* 精确）；雨量＝与 RIM 同一驱动（D-T 用 CMORPH 中心像元，D-B 闭合比另用雨量计）；小时值＝两个半步末值之均（K5）。
  M7 Obukhov：L=u*_w³/(κ·g·β_S·S0·P)，g=9.81、β_S=7.5×10⁻⁴（g/kg）⁻¹（常数，只用于排序与三分位，不改秩），
     u*_w 取 [t0,t0+6 h) 均值，P＝雨量计 0–6 h 平均雨强（负值置 0）；P=0 时 ζ=0。忽略雨致降温与日射加热的浮力。
  M8 ρ_e＝ΔS_LoW,e(z1)/ΔS_RIM,e(z1)（均 CMORPH）；只对 ΔS_RIM,e(z1) ≤ −0.005 psu 的事件定义（排序键，防分母近 0）。
  M9 收支积分（观测与模型同一套）：深度 {1,5,10} m，0–1 m 取 1 m 值，梯形 → 权重 3/4.5/2.5；淡水高 h＝Σw·(−ΔS_z/S0_z)。
     Q10＝Σh_obs/Σh_RIM；R1 与 φ 在同一事件集、同一次重抽中算。
  M10 雨量时间匹配：观测小时箱 k（k=0..5）对应累计雨量 Σ_{j<k}r_j＋r_k/2（箱中点）；模型/壁层解按其半步末时刻累计。
  M11 λ＝ln((G6+0.1)/(C6+0.1))，G6 雨量计 [t0,t0+6 h) 累计（负值置 0），C6 CMORPH 中心像元同窗累计（mm）。
  M12 三分位：在该诊断分析集的事件上按 numpy.quantile(1/3, 2/3) 切分（≤q1、≤q2、其余），切点固定、不随重抽变。
  M13 推断：站×季整簇 bootstrap（p2.cluster_boot：簇键排序、PCG64(20260926)、B=10000、Σnum/Σden、百分位）；同一诊断内
     所有量用同一次抽样（对比、比值之比、φ 逐次计算）；另报以站为簇的 bootstrap（少簇替代）与逐站剔除 jackknife t 区间。
     「弱通过」：主 CI 满足而少簇替代不满足——三分位类用站簇 bootstrap；P3（只约 5 站）用逐站剔除 φ 的最小值 >0。
     单格可评条件：n≥30 事件且 ≥8 簇；D-B 可评：n≥40 且 ≥8 簇。
  M14 时间后移：模型半步序列前补 k=1、2 个半步（CMORPH 缓存有该小时即用，否则置 0，计数），小时值改用前移 k 个半步的输出。
  M15 配对采样：同一事件前 6 h 与 0–6 h 窗内，只取 0.5 m 与 1 m 同时有效的小时，按 I8 求两层 ΔS；模型同样只取这些小时。
  M16 同传感器样（小时级）：|SSS_MAPCO2 − S_1m层| ≤ 0.003 的同小时比例，逐站、全时段与事件窗 [t0−6, t0+6) 各报。
  M17 判定、出口与可评规则全部写在 decide()；不做多重比较校正（总判据为交集-并检验，各支柱 α=0.05）。

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

VERSION = "p4mech-2026-09-27a"
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
P2_DIR_DEFAULT = _rp.upstream("p2_dir")  # [repro] 读 p2 阶段输出目录
NAN = float("nan")

# ---- 物理常数（M5、M7） ----
KAPPA = 0.4
RHO_A, RHO_W = 1.22, 1025.0
G_ACC = 9.81
BETA_S = 7.5e-4
CD_UMIN = 0.5
HALF_S = 1800.0

# ---- 判据（运行前写定） ----
SESOI_LO = 0.65
PHI_LAT = 0.10                                  # 否定信号：φ 的 95% CI 上端 <0.10（深层几乎无回收）
RTH_MIN = 0.5                                   # P2b：R_obs/R_th 的 95% CI 下端 ≥0.5
RHO_MIN_DEN = -0.005
CELL_MIN_N, CELL_MIN_G = 30, 8
BUDGET_MIN_N, BUDGET_MIN_G = 40, 8
V0_TOL = 1e-4

# ---- 实现层 ----
B = p2.BOOT_B
SEED = p2.BOOT_SEED
DEEP_TARGETS = (10.0, 20.0, 25.0)
W_10 = {1.0: 3.0, 5.0: 4.5, 10.0: 2.5}            # M9：0–10 m
W_5 = {1.0: 3.0, 5.0: 2.0}                        # 0–5 m（描述）
SITE2NAME = {st["site"]: st["name"] for st in p1.STATIONS}
I_T0 = -p2.WIN_LO_H * 2                           # 强迫半步序列中 t0 的下标（=60）
_DEEP = {}


def _np():
    import numpy
    return numpy


def isnum(x):
    return x is not None and isinstance(x, (int, float)) and math.isfinite(x)


def rnd(x, nd=5):
    return None if not isnum(x) else round(float(x), nd)


# ======================================================================== 物理（M5–M7）
def drag_coef(u):
    np = _np()
    uu = np.maximum(np.asarray(u, float), CD_UMIN)
    return 1e-3 * (2.7 / uu + 0.142 + 0.0764 * uu)


def ustar_w(u):
    np = _np()
    u = np.asarray(u, float)
    return np.sqrt(RHO_A * drag_coef(u) / RHO_W) * u


def obukhov_zeta(z1, us_mean, s0, p_mm_h):
    """ζ=z1/L，L=u*³/(κB)，B=g·β_S·S0·P（P m/s）。P≤0 → 0。"""
    if not (isnum(us_mean) and isnum(s0) and isnum(p_mm_h)) or us_mean <= 0:
        return NAN
    if p_mm_h <= 0:
        return 0.0
    bflux = G_ACC * BETA_S * s0 * p_mm_h / 1000.0 / 3600.0
    return z1 * KAPPA * bflux / us_mean ** 3


def low_conc(P, U, zs):
    """壁层解（M6）。P mm/h、U m/s 半步序列（长 L）。返回 {z: 淡水体积分数 c 在各半步末（长 L）}。"""
    np = _np()
    P = np.asarray(P, float)
    us = ustar_w(np.asarray(U, float))
    L = len(P)
    cum = np.concatenate([[0.0], np.cumsum(us * HALF_S)])          # cum[n+1]＝前 n 个半步积分到半步 n 末
    dep = P * 0.5 / 1000.0                                          # 每半步入水（m）
    n_idx = np.arange(L)[:, None]
    j_idx = np.arange(L)[None, :]
    lam = KAPPA * (cum[n_idx + 1] - cum[j_idx + 1] + 0.5 * us[j_idx] * HALF_S)
    mask = j_idx <= n_idx
    out = {}
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for z in zs:
            k = np.where(mask & (lam > 0), dep[None, :] / lam * np.exp(-z / lam), 0.0)
            out[z] = k.sum(axis=1)
    return out


def low_hourly(P, U, zs):
    """壁层解小时值：半步 48..83 两两平均（与 RIM 小时值 K5 同位）→ 18 个小时（t0−6…t0+11）。"""
    c = low_conc(P, U, zs)
    return {z: v[48:].reshape(-1, 2).mean(axis=1) for z, v in c.items()}


# ======================================================================== 深层盐度（M3）
def _wrap_loader(orig):
    def f(http, cache, base, kind, per_hour, log, qc_drop):
        out = orig(http, cache, base, kind, per_hour, log, qc_drop)
        if kind == "SALT":
            fname = base.rsplit("/", 1)[-1]
            site = fname.split("_")[1]
            name = SITE2NAME.get(site)
            apath = os.path.join(cache, "thredds", site, fname + ".ascii.gz")
            if name and os.path.exists(apath):
                add_deep_from_ascii(name, apath, per_hour)
        return out
    return f


def add_deep_from_ascii(name, apath, per_hour):
    arr = p1.parse_dap_ascii(apath, {"TIME", "DEPTH", "PSAL", "PSAL_QC"})
    times, depths = arr.get("TIME"), list(arr.get("DEPTH", []))
    vals, qcs = arr.get("PSAL"), arr.get("PSAL_QC")
    if not times or not depths or vals is None or qcs is None:
        return
    nt, nd = len(times), len(depths)
    if len(vals) != nt * nd or len(qcs) != nt * nd:
        return
    for target in DEEP_TARGETS:
        if target in depths:
            j = depths.index(target)
            pairs = p1.file_to_hourly(times, vals[j::nd], qcs[j::nd], per_hour, {}, 0.0, 45.0)
            acc = _DEEP.setdefault(name, {}).setdefault(f"s{int(target)}", {})
            for h, v in pairs:
                s = acc.setdefault(h, [0.0, 0])
                s[0] += v
                s[1] += 1


def attach_deep(stations):
    np = _np()
    for stn in stations:
        acc = _DEEP.get(stn["name"], {})
        for key in ("s10", "s20", "s25"):
            a = np.full(stn["n"], np.nan)
            for h, (s, c) in acc.get(key, {}).items():
                i = h - stn["h0"]
                if 0 <= i < stn["n"] and c:
                    a[i] = s / c
            stn["ser"][key] = a


# ======================================================================== 事件与派生量
class PixView:
    """cmorph_pixels.csv 只读视图；mode='c' 中心像元（v0,v1）、'm' 3×3 均值（m0,m1）。"""

    def __init__(self, data, mode="c"):
        self.data, self.mode = data, mode

    def get(self, h, st):
        r = self.data.get(h)
        if r is None:
            return None
        v = r.get(st)
        if v is None:
            return None
        return (v[0], v[1]) if self.mode == "c" else (v[2], v[3])


def rim_hourly(P, U, z, shift=0):
    """RIM 稀释因子小时值（18 个）；shift=k：P/U 前已补 k 个半步，小时值取前移 k 个半步的输出（M14）。"""
    F = p2.rim_factor(P, U, z)
    return F[:36].reshape(-1, 2).mean(axis=1) if shift else F.reshape(-1, 2).mean(axis=1)


def padded(P, U, k, stn, H, pix):
    np = _np()
    pre, n0 = [], 0
    for m in range(k, 0, -1):                          # 半步 −m：小时 H−31（k≤2 时都在该小时内）
        r = pix.get(H + p2.WIN_LO_H - 1, stn["name"])
        v = r[1] if (r and m == 1) else (r[0] if r else None)
        if v is None or not math.isfinite(v):
            v, n0 = 0.0, n0 + 1
        pre.append(v)
    return np.concatenate([pre, P]), np.concatenate([[U[0]] * k, U]), n0


def delta_at(arr_hourly, idx_pre, idx_post):
    np = _np()
    a = np.asarray(arr_hourly, float)
    pre, post = a[idx_pre], a[idx_post]
    if not len(pre) or not len(post) or not (np.all(np.isfinite(pre)) and np.all(np.isfinite(post))):
        return NAN
    return float(np.mean(post) - np.median(pre))


def cum_mid_hourly(r_hourly):
    """M10：观测小时箱 k 的累计雨量 Σ_{j<k} r_j + r_k/2，返回 6 个值的均值（mm）。"""
    tot, acc = 0.0, []
    for x in r_hourly:
        acc.append(tot + 0.5 * x)
        tot += x
    return sum(acc) / len(acc)


def cum_model_hourly(p_half_mm):
    """模型/壁层解时刻（半步末，两两平均）的累计雨量均值；p_half_mm 为 t0 起 12 个半步的入水量（mm）。"""
    tot, ends = 0.0, []
    for x in p_half_mm:
        tot += x
        ends.append(tot)
    hrs = [(ends[2 * k] + ends[2 * k + 1]) / 2 for k in range(6)]
    return sum(hrs) / 6


def event_mech(stn, e, pixc, pixm):
    """一个事件的全部 P4 派生量（模型与观测）。观测量只在正式模式调用。"""
    np = _np()
    rec = p2.event_record(stn, e, pixc, ("cmorph",), (p2.PRIMARY_WIN,), k27=True)
    p2.fill_obs(rec, stn)
    ser, i, H = stn["ser"], e["i"], e["hour"]
    z1 = rec["s1_depth_m"]
    out = {"station": stn["name"], "group": stn["group"], "season": e["season"], "onset_utc": e["onset_utc"],
           "z1": z1, "rain24": e["acc24"],
           "obs1": rec["obs"].get((1.0, p2.PRIMARY_WIN), NAN), "obs05": rec["obs"].get((0.5, p2.PRIMARY_WIN), NAN),
           "obs5": rec["obs"].get((5.0, p2.PRIMARY_WIN), NAN),
           "rim1": rec["model"][("cmorph", "rim", 1.0, p2.PRIMARY_WIN)],
           "rim05": rec["model"][("cmorph", "rim", 0.5, p2.PRIMARY_WIN)],
           "rim5": rec["model"][("cmorph", "rim", 5.0, p2.PRIMARY_WIN)],
           "rim1_k27": rec["model"].get((p2.K27_TAG, "rim", 1.0, p2.PRIMARY_WIN), NAN)}
    s0 = {0.5: rec["s0"][0.5], 1.0: rec["s0"][1.0], 5.0: rec["s0"][5.0], 10.0: p2.pre_median(ser["s10"], i),
          20.0: p2.pre_median(ser["s20"], i), 25.0: p2.pre_median(ser["s25"], i)}
    for z in (10.0, 20.0, 25.0):
        out[f"obs{int(z)}"] = p2.obs_delta(ser[f"s{int(z)}"], i, p2.PRIMARY_WIN)
    pre_i, post_i = list(range(0, 6)), list(range(6, 12))
    P, U, diag = p2.event_forcing(stn, H, "cmorph", pixc)
    if P is None:
        out["model_ok"] = False
        return out, rec
    out["model_ok"] = True
    zm = {0.0: z1, 1.0: z1, 5.0: 5.0, 10.0: 10.0, 20.0: 20.0, 25.0: 25.0}
    # RIM 其余深度（0、10、20、25 m）
    for zl, zz in ((0.0, 0.0), (10.0, 10.0), (20.0, 20.0), (25.0, 25.0)):
        s0z = s0[1.0] if zl == 0.0 else s0[zl]
        d = delta_at(rim_hourly(P, U, zz), pre_i, post_i)
        out[f"rim{int(zl)}"] = s0z * d if isnum(s0z) and isnum(d) else NAN
    # 壁层解（CMORPH 驱动）
    lh = low_hourly(P, U, (0.0, z1, 5.0, 10.0))
    for zl, zz in ((0.0, 0.0), (1.0, z1), (5.0, 5.0), (10.0, 10.0)):
        s0z = s0[1.0] if zl == 0.0 else s0[zl]
        d = delta_at(-lh[zz], pre_i, post_i)
        out[f"low{int(zl)}"] = s0z * d if isnum(s0z) and isnum(d) else NAN
    # 风、u*、雨量
    u_ev = U[I_T0:I_T0 + 12]
    out["u10"] = float(np.mean(u_ev))
    us_ev = float(np.mean(ustar_w(u_ev)))
    r6 = np.asarray(ser["rain"][i:i + 6], float)
    r6 = np.where(np.isfinite(r6), np.maximum(r6, 0.0), np.nan)
    out["g6"] = float(np.sum(r6)) if np.all(np.isfinite(r6)) else NAN
    out["c6"] = float(np.sum(P[I_T0:I_T0 + 12]) * 0.5)
    out["zeta"] = obukhov_zeta(z1, us_ev, s0[1.0], out["g6"] / 6.0 if isnum(out["g6"]) else NAN)
    out["lam"] = math.log((out["g6"] + 0.1) / (out["c6"] + 0.1)) if isnum(out["g6"]) else NAN
    out["rho"] = (out["low1"] / out["rim1"]) if (isnum(out["low1"]) and isnum(out["rim1"]) and out["rim1"] <= RHO_MIN_DEN) \
        else NAN
    # 雨量计驱动的壁层解（收支闭合比）与雨量时间匹配
    Pg, Ug, _ = p2.event_forcing(stn, H, "gauge", pixc)
    if Pg is not None:
        lg = low_hourly(Pg, Ug, (z1, 5.0, 10.0))
        for zl, zz in ((1.0, z1), (5.0, 5.0), (10.0, 10.0)):
            d = delta_at(-lg[zz], pre_i, post_i)
            out[f"lowg{int(zl)}"] = s0[zl] * d if isnum(s0[zl]) and isnum(d) else NAN
        out["pcum_gauge_model"] = cum_model_hourly(list(Pg[I_T0:I_T0 + 12] * 0.5))
    out["pcum_gauge_obs"] = cum_mid_hourly(list(r6)) if np.all(np.isfinite(r6)) else NAN
    out["pcum_cmorph_model"] = cum_model_hourly(list(P[I_T0:I_T0 + 12] * 0.5))
    for z in (1.0, 5.0, 10.0, 20.0, 25.0):
        out[f"s0_{int(z)}"] = s0[z]
    # 3×3 均值驱动（D-P 描述）
    Pm, Um, _ = p2.event_forcing(stn, H, "cmorph", pixm)
    if Pm is not None:
        d = delta_at(rim_hourly(Pm, Um, z1), pre_i, post_i)
        out["rim1_m3"] = s0[1.0] * d if isnum(s0[1.0]) and isnum(d) else NAN
    # 时间后移（M14）
    for k in (1, 2):
        Pp, Up, nz = padded(P, U, k, stn, H, pixc)
        d = delta_at(rim_hourly(Pp, Up, z1, shift=k), pre_i, post_i)
        out[f"rim1_shift{k}"] = s0[1.0] * d if isnum(s0[1.0]) and isnum(d) else NAN
        out[f"shift{k}_pad_zero"] = nz
    # 配对采样（M15）
    s05, s1 = ser["sss05"], ser["s1"]
    pre_h = [h for h in range(max(0, i - 6), i) if math.isfinite(s05[h]) and math.isfinite(s1[h])]
    post_h = [h for h in range(i, i + 6) if math.isfinite(s05[h]) and math.isfinite(s1[h])]
    if pre_h and post_h:
        out["obs05_m"] = float(np.mean(s05[post_h]) - np.median(s05[pre_h]))
        out["obs1_m"] = float(np.mean(s1[post_h]) - np.median(s1[pre_h]))
        Fh = rim_hourly(P, U, z1)
        out["rim1_m"] = s0[1.0] * delta_at(Fh, [h - i + 6 for h in pre_h], [h - i + 6 for h in post_h]) \
            if isnum(s0[1.0]) else NAN
    # E5 代理（非盲程度量化用）
    out["e5proxy"] = (out["rim1"] / out["rim05"]) if (isnum(out["rim1"]) and isnum(out["rim05"])
                                                      and out["rim05"] <= -0.005 and abs(z1 - 1.0) < 1e-9) else NAN
    return out, rec


# ======================================================================== 统计（M12、M13）
def ci95(a):
    np = _np()
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if len(a) < 100:
        return [None, None], int(len(a))
    return [rnd(np.percentile(a, 2.5)), rnd(np.percentile(a, 97.5))], int(len(a))


def boot_pairs(keys, pairs):
    """p2.cluster_boot 的封装：返回 ({名: 比值重抽数组}, 簇数)。"""
    return p2.cluster_boot(keys, pairs, B, SEED)


def jackknife_station(rows, fn):
    """逐站剔除：fn(rows)→标量。返回 {points, se, t_ci95, G}。"""
    from scipy import stats
    np = _np()
    sts = sorted({r["station"] for r in rows})
    G = len(sts)
    if G < 3:
        return {"G": G, "evaluable": False}
    th = [fn([r for r in rows if r["station"] != s]) for s in sts]
    th = np.array([x if isnum(x) else np.nan for x in th])
    ok = th[np.isfinite(th)]
    if len(ok) < 3:
        return {"G": G, "evaluable": False}
    full = fn(rows)
    se = math.sqrt((len(ok) - 1) / len(ok) * float(np.sum((ok - ok.mean()) ** 2)))
    tq = float(stats.t.ppf(0.975, len(ok) - 1))
    return {"G": G, "evaluable": True, "loo_range": [rnd(float(ok.min())), rnd(float(ok.max()))], "se": rnd(se),
            "t_ci95": [rnd(full - tq * se), rnd(full + tq * se)]}


def ratio_of(rows, num, den):
    n = sum(r[num] for r in rows)
    d = sum(r[den] for r in rows)
    return n / d if d else NAN


def pooled(rows, num, den, label):
    """合并比值＋站×季与站两种簇的 CI＋jackknife。"""
    rows = [r for r in rows if isnum(r.get(num)) and isnum(r.get(den))]
    res = {"label": label, "n": len(rows)}
    if len(rows) < 2:
        res["evaluable"] = False
        return res
    pairs = {"r": ([r[num] for r in rows], [r[den] for r in rows])}
    b1, g1 = boot_pairs([(r["station"], r["season"]) for r in rows], pairs)
    b2, g2 = boot_pairs([r["station"] for r in rows], pairs)
    res.update({"point": rnd(ratio_of(rows, num, den)), "G_station_season": g1, "G_station": g2,
                "ci95": ci95(b1["r"])[0], "ci95_station": ci95(b2["r"])[0],
                "jackknife_station": jackknife_station(rows, lambda rr: ratio_of(rr, num, den)),
                "evaluable": True, "_boot": b1["r"], "_boot_st": b2["r"]})
    return res


def terciles(rows, keyf):
    np = _np()
    v = np.array([keyf(r) for r in rows], float)
    q1, q2 = np.quantile(v, [1 / 3, 2 / 3])
    lab = np.where(v <= q1, 0, np.where(v <= q2, 1, 2))
    return lab, [float(q1), float(q2)], v


def tercile_block(rows, keyf, label, num="obs1", den="rim1", theory=None, loglog=False):
    """三分位合并比值、上－下对比（同一次抽样）、可评性；theory＝(num,den) 另报理论比值（如 low1/rim1）。"""
    np = _np()
    rows = [r for r in rows if isnum(r.get(num)) and isnum(r.get(den)) and isnum(keyf(r))
            and (theory is None or (isnum(r.get(theory[0])) and isnum(r.get(theory[1]))))]
    res = {"label": label, "n": len(rows)}
    if len(rows) < 9:
        res["evaluable"] = False
        return res
    lab, cuts, v = terciles(rows, keyf)
    res["cuts"] = [rnd(c, 6) for c in cuts]
    out = {}
    for kk_name, keys in (("ss", [(r["station"], r["season"]) for r in rows]), ("st", [r["station"] for r in rows])):
        pairs = {}
        for k in range(3):
            m = lab == k
            pairs[f"t{k}"] = ([r[num] if m[j] else 0.0 for j, r in enumerate(rows)],
                              [r[den] if m[j] else 0.0 for j, r in enumerate(rows)])
            if theory:
                pairs[f"th{k}"] = ([r[theory[0]] if m[j] else 0.0 for j, r in enumerate(rows)],
                                   [r[theory[1]] if m[j] else 0.0 for j, r in enumerate(rows)])
        out[kk_name] = boot_pairs(keys, pairs)
    cells = []
    for k in range(3):
        sub = [r for j, r in enumerate(rows) if lab[j] == k]
        g = len({(r["station"], r["season"]) for r in sub})
        cell = {"k": k, "n": len(sub), "G_station_season": g, "G_station": len({r["station"] for r in sub}),
                "key_median": rnd(float(np.median(v[lab == k])), 6) if sub else None, "point": rnd(ratio_of(sub, num, den)),
                "ci95": ci95(out["ss"][0][f"t{k}"])[0], "ci95_station": ci95(out["st"][0][f"t{k}"])[0],
                "evaluable": len(sub) >= CELL_MIN_N and g >= CELL_MIN_G}
        if theory:
            cell["theory_point"] = rnd(ratio_of(sub, theory[0], theory[1]))
        cells.append(cell)
    res["cells"] = cells
    for nm, (bb, _g) in (("ss", out["ss"]), ("st", out["st"])):
        d = bb["t2"] - bb["t0"]
        res[f"contrast_top_minus_bottom_ci95_{nm}"] = ci95(d)[0]
    res["contrast_top_minus_bottom_point"] = rnd(cells[2]["point"] - cells[0]["point"]) \
        if isnum(cells[2]["point"]) and isnum(cells[0]["point"]) else None
    if loglog and all(isnum(c["key_median"]) and c["key_median"] > 0 for c in cells):
        lx = np.log([c["key_median"] for c in cells])
        bb = out["ss"][0]
        with np.errstate(divide="ignore", invalid="ignore"):
            ly = np.log(np.stack([bb["t0"], bb["t1"], bb["t2"]], axis=1))
        xm = lx - lx.mean()
        slope = (ly - ly.mean(axis=1, keepdims=True)) @ xm / float(xm @ xm)
        pts = [c["point"] for c in cells]
        res["loglog_slope"] = {"point": rnd(float(np.log(pts) @ xm / float(xm @ xm))) if all(isnum(p) and p > 0 for p in pts)
                               else None, "ci95": ci95(slope)[0]}
        if theory:
            tp = [c["theory_point"] for c in cells]
            res["loglog_slope_theory"] = rnd(float(np.log(tp) @ xm / float(xm @ xm))) \
                if all(isnum(p) and p > 0 for p in tp) else None
    res["evaluable"] = True
    return res


def budget_block(rows):
    """D-B（M9）：Q10、R1、φ 同一次抽样；闭合比与 0–5 m、20/25 m 版描述。"""
    np = _np()

    def h(r, pref, wts):
        s = 0.0
        for z, w in wts.items():
            dz, s0z = r.get(f"{pref}{int(z)}"), r.get(f"s0_{int(z)}")
            if not (isnum(dz) and isnum(s0z)) or s0z <= 0:
                return NAN
            s += w * (-dz / s0z)
        return s

    A = []
    for r in rows:
        if r["group"] != p2.GROUP_A or not r.get("model_ok"):
            continue
        x = dict(r)
        x["h_obs"], x["h_rim"] = h(r, "obs", W_10), h(r, "rim", W_10)
        x["h_lowg"] = h(r, "lowg", W_10)
        if isnum(x["h_obs"]) and isnum(x["h_rim"]) and isnum(r["obs1"]) and isnum(r["rim1"]):
            A.append(x)
    res = {"n": len(A), "G_station_season": len({(r["station"], r["season"]) for r in A}),
           "G_station": len({r["station"] for r in A}), "stations": sorted({r["station"] for r in A})}
    res["evaluable"] = res["n"] >= BUDGET_MIN_N and res["G_station_season"] >= BUDGET_MIN_G
    if res["n"] < 5:
        return res
    for nm, keys in (("ss", [(r["station"], r["season"]) for r in A]), ("st", [r["station"] for r in A])):
        bb, _ = boot_pairs(keys, {"Q": ([r["h_obs"] for r in A], [r["h_rim"] for r in A]),
                                  "R1": ([r["obs1"] for r in A], [r["rim1"] for r in A])})
        with np.errstate(divide="ignore", invalid="ignore"):
            phi = (bb["Q"] - bb["R1"]) / (1.0 - bb["R1"])
        res[f"Q10_ci95_{nm}"] = ci95(bb["Q"])[0]
        res[f"R1_ci95_{nm}"] = ci95(bb["R1"])[0]
        res[f"phi_ci95_{nm}"], res[f"phi_nboot_{nm}"] = ci95(phi)
    Q = ratio_of(A, "h_obs", "h_rim")
    R1 = ratio_of(A, "obs1", "rim1")
    res.update({"Q10": rnd(Q), "R1": rnd(R1), "phi": rnd((Q - R1) / (1 - R1)) if isnum(Q) and isnum(R1) and R1 != 1 else None,
                "jackknife_phi_station": jackknife_station(
                    A, lambda rr: (ratio_of(rr, "h_obs", "h_rim") - ratio_of(rr, "obs1", "rim1"))
                    / (1 - ratio_of(rr, "obs1", "rim1")))})
    # 描述：形状比 Q/R1（观测）与壁层解同口径（CMORPH 驱动 low1/5/10 当作「观测」）的理论参照
    res["shape_obs_Q_over_R1"] = rnd(Q / R1) if isnum(Q) and isnum(R1) and R1 else None
    lw = [dict(r, h_low=h(r, "low", W_10)) for r in A]
    lw = [r for r in lw if isnum(r["h_low"]) and isnum(r.get("low1"))]
    if lw:
        ql, rl = ratio_of(lw, "h_low", "h_rim"), ratio_of(lw, "low1", "rim1")
        res["shape_low_Q_over_R1"] = rnd(ql / rl) if isnum(ql) and isnum(rl) and rl else None
        res["phi_if_obs_had_low_shape"] = rnd(R1 * (ql / rl - 1) / (1 - R1)) if isnum(ql) and isnum(rl) and rl and R1 != 1 else None
    # 描述：闭合比（mm 对 mm）
    gg = [r for r in A if isnum(r.get("pcum_gauge_obs")) and r["pcum_gauge_obs"] > 0]
    res["closure_obs_over_gauge"] = rnd(sum(r["h_obs"] for r in gg) * 1000 / sum(r["pcum_gauge_obs"] for r in gg)) if gg else None
    gm = [r for r in A if isnum(r.get("h_lowg")) and isnum(r.get("pcum_gauge_model")) and r["pcum_gauge_model"] > 0]
    res["closure_low_gauge"] = rnd(sum(r["h_lowg"] for r in gm) * 1000 / sum(r["pcum_gauge_model"] for r in gm)) if gm else None
    cm = [r for r in A if isnum(r.get("pcum_cmorph_model")) and r["pcum_cmorph_model"] > 0]
    res["closure_rim_cmorph"] = rnd(sum(r["h_rim"] for r in cm) * 1000 / sum(r["pcum_cmorph_model"] for r in cm)) if cm else None
    if gg:
        bb, _ = boot_pairs([(r["station"], r["season"]) for r in gg],
                           {"C": ([r["h_obs"] * 1000 for r in gg], [r["pcum_gauge_obs"] for r in gg])})
        res["closure_obs_over_gauge_ci95"] = ci95(bb["C"])[0]
    # 描述：φ 只用 5–10 m 段（1–5 m 段在设计时已部分看过）
    for r in A:
        r["h_obs510"] = 2.5 * (-r["obs5"] / r["s0_5"]) + 2.5 * (-r["obs10"] / r["s0_10"])
        r["h_rim510"] = 2.5 * (-r["rim5"] / r["s0_5"]) + 2.5 * (-r["rim10"] / r["s0_10"])
    res["ratio_5_10m_segment"] = rnd(ratio_of(A, "h_obs510", "h_rim510"))
    # 描述：延伸到 20/25 m
    ext = []
    for r in A:
        for zd in (20.0, 25.0):
            if isnum(r.get(f"obs{int(zd)}")) and isnum(r.get(f"rim{int(zd)}")) and isnum(r.get(f"s0_{int(zd)}")):
                w = dict(W_10)
                w[10.0] += (zd - 10) / 2
                w[zd] = (zd - 10) / 2
                x = dict(r)
                x["h_o"], x["h_r"] = h(r, "obs", w), h(r, "rim", w)
                if isnum(x["h_o"]) and isnum(x["h_r"]):
                    ext.append(x)
                break
    res["deep_extension"] = {"n": len(ext), "Q_deep": rnd(ratio_of(ext, "h_o", "h_r")) if ext else None,
                             "R1_same_events": rnd(ratio_of(ext, "obs1", "rim1")) if ext else None}
    return res


def budget_05m_all(rows):
    """0–5 m 版收支（A＋B，描述；B 群 1 m 层为替代层时仍按其值代 0–1 m）。"""
    X = []
    for r in rows:
        if not r.get("model_ok"):
            continue
        vals = [r.get(k) for k in ("obs1", "obs5", "rim1", "rim5", "s0_1", "s0_5")]
        if all(isnum(v) for v in vals):
            X.append({"station": r["station"], "season": r["season"],
                      "ho": 3 * (-r["obs1"] / r["s0_1"]) + 2 * (-r["obs5"] / r["s0_5"]),
                      "hr": 3 * (-r["rim1"] / r["s0_1"]) + 2 * (-r["rim5"] / r["s0_5"]),
                      "obs1": r["obs1"], "rim1": r["rim1"]})
    return {"n": len(X), "Q5": rnd(ratio_of(X, "ho", "hr")) if X else None,
            "R1_same_events": rnd(ratio_of(X, "obs1", "rim1")) if X else None}


def spearman(xs, ys):
    from scipy import stats
    pts = [(x, y) for x, y in zip(xs, ys) if isnum(x) and isnum(y)]
    if len(pts) < 10:
        return {"n": len(pts), "rho": None}
    r = stats.spearmanr([p[0] for p in pts], [p[1] for p in pts])
    return {"n": len(pts), "rho": rnd(float(r.statistic), 4)}


def same_sensor_hourly(stations, events):
    np = _np()
    out = {}
    for stn in stations:
        a, b = stn["ser"]["sss05"], stn["ser"]["s1"]
        ok = np.isfinite(a) & np.isfinite(b)
        win = np.zeros(stn["n"], bool)
        for e in events:
            if e["_stn"] == stn["name"]:
                win[max(0, e["i"] - 6):e["i"] + 6] = True
        d_all = np.abs(a[ok] - b[ok])
        d_win = np.abs(a[ok & win] - b[ok & win])
        out[stn["name"]] = {"group": stn["group"], "n_hours": int(ok.sum()),
                            "frac_le_0.003_all": rnd(float(np.mean(d_all <= 0.003)), 4) if len(d_all) else None,
                            "n_hours_event_windows": int((ok & win).sum()),
                            "frac_le_0.003_event_windows": rnd(float(np.mean(d_win <= 0.003)), 4) if len(d_win) else None}
    return out


# ======================================================================== 判定（M17）
def _cell(block, k):
    return block["cells"][k] if block.get("evaluable") and block.get("cells") else None


def decide(S):
    """S：诊断汇总。返回 verdict dict。每支柱：pass／weak（主 CI 过、站簇 CI 不过）／fail／not_evaluable。"""
    pil = {}

    def upper_below(cell, thr):
        if cell is None or not cell.get("evaluable"):
            return "not_evaluable"
        hi, hs = (cell.get("ci95") or [None, None])[1], (cell.get("ci95_station") or [None, None])[1]
        if hi is None:
            return "not_evaluable"
        if hi < thr:
            return "pass" if (hs is not None and hs < thr) else "weak"
        return "fail"

    pil["P1a"] = upper_below(_cell(S["D_U"], 2), SESOI_LO)
    pil["P1b"] = upper_below(_cell(S["D_Z"], 0), SESOI_LO)
    pil["P4"] = upper_below(_cell(S["D_P"], 2), SESOI_LO)
    t = S["D_T_order"]
    if not t.get("evaluable") or not all(c["evaluable"] for c in (t["cells"][0], t["cells"][2])):
        pil["P2a"] = "not_evaluable"
    else:
        lo, los = t["contrast_top_minus_bottom_ci95_ss"][0], t["contrast_top_minus_bottom_ci95_st"][0]
        pil["P2a"] = "fail" if lo is None or lo <= 0 else ("pass" if los is not None and los > 0 else "weak")
    m = S["D_T_level"]
    if not m.get("evaluable") or m.get("ci95") is None or m["ci95"][0] is None:
        pil["P2b"] = "not_evaluable"
    else:
        lo, los = m["ci95"][0], (m.get("ci95_station") or [None])[0]
        pil["P2b"] = "fail" if lo < RTH_MIN else ("pass" if los is not None and los >= RTH_MIN else "weak")
    b = S["D_B"]
    if not b.get("evaluable") or b.get("phi") is None or b.get("phi_ci95_ss", [None])[0] is None:
        pil["P3"] = "not_evaluable"
    else:
        ok = b["phi_ci95_ss"][0] > 0
        jk = b.get("jackknife_phi_station") or {}
        oks = bool(jk.get("evaluable")) and jk["loo_range"][0] is not None and jk["loo_range"][0] > 0   # 站少（≈5），稳健性用逐站剔除
        pil["P3"] = "fail" if not ok else ("pass" if oks else "weak")
    # 否定信号
    sig = {}
    dz = S["D_Z"]
    c0 = _cell(dz, 0)
    sig["trap_S"] = bool(dz.get("evaluable") and c0 and c0.get("ci95") and c0["ci95"][1] is not None
                         and c0["ci95"][1] >= SESOI_LO and dz["contrast_top_minus_bottom_ci95_ss"][1] is not None
                         and dz["contrast_top_minus_bottom_ci95_ss"][1] < 0)
    du = S["D_U"]
    c2 = _cell(du, 2)
    sig["trap_S_wind"] = bool(du.get("evaluable") and c2 and c2.get("ci95") and c2["ci95"][1] is not None
                              and c2["ci95"][1] >= SESOI_LO and du["contrast_top_minus_bottom_ci95_ss"][0] is not None
                              and du["contrast_top_minus_bottom_ci95_ss"][0] > 0)
    sig["lateral_L"] = bool(b.get("evaluable") and b.get("phi_ci95_ss") and b["phi_ci95_ss"][1] is not None
                            and 0 <= b["phi_ci95_ss"][1] < PHI_LAT)
    sig["shallow_S_budget"] = bool(b.get("evaluable") and b.get("phi_ci95_ss") and b["phi_ci95_ss"][1] is not None
                                   and b["phi_ci95_ss"][1] < 0)
    sig["mismatch_M"] = pil["P4"] == "fail"
    vals = list(pil.values())
    if all(v == "pass" for v in vals):
        verdict = "positive"
    elif (sig["trap_S"] or sig["trap_S_wind"] or sig["shallow_S_budget"]) and (sig["lateral_L"] or sig["mismatch_M"]):
        verdict = "negative_multiple"
    elif sig["trap_S"] or sig["trap_S_wind"] or sig["shallow_S_budget"]:
        verdict = "negative_S"
    elif sig["lateral_L"] or sig["mismatch_M"]:
        verdict = "negative_L"
    else:
        verdict = "mixed"
    return {"pillars": pil, "negative_signals": sig, "verdict": verdict,
            "rule": "positive＝P1a、P1b、P2a、P2b、P3、P4 全部 pass（weak 或 not_evaluable 均不算）；否则按否定信号分出口，"
                    "无否定信号即 mixed。"}


# ======================================================================== 主流程
def build(args, out_dir, log):
    orig = p1.load_gtmba_file
    p1.load_gtmba_file = _wrap_loader(orig)
    try:
        stations, ev, val, prim, extra, plan = p2.build_all(args, out_dir, log)
    finally:
        p1.load_gtmba_file = orig
    attach_deep(stations)
    return stations, ev[p2.PRIMARY_THRESHOLD], val, prim, plan


def load_pixels(p2_dir, hours):
    data, missing = p2.PixelCache.parse(os.path.join(p2_dir, "cmorph_pixels.csv"))
    need = set(hours)
    return data, {"hours_needed": len(need), "hours_present": sum(1 for h in need if h in data),
                  "hours_missing_marked": sum(1 for h in need if h in missing)}


def availability(stations, events):
    """--plan：只数有效小时，不算任何变化量。"""
    np = _np()
    sbn = {s["name"]: s for s in stations}
    cnt = {}
    for e in events:
        stn = sbn[e["_stn"]]
        i = e["i"]
        c = cnt.setdefault(stn["name"], {"events": 0, "with_s5": 0, "with_s10": 0, "with_s5_s10": 0, "with_s20_or_25": 0})
        c["events"] += 1

        def has(k):
            a = stn["ser"][k]
            return bool(np.isfinite(a[max(0, i - 6):i]).any() and np.isfinite(a[i:i + 6]).any())
        h5, h10 = has("s5"), has("s10")
        c["with_s5"] += h5
        c["with_s10"] += h10
        c["with_s5_s10"] += h5 and h10
        c["with_s20_or_25"] += has("s20") or has("s25")
    tot = {k: sum(v[k] for v in cnt.values()) for k in ("events", "with_s5", "with_s10", "with_s5_s10", "with_s20_or_25")}
    return {"by_station": cnt, "total": tot}


def analyze(rows, stations, events, ref, out_dir=None):
    """逐事件派生量 → 各诊断汇总＋判定。ref＝P2 p2_summary.primary.R_RIM_1m（V0；自测传 None 跳过）。"""
    S = {}
    # V0：复算主结果
    main = pooled([r for r in rows if r.get("model_ok")], "obs1", "rim1", "V0 R_RIM(1 m)")
    if ref is None:
        v0_ok = True
    else:
        v0_ok = (abs(main["point"] - ref["point"]) <= V0_TOL and all(abs(a - b) <= V0_TOL for a, b in zip(main["ci95"], ref["ci95"]))
                 and main["n"] == ref["n_events"])
    S["V0"] = {"this": {k: main[k] for k in ("n", "point", "ci95", "G_station_season")},
               "p2": {k: ref[k] for k in ("n_events", "point", "ci95")} if ref else None, "pass": bool(v0_ok)}
    if not v0_ok:
        if out_dir:
            p2.jdump(S, os.path.join(out_dir, "p4_mech_abort.json"))
        raise RuntimeError(f"V0 复算不一致：{S['V0']}")
    ok = [r for r in rows if r.get("model_ok")]
    # D-U
    S["D_U"] = tercile_block(ok, lambda r: r["u10"], "R_RIM(1 m) by U10 tercile", theory=("low1", "rim1"), loglog=True)
    S["D_U_k27"] = tercile_block(ok, lambda r: r["u10"], "K27 变体 by U10 tercile", den="rim1_k27")
    cross = {}
    lab_u, _, _ = terciles([r for r in ok if isnum(r["u10"])], lambda r: r["u10"])
    lab_r, _, _ = terciles([r for r in ok if isnum(r["u10"])], lambda r: r["rain24"])
    su = [r for r in ok if isnum(r["u10"])]
    for a in range(3):
        for b_ in range(3):
            sub = [r for j, r in enumerate(su) if lab_u[j] == a and lab_r[j] == b_]
            cross[f"U{a}_rain{b_}"] = {k: v for k, v in pooled(sub, "obs1", "rim1", "").items() if not k.startswith("_")} \
                if len(sub) >= 10 else {"n": len(sub)}
    S["D_U_cross_rain"] = cross
    # D-Z
    S["D_Z"] = tercile_block(ok, lambda r: r["zeta"], "R_RIM(1 m) by zeta=z1/L tercile (k=0 最不稳定)")
    # D-T
    S["D_T_order"] = tercile_block(ok, lambda r: r["rho"], "R_RIM(1 m) by LoW 预测比 rho tercile", theory=("low1", "rim1"))
    lv = pooled(ok, "obs1", "low1", "R_obs/R_th＝Σobs1/ΣLoW1")
    S["D_T_level"] = {k: v for k, v in lv.items() if not k.startswith("_")}
    S["D_T_level"]["R_th_pooled"] = rnd(ratio_of([r for r in ok if isnum(r["low1"]) and isnum(r["rim1"])], "low1", "rim1"))
    t0v = [r for r in ok if all(isnum(r.get(k)) for k in ("low0", "rim0", "low1", "rim1"))]
    S["D_T_surface"] = {"n": len(t0v), "s_th_surface": rnd(ratio_of(t0v, "low0", "rim0")),
                        "s_th_1m": rnd(ratio_of(t0v, "low1", "rim1")),
                        "note": "壁层解隐含的海面比／1 m 比；≈1 表示在该理论下 1 m 比值可外推到海面（描述）"}
    # D-B
    S["D_B"] = budget_block(rows)
    S["D_B_05m_all"] = budget_05m_all(rows)
    # D-P
    S["D_P"] = tercile_block(ok, lambda r: r["lam"], "R_RIM(1 m) by λ=ln(G6/C6) tercile（k=2 雨量计相对最多）")
    m3 = pooled([r for r in ok if isnum(r.get("rim1_m3"))], "obs1", "rim1_m3", "3×3 均值驱动 R_RIM(1 m)")
    S["D_P_3x3"] = {k: v for k, v in m3.items() if not k.startswith("_")}
    # D-S
    ds = {}
    for g in (None, p2.GROUP_A, p2.GROUP_B):
        sub = [r for r in ok if g is None or r["group"] == g]
        tag = "all" if g is None else ("A" if g == p2.GROUP_A else "B")
        ds[f"inversion_matched_{tag}"] = {k: v for k, v in pooled(sub, "obs05_m", "obs1_m", "Σ0.5/Σ1 同小时").items()
                                          if not k.startswith("_")}
        ds[f"inversion_unmatched_{tag}"] = rnd(ratio_of([r for r in sub if isnum(r["obs05"]) and isnum(r["obs1"])],
                                                        "obs05", "obs1"))
    ds["R1_obs_matched_model_full"] = {k: v for k, v in pooled(ok, "obs1_m", "rim1", "").items() if not k.startswith("_")}
    ds["R1_obs_matched_model_matched"] = {k: v for k, v in pooled(ok, "obs1_m", "rim1_m", "").items() if not k.startswith("_")}
    for k in (1, 2):
        ds[f"R1_model_shift_{30 * k}min"] = {kk: v for kk, v in pooled(ok, "obs1", f"rim1_shift{k}", "").items()
                                             if not kk.startswith("_")}
        ds[f"shift{k}_pad_zero_events"] = sum(int(r.get(f"shift{k}_pad_zero", 0) > 0) for r in ok)
    ds["same_sensor_hourly"] = same_sensor_hourly(stations, events)
    S["D_S"] = ds
    # 非盲程度量化
    S["nonblind_quantification"] = {
        "spearman_u10_vs_E5proxy": spearman([r["u10"] for r in ok], [r["e5proxy"] for r in ok]),
        "spearman_rho_vs_rain24": spearman([r["rho"] for r in ok], [r["rain24"] for r in ok]),
        "spearman_rho_vs_E5proxy": spearman([r["rho"] for r in ok], [r["e5proxy"] for r in ok]),
        "spearman_zeta_vs_rain24": spearman([r["zeta"] for r in ok], [r["rain24"] for r in ok]),
        "spearman_lam_vs_rain24": spearman([r["lam"] for r in ok], [r["rain24"] for r in ok])}
    for key in ("D_U", "D_U_k27", "D_Z", "D_T_order", "D_P"):
        S[key] = {k: v for k, v in S[key].items() if not k.startswith("_")}
    S["verdict"] = decide(S)
    return S


def run(args, out_dir, log, plan_only):
    t0 = time.monotonic()
    stations, events, val, prim, plan = build(args, out_dir, log)
    data, cov = load_pixels(args.p2_dir, prim)
    info = {"script": "p4_mech.py", "version": VERSION,
            "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "n_events": len(events), "rebuild_validation": {k: val[k] for k in ("n_ref", "n_rebuilt", "d4_count")},
            "cmorph_pixel_cache": cov, "deep_salinity_stations": sorted(_DEEP),
            "downloads": {"cache_http_requests": plan.get("cache_http_requests"), "cache_links": plan.get("cache_links"),
                          "new_data_bytes_expected": 0}}
    if plan_only:
        info["availability_counts_only"] = availability(stations, events)
        info["runtime_s"] = round(time.monotonic() - t0, 1)
        p2.jdump(info, os.path.join(out_dir, "p4_mech_plan.json"))
        print(json.dumps(info, ensure_ascii=False)[:3000], flush=True)
        return 0
    sbn = {s["name"]: s for s in stations}
    pixc, pixm = PixView(data, "c"), PixView(data, "m")
    rows = [event_mech(sbn[e["_stn"]], e, pixc, pixm)[0] for e in events]
    ref = json.load(open(os.path.join(args.p2_dir, "p2_summary.json"), encoding="utf-8"))["primary"]["R_RIM_1m"]
    S = dict(info)
    S.update(analyze(rows, stations, events, ref, out_dir))
    S["runtime_s"] = round(time.monotonic() - t0, 1)
    p2.jdump(S, os.path.join(out_dir, "p4_mech_summary.json"))
    cols = sorted({k for r in rows for k in r})
    with open(os.path.join(out_dir, "p4_mech_events.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (rnd(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
    print(json.dumps(S["verdict"], ensure_ascii=False), flush=True)
    return 0


# ======================================================================== 自测（合成数据，无网络）
def _fd_wall(a_of_t, F, t_end, dt=5.0, H=200.0, dz=0.05):
    """∂c/∂t=∂z(a(t)·z·∂z c) 的隐式有限差分（零通量上下界），初始 F 置于首格；返回 (z 格心, c)。"""
    np = _np()
    from scipy.linalg import solve_banded
    n = int(H / dz)
    zc = (np.arange(n) + 0.5) * dz
    zf = np.arange(1, n) * dz
    c = np.zeros(n)
    c[0] = F / dz
    t = 0.0
    while t < t_end - 1e-9:
        a = a_of_t(t + dt / 2)
        k = a * zf / dz ** 2 * dt
        ab = np.zeros((3, n))
        ab[1] = 1.0
        ab[1, :-1] += k
        ab[1, 1:] += k
        ab[0, 1:] = -k
        ab[2, :-1] = -k
        c = solve_banded((1, 1), ab, c)
        t += dt
    return zc, c


def selftest(out_dir=None):
    np = _np()
    res, fails = [], []

    def check(cond, msg):
        res.append(("PASS" if cond else "FAIL") + " " + msg)
        if not cond:
            fails.append(msg)

    # 1 壁层解：解析式对有限差分（常 u*、时变 u*）、质量守恒
    us = 7e-3
    a = KAPPA * us
    F = 0.01
    zc, c = _fd_wall(lambda t: a, F, 3600.0)
    lam = a * 3600.0
    for z in (1.0, 5.0, 10.0):
        ana = F / lam * math.exp(-z / lam)
        num = float(np.interp(z, zc, c))
        check(abs(num / ana - 1) < 0.03, f"壁层解析 vs FD（常 u*）z={z}: {num:.4e}/{ana:.4e}")
    check(abs(float(c.sum() * 0.05) / F - 1) < 1e-6, "FD 质量守恒")
    zc2, c2 = _fd_wall(lambda t: KAPPA * (3e-3 if t < 1800 else 9e-3), F, 3600.0)
    lam2 = KAPPA * (3e-3 * 1800 + 9e-3 * 1800)
    for z in (1.0, 5.0):
        ana = F / lam2 * math.exp(-z / lam2)
        check(abs(float(np.interp(z, zc2, c2)) / ana - 1) < 0.03, f"壁层解析 vs FD（时变 u*）z={z}")
    # 2 low_conc 离散实现：单半步入水、常风 → 与解析式一致
    L = 84
    P = np.zeros(L)
    P[60] = 10.0                                      # 10 mm/h × 0.5 h = 5 mm
    U = np.full(L, 6.0)
    usw = float(ustar_w(np.array([6.0]))[0])
    cz = low_conc(P, U, (1.0,))[1.0]
    n = 65
    lam_n = KAPPA * usw * HALF_S * (n - 60 + 0.5)
    check(abs(cz[n] / (0.005 / lam_n * math.exp(-1 / lam_n)) - 1) < 1e-9, "low_conc 与解析式逐值一致")
    check(np.all(cz[:60] == 0), "low_conc 因果性（入水前为 0）")
    # 3 Cd、u*、Obukhov 手算
    check(abs(float(drag_coef(6.0)) - 1e-3 * (0.45 + 0.142 + 0.4584)) < 1e-12, "Cd(6 m/s) 手算")
    check(abs(float(drag_coef(0.1)) - float(drag_coef(0.5))) < 1e-15, "Cd 低风下限 0.5 m/s")
    z = obukhov_zeta(1.0, 7e-3, 35.0, 2.0)
    Bf = 9.81 * 7.5e-4 * 35.0 * 2.0 / 3.6e6
    check(abs(z - KAPPA * Bf / 7e-3 ** 3) < 1e-12 and obukhov_zeta(1.0, 7e-3, 35.0, 0.0) == 0.0, "Obukhov ζ 手算、无雨为 0")
    # 4 雨量时间匹配
    check(abs(cum_mid_hourly([1, 1, 1, 1, 1, 1]) - 3.0) < 1e-12, "箱中点累计（均匀 1 mm/h → 3 mm）")
    check(abs(cum_model_hourly([0.5] * 12) - (sum((k + 0.75) for k in range(6)) / 6)) < 1e-12, "模型时刻累计")
    # 5 RIM 时间后移：shift=0 等于 p2.model_hourly；shift=1 等于手工前移
    rng = np.random.default_rng(1)
    Pm = np.where(rng.random(84) < 0.2, rng.gamma(2, 3, 84), 0.0)
    Um = rng.uniform(2, 9, 84)
    h0v, _ = p2.model_hourly(Pm, Um, 35.0, {0.5: 0.5, 1.0: 1.0, 5.0: 5.0})
    check(np.allclose(rim_hourly(Pm, Um, 1.0), h0v[1.0], equal_nan=True), "rim_hourly(shift 0)＝p2.model_hourly")
    Pp = np.concatenate([[0.0], Pm])
    Up = np.concatenate([[Um[0]], Um])
    Ffull = p2.rim_factor(Pp, Up, 1.0)                 # 半步 47..83（原下标）
    man = np.array([(Ffull[2 * q] + Ffull[2 * q + 1]) / 2 for q in range(18)])
    check(np.allclose(rim_hourly(Pp, Up, 1.0, shift=1), man, equal_nan=True), "rim_hourly(shift 1) 手工一致")
    # 6 收支权重
    check(abs(sum(W_10.values()) - 10) < 1e-12 and abs(sum(W_5.values()) - 5) < 1e-12, "梯形权重和＝积分深度")
    # 7 统计：三分位对比、φ、判定
    rng = np.random.default_rng(7)
    rows = []
    for st in range(12):
        for se in range(4):
            for _ in range(15):
                u = rng.uniform(1, 11)
                den = -rng.gamma(2, 0.05)
                true_r = 0.2 + 0.05 * u                   # R 随 U 上升（0.25→0.75）
                rows.append({"station": f"S{st}", "season": se, "u10": u, "rim1": den,
                             "obs1": true_r * den + rng.normal(0, 0.01), "low1": 0.5 * den})
    tb = tercile_block(rows, lambda r: r["u10"], "syn", theory=("low1", "rim1"), loglog=True)
    c = tb["cells"]
    check(tb["evaluable"] and c[0]["point"] < c[1]["point"] < c[2]["point"], "三分位单调（合成 R 随 U 上升）")
    check(tb["contrast_top_minus_bottom_ci95_ss"][0] > 0, "上－下对比 CI>0（合成）")
    check(abs(c[2]["theory_point"] - 0.5) < 1e-9, "理论比值按格合并")
    check(tb["loglog_slope"]["point"] > 0 and abs(tb["loglog_slope_theory"]) < 1e-9, "对数斜率（观测>0、理论 0）")
    for r in rows:
        r["obs1n"] = 0.3 * r["rim1"] + rng.normal(0, 0.01)
    tn = tercile_block(rows, lambda r: r["u10"], "null", num="obs1n")
    lo, hi = tn["contrast_top_minus_bottom_ci95_ss"]
    check(lo < 0 < hi, "零效应对比 CI 含 0（合成）")
    # φ：观测剖面更深
    brow = []
    for st in range(6):
        for se in range(4):
            for _ in range(5):
                m1, m5, m10 = -0.2, -0.06, -0.005
                o = rng.normal(0, 0.002, 3)
                brow.append({"station": f"S{st}", "season": se, "group": p2.GROUP_A, "model_ok": True,
                             "obs1": 0.3 * m1 + o[0], "obs5": 0.8 * 0.3 * m1 + o[1], "obs10": 0.6 * 0.3 * m1 + o[2],
                             "rim1": m1, "rim5": m5, "rim10": m10, "lowg1": NAN, "lowg5": NAN, "lowg10": NAN,
                             "s0_1": 35.0, "s0_5": 35.0, "s0_10": 35.0, "pcum_gauge_obs": 5.0,
                             "pcum_gauge_model": 5.0, "pcum_cmorph_model": 5.0})
    bb = budget_block(brow)
    hq = (3 * 0.06 + 4.5 * 0.048 + 2.5 * 0.036) / (3 * 0.2 + 4.5 * 0.06 + 2.5 * 0.005)
    phi_true = (hq - 0.3) / 0.7
    check(bb["evaluable"] and abs(bb["phi"] - phi_true) < 0.02 and bb["phi_ci95_ss"][0] > 0,
          f"φ 合成（真值 {phi_true:.3f}，得 {bb['phi']}）")
    # 判定逻辑：构造各分支
    def cell(hi, his=None, ev=True, pt=0.3):
        return {"evaluable": ev, "ci95": [0.1, hi], "ci95_station": [0.1, his if his is not None else hi], "point": pt}

    def blk(cells, con=(0.05, 0.3), cons=None):
        return {"evaluable": True, "cells": cells, "contrast_top_minus_bottom_ci95_ss": list(con),
                "contrast_top_minus_bottom_ci95_st": list(cons or con)}
    good = {"D_U": blk([cell(0.4), cell(0.4), cell(0.5)]),
            "D_Z": blk([cell(0.5), cell(0.4), cell(0.4)], con=(-0.1, 0.1)),
            "D_P": blk([cell(0.4), cell(0.4), cell(0.5)]),
            "D_T_order": blk([cell(0.4), cell(0.4), cell(0.5)], con=(0.05, 0.3)),
            "D_T_level": {"evaluable": True, "ci95": [0.8, 1.6], "ci95_station": [0.7, 1.8]},
            "D_B": {"evaluable": True, "phi": 0.4, "phi_ci95_ss": [0.1, 0.7], "phi_ci95_st": [0.05, 0.8],
                    "jackknife_phi_station": {"evaluable": True, "loo_range": [0.2, 0.5]}}}
    check(decide(good)["verdict"] == "positive", "判定：全过 → positive")
    import copy
    g2 = copy.deepcopy(good)
    g2["D_U"]["cells"][2] = cell(0.6, his=0.7)
    check(decide(g2)["verdict"] == "mixed" and decide(g2)["pillars"]["P1a"] == "weak", "判定：站簇不过 → weak → mixed")
    g3 = copy.deepcopy(good)
    g3["D_Z"]["cells"][0] = cell(0.9)
    g3["D_Z"]["contrast_top_minus_bottom_ci95_ss"] = [-0.5, -0.1]
    check(decide(g3)["verdict"] == "negative_S", "判定：稳定度困住信号 → negative_S")
    g4 = copy.deepcopy(good)
    g4["D_B"] = {"evaluable": True, "phi": 0.01, "phi_ci95_ss": [-0.05, 0.08], "phi_ci95_st": [-0.1, 0.1]}
    check(decide(g4)["verdict"] == "negative_L", "判定：深层无回收 → negative_L")
    g4b = copy.deepcopy(good)
    g4b["D_B"] = {"evaluable": True, "phi": -0.2, "phi_ci95_ss": [-0.35, -0.05]}
    check(decide(g4b)["verdict"] == "negative_S", "判定：观测淡水比模型更浅 → negative_S")
    g4c = copy.deepcopy(good)
    g4c["D_B"] = {"evaluable": True, "phi": 0.05, "phi_ci95_ss": [-0.1, 0.2]}
    check(decide(g4c)["verdict"] == "mixed", "判定：φ 区间跨 0 且上端 ≥0.10 → mixed")
    g5 = copy.deepcopy(good)
    g5["D_B"] = {"evaluable": False, "phi": None}
    check(decide(g5)["verdict"] == "mixed" and decide(g5)["pillars"]["P3"] == "not_evaluable", "判定：P3 不可评 → mixed")
    g7 = copy.deepcopy(good)
    g7["D_B"]["jackknife_phi_station"]["loo_range"] = [-0.05, 0.5]
    check(decide(g7)["pillars"]["P3"] == "weak", "判定：P3 逐站剔除有负值 → weak")
    g6 = copy.deepcopy(good)
    g6["D_T_level"] = {"evaluable": True, "ci95": [0.2, 0.4], "ci95_station": [0.2, 0.45]}
    check(decide(g6)["pillars"]["P2b"] == "fail" and decide(g6)["verdict"] == "mixed", "判定：量级不符 → P2b fail → mixed")
    # 8 深层盐度解析：合成 DAP ascii
    import gzip
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        pth = os.path.join(td, "x.ascii.gz")
        t = [20000.0 + k / 24 for k in range(4)]
        with gzip.open(pth, "wt") as f:
            f.write("Dataset {\n} x;\n---------------------------------------------\n")
            f.write("TIME[4]\n" + ", ".join(str(x) for x in t) + "\n\nDEPTH[3]\n1.0, 5.0, 10.0\n\n")
            f.write("PSAL.PSAL[4][3]\n[0], 35.0, 35.1, 35.2\n[1], 34.9, 35.1, 35.21\n[2], 34.8, 35.0, 35.22\n"
                    "[3], 34.7, 35.0, 35.23\n\n")
            f.write("PSAL_QC.PSAL_QC[4][3]\n[0], 1, 1, 1\n[1], 1, 1, 4\n[2], 1, 1, 2\n[3], 1, 1, 1\n\n")
        _DEEP.clear()
        add_deep_from_ascii("SYN", pth, 1)
        got = _DEEP.get("SYN", {}).get("s10", {})
        vals = sorted((h, s / c_) for h, (s, c_) in got.items())
        check(len(vals) == 3 and abs(vals[0][1] - 35.2) < 1e-9 and abs(vals[-1][1] - 35.23) < 1e-9,
              "深层盐度解析：10 m 取值、QC=4 剔除")
        _DEEP.clear()
    # 9 配对采样：0.5 m 与 1 m 同源时倒置比＝1
    x = {"station": "S", "season": 0}
    rr = [dict(x, station=f"S{k % 12}", season=k % 4, obs05_m=v, obs1_m=v) for k, v in enumerate(rng.normal(-0.05, 0.02, 200))]
    pr = pooled(rr, "obs05_m", "obs1_m", "syn")
    check(abs(pr["point"] - 1) < 1e-12 and pr["ci95"][0] == 1.0, "配对采样同源 → 倒置比 1")
    # 10 端到端：p2 合成两站夹具（A 真 1 m、B 替代层）＋合成 10 m → event_mech → analyze（真值比 0.5）
    with tempfile.TemporaryDirectory() as td:
        fx = p2._e2e_fixture(td)
        pixd = {h: {nm: (v[0], v[1], v[0], v[1], 9, 9) for nm, v in d.items()} for h, d in fx["cm"].items()}
        pixc, pixm = PixView(pixd, "c"), PixView(pixd, "m")
        for stn in fx["stations"]:
            n_ = stn["n"]
            Pfull = np.zeros(2 * n_)
            for hh in range(n_):
                Pfull[2 * hh:2 * hh + 2] = pixd[stn["h0"] + hh][stn["name"]][:2]
            Uh = np.repeat(stn["ser"]["wind"] * p2.WIND_FACTOR, 2)
            for key, zt in (("s10", 10.0), ("s20", 20.0), ("s25", 25.0)):
                if stn["group"] == p2.GROUP_A and key != "s25":
                    F = np.ones(2 * n_)
                    F[48:] = p2.rim_factor(Pfull, Uh, zt)
                    stn["ser"][key] = 35.0 * (1 - 0.5 * (1 - F.reshape(-1, 2).mean(axis=1)))
                else:
                    stn["ser"][key] = np.full(n_, np.nan)
        sbn = {s["name"]: s for s in fx["stations"]}
        evs = fx["ev"][p2.PRIMARY_THRESHOLD]
        rows_e = [event_mech(sbn[e["_stn"]], e, pixc, pixm)[0] for e in evs]
        try:
            Se = analyze(rows_e, fx["stations"], evs, None)
            err = None
        except Exception as ex:
            Se, err = None, repr(ex)
        check(err is None, f"端到端 analyze 无异常（{err}）")
        if Se:
            check(abs(Se["V0"]["this"]["point"] - 0.5) < 0.03, f"端到端 R_RIM(1 m)≈0.5（得 {Se['V0']['this']['point']}）")
            check(Se["D_B"]["n"] > 0 and Se["D_B"]["phi"] is not None and abs(Se["D_B"]["phi"]) < 0.1,
                  f"端到端 φ≈0（各深度同比，得 {Se['D_B'].get('phi')}）")
            check(abs(Se["D_T_level"]["R_th_pooled"] or 0) > 0 and Se["verdict"]["verdict"] in
                  ("positive", "mixed", "negative_S", "negative_L", "negative_multiple"), "端到端 判定字段齐")
            check(abs((Se["D_P_3x3"]["point"] or 0) - Se["V0"]["this"]["point"]) < 1e-9, "端到端 3×3＝中心（夹具同值）")
    out = {"script": "p4_mech.py", "version": VERSION, "n_checks": len(res), "n_fail": len(fails), "checks": res,
           "pass": not fails}
    if out_dir:
        p2.jdump(out, os.path.join(out_dir, "p4_mech_selftest.json"))
    for line in res:
        print(line)
    print(f"自测 {len(res) - len(fails)}/{len(res)} 通过", flush=True)
    return 0 if not fails else 4


def main(argv=None):
    ap = argparse.ArgumentParser(description="P4 机制诊断")
    ap.add_argument("--out")
    ap.add_argument("--p1-events", default=p2.P1_EVENTS_DEFAULT)
    ap.add_argument("--p1b-dir", default=p2.P1B_DIR_DEFAULT)
    ap.add_argument("--p2-dir", default=P2_DIR_DEFAULT)
    ap.add_argument("--plan", action="store_true")
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
    if args.plan:
        rc = selftest(out_dir)
        if rc:
            return rc
    log = p1.Log(os.path.join(out_dir, "p4_mech_log.txt"))
    log.log(f"=== start {VERSION} plan={args.plan} p2_dir={args.p2_dir} out={out_dir}", echo=True)
    try:
        rc = run(args, out_dir, log, args.plan)
    except Exception:
        log.log("FATAL：\n" + traceback.format_exc(), echo=True)
        return 3
    log.log("=== done")
    log.close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
