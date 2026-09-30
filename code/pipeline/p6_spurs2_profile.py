#!/usr/bin/env python3
"""p6_spurs2_profile.py — P6：SPURS-2 船载／浮标近表层剖面检验。用不依赖 L 波段反演的点观测，
比较雨后「海面（盐度 snake，0.01–0.02 m）」与「1 m」的观测淡化和 RIM-3 同层预测之比，判别
R1（RIM-3 幅度整体偏大，海面也偏大）与 R2（海面幅度对、剖面画得太深）；另检验 P5「雨中过境超出」疑点。

判读在 judge()／judge_rain()。性质：探索性；方法与判读在读剖面盐度数值之前写定。

用法：
  python p6_spurs2_profile.py --selftest [--out DIR]     合成数据自测（无网络、不读原件）
  python p6_spurs2_profile.py --full --root RAW --cache CACHE --out OUT
        自测 → 读原件（只读）→ 小时化（缓存 npz）→ 事件 → 模型 → 先写 p6_power.json → 再写 p6_summary.json、p6_events.csv
依赖：numpy、scipy、netCDF4＋同目录 p2_rim_test.py（RIM-3、风插补、簇 bootstrap；只 import 不改，
  其 import 的 p1_events.py、p1b_extend.py 也须在同目录）。
数据（raw 原始档，只读，不改不删不移动）：`<root>/SPURS2_{METEO,SALINITYSNAKE,SSP,MOORING_CENTRAL}/…nc`。
产物（<out>/）：p6_selftest.json、p6_power.json（先写）、p6_summary.json、p6_events.csv、p6_log.txt。缓存（<cache>/）：hourly_*.npz。
退出码：0 跑完（无论判读）；3 其他异常；4 自测不过。

实现选择（W 条）：
  W1 时间：各文件 time「days since 1950-01-01」→ Unix 秒；UTC 小时箱 [HH:00, HH+1:00)（P1 I4）。
  W2 小时化（盐度）：小时内有效样本均值，需 ≥25% 的期望样本数（期望＝3600/该变量中位采样间隔）；否则该小时缺。
  W3 雨：船（Revelle）用 METEO `precipitation_rate`（光学雨量计，mm/h，1 min）；浮标用 WHOI_Met `rainfall_rate`（mm/h，1 min）。
     事件用小时雨（≥50% 分钟有效，P1 I4）；RIM 用 30 min 半步均值（≥50% 分钟有效，否则置 0 并计数，同 P2 K21 雨量计口径）；负值置 0。
  W4 事件（主）：雨小时＝小时雨 ≥0.4 mm/h（P1 I5）；onset＝有效雨小时，其前 6 h 全有效且无雨；[t0, t0+24h) 全有效、雨小时累积 ≥10 mm
     （P1 I6 的 10 mm，干期由 24 h 放宽到 6 h，以保留足够事件）；另要求 [t0−30 h, t0+12 h) 在雨量记录内（RIM 24 h 雨史）。
  W5 风：船用 METEO `neutral_wind_speed_relative_to_earth`（已折算 10 m 中性）；浮标用 √(ewind²+nwind²)，按 P2 K9 同式由
     3.23 m（文件属性 wind_sensor_height）换算 10 m（×ln(10/z0)/ln(3.23/z0)，z0=1e-4）。小时均值，p2.fill_wind ≤6 h 线性插补，
     下限 0.1 m/s，重复到两个半步；[t0−30, t0+12) 内仍缺 → 事件模型无效（计数）。
  W6 RIM-3：p2.rim_factor（K2–K4：d0 7×7 双线性、48 历史半步＋当前半步、当前项 t=1 s）逐层算，84 个半步 [t0−30 h, t0+12 h)
     → 输出半步 48..83；模型小时值＝该小时两个半步输出的平均（P2 K5），对应小时 t0−6…t0+11。自测逐值对 p2.model_hourly。
  W7 层与模型深度：海面＝snake `salinity`（Salinity_Snake 文件，1 s；0.01–0.02 m）→ 模型 z=0（与 Witte 全球、P5 同）；
     1 m＝SSP `sea_water_salinity_at_1p10m` → z=1.10；剖面描述层 SSP 0.05／0.12／0.23／0.54 m、METEO `salinity_at_2m`／`_3m`（USPS）、
     `salinity_at_5m`（TSG）→ 各自名义深度；浮标 1 m＝CentralMooring_SubsurfaceTemperatureSalinity `salinity` 深度 1.0 m 那一列 → z=1.0。
  W8 ΔS（P2 I8／K6 同式，但模型与观测按小时配对）：ΔS_obs＝mean(观测, 雨后窗内有值小时) − median(观测, [t0−6, t0) 有值小时)；
     ΔS_RIM＝S0·[mean F(同一批雨后小时) − median F(同一批雨前小时)]，S0＝该层观测雨前中位数。雨前、雨后窗各至少 1 个有值小时，否则该层无效。
     主窗 [t0, t0+6 h)。（P2 的 1 m 观测几乎逐时齐全，配对与否等价；船载层断续，配对才对等。）
  W9 主估计量 R_sfc＝Σ ΔS_obs(snake)/Σ ΔS_RIM(0)，Revelle 两航次主事件集；R_1m,SSP、R_1m,moor 同式（一致性检查）。
  W10 簇：船＝航次×ISO 周（主 CI，p2.cluster_boot，B=10000，seed 20260926，百分位 95%）；另报事件级 bootstrap 与逐航次留一。
      浮标＝日历月簇。
  W11 可评：n≥15、簇≥4、两航次各 ≥3 事件、bootstrap 中 Σ ΔS_RIM ≥0 的比例 ≤1%。
  W12 判读 judge()：按序（不可评 → 支持 R1 → 支持 R2 → 不能区分），稳健标注＝事件级 CI 与逐航次留一同侧；
      平流标注＝R_rel5（snake−5 m TSG 差分）类别是否与主判相同。
  W13 雨中／雨停：[t0, t0+12 h) 内 snake 有值的小时按该小时 METEO 雨 ≥0.4 mm/h 分「雨中」「雨停」，各自按 W8 同式
      （只取该组小时）算 ΔS_obs、ΔS_RIM(0)，R_rain、R_dry 与 ln(R_rain/R_dry)，同一簇抽样；judge_rain() 判读。
  W14 R_rel5：逐小时 d_obs＝S_snake−S_5m（两者同小时都有值），ΔD_obs 按 W8；d_RIM＝S0_sfc·F(0)−S0_5·F(5)，同小时配对。
  W15 层比与 Q（SSP 配对子集＝snake 与 SSP 1.10 m 都有效的事件）：观测层比 ΣΔS_sfc/ΣΔS_1m、RIM 层比、Q＝R_sfc/R_1m（同一抽样）；
      SSP 内部层比（0.12 m 对 1.10 m）同法。
  W16 敏感性（只报，不改判读）：S1 干期 24 h（P1 原样）；S2 干期 12 h；S3 阈值 5 mm；S4 窗 [0,3) h；S5 窗 [0,12) h；
      S6 航迹位移 ≤30 km（雨前 6 h 与雨后 6 h 船位中位点距离，METEO 经纬度）；S7 海面改用 METEO `sea_surface_salinity`（snake 1 min 副本）；
      S8 海面模型深度 z=0.015 m；S9 当前项深度因子 t=1800 s（P2 K27 同式）；S10 雨前、雨后窗各需 ≥4 个有值小时。
  W17 功效（p6_power.json，先于任何观测 ΔS 合计写盘）：主事件集模型侧 x 的 n、Σx、均值、RMS、簇数；无雨对照窗（Revelle，
      [c−12, c+6) 无雨小时、c 每 6 h 取一个）上 snake 的 ΔS 标准差 σ_ctrl（观测噪声代理）与由此推的 SE(R)≈σ√n/|Σx|；
      另给 σ＝0.1／0.2／0.4 psu 三档假设的 SE。只作信息，不改判读。
Change Log：
  2026-09-27 初版（P6）。
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

import numpy as np

import p2_rim_test as p2

VERSION = "p6-2026-09-27a"
NAN = float("nan")
THRESH = 0.65
DEN_NONNEG_MAX = 0.01
RAIN_HOUR_MM = 0.4
PRE_H = 6
HIST_H = 30               # [t0−30 h, t0+12 h) 雨史
POST_MAX_H = 12
MIN_FRAC_SAL = 0.25
MIN_FRAC_RAIN = 0.5
WIND_Z_MOOR = 3.23
WIND_FACTOR_MOOR = math.log(10.0 / p2.WIND_Z0_M) / math.log(WIND_Z_MOOR / p2.WIND_Z0_M)
B = p2.BOOT_B
SEED = p2.BOOT_SEED
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
ROOT_DEFAULT = _rp.raw_sub("spurs2")  # [repro] 路径来自集中配置
CACHE_DEFAULT = _rp.path("p6_cache")  # [repro] 路径来自集中配置

CRUISES = {
    "RR1610": {"met": "SPURS2_METEO/SPURS2_RR1610_WHOI_Met.nc",
               "snake": "SPURS2_SALINITYSNAKE/SPURS2_RR1610_Salinity_Snake.nc",
               "ssp": "SPURS2_SSP/SPURS2_RR1610_SSP.nc"},
    "RR1720": {"met": "SPURS2_METEO/SPURS2_RR1720_WHOI_Met.nc",
               "snake": "SPURS2_SALINITYSNAKE/SPURS2_RR1720_Salinity_Snake.nc",
               "ssp": "SPURS2_SSP/SPURS2_RR1720_SSP.nc"},
}
# 层名 → (文件键, 变量名, 模型深度 m)
LAYERS = {
    "snake": ("snake", "salinity", 0.0),
    "snake_met": ("met", "sea_surface_salinity", 0.0),
    "ssp005": ("ssp", "sea_water_salinity_at_0p05m", 0.05),
    "ssp012": ("ssp", "sea_water_salinity_at_0p12m", 0.12),
    "ssp023": ("ssp", "sea_water_salinity_at_0p23m", 0.23),
    "ssp054": ("ssp", "sea_water_salinity_at_0p54m", 0.54),
    "ssp110": ("ssp", "sea_water_salinity_at_1p10m", 1.10),
    "usps2": ("met", "salinity_at_2m", 2.0),
    "usps3": ("met", "salinity_at_3m", 3.0),
    "tsg5": ("met", "salinity_at_5m", 5.0),
}
PROFILE_LAYERS = ["snake", "ssp005", "ssp012", "ssp023", "ssp054", "ssp110", "usps2", "usps3", "tsg5"]
MOOR = {"sal": "SPURS2_MOORING_CENTRAL/SPURS2_CentralMooring_SubsurfaceTemperatureSalinity.nc",
        "met": "SPURS2_MOORING_CENTRAL/SPURS2_WHOI_Met.nc"}
PRIMARY = {"name": "primary", "dry_h": 6, "thr": 10.0, "win": (0, 6), "sfc": "snake", "z_sfc": 0.0,
           "t_cur": p2.T_CURRENT_S, "disp_max": None, "min_hours": 1}
SENS = [
    dict(PRIMARY, name="S1_dry24", dry_h=24),
    dict(PRIMARY, name="S2_dry12", dry_h=12),
    dict(PRIMARY, name="S3_thr5", thr=5.0),
    dict(PRIMARY, name="S4_win0_3", win=(0, 3)),
    dict(PRIMARY, name="S5_win0_12", win=(0, 12)),
    dict(PRIMARY, name="S6_disp30", disp_max=30.0),
    dict(PRIMARY, name="S7_snake_met", sfc="snake_met"),
    dict(PRIMARY, name="S8_z0015", z_sfc=0.015),
    dict(PRIMARY, name="S9_tcur1800", t_cur=1800.0),
    dict(PRIMARY, name="S10_min4h", min_hours=4),
]
LOG = []


# ======================================================================== 小工具
def log(msg):
    line = f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} {msg}"
    LOG.append(line)
    print(line, flush=True)


def isnum(x):
    return isinstance(x, (int, float, np.floating)) and math.isfinite(float(x))


def rnd(x, nd=5):
    return round(float(x), nd) if isnum(x) else None


def iso(h):
    return datetime.fromtimestamp(int(h) * 3600, timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha256_file(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jdump(obj, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1, default=str)
    os.replace(tmp, path)


def haversine(la1, lo1, la2, lo2):
    p = math.pi / 180
    a = math.sin((la2 - la1) * p / 2) ** 2 + math.cos(la1 * p) * math.cos(la2 * p) * math.sin((lo2 - lo1) * p / 2) ** 2
    return 2 * 6371.0 * math.asin(min(1.0, math.sqrt(a)))


# ======================================================================== 读原件（只读）与小时化（W1–W3、W5）
def nc_time(ds, name="time"):
    import netCDF4
    tv = ds.variables[name]
    raw = np.asarray(np.ma.filled(tv[:].astype(float), np.nan), float)
    b0 = netCDF4.num2date(0.0, tv.units, getattr(tv, "calendar", "standard"),
                          only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    b1 = netCDF4.num2date(1.0, tv.units, getattr(tv, "calendar", "standard"),
                          only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    t0 = b0.replace(tzinfo=timezone.utc).timestamp()
    sc = b1.replace(tzinfo=timezone.utc).timestamp() - t0
    return t0 + raw * sc


def nc_var(ds, name):
    v = ds.variables[name]
    v.set_auto_maskandscale(True)
    return np.asarray(np.ma.filled(v[:].astype(float), np.nan), float)


def bin_mean(t, x, step_s, min_frac):
    """按 step_s 秒箱求均值；返回 (k0, 均值数组)，k＝⌊t/step⌋。有效样本 <min_frac×期望数的箱为 NaN。"""
    ok = np.isfinite(t) & np.isfinite(x)
    t, x = t[ok], x[ok]
    if len(t) < 2:
        return None, np.array([])
    dt = float(np.median(np.diff(np.sort(t))))
    exp_n = max(1.0, step_s / dt) if dt > 0 else 1.0
    k = np.floor(t / step_s).astype(np.int64)
    k0 = int(k.min())
    n = int(k.max()) - k0 + 1
    cnt = np.bincount(k - k0, minlength=n)
    sm = np.bincount(k - k0, weights=x, minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = sm / cnt
    m[cnt < min_frac * exp_n] = np.nan
    return k0, m


class Series:
    """整数时间索引（小时或半步）上的序列，越界读为 NaN。"""

    def __init__(self, k0, arr):
        self.k0 = int(k0) if k0 is not None else 0
        self.a = np.asarray(arr, float)

    def get(self, k_lo, k_hi):
        out = np.full(k_hi - k_lo, np.nan)
        a, b = max(k_lo, self.k0), min(k_hi, self.k0 + len(self.a))
        if b > a:
            out[a - k_lo:b - k_lo] = self.a[a - self.k0:b - self.k0]
        return out

    def span(self):
        return self.k0, self.k0 + len(self.a)


def load_all(root, cache, log_fn):
    """读原件并小时化；结果缓存到 cache/hourly_all.npz（键＝平台|量）。返回 {键: Series}。"""
    import netCDF4
    os.makedirs(cache, exist_ok=True)
    cp = os.path.join(cache, "hourly_all.npz")
    src_sig = {rel: os.path.getsize(os.path.join(root, rel)) for c in CRUISES.values() for rel in c.values()}
    src_sig.update({rel: os.path.getsize(os.path.join(root, rel)) for rel in MOOR.values()})
    if os.path.exists(cp):
        z = np.load(cp, allow_pickle=False)
        meta = json.loads(str(z["__meta__"]))
        if meta.get("src_sig") == src_sig and meta.get("version") == VERSION:
            log_fn(f"cache hit {cp}")
            return {k: Series(int(z[k + "::k0"]), z[k]) for k in meta["keys"]}, meta
    S, info = {}, {"files": {}}
    for cr, fs in CRUISES.items():
        dsm = netCDF4.Dataset(os.path.join(root, fs["met"]))
        tm = nc_time(dsm)
        r = np.maximum(nc_var(dsm, "precipitation_rate"), 0.0)
        k0, hr = bin_mean(tm, r, 3600.0, MIN_FRAC_RAIN)
        S[f"{cr}|rain_h"] = Series(k0, hr)
        k0, hs = bin_mean(tm, r, 1800.0, MIN_FRAC_RAIN)
        S[f"{cr}|rain_hs"] = Series(k0, hs)
        k0, w = bin_mean(tm, nc_var(dsm, "neutral_wind_speed_relative_to_earth"), 3600.0, MIN_FRAC_RAIN)
        S[f"{cr}|wind_h"] = Series(k0, w)
        k0, la = bin_mean(tm, nc_var(dsm, "latitude"), 3600.0, MIN_FRAC_RAIN)
        S[f"{cr}|lat_h"] = Series(k0, la)
        k0, lo = bin_mean(tm, nc_var(dsm, "longitude"), 3600.0, MIN_FRAC_RAIN)
        S[f"{cr}|lon_h"] = Series(k0, lo)
        opened = {"met": (dsm, tm)}
        for fk in ("snake", "ssp"):
            ds = netCDF4.Dataset(os.path.join(root, fs[fk]))
            opened[fk] = (ds, nc_time(ds))
        for lname, (fk, vname, _z) in LAYERS.items():
            ds, tt = opened[fk]
            k0, m = bin_mean(tt, nc_var(ds, vname), 3600.0, MIN_FRAC_SAL)
            S[f"{cr}|{lname}"] = Series(k0, m)
        for ds, _ in opened.values():
            ds.close()
        info["files"][cr] = fs
    ds = netCDF4.Dataset(os.path.join(root, MOOR["met"]))
    tm = nc_time(ds)
    r = np.maximum(nc_var(ds, "rainfall_rate"), 0.0)
    k0, hr = bin_mean(tm, r, 3600.0, MIN_FRAC_RAIN)
    S["MOOR|rain_h"] = Series(k0, hr)
    k0, hs = bin_mean(tm, r, 1800.0, MIN_FRAC_RAIN)
    S["MOOR|rain_hs"] = Series(k0, hs)
    spd = np.hypot(nc_var(ds, "ewind"), nc_var(ds, "nwind"))
    k0, w = bin_mean(tm, spd, 3600.0, MIN_FRAC_RAIN)
    S["MOOR|wind_h"] = Series(k0, w * WIND_FACTOR_MOOR)
    ds.close()
    ds = netCDF4.Dataset(os.path.join(root, MOOR["sal"]))
    tt = nc_time(ds)
    depth = nc_var(ds, "depth")
    j = int(np.argmin(np.abs(depth - 1.0)))
    if abs(depth[j] - 1.0) > 1e-6:
        raise RuntimeError(f"浮标无 1.0 m 层（最近 {depth[j]}）")
    sal = ds.variables["salinity"]
    dims = sal.dimensions
    arr = nc_var(ds, "salinity")
    col = arr[j, :] if dims[0] == "depth" else arr[:, j]
    k0, m = bin_mean(tt, col, 3600.0, MIN_FRAC_SAL)
    S["MOOR|s1"] = Series(k0, m)
    info["moor_depth_index"] = j
    ds.close()
    meta = {"version": VERSION, "src_sig": src_sig, "keys": sorted(S), "info": info}
    np.savez(cp, __meta__=json.dumps(meta), **{k: s.a for k, s in S.items()}, **{k + "::k0": np.int64(s.k0) for k, s in S.items()})
    log_fn(f"cache written {cp}")
    return S, meta


# ======================================================================== 事件（W4）
def find_events(rain_h, dry_h, thr):
    """rain_h：Series（小时雨 mm/h）。返回 onset 小时列表与累积量（P1 I6，干期 dry_h）。"""
    k0, rr = rain_h.k0, rain_h.a
    valid = np.isfinite(rr)
    wet = valid & (rr >= RAIN_HOUR_MM)
    out = []
    for i in range(dry_h, len(rr) - 24 + 1):
        if not wet[i]:
            continue
        if not (valid[i - dry_h:i].all() and not wet[i - dry_h:i].any()):
            continue
        if not valid[i:i + 24].all():
            continue
        acc = float(np.where(wet[i:i + 24], rr[i:i + 24], 0.0).sum())
        if acc >= thr:
            out.append((k0 + i, acc))
    return out


# ======================================================================== 模型（W5、W6）
def forcing(S, plat, H, wind_factor=1.0):
    """[H−30, H+12) 的半步雨 P 与半步风 U（长 84）；返回 (P, U, diag)。"""
    lo, hi = H - HIST_H, H + POST_MAX_H
    rs = S[f"{plat}|rain_hs"]
    a, b = rs.span()
    if 2 * lo < a or 2 * hi > b:
        return None, None, {"reason": "rain_record_short"}
    P = rs.get(2 * lo, 2 * hi)
    n_fill = int((~np.isfinite(P)).sum())
    P = np.where(np.isfinite(P), P, 0.0)
    pad = p2.WIND_GAP_MAX_H
    w = S[f"{plat}|wind_h"].get(lo - pad, hi + pad)
    w = p2.fill_wind(w, p2.WIND_GAP_MAX_H)[pad:pad + (hi - lo)] * wind_factor
    if not np.all(np.isfinite(w)):
        return None, None, {"reason": "wind_gap"}
    n_floor = int((w < p2.WIND_FLOOR).sum())
    w = np.maximum(w, p2.WIND_FLOOR)
    return P, np.repeat(w, 2), {"reason": None, "rain_halfsteps_filled0": n_fill, "wind_floor_hours": n_floor,
                                "p_gt_d0_range": int((P > p2.R_GRID[-1]).sum())}


def model_hourly_F(P, U, z, t_cur=p2.T_CURRENT_S):
    """W6：小时 t0−6…t0+11 的 F(z)（18 个）。"""
    F = p2.rim_factor(P, U, z, t_cur_depth=t_cur)
    return F.reshape(-1, 2).mean(axis=1)


# ======================================================================== ΔS（W8、W13、W14）
def hours_rel(lo, hi):
    """相对 t0 的小时 [lo, hi) → 模型小时数组下标（模型下标 0 ↔ t0−6）。"""
    return np.arange(lo + PRE_H, hi + PRE_H)


def delta_pair(obs18, F18, win, min_hours=1, mask_post=None):
    """obs18、F18：小时 t0−6…t0+11。返回 (ΔS_obs, ΔS_RIM, S0, n_pre, n_post)；无效返回 NaN。"""
    pre_i = np.arange(0, PRE_H)
    post_i = hours_rel(*win)
    pre_ok = pre_i[np.isfinite(obs18[pre_i]) & np.isfinite(F18[pre_i])]
    post_m = np.isfinite(obs18[post_i]) & np.isfinite(F18[post_i])
    if mask_post is not None:
        post_m &= mask_post
    post_ok = post_i[post_m]
    if len(pre_ok) < min_hours or len(post_ok) < min_hours:
        return NAN, NAN, NAN, len(pre_ok), len(post_ok)
    s0 = float(np.median(obs18[pre_ok]))
    dso = float(obs18[post_ok].mean() - s0)
    dsm = s0 * float(F18[post_ok].mean() - np.median(F18[pre_ok]))
    return dso, dsm, s0, len(pre_ok), len(post_ok)


def event_block(S, plat, H, layers, cfg, rain_mask12=None, wind_factor=1.0):
    """一个事件、一组层的观测与模型 ΔS。layers：{层名: (序列键, 模型深度)}。"""
    P, U, diag = forcing(S, plat, H, wind_factor)
    rec = {"diag": diag, "layers": {}}
    if P is None:
        return rec
    Fcache = {}
    for lname, (skey, z) in layers.items():
        zz = cfg["z_sfc"] if lname in ("snake", "snake_met") else z
        tc = cfg["t_cur"]
        if (zz, tc) not in Fcache:
            Fcache[(zz, tc)] = model_hourly_F(P, U, zz, tc)
        F18 = Fcache[(zz, tc)]
        obs18 = S[skey].get(H - PRE_H, H + POST_MAX_H)
        dso, dsm, s0, npre, npost = delta_pair(obs18, F18, cfg["win"], cfg["min_hours"])
        rec["layers"][lname] = {"obs": dso, "rim": dsm, "s0": s0, "n_pre": npre, "n_post": npost,
                                "obs18": obs18, "F18": F18}
    return rec


# ======================================================================== 统计（W10、W11）
def boot(recs, pairs, key_f):
    """pairs：{名: (num 函数, den 函数)}；返回 (bootstrap 数组 dict, 簇数)。"""
    keys = [key_f(r) for r in recs]
    arrs = {k: (np.array([f(r) for r in recs], float), np.array([g(r) for r in recs], float)) for k, (f, g) in pairs.items()}
    return p2.cluster_boot(keys, arrs)


def pct(a):
    a = a[np.isfinite(a)]
    if len(a) < 100:
        return [None, None]
    return [rnd(np.percentile(a, 2.5)), rnd(np.percentile(a, 97.5))]


def ratio_block(recs, num_f, den_f, clus_f, extra=None):
    n = len(recs)
    out = {"n": n}
    if n == 0:
        return dict(out, evaluable=False, point=None)
    num = sum(num_f(r) for r in recs)
    den = sum(den_f(r) for r in recs)
    out.update({"sum_obs": rnd(num), "sum_rim": rnd(den), "point": rnd(num / den) if den else None})
    pairs = {"R": (num_f, den_f)}
    if extra:
        pairs.update(extra)
    bc, K = boot(recs, pairs, clus_f)
    be, _ = boot(recs, {"R": (num_f, den_f)}, lambda r: r["uid"])
    dsum = np.array([den_f(r) for r in recs], float)
    idx = {}
    for r in recs:
        idx.setdefault(clus_f(r), []).append(r)
    # 分母 ≥0 的比例：用同一簇抽样重算
    uk = sorted(idx)
    kid = {k: j for j, k in enumerate(uk)}
    cid = np.array([kid[clus_f(r)] for r in recs])
    dk = np.bincount(cid, weights=dsum, minlength=len(uk))
    rng = np.random.Generator(np.random.PCG64(SEED))
    draw = rng.integers(0, len(uk), size=(B, len(uk)))
    frac_nonneg = float((dk[draw].sum(axis=1) >= 0).mean())
    out.update({"clusters": K, "ci95_cluster": pct(bc["R"]), "ci95_event": pct(be["R"]),
                "frac_boot_den_nonneg": rnd(frac_nonneg, 4), "_boot": bc})
    return out


def loo_cruise(recs, num_f, den_f):
    out = {}
    for cr in sorted({r["cruise"] for r in recs}):
        rs = [r for r in recs if r["cruise"] != cr]
        d = sum(den_f(r) for r in rs)
        out[f"drop_{cr}"] = rnd(sum(num_f(r) for r in rs) / d) if rs and d else None
    return out


def category(lo, hi):
    if lo is None or hi is None:
        return "不可评"
    if hi < THRESH:
        return "海面未见淡化（异常）：不支持 R2" if hi < 0 else "支持 R1（海面也显著小于 RIM-3）"
    if lo > THRESH:
        return "支持 R2（海面与 RIM-3 相符或更大）"
    return "不能区分"


def evaluable(blk, recs):
    per = {}
    for r in recs:
        per[r["cruise"]] = per.get(r["cruise"], 0) + 1
    conds = {"n>=15": blk["n"] >= 15, "clusters>=4": (blk.get("clusters") or 0) >= 4,
             "each_cruise>=3": len(per) == len(CRUISES) and min(per.values()) >= 3,
             "frac_den_nonneg<=0.01": blk.get("frac_boot_den_nonneg") is not None and blk["frac_boot_den_nonneg"] <= DEN_NONNEG_MAX}
    return all(conds.values()), conds, per


def judge(blk, ok, loo, rel_cat):
    """W12：按序取第一条成立者；稳健标注；平流标注。"""
    if not ok:
        return {"category": "不可评", "rule": "第六节第 1 条", "robust": None, "advection": None}
    lo, hi = blk["ci95_cluster"]
    if lo is None:
        return {"category": "不可评", "rule": "第六节第 1 条（bootstrap 不可算）", "robust": None, "advection": None}
    cat = category(lo, hi)
    rule = {"支持 R1（海面也显著小于 RIM-3）": "第六节第 2 条", "海面未见淡化（异常）：不支持 R2": "第六节第 2 条（上端 <0）",
            "支持 R2（海面与 RIM-3 相符或更大）": "第六节第 3 条", "不能区分": "第六节第 4 条"}[cat]
    robust = None
    if cat != "不能区分":
        elo, ehi = blk["ci95_event"]
        lv = [v for v in loo.values() if v is not None]
        if cat.startswith("支持 R2"):
            robust = elo is not None and elo > THRESH and len(lv) == len(loo) and min(lv) > THRESH
        else:
            robust = ehi is not None and ehi < THRESH and len(lv) == len(loo) and max(lv) < THRESH
        robust = "稳健" if robust else "依赖簇定义"
    adv = None if rel_cat is None else ("平流订正同类" if rel_cat == cat else f"依赖平流订正（R_rel5 为「{rel_cat}」）")
    return {"category": cat, "rule": rule, "robust": robust, "advection": adv}


def judge_rain(blk):
    """W13：雨中组 R_sfc,rain 与 1 比较。"""
    lo, hi = blk.get("ci95_cluster") or [None, None]
    if blk.get("n", 0) < 8 or lo is None:
        return "不可评（雨中组 n<8 或 CI 不可算）"
    if hi < 1.0:
        return "G1：现场雨中海面淡化小于 RIM-3，P5 雨中过境比值 1.98 不被现场复现"
    if lo > 1.0:
        return "G2：现场雨中海面淡化也超过 RIM-3，P5 雨中超出可以是真实海面信号"
    return "G3：不能区分"


# ======================================================================== 事件集构建
def build_ship(S, cfg):
    recs, drop = [], {}
    uid = 0
    for cr in CRUISES:
        evs = find_events(S[f"{cr}|rain_h"], cfg["dry_h"], cfg["thr"])
        for H, acc in evs:
            layers = {ln: (f"{cr}|{ln}", LAYERS[ln][2]) for ln in PROFILE_LAYERS + ["snake_met"]}
            b = event_block(S, cr, H, layers, cfg)
            if b["diag"].get("reason"):
                drop[b["diag"]["reason"]] = drop.get(b["diag"]["reason"], 0) + 1
                continue
            la1, lo1 = np.nanmedian(S[f"{cr}|lat_h"].get(H - 6, H)), np.nanmedian(S[f"{cr}|lon_h"].get(H - 6, H))
            la2, lo2 = np.nanmedian(S[f"{cr}|lat_h"].get(H, H + 6)), np.nanmedian(S[f"{cr}|lon_h"].get(H, H + 6))
            disp = haversine(la1, lo1, la2, lo2) if all(np.isfinite([la1, lo1, la2, lo2])) else NAN
            if cfg["disp_max"] is not None and not (np.isfinite(disp) and disp <= cfg["disp_max"]):
                drop["disp"] = drop.get("disp", 0) + 1
                continue
            iy, iw, _ = datetime.fromtimestamp(H * 3600, timezone.utc).isocalendar()
            rain12 = S[f"{cr}|rain_h"].get(H, H + POST_MAX_H)
            uid += 1
            recs.append({"uid": uid, "cruise": cr, "H": H, "onset": iso(H), "acc24": acc, "cluster": f"{cr}-W{iw:02d}",
                         "disp_km": disp, "diag": b["diag"], "L": b["layers"], "rain12": rain12})
    return recs, drop


def build_moor(S, cfg):
    recs, drop = [], {}
    evs = find_events(S["MOOR|rain_h"], cfg["dry_h"], cfg["thr"])
    for i, (H, acc) in enumerate(evs):
        b = event_block(S, "MOOR", H, {"s1": ("MOOR|s1", 1.0)}, dict(cfg, z_sfc=0.0))
        if b["diag"].get("reason"):
            drop[b["diag"]["reason"]] = drop.get(b["diag"]["reason"], 0) + 1
            continue
        d = datetime.fromtimestamp(H * 3600, timezone.utc)
        recs.append({"uid": i + 1, "cruise": "MOOR", "H": H, "onset": iso(H), "acc24": acc,
                     "cluster": f"{d.year}-{d.month:02d}", "L": b["layers"]})
    return recs, drop


def lay_ok(r, ln):
    L = r["L"].get(ln)
    return L is not None and isnum(L["obs"]) and isnum(L["rim"])


def f_obs(ln):
    return lambda r: r["L"][ln]["obs"]


def f_rim(ln):
    return lambda r: r["L"][ln]["rim"]


def clus(r):
    return r["cluster"]


def strip(b):
    return {k: v for k, v in b.items() if not k.startswith("_")}


# ======================================================================== 专项量
def rel5_block(recs, cfg):
    """W14：snake − 5 m 差分。"""
    rows = []
    for r in recs:
        Ls, L5 = r["L"].get(cfg["sfc"]), r["L"].get("tsg5")
        if Ls is None or L5 is None:
            continue
        o = Ls["obs18"] - L5["obs18"]
        s0s = np.nanmedian(Ls["obs18"][:PRE_H]) if np.isfinite(Ls["obs18"][:PRE_H]).any() else NAN
        s05 = np.nanmedian(L5["obs18"][:PRE_H]) if np.isfinite(L5["obs18"][:PRE_H]).any() else NAN
        if not (isnum(s0s) and isnum(s05)):
            continue
        m = s0s * Ls["F18"] - s05 * L5["F18"]
        pre_i = np.arange(0, PRE_H)
        post_i = hours_rel(*cfg["win"])
        pre_ok = pre_i[np.isfinite(o[pre_i])]
        post_ok = post_i[np.isfinite(o[post_i])]
        if len(pre_ok) < cfg["min_hours"] or len(post_ok) < cfg["min_hours"]:
            continue
        rows.append(dict(r, d_obs=float(o[post_ok].mean() - np.median(o[pre_ok])),
                         d_rim=float(m[post_ok].mean() - np.median(m[pre_ok]))))
    blk = ratio_block(rows, lambda r: r["d_obs"], lambda r: r["d_rim"], clus)
    blk["category"] = category(*blk["ci95_cluster"]) if blk.get("ci95_cluster") else "不可评"
    return blk


def rain_split(recs, cfg):
    """W13：[t0, t0+12) 内 snake 有值小时按雨中／雨停分组。"""
    out = {}
    for grp in ("rain", "dry"):
        rows = []
        for r in recs:
            L = r["L"].get(cfg["sfc"])
            if L is None:
                continue
            wet = np.isfinite(r["rain12"]) & (r["rain12"] >= RAIN_HOUR_MM)
            mask = wet if grp == "rain" else (np.isfinite(r["rain12"]) & ~wet)
            dso, dsm, s0, npre, npost = delta_pair(L["obs18"], L["F18"], (0, POST_MAX_H), 1, mask_post=mask)
            if isnum(dso) and isnum(dsm):
                rows.append(dict(r, g_obs=dso, g_rim=dsm, g_hours=npost))
        out[grp] = rows
    res = {}
    for grp, rows in out.items():
        blk = ratio_block(rows, lambda r: r["g_obs"], lambda r: r["g_rim"], clus)
        blk["hours"] = int(sum(r["g_hours"] for r in rows))
        res[grp] = blk
    both = [r for r in out["rain"]]
    dry_by = {r["uid"]: r for r in out["dry"]}
    paired = [dict(r, d_obs=dry_by[r["uid"]]["g_obs"], d_rim=dry_by[r["uid"]]["g_rim"]) for r in both if r["uid"] in dry_by]
    if paired:
        bc, K = boot(paired, {"rain": (lambda r: r["g_obs"], lambda r: r["g_rim"]),
                              "dry": (lambda r: r["d_obs"], lambda r: r["d_rim"])}, clus)
        with np.errstate(invalid="ignore", divide="ignore"):
            lr = np.log(bc["rain"] / bc["dry"])
        num_r = sum(r["g_obs"] for r in paired) / sum(r["g_rim"] for r in paired)
        num_d = sum(r["d_obs"] for r in paired) / sum(r["d_rim"] for r in paired)
        res["paired"] = {"n": len(paired), "clusters": K, "R_rain": rnd(num_r), "R_dry": rnd(num_d),
                         "ln_ratio": rnd(math.log(num_r / num_d)) if num_r > 0 and num_d > 0 else None,
                         "ln_ratio_ci95": pct(lr)}
    res["reading"] = judge_rain(res["rain"])
    return {k: strip(v) if isinstance(v, dict) else v for k, v in res.items()}


def layer_q_block(recs, sfc):
    """W15：snake 与 SSP 1.10 m 都有效的事件；同一抽样下 R_sfc、R_1m、Q、观测层比与 RIM 层比。"""
    rows = [r for r in recs if lay_ok(r, sfc) and lay_ok(r, "ssp110")]
    out = {"n": len(rows)}
    if len(rows) < 3:
        return dict(out, note="事件 <3，不算")
    so, sm = sum(r["L"][sfc]["obs"] for r in rows), sum(r["L"][sfc]["rim"] for r in rows)
    oo, om = sum(r["L"]["ssp110"]["obs"] for r in rows), sum(r["L"]["ssp110"]["rim"] for r in rows)
    bc, K = boot(rows, {"Rs": (f_obs(sfc), f_rim(sfc)), "R1": (f_obs("ssp110"), f_rim("ssp110")),
                        "LRo": (f_obs(sfc), f_obs("ssp110")), "LRm": (f_rim(sfc), f_rim("ssp110"))}, clus)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = bc["Rs"] / bc["R1"]
    out.update({"clusters": K, "R_sfc": rnd(so / sm), "R_sfc_ci95": pct(bc["Rs"]), "R_1m": rnd(oo / om), "R_1m_ci95": pct(bc["R1"]),
                "Q": rnd((so / sm) / (oo / om)) if oo and om else None, "Q_ci95": pct(q),
                "layer_ratio_obs": rnd(so / oo) if oo else None, "layer_ratio_obs_ci95": pct(bc["LRo"]),
                "layer_ratio_rim": rnd(sm / om) if om else None, "layer_ratio_rim_ci95": pct(bc["LRm"])})
    rows2 = [r for r in recs if lay_ok(r, "ssp012") and lay_ok(r, "ssp110")]
    if len(rows2) >= 3:
        b2, K2 = boot(rows2, {"Ro": (f_obs("ssp012"), f_obs("ssp110")), "Rm": (f_rim("ssp012"), f_rim("ssp110")),
                              "R012": (f_obs("ssp012"), f_rim("ssp012"))}, clus)
        a_o = sum(r["L"]["ssp012"]["obs"] for r in rows2) / sum(r["L"]["ssp110"]["obs"] for r in rows2)
        a_m = sum(r["L"]["ssp012"]["rim"] for r in rows2) / sum(r["L"]["ssp110"]["rim"] for r in rows2)
        out["ssp_internal_012_vs_110"] = {"n": len(rows2), "clusters": K2, "layer_ratio_obs": rnd(a_o), "ci95": pct(b2["Ro"]),
                                          "layer_ratio_rim": rnd(a_m), "rim_ci95": pct(b2["Rm"]),
                                          "R_012": rnd(sum(r["L"]["ssp012"]["obs"] for r in rows2) / sum(r["L"]["ssp012"]["rim"] for r in rows2)),
                                          "R_012_ci95": pct(b2["R012"])}
    return out


def power_block(recs, S, cfg):
    """W17：只用模型 x 与无雨对照窗观测噪声；在任何事件观测 ΔS 合计之前写盘。"""
    xs = np.array([r["L"][cfg["sfc"]]["rim"] for r in recs if lay_ok(r, cfg["sfc"])], float)
    n = len(xs)
    K = len({r["cluster"] for r in recs if lay_ok(r, cfg["sfc"])})
    ctrl = []
    for cr in CRUISES:
        rh = S[f"{cr}|rain_h"]
        a, b = rh.span()
        for c in range(a + 12, b - 6, 6):
            seg = rh.get(c - 12, c + 6)
            if not (np.isfinite(seg).all() and (seg < RAIN_HOUR_MM).all()):
                continue
            obs18 = S[f"{cr}|{cfg['sfc']}"].get(c - PRE_H, c + POST_MAX_H)
            dso, _, _, _, _ = delta_pair(obs18, np.zeros(18), cfg["win"], 1)
            if isnum(dso):
                ctrl.append(dso)
    ctrl = np.array(ctrl, float)
    sx = float(xs.sum()) if n else NAN
    out = {"n": n, "clusters": K, "sum_x": rnd(sx), "mean_x": rnd(xs.mean()) if n else None,
           "rms_x": rnd(np.sqrt((xs ** 2).mean())) if n else None,
           "n_ctrl_windows": int(len(ctrl)), "sd_ctrl_dS": rnd(ctrl.std(ddof=1)) if len(ctrl) > 2 else None,
           "note": "只含模型侧 x 与无雨对照窗 ΔS 的离散度（观测噪声代理，含船动带来的水平梯度）；只作信息，不改判读"}
    if n and abs(sx) > 0:
        for sg in (0.1, 0.2, 0.4):
            out[f"se_R_if_sigma_{sg}"] = rnd(sg * math.sqrt(n) / abs(sx))
        if out["sd_ctrl_dS"]:
            se = out["sd_ctrl_dS"] * math.sqrt(n) / abs(sx)
            out["se_R_ctrl"] = rnd(se)
            out["ci_halfwidth_ctrl"] = rnd(1.96 * se)
            out["separable_0p3_vs_1_at_0p65"] = bool(1.96 * se < 0.35)
    return out


# ======================================================================== 主流程
def run_cfg(S, cfg, full=True):
    recs, drop = build_ship(S, cfg)
    sfc = cfg["sfc"]
    rows = [r for r in recs if lay_ok(r, sfc)]
    blk = ratio_block(rows, f_obs(sfc), f_rim(sfc), clus)
    ok, conds, per = evaluable(blk, rows) if rows else (False, {}, {})
    loo = loo_cruise(rows, f_obs(sfc), f_rim(sfc))
    rel = rel5_block(recs, cfg)
    j = judge(blk, ok, loo, rel.get("category"))
    res = {"cfg": {k: v for k, v in cfg.items()}, "n_events_all": len(recs), "dropped": drop,
           "R_sfc": dict(strip(blk), evaluable=ok, evaluable_conditions=conds, per_cruise=per, loo_cruise=loo),
           "R_rel5": strip(rel), "judgement": j}
    if full:
        return res, recs
    return res, None


def summarize_layers(recs, cfg):
    prof = {}
    for ln in PROFILE_LAYERS:
        rows = [r for r in recs if lay_ok(r, ln)]
        if len(rows) >= 3:
            b = ratio_block(rows, f_obs(ln), f_rim(ln), clus)
            prof[ln] = dict(strip(b), z_model=LAYERS[ln][2] if ln != "snake" else cfg["z_sfc"])
        else:
            prof[ln] = {"n": len(rows), "note": "事件 <3，不算"}
    return prof


def events_csv(recs, path):
    cols = ["cruise", "onset", "acc24", "cluster", "disp_km"]
    lays = PROFILE_LAYERS + ["snake_met"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols + [f"{ln}_{k}" for ln in lays for k in ("obs", "rim", "s0", "npre", "npost")] + ["rain_hours_0_12"])
        for r in recs:
            row = [r["cruise"], r["onset"], rnd(r["acc24"], 2), r["cluster"], rnd(r["disp_km"], 1)]
            for ln in lays:
                L = r["L"].get(ln, {})
                row += [rnd(L.get("obs")), rnd(L.get("rim")), rnd(L.get("s0"), 4), L.get("n_pre"), L.get("n_post")]
            row.append(int(np.nansum(r["rain12"] >= RAIN_HOUR_MM)))
            w.writerow(row)


def run_full(a):
    t_start = time.time()
    out = a.out
    os.makedirs(out, exist_ok=True)
    st = selftest()
    jdump(st, os.path.join(out, "p6_selftest.json"))
    if not st["all_ok"]:
        log("自测不过，退出 4")
        return 4
    here = os.path.dirname(os.path.abspath(__file__))
    code_sha = {fn: sha256_file(os.path.join(here, fn)) for fn in ("p6_spurs2_profile.py", "p2_rim_test.py", "p1_events.py", "p1b_extend.py")}
    log(f"code sha {json.dumps(code_sha)}")
    S, meta = load_all(a.root, a.cache, log)
    log(f"series {len(S)}")
    # ---- 主事件集与功效（先写）
    recs, drop = build_ship(S, PRIMARY)
    log(f"primary events {len(recs)} dropped {drop}")
    pw = power_block(recs, S, PRIMARY)
    pw.update({"version": VERSION, "code_sha256": code_sha, "written": datetime.now(timezone.utc).isoformat()})
    jdump(pw, os.path.join(out, "p6_power.json"))
    log(f"power written {json.dumps({k: v for k, v in pw.items() if k != 'code_sha256'}, ensure_ascii=False)}")
    time.sleep(1.0)
    # ---- 主结果与判读
    main_res, _ = run_cfg(S, PRIMARY)
    prof = summarize_layers(recs, PRIMARY)
    rain = rain_split(recs, PRIMARY)
    lq = layer_q_block(recs, "snake")
    rows1 = [r for r in recs if lay_ok(r, "ssp110")]
    r1_ssp = strip(ratio_block(rows1, f_obs("ssp110"), f_rim("ssp110"), clus)) if rows1 else {"n": 0}
    moor_cfg = dict(PRIMARY, dry_h=24)          # 浮标：P1 事件定义原样（干期 24 h）
    mrecs, mdrop = build_moor(S, moor_cfg)
    mrows = [r for r in mrecs if lay_ok(r, "s1")]
    r1_moor = strip(ratio_block(mrows, f_obs("s1"), f_rim("s1"), clus)) if mrows else {"n": 0}
    r1_moor.update({"events_all": len(mrecs), "dropped": mdrop, "def": "P1 原样：10 mm、干期 24 h、0.4 mm/h；窗 [0,6) h；簇＝日历月"})
    r1_moor["reading"] = ("与 P2 同侧（CI 上端 <0.65）" if r1_moor.get("ci95_cluster") and r1_moor["ci95_cluster"][1] is not None
                          and r1_moor["ci95_cluster"][1] < THRESH else "未落在 P2 同侧（CI 上端 ≥0.65 或不可算）：本海区 1 m 与 P2 不一致，须在结果中说明")
    mrows6 = [r for r in build_moor(S, PRIMARY)[0] if lay_ok(r, "s1")]
    r1_moor_dry6 = strip(ratio_block(mrows6, f_obs("s1"), f_rim("s1"), clus)) if mrows6 else {"n": 0}
    sens = {}
    for c in SENS:
        res, _ = run_cfg(S, c, full=False)
        j = res["judgement"]
        sens[c["name"]] = {"n": res["R_sfc"]["n"], "clusters": res["R_sfc"].get("clusters"), "point": res["R_sfc"].get("point"),
                           "ci95_cluster": res["R_sfc"].get("ci95_cluster"), "evaluable": res["R_sfc"]["evaluable"],
                           "category": j["category"], "same_as_primary": j["category"] == main_res["judgement"]["category"],
                           "R_rel5": {"n": res["R_rel5"].get("n"), "point": res["R_rel5"].get("point"),
                                      "ci95_cluster": res["R_rel5"].get("ci95_cluster"), "category": res["R_rel5"].get("category")}}
    summ = {"version": VERSION, "code_sha256": code_sha, "data_meta": meta.get("info"),
            "primary": main_res, "verdict": main_res["judgement"], "profile_R_by_layer": prof,
            "rain_vs_dry": rain, "layer_ratio_and_Q_ssp_subset": lq, "R_1m_ssp": r1_ssp,
            "R_1m_mooring": r1_moor, "R_1m_mooring_dry6_info": r1_moor_dry6, "sensitivity": sens,
            "disp_km_primary": {"median": rnd(np.nanmedian([r["disp_km"] for r in recs])) if recs else None,
                                "max": rnd(np.nanmax([r["disp_km"] for r in recs])) if recs else None},
            "elapsed_s": round(time.time() - t_start, 1)}
    jdump(summ, os.path.join(out, "p6_summary.json"))
    events_csv(recs, os.path.join(out, "p6_events.csv"))
    log(f"verdict {json.dumps(main_res['judgement'], ensure_ascii=False)}")
    log("done rc=0")
    return 0


# ======================================================================== 自测
def selftest():
    res = {}
    rng = np.random.default_rng(7)
    # T1：W6 模型小时 F 与 p2.model_hourly 逐值一致（z 0.5/1/5）
    P = np.where(rng.random(84) < 0.3, rng.gamma(1.5, 6.0, 84), 0.0)
    U = np.repeat(rng.uniform(1, 11, 42), 2)
    ref, _ = p2.model_hourly(P, U, 34.0, {0.5: 0.5, 1.0: 1.0, 5.0: 5.0})
    mx = max(float(np.nanmax(np.abs(model_hourly_F(P, U, z) - ref[z]))) for z in (0.5, 1.0, 5.0))
    res["T1_model_equals_p2_model_hourly"] = mx == 0.0
    # T2：z=0 与 p2._rim_literal（Witte 代码逐行移植）一致
    try:
        lit = p2._rim_literal(P, U, np.full(len(P), 34.0))
        F0 = p2.rim_factor(P, U, 0.0)
        res["T2_z0_equals_literal"] = bool(np.nanmax(np.abs(np.asarray(lit) / 34.0 - F0)) < 1e-12)
    except Exception as e:
        res["T2_z0_equals_literal"] = f"error {e!r}"
    # T3：bin_mean 与覆盖门限
    t = np.arange(0, 7200, 1.0)
    x = np.ones_like(t)
    x[3600:3600 + 2600] = np.nan                      # 第二小时剩 1000/3600＝28% ≥25% → 有效
    k0, m = bin_mean(t, x, 3600.0, MIN_FRAC_SAL)
    res["T3_bin_mean"] = k0 == 0 and m[0] == 1.0 and m[1] == 1.0
    x[3600:3600 + 2800] = np.nan                      # 剩 800/3600＝22% <25% → 缺
    k0, m = bin_mean(t, x, 3600.0, MIN_FRAC_SAL)
    res["T3b_bin_cov"] = bool(np.isnan(m[1]))
    # T4：事件：干 6 h 后下 12 mm
    rr = np.zeros(80)
    rr[40:43] = [5.0, 5.0, 2.0]
    ev = find_events(Series(100, rr), 6, 10.0)
    res["T4_events"] = ev == [(140, 12.0)]
    rr[36] = 0.5                                      # 40 不再是 onset（前 6 h 有雨）；36 成 onset，累积 12.5
    res["T4b_dry"] = find_events(Series(100, rr), 6, 10.0) == [(136, 12.5)]
    # T5：delta_pair 配对与已知值
    obs = np.array([34.0] * 6 + [33.0] * 6 + [np.nan] * 6)
    F = np.array([1.0] * 6 + [0.99] * 6 + [1.0] * 6)
    d = delta_pair(obs, F, (0, 6))
    res["T5_delta"] = abs(d[0] + 1.0) < 1e-12 and abs(d[1] + 0.34) < 1e-12 and d[3] == 6 and d[4] == 6
    obs[6:12] = np.nan
    res["T5b_invalid"] = not isnum(delta_pair(obs, F, (0, 6))[0])
    # T6：端到端合成：观测＝ρ×RIM（海面），1 m＝ρ1×RIM；R 应复现 ρ
    S = {}
    H0 = 400000
    n = 24 * 60
    rain_h = np.zeros(n)
    for i, k in enumerate(range(60, n - 40, 48)):
        rain_h[k:k + 3] = [8.0, 6.0, 4.0]
    for cr in CRUISES:
        S[f"{cr}|rain_h"] = Series(H0, rain_h)
        S[f"{cr}|rain_hs"] = Series(2 * H0, np.repeat(rain_h, 2))
        S[f"{cr}|wind_h"] = Series(H0, np.full(n, 5.0))
        S[f"{cr}|lat_h"] = Series(H0, np.full(n, 10.0))
        S[f"{cr}|lon_h"] = Series(H0, np.full(n, -125.0))
    rho = {"snake": 1.0, "snake_met": 1.0, "ssp005": 0.9, "ssp012": 0.8, "ssp023": 0.6, "ssp054": 0.45, "ssp110": 0.3,
           "usps2": 0.3, "usps3": 0.3, "tsg5": 0.3}
    for cr in CRUISES:
        U = np.full(2 * n, 5.0)
        P = np.repeat(rain_h, 2)
        for ln, (_, _, z) in LAYERS.items():
            F = np.full(n, np.nan)
            Fz = p2.rim_factor(P, U, z).reshape(-1, 2).mean(axis=1)    # 小时 24..n−1
            F[24:] = Fz
            S[f"{cr}|{ln}"] = Series(H0, 34.0 * (1.0 + rho[ln] * (F - 1.0)))
    recs, _ = build_ship(S, PRIMARY)
    ok = len(recs) >= 20
    for ln in ("snake", "ssp110", "tsg5"):
        rows = [r for r in recs if lay_ok(r, ln)]
        num = sum(r["L"][ln]["obs"] for r in rows)
        den = sum(r["L"][ln]["rim"] for r in rows)
        ok = ok and abs(num / den - rho[ln]) < 0.02
    res["T6_e2e_ratio_recovers_rho"] = bool(ok)
    main_res, _ = run_cfg(S, PRIMARY)
    res["T6b_e2e_judge_R2"] = main_res["judgement"]["category"].startswith("支持 R2")
    rel = main_res["R_rel5"]
    res["T6c_rel5_gt1"] = rel.get("point") is not None and rel["point"] > 1.0
    lq = layer_q_block(recs, "snake")
    res["T6d_Q"] = lq.get("Q") is not None and abs(lq["Q"] - 1.0 / 0.3) < 0.1
    # T7：R1 情形（全层 ρ=0.3）应判 R1
    for cr in CRUISES:
        for ln, (_, _, z) in LAYERS.items():
            a = S[f"{cr}|{ln}"].a
            S[f"{cr}|{ln}"] = Series(H0, 34.0 + (a - 34.0) * 0.3 / rho[ln])
    main1, _ = run_cfg(S, PRIMARY)
    res["T7_e2e_judge_R1"] = main1["judgement"]["category"].startswith("支持 R1")
    # T8：judge 各分支
    res["T8_category"] = (category(0.1, 0.5).startswith("支持 R1") and category(0.7, 2).startswith("支持 R2")
                          and category(0.5, 0.9) == "不能区分" and category(-1, -0.1).startswith("海面未见"))
    res = {k: (bool(v) if isinstance(v, (bool, np.bool_)) else v) for k, v in res.items()}
    res["all_ok"] = all(v is True for v in res.values())
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--root", default=ROOT_DEFAULT)
    ap.add_argument("--cache", default=CACHE_DEFAULT)
    ap.add_argument("--out", default=os.environ.get("REPRO_OUTPUT_DIR"))
    a = ap.parse_args(argv)
    rc = 0
    try:
        if a.selftest:
            st = selftest()
            print(json.dumps(st, ensure_ascii=False, indent=1))
            if a.out:
                os.makedirs(a.out, exist_ok=True)
                jdump(st, os.path.join(a.out, "p6_selftest.json"))
            rc = 0 if st["all_ok"] else 4
        elif a.full:
            rc = run_full(a)
        else:
            ap.error("need --selftest or --full")
    except Exception:
        log("异常：" + traceback.format_exc())
        rc = 3
    if a.out and a.full:
        with open(os.path.join(a.out, "p6_log.txt"), "a") as f:
            f.write("\n".join(LOG) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
