#!/usr/bin/env python3
"""p7e_piston_spurs1.py — P7e：第二个独立现场检验。PISTON 2018/19 船载 SurfOtter 近表层（最浅盐度 0.35–0.36 m）＋
船载光学雨量计，按 P6 同构的定义、估计量、推断与判读，看雨后近表层观测淡化与 RIM-3 同层预测之比；SPURS-1 WHOI 中心浮标
（0.865 m SBE37＋同址雨量计）只描述 R_1m_obs 是否与 0.30 相容（不判 R1/R2）。

性质：探索性；方法与判读在读盐度数值之前写定。

用法：
  python p7e_piston_spurs1.py --selftest [--out DIR]           合成数据自测（不读原件）
  python p7e_piston_spurs1.py --load-check --piston-root P --spurs1-root S --out DIR
        装载核对：只输出各小时序列的有值小时数／起止与每个事件各层的有值小时计数（不输出任何盐度数值或淡化量）
  python p7e_piston_spurs1.py --full --piston-root P --spurs1-root S --out DIR
        自测 → 读原件（只读）→ 小时化 → 事件 → 模型 → 先写 p7e_power.json → 再写 p7e_summary.json、p7e_piston_events.csv、p7e_spurs1_events.csv
依赖：numpy、netCDF4＋同目录 p6_spurs2_profile.py（事件、小时化、RIM-3 调用、ΔS 配对、簇 bootstrap、判读函数；只 import 不改）
  → p2_rim_test.py（RIM-3 `rim_factor`、风插补、`cluster_boot`）→ p1_events.py、p1b_extend.py。全部路径由命令行给，代码内无路径。
产物（<out>/）：p7e_selftest.json、p7e_power.json（先写）、p7e_summary.json、p7e_piston_events.csv、p7e_spurs1_events.csv、p7e_log.txt；
  --load-check 只写 p7e_loadcheck.json。退出码：0 跑完（无论判读）；3 异常；4 自测不过。

实现选择（X 条）：
  X1 PISTON 两航次：TN2018（R/V Thompson，2018-08-19—10-13）、SR2019（R/V Sally Ride，2019-09-06—09-25）。
     雨：nav-met-sea-1min `prate`（NOAA 光学雨量计，mm/hr）；负值置 0；小时与 30 min 半步按 P6 W3（≥50% 分钟有效）。
  X2 风：`wspd_10N`（10 m 中性、对水）小时均值；该小时缺时用 `wspd_sfc`（传感器高度真风、对水；TN2018 18.00 m、SR2019 14.75 m，
     运行时核对 long_name）的小时均值 ×ln(10/z0)/ln(z_s/z0)（z0＝p2.WIND_Z0_M，P2 K9 同式）补；之后 p6.forcing 原样（≤6 h 插补、下限 0.1）。
  X3 层：SurfOtter `S`（time×depth_CTD，1 min；time 为 Unix 秒）逐文件小时化（P6 W2：≥25% 期望样本）；
     so_top＝该文件 depth_CTD 最浅且 ≤0.5 m 的一列（0.35／0.36 m；TN2018 08-30 文件最浅 1.02 m，不进 so_top）；
     so_1m＝depth_CTD 落在 [0.8, 1.3] m 内离 1.0 m 最近的一列（无则缺）；各文件小时序列拼接（时间不重叠，重叠即报错）；
     每层另存逐小时名义深度序列。船 TSG `ssea_ship`（TN2018 5.00 m、SR2019 3.00 m，long_name 核对）作深层，代码键名沿用 P6 的 "tsg5"
     （以复用 p6.rel5_block；PISTON 的 "tsg5"＝船 TSG，深度按航次）。平台经纬度＝nav-met `lat`／`lon`（SurfOtter 由船拖曳）。
  X4 事件：p6.find_events（P6 W4：雨小时 ≥0.4 mm/h、主干期 6 h、24 h 累积 ≥10 mm）；[t0−30 h, t0+12 h) 须在雨量记录内（p6.forcing）。
  X5 模型深度：so_top／so_1m 取该事件 [t0−6 h, t0+12 h) 内该层有值小时名义深度的中位数；tsg5 取航次固定深度；
     RIM-3＝p6.event_block → p6.model_hourly_F → p2.rim_factor（K2–K4，当前项 t=1 s），一行未改。
  X6 ΔS：p6.delta_pair（P6 W8：雨后窗有值小时均值 − 雨前 6 h 中位数，模型按同一批小时配对），主窗 [0, 6) h。
  X7 PISTON 主估计量 R_top＝Σ ΔS_obs(so_top)/Σ ΔS_RIM(so_top)；簇＝航次×ISO 周；p6.ratio_block（B=10000、seed 20260926、百分位 95%）；
     可评（P6 W11 原样）：n≥15、簇≥4、两航次各 ≥3、bootstrap 中 Σ ΔS_RIM ≥0 的比例 ≤1%；判读 p6.judge（P6 第六节原样：不可评 →
     CI 上端 <0.65 → 下端 >0.65 → 不能区分；稳健标注＝事件级 CI 与逐航次留一同侧；平流标注＝R_rel（so_top − TSG 差分，p6.rel5_block））。
  X8 PISTON 描述：各层 R(z)（so_top、so_1m、tsg5，n≥3）；层梯度 Q_top＝R_top/R_1m 与观测／RIM 层比（so_top 与 so_1m 同一事件都有效的子集，
     同一簇抽样）；雨中／雨停 p6.rain_split（P6 6.3 原样，sfc＝so_top）。
  X9 PISTON 敏感性（只报）：S1 干期 24 h；S2 干期 12 h；S3 阈值 5 mm；S4 窗 [0,3)；S5 窗 [0,12)；S6 位移 ≤30 km；S9 当前项 t=1800 s；
     S10 雨前、雨后窗各 ≥4 有值小时。（P6 的 S7、S8 无对应，不设。）
  X10 SPURS-1：`SPURS_WHOI_1_D_M.nc`（1 min）的 RAIN（mm/hour，负值置 0）、UWND/VWND（高度取文件 HEIGHT_WND，按 X2 同一对数律换算 10 m）、
     PSAL（SBE37，深度取文件 DEPTH＝0.865 m）→ 层 s086；`SPURS_WHOI_1_D_TS.nc` 的 PSAL 在 DEPTH＝2.1 m、5.2 m 两列 → 层 s21、tsg5（5.2 m）。
     事件主定义＝P1 原样（干期 24 h、10 mm、窗 [0,6)，与 P6 浮标一致性检查同口径）；簇＝日历月。
  X11 SPURS-1 读法 reading_s1()：n<8 或簇<4 或 CI 不可算 →「不可评」；否则按 95% 簇 CI 与 0.30 的关系三分（含／整体小于／整体大于），
     并标注是否 CI 上端 <0.65（与 P2 同侧）、是否含 1。信息项：干期 6 h、窗 [0,12)、阈值 5 mm；s21、tsg5 层 R；R_rel（s086 − 5.2 m）。
  X12 功效（p7e_power.json，先于任何观测 ΔS 合计写盘）：PISTON（so_top）与 SPURS-1（s086）各自模型侧 x 的 n、Σx、均值、RMS、簇数；
     无雨对照窗（[c−12, c+6) 全有效且无雨、c 每 6 h 一个）上观测 ΔS 的标准差与推得的 SE(R)；σ＝0.1／0.2／0.4 三档 SE。只作信息。
Change Log：
  2026-09-27 初版（P7e）。
"""

import argparse
import csv
import glob
import hashlib
import json
import math
import os
import re
import sys
import time
import traceback
from datetime import datetime, timezone

import numpy as np

import p6_spurs2_profile as p6
import p2_rim_test as p2

VERSION = "p7e-2026-09-27a"
NAN = float("nan")
THRESH = p6.THRESH
R_REF = 0.30
PRE_H = p6.PRE_H
POST_MAX_H = p6.POST_MAX_H
RAIN_HOUR_MM = p6.RAIN_HOUR_MM
TOP_MAX_Z = 0.5
ONE_M_RANGE = (0.8, 1.3)
S1_MIN_N, S1_MIN_K = 8, 4

PISTON = {
    "TN2018": {"dir": "PISTON-ONR-NOAA_RVThompson_2018",
               "met": "PISTON-nav-met-sea-1min_RV-Thompson_20180819_R1_thru_20181012.nc",
               "otter": "piston-SurfOtter_RV-Thompson_*.nc", "wind_z": 18.00, "tsg_z": 5.0},
    "SR2019": {"dir": "PISTON-ONR-NOAA_RVSallyRide_2019",
               "met": "PISTON-nav-met-sea-1min_RV-Sally-Ride_20190906_R1_thru_20190925.nc",
               "otter": "piston-SurfOtter_RV-Sally-Ride_*.nc", "wind_z": 14.75, "tsg_z": 3.0},
}
SPURS1 = {"met": os.path.join("SPURS1_MOORING_WHOI", "SPURS_WHOI_1_D_M.nc"),
          "ts": os.path.join("SPURS1_MOORING_WHOI", "SPURS_WHOI_1_D_TS.nc")}
P_LAYERS = ["so_top", "so_1m", "tsg5"]
S1_LAYERS = ["s086", "s21", "tsg5"]

PRIMARY = {"name": "primary", "dry_h": 6, "thr": 10.0, "win": (0, 6), "sfc": "so_top", "z_sfc": 0.0,
           "t_cur": p2.T_CURRENT_S, "disp_max": None, "min_hours": 1}
SENS = [
    dict(PRIMARY, name="S1_dry24", dry_h=24),
    dict(PRIMARY, name="S2_dry12", dry_h=12),
    dict(PRIMARY, name="S3_thr5", thr=5.0),
    dict(PRIMARY, name="S4_win0_3", win=(0, 3)),
    dict(PRIMARY, name="S5_win0_12", win=(0, 12)),
    dict(PRIMARY, name="S6_disp30", disp_max=30.0),
    dict(PRIMARY, name="S9_tcur1800", t_cur=1800.0),
    dict(PRIMARY, name="S10_min4h", min_hours=4),
]
S1_PRIMARY = dict(PRIMARY, name="s1_primary_P1def", dry_h=24, sfc="s086")
S1_INFO = [dict(S1_PRIMARY, name="s1_dry6", dry_h=6),
           dict(S1_PRIMARY, name="s1_win0_12", win=(0, 12)),
           dict(S1_PRIMARY, name="s1_thr5", thr=5.0)]
LOG = []


def log(msg):
    line = f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} {msg}"
    LOG.append(line)
    print(line, flush=True)


rnd = p6.rnd
isnum = p6.isnum


def jdump(obj, path):
    p6.jdump(obj, path)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ======================================================================== 小工具（X2、X3）
def log_factor(z_sensor):
    return math.log(10.0 / p2.WIND_Z0_M) / math.log(z_sensor / p2.WIND_Z0_M)


def to_grid(k0, a, lo, hi):
    out = np.full(hi - lo, np.nan)
    if k0 is None or len(a) == 0:
        return out
    s, e = max(lo, k0), min(hi, k0 + len(a))
    if e > s:
        out[s - lo:e - lo] = a[s - k0:e - k0]
    return out


def wind_compose(k10, w10, ks, ws, z_sensor):
    """X2：10N 小时值优先，缺时用传感器高度真风×对数律。返回 (k0, 数组, 补的小时数)。"""
    ks_ = [k for k, a in ((k10, w10), (ks, ws)) if k is not None and len(a)]
    if not ks_:
        return None, np.array([]), 0
    lo = min(ks_)
    hi = max(k + len(a) for k, a in ((k10, w10), (ks, ws)) if k is not None and len(a))
    a10 = to_grid(k10, w10, lo, hi)
    as_ = to_grid(ks, ws, lo, hi) * log_factor(z_sensor)
    fill = ~np.isfinite(a10) & np.isfinite(as_)
    return lo, np.where(np.isfinite(a10), a10, as_), int(fill.sum())


def pick_cols(dz):
    """X3：返回 (最浅且 ≤0.5 m 的列下标或 None, [0.8,1.3] 内离 1.0 m 最近的列下标或 None)。"""
    dz = np.asarray(dz, float)
    top = int(np.argmin(dz)) if np.nanmin(dz) <= TOP_MAX_Z else None
    cand = [j for j in range(len(dz)) if ONE_M_RANGE[0] <= dz[j] <= ONE_M_RANGE[1]]
    one = min(cand, key=lambda j: abs(dz[j] - 1.0)) if cand else None
    return top, one


def merge_series(parts):
    """parts：[(k0, 数组)]，时间不重叠；返回 p6.Series。"""
    parts = [(k0, np.asarray(a, float)) for k0, a in parts if k0 is not None and len(a)]
    if not parts:
        return p6.Series(0, np.array([]))
    lo = min(k0 for k0, _ in parts)
    hi = max(k0 + len(a) for k0, a in parts)
    out = np.full(hi - lo, np.nan)
    for k0, a in parts:
        seg = out[k0 - lo:k0 - lo + len(a)]
        m = np.isfinite(a)
        if np.any(np.isfinite(seg) & m):
            raise RuntimeError("SurfOtter 文件小时重叠")
        seg[m] = a[m]
    return p6.Series(lo, out)


def nc_time_any(ds, name):
    v = ds.variables[name]
    u = str(getattr(v, "units", ""))
    if " since " in u:
        return p6.nc_time(ds, name)
    if u.strip() == "s" and "unix time" in str(getattr(v, "long_name", "")).lower():
        return p6.nc_var(ds, name)
    raise RuntimeError(f"未知时间单位 {name}: {u!r}")


def check_ln(ds, var, needle):
    ln = str(getattr(ds.variables[var], "long_name", ""))
    if needle not in ln:
        raise RuntimeError(f"{var} long_name 不含 {needle!r}：{ln!r}")


# ======================================================================== 读原件（只读）
def load_piston(root, S, info):
    import netCDF4
    for cr, c in PISTON.items():
        d = os.path.join(root, c["dir"])
        ds = netCDF4.Dataset(os.path.join(d, c["met"]))
        tm = p6.nc_time(ds, "time")
        r = np.maximum(p6.nc_var(ds, "prate"), 0.0)
        k0, hr = p6.bin_mean(tm, r, 3600.0, p6.MIN_FRAC_RAIN)
        S[f"{cr}|rain_h"] = p6.Series(k0, hr)
        k0, hs = p6.bin_mean(tm, r, 1800.0, p6.MIN_FRAC_RAIN)
        S[f"{cr}|rain_hs"] = p6.Series(k0, hs)
        check_ln(ds, "wspd_sfc", f"{c['wind_z']:.2f} m")
        k10, w10 = p6.bin_mean(tm, p6.nc_var(ds, "wspd_10N"), 3600.0, p6.MIN_FRAC_RAIN)
        ks, ws = p6.bin_mean(tm, p6.nc_var(ds, "wspd_sfc"), 3600.0, p6.MIN_FRAC_RAIN)
        kw, w, nfill = wind_compose(k10, w10, ks, ws, c["wind_z"])
        S[f"{cr}|wind_h"] = p6.Series(kw, w)
        for nm in ("lat", "lon"):
            k0, a = p6.bin_mean(tm, p6.nc_var(ds, nm), 3600.0, p6.MIN_FRAC_RAIN)
            S[f"{cr}|{nm}_h"] = p6.Series(k0, a)
        check_ln(ds, "ssea_ship", f"{c['tsg_z']:.2f} depth")
        k0, a = p6.bin_mean(tm, p6.nc_var(ds, "ssea_ship"), 3600.0, p6.MIN_FRAC_SAL)
        S[f"{cr}|tsg5"] = p6.Series(k0, a)
        ds.close()
        parts = {"so_top": [], "so_1m": [], "so_top_z": [], "so_1m_z": []}
        files = sorted(glob.glob(os.path.join(d, c["otter"])))
        finfo = []
        for fp in files:
            ds = netCDF4.Dataset(fp)
            t = nc_time_any(ds, "time")
            dz = p6.nc_var(ds, "depth_CTD")
            if tuple(ds.variables["S"].dimensions) != ("time", "depth_CTD"):
                raise RuntimeError(f"{fp} S 维序 {ds.variables['S'].dimensions}")
            sal = p6.nc_var(ds, "S")
            top, one = pick_cols(dz)
            fi = {"file": os.path.basename(fp), "depth_CTD": [float(x) for x in dz], "top_col": top, "one_col": one}
            for ln, j in (("so_top", top), ("so_1m", one)):
                if j is None:
                    continue
                k0, m = p6.bin_mean(t, sal[:, j], 3600.0, p6.MIN_FRAC_SAL)
                parts[ln].append((k0, m))
                parts[ln + "_z"].append((k0, np.where(np.isfinite(m), float(dz[j]), np.nan)))
            ds.close()
            finfo.append(fi)
        for k, v in parts.items():
            S[f"{cr}|{k}"] = merge_series(v)
        info[cr] = {"wind_hours_filled_from_sensor_height": nfill, "surfotter_files": finfo,
                    "wind_sensor_z": c["wind_z"], "tsg_z": c["tsg_z"]}


def load_spurs1(root, S, info):
    import netCDF4
    ds = netCDF4.Dataset(os.path.join(root, SPURS1["met"]))
    tm = p6.nc_time(ds, "TIME")
    r = np.maximum(p6.nc_var(ds, "RAIN"), 0.0)
    k0, hr = p6.bin_mean(tm, r, 3600.0, p6.MIN_FRAC_RAIN)
    S["S1M|rain_h"] = p6.Series(k0, hr)
    k0, hs = p6.bin_mean(tm, r, 1800.0, p6.MIN_FRAC_RAIN)
    S["S1M|rain_hs"] = p6.Series(k0, hs)
    zw = float(np.asarray(ds.variables["HEIGHT_WND"][:]).item())
    spd = np.hypot(p6.nc_var(ds, "UWND"), p6.nc_var(ds, "VWND"))
    k0, w = p6.bin_mean(tm, spd, 3600.0, p6.MIN_FRAC_RAIN)
    S["S1M|wind_h"] = p6.Series(k0, w * log_factor(zw))
    z086 = float(np.asarray(ds.variables["DEPTH"][:]).item())
    k0, a = p6.bin_mean(tm, p6.nc_var(ds, "PSAL"), 3600.0, p6.MIN_FRAC_SAL)
    S["S1M|s086"] = p6.Series(k0, a)
    ds.close()
    ds = netCDF4.Dataset(os.path.join(root, SPURS1["ts"]))
    tt = p6.nc_time(ds, "TIME")
    dep = p6.nc_var(ds, "DEPTH")
    if tuple(ds.variables["PSAL"].dimensions) != ("TIME", "DEPTH"):
        raise RuntimeError(f"TS PSAL 维序 {ds.variables['PSAL'].dimensions}")
    ps = p6.nc_var(ds, "PSAL")
    idx = {}
    for key, zz in (("s21", 2.1), ("tsg5", 5.2)):
        j = int(np.argmin(np.abs(dep - zz)))
        if abs(dep[j] - zz) > 1e-6:
            raise RuntimeError(f"TS 无 {zz} m 层（最近 {dep[j]}）")
        k0, a = p6.bin_mean(tt, ps[:, j], 3600.0, p6.MIN_FRAC_SAL)
        S[f"S1M|{key}"] = p6.Series(k0, a)
        idx[key] = j
    ds.close()
    info["S1M"] = {"wind_sensor_z": zw, "s086_depth": z086, "ts_depth_index": idx}
    return {"s086": z086, "s21": 2.1, "tsg5": 5.2}


# ======================================================================== 事件集（X4、X5、X10）
def build_piston(S, cfg):
    recs, drop = [], {}
    uid = 0
    for cr in PISTON:
        for H, acc in p6.find_events(S[f"{cr}|rain_h"], cfg["dry_h"], cfg["thr"]):
            layers = {}
            for ln in ("so_top", "so_1m"):
                zz = S[f"{cr}|{ln}_z"].get(H - PRE_H, H + POST_MAX_H)
                zz = zz[np.isfinite(zz)]
                if len(zz):
                    layers[ln] = (f"{cr}|{ln}", float(np.median(zz)))
            layers["tsg5"] = (f"{cr}|tsg5", PISTON[cr]["tsg_z"])
            b = p6.event_block(S, cr, H, layers, cfg)
            if b["diag"].get("reason"):
                drop[b["diag"]["reason"]] = drop.get(b["diag"]["reason"], 0) + 1
                continue
            la1, lo1 = np.nanmedian(S[f"{cr}|lat_h"].get(H - 6, H)), np.nanmedian(S[f"{cr}|lon_h"].get(H - 6, H))
            la2, lo2 = np.nanmedian(S[f"{cr}|lat_h"].get(H, H + 6)), np.nanmedian(S[f"{cr}|lon_h"].get(H, H + 6))
            disp = p6.haversine(la1, lo1, la2, lo2) if all(np.isfinite([la1, lo1, la2, lo2])) else NAN
            if cfg["disp_max"] is not None and not (np.isfinite(disp) and disp <= cfg["disp_max"]):
                drop["disp"] = drop.get("disp", 0) + 1
                continue
            _, iw, _ = datetime.fromtimestamp(H * 3600, timezone.utc).isocalendar()
            uid += 1
            recs.append({"uid": uid, "cruise": cr, "H": H, "onset": p6.iso(H), "acc24": acc, "cluster": f"{cr}-W{iw:02d}",
                         "disp_km": disp, "diag": b["diag"], "L": b["layers"],
                         "z": {ln: z for ln, (_, z) in layers.items()},
                         "rain12": S[f"{cr}|rain_h"].get(H, H + POST_MAX_H)})
    return recs, drop


def build_spurs1(S, cfg, zmap):
    recs, drop = [], {}
    for i, (H, acc) in enumerate(p6.find_events(S["S1M|rain_h"], cfg["dry_h"], cfg["thr"])):
        layers = {ln: (f"S1M|{ln}", zmap[ln]) for ln in S1_LAYERS}
        b = p6.event_block(S, "S1M", H, layers, cfg)
        if b["diag"].get("reason"):
            drop[b["diag"]["reason"]] = drop.get(b["diag"]["reason"], 0) + 1
            continue
        d = datetime.fromtimestamp(H * 3600, timezone.utc)
        recs.append({"uid": i + 1, "cruise": "S1M", "H": H, "onset": p6.iso(H), "acc24": acc,
                     "cluster": f"{d.year}-{d.month:02d}", "disp_km": 0.0, "diag": b["diag"], "L": b["layers"],
                     "z": dict(zmap), "rain12": S["S1M|rain_h"].get(H, H + POST_MAX_H)})
    return recs, drop


# ======================================================================== 统计与判读（X7、X8、X11）
def evaluable(blk, rows, cruises):
    per = {}
    for r in rows:
        per[r["cruise"]] = per.get(r["cruise"], 0) + 1
    conds = {"n>=15": blk["n"] >= 15, "clusters>=4": (blk.get("clusters") or 0) >= 4,
             "each_cruise>=3": len(per) == len(cruises) and min(per.values()) >= 3,
             "frac_den_nonneg<=0.01": blk.get("frac_boot_den_nonneg") is not None and blk["frac_boot_den_nonneg"] <= p6.DEN_NONNEG_MAX}
    return all(conds.values()), conds, per


def layer_block(recs, ln):
    rows = [r for r in recs if p6.lay_ok(r, ln)]
    if len(rows) < 3:
        return {"n": len(rows), "note": "事件 <3，不算"}
    b = p6.strip(p6.ratio_block(rows, p6.f_obs(ln), p6.f_rim(ln), p6.clus))
    zs = [r["z"].get(ln) for r in rows]
    b["z_model_median"] = rnd(float(np.median(zs)), 3)
    return b


def q_block(recs, a="so_top", b="so_1m"):
    rows = [r for r in recs if p6.lay_ok(r, a) and p6.lay_ok(r, b)]
    out = {"n": len(rows)}
    if len(rows) < 3:
        return dict(out, note="事件 <3，不算")
    so, sm = sum(r["L"][a]["obs"] for r in rows), sum(r["L"][a]["rim"] for r in rows)
    oo, om = sum(r["L"][b]["obs"] for r in rows), sum(r["L"][b]["rim"] for r in rows)
    bc, K = p6.boot(rows, {"Ra": (p6.f_obs(a), p6.f_rim(a)), "Rb": (p6.f_obs(b), p6.f_rim(b)),
                           "LRo": (p6.f_obs(a), p6.f_obs(b)), "LRm": (p6.f_rim(a), p6.f_rim(b))}, p6.clus)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = bc["Ra"] / bc["Rb"]
    out.update({"clusters": K, f"R_{a}": rnd(so / sm), f"R_{a}_ci95": p6.pct(bc["Ra"]),
                f"R_{b}": rnd(oo / om), f"R_{b}_ci95": p6.pct(bc["Rb"]),
                "Q": rnd((so / sm) / (oo / om)) if oo and om else None, "Q_ci95": p6.pct(q),
                "layer_ratio_obs": rnd(so / oo) if oo else None, "layer_ratio_obs_ci95": p6.pct(bc["LRo"]),
                "layer_ratio_rim": rnd(sm / om) if om else None, "layer_ratio_rim_ci95": p6.pct(bc["LRm"])})
    return out


def reading_s1(blk):
    lo, hi = blk.get("ci95_cluster") or [None, None]
    if blk.get("n", 0) < S1_MIN_N or (blk.get("clusters") or 0) < S1_MIN_K or lo is None or hi is None:
        return {"reading": "不可评（n<8 或簇<4 或 CI 不可算）", "ci_upper_lt_0p65": None, "ci_contains_1": None}
    if lo <= R_REF <= hi:
        rd = "与 0.30 相容（95% 簇 CI 含 0.30）"
    elif hi < R_REF:
        rd = "小于 0.30（CI 上端 <0.30）"
    else:
        rd = "大于 0.30（CI 下端 >0.30）"
    return {"reading": rd, "ci_upper_lt_0p65": bool(hi < THRESH), "ci_contains_1": bool(lo <= 1.0 <= hi)}


def power_block(recs, S, cfg, plats):
    """X12：只用模型 x 与无雨对照窗观测噪声。"""
    ln = cfg["sfc"]
    xs = np.array([r["L"][ln]["rim"] for r in recs if p6.lay_ok(r, ln)], float)
    n = len(xs)
    K = len({r["cluster"] for r in recs if p6.lay_ok(r, ln)})
    ctrl = []
    for pl in plats:
        rh = S[f"{pl}|rain_h"]
        a, b = rh.span()
        for c in range(a + 12, b - 6, 6):
            seg = rh.get(c - 12, c + 6)
            if not (np.isfinite(seg).all() and (seg < RAIN_HOUR_MM).all()):
                continue
            obs18 = S[f"{pl}|{ln}"].get(c - PRE_H, c + POST_MAX_H)
            dso = p6.delta_pair(obs18, np.zeros(18), cfg["win"], 1)[0]
            if isnum(dso):
                ctrl.append(dso)
    ctrl = np.array(ctrl, float)
    sx = float(xs.sum()) if n else NAN
    out = {"layer": ln, "n": n, "clusters": K, "sum_x": rnd(sx), "mean_x": rnd(xs.mean()) if n else None,
           "rms_x": rnd(np.sqrt((xs ** 2).mean())) if n else None,
           "n_ctrl_windows": int(len(ctrl)), "sd_ctrl_dS": rnd(ctrl.std(ddof=1)) if len(ctrl) > 2 else None}
    if n and abs(sx) > 0:
        for sg in (0.1, 0.2, 0.4):
            out[f"se_R_if_sigma_{sg}"] = rnd(sg * math.sqrt(n) / abs(sx))
        if out["sd_ctrl_dS"]:
            se = out["sd_ctrl_dS"] * math.sqrt(n) / abs(sx)
            out["se_R_ctrl"] = rnd(se)
            out["ci_halfwidth_ctrl"] = rnd(1.96 * se)
    return out


def run_piston(S, cfg):
    recs, drop = build_piston(S, cfg)
    ln = cfg["sfc"]
    rows = [r for r in recs if p6.lay_ok(r, ln)]
    blk = p6.ratio_block(rows, p6.f_obs(ln), p6.f_rim(ln), p6.clus)
    ok, conds, per = evaluable(blk, rows, PISTON) if rows else (False, {}, {})
    loo = p6.loo_cruise(rows, p6.f_obs(ln), p6.f_rim(ln))
    rel = p6.rel5_block(recs, cfg)
    j = p6.judge(blk, ok, loo, rel.get("category"))
    res = {"cfg": dict(cfg), "n_events_all": len(recs), "dropped": drop,
           "R_top": dict(p6.strip(blk), evaluable=ok, evaluable_conditions=conds, per_cruise=per, loo_cruise=loo),
           "R_rel_tsg": p6.strip(rel), "judgement": j}
    return res, recs


def run_spurs1(S, cfg, zmap):
    recs, drop = build_spurs1(S, cfg, zmap)
    rows = [r for r in recs if p6.lay_ok(r, "s086")]
    blk = p6.strip(p6.ratio_block(rows, p6.f_obs("s086"), p6.f_rim("s086"), p6.clus)) if rows else {"n": 0}
    blk.update(reading_s1(blk))
    blk.update({"events_all": len(recs), "dropped": drop})
    return blk, recs


def events_csv(recs, layers, path):
    cols = ["platform", "onset", "acc24", "cluster", "disp_km"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols + [f"{ln}_{k}" for ln in layers for k in ("z", "obs", "rim", "s0", "npre", "npost")] + ["rain_hours_0_12"])
        for r in recs:
            row = [r["cruise"], r["onset"], rnd(r["acc24"], 2), r["cluster"], rnd(r["disp_km"], 1)]
            for ln in layers:
                L = r["L"].get(ln, {})
                row += [rnd(r["z"].get(ln), 3), rnd(L.get("obs")), rnd(L.get("rim")), rnd(L.get("s0"), 4), L.get("n_pre"), L.get("n_post")]
            row.append(int(np.nansum(r["rain12"] >= RAIN_HOUR_MM)))
            w.writerow(row)


# ======================================================================== 主流程
def code_shas():
    here = os.path.dirname(os.path.abspath(__file__))
    return {fn: sha256_file(os.path.join(here, fn)) for fn in
            ("p7e_piston_spurs1.py", "p6_spurs2_profile.py", "p2_rim_test.py", "p1_events.py", "p1b_extend.py")}


def load_all(a):
    S, info = {}, {}
    load_piston(a.piston_root, S, info)
    zmap = load_spurs1(a.spurs1_root, S, info)
    return S, info, zmap


def run_load_check(a):
    """装载核对：只输出有值小时计数与起止、事件与各层有值小时数；不输出盐度数值、不算淡化量。"""
    os.makedirs(a.out, exist_ok=True)
    S, info, zmap = load_all(a)
    ser = {}
    for k, s in sorted(S.items()):
        fin = np.isfinite(s.a)
        ser[k] = {"k0": p6.iso(s.k0) if len(s.a) else None, "hours": int(len(s.a)), "finite_hours": int(fin.sum())}
    ev = {}
    for tag, recs in (("piston_primary", build_piston_counts(S, PRIMARY)), ("spurs1_P1def", build_s1_counts(S, S1_PRIMARY))):
        ev[tag] = recs
    jdump({"version": VERSION, "note": "只含可用性计数，不含盐度数值或淡化量", "info": info, "series": ser, "events": ev,
           "code_sha256": code_shas()}, os.path.join(a.out, "p7e_loadcheck.json"))
    log("load-check done")
    return 0


def _avail(S, key, H):
    s = S[key].get(H - PRE_H, H + POST_MAX_H)
    f = np.isfinite(s)
    return [int(f[:PRE_H].sum()), int(f[PRE_H:PRE_H + 6].sum()), int(f[PRE_H:].sum())]


def build_piston_counts(S, cfg):
    out = []
    for cr in PISTON:
        for H, acc in p6.find_events(S[f"{cr}|rain_h"], cfg["dry_h"], cfg["thr"]):
            P, U, diag = p6.forcing(S, cr, H)
            out.append({"cruise": cr, "onset": p6.iso(H), "acc24": rnd(acc, 2), "forcing": diag.get("reason") or "ok",
                        "hours_pre6_post6_post12": {ln: _avail(S, f"{cr}|{ln}", H) for ln in P_LAYERS}})
    return out


def build_s1_counts(S, cfg):
    out = []
    for H, acc in p6.find_events(S["S1M|rain_h"], cfg["dry_h"], cfg["thr"]):
        P, U, diag = p6.forcing(S, "S1M", H)
        out.append({"onset": p6.iso(H), "acc24": rnd(acc, 2), "forcing": diag.get("reason") or "ok",
                    "hours_pre6_post6_post12": {ln: _avail(S, f"S1M|{ln}", H) for ln in S1_LAYERS}})
    return out


def run_full(a):
    t_start = time.time()
    out = a.out
    os.makedirs(out, exist_ok=True)
    st = selftest()
    jdump(st, os.path.join(out, "p7e_selftest.json"))
    if not st["all_ok"]:
        log("自测不过，退出 4")
        return 4
    code_sha = code_shas()
    log(f"code sha {json.dumps(code_sha)}")
    S, info, zmap = load_all(a)
    log(f"series {len(S)} info {json.dumps(info, ensure_ascii=False)[:600]}")
    # ---- 主事件集与功效（先写）
    precs, pdrop = build_piston(S, PRIMARY)
    srecs, sdrop = build_spurs1(S, S1_PRIMARY, zmap)
    log(f"piston primary events {len(precs)} dropped {pdrop}; spurs1 P1-def events {len(srecs)} dropped {sdrop}")
    pw = {"version": VERSION, "code_sha256": code_sha, "written": datetime.now(timezone.utc).isoformat(),
          "piston_so_top": power_block(precs, S, PRIMARY, list(PISTON)),
          "spurs1_s086": power_block(srecs, S, S1_PRIMARY, ["S1M"]),
          "note": "只含模型侧 x 与无雨对照窗 ΔS 的离散度；只作信息，不改判读"}
    jdump(pw, os.path.join(out, "p7e_power.json"))
    log("power written")
    time.sleep(1.0)
    # ---- PISTON 主结果与判读
    main_res, precs = run_piston(S, PRIMARY)
    prof = {ln: layer_block(precs, ln) for ln in P_LAYERS}
    q = q_block(precs)
    rain = p6.rain_split(precs, PRIMARY)
    sens = {}
    for c in SENS:
        res, _ = run_piston(S, c)
        j = res["judgement"]
        sens[c["name"]] = {"n": res["R_top"]["n"], "clusters": res["R_top"].get("clusters"), "point": res["R_top"].get("point"),
                           "ci95_cluster": res["R_top"].get("ci95_cluster"), "evaluable": res["R_top"]["evaluable"],
                           "category": j["category"], "same_as_primary": j["category"] == main_res["judgement"]["category"],
                           "R_rel_tsg": {"n": res["R_rel_tsg"].get("n"), "point": res["R_rel_tsg"].get("point"),
                                         "ci95_cluster": res["R_rel_tsg"].get("ci95_cluster"),
                                         "category": res["R_rel_tsg"].get("category")}}
    # ---- SPURS-1 描述
    s1_main, srecs = run_spurs1(S, S1_PRIMARY, zmap)
    s1_info = {}
    for c in S1_INFO:
        b, _ = run_spurs1(S, c, zmap)
        s1_info[c["name"]] = b
    s1_layers = {ln: layer_block(srecs, ln) for ln in S1_LAYERS}
    s1_rel = p6.rel5_block(srecs, dict(S1_PRIMARY, sfc="s086"))
    summ = {"version": VERSION, "code_sha256": code_sha, "data_info": info,
            "piston": {"primary": main_res, "verdict": main_res["judgement"], "profile_R_by_layer": prof,
                       "Q_top_vs_1m": q, "rain_vs_dry": rain, "sensitivity": sens,
                       "disp_km_primary": {"median": rnd(np.nanmedian([r["disp_km"] for r in precs])) if precs else None,
                                           "max": rnd(np.nanmax([r["disp_km"] for r in precs])) if precs else None}},
            "spurs1": {"R_1m_obs_s086": s1_main, "info_variants": s1_info, "layers_P1def": s1_layers,
                       "R_rel_s086_minus_5p2m": p6.strip(s1_rel)},
            "elapsed_s": round(time.time() - t_start, 1)}
    jdump(summ, os.path.join(out, "p7e_summary.json"))
    events_csv(precs, P_LAYERS, os.path.join(out, "p7e_piston_events.csv"))
    events_csv(srecs, S1_LAYERS, os.path.join(out, "p7e_spurs1_events.csv"))
    log(f"piston verdict {json.dumps(main_res['judgement'], ensure_ascii=False)}")
    log(f"spurs1 reading {s1_main.get('reading')}")
    log("done rc=0")
    return 0


# ======================================================================== 自测
def _synth(S, plat, n, H0, rain_h, rho, zs):
    S[f"{plat}|rain_h"] = p6.Series(H0, rain_h)
    S[f"{plat}|rain_hs"] = p6.Series(2 * H0, np.repeat(rain_h, 2))
    S[f"{plat}|wind_h"] = p6.Series(H0, np.full(n, 5.0))
    S[f"{plat}|lat_h"] = p6.Series(H0, np.full(n, 10.0))
    S[f"{plat}|lon_h"] = p6.Series(H0, np.full(n, 130.0))
    U = np.full(2 * n, 5.0)
    P = np.repeat(rain_h, 2)
    for ln, z in zs.items():
        F = np.full(n, np.nan)
        F[24:] = p2.rim_factor(P, U, z).reshape(-1, 2).mean(axis=1)
        S[f"{plat}|{ln}"] = p6.Series(H0, 34.0 * (1.0 + rho[ln] * (F - 1.0)))
        S[f"{plat}|{ln}_z"] = p6.Series(H0, np.full(n, z))


def selftest():
    res = {}
    st6 = p6.selftest()
    res["T0_p6_selftest_all_ok"] = bool(st6["all_ok"])
    rng = np.random.default_rng(11)
    P = np.where(rng.random(84) < 0.3, rng.gamma(1.5, 6.0, 84), 0.0)
    U = np.repeat(rng.uniform(1, 11, 42), 2)
    mx = 0.0
    for z in (0.35, 0.36, 0.865, 1.16, 2.1, 3.0, 5.0, 5.2):
        mx = max(mx, float(np.nanmax(np.abs(p6.model_hourly_F(P, U, z) - p2.rim_factor(P, U, z).reshape(-1, 2).mean(axis=1)))))
    res["T1_model_equals_p2_rim_factor_all_depths"] = mx == 0.0
    ref, _ = p2.model_hourly(P, U, 34.0, {0.5: 0.5, 1.0: 1.0, 5.0: 5.0})
    res["T1b_model_equals_p2_model_hourly"] = max(float(np.nanmax(np.abs(p6.model_hourly_F(P, U, z) - ref[z]))) for z in (0.5, 1.0, 5.0)) == 0.0
    # T2 风合成
    k0, w, nf = wind_compose(10, np.array([5.0, np.nan, 7.0]), 10, np.array([6.0, 6.0, np.nan]), 14.75)
    res["T2_wind_compose"] = (k0 == 10 and w[0] == 5.0 and abs(w[1] - 6.0 * log_factor(14.75)) < 1e-12 and w[2] == 7.0 and nf == 1)
    # T3 列选择
    res["T3_pick_cols"] = (pick_cols([0.35, 1.16, 2.6, 3.79, 5.3]) == (0, 1) and pick_cols([1.02, 2.45, 3.88]) == (None, 0)
                           and pick_cols([0.36, 1.87, 2.83]) == (0, None) and pick_cols([0.36, 0.81, 1.16, 1.87]) == (0, 2))
    # T4 拼接
    s = merge_series([(5, [1.0, np.nan]), (7, [2.0])])
    ok4 = s.k0 == 5 and list(np.nan_to_num(s.a, nan=-9)) == [1.0, -9, 2.0]
    try:
        merge_series([(5, [1.0, 1.0]), (6, [2.0])])
        ok4 = False
    except RuntimeError:
        pass
    res["T4_merge"] = bool(ok4)
    # T5 PISTON 端到端（R1 型：各层 ρ＝0.3；R2 型：ρ＝1）
    n = 24 * 60
    H0 = 400000
    rain_h = np.zeros(n)
    for k in range(60, n - 40, 48):
        rain_h[k:k + 3] = [8.0, 6.0, 4.0]
    zs = {"so_top": 0.35, "so_1m": 1.16, "tsg5": 3.0}
    for rho_v, tag, want in ((0.3, "R1", "支持 R1"), (1.0, "R2", "支持 R2")):
        S = {}
        for cr in PISTON:
            _synth(S, cr, n, H0, rain_h, {k: rho_v for k in zs}, zs)
        r, recs = run_piston(S, PRIMARY)
        ok = len(recs) >= 20 and abs(r["R_top"]["point"] - rho_v) < 0.02 and r["judgement"]["category"].startswith(want)
        lb = layer_block(recs, "so_1m")
        ok = ok and abs(lb["point"] - rho_v) < 0.02 and abs(q_block(recs)["Q"] - 1.0) < 0.05
        res[f"T5_piston_e2e_{tag}"] = bool(ok)
    # T6 SPURS-1 端到端（ρ＝0.3 → 与 0.30 相容）
    S = {}
    zm = {"s086": 0.865, "s21": 2.1, "tsg5": 5.2}
    n1 = 24 * 200                                     # 约 6.6 个月，月簇 ≥4
    rain1 = np.zeros(n1)
    for k in range(60, n1 - 40, 48):
        rain1[k:k + 3] = [8.0, 6.0, 4.0]
    _synth(S, "S1M", n1, H0, rain1, {k: 0.3 for k in zm}, zm)
    blk, recs = run_spurs1(S, dict(S1_PRIMARY, dry_h=24), zm)
    res["T6_spurs1_e2e"] = bool(abs(blk["point"] - 0.3) < 0.02 and blk["reading"].startswith("与 0.30 相容"))
    # T7 读法分支
    res["T7_reading_s1"] = (reading_s1({"n": 10, "clusters": 5, "ci95_cluster": [0.1, 0.5]})["reading"].startswith("与 0.30")
                            and reading_s1({"n": 10, "clusters": 5, "ci95_cluster": [0.05, 0.2]})["reading"].startswith("小于")
                            and reading_s1({"n": 10, "clusters": 5, "ci95_cluster": [0.4, 1.2]})["reading"].startswith("大于")
                            and reading_s1({"n": 5, "clusters": 5, "ci95_cluster": [0.1, 0.5]})["reading"].startswith("不可评")
                            and reading_s1({"n": 10, "clusters": 5, "ci95_cluster": [0.4, 1.2]})["ci_contains_1"] is True)
    res = {k: (bool(v) if isinstance(v, (bool, np.bool_)) else v) for k, v in res.items()}
    res["all_ok"] = all(v is True for v in res.values())
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--load-check", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--piston-root")
    ap.add_argument("--spurs1-root")
    ap.add_argument("--out", default=os.environ.get("REPRO_OUTPUT_DIR"))
    a = ap.parse_args(argv)
    rc = 0
    try:
        if a.selftest:
            st = selftest()
            print(json.dumps(st, ensure_ascii=False, indent=1))
            if a.out:
                os.makedirs(a.out, exist_ok=True)
                jdump(st, os.path.join(a.out, "p7e_selftest.json"))
            rc = 0 if st["all_ok"] else 4
        elif a.load_check or a.full:
            if not (a.piston_root and a.spurs1_root and a.out):
                ap.error("--piston-root、--spurs1-root、--out 必须给")
            rc = run_load_check(a) if a.load_check else run_full(a)
        else:
            ap.error("need --selftest / --load-check / --full")
    except Exception:
        log("异常：" + traceback.format_exc())
        rc = 3
    if a.out and (a.full or a.load_check):
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, "p7e_log.txt"), "a", encoding="utf-8") as f:
            f.write("\n".join(LOG) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
