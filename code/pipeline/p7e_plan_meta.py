#!/usr/bin/env python3
"""p7e_plan_meta.py — P7e（PISTON 2018/19 船载＋SPURS-1 WHOI 浮标）的元数据盘点：只读原始 netCDF 的元数据与可用性，不读盐度数值。

Description：
  不读盐度数值的元数据盘点，逻辑沿用 `p6_plan_meta.py`（import 其工具函数，一行未改），改动三处：
    1) 递归进入 netCDF4 group（ASDC 文件可能分组）；
    2) **比 P6 更严的盲**：盐度类变量（名字／long_name／standard_name 含 sal／sss／psal／cond）的数值型属性
       （valid_min、valid_max、actual_range、*_min、*_max、*range* 等）一律替换为 "<redacted: blind>"，不输出整文件极值；
       其余属性（units、long_name、深度说明等）照常输出；
    3) 雨量单位表补充（mm/hr 的若干写法）；累积型雨量（units 为 mm 且 long_name 含 accum／cumul）只报量级与单调性，不自动换算。
  盐度类变量只计「有限值」比例与逐小时可用性；雨与风给量级统计与多种事件定义下的事件数；每个事件标各盐度变量
  在 [t0−6,t0)、[t0,t0+6)、[t0,t0+12) 的有值小时数与平台位移——都属可用时段与平台运动，不涉及盐度数值。
Usage：python p7e_plan_meta.py --roots <bulk-disk>/.../raw/piston,<bulk-disk>/.../raw/spurs1 --out DIR   （产物 DIR/piston/、DIR/spurs1/）
Dependencies：numpy、netCDF4（Python 环境）＋同目录 p6_plan_meta.py。只读原件；产物 <out>/<根名>/p7e_plan.json、<out>/p7e_plan_log.txt。
Change Log：2026-09-27 初版。
"""

import argparse
import json
import os
import re
import sys
import time

import numpy as np
import netCDF4

import p6_plan_meta as pm

VERSION = "p7e-plan-2026-09-27a"
REDACT_ATTR = re.compile(r"(valid_min|valid_max|valid_range|actual_range|actual_min|actual_max|_min$|_max$|range|"
                         r"minimum|maximum|mean|median)", re.I)
EXTRA_MMH = ("mm hr-1", "mm hr^-1", "mm/hr", "mm h-1", "mm hour-1", "mm/h", "mm per hour")


def rain_factor(units):
    f = pm.rain_to_mmh(units)
    if f is not None:
        return f
    u = units.lower().strip()
    if u in EXTRA_MMH or u.replace(" ", "") in [x.replace(" ", "") for x in EXTRA_MMH]:
        return 1.0
    return None


def walk_groups(ds, prefix=""):
    yield prefix, ds
    for gn, g in ds.groups.items():
        yield from walk_groups(g, f"{prefix}{gn}/")


def scan_file(fp, rel, sal_avail, rain_series, pos_series):
    rec = {"bytes": os.path.getsize(fp)}
    try:
        ds = netCDF4.Dataset(fp)
    except Exception as e:
        rec["error"] = repr(e)
        return rec
    rec["global_attrs"] = pm.attrs_of(ds)
    rec["groups"] = {}
    for gp, g in walk_groups(ds):
        g.set_auto_maskandscale(True)
        grec = {"dims": {k: len(d) for k, d in g.dimensions.items()}}
        tvars = pm.find_time_var(g)
        grec["time_vars"] = tvars
        tcache = {}
        for tn in tvars:
            try:
                tu = pm.time_to_unix(g.variables[tn])
                tcache[g.variables[tn].dimensions[0]] = (tn, tu)
                fin = tu[np.isfinite(tu)]
                d = np.diff(np.sort(fin)) if len(fin) > 1 else np.array([np.nan])
                grec.setdefault("time_summary", {})[tn] = {
                    "n": int(len(tu)), "start": pm.iso(fin.min()) if len(fin) else None,
                    "end": pm.iso(fin.max()) if len(fin) else None,
                    "median_dt_s": float(np.median(d)) if len(fin) > 1 else None,
                    "n_gaps_gt_1h": int((d > 3600).sum()) if len(fin) > 1 else None,
                    "monotonic": bool(np.all(np.diff(fin) >= 0)) if len(fin) > 1 else None}
            except Exception as e:
                grec.setdefault("time_summary", {})[tn] = {"error": repr(e)[:300]}
        vars_out = {}
        lat_v = lon_v = None
        for n, v in g.variables.items():
            at = pm.attrs_of(v)
            cls = pm.classify(n, at)
            if cls == "salinity":
                at = {k: ("<redacted: blind>" if REDACT_ATTR.search(k) and not isinstance(val, str) else val)
                      for k, val in at.items()}
            info = {"dims": list(v.dimensions), "shape": list(v.shape), "dtype": str(v.dtype), "attrs": at, "class": cls}
            tdim = next((dn for dn in v.dimensions if dn in tcache), None)
            info["time_dim"] = tdim
            try:
                lname = n.lower()
                sname = str(at.get("standard_name", "")).lower()
                is_coord = (lname in pm.COORD_NAMES or sname in ("latitude", "longitude", "depth", "height", "sea_water_pressure")
                            or (v.ndim == 1 and n in g.dimensions and cls != "salinity"))
                if n in tvars:
                    pass
                elif cls == "salinity":
                    arr = pm.read_float(v)
                    fin = np.isfinite(arr)
                    info["finite_frac"] = round(float(fin.mean()), 4) if arr.size else None
                    if tdim is not None:
                        _, tu = tcache[tdim]
                        ax = v.dimensions.index(tdim)
                        fin_t = np.moveaxis(fin, ax, 0)
                        cols = {"all": fin_t} if fin_t.ndim == 1 else \
                            {f"idx{j}": fin_t.reshape(fin_t.shape[0], -1)[:, j] for j in range(int(np.prod(fin_t.shape[1:])))}
                        av = {}
                        for ck, m in cols.items():
                            hrs = pm.hourly_avail(tu, m)
                            sal_avail[f"{rel}::{gp}{n}::{ck}"] = hrs
                            segs = pm.segments(hrs, gap=1)
                            av[ck] = {"n_hours": int(len(hrs)), "n_samples": int(m.sum()), "segments_first40": segs[:40],
                                      "n_segments": len(segs)}
                        info["hourly_availability"] = av
                    del arr
                elif cls == "qc":
                    arr = pm.read_float(v)
                    fin = arr[np.isfinite(arr)]
                    u, c = np.unique(fin, return_counts=True)
                    info["value_counts"] = ({str(pm.jsafe(k)): int(cc) for k, cc in zip(u, c)} if len(u) <= 30
                                            else f"{len(u)} unique")
                elif cls in ("rain", "wind") or is_coord:
                    arr = pm.read_float(v)
                    fin = arr[np.isfinite(arr)]
                    info["finite_frac"] = round(float(np.isfinite(arr).mean()), 4) if arr.size else None
                    if len(fin):
                        info["summary"] = {"min": float(fin.min()), "max": float(fin.max()), "mean": float(fin.mean()),
                                           "q": [float(x) for x in np.quantile(fin, [0.01, 0.25, 0.5, 0.75, 0.99])]}
                        if is_coord and len(np.unique(fin)) <= 60:
                            info["unique"] = [float(x) for x in np.unique(fin)]
                        if cls == "rain":
                            info["frac_negative"] = round(float((fin < 0).mean()), 4)
                            info["monotonic_nondecreasing"] = bool(np.all(np.diff(fin) >= 0)) if len(fin) > 1 else None
                    if cls == "rain" and tdim is not None and v.ndim == 1:
                        f = rain_factor(str(at.get("units", "")))
                        info["mmh_factor"] = f
                        if f is not None:
                            hr = pm.hourly_rain(tcache[tdim][1], arr, f)
                            if hr is not None:
                                rain_series[f"{rel}::{gp}{n}"] = hr[:2]
                                info["hourly_median_dt_s"] = hr[2]
                    if lname in ("lat", "latitude") and v.ndim == 1 and tdim is not None:
                        lat_v = (tdim, arr)
                    if lname in ("lon", "longitude") and v.ndim == 1 and tdim is not None:
                        lon_v = (tdim, arr)
                else:
                    if v.size <= pm.MAX_ELEM and np.issubdtype(v.dtype, np.number):
                        arr = pm.read_float(v)
                        info["finite_frac"] = round(float(np.isfinite(arr).mean()), 4) if arr.size else None
            except Exception as e:
                info["error"] = repr(e)[:300]
            vars_out[n] = info
        if lat_v and lon_v and lat_v[0] == lon_v[0]:
            pos_series[f"{rel}::{gp}"] = (tcache[lat_v[0]][1], lat_v[1], lon_v[1])
        grec["variables"] = vars_out
        rec["groups"][gp or "/"] = grec
    ds.close()
    return rec


def plan_root(root, out):
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    files = sorted(os.path.join(dp, f) for dp, _, fs in os.walk(root) for f in fs if f.endswith(".nc"))
    pm.log(f"start {VERSION} root={root} files={len(files)}")
    plan = {"version": VERSION, "root": root,
            "blind_note": "盐度类变量只报有限值比例与逐小时可用性；其数值型范围属性（valid_min/max 等）已替换为 <redacted: blind>",
            "files": {}}
    sal_avail, rain_series, pos_series = {}, {}, {}
    for fp in files:
        rel = os.path.relpath(fp, root)
        pm.log(f"open {rel}")
        plan["files"][rel] = scan_file(fp, rel, sal_avail, rain_series, pos_series)
    evs_out = {}
    for rk, (h0, rr) in rain_series.items():
        per = {}
        for dry, thr in pm.EVENT_DEFS:
            ev = pm.find_events(h0, rr, dry, thr)
            for e in ev:
                H = e["onset_h"]
                av = {}
                for sk, hrs in sal_avail.items():
                    pre = int(np.sum((hrs >= H - 6) & (hrs < H)))
                    p6 = int(np.sum((hrs >= H) & (hrs < H + 6)))
                    p12 = int(np.sum((hrs >= H) & (hrs < H + 12)))
                    if pre or p6 or p12:
                        av[sk] = [pre, p6, p12]
                e["sal_hours_pre6_post6_post12"] = av
                disp = {}
                for pk, (tu, la, lo) in pos_series.items():
                    m1 = (tu >= (H - 6) * 3600) & (tu < H * 3600) & np.isfinite(la) & np.isfinite(lo)
                    m2 = (tu >= H * 3600) & (tu < (H + 6) * 3600) & np.isfinite(la) & np.isfinite(lo)
                    if m1.sum() and m2.sum():
                        disp[pk] = round(pm.haversine(float(np.median(la[m1])), float(np.median(lo[m1])),
                                                      float(np.median(la[m2])), float(np.median(lo[m2]))), 1)
                e["disp_km_pre6_to_post6"] = disp
            per[f"dry{dry}_thr{int(thr)}"] = {"n": len(ev), "events": ev}
        per["_series"] = {"start": pm.iso(h0 * 3600), "hours": len(rr), "valid_hours": int(np.isfinite(rr).sum()),
                          "wet_hours": int(np.sum(np.isfinite(rr) & (rr >= pm.RAIN_HOUR_MM))),
                          "total_mm_wet": round(float(np.nansum(np.where(rr >= pm.RAIN_HOUR_MM, rr, 0.0))), 1)}
        evs_out[rk] = per
    plan["rain_events"] = evs_out
    plan["salinity_series_hours"] = {k: int(len(v)) for k, v in sal_avail.items()}
    plan["position_series"] = {k: {"n": int(len(v[0]))} for k, v in pos_series.items()}
    plan["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(out, "p7e_plan.json"), "w") as f:
        json.dump(plan, f, ensure_ascii=False, indent=1, default=str)
    pm.log(f"done {root} {plan['elapsed_s']} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roots", required=True, help="逗号分隔；每个根的产物写 <out>/<根目录名>/")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for root in [x for x in a.roots.split(",") if x]:
        plan_root(root, os.path.join(a.out, os.path.basename(root.rstrip("/"))))
    with open(os.path.join(a.out, "p7e_plan_log.txt"), "w") as f:
        f.write("\n".join(pm.LOG) + "\n")


if __name__ == "__main__":
    sys.exit(main())
