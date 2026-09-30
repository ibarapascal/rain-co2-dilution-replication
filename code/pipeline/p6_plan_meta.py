#!/usr/bin/env python3
"""p6_plan_meta.py — P6（SPURS-2 近表层剖面检验）的元数据盘点：只读 SPURS-2 原始 netCDF 的元数据与可用性，不读盐度数值。

Description：
  不读盐度数值的元数据盘点。对 bulk `raw/spurs2/` 下每个 .nc：全局属性、维、每个变量的名字／维／形状／dtype／属性；
  坐标类（时间、经纬度、深度、压力、高度）给范围与唯一值；时间轴给起止、中位采样间隔、>1 h 的断档数；
  **盐度类变量（名字或属性含 sal／sss／psal／cond）只计「有限值」的比例与逐小时可用性，从不输出、统计或比较任何数值**；
  雨量类（rain／precip）给量级统计、按 P1 I4/I5/I6 小时化后在多种事件定义下的事件数与 onset 列表；风类给量级统计；
  QC 类整数变量给取值计数；其余变量只给有限值比例。再对每个事件、每个盐度变量（按深度下标分列）判 [t0−6,t0)、[t0,t0+6)、
  [t0,t0+12) 是否有值，以及平台经纬度在雨前／雨后 6 h 的位移（km）——都属于「数据可用时段」与平台运动，不涉及盐度数值。
Usage：python p6_plan_meta.py --root <raw_root>/spurs2 --out DIR
Dependencies：numpy、netCDF4（Python 环境）。只读原件，不写 raw、不改不删任何文件；产物只有 <out>/p6_plan.json 与 p6_plan_log.txt。
Change Log：2026-09-27 初版。
"""

import argparse
import json
import math
import os
import re
import sys
import time
from datetime import datetime, timezone

import numpy as np
import netCDF4

VERSION = "p6-plan-2026-09-27a"
SAL_PAT = re.compile(r"(sal|sss|psal|cond)", re.I)
RAIN_PAT = re.compile(r"(rain|precip|prcp|\bprc\b)", re.I)
WIND_PAT = re.compile(r"(wind|wspd|wnd|\bws\b|\bu10|\bspd)", re.I)
COORD_NAMES = {"time", "lat", "lon", "latitude", "longitude", "depth", "pres", "pressure", "height", "z", "depth_ssp"}
RAIN_HOUR_MM = 0.4
EVENT_DEFS = [(24, 10.0), (12, 10.0), (6, 10.0), (24, 5.0), (12, 5.0), (6, 5.0)]
MAX_ELEM = 40_000_000
LOG = []


def log(msg):
    line = f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} {msg}"
    LOG.append(line)
    print(line, flush=True)


def jsafe(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, np.ndarray):
        return [jsafe(v) for v in x.tolist()] if x.size <= 64 else f"<array {x.shape}>"
    if isinstance(x, bytes):
        return x.decode("utf-8", "replace")[:400]
    if isinstance(x, str):
        return x[:600]
    if isinstance(x, (list, tuple)):
        return [jsafe(v) for v in x][:64]
    return x


def attrs_of(obj):
    return {k: jsafe(obj.getncattr(k)) for k in obj.ncattrs()}


def read_float(v):
    """整变量读成 float64，填充值→NaN；元素过多时抛错（本 plan 只读 1–2 维小变量）。"""
    if v.size > MAX_ELEM:
        raise MemoryError(f"{v.name} size {v.size} > {MAX_ELEM}")
    v.set_auto_maskandscale(True)
    a = v[:]
    if np.ma.isMaskedArray(a):
        a = np.ma.filled(a.astype(float), np.nan)
    return np.asarray(a, dtype=float)


def time_to_unix(tv):
    units = getattr(tv, "units", "")
    cal = getattr(tv, "calendar", "standard")
    raw = read_float(tv)
    base = netCDF4.num2date(0.0, units, cal, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    one = netCDF4.num2date(1.0, units, cal, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    b = base.replace(tzinfo=timezone.utc).timestamp()
    s = (one.replace(tzinfo=timezone.utc).timestamp() - b)
    return b + raw * s


def iso(ts):
    return datetime.fromtimestamp(float(ts), timezone.utc).isoformat(timespec="seconds") if np.isfinite(ts) else None


def classify(name, attrs):
    txt = " ".join([name, str(attrs.get("long_name", "")), str(attrs.get("standard_name", ""))])
    if "flag_values" in attrs or "flag_meanings" in attrs or re.search(r"(_qc$|^qc_|qc$|flag)", name, re.I):
        return "qc"
    if SAL_PAT.search(txt) and not re.search(r"(time|lat|lon)", name, re.I):
        return "salinity"
    if RAIN_PAT.search(txt):
        return "rain"
    if WIND_PAT.search(txt) and not re.search(r"(dir|stress|tau)", txt, re.I):
        return "wind"
    return "other"


def find_time_var(ds):
    cands = []
    for n, v in ds.variables.items():
        u = str(getattr(v, "units", ""))
        if " since " in u and v.ndim == 1:
            cands.append((0 if n.lower() in ("time", "t") else 1, n))
    return [n for _, n in sorted(cands)]


def hourly_avail(tunix, mask):
    """mask：与时间轴同长的布尔（该时刻有有限值）；返回有值的 UTC 小时整数集合（排序 numpy 数组）。"""
    h = np.floor(tunix[mask] / 3600.0)
    h = h[np.isfinite(h)]
    return np.unique(h.astype(np.int64))


def segments(hours, gap=1):
    if len(hours) == 0:
        return []
    out, s, p = [], hours[0], hours[0]
    for h in hours[1:]:
        if h - p > gap:
            out.append([iso(s * 3600), iso((p + 1) * 3600), int(p - s + 1)])
            s = h
        p = h
    out.append([iso(s * 3600), iso((p + 1) * 3600), int(p - s + 1)])
    return out


def rain_to_mmh(units):
    u = units.lower().replace(" ", "")
    if u in ("mm/hr", "mm/h", "mmh-1", "mm/hour", "mmhr-1", "mm.h-1", "mmh^-1", "mm/hr.", "millimeters/hour", "millimeter/hour"):
        return 1.0
    if u in ("mm/min", "mmmin-1"):
        return 60.0
    if u in ("mm/s", "mms-1", "kgm-2s-1", "kg/m2/s", "kgm-2s^-1"):
        return 3600.0
    if u in ("m/s", "ms-1"):
        return 3.6e6
    return None


def hourly_rain(tunix, r, factor):
    """P1 I4：小时箱均值，需 ≥50% 子样本有效（按该箱内中位采样间隔推算期望样本数）。"""
    ok = np.isfinite(tunix) & np.isfinite(r)
    t, v = tunix[ok], r[ok] * factor
    if len(t) < 2:
        return None
    dt = float(np.median(np.diff(np.sort(t))))
    exp_n = max(1, int(round(3600.0 / dt))) if dt > 0 else 1
    h = np.floor(t / 3600.0).astype(np.int64)
    h0, h1 = int(h.min()), int(h.max())
    n = h1 - h0 + 1
    cnt = np.bincount(h - h0, minlength=n)
    sm = np.bincount(h - h0, weights=v, minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sm / cnt
    mean[cnt < 0.5 * exp_n] = np.nan
    return h0, mean, dt


def find_events(h0, rr, dry_h, thr):
    """P1 I6：onset＝有效有雨小时（≥0.4），前 dry_h 小时全有效且无雨；[t0,t0+24h) 雨小时累积 ≥thr 且 24 h 全有效。
    onset 之间要求：前 dry_h 小时无雨已保证间隔。"""
    ev = []
    n = len(rr)
    valid = np.isfinite(rr)
    wet = valid & (rr >= RAIN_HOUR_MM)
    for i in range(dry_h, n - 24):
        if not wet[i]:
            continue
        pre = slice(i - dry_h, i)
        if not (valid[pre].all() and not wet[pre].any()):
            continue
        post = slice(i, i + 24)
        if not valid[post].all():
            continue
        acc = float(np.where(wet[post], rr[post], 0.0).sum())
        if acc >= thr:
            ev.append({"onset_h": int(h0 + i), "onset": iso((h0 + i) * 3600), "acc24_mm": round(acc, 2),
                       "peak_mmh": round(float(np.nanmax(rr[post])), 2), "wet_hours_24": int(wet[post].sum())})
    return ev


def haversine(la1, lo1, la2, lo2):
    p = math.pi / 180
    a = math.sin((la2 - la1) * p / 2) ** 2 + math.cos(la1 * p) * math.cos(la2 * p) * math.sin((lo2 - lo1) * p / 2) ** 2
    return 2 * 6371.0 * math.asin(min(1.0, math.sqrt(a)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=__import__("repro_paths").raw_sub("spurs2"))  # 路径来自集中配置
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t_start = time.time()
    files = []
    for dp, _, fs in os.walk(a.root):
        for f in sorted(fs):
            if f.endswith(".nc"):
                files.append(os.path.join(dp, f))
    files.sort()
    log(f"start {VERSION} files={len(files)}")
    plan = {"version": VERSION, "root": a.root, "blind_note": "盐度类变量只报有限值比例与逐小时可用性，不含任何数值", "files": {}}
    sal_avail = {}      # key -> sorted hours
    rain_series = {}    # key -> (h0, hourly mm/h)
    pos_series = {}     # file -> (tunix, lat, lon)
    for fp in files:
        rel = os.path.relpath(fp, a.root)
        log(f"open {rel}")
        rec = {"bytes": os.path.getsize(fp)}
        try:
            ds = netCDF4.Dataset(fp)
        except Exception as e:
            rec["error"] = repr(e)
            plan["files"][rel] = rec
            continue
        ds.set_auto_maskandscale(True)
        rec["global_attrs"] = attrs_of(ds)
        rec["groups"] = list(ds.groups.keys())
        rec["dims"] = {k: len(d) for k, d in ds.dimensions.items()}
        tvars = find_time_var(ds)
        rec["time_vars"] = tvars
        tcache = {}
        for tn in tvars:
            try:
                tu = time_to_unix(ds.variables[tn])
                tcache[ds.variables[tn].dimensions[0]] = (tn, tu)
                fin = tu[np.isfinite(tu)]
                d = np.diff(np.sort(fin)) if len(fin) > 1 else np.array([np.nan])
                rec.setdefault("time_summary", {})[tn] = {
                    "n": int(len(tu)), "start": iso(fin.min()) if len(fin) else None, "end": iso(fin.max()) if len(fin) else None,
                    "median_dt_s": float(np.median(d)) if len(fin) > 1 else None,
                    "n_gaps_gt_1h": int((d > 3600).sum()) if len(fin) > 1 else None,
                    "monotonic": bool(np.all(np.diff(fin) >= 0)) if len(fin) > 1 else None}
            except Exception as e:
                rec.setdefault("time_summary", {})[tn] = {"error": repr(e)}
        vars_out = {}
        lat_v = lon_v = None
        for n, v in ds.variables.items():
            at = attrs_of(v)
            info = {"dims": list(v.dimensions), "shape": list(v.shape), "dtype": str(v.dtype), "attrs": at}
            cls = classify(n, at)
            info["class"] = cls
            tdim = next((dname for dname in v.dimensions if dname in tcache), None)
            info["time_dim"] = tdim
            try:
                lname = n.lower()
                sname = str(at.get("standard_name", "")).lower()
                is_coord = (lname in COORD_NAMES or sname in ("latitude", "longitude", "depth", "height", "sea_water_pressure")
                            or (v.ndim == 1 and n in ds.dimensions and cls not in ("salinity",)))
                if n in tvars:
                    pass
                elif cls == "salinity":
                    arr = read_float(v)
                    fin = np.isfinite(arr)
                    info["finite_frac"] = round(float(fin.mean()), 4) if arr.size else None
                    if tdim is not None:
                        tn, tu = tcache[tdim]
                        ax = v.dimensions.index(tdim)
                        fin_t = np.moveaxis(fin, ax, 0)
                        if fin_t.ndim == 1:
                            cols = {"all": fin_t}
                        else:
                            flat = fin_t.reshape(fin_t.shape[0], -1)
                            cols = {f"idx{j}": flat[:, j] for j in range(flat.shape[1])}
                        av = {}
                        for ck, m in cols.items():
                            hrs = hourly_avail(tu, m)
                            key = f"{rel}::{n}::{ck}"
                            sal_avail[key] = hrs
                            av[ck] = {"n_hours": int(len(hrs)), "n_samples": int(m.sum()), "segments_first40": segments(hrs, gap=1)[:40],
                                      "n_segments": len(segments(hrs, gap=1))}
                        info["hourly_availability"] = av
                elif cls == "qc":
                    arr = read_float(v)
                    fin = arr[np.isfinite(arr)]
                    u, c = np.unique(fin, return_counts=True)
                    info["value_counts"] = {str(jsafe(k)): int(cc) for k, cc in zip(u[:30], c[:30])} if len(u) <= 30 else f"{len(u)} unique"
                elif cls in ("rain", "wind") or is_coord:
                    arr = read_float(v)
                    fin = arr[np.isfinite(arr)]
                    info["finite_frac"] = round(float(np.isfinite(arr).mean()), 4) if arr.size else None
                    if len(fin):
                        info["summary"] = {"min": float(fin.min()), "max": float(fin.max()), "mean": float(fin.mean()),
                                           "q": [float(x) for x in np.quantile(fin, [0.01, 0.25, 0.5, 0.75, 0.99])]}
                        if is_coord and len(np.unique(fin)) <= 60:
                            info["unique"] = [float(x) for x in np.unique(fin)]
                        if cls == "rain":
                            info["frac_negative"] = round(float((fin < 0).mean()), 4)
                    if cls == "rain" and tdim is not None and v.ndim == 1:
                        f = rain_to_mmh(str(at.get("units", "")))
                        info["mmh_factor"] = f
                        if f is not None:
                            hr = hourly_rain(tcache[tdim][1], arr, f)
                            if hr is not None:
                                rain_series[f"{rel}::{n}"] = hr[:2]
                                info["hourly_median_dt_s"] = hr[2]
                    if lname in ("lat", "latitude") and v.ndim == 1 and tdim is not None:
                        lat_v = (tdim, arr)
                    if lname in ("lon", "longitude") and v.ndim == 1 and tdim is not None:
                        lon_v = (tdim, arr)
                else:
                    if v.size <= MAX_ELEM:
                        arr = read_float(v) if np.issubdtype(v.dtype, np.number) else None
                        if arr is not None:
                            info["finite_frac"] = round(float(np.isfinite(arr).mean()), 4) if arr.size else None
            except Exception as e:
                info["error"] = repr(e)[:300]
            vars_out[n] = info
        if lat_v and lon_v and lat_v[0] == lon_v[0]:
            pos_series[rel] = (tcache[lat_v[0]][1], lat_v[1], lon_v[1])
        rec["variables"] = vars_out
        ds.close()
        plan["files"][rel] = rec
    # ---- 雨事件（多种定义），并对每个事件标各盐度变量的可用性与平台位移
    evs_out = {}
    for rk, (h0, rr) in rain_series.items():
        per = {}
        for dry, thr in EVENT_DEFS:
            ev = find_events(h0, rr, dry, thr)
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
                        disp[pk] = round(haversine(float(np.median(la[m1])), float(np.median(lo[m1])),
                                                   float(np.median(la[m2])), float(np.median(lo[m2]))), 1)
                e["disp_km_pre6_to_post6"] = disp
            per[f"dry{dry}_thr{int(thr)}"] = {"n": len(ev), "events": ev}
        valid_h = int(np.isfinite(rr).sum())
        per["_series"] = {"start": iso(h0 * 3600), "hours": len(rr), "valid_hours": valid_h,
                          "wet_hours": int(np.sum(np.isfinite(rr) & (rr >= RAIN_HOUR_MM))),
                          "total_mm_wet": round(float(np.nansum(np.where(rr >= RAIN_HOUR_MM, rr, 0.0))), 1)}
        evs_out[rk] = per
    plan["rain_events"] = evs_out
    plan["position_series"] = {k: {"n": int(len(v[0]))} for k, v in pos_series.items()}
    plan["elapsed_s"] = round(time.time() - t_start, 1)
    with open(os.path.join(a.out, "p6_plan.json"), "w") as f:
        json.dump(plan, f, ensure_ascii=False, indent=1, default=str)
    with open(os.path.join(a.out, "p6_plan_log.txt"), "w") as f:
        f.write("\n".join(LOG) + "\n")
    log(f"done {plan['elapsed_s']} s")


if __name__ == "__main__":
    sys.exit(main())
