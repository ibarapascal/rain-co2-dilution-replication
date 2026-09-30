#!/usr/bin/env python3
"""p8b_plan_meta.py — P8b（SPURS-2 Lady Amber 同平台「海面／1 m」层比）的元数据盘点：只读元数据与可用性，不读盐度数值。

Description：
  不读盐度数值的元数据盘点。对 Lady Amber 6 个 .nc：调用 `p7e_plan_meta.scan_file`（一行未改；其内部复用
  `p6_plan_meta`）得到变量名／维／属性／时间覆盖／QC 取值计数／雨与风量级／盐度类变量逐小时可用性；在此之上：
    1) **更严的盲**：除 p7e 已屏蔽的盐度类变量数值型范围属性外，名字或 long_name 含 sal／sss／psal／cond／dens／sigma
       的任何变量，及全局属性里键名含 sal／sss 的数值属性，一律替换为 "<redacted: blind>"；
    2) 雨量变体：P6 plan 发现 Lady Amber 雨量率「无雨时存为缺测」，严格口径（小时 ≥50% 有效）会漏掉干期。另给
       nan0 口径（时间戳存在的样本里雨量 NaN 记 0，小时 ≥50% 期望样本数即有效）与累积量差分口径（units 为 mm 的累积量：
       小时内正增量之和，单步跳变 >50 mm 视为复位不计）的小时序列、NaN／零／正值占比；
    3) 借用雨量：中心浮标 `SPURS2_WHOI_Met.nc` rainfall_rate（10.05°N、125.03°W）与 Revelle 两航次 METEO precipitation_rate，
       按 Lady Amber 逐小时船位算离浮标、离 Revelle 的距离，统计 25／50／100 km 内的小时数；
    4) 事件：各雨源 × 6 种事件定义（p6_plan_meta.EVENT_DEFS）下的 onset，逐事件给各盐度序列在 [t0−6,t0)、[t0,t0+6)、
       [t0,t0+12) 的有值小时数、船位位移、借用雨源时 onset 前后 6 h 船离雨源的最大距离；并给同一文件内任两个盐度序列
       「前后窗都有值」的事件数（只数小时，不涉及数值）。
  盐度类变量只算有限值比例与逐小时可用性，**不输出、不统计、不比较任何盐度数值，不算任何淡化量**。
Usage：python p8b_plan_meta.py --la-root RAW/spurs2/SPURS2_LADYAMBER --spurs2-root RAW/spurs2 --out DIR
Dependencies：numpy、netCDF4（Python 环境）＋同目录 p7e_plan_meta.py、p6_plan_meta.py（只 import 不改）。
  只读原件；产物 <out>/p8b_plan.json、<out>/p8b_plan_summary.json、<out>/p8b_plan_log.txt。
Change Log：2026-09-27 初版。
"""

import argparse
import json
import math
import os
import re
import sys
import time
from itertools import combinations

import numpy as np
import netCDF4

import p6_plan_meta as pm
import p7e_plan_meta as p7

VERSION = "p8b-plan-2026-09-27a"
BLIND_NAME = re.compile(r"(sal|sss|psal|cond|dens|sigma)", re.I)
MOOR_LATLON = (10.05, -125.03)          # 名义；若文件里有 lat/lon 变量或属性则以文件为准
DIST_BINS = (25.0, 50.0, 100.0)
BORROW_MAX_KM = 50.0


def redact_more(rec):
    for gk, g in rec.get("groups", {}).items():
        for n, info in g.get("variables", {}).items():
            at = info.get("attrs", {})
            txt = " ".join([n, str(at.get("long_name", "")), str(at.get("standard_name", ""))])
            if BLIND_NAME.search(txt):
                info["attrs"] = {k: ("<redacted: blind>" if p7.REDACT_ATTR.search(k) and not isinstance(v, str) else v)
                                 for k, v in at.items()}
    ga = rec.get("global_attrs", {})
    rec["global_attrs"] = {k: ("<redacted: blind>" if re.search(r"(sal|sss)", k, re.I) and not isinstance(v, str) else v)
                           for k, v in ga.items()}
    return rec


def hourly_by_stamp(tu, r, factor, accum=False):
    """nan0／累积口径：以时间戳存在的样本数判小时有效（≥50% 期望数）。返回 (h0, 小时 mm/h, 诊断)。"""
    ok_t = np.isfinite(tu)
    t = tu[ok_t]
    v = r[ok_t]
    if len(t) < 2:
        return None
    order = np.argsort(t)
    t, v = t[order], v[order]
    dt = float(np.median(np.diff(t)))
    exp_n = max(1, int(round(3600.0 / dt))) if dt > 0 else 1
    h = np.floor(t / 3600.0).astype(np.int64)
    h0, h1 = int(h.min()), int(h.max())
    n = h1 - h0 + 1
    cnt = np.bincount(h - h0, minlength=n)
    diag = {"median_dt_s": dt, "frac_nan": round(float(np.isnan(v).mean()), 4)}
    fin = v[np.isfinite(v)]
    diag["frac_zero_of_finite"] = round(float((fin == 0).mean()), 4) if len(fin) else None
    diag["frac_pos_of_finite"] = round(float((fin > 0).mean()), 4) if len(fin) else None
    if accum:
        d = np.diff(v)
        d = np.where(np.isfinite(d) & (d > 0) & (d <= 50.0), d, 0.0)
        diag["n_resets_or_jumps"] = int(np.sum(np.isfinite(np.diff(v)) & ((np.diff(v) < 0) | (np.diff(v) > 50.0))))
        amt = np.bincount(h[1:] - h0, weights=d, minlength=n)
        rate = amt.astype(float)
    else:
        vv = np.where(np.isfinite(v), v, 0.0) * factor
        sm = np.bincount(h - h0, weights=vv, minlength=n)
        with np.errstate(invalid="ignore", divide="ignore"):
            rate = sm / cnt
    rate = np.asarray(rate, float)
    rate[cnt < 0.5 * exp_n] = np.nan
    diag["valid_hours"] = int(np.isfinite(rate).sum())
    diag["wet_hours"] = int(np.sum(np.isfinite(rate) & (rate >= pm.RAIN_HOUR_MM)))
    return h0, rate, diag


def hourly_pos(tu, la, lo):
    ok = np.isfinite(tu) & np.isfinite(la) & np.isfinite(lo)
    h = np.floor(tu[ok] / 3600.0).astype(np.int64)
    if len(h) == 0:
        return {}
    out = {}
    hs, idx = np.unique(h, return_index=True)
    las, los = la[ok], lo[ok]
    for j, hh in enumerate(hs):
        a = idx[j]
        b = idx[j + 1] if j + 1 < len(idx) else len(h)
        out[int(hh)] = (float(np.median(las[a:b])), float(np.median(los[a:b])))
    return out


def read_1d(ds, name):
    v = ds.variables[name]
    v.set_auto_maskandscale(True)
    a = v[:]
    return np.asarray(np.ma.filled(a.astype(float), np.nan), float)


def la_extra(fp, rel):
    """Lady Amber 单文件：雨量变体小时序列与逐小时船位。"""
    ds = netCDF4.Dataset(fp)
    ds.set_auto_maskandscale(True)
    tv = pm.find_time_var(ds)
    tcache = {}
    for tn in tv:
        try:
            tcache[ds.variables[tn].dimensions[0]] = pm.time_to_unix(ds.variables[tn])
        except Exception:
            pass
    rains, diags, pos = {}, {}, {}
    lat_n = lon_n = None
    for n, v in ds.variables.items():
        at = pm.attrs_of(v)
        cls = pm.classify(n, at)
        if n.lower() in ("lat", "latitude"):
            lat_n = n
        if n.lower() in ("lon", "longitude"):
            lon_n = n
        if cls != "rain" or v.ndim != 1 or v.dimensions[0] not in tcache:
            continue
        units = str(at.get("units", ""))
        tu = tcache[v.dimensions[0]]
        r = read_1d(ds, n)
        f = p7.rain_factor(units)
        u = units.lower().replace(" ", "")
        if f is not None:
            hr = hourly_by_stamp(tu, r, f)
            if hr:
                rains[f"{rel}::{n}::nan0"] = hr[:2]
                diags[f"{rel}::{n}::nan0"] = dict(hr[2], units=units)
        elif u in ("mm", "millimeters", "millimeter", "kgm-2"):
            hr = hourly_by_stamp(tu, r, 1.0, accum=True)
            if hr:
                rains[f"{rel}::{n}::accum"] = hr[:2]
                diags[f"{rel}::{n}::accum"] = dict(hr[2], units=units)
        else:
            diags[f"{rel}::{n}::unhandled"] = {"units": units}
    if lat_n and lon_n and ds.variables[lat_n].ndim == 1 and ds.variables[lat_n].dimensions[0] in tcache:
        pos = hourly_pos(tcache[ds.variables[lat_n].dimensions[0]], read_1d(ds, lat_n), read_1d(ds, lon_n))
    ds.close()
    return rains, diags, pos


def ext_source(fp, rain_name, latlon_fixed=None):
    """借用雨源：严格小时化（P1 I4，同 P6）＋逐小时位置（船）或固定位置（浮标）。"""
    ds = netCDF4.Dataset(fp)
    ds.set_auto_maskandscale(True)
    tn = pm.find_time_var(ds)[0]
    tu = pm.time_to_unix(ds.variables[tn])
    r = np.maximum(read_1d(ds, rain_name), 0.0)
    f = p7.rain_factor(str(getattr(ds.variables[rain_name], "units", "")))
    hr = pm.hourly_rain(tu, r, f if f else 1.0)
    pos = {}
    fixed = latlon_fixed
    names = {n.lower(): n for n in ds.variables}
    if "latitude" in names or "lat" in names:
        la = read_1d(ds, names.get("latitude", names.get("lat")))
        lo = read_1d(ds, names.get("longitude", names.get("lon")))
        if la.size == tu.size and la.size > 1:
            pos = hourly_pos(tu, la, lo)
            fixed = None
        elif la.size >= 1:
            fixed = (float(np.nanmedian(la)), float(np.nanmedian(lo)))
    ds.close()
    return hr[0], hr[1], pos, fixed, {"units_factor": f}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--la-root", required=True)
    ap.add_argument("--spurs2-root", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t_start = time.time()
    files = sorted(os.path.join(dp, f) for dp, _, fs in os.walk(a.la_root) for f in fs if f.endswith(".nc"))
    pm.log(f"start {VERSION} files={len(files)}")
    plan = {"version": VERSION, "la_root": a.la_root,
            "blind_note": "盐度类变量只报有限值比例与逐小时可用性；盐度／电导／密度类变量与全局属性的数值型范围属性已替换为 <redacted: blind>",
            "files": {}}
    sal_avail, rain_strict, pos_series = {}, {}, {}
    rain_var, rain_diag, pos_h = {}, {}, {}
    for fp in files:
        rel = os.path.relpath(fp, a.la_root).replace("\\", "/")
        pm.log(f"open {rel}")
        rec = p7.scan_file(fp, rel, sal_avail, rain_strict, pos_series)
        plan["files"][rel] = redact_more(rec)
        rs, dg, ph = la_extra(fp, rel)
        rain_var.update(rs)
        rain_diag.update(dg)
        pos_h[rel] = ph
    for k, (h0, rr) in rain_strict.items():
        rain_var[k + "::strict"] = (h0, rr)
    plan["rain_variant_diag"] = rain_diag
    # ---- 借用雨源
    src = {}
    try:
        h0, rr, pos, fixed, dg = ext_source(os.path.join(a.spurs2_root, "SPURS2_MOORING_CENTRAL", "SPURS2_WHOI_Met.nc"),
                                            "rainfall_rate", MOOR_LATLON)
        src["MOOR"] = {"h0": h0, "rr": rr, "pos": pos, "fixed": fixed or MOOR_LATLON}
    except Exception as e:
        plan["moor_error"] = repr(e)[:300]
    for cr in ("RR1610", "RR1720"):
        try:
            h0, rr, pos, fixed, dg = ext_source(os.path.join(a.spurs2_root, "SPURS2_METEO", f"SPURS2_{cr}_WHOI_Met.nc"),
                                                "precipitation_rate")
            src[cr] = {"h0": h0, "rr": rr, "pos": pos, "fixed": fixed}
        except Exception as e:
            plan[f"{cr}_error"] = repr(e)[:300]

    def src_dist(s, H, rel):
        p = pos_h[rel].get(H)
        if p is None:
            return None
        q = s["fixed"] if s["fixed"] else s["pos"].get(H)
        if q is None:
            return None
        return pm.haversine(p[0], p[1], q[0], q[1])

    prox = {}
    for rel in pos_h:
        hrs = sorted(pos_h[rel])
        prox[rel] = {"hours_with_position": len(hrs),
                     "lat_range": [min(v[0] for v in pos_h[rel].values()), max(v[0] for v in pos_h[rel].values())] if hrs else None,
                     "lon_range": [min(v[1] for v in pos_h[rel].values()), max(v[1] for v in pos_h[rel].values())] if hrs else None,
                     "start": pm.iso(hrs[0] * 3600) if hrs else None, "end": pm.iso((hrs[-1] + 1) * 3600) if hrs else None}
        for sk, s in src.items():
            d = [src_dist(s, H, rel) for H in hrs]
            d = np.array([x for x in d if x is not None], float)
            prox[rel][sk] = {f"hours_within_{int(b)}km": int((d <= b).sum()) for b in DIST_BINS}
            prox[rel][sk]["n_hours_dist_defined"] = int(len(d))
            prox[rel][sk]["min_km"] = round(float(d.min()), 1) if len(d) else None
    plan["proximity"] = prox
    for sk, s in src.items():
        rain_var[f"BORROW_{sk}"] = (s["h0"], s["rr"])

    # ---- 事件与盐度可用性
    sal_keys = sorted(sal_avail)
    evs_out = {}
    for rk, (h0, rr) in rain_var.items():
        borrowed = rk.startswith("BORROW_")
        la_file = None if borrowed else rk.split("::")[0]
        per = {}
        for dry, thr in pm.EVENT_DEFS:
            ev = pm.find_events(h0, rr, dry, thr)
            kept = []
            for e in ev:
                H = e["onset_h"]
                files_here = [la_file] if la_file else list(pos_h)
                if borrowed:
                    s = src[rk.replace("BORROW_", "")]
                    best = None
                    for rel in files_here:
                        ds_ = [src_dist(s, h, rel) for h in range(H - 6, H + 6)]
                        if all(x is not None for x in ds_):
                            mx = max(ds_)
                            if best is None or mx < best[1]:
                                best = (rel, mx)
                    if best is None:
                        continue
                    e["la_file"] = best[0]
                    e["max_km_to_source_pm6h"] = round(best[1], 1)
                    files_here = [best[0]]
                av = {}
                for sk in sal_keys:
                    if sk.split("::")[0] not in files_here:
                        continue
                    hrs = sal_avail[sk]
                    pre = int(np.sum((hrs >= H - 6) & (hrs < H)))
                    p6 = int(np.sum((hrs >= H) & (hrs < H + 6)))
                    p12 = int(np.sum((hrs >= H) & (hrs < H + 12)))
                    if pre or p6 or p12:
                        av[sk] = [pre, p6, p12]
                e["sal_hours_pre6_post6_post12"] = av
                rel = files_here[0]
                p1 = [pos_h[rel][h] for h in range(H - 6, H) if h in pos_h[rel]]
                p2_ = [pos_h[rel][h] for h in range(H, H + 6) if h in pos_h[rel]]
                if p1 and p2_:
                    e["disp_km_pre6_to_post6"] = round(pm.haversine(float(np.median([x[0] for x in p1])), float(np.median([x[1] for x in p1])),
                                                                    float(np.median([x[0] for x in p2_])), float(np.median([x[1] for x in p2_]))), 1)
                kept.append(e)
            single = {sk: sum(1 for e in kept if sk in e["sal_hours_pre6_post6_post12"]
                              and e["sal_hours_pre6_post6_post12"][sk][0] >= 1 and e["sal_hours_pre6_post6_post12"][sk][1] >= 1)
                      for sk in sal_keys}
            pair = {}
            for k1, k2 in combinations(sal_keys, 2):
                if k1.split("::")[0] != k2.split("::")[0]:
                    continue
                c = 0
                for e in kept:
                    A = e["sal_hours_pre6_post6_post12"]
                    if k1 in A and k2 in A and A[k1][0] >= 1 and A[k1][1] >= 1 and A[k2][0] >= 1 and A[k2][1] >= 1:
                        c += 1
                if c:
                    pair[f"{k1} | {k2}"] = c
            per[f"dry{dry}_thr{int(thr)}"] = {"n": len(kept), "n_raw_onsets": len(ev),
                                              "n_with_layer_pre_and_post6": {k: v for k, v in single.items() if v},
                                              "n_with_pair_pre_and_post6": pair, "events": kept}
        per["_series"] = {"start": pm.iso(h0 * 3600), "hours": len(rr), "valid_hours": int(np.isfinite(rr).sum()),
                          "wet_hours": int(np.sum(np.isfinite(rr) & (rr >= pm.RAIN_HOUR_MM)))}
        evs_out[rk] = per
    plan["rain_events"] = evs_out
    plan["salinity_series_hours"] = {k: int(len(v)) for k, v in sal_avail.items()}
    plan["elapsed_s"] = round(time.time() - t_start, 1)
    with open(os.path.join(a.out, "p8b_plan.json"), "w", encoding="utf-8") as f:
        json.dump(plan, f, ensure_ascii=False, indent=1, default=str)
    # 紧凑摘要（不含逐事件列表）
    summ = {"version": VERSION, "salinity_series_hours": plan["salinity_series_hours"], "proximity": prox,
            "rain_variant_diag": rain_diag,
            "events": {rk: {dk: ({kk: vv for kk, vv in dv.items() if kk != "events"} if isinstance(dv, dict) else dv)
                            for dk, dv in per.items()} for rk, per in evs_out.items()}}
    with open(os.path.join(a.out, "p8b_plan_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1, default=str)
    pm.log(f"done {plan['elapsed_s']} s")
    with open(os.path.join(a.out, "p8b_plan_log.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(pm.LOG) + "\n")


if __name__ == "__main__":
    sys.exit(main())
