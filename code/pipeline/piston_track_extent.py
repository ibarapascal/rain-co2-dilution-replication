"""PISTON 2018/2019 船载航迹范围提取（逐小时船位与范围，供站点地图用；不参与任何估计量）。

只读 bulk raw/piston 两个航次的 nav-met-sea-flux-60min 原件，取逐小时船位，
输出每航次经纬度范围（min/max 与 5/95 百分位）与逐小时航迹 CSV。只读原件、不改不删。
用法：python piston_track_extent.py --root <raw/piston> --out <输出目录>
"""
import argparse
import csv
import glob
import hashlib
import json
import os

import numpy as np
from netCDF4 import Dataset

LAT_NAMES = ("lat", "latitude", "LAT", "Latitude", "gps_lat", "lat_gps")
LON_NAMES = ("lon", "longitude", "LON", "Longitude", "gps_lon", "lon_gps")


def find_var(ds, names):
    groups = [ds] + list(ds.groups.values())
    for g in groups:
        for n in names:
            if n in g.variables:
                return g.variables[n]
    for g in groups:
        for n, v in g.variables.items():
            sn = str(getattr(v, "standard_name", "")).lower()
            if (names is LAT_NAMES and sn == "latitude") or (names is LON_NAMES and sn == "longitude"):
                return v
    raise KeyError(names)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    summary = {}
    with open(os.path.join(a.out, "piston_track_hourly.csv"), "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["cruise", "i", "lat", "lon"])
        for path in sorted(glob.glob(os.path.join(a.root, "*", "PISTON-nav-met-sea-flux-60min_*.nc"))):
            cruise = os.path.basename(os.path.dirname(path))
            with Dataset(path) as ds:
                lat = np.ma.filled(find_var(ds, LAT_NAMES)[:].astype(float), np.nan).ravel()
                lon = np.ma.filled(find_var(ds, LON_NAMES)[:].astype(float), np.nan).ravel()
            ok = np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90) & (np.abs(lon) <= 360)
            lat, lon = lat[ok], lon[ok]
            for i, (y, x) in enumerate(zip(lat, lon)):
                w.writerow([cruise, i, f"{y:.4f}", f"{x:.4f}"])
            summary[cruise] = {
                "file": os.path.basename(path), "sha256": sha256(path), "n_hours": int(ok.sum()),
                "lat_min": float(lat.min()), "lat_max": float(lat.max()),
                "lon_min": float(lon.min()), "lon_max": float(lon.max()),
                "lat_p5_p95": [float(np.percentile(lat, 5)), float(np.percentile(lat, 95))],
                "lon_p5_p95": [float(np.percentile(lon, 5)), float(np.percentile(lon, 95))],
            }
    with open(os.path.join(a.out, "piston_track_extent.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
