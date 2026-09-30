"""P3 补充下载（主下载 p3-fetch 不含的两样；只下载与格式转换，不做 P3 计算）。

用途：P3 全球重算的输入；函数全部 import 自 p3_fetch.py（同目录），不复制实现。
  1. era5_ws10/200101.nc：ERA5 u10/v10 2001-01（NCAR S3 镜像，URL 规律同主下载）→ 与主下载完全相同的 derive_ws10
     （60°S–60°N、ws10 int16×0.01、纬度降序、zlib）先转整月到临时文件，再原样截取 2001-01-01 当天 24 个时次
     （int16 原值与属性照抄，不重新量化）；供 2000-12-31 23:30 半步用 2001-01-01 00Z（L32）。
     HDF5 分块经 HTTP Range 只取当天需要解析 HDF5 元数据，环境里没有 fsspec，故下整月（约 3 GB）再截，原文件与整月临时文件用完即删。
  2. woa09/basin.msk：WOA09 1° 分洋盆掩膜（NCEI 公开免登录，17,534,880 B，FORTRAN 文本 10F8.0，按层排列，
     每层 64,800 值＝纬度 −89.5→89.5 × 经度 0.5→359.5，经度最快；第一层＝表层）→ 另写 woa09/WOA09_basin_surface_1deg.csv
     （latitude,longitude,basin_code；−100＝陆地），供 p3_global.load_basins（分洋盆累计）读取。
     注意：load_basins 取 sorted(glob('woa09/**/*basin*'))[0]；CSV 名以大写 W 开头，按字节序排在 basin.msk 之前。
清单：写独立的 manifest_extra.jsonl（主下载可能同时在追加 manifest.jsonl，不共写一个文件）；字段同主清单。
输出：REPRO_OUTPUT_DIR/fetch_extra_summary.json（完成契约）或 fetch_extra_status.json（未完成）；日志 fetch_extra_log.txt。
退出码：0 全部完成；2 有失败（可原样重提续跑，manifest_extra 断点续传）；3 前置条件（磁盘、参数）。
用法：p3_fetch_extra.py --cache <fast_root>/p3-cache
版本：p3-fetch-extra-2026-09-26a
"""
import argparse, collections, concurrent.futures as cf, json, os, shutil, sys, threading, time, traceback

import p3_fetch as F

VERSION = "p3-fetch-extra-2026-09-26a"
ERA5_YM = 200101
KEEP_STEPS = 24  # 2001-01-01 00Z–23Z
WOA09_URL = "https://www.ncei.noaa.gov/data/oceans/woa/WOA09/MASKS/basin.msk"
WOA09_REL = "woa09/basin.msk"
WOA09_CSV = "woa09/WOA09_basin_surface_1deg.csv"
WOA09_BYTES = 17534880
NLAT, NLON = 180, 360


def extra_manifest(root):
    """F.Manifest 的同款对象，只是落在 manifest_extra.jsonl。"""
    m = F.Manifest.__new__(F.Manifest)
    m.root, m.path, m.recs, m.lock = root, os.path.join(root, "manifest_extra.jsonl"), {}, threading.Lock()
    if os.path.exists(m.path):
        with open(m.path, encoding="utf-8") as f:
            for line in f:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                m.recs[r["path"]] = r
    return m


# ======================================================================== ERA5 200101

def subset_day(full_path, out_path, nkeep):
    """整月 ws10 → 前 nkeep 个时次；int16 原值、变量与全局属性、压缩与分块照抄。返回校验摘要。"""
    import numpy as np, netCDF4
    tmp = out_path + ".tmp"
    with F.NC_LOCK:
        src = netCDF4.Dataset(full_path)
        try:
            dst = netCDF4.Dataset(tmp, "w", format="NETCDF4")
            for name, dim in src.dimensions.items():
                dst.createDimension(name, nkeep if name == "time" else len(dim))
            for name, v in src.variables.items():
                v.set_auto_maskandscale(False)
                filt = v.filters() or {}
                ch = v.chunking()
                o = dst.createVariable(name, v.dtype, v.dimensions, zlib=bool(filt.get("zlib")),
                                       complevel=filt.get("complevel") or 4, shuffle=bool(filt.get("shuffle")),
                                       chunksizes=ch if isinstance(ch, list) else None,
                                       fill_value=getattr(v, "_FillValue", None))
                o.set_auto_maskandscale(False)
                o.setncatts({k: v.getncattr(k) for k in v.ncattrs() if k != "_FillValue"})
                o[:] = v[:nkeep] if v.dimensions and v.dimensions[0] == "time" else v[:]
            dst.setncatts({**{k: src.getncattr(k) for k in src.ncattrs()},
                           "subset": f"time[0:{nkeep}] of the full-month ws10 derived by p3_fetch.derive_ws10 "
                                     f"(int16 values copied verbatim)", "subset_by": VERSION})
            dst.close()
            chk = {}
            with netCDF4.Dataset(tmp) as d:
                for name in ("time", "latitude", "longitude", "ws10"):
                    a, b = src.variables[name], d.variables[name]
                    a.set_auto_maskandscale(False)
                    b.set_auto_maskandscale(False)
                    sa = a[:nkeep] if a.dimensions[0] == "time" else a[:]
                    assert np.array_equal(np.asarray(sa), np.asarray(b[:])), f"{name} 截取后不一致"
                    assert {k: str(a.getncattr(k)) for k in a.ncattrs()} == {k: str(b.getncattr(k)) for k in b.ncattrs()}, name
                chk = {"subset_steps": nkeep, "full_steps": len(src.dimensions["time"]), "subset_equal_raw": True,
                       "ws10_filters": d.variables["ws10"].filters(), "ws10_chunking": d.variables["ws10"].chunking()}
        finally:
            src.close()
    os.replace(tmp, out_path)
    return chk


def verify_against_uv(path, u_path, v_path, nt):
    """用 p3_fetch.verify_ws10 把截取后的文件（第 0、nt//2、nt−1 步）与 u/v 原文件重算比对。"""
    import numpy as np, netCDF4
    with F.NC_LOCK:
        du, dv = netCDF4.Dataset(u_path), netCDF4.Dataset(v_path)
    try:
        with F.NC_LOCK:
            uvar, vvar = F.find_var(du, "10U"), F.find_var(dv, "10V")
            lat = du.variables[uvar.dimensions[1]][:]
            for v in (uvar, vvar):
                v.set_auto_mask(False)
        sel = np.where((lat >= -60 - 1e-6) & (lat <= 60 + 1e-6))[0]
        return F.verify_ws10(path, uvar, vvar, int(sel[0]), int(sel[-1]) + 1, nt)
    finally:
        with F.NC_LOCK:
            du.close()
            dv.close()


def do_era5(root, man, stats, pool):
    et = F.era5_task(ERA5_YM)
    if man.have(et["rel"]):
        print(F.log(f"[ERA5] {et['rel']} 已在 manifest_extra，跳过"), flush=True)
        return True
    for k in ("u", "v"):
        t = et["raw"][k]
        if not man.have(t["rel"]):
            F.s3_download(t, root, man, stats, pool)
    raw = {k: man.recs[et["raw"][k]["rel"]] for k in ("u", "v")}
    upath, vpath = (os.path.join(root, et["raw"][k]["rel"]) for k in ("u", "v"))
    full = os.path.join(root, "era5_raw", f"{ERA5_YM}_ws10_fullmonth.nc")
    inputs = [{"url": r["url"], "path": r["path"], "bytes": r["bytes"], "sha256": r["sha256"]} for r in raw.values()]
    t0 = time.monotonic()
    chk = F.derive_ws10(upath, vpath, full, expect_nt=24 * len(F.days_of_month(*divmod(ERA5_YM, 100))),
                        meta={"inputs": [{k: x[k] for k in ("url", "sha256", "bytes")} for x in inputs]})
    out = os.path.join(root, et["rel"])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    chk["subset"] = subset_day(full, out, KEEP_STEPS)
    chk["subset_verify"] = verify_against_uv(out, upath, vpath, KEEP_STEPS)
    man.add(source="era5_ws10", path=et["rel"], bytes=os.path.getsize(out), sha256=F.sha256_file(out), status="ok",
            inputs=inputs, check=chk, note=f"only 2001-01-01 00-23Z ({KEEP_STEPS} steps) kept; clarify L32")
    os.remove(full)
    # [repro] 上一行删的是整月 ws10 临时文件（派生物，非原件），照原逻辑；ERA5 u/v 原件保留（原为逐个 os.remove 并记
    # [repro] deleted_after_derive）；raw_sync.py A6 复制进 bulk raw 档
    print(F.log(f"[ERA5] {ERA5_YM} 转换＋截取完成 {time.monotonic() - t0:.0f} s，{json.dumps(chk, ensure_ascii=False)}"),
          flush=True)
    return True


# ======================================================================== WOA09

def woa09_surface(msk_path):
    """basin.msk 第一层（表层）→ [(lat, lon, code)]；按布局与地理抽查断言。"""
    vals = []
    with open(msk_path, encoding="ascii") as f:
        for line in f:
            vals.extend(float(x) for x in line.split())
            if len(vals) >= NLAT * NLON:
                break
    assert len(vals) >= NLAT * NLON, f"值不足一层：{len(vals)}"
    s = vals[:NLAT * NLON]
    at = lambda la, lo: s[int(la + 89.5) * NLON + int(lo - 0.5)]
    assert all(x == -100.0 for x in s[:NLON]), "−89.5° 纬圈应全为陆地（−100）"
    # 地理抽查（2026-09-26 从文件头 1.1 MB 核过）：赤道太平洋 2、大西洋 1、印度洋 3；40.5°N 日本海 12
    for la, lo, want in ((0.5, 154.5, 2), (0.5, 330.5, 1), (0.5, 89.5, 3), (0.5, 29.5, -100), (40.5, 139.5, 12)):
        assert at(la, lo) == want, f"抽查 ({la},{lo}) = {at(la, lo)} ≠ {want}"
    return [(-89.5 + i // NLON, 0.5 + i % NLON, int(c)) for i, c in enumerate(s)]


def do_woa09(root, man, stats):
    task = {"source": "woa09", "url": WOA09_URL, "rel": WOA09_REL}
    if not man.have(WOA09_REL):
        rec = F.get_whole(task, root, man, stats)
        if rec["bytes"] != WOA09_BYTES:
            print(F.log(f"[WOA09] 体量 {rec['bytes']} ≠ HEAD 所见 {WOA09_BYTES}（仅提示）"), flush=True)
    if man.have(WOA09_CSV):
        return True
    rows = woa09_surface(os.path.join(root, WOA09_REL))
    out = os.path.join(root, WOA09_CSV)
    with open(out + ".tmp", "w", encoding="utf-8") as f:
        f.write("latitude,longitude,basin_code\n")
        for la, lo, c in rows:
            f.write(f"{la:g},{lo:g},{c}\n")
    os.replace(out + ".tmp", out)
    counts = collections.Counter(c for _, _, c in rows)
    man.add(source="woa09_derived", path=WOA09_CSV, bytes=os.path.getsize(out), sha256=F.sha256_file(out), status="ok",
            derived_from={"path": WOA09_REL, "sha256": man.recs[WOA09_REL]["sha256"]},
            layout="level 1 (surface) of basin.msk; lat -89.5..89.5 x lon 0.5..359.5, lon fastest; -100=land",
            code_counts={str(k): v for k, v in sorted(counts.items())})
    print(F.log(f"[WOA09] 表层掩膜 {len(rows)} 格，分区码计数 {dict(sorted(counts.items()))}"), flush=True)
    return True


# ======================================================================== 编排

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", help="p3-cache 根目录（绝对路径）")
    ap.add_argument("--min-free-gb", type=float, default=10.0)
    a = ap.parse_args()
    if not a.cache or not os.path.isabs(a.cache) or not os.path.isdir(a.cache):
        print("需要已存在的绝对路径 --cache", flush=True)
        return 3
    free = shutil.disk_usage(a.cache).free / 1e9
    if free < a.min_free_gb:
        print(f"磁盘空闲 {free:.0f} GB < {a.min_free_gb} GB，停止", flush=True)
        return 3
    out = os.environ.get("REPRO_OUTPUT_DIR")
    if out:
        F.LOGF = open(os.path.join(out, "fetch_extra_log.txt"), "a", encoding="utf-8")
    print(F.log(f"[P3 extra] {VERSION}（p3_fetch {F.VERSION}），空闲 {free:.0f} GB，根 {a.cache}"), flush=True)
    man = extra_manifest(a.cache)
    stats = F.Stats(2 * F.EST["era5_raw"] + WOA09_BYTES)
    stats.total_files = 5
    pool = cf.ThreadPoolExecutor(max_workers=F.N_S3, thread_name_prefix="s3")
    t0 = time.monotonic()
    status, failures = {}, []
    for name, fn in (("woa09", lambda: do_woa09(a.cache, man, stats)),
                     ("era5_200101", lambda: do_era5(a.cache, man, stats, pool))):
        try:
            status[name] = "ok" if fn() else "failed"
        except Exception as e:
            traceback.print_exc()
            status[name] = "failed"
            failures.append((name, f"{type(e).__name__}: {e}"))
            print(F.log(f"失败 {name}：{type(e).__name__}: {e}"), flush=True)
    pool.shutdown(wait=False)
    ok = all(v == "ok" for v in status.values())
    summary = {"version": VERSION, "p3_fetch_version": F.VERSION, "status": status, "failures": failures,
               "seconds": round(time.monotonic() - t0, 1), "bytes_this_run": stats.bytes,
               "files": {k: {x: r.get(x) for x in ("bytes", "sha256", "status")} for k, r in man.recs.items()}}
    if out:
        with open(os.path.join(out, "fetch_extra_summary.json" if ok else "fetch_extra_status.json"), "w",
                  encoding="utf-8") as f:
            json.dump(summary, f, ensure_ascii=False, indent=1)
    print(F.log(f"[P3 extra] {'全部完成' if ok else '未完成'}：{json.dumps(status, ensure_ascii=False)}，"
                f"{stats.bytes / 1e9:.2f} GB，{summary['seconds']} s"), flush=True)
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
