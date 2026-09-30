#!/usr/bin/env python3
"""raw_sync.py — 把各阶段自带下载器落在缓存／输出目录里的原件复制进 bulk 原始档（只复制、不删源、逐文件 sha256）。

为什么需要：P1／P1b（系泊缓存）与 P3 下载器（p3_fetch／p3_fetch_extra）有自己的落盘布局，不经 raw_store；整合版已去掉它们的
「派生后删原件」，但原件仍在缓存目录里。本脚本在这些阶段之后跑一次，把原件按 raw_store 目录约定登记进 `<raw_root>/manifest.jsonl`。
P2／P5／P7a／raw_fetch_list／p7e_fetch 的下载本来就直接落 raw 档，不需要本脚本。

阶段（复用 raw_archive.py 的 Copier 与 stage_a1／stage_a2，一行分析逻辑都没有）：
  A1  <p3_cache>/cmorph → raw/cmorph/YYYY/MM/DD/（核 P3 manifest sha256）
  A2  <p3_cache>/{hycom_sss,oisst,glodap,watson,woa09}、manifest*.jsonl、fetch_done.json、<out_root>/p3-fetch* 日志 → raw/p3-cache/…
  A3  P1 缓存（upstream p1_cache）→ raw/p1/cache/
  A4  P1b 缓存（upstream p1b_cache）→ raw/p1b/cache/
  A5  P5 缓存 smap／cmr（npz 为子集解析后存盘，非字节原样；原样子集由 p5_sss_sat 另存 raw/smap_opendap_subset/）→ raw/p5-cache/
  A6  （新增）<p3_cache>/era5_raw 下 ERA5 u10/v10 原件（P3 manifest 与 manifest_extra 的 era5_raw ok 行）→ raw/era5/<S3 路径>，
      核 P3 manifest sha256；整合版 p3_fetch 不再删这些原件（历史版删了，历史原件由 raw_archive.py B4 重下补档）。

用法：python raw_sync.py --stages A1,A2,A3,A4,A6 [--raw-root DIR] [--out DIR] [--tag T]
退出码：0 完成且无 sha 不符；4 有 sha 不符或源缺失（文件照样复制保留，记 manifest）；3 未预期异常／挂载哨兵不过。
Change Log：
  2026-09-27 初版。
"""

import argparse
import json
import os
import sys
import traceback

import raw_archive as ra
import raw_store as rs
import repro_paths as _rp

STAGES = ("A1", "A2", "A3", "A4", "A5", "A6")


def stage_a6(cp, rows):
    st = {"files": 0, "bytes": 0, "copied": 0, "exists": 0, "sha_mismatch": 0, "missing_src": 0}
    for r in ra.era5_list(rows):
        st["files"] += 1
        src = os.path.join(ra.P3_CACHE, r["path"])
        if not os.path.exists(src):
            st["missing_src"] += 1
            continue
        s, size = cp.copy(src, ra.era5_rel(r["url"]), "era5_raw", r.get("sha256"), r.get("url"))
        st["bytes"] += size
        st["copied" if s.startswith("copied") else "exists"] += 1
        st["sha_mismatch"] += s == "copied_sha_mismatch_vs_source_manifest"
    st["verified"] = st["missing_src"] == 0 and st["sha_mismatch"] == 0
    cp.log.log(f"[A6] ERA5 u/v 原件：{json.dumps(st, ensure_ascii=False)}")
    return st


def main(argv=None):
    ap = argparse.ArgumentParser(description="原件同步进 bulk raw 档")
    ap.add_argument("--stages", default="A1,A2,A3,A4,A6")
    ap.add_argument("--raw-root", default=rs.RAW_ROOT_DEFAULT)
    ap.add_argument("--out")
    ap.add_argument("--tag", default=_rp.run_tag("raw-sync"))
    a = ap.parse_args(argv)
    out = a.out or os.environ.get("REPRO_OUTPUT_DIR")
    if not out:
        print("需要 --out 或 REPRO_OUTPUT_DIR", file=sys.stderr)
        return 3
    os.makedirs(out, exist_ok=True)
    log = ra.Log(os.path.join(out, "raw_sync_log.txt"))
    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    bad_st = [s for s in stages if s not in STAGES]
    if bad_st:
        print(f"未知阶段 {bad_st}", file=sys.stderr)
        return 3
    status = {"stages": stages, "raw_root": a.raw_root, "start": rs.now_iso()}
    try:
        root = rs.ensure_root(a.raw_root)
        cp = ra.Copier(root, rs.Manifest(root), log, a.tag)
        rows = ra.p3_manifest_rows()
        if "A1" in stages:
            status["A1"] = ra.stage_a1(cp, rows)
        if "A2" in stages:
            status["A2"] = ra.stage_a2(cp, rows)
        if "A3" in stages:
            status["A3"] = cp.tree(ra.P1_CACHE, "p1/cache", "p1_mooring_cache")
        if "A4" in stages:
            status["A4"] = cp.tree(ra.P1B_CACHE, "p1b/cache", "p1b_mooring_cache")
        if "A5" in stages:
            status["A5"] = {s: cp.tree(os.path.join(ra.P5_CACHE, s), f"p5-cache/{s}",
                                       "p5_smap_subset_npz" if s == "smap" else "p5_cmr_json") for s in ("smap", "cmr")}
        if "A6" in stages:
            status["A6"] = stage_a6(cp, rows)
        rc = 0

        def bad(x):
            if isinstance(x, dict):
                if x.get("sha_mismatch") or x.get("missing_src") or x.get("verified") is False \
                        or x.get("verified_count_bytes") is False:
                    return True
                return any(bad(v) for v in x.values() if isinstance(v, dict))
            return False
        if any(bad(status.get(s)) for s in stages):
            rc = 4
    except Exception:
        log.log("FATAL 未预期异常：\n" + traceback.format_exc())
        rc = 3
    status["end"] = rs.now_iso()
    status["rc"] = rc
    with open(os.path.join(out, "raw_sync_status.json"), "w", encoding="utf-8") as f:
        json.dump(status, f, ensure_ascii=False, indent=1)
    log.log(f"=== done rc={rc}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
