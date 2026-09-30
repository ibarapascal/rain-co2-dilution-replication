#!/usr/bin/env python3
"""raw_archive.py — 原始数据补档：把各阶段用过的全部原始下载文件在 bulk 原始档留完整原件＋文件级索引。

Description：
  原始数据一律保留（含 P5 下载的原件），合计限速 4 MB/s。
  阶段 A 把盘上已有的原件复制到 bulk（只复制、不删源、逐文件 sha256）；阶段 B 按限速重下已被旧脚本删掉的原件。
  共用 I/O（挂载哨兵、令牌桶、manifest、可续传下载）在 `raw_store.py`；CMORPH 路径与 URL 规则与 `p2_rim_test.hour_url`、
  `p3_fetch.CMORPH_BASE` 相同，ERA5 URL/sha256/etag 取自 P3 的 `manifest*.jsonl`。

阶段（按序；每个文件落盘后算 sha256 并追加 `<raw>/manifest.jsonl`；已在 raw 档的文件跳过＝断点续传）：
  A1  P3 CMORPH 2000 年原件（fast p3-cache/cmorph/YYYY/MM）→ raw/cmorph/YYYY/MM/DD/（与 NCEI 层级一致），sha256 对 P3 manifest
  A2  P3 其余原件 hycom_sss、oisst、glodap、watson、woa09 ＋ manifest.jsonl、manifest_extra.jsonl、fetch_done.json、
      runs/p3-fetch*（下载日志与 fetch_summary）→ raw/p3-cache/…；有 P3 manifest 记录的逐个核 sha256
  A3  P1 系泊原始缓存（P1 输出目录下 cache/）→ raw/p1/cache/
  A4  P1b 系泊原始缓存（P1b 输出目录下 cache/，下游 2.5 万个符号链接指向它：只读、只复制）→ raw/p1b/cache/
  A5  P5 已取的 SMAP OPeNDAP 子集（p5-cache/smap，npz＝子集解析后存盘，非字节原样）与 CMR 查询结果（p5-cache/cmr）→ raw/p5-cache/
  B1  P5 取消前已抽窗口、原件被删的 CMORPH 小时（p5-cache/cmorph/<站>/<hour>.npy 的并集）
  B2  P5 用到的 SMAP L2C V6 完整颗粒（p5-cache/smap 下 787 个 granule；受保护，需 Earthdata token）＋公开 .md5，逐个核 md5
  B3  P2 主集与 P2-sens 的 CMORPH 小时（两份 cmorph_pixels.csv 的小时并集；与 A1、P5 已归档者按文件存在去重）
  B4  ERA5 u10/v10 原件（P3 manifest 与 manifest_extra 的 era5_raw 行，去重后 28 个），核 sha256，不符记 sha_mismatch 并保留

限速（所有下载合计 ≤4 MB/s，1 MB＝10^6 B）：本进程一个令牌桶；每 15 s 看 `<raw>/.p5_fetch_heartbeat`——
  mtime 在 180 s 内且内容不是 done ⇒ P5 正在下载（P5 a6 自己限 2 MB/s）⇒ 本任务限 --rate-shared（2）；否则 --rate-full（4）。
  P5 下载段结束会把心跳写成 done；P5 崩溃则心跳 180 s 后过期，本任务自动升到 4。P5 还没开始时本任务先按 4 跑，P5 一开跑
  （首个心跳）15 s 内降到 2。

认证：只有 B2 需要。token 只从配置的 token 文件 `<token-file>` 读进内存（沿用 p5_sss_sat.read_token），
  只放进发往 Earthdata 主机（p5_sss_sat.auth_host 判定）的 Authorization 头；重定向到 S3 等其他主机时不带；不进 argv、env、
  日志、manifest、仓库；异常文本经 redact()。

Usage：
  python raw_archive.py --plan          只统计各阶段待处理文件数与字节（本地列表＋已知均值，不联网）→ <out>/raw_plan.json
  python raw_archive.py --selftest-net  小自测（<100 MB）：写 bulk 下 raw-selftest/（与正式 raw 档分开），见 selftest_net()
  python raw_archive.py --run           正式补档（A1–A5、B1–B4）；--stages 可只跑部分，如 --stages B2,B4
  --raw-root DIR（默认 raw_store.RAW_ROOT_DEFAULT）  --out DIR（日志与状态；默认 $REPRO_OUTPUT_DIR）
Dependencies：Python 环境（标准库；自测里用 netCDF4 打开一个 CMORPH 文件核 time 轴）；同目录 raw_store.py、p5_sss_sat.py（只取认证函数）。
退出码：0 全部完成；2 有文件多次失败（重跑即续传）；3 未预期异常／挂载哨兵不过；5 Earthdata 认证失败（其余阶段照跑）。

Change Log：
  2026-09-27 初版。
"""

import argparse
import concurrent.futures as cf
import csv
import glob
import json
import os
import re
import shutil
import sys
import threading
import time
import traceback

import raw_store as rs

import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
# [repro] 以下路径来自集中配置（各阶段输出目录与 token 文件）
PAPER = _rp.path("fast_root")
P3_CACHE = _rp.path("p3_cache")
P5_CACHE = _rp.path("p5_cache")
OUT_ROOT = _rp.path("out_root")
P1_CACHE = _rp.upstream("p1_cache")
P1B_CACHE = _rp.upstream("p1b_cache")
P2_PIX = _rp.upstream("p2_pixels")
P2S_PIX = _rp.upstream("p2sens_pixels")
TOKEN_FILE = _rp.path("token_file")  # 只取路径，内容只进 Earthdata 主机请求头
SMAP_PROT = "https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/SMAP_RSS_L2_SSS_V6/"
SMAP_PUB = "https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-public/SMAP_RSS_L2_SSS_V6/"
SELFTEST_ROOT = _rp.path("raw_selftest_root")   # 每次自测一个时间戳子目录  [repro] 路径来自集中配置
STAGES = ("A1", "A2", "A3", "A4", "A5", "B1", "B2", "B3", "B4")
EST = {"cmorph": 1.55e6, "smap": 66.3e6}          # 均值：CMORPH 取 P5 实测 4.51 GB/3000 与 P2 HEAD 1.64 MB 之间；SMAP 取 CMR 样本
CONC = {"cmorph": 4, "smap": 2, "era5": 1}
MAX_CONSEC_FAIL = 40
VERSION = "raw-archive-2026-09-27a"


class Log:
    def __init__(self, path):
        self.f = open(path, "a", encoding="utf-8")
        self.lock = threading.Lock()

    def log(self, msg, echo=True):
        line = f"{rs.now_iso()} {msg}"
        with self.lock:
            self.f.write(line + "\n")
            self.f.flush()
        if echo:
            print(line, flush=True)


# ================================================================ 清单
def p3_manifest_rows():
    rows = []
    for fn in ("manifest.jsonl", "manifest_extra.jsonl"):
        p = os.path.join(P3_CACHE, fn)
        if os.path.exists(p):
            with open(p, encoding="utf-8") as f:
                rows += [json.loads(ln) for ln in f if ln.strip()]
    return rows


def p3_sha_by_path(rows):
    return {r["path"]: r["sha256"] for r in rows if r.get("status") == "ok" and r.get("sha256") and r.get("path")}


def p3_url_by_path(rows):
    return {r["path"]: r["url"] for r in rows if r.get("status") == "ok" and r.get("url") and r.get("path")}


def era5_list(rows):
    seen = {}
    for r in rows:
        if r.get("source") == "era5_raw" and r.get("status") == "ok" and r.get("url"):
            seen.setdefault(r["url"], r)
    return sorted(seen.values(), key=lambda r: r["url"])


def era5_rel(url):
    return "era5/" + url.split(".amazonaws.com/", 1)[1]


def p5_cached_hours():
    hs = set()
    for d in glob.glob(os.path.join(P5_CACHE, "cmorph", "*", "")):
        for fn in os.listdir(d):
            if fn.endswith(".npy") and not fn.endswith(".part.npy"):
                try:
                    hs.add(int(fn[:-4]))
                except ValueError:
                    pass
    return hs


def smap_granules():
    g = set()
    for d in glob.glob(os.path.join(P5_CACHE, "smap", "*", "")):
        for fn in os.listdir(d):
            m = re.match(r"^(RSS_SMAP_SSS_L2C_.+_V06\.0)\.npz(\.missing|\.nocover)?$", fn)   # P5 缓存：<gid>.npz[.missing|.nocover]
            if m:
                g.add(m.group(1))
    return sorted(g)


def smap_rel(gid):
    m = re.search(r"_(\d{4})\d{4}T\d{6}_", gid)
    return f"smap_l2c_v6/{m.group(1) if m else 'unknown'}/{gid}.nc"


def p2_hours():
    hs = set()
    for p in (P2_PIX, P2S_PIX):
        with open(p, encoding="utf-8") as f:
            for row in csv.reader(f):
                if row and row[0] != "hour":
                    try:
                        hs.add(int(row[0]))
                    except ValueError:
                        pass
    return hs


def walk_files(src):
    out, links = [], 0
    for dp, dns, fns in os.walk(src, followlinks=False):
        dns.sort()
        for fn in sorted(fns):
            p = os.path.join(dp, fn)
            if os.path.islink(p):
                links += 1
                continue
            out.append(p)
    return out, links


# ================================================================ 阶段 A：本地复制
class Copier:
    def __init__(self, root, man, log, tag):
        self.root, self.man, self.log, self.tag = root, man, log, tag

    def copy(self, src, rel, dataset, expect_sha=None, url=None):
        dst = os.path.join(self.root, rel)
        if os.path.exists(dst) and os.path.getsize(dst) == os.path.getsize(src):
            return "exists", os.path.getsize(src)
        import repro_io as _rio  # [repro] raw archive is append-only: an existing file of different size is not replaced
        _rio.refuse_overwrite(dst)  # [repro] (os.replace below is reached only for new files)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        part = os.path.join(self.root, ".tmp", rel.replace("/", "__") + ".copy.part")
        shutil.copyfile(src, part)
        shutil.copystat(src, part)
        size = os.path.getsize(part)
        if size != os.path.getsize(src):
            raise IOError(f"复制后长度不符：{src}")
        sha = rs.sha256_file(part)
        status = "copied"
        if expect_sha and sha != expect_sha:
            status = "copied_sha_mismatch_vs_source_manifest"
        os.replace(part, dst)
        rec = dict(dataset=dataset, url=url, path=rel, bytes=size, sha256=sha, status=status, src=src,
                   source_job=self.tag)
        if expect_sha:
            rec["expect_sha256"] = expect_sha
        self.man.add(**rec)
        return status, size

    def tree(self, src, rel_prefix, dataset, sha_by_rel=None, sha_prefix="", url_by_rel=None):
        files, links = walk_files(src)
        stat = {"src": src, "dst": os.path.join(self.root, rel_prefix), "files": len(files), "symlinks_skipped": links,
                "bytes": 0, "copied": 0, "exists": 0, "sha_mismatch": 0}
        for p in files:
            r = os.path.relpath(p, src)
            exp = (sha_by_rel or {}).get(sha_prefix + r)
            st, size = self.copy(p, f"{rel_prefix}/{r}", dataset, expect_sha=exp, url=(url_by_rel or {}).get(sha_prefix + r))
            stat["bytes"] += size
            stat["copied" if st.startswith("copied") else "exists"] += 1
            if st == "copied_sha_mismatch_vs_source_manifest":
                stat["sha_mismatch"] += 1
        dfiles, _ = walk_files(os.path.join(self.root, rel_prefix))
        stat["dst_files"] = len(dfiles)
        stat["dst_bytes"] = sum(os.path.getsize(x) for x in dfiles)
        stat["verified_count_bytes"] = (stat["dst_files"] == stat["files"] and stat["dst_bytes"] == stat["bytes"])
        self.log.log(f"[A] {dataset} {src} → {rel_prefix}：{json.dumps(stat, ensure_ascii=False)}")
        return stat


def stage_a1(cp, rows):
    cm = [r for r in rows if r.get("source") == "cmorph" and r.get("status") == "ok"]
    st = {"files": len(cm), "bytes": 0, "copied": 0, "exists": 0, "sha_mismatch": 0, "missing_src": 0}
    for r in cm:
        src = os.path.join(P3_CACHE, r["path"])
        if not os.path.exists(src):
            st["missing_src"] += 1
            continue
        s, size = cp.copy(src, rs.cmorph_rel_from_name(os.path.basename(r["path"])), "cmorph", r.get("sha256"), r.get("url"))
        st["bytes"] += size
        st["copied" if s.startswith("copied") else "exists"] += 1
        st["sha_mismatch"] += s == "copied_sha_mismatch_vs_source_manifest"
    st["verified"] = st["missing_src"] == 0 and st["sha_mismatch"] == 0 and st["copied"] + st["exists"] == st["files"]
    cp.log.log(f"[A1] P3 CMORPH：{json.dumps(st, ensure_ascii=False)}")
    return st


def stage_a2(cp, rows):
    sha, url = p3_sha_by_path(rows), p3_url_by_path(rows)
    out = {}
    for sub in ("hycom_sss", "oisst", "glodap", "watson", "woa09"):
        out[sub] = cp.tree(os.path.join(P3_CACHE, sub), f"p3-cache/{sub}", f"p3_{sub}", sha, sha_prefix=f"{sub}/",
                           url_by_rel=url)
    for fn in ("manifest.jsonl", "manifest_extra.jsonl", "fetch_done.json"):
        s, size = cp.copy(os.path.join(P3_CACHE, fn), f"p3-cache/{fn}", "p3_meta")
        out[fn] = {"status": s, "bytes": size}
    for sub in ("p3-fetch", "p3-fetch-extra"):
        d = os.path.join(OUT_ROOT, sub)  # [repro] 各阶段输出在 out_root/<阶段>
        if os.path.isdir(d):
            out["runs_" + sub] = cp.tree(d, f"p3-cache/_fetch_runs/{sub}", "p3_fetch_run_log")
    return out


# ================================================================ 阶段 B：网络补下
class RateCtl(threading.Thread):
    def __init__(self, bucket, root, shared, full, log):
        super().__init__(daemon=True)
        self.bucket, self.root, self.shared, self.full, self.log = bucket, root, shared, full, log
        self.cur = None
        self.stop = threading.Event()
        self.switches = []

    def tick(self):
        want = self.shared if rs.p5_active(self.root) else self.full
        if want != self.cur:
            self.bucket.set_rate(want * rs.MB)
            self.switches.append([rs.now_iso(), want])
            self.log.log(f"[rate] 限速 → {want} MB/s（P5 心跳{'在' if want == self.shared else '不在'}）")
            self.cur = want

    def run(self):
        while not self.stop.is_set():
            try:
                self.tick()
            except Exception as e:
                self.log.log(f"[rate] 检查心跳出错（保持当前限速）：{e}")
            self.stop.wait(15)


def run_pool(name, items, fn, conc, log, status, key):
    """items：任务列表；fn(item) → (kind, info)。kind ∈ ok/exists/missing/mismatch。多次失败的收集后末轮再试一次。"""
    st = status.setdefault(key, {"todo": len(items), "ok": 0, "exists": 0, "missing": 0, "mismatch": 0, "failed": 0,
                                 "bytes": 0})
    failed, consec, t0 = [], [0], time.monotonic()
    lock = threading.Lock()

    def one(it):
        try:
            kind, size = fn(it)
            with lock:
                st[kind] += 1
                st["bytes"] += size or 0
                consec[0] = 0
        except rs.AuthFail:
            raise
        except Exception as e:
            with lock:
                failed.append(it)
                consec[0] += 1
            log.log(f"[{key}] 失败 {it if not isinstance(it, dict) else it.get('url')}：{e}", echo=False)
        n = st["ok"] + st["exists"] + st["missing"] + st["mismatch"]
        if n and n % 250 == 0:
            el = time.monotonic() - t0
            log.log(f"[{key}] {n}/{len(items)}，本段 {st['bytes'] / 1e9:.2f} GB，{st['bytes'] / 1e6 / max(el, 1):.2f} MB/s")

    with cf.ThreadPoolExecutor(max_workers=conc) as ex:
        it = iter(items)
        pend = set()
        for x in it:
            pend.add(ex.submit(one, x))
            if len(pend) >= conc * 2:
                break
        while pend:
            done, pend = cf.wait(pend, return_when=cf.FIRST_COMPLETED)
            for f in done:
                f.result()          # AuthFail 在此上抛
            if consec[0] >= MAX_CONSEC_FAIL:
                for f in pend:
                    f.cancel()
                log.log(f"[{key}] 连续 {consec[0]} 个失败，本段中止（重跑即续传）")
                st["aborted"] = True
                break
            for x in it:
                pend.add(ex.submit(one, x))
                if len(pend) >= conc * 2:
                    break
    still = []
    if not st.get("aborted"):
        for x in failed:
            try:
                kind, size = fn(x)
                st[kind] += 1
                st["bytes"] += size or 0
            except rs.AuthFail:
                raise
            except Exception as e:
                still.append(x)
                log.log(f"[{key}] 末轮仍失败：{e}", echo=False)
    else:
        still = failed
    st["failed"] = len(still)
    st["seconds"] = round(time.monotonic() - t0, 1)
    log.log(f"[{key}] 完成：{json.dumps(st, ensure_ascii=False)}")
    return st


def make_cmorph_fn(dl):
    def fn(h):
        kind, val, rec = dl.get(rs.cmorph_url(h), rs.cmorph_rel(h), "cmorph")
        if kind == "exists":
            return "exists", 0
        if kind == "missing":
            return "missing", 0
        return "ok", rec["bytes"]
    return fn


def make_smap_fn(dl):
    def fn(gid):
        rel = smap_rel(gid)
        md5_rel = rel + ".md5"
        k0, md5_path, _ = dl.get(SMAP_PUB + gid + ".nc.md5", md5_rel, "smap_l2c_v6_md5")
        exp = None
        if k0 in ("ok", "exists"):
            with open(md5_path, encoding="utf-8", errors="replace") as f:
                tok = f.read().split()
            exp = tok[0].lower() if tok and re.fullmatch(r"[0-9a-fA-F]{32}", tok[0]) else None
        kind, val, rec = dl.get(SMAP_PROT + gid + ".nc", rel, "smap_l2c_v6", expect_md5=exp,
                                extra={"granule": gid, "collection": "C2832221740-POCLOUD"})
        if kind == "exists":
            return "exists", 0
        if kind == "missing":
            return "missing", 0
        return ("ok" if rec["status"] == "ok" else "mismatch"), rec["bytes"]
    return fn


def make_era5_fn(dl):
    def fn(r):
        kind, val, rec = dl.get(r["url"], era5_rel(r["url"]), "era5_raw", expect_bytes=r.get("bytes"),
                                expect_sha256=r.get("sha256"),
                                extra={"expect_etag": r.get("etag"), "p3_manifest_t": r.get("t")})
        if kind == "exists":
            return "exists", 0
        if kind == "missing":
            return "missing", 0
        return ("ok" if rec["status"] == "ok" else "mismatch"), rec["bytes"]
    return fn


def auth_bits():
    import p5_sss_sat as p5
    return p5.read_token, p5.auth_host, p5.redact, p5.AuthError


# ================================================================ plan / run / selftest
def plan(root):
    rows = p3_manifest_rows()
    have = lambda rel: os.path.exists(os.path.join(root, rel))
    a1 = [r for r in rows if r.get("source") == "cmorph" and r.get("status") == "ok"]
    p5h = p5_cached_hours()
    gids = smap_granules()
    p2h = p2_hours()
    p3h_rels = {rs.cmorph_rel_from_name(os.path.basename(r["path"])) for r in a1}
    b3 = sorted(h for h in p2h if h not in p5h and rs.cmorph_rel(h) not in p3h_rels)
    e5 = era5_list(rows)

    def tree_size(d):
        fs, links = walk_files(d)
        return {"files": len(fs), "bytes": sum(os.path.getsize(x) for x in fs), "symlinks": links}

    out = {"version": VERSION, "raw_root": root, "t": rs.now_iso(),
           "A1_p3_cmorph": {"files": len(a1), "bytes": sum(r["bytes"] for r in a1),
                            "already": sum(have(rs.cmorph_rel_from_name(os.path.basename(r["path"]))) for r in a1)},
           "A2_p3_other": {s: tree_size(os.path.join(P3_CACHE, s)) for s in ("hycom_sss", "oisst", "glodap", "watson", "woa09")},
           "A3_p1_cache": tree_size(P1_CACHE), "A4_p1b_cache": tree_size(P1B_CACHE),
           "A5_p5_cache": {s: tree_size(os.path.join(P5_CACHE, s)) for s in ("smap", "cmr")},
           "B1_p5_deleted_cmorph": {"hours": len(p5h), "already": sum(have(rs.cmorph_rel(h)) for h in p5h),
                                    "est_bytes": int(len(p5h) * EST["cmorph"])},
           "B2_smap_granules": {"granules": len(gids), "already": sum(have(smap_rel(g)) for g in gids),
                                "est_bytes": int(len(gids) * EST["smap"])},
           "B3_p2_cmorph": {"p2_union_hours": len(p2h), "overlap_p5_cached": len(p2h & p5h),
                            "to_fetch_upper": len(b3), "est_bytes_upper": int(len(b3) * EST["cmorph"]),
                            "note": "P5 余下约 8.2k 小时中与 P2 重叠的部分由 P5 先落档，届时按存在跳过，实下会更少"},
           "B4_era5": {"files": len(e5), "bytes": sum(r["bytes"] for r in e5), "already": sum(have(era5_rel(r["url"])) for r in e5)}}
    net = out["B1_p5_deleted_cmorph"]["est_bytes"] + out["B2_smap_granules"]["est_bytes"] + \
        out["B3_p2_cmorph"]["est_bytes_upper"] + out["B4_era5"]["bytes"]
    out["network_bytes_upper"] = net
    out["hours_at_3MBps"] = round(net / 3e6 / 3600, 1)
    return out


def run(args, out, log):
    root = rs.ensure_root(args.raw_root)
    man = rs.Manifest(root)
    tag = args.tag
    stages = [s.strip() for s in args.stages.split(",")] if args.stages else list(STAGES)
    status = {"version": VERSION, "start": rs.now_iso(), "raw_root": root, "stages": stages}
    spath = [os.path.join(out, "raw_archive_status.json"), os.path.join(root, "archive_status.json")]

    def dump():
        status["updated"] = rs.now_iso()
        for p in spath:
            tmp = p + ".part"
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(status, f, ensure_ascii=False, indent=1)
            os.replace(tmp, p)

    status["plan"] = plan(root)
    dump()
    rows = p3_manifest_rows()
    cp = Copier(root, man, log, tag)
    if "A1" in stages:
        status["A1"] = stage_a1(cp, rows); dump()
    if "A2" in stages:
        status["A2"] = stage_a2(cp, rows); dump()
    if "A3" in stages:
        status["A3"] = cp.tree(P1_CACHE, "p1/cache", "p1_mooring_cache"); dump()
    if "A4" in stages:
        status["A4"] = cp.tree(P1B_CACHE, "p1b/cache", "p1b_mooring_cache"); dump()
    if "A5" in stages:
        status["A5"] = {s: cp.tree(os.path.join(P5_CACHE, s), f"p5-cache/{s}",
                                   "p5_smap_subset_npz" if s == "smap" else "p5_cmr_json") for s in ("smap", "cmr")}
        dump()

    bucket = rs.TokenBucket(args.rate_full * rs.MB)
    ctl = RateCtl(bucket, root, args.rate_shared, args.rate_full, log)
    ctl.tick()
    ctl.start()
    status["rate"] = {"shared_mbps": args.rate_shared, "full_mbps": args.rate_full, "switches": ctl.switches}
    rc = 0
    ticker_stop = threading.Event()

    def ticker():
        while not ticker_stop.wait(60):
            try:
                dump()
            except Exception:
                pass

    threading.Thread(target=ticker, daemon=True).start()
    dl_anon = rs.Downloader(root, bucket, man, "arc", tag, log=log)
    dls = [dl_anon]
    try:
        if "B1" in stages:
            hs = sorted(p5_cached_hours())
            log.log(f"[B1] P5 已抽窗口小时 {len(hs)}")
            run_pool("B1", hs, make_cmorph_fn(dl_anon), CONC["cmorph"], log, status, "B1"); dump()
        if "B2" in stages:
            read_token, auth_host, redact, AuthError = auth_bits()
            try:
                token = read_token(TOKEN_FILE)
                dl_auth = rs.Downloader(root, bucket, man, "arc", tag, token=token, auth_host=auth_host,
                                        redact=redact, log=log)
                dls.append(dl_auth)
                gids = smap_granules()
                log.log(f"[B2] SMAP 颗粒 {len(gids)}")
                run_pool("B2", gids, make_smap_fn(dl_auth), CONC["smap"], log, status, "B2")
            except (rs.AuthFail, AuthError) as e:
                log.log(f"[B2] 认证失败，跳过 SMAP（其余阶段照跑）：{redact(e)}")
                status["B2_auth_error"] = redact(str(e))
                rc = 5
            dump()
        if "B3" in stages:
            hs = sorted(p2_hours())
            log.log(f"[B3] P2＋P2-sens 小时并集 {len(hs)}（已在 raw 档者跳过）")
            run_pool("B3", hs, make_cmorph_fn(dl_anon), CONC["cmorph"], log, status, "B3"); dump()
        if "B4" in stages:
            e5 = era5_list(rows)
            log.log(f"[B4] ERA5 原件 {len(e5)}")
            run_pool("B4", e5, make_era5_fn(dl_anon), CONC["era5"], log, status, "B4"); dump()
    finally:
        ctl.stop.set()
        ticker_stop.set()
        status["bytes_downloaded_total"] = sum(d.bytes for d in dls)
        status["end"] = rs.now_iso()
        dump()
    bad = any(isinstance(status.get(k), dict) and (status[k].get("failed") or status[k].get("aborted"))
              for k in ("B1", "B2", "B3", "B4"))
    if bad and rc == 0:
        rc = 2
    return rc


def selftest_net(args, out, log):
    """小自测（<100 MB，写 SELFTEST_ROOT）：挂载哨兵、令牌桶、本地复制核 sha、CMORPH 下载＋复用＋time 轴、
    SMAP 颗粒断线续传＋md5、ERA5 HEAD 对 P3 manifest、心跳判定。"""
    res = []

    def chk(name, ok, detail=""):
        res.append({"name": name, "ok": bool(ok), "detail": detail})
        log.log(f"  {'ok ' if ok else 'BAD'}  {name}  {detail}")

    root = rs.ensure_root(os.path.join(SELFTEST_ROOT, time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())))
    chk("挂载哨兵：bulk 为挂载点", os.path.ismount(_rp.volume_of(args.raw_root)))
    try:
        rs.ensure_root("/Volumes/__no_such_volume__/x")
        chk("挂载哨兵：未挂载卷拒写", False)
    except RuntimeError:
        chk("挂载哨兵：未挂载卷拒写", True)
    man = rs.Manifest(root)
    tb = rs.TokenBucket(4 * rs.MB)
    t0 = time.monotonic()
    for _ in range(32):
        tb.consume(rs.CHUNK)
    el = time.monotonic() - t0
    chk("令牌桶 4 MB/s：2.1 MB 约 0.5 s", 0.35 < el < 0.9, f"{el:.2f}s")
    # 本地复制
    rows = p3_manifest_rows()
    cp = Copier(root, man, log, "selftest")
    cm = [r for r in rows if r.get("source") == "cmorph" and r.get("status") == "ok"][:3]
    ok = True
    for r in cm:
        s, _ = cp.copy(os.path.join(P3_CACHE, r["path"]), rs.cmorph_rel_from_name(os.path.basename(r["path"])),
                       "cmorph", r["sha256"], r["url"])
        ok &= s == "copied"
    chk("A1 复制 3 个 P3 CMORPH 并核 sha256", ok)
    # CMORPH 下载
    bucket = rs.TokenBucket(2 * rs.MB)
    dl = rs.Downloader(root, bucket, man, "arc", "selftest", log=log)
    hs = sorted(p2_hours())[:3]
    t0 = time.monotonic()
    got = [dl.get(rs.cmorph_url(h), rs.cmorph_rel(h), "cmorph") for h in hs]
    el = time.monotonic() - t0
    rate = dl.bytes / max(el, 1e-6) / 1e6
    chk("CMORPH 3 小时下载", all(k == "ok" for k, _, _ in got), f"{dl.bytes} B，{rate:.2f} MB/s")
    chk("CMORPH 限速 ≤2.2 MB/s", rate <= 2.2, f"{rate:.2f}")
    again = [dl.get(rs.cmorph_url(h), rs.cmorph_rel(h), "cmorph")[0] for h in hs]
    chk("CMORPH 已归档即复用不重下", all(k == "exists" for k in again))
    cr = rs.CmorphRaw(root, 2.0, owner="p5", source_job="selftest", bucket=bucket, manifest=man)
    k, pth = cr.fetch(hs[0])
    chk("CmorphRaw（P5 用）命中 raw 档", k == "ok" and cr.reused == 1 and cr.downloaded == 0)
    try:
        import netCDF4
        with netCDF4.Dataset(pth) as ds:
            t = [int(x) for x in ds["time"][:]]
        chk("原件 time 轴与小时一致", t == [hs[0] * 3600, hs[0] * 3600 + 1800], str(t))
    except Exception as e:
        chk("原件 time 轴与小时一致", False, repr(e))
    # SMAP 断线续传＋md5
    read_token, auth_host, redact, AuthError = auth_bits()
    try:
        token = read_token(TOKEN_FILE)
        dla = rs.Downloader(root, rs.TokenBucket(4 * rs.MB), man, "arc", "selftest", token=token,
                            auth_host=auth_host, redact=redact, log=log)
        gid = smap_granules()[0]
        dla.abort_after = 20 * 1024 * 1024
        try:
            make_smap_fn(dla)(gid)
            chk("SMAP 模拟断线", False, "没有触发")
        except Exception as e:
            dla._drop()
            part = dla.part_path(smap_rel(gid))
            chk("SMAP 模拟断线留下 part", os.path.exists(part) and os.path.getsize(part) >= 20 * 1024 * 1024,
                f"{type(e).__name__} part={os.path.getsize(part) if os.path.exists(part) else None}")
        b0, t0 = dla.bytes, time.monotonic()
        kind, size = make_smap_fn(dla)(gid)
        r4 = (dla.bytes - b0) / max(time.monotonic() - t0, 1e-6) / 1e6
        chk("限速 4 MB/s 下实测不超（续传段）", r4 <= 4.3, f"{r4:.2f} MB/s（低于 4 属网络或服务器限制，只作参考）")
        recs = [r for r in man.load() if r.get("dataset") == "smap_l2c_v6"]
        chk("SMAP Range 续传完成且 md5 一致", kind == "ok" and recs and recs[-1].get("resumed") and
            recs[-1].get("md5") == recs[-1].get("expect_md5"), f"{size} B")
        chk("manifest 不含 token", token not in open(man.path, encoding="utf-8").read())
    except AuthError as e:
        chk("SMAP 认证", False, redact(e))
    # ERA5 HEAD
    e5 = era5_list(rows)
    import http.client
    import urllib.parse
    okn = 0
    for r in e5[:2]:
        u = urllib.parse.urlsplit(r["url"])
        c = http.client.HTTPSConnection(u.hostname, timeout=60)
        c.request("HEAD", u.path, headers={"User-Agent": rs.UA})
        h = c.getresponse()
        h.read()
        okn += int(h.status == 200 and int(h.headers.get("Content-Length", -1)) == r["bytes"]
                   and h.headers.get("ETag", "").strip('"') == r.get("etag"))
        c.close()
    chk("ERA5 HEAD 长度与 etag 对 P3 manifest", okn == 2, f"{okn}/2（共 {len(e5)} 个原件）")
    # 心跳
    rs.heartbeat_write(root, "active")
    a = rs.p5_active(root)
    rs.heartbeat_write(root, "done")
    b = rs.p5_active(root)
    chk("心跳 active→限 2，done→升 4", a and not b)
    tot = sum(r.get("bytes") or 0 for r in man.load() if r.get("status") in ("ok", "copied"))
    chk("自测总字节 <100 MB", tot < 100e6, f"{tot / 1e6:.1f} MB")
    summary = {"version": VERSION, "root": root, "passed": sum(r["ok"] for r in res), "total": len(res),
               "all_ok": all(r["ok"] for r in res), "checks": res, "plan": plan(args.raw_root)}
    with open(os.path.join(out, "raw_selftest.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    log.log(f"[selftest] {summary['passed']}/{summary['total']}")
    return 0 if summary["all_ok"] else 4


def main(argv=None):
    ap = argparse.ArgumentParser(description="原始数据补档（bulk）")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--plan", action="store_true")
    g.add_argument("--selftest-net", action="store_true")
    g.add_argument("--run", action="store_true")
    ap.add_argument("--raw-root", default=rs.RAW_ROOT_DEFAULT)
    ap.add_argument("--out")
    ap.add_argument("--stages", default="")
    ap.add_argument("--rate-full", type=float, default=4.0)
    ap.add_argument("--rate-shared", type=float, default=2.0)
    ap.add_argument("--tag", default="raw-archive")
    args = ap.parse_args(argv)
    if args.rate_full > 4.0 or args.rate_shared > args.rate_full:
        print("限速参数越界（总上限 4 MB/s）", file=sys.stderr)
        return 3
    out = args.out or os.environ.get("REPRO_OUTPUT_DIR")
    if not out:
        print("需要 --out 或 REPRO_OUTPUT_DIR", file=sys.stderr)
        return 3
    os.makedirs(out, exist_ok=True)
    log = Log(os.path.join(out, "raw_archive_log.txt"))
    log.log(f"=== start {VERSION} mode={'plan' if args.plan else 'selftest' if args.selftest_net else 'run'} "
            f"raw={args.raw_root} out={out}")
    try:
        if args.plan:
            p = plan(args.raw_root)
            with open(os.path.join(out, "raw_plan.json"), "w", encoding="utf-8") as f:
                json.dump(p, f, ensure_ascii=False, indent=1)
            log.log(json.dumps(p, ensure_ascii=False))
            rc = 0
        elif args.selftest_net:
            rc = selftest_net(args, out, log)
        else:
            rc = run(args, out, log)
    except Exception:
        log.log("FATAL 未预期异常：\n" + traceback.format_exc())
        rc = 3
    log.log(f"=== done rc={rc}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
