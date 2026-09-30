#!/usr/bin/env python3
"""raw_fetch_list.py — 按清单把 PO.DAAC 原始颗粒下载进 bulk 原始档（只下载，不分析）；首用：SPURS-2 0–1 m 剖面相关 8 个 collection。

Description：
  原始数据一律保留；0–1 m 剖面检验的数据先下好，再分析。
  清单由 `--cmr-list` 从 CMR（公开元数据，无 token、无个人信息）生成，存 `pipeline/raw_archive_lists/<name>.json`；
  `--run` 按清单下载到 `<raw>/<rel>`（rel＝`spurs2/<ShortName>/<PO.DAAC 原文件名>`），逐文件核 CMR 给的字节数与 MD5，
  同时下公开 `.md5` 旁注文件并记它与 CMR 值是否一致；每个落盘文件追加 `<raw>/manifest.jsonl`（dataset＝`spurs2`）。
  共用 I/O 全部复用 `raw_store.py`（挂载哨兵、令牌桶、manifest、可续传下载），线程池与日志复用 `raw_archive.py`（run_pool、Log、auth_bits）；
  两个文件一行未改（补档与 P5 下载共用它们）。

排队与限速（所有下载合计 ≤4 MB/s）：本任务接替 P5 a6 的 2 MB/s 名额。
  1) 等待：每 15 s 看 `<raw>/.p5_fetch_heartbeat`，直到 `raw_store.p5_active()` 为假（内容 done，或 180 s 无更新）；
     等待总时限 `--wait-max-s`（默认 6 h），超时退出码 6、不下载。
  2) 接替：写心跳 `active <时刻> pid=<本进程>`（raw_store.heartbeat_write 原格式），再等 `--settle-s`（默认 20 s，大于补档任务
     15 s 的查心跳周期）让补档任务降回 2 MB/s，然后本任务以 `--rate`（默认 2 MB/s）下载；下载期间每 20 s 重写心跳。
  3) 让位：重写心跳前先读文件，若是别的进程（pid≠本进程）在 60 s 内写的 active（例如 P5 恢复下载或重提），本任务暂停下载、
     不再写心跳，等它 done／过期后再按第 2 步接替。
  4) 结束（成功、失败、异常）时若心跳文件仍是本进程写的，写 done（补档任务随即升回 4 MB/s）。

认证：PO.DAAC 受保护桶需 Earthdata token。只从 workstation 既有绝对路径 `raw_archive.TOKEN_FILE` 读进内存（沿用 p5_sss_sat.read_token），
  只放进发往 Earthdata 主机（p5_sss_sat.auth_host）的 Authorization 头；重定向到 S3 预签名地址时不带；不进 argv、env、日志、manifest、
  仓库；异常文本经 redact()。不发送任何用户邮箱或个人信息（UA 用 raw_store.UA）。

Usage：
  laptop（秒级，只读 CMR）：python raw_fetch_list.py --cmr-list SPURS2_SSP=C2491772351-POCLOUD,... --list-out raw_archive_lists/spurs2.json
  workstation：python raw_fetch_list.py --run --list raw_archive_lists/spurs2.json --raw-root <bulk-disk>/.../raw --rate 2 \
          --tag raw-spurs2:$REPRO_ATTEMPT_ID --out <fast-disk>/.../runs/raw-spurs2/$REPRO_ATTEMPT_ID
  离线自测（laptop 或 workstation，约 10 s，临时目录，不联网）：python raw_fetch_list.py --selftest
Dependencies：标准库；同目录 raw_store.py、raw_archive.py（import 时无副作用）、p5_sss_sat.py（仅取认证函数）。
退出码：0 全部完成且核对一致；2 有文件多次失败（同一 JobSpec 重提即续传，已在档者跳过）；3 未预期异常／挂载哨兵不过／清单不合法；
  4 有文件字节或 MD5 与 CMR 不符（文件照样保留，manifest 记 md5_mismatch）；5 Earthdata 认证失败；6 等待 P5 超时。

Change Log：
  2026-09-27 初版（SPURS-2 八个 collection：SSP、SALINITYSNAKE、USPS、LADYAMBER、MOORING_CENTRAL、WAVEGLIDER、SAILDRONE、METEO）。
"""

import argparse
import json
import os
import re
import sys
import tempfile
import threading
import time
import traceback
import urllib.parse
import urllib.request

import raw_store as rs

VERSION = "raw-fetch-list-2026-09-27a"
CMR = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"
FOREIGN_FRESH_S = 60


# ================================================================ 清单（laptop，公开 CMR）
def cmr_list(spec, name, prefix):
    """spec：ShortName=ConceptId,…；返回清单 dict。只取 archive.podaac 的 https 下载链接（GET DATA 与 .md5 旁注）。"""
    files, colls = [], {}
    for part in spec.split(","):
        short, cid = part.strip().split("=")
        q = urllib.parse.urlencode({"collection_concept_id": cid, "page_size": 2000})
        req = urllib.request.Request(f"{CMR}?{q}", headers={"User-Agent": rs.UA})
        with urllib.request.urlopen(req, timeout=60) as r:
            hits = int(r.headers.get("CMR-Hits", "-1"))
            d = json.load(r)
        items = d.get("items", [])
        if hits != len(items):
            raise RuntimeError(f"{cid}：CMR-Hits {hits} ≠ 返回 {len(items)}（需翻页）")
        colls[short] = {"concept_id": cid, "granules": len(items), "files": 0, "bytes": 0}
        for it in items:
            u = it["umm"]
            info = {a["Name"]: a for a in u["DataGranule"]["ArchiveAndDistributionInformation"]}
            for ru in u.get("RelatedUrls", []):
                url = ru["URL"]
                if not url.startswith("https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-"):
                    continue
                fn = url.rsplit("/", 1)[1]
                a = info.get(fn)
                if a is None:
                    raise RuntimeError(f"{cid} {u['GranuleUR']}：{fn} 在 ArchiveAndDistributionInformation 里没有大小")
                ck = a.get("Checksum") or {}
                files.append({"short_name": short, "concept_id": cid, "granule": u["GranuleUR"],
                              "kind": "md5_sidecar" if fn.endswith(".md5") else "data",
                              "url": url, "rel": f"{prefix}/{short}/{fn}", "bytes": int(a["SizeInBytes"]),
                              "md5": ck.get("Value", "").lower() if ck.get("Algorithm", "").upper() == "MD5" else None,
                              "production": u["DataGranule"].get("ProductionDateTime")})
                colls[short]["files"] += 1
                colls[short]["bytes"] += int(a["SizeInBytes"])
    files.sort(key=lambda f: (f["short_name"], f["granule"], f["kind"] != "md5_sidecar"))
    return {"name": name, "built": rs.now_iso(), "source": "CMR granules.umm_json（公开元数据，无 token）",
            "rel_prefix": prefix, "collections": colls, "n_files": len(files),
            "total_bytes": sum(f["bytes"] for f in files), "files": files}


def load_list(path):
    with open(path, encoding="utf-8") as f:
        L = json.load(f)
    rels = [x["rel"] for x in L["files"]]
    if len(set(rels)) != len(rels):
        raise ValueError("清单 rel 有重复")
    for x in L["files"]:
        if ".." in x["rel"].split("/") or x["rel"].startswith("/"):
            raise ValueError(f"清单 rel 不合法：{x['rel']}")
        if not x["url"].startswith("https://"):
            raise ValueError(f"清单 url 非 https：{x['url']}")
    if L["n_files"] != len(L["files"]) or L["total_bytes"] != sum(x["bytes"] for x in L["files"]):
        raise ValueError("清单 n_files／total_bytes 与条目不符")
    return L


# ================================================================ 心跳接力
def hb_read(root):
    p = os.path.join(root, rs.HEARTBEAT_NAME)
    try:
        st = os.stat(p)
        with open(p) as f:
            s = f.read(128).strip()
    except OSError:
        return None
    m = re.search(r"pid=(\d+)", s)
    return {"state": s.split(" ", 1)[0] if s else "", "pid": int(m.group(1)) if m else None,
            "age": time.time() - st.st_mtime, "raw": s}


class Relay(threading.Thread):
    """持有 P5 的 2 MB/s 名额：写心跳、发现别的写者在下载就暂停让位（令牌桶置 0）。"""

    def __init__(self, root, bucket, rate_bps, log, every=rs.HEARTBEAT_EVERY_S, settle=20.0):
        super().__init__(daemon=True)
        self.root, self.bucket, self.rate_bps, self.log = root, bucket, rate_bps, log
        self.every, self.settle = every, settle
        self.stop = threading.Event()
        self.paused = False
        self.events = []

    def _ev(self, what):
        self.events.append([rs.now_iso(), what])
        self.log.log(f"[relay] {what}")

    def take(self):
        rs.heartbeat_write(self.root, "active")
        self.stop.wait(self.settle)
        self.bucket.set_rate(self.rate_bps)
        self.paused = False
        self._ev(f"接替：心跳已写，{self.settle:.0f} s 后以 {self.rate_bps / rs.MB:g} MB/s 下载")

    def step(self):
        h = hb_read(self.root)
        foreign = (h is not None and h["pid"] != os.getpid() and h["state"] != "done" and h["age"] < FOREIGN_FRESH_S)
        if foreign:
            if not self.paused:
                self.bucket.set_rate(0)
                self.paused = True
                self._ev(f"让位：发现别的进程在写心跳（{h['raw']}），暂停下载")
            return
        if self.paused:
            if not rs.p5_active(self.root):
                self.take()
            return
        rs.heartbeat_write(self.root, "active")

    def run(self):
        while not self.stop.wait(self.every):
            try:
                self.step()
            except Exception as e:
                self.log.log(f"[relay] 心跳出错（保持当前状态）：{e}")

    def release(self):
        self.stop.set()
        h = hb_read(self.root)
        if h is not None and h["pid"] == os.getpid():
            rs.heartbeat_write(self.root, "done")
            self._ev("结束：心跳写 done")


def wait_turn(root, max_s, log, poll=15.0):
    t0 = time.monotonic()
    last = -1
    while rs.p5_active(root):
        el = time.monotonic() - t0
        if el > max_s:
            log.log(f"[wait] 等待 P5 心跳超过 {max_s:.0f} s，退出（不下载）：{(hb_read(root) or {}).get('raw')}")
            return False, el
        if int(el // 600) != last:
            last = int(el // 600)
            log.log(f"[wait] P5 仍在下载（{(hb_read(root) or {}).get('raw')}），已等 {el / 60:.0f} min")
        time.sleep(poll)
    el = time.monotonic() - t0
    log.log(f"[wait] P5 心跳已 done 或过期（{(hb_read(root) or {}).get('raw')}），等了 {el / 60:.1f} min")
    return True, el


# ================================================================ 下载
def make_granule_fn(dl, by_gran, dataset):
    """每个颗粒：先 .md5 旁注（公开），再数据文件；都核 CMR 字节与 MD5。"""
    def fn(gkey):
        fs = by_gran[gkey]
        side = next((x for x in fs if x["kind"] == "md5_sidecar"), None)
        sidecar_md5 = None
        kinds, total = [], 0
        for x in ([side] if side else []) + [x for x in fs if x["kind"] != "md5_sidecar"]:
            extra = {"collection": x["concept_id"], "short_name": x["short_name"], "granule": x["granule"],
                     "file_kind": x["kind"], "expect_bytes": x["bytes"]}
            if x["kind"] != "md5_sidecar" and sidecar_md5:
                extra["sidecar_md5"] = sidecar_md5
                extra["sidecar_agrees_cmr"] = (sidecar_md5 == x["md5"])
            kind, val, rec = dl.get(x["url"], x["rel"], dataset, extra=extra, expect_bytes=x["bytes"],
                                    expect_md5=x["md5"])
            if x["kind"] == "md5_sidecar" and kind in ("ok", "exists"):
                with open(val, encoding="utf-8", errors="replace") as f:
                    tok = f.read().split()
                sidecar_md5 = tok[0].lower() if tok and re.fullmatch(r"[0-9a-fA-F]{32}", tok[0]) else None
            if kind == "ok":
                kinds.append("ok" if rec["status"] == "ok" else "mismatch")
                total += rec["bytes"]
            else:
                kinds.append(kind)
        for k in ("mismatch", "missing", "ok"):
            if k in kinds:
                return k, total
        return "exists", 0
    return fn


def verify_on_disk(root, L, man_rows):
    """落盘后复核：每个清单文件在档、长度等于 CMR；manifest 里本清单的 md5 记录。"""
    have = {}
    for r in man_rows:
        if r.get("path"):
            have[r["path"]] = r
    out = {"files": len(L["files"]), "present": 0, "bytes_ok": 0, "md5_ok": 0, "md5_bad": [], "absent": [],
           "sidecar_disagree": []}
    for x in L["files"]:
        p = os.path.join(root, x["rel"])
        if not os.path.exists(p):
            out["absent"].append(x["rel"])
            continue
        out["present"] += 1
        out["bytes_ok"] += os.path.getsize(p) == x["bytes"]
        r = have.get(x["rel"])
        if r and r.get("md5") and r.get("md5") == x["md5"]:
            out["md5_ok"] += 1
        elif r and r.get("md5"):
            out["md5_bad"].append(x["rel"])
        if r and r.get("sidecar_agrees_cmr") is False:
            out["sidecar_disagree"].append(x["rel"])
    out["all_ok"] = (out["present"] == out["files"] == out["bytes_ok"] == out["md5_ok"])
    return out


def run(args, out, log):
    import raw_archive as ra
    L = load_list(args.list)
    root = rs.ensure_root(args.raw_root)
    man = rs.Manifest(root)
    status = {"version": VERSION, "start": rs.now_iso(), "raw_root": root, "list": args.list, "list_name": L["name"],
              "n_files": L["n_files"], "total_bytes": L["total_bytes"], "rate_mbps": args.rate}
    spath = os.path.join(out, "raw_fetch_status.json")

    def dump():
        status["updated"] = rs.now_iso()
        tmp = spath + ".part"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(status, f, ensure_ascii=False, indent=1)
        os.replace(tmp, spath)

    dump()
    read_token, auth_host, redact, AuthError = ra.auth_bits()
    try:
        token = read_token(ra.TOKEN_FILE)
    except AuthError as e:
        log.log(f"认证：{redact(e)}")
        status["auth_error"] = redact(str(e))
        dump()
        return 5
    ok, waited = wait_turn(root, args.wait_max_s, log)
    status["waited_s"] = round(waited, 1)
    if not ok:
        status["result"] = "wait_timeout"
        dump()
        return 6
    bucket = rs.TokenBucket(args.rate * rs.MB)
    bucket.set_rate(0)
    relay = Relay(root, bucket, args.rate * rs.MB, log, settle=args.settle_s)
    status["relay_events"] = relay.events
    relay.take()
    relay.start()
    stop_tick = threading.Event()

    def ticker():
        while not stop_tick.wait(60):
            try:
                dump()
            except Exception:
                pass

    threading.Thread(target=ticker, daemon=True).start()
    dl = rs.Downloader(root, bucket, man, "rfl", args.tag, token=token, auth_host=auth_host, redact=redact, log=log)
    by_gran = {}
    for x in L["files"]:
        by_gran.setdefault(f"{x['short_name']}/{x['granule']}", []).append(x)
    rc = 0
    try:
        log.log(f"[{L['name']}] {len(by_gran)} 颗粒 {L['n_files']} 文件 {L['total_bytes'] / 1e9:.3f} GB")
        st = ra.run_pool(L["name"], sorted(by_gran), make_granule_fn(dl, by_gran, L["name"]), args.conc, log, status,
                         "fetch")
        if st.get("failed") or st.get("aborted"):
            rc = 2
    except rs.AuthFail as e:
        log.log(f"Earthdata 认证失败：{redact(e)}")
        status["auth_error"] = redact(str(e))
        rc = 5
    finally:
        relay.release()
        stop_tick.set()
        status["bytes_downloaded"] = dl.bytes
        status["http_requests"] = dl.n_req
        status["verify"] = v = verify_on_disk(root, L, [r for r in man.load() if r.get("dataset") == L["name"]])
        status["end"] = rs.now_iso()
        dump()
    log.log(f"[verify] {json.dumps({k: v[k] for k in ('files', 'present', 'bytes_ok', 'md5_ok', 'all_ok')})}"
            f" md5_bad={len(v['md5_bad'])} absent={len(v['absent'])} sidecar_disagree={len(v['sidecar_disagree'])}")
    if rc == 0 and v["md5_bad"]:
        rc = 4
    if rc == 0 and not v["all_ok"]:
        rc = 2
    return rc


# ================================================================ 离线自测（心跳接力与清单校验；不联网）
def selftest(log):
    res = []

    def chk(name, ok, detail=""):
        res.append({"name": name, "ok": bool(ok), "detail": detail})
        log.log(f"  {'ok ' if ok else 'BAD'}  {name}  {detail}")

    d = tempfile.mkdtemp(prefix="rfl-selftest-")
    hb = os.path.join(d, rs.HEARTBEAT_NAME)

    def foreign(state, age=0.0):
        with open(hb, "w") as f:
            f.write(f"{state} {rs.now_iso()} pid=1\n")
        t = time.time() - age
        os.utime(hb, (t, t))

    foreign("active")
    ok, el = wait_turn(d, 1.0, log, poll=0.3)
    chk("等待：别的进程 active 且新鲜 → 超时返回假", not ok and el > 1.0, f"{el:.1f}s")
    foreign("active", age=rs.HEARTBEAT_FRESH_S + 5)
    ok, _ = wait_turn(d, 5.0, log, poll=0.3)
    chk("等待：心跳 180 s 无更新 → 放行", ok)
    foreign("done")
    ok, _ = wait_turn(d, 5.0, log, poll=0.3)
    chk("等待：心跳 done → 放行", ok)
    b = rs.TokenBucket(2 * rs.MB)
    r = Relay(d, b, 2 * rs.MB, log, every=0.2, settle=0.3)
    r.take()
    h = hb_read(d)
    chk("接替：心跳写成 active、pid＝本进程、格式同 raw_store", h["state"] == "active" and h["pid"] == os.getpid()
        and rs.p5_active(d), h["raw"])
    chk("接替后令牌桶＝2 MB/s", b.rate == 2 * rs.MB)
    foreign("active")
    r.step()
    chk("让位：别的进程写 active → 暂停（桶置 0）且不覆盖心跳", r.paused and b.rate == 0 and hb_read(d)["pid"] == 1)
    r.step()
    chk("让位中：别的进程仍新鲜 → 继续暂停", r.paused)
    foreign("done")
    r.step()
    chk("别的进程 done → 重新接替", (not r.paused) and b.rate == 2 * rs.MB and hb_read(d)["pid"] == os.getpid())
    r.step()
    chk("正常：重写本进程心跳", hb_read(d)["pid"] == os.getpid() and hb_read(d)["state"] == "active")
    r.release()
    chk("结束：本进程心跳 → done（p5_active 为假）", hb_read(d)["state"] == "done" and not rs.p5_active(d))
    foreign("active")
    r2 = Relay(d, b, 2 * rs.MB, log)
    r2.release()
    chk("结束：心跳不是本进程写的 → 不写 done", hb_read(d)["pid"] == 1 and hb_read(d)["state"] == "active")
    # 清单校验
    here = os.path.dirname(os.path.abspath(__file__))
    lp = os.path.join(here, "raw_archive_lists", "spurs2.json")
    if os.path.exists(lp):
        L = load_list(lp)
        chk("清单 spurs2.json 合法（rel 唯一、https、合计自洽）", True, f"{L['n_files']} 文件 {L['total_bytes']} B")
        chk("清单每个文件都有 CMR MD5", all(x["md5"] and re.fullmatch(r"[0-9a-f]{32}", x["md5"]) for x in L["files"]))
    bad = os.path.join(d, "bad.json")
    with open(bad, "w") as f:
        json.dump({"name": "x", "n_files": 1, "total_bytes": 1,
                   "files": [{"rel": "../x", "url": "https://a/b", "bytes": 1}]}, f)
    try:
        load_list(bad)
        chk("清单校验：rel 含 .. 拒绝", False)
    except ValueError:
        chk("清单校验：rel 含 .. 拒绝", True)
    n_ok = sum(x["ok"] for x in res)
    log.log(f"selftest {n_ok}/{len(res)}（临时目录 {d} 留作查看，不含数据）")
    return 0 if n_ok == len(res) else 4


# ================================================================ main
class _Log:
    def __init__(self, path=None):
        self.f = open(path, "a", encoding="utf-8") if path else None
        self.lock = threading.Lock()

    def log(self, msg, echo=True):
        line = f"{rs.now_iso()} {msg}"
        with self.lock:
            if self.f:
                self.f.write(line + "\n")
                self.f.flush()
        if echo:
            print(line, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="按清单下载 PO.DAAC 原始颗粒进 bulk 原始档")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--cmr-list", metavar="SHORT=CID,…")
    g.add_argument("--run", action="store_true")
    g.add_argument("--selftest", action="store_true")
    ap.add_argument("--list-out")
    ap.add_argument("--list-name", default="spurs2")
    ap.add_argument("--list")
    ap.add_argument("--raw-root", default=rs.RAW_ROOT_DEFAULT)
    ap.add_argument("--out")
    ap.add_argument("--rate", type=float, default=2.0)
    ap.add_argument("--conc", type=int, default=2)
    ap.add_argument("--wait-max-s", type=float, default=6 * 3600)
    ap.add_argument("--settle-s", type=float, default=20.0)
    ap.add_argument("--tag", default="raw-fetch-list")
    args = ap.parse_args(argv)
    if args.cmr_list:
        L = cmr_list(args.cmr_list, args.list_name, args.list_name)
        with open(args.list_out, "w", encoding="utf-8") as f:
            json.dump(L, f, ensure_ascii=False, indent=1)
            f.write("\n")
        print(json.dumps({k: L[k] for k in ("collections", "n_files", "total_bytes")}, ensure_ascii=False, indent=1))
        return 0
    if args.selftest:
        return selftest(_Log())
    if args.rate > 2.0 or args.rate <= 0:
        print("限速越界：本任务只接替 P5 的 2 MB/s 名额", file=sys.stderr)
        return 3
    out = args.out or os.environ.get("REPRO_OUTPUT_DIR")
    if not out or not args.list:
        print("需要 --list 与 --out（或 REPRO_OUTPUT_DIR）", file=sys.stderr)
        return 3
    os.makedirs(out, exist_ok=True)
    log = _Log(os.path.join(out, "raw_fetch_log.txt"))
    log.log(f"=== start {VERSION} list={args.list} raw={args.raw_root} out={out} rate={args.rate} MB/s")
    try:
        rc = run(args, out, log)
    except Exception:
        log.log("FATAL 未预期异常：\n" + traceback.format_exc())
        rc = 3
    log.log(f"=== done rc={rc}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
