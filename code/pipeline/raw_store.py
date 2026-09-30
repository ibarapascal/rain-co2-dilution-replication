#!/usr/bin/env python3
"""raw_store.py — 原始下载文件归档（bulk 原始档）的共用 I/O 层：挂载哨兵、限速令牌桶、原件落盘路径、
raw/manifest.jsonl 追加、可续传的限速 HTTPS 下载（CMORPH 小时文件、Earthdata 受保护文件、S3 大文件）、P5 下载心跳。

只做 I/O，不含任何分析逻辑；被 `p5_sss_sat.py`（a6 起）与 `raw_archive.py` 共用。
规则：原始数据一律保留——原件一律落盘保留，本模块没有任何删除原件的代码路径
（下载中的 `.part` 只在 `<raw>/.tmp/` 里续写或覆盖重下，完成后原子改名到正式位置）。

目录约定（`RAW_ROOT_DEFAULT` 下）：
  cmorph/<YYYY>/<MM>/<DD>/CMORPH_V1.0_ADJ_8km-30min_<YYYYMMDDHH>.nc   与 NCEI 路径层级相同
  smap_l2c_v6/<YYYY>/<granule>.nc（＋ .nc.md5）                      PO.DAAC SMAP_RSS_L2_SSS_V6 完整颗粒
  era5/e5.oper.an.sfc/<YYYYMM>/<原文件名>                              NSF NCAR ERA5 镜像原件
  smap_opendap_subset/<站>/<granule>.<what>.nc4                         P5 OPeNDAP 子集的原样返回（a6 起新取的才有）
  p1/cache/、p1b/cache/、p3-cache/…、p5-cache/…                         本地已有原件的复制件（raw_archive 阶段 A）
  manifest.jsonl                                                        一行一个文件：url、path（相对 raw 根）、bytes、sha256、t、status、source_job …
  .tmp/                                                                 下载中的 .part（续传用）
  .p5_fetch_heartbeat                                                   P5 下载段心跳（raw_archive 据此切换 2↔4 MB/s）

Change Log：
  2026-09-27 初版（随 P5 a6 与 raw_archive.py 同批）。
"""

import fcntl
import hashlib
import http.client
import json
import os
import threading
import time
import urllib.parse
from datetime import datetime, timezone

import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
RAW_ROOT_DEFAULT = _rp.path("raw_root")  # [repro] 路径来自集中配置
MB = 1_000_000
CHUNK = 64 * 1024
HEARTBEAT_NAME = ".p5_fetch_heartbeat"
HEARTBEAT_FRESH_S = 180          # 心跳 mtime 在此秒数内且内容不是 done ⇒ P5 仍在下载
HEARTBEAT_EVERY_S = 20
RETRYABLE = (408, 429, 500, 502, 503, 504)
MAX_ATTEMPTS = 5
UA = "rain-co2-dilution-replication-raw/1.0 (research archive; python)"


class DownloadError(Exception):
    """重试用尽或不可重试的 HTTP 状态。"""


class NotFound404(Exception):
    pass


class AuthFail(Exception):
    pass


class _TestAbort(Exception):
    pass


def now_iso():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ---------------------------------------------------------------- 挂载哨兵
def ensure_root(root):
    """root 在 <mount-root>/<卷>/ 下时，要求该卷确为挂载点（防止卷没挂上时写进启动盘的同名空目录）。"""
    root = os.path.abspath(root)
    parts = root.split(os.sep)
    if len(parts) > 2 and parts[1] == "Volumes":
        vol = os.sep + os.path.join(parts[1], parts[2])
        if not os.path.ismount(vol):
            raise RuntimeError(f"挂载哨兵：{vol} 不是挂载点，拒绝写入 {root}")
    os.makedirs(os.path.join(root, ".tmp"), exist_ok=True)
    return root


# ---------------------------------------------------------------- 限速
class TokenBucket:
    """线程安全令牌桶；rate 单位 字节/秒；同一进程内所有下载线程共用一个。"""

    # 桶容量取 max(1 MiB, 0.5 s×rate)：workstation 上 LaunchAgent 后台进程的 sleep 会被系统合并拉长（自测实测
    # 256 KiB 容量时只跑到设定的约一半），容量够大才能在醒来时补足欠的额度；长期平均仍严格等于 rate。
    def __init__(self, rate_bps, burst=None):
        self.fixed_burst = burst
        self.rate = float(rate_bps)
        self.burst = self._burst()
        self.tokens = 0.0
        self.last = time.monotonic()
        self.lock = threading.Lock()

    def _burst(self):
        return float(self.fixed_burst) if self.fixed_burst else max(1024 * 1024, 0.5 * self.rate)

    def set_rate(self, rate_bps):
        with self.lock:
            self.rate = float(rate_bps)
            self.burst = self._burst()

    def consume(self, n):
        while True:
            with self.lock:
                t = time.monotonic()
                self.tokens = min(self.burst, self.tokens + (t - self.last) * self.rate)
                self.last = t
                if self.tokens >= n:
                    self.tokens -= n
                    return
                wait = (n - self.tokens) / max(self.rate, 1.0)
            time.sleep(min(wait, 0.5))


# ---------------------------------------------------------------- manifest
class Manifest:
    def __init__(self, root):
        self.path = os.path.join(root, "manifest.jsonl")
        self.lock = threading.Lock()

    def add(self, **rec):
        rec.setdefault("t", now_iso())
        line = (json.dumps(rec, ensure_ascii=False, sort_keys=False) + "\n").encode("utf-8")
        with self.lock:
            fd = os.open(self.path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
            try:
                fcntl.flock(fd, fcntl.LOCK_EX)       # P5 与 raw_archive 两个进程同时追加
                os.write(fd, line)
                fcntl.flock(fd, fcntl.LOCK_UN)
            finally:
                os.close(fd)
        return rec

    def load(self):
        out = []
        if os.path.exists(self.path):
            with open(self.path, encoding="utf-8") as f:
                for ln in f:
                    try:
                        out.append(json.loads(ln))
                    except ValueError:
                        continue
        return out


# ---------------------------------------------------------------- 路径
def cmorph_path(h):
    d = datetime.fromtimestamp(h * 3600, tz=timezone.utc)
    return (f"/data/cmorph-high-resolution-global-precipitation-estimates/access/30min/8km/"
            f"{d:%Y/%m/%d}/CMORPH_V1.0_ADJ_8km-30min_{d:%Y%m%d%H}.nc")


CMORPH_HOST = "www.ncei.noaa.gov"


def cmorph_url(h):
    return "https://" + CMORPH_HOST + cmorph_path(h)


def cmorph_rel(h):
    d = datetime.fromtimestamp(h * 3600, tz=timezone.utc)
    return f"cmorph/{d:%Y/%m/%d}/CMORPH_V1.0_ADJ_8km-30min_{d:%Y%m%d%H}.nc"


def cmorph_rel_from_name(name):
    """CMORPH_V1.0_ADJ_8km-30min_YYYYMMDDHH.nc → cmorph/YYYY/MM/DD/<name>。"""
    s = name.rsplit("_", 1)[1][:10]
    return f"cmorph/{s[0:4]}/{s[4:6]}/{s[6:8]}/{name}"


# ---------------------------------------------------------------- 心跳
def heartbeat_write(root, state="active"):
    p = os.path.join(root, HEARTBEAT_NAME)
    tmp = p + f".{os.getpid()}.tmp"
    with open(tmp, "w") as f:
        f.write(f"{state} {now_iso()} pid={os.getpid()}\n")
    os.replace(tmp, p)


def p5_active(root):
    p = os.path.join(root, HEARTBEAT_NAME)
    try:
        st = os.stat(p)
        with open(p) as f:
            s = f.read(64)
    except OSError:
        return False
    return (not s.startswith("done")) and (time.time() - st.st_mtime) < HEARTBEAT_FRESH_S


# ---------------------------------------------------------------- 下载
class Downloader:
    """限速、可续传（Range）、逐跳处理重定向的 HTTPS GET。

    - 连接按 (线程, host) 复用；重定向逐跳判断，只有 auth_host(host) 为真时才附 Authorization（token 由调用方给，不落任何文件）。
    - 落盘：<root>/.tmp/<rel 扁平化>.<owner>.part → 核长度 → sha256（流式）→ os.replace 到 <root>/<rel>。
    - 已存在的 <root>/<rel> 视为已归档，直接返回，不重下。
    """

    def __init__(self, root, bucket, manifest, owner, source_job, token=None, auth_host=None, redact=None, log=None):
        self.root = root
        self.bucket = bucket
        self.man = manifest
        self.owner = owner
        self.source_job = source_job
        self.token = token
        self.auth_host = auth_host or (lambda h: False)
        self.redact = redact or (lambda s: str(s))
        self.log = log
        self.local = threading.local()
        self.lock = threading.Lock()
        self.bytes = 0
        self.n_req = 0
        self.hb_root = None
        self.hb_last = 0.0
        self.abort_after = None      # 仅自测：读到该字节数后模拟断线，用来验证 Range 续传

    # 连接
    def _conn(self, host, fresh=False):
        cs = getattr(self.local, "cs", None)
        if cs is None:
            cs = self.local.cs = {}
        c = cs.get(host)
        if fresh and c is not None:
            try:
                c.close()
            except Exception:
                pass
            c = None
        if c is None:
            c = http.client.HTTPSConnection(host, timeout=180)
            cs[host] = c
        return c

    def _drop(self, host=None):
        """丢弃本线程的连接（host=None 时全部），读写出错后用。"""
        cs = getattr(self.local, "cs", None) or {}
        for k in ([host] if host else list(cs)):
            c = cs.pop(k, None)
            if c is not None:
                try:
                    c.close()
                except Exception:
                    pass

    def _beat(self):
        if self.hb_root and time.monotonic() - self.hb_last > HEARTBEAT_EVERY_S:
            self.hb_last = time.monotonic()
            try:
                heartbeat_write(self.hb_root, "active")
            except OSError:
                pass

    def final_path(self, rel):
        return os.path.join(self.root, rel)

    def part_path(self, rel):
        return os.path.join(self.root, ".tmp", rel.replace("/", "__") + f".{self.owner}.part")

    def _open(self, url, start):
        """返回 (response, host)，已跟完重定向；response.status 为 200 或 206。"""
        cur = url
        for _hop in range(8):
            u = urllib.parse.urlsplit(cur)
            if u.scheme != "https":
                raise DownloadError(f"拒绝非 https：{u.scheme}://{u.netloc}")
            host = (u.hostname or "").lower()
            if host == "urs.earthdata.nasa.gov":
                raise AuthFail("被重定向到 Earthdata 登录页（token 未被接受）")
            path = u.path + ("?" + u.query if u.query else "")
            h = {"User-Agent": UA, "Connection": "keep-alive"}
            if start > 0:
                h["Range"] = f"bytes={start}-"
            if self.token and self.auth_host(host):
                h["Authorization"] = "Bearer " + self.token
            c = self._conn(host)
            try:
                c.request("GET", path, headers=h)
                r = c.getresponse()
            except Exception:
                self._drop()
                raise
            with self.lock:
                self.n_req += 1
            if r.status in (301, 302, 303, 307, 308):
                loc = r.headers.get("Location")
                r.read()
                if not loc:
                    raise DownloadError(f"HTTP {r.status} 无 Location @ {host}")
                cur = urllib.parse.urljoin(cur, loc)
                continue
            if r.status in (200, 206):
                return r, host
            r.read()
            if r.status == 404:
                raise NotFound404(f"HTTP 404 @ {host}")
            if r.status in (401, 403) and self.auth_host(host):
                raise AuthFail(f"HTTP {r.status} @ {host}")
            if r.status == 416:
                return r, host
            raise DownloadError(f"HTTP {r.status} @ {host}", r.status)
        raise DownloadError("重定向超过 8 次")

    def get(self, url, rel, dataset, extra=None, expect_bytes=None, expect_sha256=None, expect_md5=None,
            n404_missing=2, attempts=MAX_ATTEMPTS):
        """下载 url 到 <root>/<rel>。返回 ("exists"|"ok", 绝对路径, rec) 或 ("missing", 原因, rec)。
        给了 expect_sha256／expect_md5 时核对；不符记 status=sha_mismatch／md5_mismatch，文件照样落到正式位置保留。"""
        dest = self.final_path(rel)
        if os.path.exists(dest):
            return "exists", dest, None
        part = self.part_path(rel)
        os.makedirs(os.path.dirname(part), exist_ok=True)
        last, n404 = None, 0
        for attempt in range(1, attempts + 1):
            start = os.path.getsize(part) if os.path.exists(part) else 0
            if expect_bytes is not None and start > expect_bytes:
                start = 0
            try:
                r, host = self._open(url, start)
                if r.status == 416 or (r.status == 200 and start > 0):
                    # 服务器不认 Range 或 part 已满：从头重下（覆盖 part，不删）
                    if r.status == 416:
                        r.read()
                        start = 0
                        r, host = self._open(url, 0)
                    else:
                        start = 0
                total = None
                cr = r.headers.get("Content-Range")
                if r.status == 206 and cr and "/" in cr:
                    try:
                        total = int(cr.rsplit("/", 1)[1])
                    except ValueError:
                        total = None
                elif r.status == 200:
                    cl = r.headers.get("Content-Length")
                    total = int(cl) if cl else None
                got = 0
                with open(part, "r+b" if start > 0 else "wb") as f:
                    if start > 0:
                        f.seek(start)
                        f.truncate()
                    while True:
                        self.bucket.consume(CHUNK)
                        b = r.read(CHUNK)
                        if not b:
                            break
                        f.write(b)
                        got += len(b)
                        with self.lock:
                            self.bytes += len(b)
                        self._beat()
                        if self.abort_after is not None and got >= self.abort_after:
                            self.abort_after = None
                            raise _TestAbort("自测：模拟中途断线")
                size = start + got
                if total is not None and size != total:
                    raise IOError(f"长度 {size} ≠ {total}")
                if expect_bytes is not None and size != expect_bytes:
                    raise IOError(f"长度 {size} ≠ 预期 {expect_bytes}")
                sha = sha256_file(part)
                status = "ok"
                rec = dict(dataset=dataset, url=url, path=rel, bytes=size, sha256=sha, status=status,
                           source_job=self.source_job, resumed=start > 0)
                if expect_sha256:
                    rec["expect_sha256"] = expect_sha256
                    if sha != expect_sha256:
                        rec["status"] = "sha_mismatch"
                if expect_md5:
                    md5 = md5_file(part)
                    rec["md5"], rec["expect_md5"] = md5, expect_md5
                    if md5 != expect_md5 and rec["status"] == "ok":
                        rec["status"] = "md5_mismatch"
                et = (r.headers.get("ETag") or "").strip('"')
                if et:
                    rec["etag"] = et
                if extra:
                    rec.update(extra)
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                os.replace(part, dest)
                self.man.add(**rec)
                return "ok", dest, rec
            except NotFound404:
                n404 += 1
                last = "HTTP 404"
                if n404 >= n404_missing:
                    rec = self.man.add(dataset=dataset, url=url, path=rel, bytes=0, sha256=None,
                                       status="missing_404", source_job=self.source_job)
                    return "missing", "HTTP404", rec
            except (AuthFail, _TestAbort):
                raise
            except DownloadError as e:
                last = self.redact(e)
                code = e.args[1] if len(e.args) > 1 else None
                if code is not None and code not in RETRYABLE:
                    raise DownloadError(f"{url} → {last}")
            except Exception as e:
                last = self.redact(f"{type(e).__name__}: {e}")
                self._drop()
            if attempt < attempts:
                time.sleep(min(40, 5 * 2 ** (attempt - 1)))
        raise DownloadError(f"{url} → {attempts} 次失败：{last}")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class CmorphRaw:
    """P5（及 raw_archive）用的 CMORPH 小时取数：raw 档里有就直接用，没有就限速下载进 raw 档；接口同 p2.Fetcher.fetch。

    fetch(h) → ("ok", 原件绝对路径) 或 ("missing", 原因)；重试用尽抛 exc_cls（P5 传 p2.DataSourceError，行为与 p2.Fetcher 一致）。
    原件永远留在 raw 档，调用方不得删除。
    """

    def __init__(self, root, rate_mbps, owner, source_job, exc_cls=DownloadError, log=None, heartbeat=False,
                 bucket=None, manifest=None):
        self.root = ensure_root(root)
        self.bucket = bucket or TokenBucket(rate_mbps * MB)
        self.man = manifest or Manifest(self.root)
        self.dl = Downloader(self.root, self.bucket, self.man, owner, source_job, log=log)
        if heartbeat:
            self.dl.hb_root = self.root
        self.exc_cls = exc_cls
        self.reused = 0
        self.downloaded = 0

    @property
    def bytes(self):
        return self.dl.bytes

    @property
    def n_req(self):
        return self.dl.n_req

    def fetch(self, h):
        try:
            kind, val, _rec = self.dl.get(cmorph_url(h), cmorph_rel(h), "cmorph")
        except (DownloadError, AuthFail) as e:
            raise self.exc_cls(str(e))
        if kind == "exists":
            with self.dl.lock:
                self.reused += 1
            return "ok", val
        if kind == "ok":
            with self.dl.lock:
                self.downloaded += 1
            return "ok", val
        return "missing", val
