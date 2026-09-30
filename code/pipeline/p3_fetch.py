"""P3 下载器（只下载、不做任何 P3 计算）：把 P3 全部输入落到 fast_root 下的 p3-cache/。

用途：P3（全球重算）的数据准备（只用标准库＋numpy＋netCDF4）。
目录（根＝--cache）：
  cmorph/YYYY/MM/<原文件名>            CMORPH CDR V1.0 ADJ 8 km/30 min 逐小时文件（原样）
  era5_ws10/YYYYMM.nc                  ERA5 u10/v10（NCAR S3 镜像）→ 60°S–60°N → ws10 int16×0.01 m/s，(time, latitude, longitude)
  hycom_sss/YYYYMMDD.nc                HYCOM GLBv0.08 expt_53.X 表层盐度，NCSS 每天 12Z（原样）
  oisst/<原文件名>                     OISST v2.1 逐日（原样）
  watson/<tar 名>＋tar 内全部成员       Watson 等 2020 RECCAP2 表层 fCO2（Zenodo）
  glodap/<tar.gz 名>＋tar 内全部 .nc    GLODAPv2.2016b mapped
  month_ready/YYYYMM.json              该月（200001 含 1999-12-31 起转与静态源）全部源齐了才写
  manifest.jsonl                       每个落盘文件一行：source、url、path（相对根）、bytes、sha256（只追加，后写覆盖前写）
  fetch_done.json                      全部完成时写：各源文件数、总字节、耗时
顺序：先 200001 组（1999-12-31 起转＋1 月＋Watson＋GLODAP＋ERA5 199912/200001），齐了写 month_ready/200001.json，再逐月推进；
  每月下载完即开下月，ERA5 转换在后台单线程做，转完才写该月 month_ready。
并发：NCEI（CMORPH＋OISST）4 路持久连接；S3 单文件分 32 MiB 段 8 路；HYCOM NCSS 2 路；Zenodo/GLODAP 各 1 路。
断点续传：manifest 里 ok 且文件大小一致就跳过（--verify 时再核 sha256）；S3 大文件按段记进度（.part.json）；
  单流下载的 .part 用 Range 续传；ERA5 原文件在 manifest 里记 sha，转换校验通过后才删。
退出码：0 全部完成；2 有文件失败（可原样重提续跑）；3 前置条件不满足（磁盘、参数）。
用法：p3_fetch.py --cache <fast_root>/p3-cache [--verify] [--selftest]
  --selftest：不联网，用合成 u/v 文件测 ERA5 转换与校验、manifest 跳过逻辑、任务清单计数，并报峰值 RSS。
版本：p3-fetch-2026-09-26a
"""
import argparse, concurrent.futures as cf, datetime as dt, hashlib, http.client, json, math, os, resource
import shutil, sys, tarfile, threading, time, traceback
from urllib.parse import urlsplit, urljoin

VERSION = "p3-fetch-2026-09-26b"  # b：derive_ws10 先建输出目录（a 版漏建致 ERA5 转换全部失败）
UA = "rain-co2-dilution-replication-p3/1.0 (research; python)"
YEAR = 2000
SPINUP_DAY = dt.date(1999, 12, 31)
CMORPH_BASE = "https://www.ncei.noaa.gov/data/cmorph-high-resolution-global-precipitation-estimates/access/30min/8km"
OISST_BASE = "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr"
ERA5_BASE = "https://nsf-ncar-era5.s3.amazonaws.com/e5.oper.an.sfc"
ERA5_VARS = {"u": "128_165_10u", "v": "128_166_10v"}
HYCOM_NCSS = "https://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_53.X/data"
WATSON_URL = "https://zenodo.org/records/7990823/files/surface_co2__UOEX_Wat20_1985_2019_v20211204.tar?download=1"
WATSON_NAME = "surface_co2__UOEX_Wat20_1985_2019_v20211204.tar"
GLODAP_URL = "https://glodap.info/glodap_files/v2.2023/GLODAPv2.2016b.MappedProduct.tar.gz"
GLODAP_NAME = "GLODAPv2.2016b.MappedProduct.tar.gz"
N_NCEI, N_S3, N_HYCOM, N_MISC = 4, 8, 2, 2
S3_CHUNK = 32 << 20
MAX_ATTEMPTS = 6
MAX_TOTAL_FAIL = 200
PRINT_EVERY = 250
HEARTBEAT_S = 900
EST = {"cmorph": 1.64e6, "oisst": 1.69e6, "hycom": 30e6, "era5_raw": 1.52e9, "watson": 675.8e6, "glodap": 211.5e6}
WS_SCALE, WS_FILL = 0.01, -32767
NC_LOCK = threading.Lock()  # HDF5 非线程安全：所有 netCDF4 调用串行
RETRYABLE = (408, 429, 500, 502, 503, 504)


class Fatal(Exception):
    """不可重试（如 403/非预期 404）。"""


# ======================================================================== 任务清单

def days_of_month(y, m):
    d0 = dt.date(y, m, 1)
    d1 = dt.date(y + (m == 12), m % 12 + 1, 1)
    return [d0 + dt.timedelta(i) for i in range((d1 - d0).days)]


def cmorph_task(day, hh):
    name = f"CMORPH_V1.0_ADJ_8km-30min_{day:%Y%m%d}{hh:02d}.nc"
    return {"source": "cmorph", "url": f"{CMORPH_BASE}/{day:%Y/%m/%d}/{name}", "rel": f"cmorph/{day:%Y/%m}/{name}"}


def oisst_task(day):
    name = f"oisst-avhrr-v02r01.{day:%Y%m%d}.nc"
    return {"source": "oisst", "url": f"{OISST_BASE}/{day:%Y%m}/{name}", "rel": f"oisst/{name}"}


def hycom_task(day):
    q = (f"?var=salinity&north=60&south=-60&west=-180&east=180&horizStride=1"
         f"&time={day:%Y-%m-%d}T12:00:00Z&vertCoord=0&accept=netcdf4")
    return {"source": "hycom", "url": f"{HYCOM_NCSS}/{day.year}{q}", "rel": f"hycom_sss/{day:%Y%m%d}.nc"}


def era5_names(ym):
    y, m = divmod(ym, 100)
    nd = len(days_of_month(y, m))
    return {k: f"e5.oper.an.sfc.{code}.ll025sc.{ym}0100_{ym}{nd:02d}23.nc" for k, code in ERA5_VARS.items()}


def era5_task(ym):
    n = era5_names(ym)
    return {"source": "era5_ws10", "ym": ym, "rel": f"era5_ws10/{ym}.nc",
            "raw": {k: {"source": "era5_raw", "url": f"{ERA5_BASE}/{ym}/{v}", "rel": f"era5_raw/{v}"} for k, v in n.items()}}


def build_groups():
    """[(月标签 YYYYMM, days, era5 月, 静态源)]；200001 组含起转日、ERA5 199912 与 Watson/GLODAP。"""
    groups = []
    for m in range(1, 13):
        days = days_of_month(YEAR, m)
        era = [YEAR * 100 + m]
        static = []
        if m == 1:
            days = [SPINUP_DAY] + days
            era = [199912] + era
            static = [{"source": "watson", "url": WATSON_URL, "rel": f"watson/{WATSON_NAME}", "extract": "all"},
                      {"source": "glodap", "url": GLODAP_URL, "rel": f"glodap/{GLODAP_NAME}", "extract": "nc"}]
        groups.append((YEAR * 100 + m, days, era, static))
    return groups


# ======================================================================== manifest 与进度

class Manifest:
    def __init__(self, root):
        self.root = root
        self.path = os.path.join(root, "manifest.jsonl")
        self.recs = {}
        self.lock = threading.Lock()
        if os.path.exists(self.path):
            with open(self.path, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        r = json.loads(line)
                    except ValueError:
                        continue  # 崩溃时最后一行可能半截
                    self.recs[r["path"]] = r

    def add(self, **rec):
        rec = {"t": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), **rec}
        line = json.dumps(rec, ensure_ascii=False) + "\n"
        with self.lock:
            with open(self.path, "a", encoding="utf-8") as f:
                f.write(line)
                f.flush()
                os.fsync(f.fileno())
            self.recs[rec["path"]] = rec
        return rec

    def have(self, rel, verify=False):
        r = self.recs.get(rel)
        if not r:
            return False
        if r.get("status") == "missing":
            return True
        if r.get("status") == "done":  # tar 解包标记：成员都在才算
            return all(self.have(x, verify) for x in r.get("members", []))
        if r.get("status") != "ok":
            return False
        p = os.path.join(self.root, rel)
        if not os.path.isfile(p) or os.path.getsize(p) != r["bytes"]:
            return False
        return not verify or sha256_file(p) == r["sha256"]


class Stats:
    def __init__(self, remaining_est):
        self.lock = threading.RLock()
        self.t0 = time.monotonic()
        self.bytes = 0
        self.files = 0
        self.remaining = remaining_est
        self.win = [(self.t0, 0)]
        self.group = "-"
        self.total_files = 0
        self.last_print = self.t0

    def add_bytes(self, n):
        with self.lock:
            self.bytes += n
            self.remaining = max(0.0, self.remaining - n)

    def line(self, tag):
        with self.lock:
            return self._line(tag)

    def _line(self, tag):
        now = time.monotonic()
        el = now - self.t0
        rate = self.bytes / 1e6 / max(el, 1e-6)
        tw, bw = self.win[0]
        recent = (self.bytes - bw) / 1e6 / max(now - tw, 1e-6)
        self.win.append((now, self.bytes))
        self.win = [x for x in self.win if now - x[0] <= 1800] or [(now, self.bytes)]
        eta = self.remaining / 1e6 / max(recent if recent > 0 else rate, 1e-6) / 3600
        self.last_print = now
        return (f"[P3 {tag}] {self.files}/{self.total_files} 文件，{self.bytes / 1e9:.2f} GB，平均 {rate:.2f} MB/s，"
                f"近 30 min {recent:.2f} MB/s，当前组 {self.group}，预计剩余 {eta:.2f} h（约 {self.remaining / 1e9:.1f} GB）")

    def file_done(self):
        with self.lock:
            self.files += 1
            if self.files % PRINT_EVERY == 0:
                print(self._line("进度"), flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(8 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def log(msg):
    ts = dt.datetime.now(dt.timezone.utc).strftime("%H:%M:%SZ")
    line = f"{ts} {msg}"
    if LOGF:
        with LOG_LOCK:
            LOGF.write(line + "\n")
            LOGF.flush()
    return line


LOGF = None
LOG_LOCK = threading.Lock()


# ======================================================================== HTTP

class Conns:
    """每线程每 host 一条持久 HTTPS 连接（P2 Fetcher 同款）。"""

    def __init__(self):
        self.local = threading.local()

    def get(self, host, fresh=False):
        d = getattr(self.local, "d", None)
        if d is None:
            d = self.local.d = {}
        c = d.get(host)
        if fresh and c is not None:
            try:
                c.close()
            except Exception:
                pass
            c = None
        if c is None:
            c = http.client.HTTPSConnection(host, timeout=900 if "hycom" in host else 180)
            d[host] = c
        return c


CONNS = Conns()


def backoff(attempt):
    time.sleep(min(60, 5 * 2 ** (attempt - 1)))


def get_whole(task, root, man, stats, allow_404=False, validate=None):
    """单流下载到 <rel>.part（支持 Range 续传）→ 校验 → 改名 → 记 manifest。返回记录。"""
    url0, rel, source = task["url"], task["rel"], task["source"]
    dest = os.path.join(root, rel)
    part = dest + ".part"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    last, n404 = None, 0
    for attempt in range(1, MAX_ATTEMPTS + 1):
        url = url0
        try:
            for _hop in range(6):
                u = urlsplit(url)
                c = CONNS.get(u.netloc, fresh=attempt > 1)
                have = os.path.getsize(part) if os.path.exists(part) else 0
                hdr = {"User-Agent": UA, "Connection": "keep-alive"}
                if have:
                    hdr["Range"] = f"bytes={have}-"
                c.request("GET", u.path + ("?" + u.query if u.query else ""), headers=hdr)
                r = c.getresponse()
                if r.status in (301, 302, 303, 307, 308):
                    r.read()
                    url = urljoin(url, r.headers["Location"])
                    continue
                break
            if r.status in (200, 206):
                h = hashlib.sha256()
                if r.status == 206 and have:
                    total = int(r.headers.get("Content-Range", "*/-1").split("/")[-1])
                    with open(part, "rb") as f:
                        while True:
                            b = f.read(8 << 20)
                            if not b:
                                break
                            h.update(b)
                    mode, got = "ab", have
                else:
                    clen = int(r.headers.get("Content-Length") or -1)
                    total = clen
                    mode, got = "wb", 0
                with open(part, mode) as f:
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        f.write(b)
                        h.update(b)
                        got += len(b)
                        stats.add_bytes(len(b))
                if total >= 0 and got != total:
                    raise IOError(f"长度 {got} ≠ {total}")
                try:
                    extra = validate(part) if validate else {}
                except Exception as e:
                    os.remove(part)
                    raise IOError(f"校验不过：{type(e).__name__}: {e}")
                os.replace(part, dest)
                rec = man.add(source=source, url=url0, path=rel, bytes=got, sha256=h.hexdigest(), status="ok", **extra)
                stats.file_done()
                return rec
            r.read()
            if r.status == 404 and allow_404:
                n404 += 1
                if n404 >= 2:
                    rec = man.add(source=source, url=url0, path=rel, status="missing", reason="HTTP404")
                    log(f"{rel} 缺测（HTTP404×2）")
                    stats.file_done()
                    return rec
                last = "HTTP 404"
            elif r.status == 416:
                os.remove(part)
                last = "HTTP 416（丢弃 .part 重下）"
            else:
                last = f"HTTP {r.status}"
                if r.status not in RETRYABLE and r.status != 404:
                    raise Fatal(f"{url} → {last}")
        except Fatal:
            raise
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
            CONNS.get(urlsplit(url).netloc, fresh=True)
        if attempt < MAX_ATTEMPTS:
            backoff(attempt)
    raise IOError(f"{rel}：{MAX_ATTEMPTS} 次失败：{last}")


def s3_head(url):
    u = urlsplit(url)
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            c = CONNS.get(u.netloc, fresh=attempt > 1)
            c.request("HEAD", u.path, headers={"User-Agent": UA, "Connection": "keep-alive"})
            r = c.getresponse()
            r.read()
            if r.status == 200:
                return int(r.headers["Content-Length"]), r.headers.get("ETag", "").strip('"')
            if r.status not in RETRYABLE:
                raise Fatal(f"HEAD {url} → HTTP {r.status}")
        except Fatal:
            raise
        except Exception:
            CONNS.get(u.netloc, fresh=True)
        backoff(attempt)
    raise IOError(f"HEAD {url} 失败")


def s3_chunk(url, fd, i, size, stats):
    a = i * S3_CHUNK
    b = min(size, a + S3_CHUNK) - 1
    u = urlsplit(url)
    last = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            c = CONNS.get(u.netloc, fresh=attempt > 1)
            c.request("GET", u.path, headers={"User-Agent": UA, "Connection": "keep-alive", "Range": f"bytes={a}-{b}"})
            r = c.getresponse()
            if r.status != 206:
                r.read()
                last = f"HTTP {r.status}"
                if r.status not in RETRYABLE:
                    raise Fatal(f"{url} 段 {i} → {last}")
            else:
                off, got = a, 0
                while True:
                    buf = r.read(1 << 20)
                    if not buf:
                        break
                    os.pwrite(fd, buf, off)
                    off += len(buf)
                    got += len(buf)
                    stats.add_bytes(len(buf))
                if got == b - a + 1:
                    return i
                last = f"段长 {got} ≠ {b - a + 1}"
        except Fatal:
            raise
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
            CONNS.get(u.netloc, fresh=True)
        if attempt < MAX_ATTEMPTS:
            backoff(attempt)
    raise IOError(f"{url} 段 {i}：{last}")


def s3_download(task, root, man, stats, pool):
    """大文件 8 路分段并发；.part.json 记已完成段，可跨任务续传。"""
    url, rel = task["url"], task["rel"]
    dest = os.path.join(root, rel)
    part, state = dest + ".part", dest + ".part.json"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    size, etag = s3_head(url)
    n = math.ceil(size / S3_CHUNK)
    done = set()
    if os.path.exists(part) and os.path.exists(state):
        try:
            st = json.load(open(state))
            if st.get("etag") == etag and st.get("size") == size and os.path.getsize(part) == size:
                done = set(st["done"])
        except ValueError:
            pass
    if not done:
        with open(part, "wb") as f:
            f.truncate(size)
    lock = threading.Lock()
    fd = os.open(part, os.O_RDWR)
    try:
        futs = [pool.submit(s3_chunk, url, fd, i, size, stats) for i in range(n) if i not in done]
        err = None
        for fut in cf.as_completed(futs):
            try:
                i = fut.result()
            except Exception as e:
                err = err or e
                continue
            with lock:
                done.add(i)
                tmp = state + ".tmp"
                with open(tmp, "w") as f:
                    json.dump({"url": url, "etag": etag, "size": size, "chunk": S3_CHUNK, "done": sorted(done)}, f)
                os.replace(tmp, state)
        if err:
            raise err
        os.fsync(fd)
    finally:
        os.close(fd)
    sha = sha256_file(part)
    os.replace(part, dest)
    os.remove(state)
    rec = man.add(source=task["source"], url=url, path=rel, bytes=size, sha256=sha, etag=etag, status="ok")
    stats.file_done()
    return rec


# ======================================================================== 校验与派生

def validate_hycom(path):
    with open(path, "rb") as f:
        if f.read(8) != b"\x89HDF\r\n\x1a\n":
            raise ValueError("HYCOM 返回的不是 netCDF4")
    import netCDF4
    with NC_LOCK:
        with netCDF4.Dataset(path) as ds:
            v = ds.variables["salinity"]
            lat = ds.variables["lat" if "lat" in ds.variables else "latitude"][:]
            t = ds.variables["time"]
            tval = netCDF4.num2date(t[:], t.units, getattr(t, "calendar", "standard"))
            shape = list(v.shape)
    if lat.min() > -59.9 or lat.max() < 59.9:
        raise ValueError(f"HYCOM 纬度范围不足 {lat.min()}..{lat.max()}")
    return {"hycom_time": str(tval[0]), "shape": shape}


def find_var(ds, key):
    for name, v in ds.variables.items():
        if v.ndim == 3 and key in name.upper():
            return v
    cands = [v for v in ds.variables.values() if v.ndim == 3]
    if len(cands) == 1:
        return cands[0]
    raise KeyError(f"找不到 {key} 变量：{list(ds.variables)}")


def derive_ws10(u_path, v_path, out_path, expect_nt=None, meta=None):
    """u10/v10 → 60°S–60°N 风速 int16（scale 0.01，offset 0，fill −32767），逐时间块处理；返回校验摘要。"""
    import numpy as np, netCDF4
    tmp = out_path + ".tmp"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with NC_LOCK:
        du, dv = netCDF4.Dataset(u_path), netCDF4.Dataset(v_path)
    try:
        with NC_LOCK:
            uvar, vvar = find_var(du, "10U"), find_var(dv, "10V")
            assert uvar.dimensions == vvar.dimensions and uvar.shape == vvar.shape
            tdim, ydim, xdim = uvar.dimensions
            lat = du.variables[ydim][:]
            lon = du.variables[xdim][:]
            tim = du.variables[tdim]
            assert np.array_equal(lat, dv.variables[ydim][:]) and np.array_equal(tim[:], dv.variables[tdim][:])
            assert lat[0] > lat[-1], "纬度应为降序"
            sel = np.where((lat >= -60 - 1e-6) & (lat <= 60 + 1e-6))[0]
            i0, i1 = int(sel[0]), int(sel[-1]) + 1
            nt, ny, nx = uvar.shape[0], i1 - i0, len(lon)
            if expect_nt is not None:
                assert nt == expect_nt, f"时间步 {nt} ≠ {expect_nt}"
            ch = uvar.chunking()
            tb = ch[0] if isinstance(ch, list) and 0 < ch[0] <= 96 else 24
            for v in (uvar, vvar):
                v.set_auto_mask(False)
            fills = [getattr(v, "_FillValue", None) for v in (uvar, vvar)]
            out = netCDF4.Dataset(tmp, "w", format="NETCDF4")
            out.createDimension("time", nt)
            out.createDimension("latitude", ny)
            out.createDimension("longitude", nx)
            for name, src, sl in (("time", tim, slice(None)), ("latitude", du.variables[ydim], slice(i0, i1)),
                                  ("longitude", du.variables[xdim], slice(None))):
                dim = {"time": "time", "latitude": "latitude", "longitude": "longitude"}[name]
                o = out.createVariable(name, src.dtype, (dim,))
                o.setncatts({k: src.getncattr(k) for k in src.ncattrs() if k != "_FillValue"})
                o[:] = src[sl]
            ws = out.createVariable("ws10", "i2", ("time", "latitude", "longitude"), zlib=True, complevel=4,
                                    shuffle=True, chunksizes=(1, ny, nx), fill_value=WS_FILL)
            ws.set_auto_maskandscale(False)
            ws.setncatts({"scale_factor": WS_SCALE, "add_offset": 0.0, "units": "m s-1",
                          "long_name": "10 m wind speed sqrt(u10^2+v10^2) from ERA5 (NCAR ds633.0 mirror)"})
            out.setncatts({"title": "ERA5 10 m wind speed, 60S-60N, for the rain CO2 dilution stage P3",
                           "created_by": VERSION, **({k: json.dumps(v) for k, v in (meta or {}).items()})})
        n_fill = 0
        for t0 in range(0, nt, tb):
            t1 = min(nt, t0 + tb)
            with NC_LOCK:
                u = np.asarray(uvar[t0:t1, i0:i1, :], dtype=np.float32)
                v = np.asarray(vvar[t0:t1, i0:i1, :], dtype=np.float32)
            bad = ~np.isfinite(u) | ~np.isfinite(v)
            for a, fv in ((u, fills[0]), (v, fills[1])):
                if fv is not None:
                    bad |= a == np.float32(fv)
            s = np.hypot(u, v)
            del u, v
            p = np.rint(s * np.float32(1.0 / WS_SCALE))
            del s
            np.clip(p, 0, 32766, out=p)
            p = p.astype(np.int16)
            p[bad] = WS_FILL
            n_fill += int(bad.sum())
            with NC_LOCK:
                ws[t0:t1] = p
            del p, bad
        with NC_LOCK:
            out.close()
        chk = verify_ws10(tmp, uvar, vvar, i0, i1, nt)
        chk.update({"nt": nt, "ny": ny, "nx": nx, "lat_first": float(lat[i0]), "lat_last": float(lat[i1 - 1]),
                    "time_block": tb, "src_chunking": ch if isinstance(ch, list) else str(ch), "n_fill": n_fill})
    finally:
        with NC_LOCK:
            du.close()
            dv.close()
    os.replace(tmp, out_path)
    return chk


def verify_ws10(path, uvar, vvar, i0, i1, nt):
    import numpy as np, netCDF4
    worst = 0.0
    with NC_LOCK:
        ds = netCDF4.Dataset(path)
    try:
        for t in sorted({0, nt // 2, nt - 1}):
            with NC_LOCK:
                w = ds.variables["ws10"]
                assert w.dtype == np.int16 and abs(w.scale_factor - WS_SCALE) < 1e-12 and w.add_offset == 0
                dec = np.ma.filled(np.ma.asarray(w[t]).astype(np.float64), np.nan)  # auto scale → m/s；填充值→NaN
                u = np.asarray(uvar[t, i0:i1, :], dtype=np.float64)
                v = np.asarray(vvar[t, i0:i1, :], dtype=np.float64)
            ref = np.hypot(u, v)
            d = np.abs(dec - ref)
            worst = max(worst, float(np.nanmax(d)))
        with NC_LOCK:
            lat = ds.variables["latitude"][:]
            dims = ds.variables["ws10"].dimensions
    finally:
        with NC_LOCK:
            ds.close()
    assert dims == ("time", "latitude", "longitude"), dims
    assert lat[0] > lat[-1] and abs(lat[0] - 60) < 1e-6 and abs(lat[-1] + 60) < 1e-6, (lat[0], lat[-1])
    assert worst <= WS_SCALE / 2 + 1e-4, f"解码误差 {worst}"
    return {"verify_max_abs_err": round(worst, 6), "verify_steps": 3}


def extract_tar(task, root, man, stats):
    """保留原 tar，另把成员按原文件名（basename）解到同目录；成员逐个记 manifest。"""
    src = os.path.join(root, task["rel"])
    d = os.path.dirname(src)
    mark = task["rel"] + "#extracted"
    if man.have(mark):
        return
    members = []
    with tarfile.open(src, "r:*") as tf:
        for m in tf:
            if not m.isfile():
                continue
            base = os.path.basename(m.name)
            if not base or base.startswith(".") or (task["extract"] == "nc" and not base.endswith(".nc")):
                continue
            rel = f"{os.path.relpath(d, root)}/{base}"
            if rel in members:
                raise ValueError(f"tar 内重名 {base}")
            dst = os.path.join(root, rel)
            h = hashlib.sha256()
            with tf.extractfile(m) as fi, open(dst + ".part", "wb") as fo:
                while True:
                    b = fi.read(8 << 20)
                    if not b:
                        break
                    fo.write(b)
                    h.update(b)
            os.replace(dst + ".part", dst)
            man.add(source=task["source"] + "_member", url=task["url"], path=rel, member=m.name, bytes=m.size,
                    sha256=h.hexdigest(), status="ok")
            members.append(rel)
    man.add(source=task["source"] + "_extract", url=task["url"], path=mark, status="done", members=members)
    print(log(f"{task['rel']} 解出 {len(members)} 个成员"), flush=True)


# ======================================================================== 编排

class Runner:
    def __init__(self, root, verify):
        self.root, self.verify = root, verify
        self.man = Manifest(root)
        self.pools = {k: cf.ThreadPoolExecutor(max_workers=n, thread_name_prefix=k)
                      for k, n in (("ncei", N_NCEI), ("s3", N_S3), ("hycom", N_HYCOM), ("misc", N_MISC),
                                   ("era5", 1), ("conv", 1))}
        self.failures = []
        self.flock = threading.Lock()
        self.ready_status = {}

    def fail(self, what, e):
        with self.flock:
            self.failures.append((what, f"{type(e).__name__}: {e}"))
            n = len(self.failures)
        print(log(f"失败 {what}：{type(e).__name__}: {e}"), flush=True)
        if n >= MAX_TOTAL_FAIL:
            print(log(f"累计失败 {n} ≥ {MAX_TOTAL_FAIL}，放弃（可原样重提续跑）"), flush=True)
            os._exit(2)

    def run_task(self, task):
        try:
            if task["source"] in ("cmorph", "oisst"):
                return get_whole(task, self.root, self.man, self.stats, allow_404=task["source"] == "cmorph")
            if task["source"] == "hycom":
                return get_whole(task, self.root, self.man, self.stats, validate=validate_hycom)
            if task["source"] in ("watson", "glodap"):
                if not self.man.have(task["rel"], self.verify):
                    get_whole(task, self.root, self.man, self.stats)
                extract_tar(task, self.root, self.man, self.stats)
                return True
        except Exception as e:
            self.fail(task["rel"], e)
            return None

    def era5_download(self, et):
        """下 u/v 原文件（已在 manifest 就跳过），返回是否齐。"""
        ok = True
        for k in ("u", "v"):
            t = et["raw"][k]
            if self.man.have(t["rel"], self.verify):
                continue
            try:
                s3_download(t, self.root, self.man, self.stats, self.pools["s3"])
            except Exception as e:
                self.fail(t["rel"], e)
                ok = False
        return ok

    def era5_convert(self, et):
        ym = et["ym"]
        try:
            raw = {k: self.man.recs[et["raw"][k]["rel"]] for k in ("u", "v")}
            y, m = divmod(ym, 100)
            t0 = time.monotonic()
            chk = derive_ws10(*(os.path.join(self.root, et["raw"][k]["rel"]) for k in ("u", "v")),
                              os.path.join(self.root, et["rel"]), expect_nt=24 * len(days_of_month(y, m)),
                              meta={"inputs": [{"url": r["url"], "sha256": r["sha256"], "bytes": r["bytes"]} for r in raw.values()]})
            out = os.path.join(self.root, et["rel"])
            self.man.add(source="era5_ws10", path=et["rel"], bytes=os.path.getsize(out), sha256=sha256_file(out),
                         status="ok", inputs=[{"url": r["url"], "path": r["path"], "bytes": r["bytes"], "sha256": r["sha256"]}
                                              for r in raw.values()], check=chk)
            # [repro] ERA5 u/v 原件保留在 era5_raw/（原始数据一律保留）：原为逐个 os.remove 并记 deleted_after_derive；
            # [repro] 现不删、不记删除行；raw_sync.py A6 把它们复制进 bulk raw 档
            self.stats.file_done()
            print(log(f"[ERA5] {ym} 转换完成 {time.monotonic() - t0:.0f} s，{chk}"), flush=True)
            return True
        except Exception as e:
            traceback.print_exc()
            self.fail(et["rel"], e)
            return False

    def write_ready(self, label, tasks, eras, statics):
        rel = f"month_ready/{label}.json"
        p = os.path.join(self.root, rel)
        missing, notok = [], []
        for t in tasks:
            r = self.man.recs.get(t["rel"])
            if not r or not self.man.have(t["rel"]):
                notok.append(t["rel"])
            elif r.get("status") == "missing":
                missing.append(t["rel"])
        for et in eras:
            if not self.man.have(et["rel"]):
                notok.append(et["rel"])
        for s in statics:
            if not self.man.have(s["rel"]) or not self.man.have(s["rel"] + "#extracted"):
                notok.append(s["rel"])
        if notok:
            self.ready_status[label] = f"未齐 {len(notok)}（例 {notok[:3]}）"
            print(log(f"[月] {label} 未齐：{len(notok)} 项未完成，例 {notok[:3]}"), flush=True)
            return False
        by = {}
        for t in tasks:
            by.setdefault(t["source"], 0)
            by[t["source"]] += 1
        nxt = label + 1 if label % 100 < 12 else None
        info = {"month": label, "written_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                "version": VERSION, "counts": by, "cmorph_missing": missing,
                "era5_ws10": [et["rel"] for et in eras], "static": [s["rel"] for s in statics],
                "spinup_day": str(SPINUP_DAY) if label == YEAR * 100 + 1 else None,
                "era5_next_month_file": f"era5_ws10/{nxt}.nc" if nxt else None,
                "note": "月末最后一个 HH:30 需要下月第一个整点风：读 era5_next_month_file（其 month_ready 未写时文件可能尚不存在）；"
                        "2000-12 的下月（200101）不在下载范围"}
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p + ".tmp", "w", encoding="utf-8") as f:
            json.dump(info, f, ensure_ascii=False, indent=1)
        os.replace(p + ".tmp", p)
        self.ready_status[label] = "ok"
        print(self.stats.line(f"月 {label} 齐"), flush=True)
        return True

    def main(self):
        groups = build_groups()
        plan = []
        for label, days, eras, statics in groups:
            tasks = [cmorph_task(d, hh) for d in days for hh in range(24)]
            tasks += [oisst_task(d) for d in days] + [hycom_task(d) for d in days]
            plan.append((label, tasks, [era5_task(ym) for ym in eras], statics))
        todo_est, n_all, n_skip = 0.0, 0, 0
        for label, tasks, eras, statics in plan:
            for t in tasks + statics:
                n_all += 1
                if self.man.have(t["rel"], self.verify):
                    n_skip += 1
                else:
                    todo_est += EST[t["source"]]
            for et in eras:  # 每月计 3 个文件：u、v 原文件＋ws10
                n_all += 3
                if self.man.have(et["rel"]):
                    n_skip += 3
                else:
                    raw_todo = [k for k in ("u", "v") if not self.man.have(et["raw"][k]["rel"])]
                    n_skip += 2 - len(raw_todo)
                    todo_est += EST["era5_raw"] * len(raw_todo)
        self.stats = Stats(todo_est)
        self.stats.total_files = n_all - n_skip
        print(log(f"[P3] {VERSION}：共 {n_all} 个文件（ERA5 每月计 u、v、ws10 三个），已有 {n_skip}，待下约 {todo_est / 1e9:.1f} GB；"
                  f"根 {self.root}"), flush=True)
        stop = threading.Event()

        def heartbeat():
            while not stop.wait(60):
                if time.monotonic() - self.stats.last_print >= HEARTBEAT_S:
                    print(self.stats.line("心跳"), flush=True)
        threading.Thread(target=heartbeat, daemon=True).start()
        t0 = time.monotonic()
        conv_futs = []
        for label, tasks, eras, statics in plan:
            self.stats.group = str(label)
            futs = []
            for t in statics:
                futs.append(self.pools["misc"].submit(self.run_task, t))
            for t in tasks:
                if self.man.have(t["rel"], self.verify):
                    continue
                pool = self.pools["hycom" if t["source"] == "hycom" else "ncei"]
                futs.append(pool.submit(self.run_task, t))
            era_futs = [(et, self.pools["era5"].submit(self.era5_download, et)) for et in eras if not self.man.have(et["rel"])]
            cf.wait(futs + [f for _, f in era_futs])
            for et, f in era_futs:
                if f.result():
                    conv_futs.append(self.pools["conv"].submit(self.era5_convert, et))
            conv_futs.append(self.pools["conv"].submit(self.write_ready, label, tasks, eras, statics))
        cf.wait(conv_futs)
        # 末轮：串行补试失败项（与 P2 同）
        if self.failures:
            print(log(f"[P3] 首轮失败 {len(self.failures)} 项，串行补试"), flush=True)
            self.failures = []
            for label, tasks, eras, statics in plan:
                if self.ready_status.get(label) == "ok":
                    continue
                for t in statics + [t for t in tasks if not self.man.have(t["rel"])]:
                    self.run_task(t)
                for et in eras:
                    if not self.man.have(et["rel"]) and self.era5_download(et):
                        self.era5_convert(et)
                self.write_ready(label, tasks, eras, statics)
        stop.set()
        el = time.monotonic() - t0
        summary = self.summarize(plan, el)
        out = os.environ.get("REPRO_OUTPUT_DIR")
        ok = all(v == "ok" for v in self.ready_status.values()) and len(self.ready_status) == len(plan)
        if ok:
            with open(os.path.join(self.root, "fetch_done.json"), "w", encoding="utf-8") as f:
                json.dump(summary, f, ensure_ascii=False, indent=1)
        if out:
            with open(os.path.join(out, "fetch_summary.json" if ok else "fetch_status.json"), "w", encoding="utf-8") as f:
                json.dump({**summary, "ready_status": self.ready_status, "failures": self.failures}, f,
                          ensure_ascii=False, indent=1)
        print(self.stats.line("结束"), flush=True)
        print(log(f"[P3] {'全部完成' if ok else '未完成'}：{json.dumps(self.ready_status, ensure_ascii=False)}"), flush=True)
        for p in self.pools.values():
            p.shutdown(wait=False)
        return 0 if ok else 2

    def summarize(self, plan, el):
        by = {}
        for rel, r in self.man.recs.items():
            s = r.get("source", "?")
            x = by.setdefault(s, {"files": 0, "bytes": 0, "missing": 0, "other": 0})
            if r.get("status") == "ok":
                x["files"] += 1
                x["bytes"] += r.get("bytes", 0)
            elif r.get("status") == "missing":
                x["missing"] += 1
            else:
                x["other"] += 1
        return {"version": VERSION, "finished_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                "seconds_this_run": round(el, 1), "bytes_this_run": self.stats.bytes, "by_source": by,
                "total_bytes_on_disk": sum(v["bytes"] for k, v in by.items() if k != "era5_raw"),
                "months": [p[0] for p in plan]}


# ======================================================================== 自测（不联网）

def selftest(workdir):
    import numpy as np, netCDF4
    os.makedirs(workdir, exist_ok=True)
    res = {}
    groups = build_groups()
    res["counts"] = {"cmorph": sum(len(g[1]) * 24 for g in groups), "oisst": sum(len(g[1]) for g in groups),
                     "hycom": sum(len(g[1]) for g in groups), "era5_months": sum(len(g[2]) for g in groups)}
    assert res["counts"] == {"cmorph": 8808, "oisst": 367, "hycom": 367, "era5_months": 13}, res["counts"]
    res["urls"] = [cmorph_task(SPINUP_DAY, 23)["url"], oisst_task(dt.date(2000, 2, 29))["url"],
                   hycom_task(dt.date(2000, 1, 1))["url"], era5_task(200002)["raw"]["v"]["url"],
                   era5_task(199912)["raw"]["u"]["url"]]
    assert res["urls"][3].endswith("2000020100_2000022923.nc") and res["urls"][4].endswith("1999120100_1999123123.nc")
    # manifest 跳过逻辑
    man = Manifest(workdir)
    p = os.path.join(workdir, "x/a.bin")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "wb").write(b"abc")
    man.add(source="t", url="u", path="x/a.bin", bytes=3, sha256=hashlib.sha256(b"abc").hexdigest(), status="ok")
    man2 = Manifest(workdir)
    assert man2.have("x/a.bin", verify=True)
    open(p, "wb").write(b"abcd")
    assert not man2.have("x/a.bin")
    # 合成 ERA5（维度、降序纬度、分块仿 NCAR 文件）
    nt = 54
    lat = np.linspace(90, -90, 721)
    lon = np.arange(1440) * 0.25
    rng = np.random.default_rng(0)
    paths = {}
    for k, name in (("u", "VAR_10U"), ("v", "VAR_10V")):
        paths[k] = os.path.join(workdir, f"syn_{k}.nc")
        with netCDF4.Dataset(paths[k], "w") as ds:
            ds.createDimension("time", nt)
            ds.createDimension("latitude", 721)
            ds.createDimension("longitude", 1440)
            t = ds.createVariable("time", "i4", ("time",))
            t.units = "hours since 1900-01-01 00:00:00"
            t.calendar = "gregorian"
            t[:] = 876576 + np.arange(nt)
            y = ds.createVariable("latitude", "f8", ("latitude",))
            y.units = "degrees_north"
            y[:] = lat
            x = ds.createVariable("longitude", "f8", ("longitude",))
            x[:] = lon
            v = ds.createVariable(name, "f4", ("time", "latitude", "longitude"), zlib=True, complevel=1,
                                  chunksizes=(27, 139, 277))
            for b in range(0, nt, 27):
                v[b:b + 27] = rng.standard_normal((min(27, nt - b), 721, 1440), dtype=np.float32) * np.float32(7)
    out = os.path.join(workdir, "syn_ws10.nc")
    res["rss_before_derive_mib"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 1)
    chk = derive_ws10(paths["u"], paths["v"], out, expect_nt=nt, meta={"inputs": ["synthetic"]})
    with netCDF4.Dataset(out) as ds:
        w = ds.variables["ws10"]
        res["ws10"] = {"dims": list(w.dimensions), "shape": list(w.shape), "dtype": str(w.dtype),
                       "scale_factor": float(w.scale_factor), "add_offset": float(w.add_offset), "units": w.units,
                       "lat0": float(ds.variables["latitude"][0]), "lat_last": float(ds.variables["latitude"][-1]),
                       "time_units": ds.variables["time"].units, "bytes": os.path.getsize(out)}
    assert res["ws10"]["shape"] == [nt, 481, 1440] and res["ws10"]["dtype"] == "int16"
    res["derive_check"] = chk
    res["peak_rss_mib"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 1)  # macOS 为字节
    for f in os.listdir(workdir):
        fp = os.path.join(workdir, f)
        shutil.rmtree(fp) if os.path.isdir(fp) else os.remove(fp)
    print(json.dumps(res, ensure_ascii=False, indent=1), flush=True)
    return res


def main():
    global LOGF
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", help="p3-cache 根目录（绝对路径）")
    ap.add_argument("--verify", action="store_true", help="跳过前对已有文件重算 sha256")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--min-free-gb", type=float, default=80.0)
    a = ap.parse_args()
    out = os.environ.get("REPRO_OUTPUT_DIR")
    if a.selftest:
        res = selftest(os.path.join(out or ".", "selftest_tmp"))
        if out:
            with open(os.path.join(out, "p3_fetch_selftest.json"), "w", encoding="utf-8") as f:
                json.dump(res, f, ensure_ascii=False, indent=1)
        return 0
    if not a.cache or not os.path.isabs(a.cache) or not os.path.isdir(a.cache):
        print("需要已存在的绝对路径 --cache", flush=True)
        return 3
    free = shutil.disk_usage(a.cache).free / 1e9
    if free < a.min_free_gb:
        print(f"磁盘空闲 {free:.0f} GB < {a.min_free_gb} GB，停止", flush=True)
        return 3
    if out:
        LOGF = open(os.path.join(out, "fetch_log.txt"), "a", encoding="utf-8")
    print(log(f"[P3] 空闲 {free:.0f} GB，python {sys.version.split()[0]}，argv {sys.argv[1:]}"), flush=True)
    return Runner(a.cache, a.verify).main()


if __name__ == "__main__":
    sys.exit(main())
