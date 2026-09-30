#!/usr/bin/env python3
"""P1：系泊雨事件计数、事件/对照 ΔS(1 m) 与 0.5 m/1 m 一致性，判 D1/D2/D5。

设计：事件定义、P1 步骤与 D1/D2/D5 阈值在运行前写定（D4 仅作信息量输出）。
  - D3、R_RIM/R_S20/P5/β、CMORPH、敏感性分析均属 P2，本脚本不做。
  - 「只数不看效应」：不计算 ΔpCO₂、ΔS(5 m)、任何模型比值；5 m 与 pCO₂ 只记可用性。

用法：
  python3 p1_events.py                      # 全量；输出到 $REPRO_OUTPUT_DIR
  python3 p1_events.py --out DIR            # 本地覆盖输出目录
  python3 p1_events.py --smoke --out DIR    # 冒烟：TAO165E 单站 2016-12 一个月，只验证端到端产物格式
  python3 p1_events.py --stations TAO165E,BOBOA   # 调试用子集（门判定仍照常输出，但不代表全量）

依赖：仅 Python 3.12 标准库。数据全部匿名 HTTP：
  - MAPCO2 3 小时序列：PMEL ERDDAP tabledap `pmel_co2_moorings_*`（CSV，按站×年分块）
  - GTMBA 高分辨率雨量/风/温盐：NDBC THREDDS OPeNDAP ASCII（OceanSITES，按部署文件分块）
      TAO 站：dods.ndbc.noaa.gov/thredds/dodsC/oceansites/DATA/<site>/OS_<site>_<dep>_<D|M>_{RAIN_10min,WIND_10min,SALT_hourly|SALT_10min}.nc
      BOBOA（RAMA 15N90E）：.../DATA_GRIDDED/RAMA/OS_15n90e_<yyyymm>_<D|M>_{TM_10m,SBP_hr,S_hr}.nc
  缓存：<out>/cache/（gzip 存盘，已下载跳过，.part 原子改名，支持断点续跑）

产物（<out>/）：p1_summary.json、p1_events.csv、p1_log.txt。退出码：0＝跑完（无论门过不过）；2＝数据源故障；3＝其他异常。

实现选择清单（设计未规定或有歧义之处的实现层决定；均未放宽任何事先写定的门）：
  I1 雨量源：GTMBA 10 分钟雨量取 NDBC OceanSITES 延时（D）/混合（M）模式文件；实时（R）文件默认不用
     （`--include-realtime` 可加入作敏感性，非默认、不作主判）。NDBC 目录对部分年份缺文件
     （如 T0N165E 2010–2014 无 RAIN），覆盖可能低于按日值估计的计数。
  I2 质量码：设计「GTMBA 质量码 1–2」按 OceanSITES QC 字面值 {1,2} 实现（good / probably good）；
     OceanSITES 的 5（adjusted / value_changed）一律剔除，剔除行数写进日志与 summary。
  I3 MAPCO2 质量标志：ERDDAP 数据集无 QC 列，把非 NaN 值视为 PMEL 终版 QC（=标志 2）后的数据【推测，未验证】。
  I4 小时化：UTC 小时箱 [HH:00, HH+1:00)，时间戳向下取整；10 分钟量取小时均值，需 ≥50% 子样本有效（≥3/6）；
     小时量（温盐）取该小时点值；多个部署文件同一小时重叠时取平均。MAPCO2 HH:17 样本归入 HH 箱。
  I5 「有雨/无雨」：小时均雨强 ≥0.4 mm/h 为有雨小时，<0.4 为无雨（GTMBA 官方雨量计精度 ±0.4 mm/h，Serra 2001）；
     设计未给无雨门限，10 分钟原始值噪声含大量负值，逐字「=0」不可实现。
  I6 事件起点 t0：一个有效有雨小时，且其前 24 h 全部有效且无雨（实现「与前一事件至少 24 h 无雨」，按更严格的
     「前 24 h 完全无雨」执行）；累积量＝[t0, t0+24h) 内有雨小时之和 ≥10 mm，且这 24 h 雨量全部有效。
  I7 数据要求「[−6 h, +24 h] 内有 …」：对 0.5 m SSS、0.5 m pCO₂、1 m S、风，各自在 [t0−6,t0)、[t0,t0+6)、
     [t0+6,t0+24) 三段中每段至少 1 个有效值。5 m S 可选，只记可用性。
  I8 ΔS(z, 0–6 h)＝[t0, t0+6h) 内有效值均值 − [t0−6h, t0) 内有效值中位数。
  I9 对照窗：候选起点 t0±k·24h（k=1..30，同 UTC 小时即同当地时刻），需满足：前 24 h 全有效且无雨、与事件相同的
     数据要求、其 [c−6, c+24) 不与本站任何 ≥10 mm 雨事件的 [e−6, e+24) 重叠；按 (|k|, k) 取最近 5 个。
     对照可被不同事件共用。D2 主判只用配满 5 个对照的事件；「≥1 个对照」口径只作敏感性输出。
  I10 D2 的 95% CI：按「站×季节(DJF/MAM/JJA/SON，不分年)」整簇 bootstrap（事件与其对照同簇重抽），
     B=10000，seed=20260926，百分位区间（与 P2 的 cluster 口径一致）；Welch 近似 CI 只作信息量。
  I11 D5：每事件的「雨前中位差」＝median(SSS0.5, 前 6 h) − median(S1, 前 6 h)，取全部事件的中位数；
     「事件异常相关」＝各事件 ΔS(0.5 m) 与 ΔS(1 m) 的 Pearson 相关（全站合并）。逐站值只作信息量。
     设计冲突：一处说 D5 不过→弃用 0.5 m；另一处（出口 3）说 D5 不过→停止。脚本只报数，不在此裁定。
     为「弃用 0.5 m」分支另算「不要求 0.5 m」的事件集（只要求 1 m S＋风），仅作信息量，不参与主判。
  I12 型态归类：TAO165E＝西太暖池；TAO8S165E＝SPCZ；BOBOA＝孟加拉湾季风；TAO170W/155W/140W/125W/110W＝中东太平洋。
  I13 站点配对：MAPCO2「110w0」（ERDDAP 标称 1°N,109°W）配 GTMBA T0N110W；BOBOA 配 RAMA 15N90E。
  I14 部署文件筛选：TAO 部署起始日 ∈ [MAPCO2 起 − 550 天, MAPCO2 止]；RAMA 网格化文件为多年合并，只要求起始
     ≤ MAPCO2 止；站时间轴只保留 MAPCO2 覆盖期 ±3 天。
  I17 诊断（信息量）：同一小时 MAPCO2 SSS 与 TAO 1 m S 之差的分布。冒烟中 TAO165E 两者几乎逐值相同
     （差 ≤0.003），MAPCO2「0.5 m SSS」可能直接取自浮标 1 m CTD【推测】——若成立，D5 会平凡通过、0.5 m 层不独立。
  I15 D4（ΔS(0.5 m) ≤ −0.2 psu 事件数）设计未列入 P1 判定，只作信息量输出。
  I16 统计不可算（事件/簇过少）时，该门记为不过（not_evaluable），不视为通过。

Change Log:
  2026-09-26 初版。
  2026-09-26 build_events 增可选参数 ctrl_ok（默认 None，P1 行为与产物不变），供 p1b_extend.py 的 R 组排除触及 P1 时段的对照。
"""

import argparse
import csv
import gzip
import http.client
import json
import math
import os
import random
import re
import shutil
import statistics
import sys
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
from array import array
from bisect import bisect_left
from datetime import datetime, timezone

VERSION = "p1-2026-09-26a"

# ---- 事先写定的参数（不得改） ----
RAIN_EVENT_MM = 10.0
DRY_BEFORE_H = 24
PRE_H = 6
POST_H = 24
DS_WIN_H = 6
CTRL_PER_EVENT = 5
CTRL_MAX_DAYS = 30
D1_MIN_TOTAL = 60
D1_MIN_EVENTS_PER_STATION = 10
D1_MIN_STATIONS = 3
D2_MIN_ABS_DIFF = 0.03
D4_DS_THRESH = -0.2
D4_MIN_EVENTS = 30
D5_MAX_ABS_MEDDIFF = 0.05
D5_MIN_CORR = 0.6
QC_OK = {1, 2}

# ---- 实现层参数（见 docstring 实现选择） ----
RAIN_HOUR_MM = 0.4          # I5
SUBHOURLY_MIN_FRAC = 0.5    # I4
BOOT_B = 10000              # I10
BOOT_SEED = 20260926
DEPLOY_LOOKBACK_DAYS = 550  # I14
PAD_H = 72                  # I14

ERDDAP = "https://data.pmel.noaa.gov/pmel/erddap"
THREDDS = "https://dods.ndbc.noaa.gov/thredds"
UA = "rain-co2-dilution-replication-p1/1.0 (research; python-urllib)"

STATIONS = [
    {"name": "TAO165E", "regime": "西太暖池", "erddap": "pmel_co2_moorings_c993_7f05_7d40",
     "src": "tao", "site": "T0N165E", "lon": 165.0},
    {"name": "TAO8S165E", "regime": "SPCZ", "erddap": "pmel_co2_moorings_8a73_1677_822e",
     "src": "tao", "site": "T8S165E", "lon": 165.0},
    {"name": "BOBOA", "regime": "孟加拉湾季风", "erddap": "pmel_co2_moorings_ee12_e5fb_4537",
     "src": "rama", "site": "15n90e", "lon": 90.0},
    {"name": "TAO170W", "regime": "中东太平洋", "erddap": "pmel_co2_moorings_b1d3_602a_9da4",
     "src": "tao", "site": "T0N170W", "lon": -170.0},
    {"name": "TAO155W", "regime": "中东太平洋", "erddap": "pmel_co2_moorings_bf98_0a93_1a18",
     "src": "tao", "site": "T0N155W", "lon": -155.0},
    {"name": "TAO140W", "regime": "中东太平洋", "erddap": "pmel_co2_moorings_509b_0e87_dbc3",
     "src": "tao", "site": "T0N140W", "lon": -140.0},
    {"name": "TAO125W", "regime": "中东太平洋", "erddap": "pmel_co2_moorings_af32_4443_213a",
     "src": "tao", "site": "T0N125W", "lon": -125.0},
    {"name": "TAO110W", "regime": "中东太平洋", "erddap": "pmel_co2_moorings_ca69_976d_4677",
     "src": "tao", "site": "T0N110W", "lon": -110.0},
]

SMOKE_STATION = "TAO165E"
SMOKE_MONTH = "2016-12"

NAN = float("nan")
EPOCH1950_H = 7305 * 24  # 1950-01-01 → 1970-01-01 的小时数

TAO_FILE_RE = re.compile(
    r"^OS_(?P<site>T\w+?)_(?P<dep>[A-Z]{2}\d{3}[A-Z])-(?P<date>\d{8})_(?P<mode>[DMR])_"
    r"(?P<var>RAIN|WIND|SALT)_(?P<res>10min|hourly)\.nc$")
RAMA_FILE_RE = re.compile(
    r"^OS_(?P<site>[0-9.]+[ns][0-9.]+[ew])_(?P<ym>\d{6})_(?P<mode>[DMR])_(?P<kind>TM_10m|SBP_hr|S_hr)\.nc$")


class FetchError(Exception):
    pass


# ======================================================================== 日志
class Log:
    def __init__(self, path):
        self.f = open(path, "a", encoding="utf-8")

    def log(self, msg, echo=False):
        ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        self.f.write(f"{ts} {msg}\n")
        self.f.flush()
        if echo:
            print(msg, flush=True)

    def close(self):
        self.f.close()


# ======================================================================== HTTP
class Http:
    def __init__(self, log, min_interval=1.0, timeout=180, max_attempts=6):
        self.log = log
        self.min_interval = min_interval
        self.timeout = timeout
        self.max_attempts = max_attempts
        self._last = 0.0
        self.n_requests = 0
        self.n_cached = 0
        self.bytes = 0

    def _pace(self):
        wait = self.min_interval - (time.monotonic() - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.monotonic()

    def get(self, url, dest, empty_ok=False):
        """下载到 dest（gzip 存盘）。返回 True＝有数据，False＝ERDDAP 明确返回无结果。"""
        if os.path.exists(dest):
            self.n_cached += 1
            return True
        if os.path.exists(dest + ".empty"):
            self.n_cached += 1
            return False
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        tmp = dest + ".part"
        last_err = None
        for attempt in range(1, self.max_attempts + 1):
            self._pace()
            self.n_requests += 1
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept-Encoding": "gzip"})
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    enc = (r.headers.get("Content-Encoding") or "").lower()
                    with open(tmp, "wb") as f:
                        if enc == "gzip":
                            shutil.copyfileobj(r, f, 1 << 20)
                        else:
                            with gzip.GzipFile(fileobj=f, mode="wb") as gz:
                                shutil.copyfileobj(r, gz, 1 << 20)
                with gzip.open(tmp, "rb") as g:  # 完整性校验（截断的 gzip 会抛异常）
                    while g.read(1 << 20):
                        pass
                self.bytes += os.path.getsize(tmp)
                os.replace(tmp, dest)
                return True
            except urllib.error.HTTPError as e:
                body = _read_err_body(e)
                if empty_ok and e.code == 404 and "no matching results" in body:
                    open(dest + ".empty", "w").close()
                    return False
                last_err = f"HTTP {e.code}: {body[:200]!r}"
                if e.code not in (408, 429, 500, 502, 503, 504):
                    raise FetchError(f"{url} → {last_err}")
            except (urllib.error.URLError, http.client.HTTPException, TimeoutError,
                    ConnectionError, OSError, EOFError) as e:
                last_err = f"{type(e).__name__}: {e}"
            if os.path.exists(tmp):
                os.remove(tmp)
            backoff = min(300, 5 * 2 ** (attempt - 1))
            self.log.log(f"retry {attempt}/{self.max_attempts} in {backoff}s: {url} ({last_err})")
            time.sleep(backoff)
        raise FetchError(f"{url} → 重试 {self.max_attempts} 次仍失败：{last_err}")


def _read_err_body(e):
    try:
        raw = e.read()
    except Exception:
        return ""
    try:
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        return raw.decode("utf-8", "replace")
    except Exception:
        return ""


def read_gz_text(path):
    with gzip.open(path, "rt", encoding="latin-1") as f:
        return f.read()


# ======================================================================== 小时序列
class Hourly:
    """站时间轴上的小时值：sum/count 累加（多部署重叠取平均）。"""

    def __init__(self, h0, n):
        self.h0 = h0
        self.n = n
        self.s = array("d", bytes(8 * n))
        self.c = array("H", bytes(2 * n))

    def add(self, h, v):
        i = h - self.h0
        if 0 <= i < self.n:
            self.s[i] += v
            if self.c[i] < 65535:
                self.c[i] += 1

    def values(self):
        out = array("d", [NAN]) * self.n
        s, c = self.s, self.c
        for i in range(self.n):
            if c[i]:
                out[i] = s[i] / c[i]
        return out

    def n_valid(self):
        return sum(1 for x in self.c if x)


def days1950_to_hour(t):
    return math.floor(t * 24.0 + 1e-6) - EPOCH1950_H


def iso_to_hour(s):
    dt = datetime.strptime(s[:19], "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc)
    return int(dt.timestamp() // 3600)


def hour_to_iso(h):
    return datetime.fromtimestamp(h * 3600, tz=timezone.utc).strftime("%Y-%m-%dT%H:%MZ")


def hour_month(h):
    return datetime.fromtimestamp(h * 3600, tz=timezone.utc).month


def season_of(h):
    m = hour_month(h)
    return {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
            6: "JJA", 7: "JJA", 8: "JJA"}.get(m, "SON")


# ======================================================================== OPeNDAP ASCII
DDS_VAR_RE = re.compile(r"\b(?:Float32|Float64|Byte|Int16|Int32|UInt16|UInt32)\s+(\w+)\[")
ASCII_HDR_RE = re.compile(r"^([A-Za-z_][\w.]*)((?:\[\d+\])+)\s*$")


def parse_dap_ascii(path, want):
    """流式解析 OPeNDAP .ascii（gzip 缓存）。返回 {短变量名: array('d')}，按 C 序展平。"""
    out = {}
    cur = None
    in_data = False
    with gzip.open(path, "rt", encoding="latin-1") as f:
        for line in f:
            if not in_data:
                if line.startswith("-----"):
                    in_data = True
                continue
            line = line.rstrip("\r\n")
            if not line.strip():
                cur = None
                continue
            m = ASCII_HDR_RE.match(line)
            if m:
                name = m.group(1).split(".")[-1]
                if name in want and name not in out:
                    cur = array("d")
                    out[name] = cur
                else:
                    cur = None
                continue
            if cur is None:
                continue
            if line.startswith("["):
                line = line[line.rfind("]") + 1:]
            for tok in line.split(","):
                tok = tok.strip()
                if tok:
                    cur.append(float(tok))
    return out


def file_to_hourly(times, vals, qcs, per_hour, qc_drop, lo, hi):
    """单文件：把 (time_days1950, value, qc) 聚成小时均值。per_hour＝每小时应有样本数。"""
    need = max(1, math.ceil(per_hour * SUBHOURLY_MIN_FRAC))
    res = []
    cur_h = None
    acc = 0.0
    cnt = 0
    for t, v, q in zip(times, vals, qcs):
        h = days1950_to_hour(t)
        if h != cur_h:
            if cur_h is not None and cnt >= need:
                res.append((cur_h, acc / cnt))
            cur_h, acc, cnt = h, 0.0, 0
        qi = int(q) if q == q else -1
        if qi not in QC_OK:
            if qi in qc_drop:
                qc_drop[qi] += 1
            else:
                qc_drop[qi] = 1
            continue
        if not (lo <= v <= hi):
            continue
        acc += v
        cnt += 1
    if cur_h is not None and cnt >= need:
        res.append((cur_h, acc / cnt))
    return res


# ======================================================================== 数据获取
def fetch_mapco2(http, st, cache, log, smoke_month=None):
    """返回 (h_start, h_end, rows[(hour, SSS, pCO2)])；按年分块下载。"""
    info_path = os.path.join(cache, "erddap", st["erddap"] + "_info.csv.gz")
    http.get(f"{ERDDAP}/info/{st['erddap']}/index.csv", info_path)
    tstart = tend = None
    for row in csv.reader(read_gz_text(info_path).splitlines()):
        if len(row) >= 5 and row[0] == "attribute" and row[1] == "NC_GLOBAL":
            if row[2] == "time_coverage_start":
                tstart = row[4]
            elif row[2] == "time_coverage_end":
                tend = row[4]
    if not tstart or not tend:
        raise FetchError(f"{st['erddap']} 无 time_coverage 元数据")
    if smoke_month:
        y, m = map(int, smoke_month.split("-"))
        ny, nm = (y + 1, 1) if m == 12 else (y, m + 1)
        chunks = [(f"{y:04d}-{m:02d}-01T00:00:00Z", f"{ny:04d}-{nm:02d}-01T00:00:00Z", f"{smoke_month}")]
    else:
        y0, y1 = int(tstart[:4]), int(tend[:4])
        chunks = [(f"{y}-01-01T00:00:00Z", f"{y + 1}-01-01T00:00:00Z", str(y)) for y in range(y0, y1 + 1)]
    rows = []
    for a, b, tag in chunks:
        q = "time,SST,SSS,pCO2_sw,pCO2_air&time%3E=" + a + "&time%3C" + b
        path = os.path.join(cache, "erddap", f"{st['erddap']}_{tag}.csv.gz")
        if not http.get(f"{ERDDAP}/tabledap/{st['erddap']}.csv?{q}", path, empty_ok=True):
            continue
        with gzip.open(path, "rt", encoding="utf-8") as f:
            rd = csv.reader(f)
            hdr = next(rd)
            next(rd)  # units 行
            it, isss, ipc = hdr.index("time"), hdr.index("SSS"), hdr.index("pCO2_sw")
            for r in rd:
                h = iso_to_hour(r[it])
                sss = _num(r[isss])
                pc = _num(r[ipc])
                if sss is not None and not (0.0 < sss < 45.0):
                    sss = None
                if pc is not None and not (50.0 < pc < 2000.0):
                    pc = None
                rows.append((h, sss, pc))
    if smoke_month:
        hs = iso_to_hour(chunks[0][0])
        he = iso_to_hour(chunks[0][1]) - 1
    else:
        hs, he = iso_to_hour(tstart), iso_to_hour(tend)
    log.log(f"[{st['name']}] MAPCO2 {tstart}..{tend}：{len(rows)} 行（{len(chunks)} 块）")
    return hs, he, rows


def _num(s):
    try:
        v = float(s)
    except ValueError:
        return None
    return v if v == v else None


def list_gtmba_files(http, st, cache, hs, he, smoke, include_rt=False):
    """按站列出要用的 OceanSITES 文件：[(url_base, kind, per_hour_hint)]。"""
    if st["src"] == "tao":
        cat_url = f"{THREDDS}/catalog/oceansites/DATA/{st['site']}/catalog.xml"
        dods = f"{THREDDS}/dodsC/oceansites/DATA/{st['site']}/"
    else:
        cat_url = f"{THREDDS}/catalog/oceansites/DATA_GRIDDED/RAMA/catalog.xml"
        dods = f"{THREDDS}/dodsC/oceansites/DATA_GRIDDED/RAMA/"
    cpath = os.path.join(cache, "thredds", f"catalog_{st['src']}_{st['site']}.xml.gz")
    http.get(cat_url, cpath)
    names = sorted(set(re.findall(r'name="(OS_[^"]+\.nc)"', read_gz_text(cpath))))
    lo_h = hs - DEPLOY_LOOKBACK_DAYS * 24
    picked = {}  # (deployment, kind) -> (mode_rank, name, per_hour, start_h)
    for n in names:
        if st["src"] == "tao":
            m = TAO_FILE_RE.match(n)
            if not m or m.group("site") != st["site"]:
                continue
            start_h = iso_to_hour(f"{m['date'][:4]}-{m['date'][4:6]}-{m['date'][6:]}T00:00:00")
            dep = m.group("dep") + "-" + m.group("date")
            kind = m.group("var")
            per_hour = 6 if m.group("res") == "10min" else 1
            res_rank = 0 if (kind != "SALT" or m.group("res") == "hourly") else 1
        else:
            m = RAMA_FILE_RE.match(n)
            if not m or m.group("site") != st["site"]:
                continue
            start_h = iso_to_hour(f"{m['ym'][:4]}-{m['ym'][4:]}-01T00:00:00")
            dep = m.group("ym")
            kind = {"TM_10m": "TM", "SBP_hr": "SALT", "S_hr": "SALT"}[m.group("kind")]
            per_hour = 6 if kind == "TM" else 1
            res_rank = 0 if m.group("kind") == "SBP_hr" else 1
        if start_h > he or (st["src"] == "tao" and start_h < lo_h):
            continue
        if m.group("mode") == "R" and not include_rt:
            continue
        rank = ({"D": 0, "M": 1, "R": 2}[m.group("mode")], res_rank)
        key = (dep, kind)
        if key not in picked or rank < picked[key][0]:
            picked[key] = (rank, n, per_hour, start_h)
    files = [(dods + v[1], k[1], v[2], v[3]) for k, v in picked.items()]
    files.sort(key=lambda x: (x[3], x[0]))
    if smoke:  # 冒烟：每类只取覆盖该月的最近一个部署
        keep = {}
        for f in files:
            if f[3] <= hs:
                keep[f[1]] = f
        files = list(keep.values())
    return files


def load_gtmba_file(http, cache, base, kind, per_hour, log, qc_drop):
    """下载并解析一个 OceanSITES 文件，返回 {series_key: [(hour, value)]}。"""
    fname = base.rsplit("/", 1)[-1]
    sub = os.path.join(cache, "thredds", fname.split("_")[1])
    dds_path = os.path.join(sub, fname + ".dds.gz")
    http.get(base + ".dds", dds_path)
    present = set(DDS_VAR_RE.findall(read_gz_text(dds_path)))
    wanted = [v for v in ("RAIN", "WSPD", "PSAL") if v in present]
    if kind == "WIND":
        wanted = [v for v in wanted if v == "WSPD"]
    elif kind == "RAIN":
        wanted = [v for v in wanted if v == "RAIN"]
    elif kind == "SALT":
        wanted = [v for v in wanted if v == "PSAL"]
    if not wanted:
        log.log(f"  跳过 {fname}：DDS 中无所需变量")
        return {}
    proj = ["TIME"] + (["DEPTH"] if "PSAL" in wanted else [])
    for v in wanted:
        if v + "_QC" not in present:
            raise FetchError(f"{fname} 缺 {v}_QC")
        proj += [f"{v}.{v}", f"{v}_QC.{v}_QC"]
    apath = os.path.join(sub, fname + ".ascii.gz")
    http.get(base + ".ascii?" + ",".join(proj), apath)
    arr = parse_dap_ascii(apath, set(["TIME", "DEPTH"] + wanted + [v + "_QC" for v in wanted]))
    times = arr.get("TIME")
    if not times:
        raise FetchError(f"{fname} 解析不到 TIME")
    nt = len(times)
    out = {}
    for v in wanted:
        vals, qcs = arr[v], arr[v + "_QC"]
        if v == "PSAL":
            depths = list(arr.get("DEPTH", []))
            nd = len(depths)
            if nd == 0 or len(vals) != nt * nd or len(qcs) != nt * nd:
                raise FetchError(f"{fname} PSAL 形状不符（nt={nt}, nd={nd}, n={len(vals)}）")
            for target, key in ((1.0, "s1"), (5.0, "s5")):
                if target in depths:
                    j = depths.index(target)
                    out[key] = file_to_hourly(times, vals[j::nd], qcs[j::nd], per_hour,
                                              qc_drop.setdefault(key, {}), 0.0, 45.0)
        else:
            if len(vals) != nt or len(qcs) != nt:
                raise FetchError(f"{fname} {v} 长度不符（{len(vals)} vs {nt}）")
            key, lo, hi = ("rain", -50.0, 500.0) if v == "RAIN" else ("wind", 0.0, 60.0)
            out[key] = file_to_hourly(times, vals, qcs, per_hour, qc_drop.setdefault(key, {}), lo, hi)
    log.log(f"  {fname}: nt={nt} → " + ", ".join(f"{k}:{len(v)}h" for k, v in out.items()))
    return out


# ======================================================================== 事件与对照
def nanmedian(xs):
    xs = [x for x in xs if x == x]
    return statistics.median(xs) if xs else NAN


def nanmean(xs):
    xs = [x for x in xs if x == x]
    return sum(xs) / len(xs) if xs else NAN


def seg_has(arr, i0, i1):
    for i in range(max(0, i0), min(len(arr), i1)):
        if arr[i] == arr[i]:
            return True
    return False


def window_ok(series, req, i):
    for k in req:
        a = series[k]
        if not (seg_has(a, i - PRE_H, i) and seg_has(a, i, i + DS_WIN_H) and seg_has(a, i + DS_WIN_H, i + POST_H)):
            return False
    return True


def delta_s(a, i):
    pre = nanmedian(a[max(0, i - PRE_H):i])
    post = nanmean(a[i:i + DS_WIN_H])
    return post - pre, pre


def find_events(series, n):
    """返回 (dry_run, onsets[(i, acc24, first_hour_mm, n_rain_h)], stats)。只用雨量定义（I5/I6）。"""
    rain = series["rain"]
    dry_run = array("i", bytes(4 * n))
    run = 0
    for i in range(n):
        r = rain[i]
        run = run + 1 if (r == r and r < RAIN_HOUR_MM) else 0
        dry_run[i] = run
    onsets = []
    st = {"onset_after_24h_dry": 0, "rain24_incomplete": 0, "below_10mm": 0}
    for i in range(DRY_BEFORE_H, n - POST_H):
        r = rain[i]
        if not (r == r and r >= RAIN_HOUR_MM) or dry_run[i - 1] < DRY_BEFORE_H:
            continue
        st["onset_after_24h_dry"] += 1
        post = rain[i:i + POST_H]
        if any(x != x for x in post):
            st["rain24_incomplete"] += 1
            continue
        wet = [x for x in post if x >= RAIN_HOUR_MM]
        acc = sum(wet)
        if acc < RAIN_EVENT_MM:
            st["below_10mm"] += 1
            continue
        onsets.append((i, acc, r, len(wet)))
    return dry_run, onsets, st


def build_events(series, n, dry_run, onsets, req, h0, ctrl_ok=None):
    """在雨事件上施加数据要求（I7）并配对照（I9）。ctrl_ok：可选的对照起点附加筛选（P1 不传，行为不变）。"""
    blocked = sorted(i for i, *_ in onsets)  # 所有 ≥10 mm 雨事件起点
    s1 = series["s1"]

    def overlaps_event(c):
        # 对照窗 [c-6, c+24) 与任一事件窗 [e-6, e+24) 相交 ⇔ |c-e| < 30
        span = PRE_H + POST_H
        j = bisect_left(blocked, c - span + 1)
        return j < len(blocked) and blocked[j] < c + span

    events, n_drop = [], 0
    for i, acc, r0, nwet in onsets:
        if not window_ok(series, req, i):
            n_drop += 1
            continue
        ctrls = []
        for k in sorted([k for k in range(-CTRL_MAX_DAYS, CTRL_MAX_DAYS + 1) if k], key=lambda k: (abs(k), k)):
            c = i + 24 * k
            if c - DRY_BEFORE_H < 0 or c + POST_H > n:
                continue
            if dry_run[c - 1] < DRY_BEFORE_H or overlaps_event(c) or not window_ok(series, req, c):
                continue
            if ctrl_ok is not None and not ctrl_ok(c):
                continue
            ctrls.append(c)
            if len(ctrls) == CTRL_PER_EVENT:
                break
        ds1, pre1 = delta_s(s1, i)
        ev = {"i": i, "hour": h0 + i, "acc24": acc, "first_mm": r0, "n_wet": nwet,
              "ds1": ds1, "pre_s1": pre1, "ctrl": ctrls,
              "ctrl_ds1": [delta_s(s1, c)[0] for c in ctrls]}
        if "sss05" in series and "sss05" in req:
            ds05, pre05 = delta_s(series["sss05"], i)
            ev["ds05"], ev["pre_s05"] = ds05, pre05
        events.append(ev)
    return events, n_drop


# ======================================================================== 统计
def describe(xs):
    xs = [x for x in xs if x == x]
    if not xs:
        return {"n": 0}
    d = {"n": len(xs), "mean": sum(xs) / len(xs), "median": statistics.median(xs)}
    if len(xs) >= 2:
        d["sd"] = statistics.stdev(xs)
        q = statistics.quantiles(xs, n=20, method="inclusive")
        d.update({"p05": q[0], "p25": q[4], "p75": q[14], "p95": q[18]})
    return {k: round(v, 5) if isinstance(v, float) else v for k, v in d.items()}


def d2_test(evs, label):
    """事件均值 − 对照均值；站×季整簇 bootstrap（I10）。evs: [(station, season, ds1, [ctrl_ds1])]"""
    evs = [e for e in evs if e[2] == e[2] and any(c == c for c in e[3])]
    ev_vals = [e[2] for e in evs]
    ct_vals = [c for e in evs for c in e[3] if c == c]
    res = {"label": label, "n_events": len(ev_vals), "n_controls": len(ct_vals)}
    if len(ev_vals) < 2 or len(ct_vals) < 2:
        res.update({"evaluable": False, "pass": False, "reason": "事件或对照不足（I16）"})
        return res
    diff = sum(ev_vals) / len(ev_vals) - sum(ct_vals) / len(ct_vals)
    clusters = {}
    for stn, sea, ds, cds in evs:
        c = clusters.setdefault((stn, sea), [0.0, 0, 0.0, 0])
        c[0] += ds
        c[1] += 1
        for x in cds:
            if x == x:
                c[2] += x
                c[3] += 1
    keys = list(clusters)
    res["n_clusters"] = len(keys)
    w_se = math.sqrt(statistics.variance(ev_vals) / len(ev_vals) + statistics.variance(ct_vals) / len(ct_vals))
    res["welch_ci95_info"] = [round(diff - 1.96 * w_se, 5), round(diff + 1.96 * w_se, 5)]
    res["diff_psu"] = round(diff, 5)
    if len(keys) < 2:
        res.update({"evaluable": False, "pass": False, "reason": "簇数 <2，bootstrap 不可算（I16）"})
        return res
    rng = random.Random(BOOT_SEED)
    cl = [clusters[k] for k in keys]
    nk = len(cl)
    boots = []
    for _ in range(BOOT_B):
        se = ne = sc = nc = 0.0
        for _j in range(nk):
            c = cl[rng.randrange(nk)]
            se += c[0]; ne += c[1]; sc += c[2]; nc += c[3]
        if ne and nc:
            boots.append(se / ne - sc / nc)
    boots.sort()
    lo = boots[int(0.025 * (len(boots) - 1))]
    hi = boots[int(math.ceil(0.975 * (len(boots) - 1)))]
    ci_has_zero = lo <= 0.0 <= hi
    small = abs(diff) < D2_MIN_ABS_DIFF
    res.update({"evaluable": True, "ci95_cluster_bootstrap": [round(lo, 5), round(hi, 5)],
                "bootstrap_B": BOOT_B, "bootstrap_seed": BOOT_SEED,
                "abs_diff_below_0.03": small, "ci_contains_0": ci_has_zero,
                "pass": (not small) and (not ci_has_zero)})
    return res


def pearson(xs, ys):
    pts = [(x, y) for x, y in zip(xs, ys) if x == x and y == y]
    if len(pts) < 3:
        return None, len(pts)
    try:
        return statistics.correlation([p[0] for p in pts], [p[1] for p in pts]), len(pts)
    except statistics.StatisticsError:
        return None, len(pts)


def d5_metrics(evs):
    diffs = [e["pre_s05"] - e["pre_s1"] for e in evs
             if e.get("pre_s05", NAN) == e.get("pre_s05", NAN) and e["pre_s1"] == e["pre_s1"]]
    r, nr = pearson([e.get("ds05", NAN) for e in evs], [e["ds1"] for e in evs])
    med = statistics.median(diffs) if diffs else None
    return {"n_pre_pairs": len(diffs), "median_pre_diff_psu": _r(med),
            "corr_event_anom": _r(r), "n_corr": nr}


def _r(x, nd=5):
    return None if x is None or x != x else round(x, nd)


# ======================================================================== 主流程
def process_station(http, st, cache, log, smoke, include_rt=False):
    t_start = time.monotonic()
    hs, he, mrows = fetch_mapco2(http, st, cache, log, SMOKE_MONTH if smoke else None)
    h0 = hs - PAD_H
    n = he + PAD_H - h0 + 1
    ser = {k: Hourly(h0, n) for k in ("rain", "wind", "s1", "s5", "sss05", "pco2")}
    for h, sss, pc in mrows:
        if sss is not None:
            ser["sss05"].add(h, sss)
        if pc is not None:
            ser["pco2"].add(h, pc)
    del mrows
    files = list_gtmba_files(http, st, cache, hs, he, smoke, include_rt)
    log.log(f"[{st['name']}] GTMBA 文件 {len(files)} 个")
    qc_drop = {}
    for base, kind, per_hour, _sh in files:
        got = load_gtmba_file(http, cache, base, kind, per_hour, log, qc_drop)
        for key, pairs in got.items():
            s = ser[key]
            for h, v in pairs:
                s.add(h, v)
        del got
    series = {k: v.values() for k, v in ser.items()}
    cov = {f"valid_hours_{k}": ser[k].n_valid() for k in ser}
    del ser
    rain, s1, wind, sss05, pco2 = (series[k] for k in ("rain", "s1", "wind", "sss05", "pco2"))
    ov = sum(1 for i in range(n) if sss05[i] == sss05[i] and pco2[i] == pco2[i]
             and rain[i] == rain[i] and s1[i] == s1[i] and wind[i] == wind[i])
    dd = [sss05[i] - s1[i] for i in range(n) if sss05[i] == sss05[i] and s1[i] == s1[i]]
    cov["sss05_minus_s1_same_hour_info"] = {
        "n": len(dd), "median": _r(statistics.median(dd)) if dd else None,
        "median_abs": _r(statistics.median(abs(x) for x in dd)) if dd else None,
        "frac_abs_le_0.003": _r(sum(1 for x in dd if abs(x) <= 0.003) / len(dd), 4) if dd else None}
    del dd
    cov.update({"axis_hours": n, "mapco2_start": hour_to_iso(hs), "mapco2_end": hour_to_iso(he),
                "overlap_mapco2_samples_all5": ov,
                "overlap_station_years_approx": round(ov * 3 / 8766, 3),
                "qc_dropped_by_code": {k: {str(q): c for q, c in v.items()} for k, v in qc_drop.items()},
                "n_gtmba_files": len(files)})
    dry_run, onsets, ostats = find_events(series, n)
    req_main = ["sss05", "pco2", "s1", "wind"]
    req_alt = ["s1", "wind"]
    ev_main, drop_main = build_events(series, n, dry_run, onsets, req_main, h0)
    ev_alt, drop_alt = build_events(series, n, dry_run, onsets, req_alt, h0)
    s5 = series["s5"]
    for e in ev_main:
        i = e["i"]
        e["has_5m"] = seg_has(s5, i - PRE_H, i) and seg_has(s5, i, i + DS_WIN_H)
        e["n_pco2"] = sum(1 for x in pco2[max(0, i - PRE_H):i + POST_H] if x == x)
        e["ctrl_hours"] = [h0 + c for c in e["ctrl"]]
        e["season"] = season_of(e["hour"])
    for e in ev_alt:
        e["season"] = season_of(e["hour"])
    log.log(f"[{st['name']}] 雨起点(前24h干) {ostats['onset_after_24h_dry']}，24h 雨量不全 {ostats['rain24_incomplete']}，"
            f"<10mm {ostats['below_10mm']}，≥10mm 雨事件 {len(onsets)}；主事件 {len(ev_main)}（数据不足剔 {drop_main}）；"
            f"无0.5m分支 {len(ev_alt)}；耗时 {time.monotonic() - t_start:.0f}s")
    del series
    return {"coverage": cov, "onset_stats": ostats, "n_rain_events_10mm": len(onsets),
            "events": ev_main, "n_dropped_data": drop_main,
            "events_alt": ev_alt, "n_dropped_data_alt": drop_alt}


def summarize(results, out_dir, smoke, http, t0, log, include_rt=False):
    rows, per_station, regimes = [], {}, {}
    all_main, all_alt = [], []
    for st in STATIONS:
        if st["name"] not in results:
            continue
        r = results[st["name"]]
        evs = r["events"]
        full = [e for e in evs if len(e["ctrl"]) == CTRL_PER_EVENT]
        regimes[st["regime"]] = regimes.get(st["regime"], 0) + len(evs)
        per_station[st["name"]] = {
            "regime": st["regime"], "gtmba_site": st["site"], "mapco2_dataset": st["erddap"],
            "coverage": r["coverage"], "rain_onset_stats": r["onset_stats"],
            "n_rain_events_10mm": r["n_rain_events_10mm"], "n_events": len(evs),
            "n_dropped_insufficient_data": r["n_dropped_data"],
            "n_events_5ctrl": len(full),
            "n_events_ctrl_lt5": len(evs) - len(full),
            "n_events_has_5m": sum(1 for e in evs if e["has_5m"]),
            "dS1_events": describe([e["ds1"] for e in evs]),
            "dS1_controls": describe([c for e in evs for c in e["ctrl_ds1"]]),
            "d5_info": d5_metrics(evs),
            "n_events_no05m_branch": len(r["events_alt"]),
        }
        for e in evs:
            all_main.append((st, e))
            rows.append({
                "station": st["name"], "regime": st["regime"], "onset_utc": hour_to_iso(e["hour"]),
                "local_solar_hour": round((e["hour"] % 24 + st["lon"] / 15.0) % 24, 2),
                "season": e["season"], "rain_first_hour_mm": round(e["first_mm"], 3),
                "rain_24h_mm": round(e["acc24"], 3), "n_rain_hours_24h": e["n_wet"],
                "pre_S1_median": _r(e["pre_s1"], 4), "dS1_0_6h": _r(e["ds1"], 4),
                "pre_S05_median": _r(e.get("pre_s05"), 4), "dS05_0_6h": _r(e.get("ds05"), 4),
                "has_5m": int(e["has_5m"]), "n_pco2_in_window": e["n_pco2"], "n_ctrl": len(e["ctrl"]),
                "ctrl_dS1_mean": _r(nanmean(e["ctrl_ds1"]), 4),
                "ctrl_onsets_utc": ";".join(hour_to_iso(h) for h in e["ctrl_hours"]),
                "ctrl_dS1": ";".join("" if x != x else f"{x:.4f}" for x in e["ctrl_ds1"]),
            })
        for e in r["events_alt"]:
            all_alt.append((st, e))

    n_total = len(all_main)
    n_st10 = sum(1 for v in per_station.values() if v["n_events"] >= D1_MIN_EVENTS_PER_STATION)
    d1 = {"n_events_total": n_total, "n_stations_ge10": n_st10,
          "threshold": f"总数 ≥{D1_MIN_TOTAL} 且 ≥{D1_MIN_EVENTS_PER_STATION} 事件的站 ≥{D1_MIN_STATIONS}",
          "events_by_station": {k: v["n_events"] for k, v in per_station.items()},
          "events_by_regime": regimes,
          "pass": n_total >= D1_MIN_TOTAL and n_st10 >= D1_MIN_STATIONS}

    def d2_input(pairs, need_full):
        out = []
        for st, e in pairs:
            if need_full and len(e["ctrl"]) != CTRL_PER_EVENT:
                continue
            if not need_full and not e["ctrl"]:
                continue
            out.append((st["name"], e["season"], e["ds1"], e["ctrl_ds1"]))
        return out

    d2 = d2_test(d2_input(all_main, True), "主判：配满 5 个对照的事件")
    d2["threshold"] = f"|事件均值−对照均值| ≥{D2_MIN_ABS_DIFF} psu 且 95% CI 不含 0"
    d2["sensitivity_ge1ctrl_info"] = d2_test(d2_input(all_main, False), "敏感性：≥1 个对照（不参与判定）")

    evs = [e for _, e in all_main]
    d5m = d5_metrics(evs)
    if d5m["median_pre_diff_psu"] is None or d5m["corr_event_anom"] is None:
        d5 = dict(d5m, evaluable=False, pass_=False, reason="配对不足（I16）")
    else:
        fail = abs(d5m["median_pre_diff_psu"]) > D5_MAX_ABS_MEDDIFF or d5m["corr_event_anom"] < D5_MIN_CORR
        d5 = dict(d5m, evaluable=True, pass_=not fail)
    d5["pass"] = d5.pop("pass_")
    d5["threshold"] = f"|雨前中位差| ≤{D5_MAX_ABS_MEDDIFF} psu 且 事件异常相关 ≥{D5_MIN_CORR}"
    d5["on_fail_note"] = "设计冲突：一条规则弃用 0.5 m 只用 1/5 m，另一条（出口 3）停止；脚本不作裁定（I11）"

    n_d4 = sum(1 for e in evs if e.get("ds05", NAN) == e.get("ds05", NAN) and e["ds05"] <= D4_DS_THRESH)
    d4 = {"n_events_dS05_le_-0.2": n_d4, "threshold": f"≥{D4_MIN_EVENTS}",
          "would_pass": n_d4 >= D4_MIN_EVENTS, "note": "信息量；未列入 P1 判定（I15）"}

    alt_n = len(all_alt)
    alt_st10 = sum(1 for st in STATIONS if st["name"] in results
                   and len(results[st["name"]]["events_alt"]) >= D1_MIN_EVENTS_PER_STATION)
    alt = {"note": "D5 不过分支：事件只要求 1 m S＋风；仅信息量（I11）",
           "D1": {"n_events_total": alt_n, "n_stations_ge10": alt_st10,
                  "would_pass": alt_n >= D1_MIN_TOTAL and alt_st10 >= D1_MIN_STATIONS},
           "D2": d2_test(d2_input(all_alt, True), "无 0.5 m 分支：配满 5 个对照")}

    summary = {
        "script": "p1_events.py", "version": VERSION, "smoke": smoke, "include_realtime_nondefault": include_rt,
        "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "runtime_s": round(time.monotonic() - t0, 1),
        "params_fixed": {
            "rain_event_mm": RAIN_EVENT_MM, "dry_before_h": DRY_BEFORE_H, "window_h": [-PRE_H, POST_H],
            "dS_window_h": [0, DS_WIN_H], "dS_baseline": "pre-6h median", "controls_per_event": CTRL_PER_EVENT,
            "control_season_days": CTRL_MAX_DAYS, "gtmba_qc_ok": sorted(QC_OK), "mapco2_flag": 2},
        "params_implementation": {
            "rain_hour_threshold_mm_h": RAIN_HOUR_MM, "subhourly_min_frac": SUBHOURLY_MIN_FRAC,
            "bootstrap_B": BOOT_B, "bootstrap_seed": BOOT_SEED, "bootstrap_cluster": "station×season(DJF/MAM/JJA/SON)",
            "deploy_lookback_days": DEPLOY_LOOKBACK_DAYS,
            "gtmba_file_modes": ["D", "M", "R"] if include_rt else ["D", "M"],
            "see_docstring": "I1–I17"},
        "gates": {"D1": d1, "D2": d2, "D5": d5, "D4_info": d4,
                  "D3": "P2（需 CMORPH＋RIM-3/S20），P1 不判"},
        "gates_verdict": {"D1": d1["pass"], "D2": d2["pass"], "D5": d5["pass"]},
        "branch_no05m_info": alt,
        "stations": per_station,
        "downloads": {"http_requests": http.n_requests, "cache_hits": http.n_cached,
                      "bytes_downloaded_gz": http.bytes},
    }
    with open(os.path.join(out_dir, "p1_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2, allow_nan=False)
    fields = ["station", "regime", "onset_utc", "local_solar_hour", "season", "rain_first_hour_mm",
              "rain_24h_mm", "n_rain_hours_24h", "pre_S1_median", "dS1_0_6h", "pre_S05_median",
              "dS05_0_6h", "has_5m", "n_pco2_in_window", "n_ctrl", "ctrl_dS1_mean",
              "ctrl_onsets_utc", "ctrl_dS1"]
    with open(os.path.join(out_dir, "p1_events.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    log.log(json.dumps(summary["gates_verdict"], ensure_ascii=False))
    return summary


def main():
    ap = argparse.ArgumentParser(description="P1 系泊雨事件")
    ap.add_argument("--out", help="输出目录（默认 $REPRO_OUTPUT_DIR）")
    ap.add_argument("--smoke", action="store_true", help=f"只跑 {SMOKE_STATION} {SMOKE_MONTH}")
    ap.add_argument("--stations", help="逗号分隔站名子集（调试用）")
    ap.add_argument("--include-realtime", action="store_true",
                    help="加入 GTMBA 实时(R)文件（非默认，仅敏感性）")
    args = ap.parse_args()
    out_dir = args.out or os.environ.get("REPRO_OUTPUT_DIR")
    if not out_dir:
        print("需要 --out 或环境变量 REPRO_OUTPUT_DIR", file=sys.stderr)
        return 3
    os.makedirs(out_dir, exist_ok=True)
    cache = os.path.join(out_dir, "cache")
    os.makedirs(cache, exist_ok=True)
    log = Log(os.path.join(out_dir, "p1_log.txt"))
    t0 = time.monotonic()
    names = [SMOKE_STATION] if args.smoke else (
        args.stations.split(",") if args.stations else [s["name"] for s in STATIONS])
    unknown = set(names) - {s["name"] for s in STATIONS}
    if unknown:
        print(f"未知站名：{sorted(unknown)}", file=sys.stderr)
        return 3
    log.log(f"=== start {VERSION} smoke={args.smoke} realtime={args.include_realtime} stations={names} out={out_dir}",
            echo=True)
    http = Http(log)
    try:
        results = {}
        for st in STATIONS:
            if st["name"] not in names:
                continue
            print(f"[{st['name']}] 下载与处理中…", flush=True)
            r = process_station(http, st, cache, log, args.smoke, args.include_realtime)
            results[st["name"]] = r
            sd = r["coverage"]["sss05_minus_s1_same_hour_info"]
            print(f"[{st['name']}] 雨事件≥10mm {r['n_rain_events_10mm']}，合格事件 {len(r['events'])}，"
                  f"配满5对照 {sum(1 for e in r['events'] if len(e['ctrl']) == CTRL_PER_EVENT)}；"
                  f"SSS0.5−S1 |差|≤0.003 占比 {sd['frac_abs_le_0.003']}", flush=True)
        s = summarize(results, out_dir, args.smoke, http, t0, log, args.include_realtime)
    except FetchError as e:
        log.log(f"FATAL 数据源故障：{e}", echo=True)
        return 2
    except Exception:
        log.log("FATAL 未预期异常：\n" + traceback.format_exc(), echo=True)
        return 3
    g = s["gates"]
    print("---- P1 门判定" + ("（SMOKE，仅验证格式，结论无效）" if args.smoke else "") + " ----")
    print(f"D1 N={g['D1']['n_events_total']}，≥10事件站={g['D1']['n_stations_ge10']} → {'过' if g['D1']['pass'] else '不过'}")
    d2 = g["D2"]
    print(f"D2 diff={d2.get('diff_psu')} psu，CI={d2.get('ci95_cluster_bootstrap')} → {'过' if d2['pass'] else '不过'}"
          + ("" if d2.get("evaluable", False) else f"（{d2.get('reason')}）"))
    d5 = g["D5"]
    print(f"D5 中位差={d5.get('median_pre_diff_psu')}，相关={d5.get('corr_event_anom')} → {'过' if d5['pass'] else '不过'}")
    print(f"D4(信息量) ΔS0.5≤−0.2 事件数={g['D4_info']['n_events_dS05_le_-0.2']}")
    print(f"下载 {http.n_requests} 次请求 / {http.bytes / 1e6:.1f} MB(gz)，缓存命中 {http.n_cached}；"
          f"耗时 {time.monotonic() - t0:.0f}s；产物在 {out_dir}")
    log.log("=== done")
    log.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
