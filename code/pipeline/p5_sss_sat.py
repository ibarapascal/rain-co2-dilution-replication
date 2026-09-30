#!/usr/bin/env python3
"""p5_sss_sat.py — P5：SMAP 卫星海表盐度共址检验。用雨后过境的卫星海表淡化与 RIM-3 足印尺度海面预测之比 R_sfc，
区分 R1（RIM-3 幅度整体偏大，海面也偏大）与 R2（海面幅度对、剖面画得太深）；系泊 0.5／1 m 作第三方。

方法与判读规则在运行前写定；判读在 judge()。性质：探索性，对卫星 SSS 盲。

用法（产物写 $REPRO_OUTPUT_DIR，未设则 --out）：
  python p5_sss_sat.py --selftest      合成数据自测（无网络；给 --out 时写 p5_selftest.json）
  python p5_sss_sat.py --plan          自测 → 重建事件 → CMR 过境计数 → 认证探测（.dmr 维序＋站点格 cellat/cellon/time，不取盐度）
                                       → 下载量估计 → p5_plan.json
  python p5_sss_sat.py --full          自测 → 下载（SMAP 子集 → 精确过境 → 所需 CMORPH 小时的站点窗口）→ 分析
  python p5_sss_sat.py --fetch         只下载（断点续传）；--analyze 只读缓存分析（与 --full 的分析段逐值相同）
  --cache DIR       缓存根（默认 <fast_root>/p5-cache）
  --p1-events / --p1b-dir               同 P2（默认 workstation 上 P1／P1b 产物；系泊风与盐度从其缓存重建）
  --token-file P    Earthdata token 文件路径（默认 workstation 绝对路径 <token-file>；只读该路径，内容只进请求头）
  --max-conc N      并发（≤4）
  --raw-root DIR    原始下载文件归档根（a6 起；默认 raw_store.RAW_ROOT_DEFAULT＝workstation <raw_root>）
  --rate-mbps R     CMORPH 下载总限速（MB/s，1 MB＝10^6 B；a6 起，默认 2）
  --raw-tag S       写进 raw/manifest.jsonl 的来源任务标记（a6 起）
依赖：numpy、scipy、netCDF4＋同目录 p1_events.py、p1b_extend.py、p2_rim_test.py、p4_mech.py（只 import，不改）；
  a6 起另 import 同目录 raw_store.py（只做原件归档 I/O）。
数据：CMR（公开）；PO.DAAC Cloud OPeNDAP 的 RSS SMAP L2C V6（Earthdata token）；NCEI CMORPH 8 km／30 min（匿名）；P1/P1b 缓存。
产物（<out>/）：p5_selftest.json、p5_plan.json（plan）、p5_fetch.json（fetch）、p5_power.json（分析段先写）、p5_summary.json、
  p5_events.csv、p5_log.txt。
缓存（--cache）：cmr/<站>.json、smap/<站>/<granule>.npz（3×3×2 子集）或 .missing、cmorph/<站>/grid.json 与 <hour>.npy
  （float32 (2,ny,nx)）、cmorph/missing.jsonl、fetch_manifest.json。
退出码：0 跑完（无论判读）；2 数据源故障（重跑即续传）；3 其他异常（含维序核对不过）；4 自测不过；5 认证失败。

实现选择（V 条）：
  V1 事件：p2.build_all 重建 646 个 10 mm 事件（K1 逐条校验），取 onset ≥ 2015-04-01T00Z；t0＝H·3600 s；簇＝站×季（P1 I10）。
  V2 轨道：CMR granules.csv 点查询（公开、无 token），按站覆盖其 SMAP 期事件时段；估计过境时刻＝地方时 06 或 18 时落在该轨时段
     （±10 min）者，否则取轨道中点；只用于挑候选轨与 plan 计数，判窗用文件 time。
  V3 子集：OPeNDAP DAP4 `<granule>.dap.nc4?dap4.ce=/v[y0:1:y1][x0:1:x1][0:1:1];…`，3×3 格，y＝⌊(lat+90)/0.25⌋（plan 实测正确）；
     **x 是逐轨变化的轨道网格**（plan a4 实测：同一行 lon 随 x 每格 −0.25°、起点逐轨不同；RSS 文档要求用 cellon 定位），所以每轨先取
     第 y 行的 cellat/cellon（1×1560×2），取离站点格心（纬 −89.875+0.25y，经 0.125+0.25⌊lon/0.25⌋）最近且 ≤0.2° 的有值格为 x；
     无则记「未覆盖」（.nocover）；x∉[1,1558] 同。3×3 子集取回后再核中心格坐标 ≤0.2°，不符即退出 3。站点落在格界时取 floor 一侧。
     维序须为 (ydim_grid, xdim_grid, look)，不符即退出 3（plan 用 .dmr 核实）。
  V4 过境：视向时刻＝中心格 time＋946684800 s；过境时刻＝有效视向平均；「覆盖」＝任一视向中心格 sss_smap_40km 非缺
     （下载阶段只凭覆盖决定要哪些 CMORPH 小时，不用 QC、不读盐度数值）。
  V5 QC（qc_valid）：无效＝iqc 位 {0,1,2,3,4,16} 任一、SSS∉[2,42]、winspd>20；剔除位 {5,6,7,8,9,10,14}；保留 11、12、13、15（雨）。
  V6 卫星值：A70＝中心格有效且 3×3 有效 ≥5 格时有效格均值；A40＝中心格；过境值＝有效视向均值。
  V7 窗：雨后取 t_s−t0∈[0,12) h 内最早的有效过境；参照取 [−72,0) h 内全部有效过境（≥1）；有效＝卫星值与模型值都非 NaN。
  V8 足印核：每格 c 以该视向该格的 cellat/cellon（OI 格实际位置）为中心的 40 km FWHM 圆高斯×cos φ（像元面积），截断 60 km，归一；
     A70 核＝该视向有效格核的均值；A40＝中心格核。坐标缺的格视为无效。距离用等距近似（R＝6371 km）。模型像元超集＝站点格心
     100 km 内的 CMORPH 像元。
  V9 RIM：vec_rim0() 与 p2.rim_factor(z=0) 同式向量化（自测逐值相对差 ≤1e−12）；q＝⌊t_s/1800⌋，用半步 q−48..q；像元 F 为 NaN
     （CMORPH 缺值、P>200 mm/h）时丢像元重归一，丢弃权重 >5% 该视向模型无效；所需小时文件缺则无效。
  V10 风：站点小时风按 P2 K9（p2.fill_wind ≤6 h、×p2.WIND_FACTOR、下限 0.1），半步 q 用小时 ⌊q/2⌋；区间内有缺则该视向无效。
  V11 ΔS_RIM,foot＝S0·[F̄(雨后) − mean F̄(参照)]；ΔSSS_sat＝S(雨后) − mean S(参照)；S0＝雨后过境中心格 sss_ref（缺则参照卫星均值）。
  V12 系泊第三方：站点像元＝CMORPH 窗口中离站点最近者；p2.rim_factor 在半步 [2H−60, q] 上算 z＝s1_depth（K28）与 0.5 m，基线＝
     雨前 6 h 小时 F（两半步平均）的中位数，S0_z＝p2.pre_median；观测 s1 取 t_s 所在小时 ±1 h、sss05 ±2 h 有效均值 − p2.pre_median。
     站点像元 z=0 另按足印同一差分结构（雨后减参照过境）算 x_pix0，给尺度因子。
  V13 推断：p2.cluster_boot（站×季、站两种簇，B=10000，seed 20260926，百分位 95%）；pm.jackknife_station（逐站剔除）。
  V14 可评：n≥30、站×季簇≥8、bootstrap 中 Σx≥0 的比例 ≤1%。
  V15 判读：judge()（含「稳健」标注）。
  V16 功效：p5_power.json 在算任何 ΣΔSSS_sat 之前写盘（只含模型 x、形式不确定度、n、簇）；视向 A70 形式 σ＝有效格 unc 均值/√n_valid
     （独立近似，偏乐观），过境 σ＝√(mean σ_look²/n_look)，事件 σ＝√(σ_post²+mean σ_ref²/k)。
  V17 敏感性 S1–S9 各自完整重算事件集，只报 n、簇、点估计、CI、可评、类别与是否与主判同类。
  V18 S6 同交点：过境地方时（UTC＋lon/15）<12 h 为上午、否则下午；参照只留与雨后过境同侧者。
  V19 S7：足印平均雨量 P̄(q)＝Σw_p P_p(q)/Σw_p（像元 NaN 同 V9 丢弃规则）驱动单列 p2.rim_factor(z=0)。
  V20 CMORPH：p2.Fetcher 整小时下载（≤4 并发、5 次重试），核 time 轴，抽各站窗口（纬向 ±17 像元、经向 ±⌈17/cos φ⌉），存 float32
      (2,ny,nx)，原文件即删；404 两次或时间轴不符记缺测（missing.jsonl）。
  V21 所需小时：对「有覆盖的雨后候选（[0,24) h）」与「有覆盖的参照候选（[−72,0) h）」都有的事件，并入每个雨后候选的 [H−30, ⌊q/2⌋]
      与每个参照候选的 [⌊(q−48)/2⌋, ⌊q/2⌋]（q 取两视向的最小／最大）。
  V22 HTTP：只对 opendap.earthdata.nasa.gov、archive.podaac.earthdata.nasa.gov、*.earthdatacloud.nasa.gov 附 Bearer；不自动跟随重定向，
      逐跳判断，到 urs.earthdata.nasa.gov 即认证失败（退出 5）；401／403 同；异常文本经 redact()；token 不写任何文件。
  V23 同一过境可同时是一个事件的雨后过境与另一事件的参照，分别计入（同 P2 对照共用）；站×季整簇重抽已含相关。
  V24 描述项：过原点斜率、逐站点估计、尺度因子、IMERG 中心格雨量与 CMORPH 足印雨量、Q＝R_sfc／R_1m,matched（同一事件子集同一抽样）。
  V25 原件保留（a6，只改 I/O；取代 V20 中「原文件即删」一句，V20 其余不变）：CMORPH 小时原件经 raw_store.CmorphRaw 取——
      raw 档（--raw-root）已有同名原件就直接抽窗口、不重下；没有就按 --rate-mbps 令牌桶限速下载进 raw 档
      （cmorph/<YYYY>/<MM>/<DD>/原文件名，与 NCEI 层级相同），逐文件记 raw/manifest.jsonl（url、path、bytes、sha256、t、status、
      source_job）；抽窗口、核 time 轴、缺测记录与 V20 完全相同；**任何原件都不删**。下载段每 ≤20 s 刷新
      raw/.p5_fetch_heartbeat，下载段结束写 done（补档任务 raw_archive.py 据此把自身限速在 2↔4 MB/s 间切换）。
      新取的 SMAP OPeNDAP 子集（行坐标与 3×3）原样另存 raw/smap_opendap_subset/<站>/（a6 提交时 790 轨×站已全在缓存，预计不触发）。

Change Log：
  2026-09-27 初版。
  2026-09-27 a2：token 默认路径改绝对路径（任务 HOME 与交互环境不同时 plan 首跑因此退出 5）；无数值路径改动。
  2026-09-27 a3：plan 首次成功运行显示 9 站中心格 cellat/cellon/time 全缺——plan 探测改为报返回形状、CE 编码与否对比、
     17×17 框内有限格定位（只取坐标与时刻，不取盐度），用于定位 CE／格号问题。
  2026-09-27 a5：plan a4 实测 x 为逐轨网格——改为逐轨按行定位 x（locate_in_row）、核中心用 cellat/cellon、
     模型像元超集改为站点 100 km 内；plan 探测改为「定位＋3×3 坐标核对」；自测加定位与实测坐标核。
  2026-09-27 a4：a3 的 plan 卡在不编码 CE 的请求上撞硬上限——探测改为编码 CE 的 3×3＋站点整行／整列坐标（90 s、2 次），
     逐站落盘 p5_plan_probe_partial.json。
  2026-09-27 a6：原始数据一律保留——此前的正式运行因逐小时删 CMORPH 原件被取消（未产出
     任何分析文件）；本版只改 I/O（V25）：原件落 bulk raw 档、先查档后下载、限速参数、心跳；分析逻辑、判读、参数一行未改。
"""

import argparse
import concurrent.futures as cf
import csv
import hashlib
import http.client
import json
import math
import os
import shutil
import socket
import sys
import tempfile
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import p1_events as p1
import p2_rim_test as p2
import p4_mech as pm
import raw_store as rs

VERSION = "p5-2026-09-27a6"
NAN = float("nan")
EPOCH2000 = 946684800

# ---- 事先写定的参数 ----
ERA_START_H = int(datetime(2015, 4, 1, tzinfo=timezone.utc).timestamp() // 3600)
THRESH = 0.65
MIN_N, MIN_G = 30, 8
DEN_NONNEG_MAX = 0.01
FWHM_KM, TRUNC_KM = 40.0, 60.0
A70_MIN_VALID = 5
QC_INVALID = (0, 1, 2, 3, 4, 16)
QC_EXCLUDE = (5, 6, 7, 8, 9, 10, 14)
SSS_RANGE = (2.0, 42.0)
WIND_MAX = 20.0
DROP_MAX = 0.05
OBS_HALF_H = {"s1": 1, "sss05": 2}
PRIMARY = {"name": "primary", "agg": "A70", "post": (0.0, 12.0), "ref": (-72.0, 0.0), "qc_ex": QC_EXCLUDE,
           "same_node": False, "drive": "pixel", "excl": ()}
SENS = [
    ("S1", "A40（只中心格）", {"agg": "A40"}),
    ("S2", "雨后窗 [0,6) h", {"post": (0.0, 6.0)}),
    ("S3", "雨后窗 [0,24) h", {"post": (0.0, 24.0)}),
    ("S4", "QC 保留位 10", {"qc_ex": (5, 6, 7, 8, 9, 14)}),
    ("S5", "另剔除位 12（风>15）", {"qc_ex": QC_EXCLUDE + (12,)}),
    ("S6", "参照只用同交点", {"same_node": True}),
    ("S7", "足印平均雨量驱动 RIM", {"drive": "rainmean"}),
    ("S8", "剔除 KEO、Papa", {"excl": ("KEO", "Papa")}),
    ("S9", "参照窗 [-48,0) h", {"ref": (-48.0, 0.0)}),
]
B = p2.BOOT_B
SEED = p2.BOOT_SEED

# ---- 数据访问 ----
COLL = "C2832221740-POCLOUD"                 # SMAP_RSS_L2_SSS_V6
AQ_COLL = "C2036882456-POCLOUD"              # AQUARIUS_L2_SSS_V5（只计数）
AQ_SPAN = ("2011-08-25T00:00:00Z", "2015-06-08T00:00:00Z")
CMR_URL = "https://cmr.earthdata.nasa.gov/search/granules.csv"
OPENDAP = f"https://opendap.earthdata.nasa.gov/collections/{COLL}/granules/"
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
TOKEN_FILE_DEFAULT = _rp.path("token_file")  # [repro] 只取路径字符串，内容仍只进请求头；位置由本机配置
CACHE_DEFAULT = _rp.path("p5_cache")  # [repro] 路径来自集中配置
VARS3 = ("sss_smap_40km", "sss_smap_40km_unc", "time", "iqc_flag", "cellat", "cellon")   # [y][x][look]
VARS2 = ("rain", "sss_ref", "winspd", "surtep")                                           # [y][x]
DIMS3 = ("ydim_grid", "xdim_grid", "look")
UA = "rain-co2-dilution-replication-p5/1.0 (research; python-urllib)"
AUTH_HOSTS = ("opendap.earthdata.nasa.gov", "archive.podaac.earthdata.nasa.gov")
AUTH_SUFFIX = ".earthdatacloud.nasa.gov"
CAND_H = (-74.0, 25.0)                       # 候选轨（估计时刻）相对 t0 的范围
FETCH_POST_H = (0.0, 24.0)                   # V21
FETCH_REF_H = (-72.0, 0.0)
WIN_LAT_PIX = 17
LOCATE_MAX_DEG = 0.2
SUPERSET_KM = 100.0
MAX_CONC = 4
MAX_ATTEMPTS = 5
PCA_FILES_PER_S = p2.PCA_FILES_PER_S
R_EARTH_KM = 6371.0
KM_PER_DEG = R_EARTH_KM * math.pi / 180.0


class AuthError(Exception):
    pass


class DataError(Exception):
    pass


class NotFound(Exception):
    pass


def _np():
    import numpy
    return numpy


def isnum(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def rnd(x, nd=5):
    return round(float(x), nd) if isnum(x) else None


def iso(ts):
    return datetime.fromtimestamp(ts, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_iso(s):
    return datetime.strptime(s[:19], "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc).timestamp()


def lon180(lon):
    return ((lon + 180.0) % 360.0) - 180.0


def safe(name):
    return name.replace("/", "_")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jdump(obj, path):
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1, default=lambda o: None)
    os.replace(tmp, path)


# ======================================================================== 认证与 HTTP（V22）
_SECRETS = []


def redact(s):
    s = str(s)
    for t in _SECRETS:
        if t:
            s = s.replace(t, "***")
    return s


def read_token(path):
    try:
        with open(path, "r", encoding="utf-8") as fh:
            raw = fh.read()
    except FileNotFoundError:
        raise AuthError(f"找不到 Earthdata token 文件 {path}")
    except OSError as e:
        raise AuthError(f"读不了 token 文件 {path}（{type(e).__name__}）")
    tok = next((ln.strip() for ln in raw.splitlines() if ln.strip()), "")
    if tok.lower().startswith("bearer "):
        tok = tok[7:].strip()
    if not tok or any(c.isspace() for c in tok):
        raise AuthError(f"token 文件 {path} 为空或格式不对（应只含一行 token）")
    _SECRETS.append(tok)
    return tok


def auth_host(host):
    h = (host or "").lower()
    return h in AUTH_HOSTS or h.endswith(AUTH_SUFFIX)


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


_OPENER = urllib.request.build_opener(_NoRedirect)


def build_headers(host, token):
    h = {"User-Agent": UA}
    if token and auth_host(host):
        h["Authorization"] = "Bearer " + token
    return h


def http_get(url, token=None, timeout=120, want_headers=False, extra=None):
    cur = url
    for _ in range(10):
        parts = urllib.parse.urlsplit(cur)
        if parts.scheme != "https":
            raise DataError(f"拒绝非 https 地址 {parts.scheme}://{parts.netloc}")
        host = parts.hostname or ""
        if host.lower() == "urs.earthdata.nasa.gov":
            raise AuthError("被重定向到 Earthdata 登录页：token 未被接受（过期／无效，或该数据应用未授权）")
        headers = build_headers(host, token)
        if extra:
            headers.update(extra)
        req = urllib.request.Request(cur, headers=headers)
        try:
            with _OPENER.open(req, timeout=timeout) as r:
                body = r.read()
                return (body, dict(r.headers)) if want_headers else body
        except urllib.error.HTTPError as e:
            code = e.code
            loc = e.headers.get("Location") if e.headers else None
            e.close()
            if code in (301, 302, 303, 307, 308) and loc:
                cur = urllib.parse.urljoin(cur, loc)
                continue
            if code in (401, 403):
                raise AuthError(f"HTTP {code} @ {host}（token 无效／过期或无权限）")
            if code == 404:
                raise NotFound(f"HTTP 404 @ {host}")
            raise DataError(f"HTTP {code} @ {host}")
        except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError, http.client.HTTPException) as e:
            raise DataError(redact(f"{type(e).__name__}: {getattr(e, 'reason', e)} @ {host}"))
    raise DataError("重定向超过 10 次")


def get_retry(url, token=None, timeout=120, attempts=MAX_ATTEMPTS, want_headers=False, extra=None):
    last = None
    for a in range(1, attempts + 1):
        try:
            return http_get(url, token, timeout, want_headers, extra)
        except (AuthError, NotFound):
            raise
        except DataError as e:
            last = e
        if a < attempts:
            time.sleep(min(40, 5 * 2 ** (a - 1)))
    raise DataError(redact(f"{attempts} 次失败：{last}"))


# ======================================================================== CMR（V2）
def cmr_granules(coll, lat, lon, t0, t1):
    """返回 [(granule_ur, start_unix, end_unix)]（点查询，公开，无 token）。"""
    out, sa = [], None
    while True:
        q = urllib.parse.urlencode({"collection_concept_id": coll, "point": f"{lon180(lon):.4f},{lat:.4f}",
                                    "temporal": f"{iso(t0)},{iso(t1)}", "page_size": 2000, "sort_key": "start_date"})
        body, hdr = get_retry(f"{CMR_URL}?{q}", None, 120, want_headers=True,
                              extra={"CMR-Search-After": sa} if sa else None)
        rows = list(csv.reader(body.decode("utf-8").splitlines()))[1:]
        for x in rows:
            if len(x) >= 4 and x[2] and x[3]:
                out.append((x[0], parse_iso(x[2]), parse_iso(x[3])))
        sa = hdr.get("CMR-Search-After") or hdr.get("Cmr-Search-After")
        if len(rows) < 2000 or not sa:
            break
    return sorted(set(out), key=lambda g: g[1])


def approx_pass(start, end, lon):
    d0 = math.floor(start / 86400.0) * 86400.0 - 86400.0
    for k in range(4):
        for lst in (6.0, 18.0):
            t = d0 + k * 86400.0 + (lst - lon180(lon) / 15.0) * 3600.0
            if start - 600 <= t <= end + 600:
                return t
    return 0.5 * (start + end)


# ======================================================================== SMAP 子集（V3–V6）
def cell_yx(lat, lon):
    y = int(math.floor((lat + 90.0) / 0.25))
    x = int(math.floor((lon % 360.0) / 0.25))
    return min(max(y, 1), 718), min(max(x, 1), 1558)


def cell_center(y, x):
    return -89.875 + 0.25 * y, 0.125 + 0.25 * x


def station_cell(lat, lon):
    """站点所在 0.25° 格的名义格心（纬、经 0–360）。"""
    y, _ = cell_yx(lat, lon)
    return -89.875 + 0.25 * y, 0.125 + 0.25 * math.floor((lon % 360.0) / 0.25)


def locate_in_row(la, lo, clat, clon):
    """la、lo：(1560, 2) 第 y 行的 cellat/cellon。返回 (x, 距离°, 该 x 两侧三格的有值数) 或 (None, 最近距离, 0)。
    多个候选（390° 网格首尾重叠）取三格有值数多者，再取距离近者。"""
    np = _np()
    ok = np.isfinite(la) & np.isfinite(lo)
    if not ok.any():
        return None, None, 0
    d = np.where(ok, np.sqrt((la - clat) ** 2 + (((lo - clon + 180.0) % 360.0) - 180.0) ** 2), np.inf)
    dmin = d.min(axis=1)
    cands = [int(i) for i in np.nonzero(dmin <= LOCATE_MAX_DEG)[0]]
    if not cands:
        return None, float(np.nanmin(dmin)), 0
    def score(i):
        n3 = int(ok[max(i - 1, 0):i + 2].sum())
        return (-n3, float(dmin[i]))
    best = min(cands, key=score)
    if not (1 <= best <= 1558):
        return None, float(dmin[best]), 0
    return best, float(dmin[best]), -score(best)[0]


def row_url(gid, y):
    ce = ";".join(f"/{v}[{y}:1:{y}][0:1:1559][0:1:1]" for v in ("cellat", "cellon"))
    return OPENDAP + urllib.parse.quote(gid) + ".dap.nc4?dap4.ce=" + urllib.parse.quote(ce, safe="")


def locate_x(gid, y, clat, clon, token, keep=None):
    buf = get_retry(row_url(gid, y), token, 120)
    if keep:
        keep(buf, f"row_y{y}", row_url(gid, y))
    sub = parse_nc(buf, ("cellat", "cellon"))
    la = sub["cellat"].reshape(-1, 2)
    lo = sub["cellon"].reshape(-1, 2)
    x, dist, n3 = locate_in_row(la, lo, clat, clon)
    return x, {"row_dist_deg": rnd(dist, 4) if dist is not None else None, "row_n3": n3,
               "row_n_finite": int((_np().isfinite(la)).sum())}


def check_center(sub, clat, clon):
    """3×3 取回后核中心格坐标（任一视向有值即核）。"""
    np = _np()
    for lk in (0, 1):
        a, b = sub["cellat"][1, 1, lk], sub["cellon"][1, 1, lk]
        if np.isfinite(a) and np.isfinite(b):
            dd = math.sqrt((a - clat) ** 2 + (((b - clon + 180.0) % 360.0) - 180.0) ** 2)
            if dd > LOCATE_MAX_DEG:
                raise ValueError(f"V3 中心格坐标偏离 {dd:.3f}°")
            return dd
    return None


def subset_url(gid, y, x):
    s3 = f"[{y - 1}:1:{y + 1}][{x - 1}:1:{x + 1}][0:1:1]"
    s2 = f"[{y - 1}:1:{y + 1}][{x - 1}:1:{x + 1}]"
    ce = ";".join([f"/{v}{s3}" for v in VARS3] + [f"/{v}{s2}" for v in VARS2])
    return OPENDAP + urllib.parse.quote(gid) + ".dap.nc4?dap4.ce=" + urllib.parse.quote(ce, safe="")


def coord_url(gid, y, x, half=1, enc=True):
    s3 = f"[{y - half}:1:{y + half}][{x - half}:1:{x + half}][0:1:1]"
    ce = ";".join(f"/{v}{s3}" for v in ("cellat", "cellon", "time"))
    return OPENDAP + urllib.parse.quote(gid) + ".dap.nc4?dap4.ce=" + (urllib.parse.quote(ce, safe="") if enc else
                                                                       urllib.parse.quote(ce, safe="/[]:;"))


def probe_coords(gid, y, x, token, half, enc):
    """plan 用：只取 cellat/cellon/time（不取盐度），报返回形状、有值格数与中心格坐标／时刻。"""
    np = _np()
    buf = get_retry(coord_url(gid, y, x, half, enc), token, 90, attempts=2)
    sub = parse_nc(buf, ("cellat", "cellon", "time"))
    la, lo, tt = sub["cellat"], sub["cellon"], sub["time"]
    c = half
    return {"bytes": len(buf), "shape": list(la.shape), "dims": sub["cellat__dims"],
            "n_finite_cellat": int(np.isfinite(la).sum()),
            "center_cellat": [rnd(v, 4) for v in la[c, c, :]], "center_cellon": [rnd(v, 4) for v in lo[c, c, :]],
            "center_time_utc": [iso(v + EPOCH2000) if np.isfinite(v) else None for v in tt[c, c, :]]}


def parse_nc(buf, names):
    """DAP4 nc4 响应（字节）→ {名: float 数组}，另给 {名}__dims。"""
    import netCDF4
    np = _np()
    out = {}
    with netCDF4.Dataset("p5mem.nc", mode="r", memory=bytes(buf)) as ds:
        pool = dict(ds.variables)
        for g in ds.groups.values():
            pool.update(g.variables)
        for n in names:
            v = pool.get(n)
            if v is None:
                raise DataError(f"子集缺变量 {n}")
            v.set_auto_maskandscale(True)
            out[n] = np.ma.filled(np.ma.asarray(v[:]).astype(float), np.nan)
            out[n + "__dims"] = [d.split("/")[-1] for d in v.dimensions]
    return out


def check_dims(sub):
    for v in VARS3:
        if v in sub and tuple(sub[v + "__dims"]) != DIMS3:
            raise ValueError(f"V3 维序不符：{v} {sub[v + '__dims']} ≠ {DIMS3}")
        if v in sub and sub[v].shape != (3, 3, 2):
            raise ValueError(f"V3 形状不符：{v} {sub[v].shape}")
    for v in VARS2:
        if v in sub and tuple(sub[v + "__dims"]) != DIMS3[:2]:
            raise ValueError(f"V3 维序不符：{v} {sub[v + '__dims']}")


def dmr_dims(xml_bytes, var):
    import xml.etree.ElementTree as ET
    root = ET.fromstring(xml_bytes)
    for el in root.iter():
        if el.attrib.get("name") == var:
            return [c.attrib.get("name", "").split("/")[-1] for c in el if c.tag.endswith("Dim")]
    return None


def save_subset(path, sub, meta):
    np = _np()
    arrs = {k: v for k, v in sub.items() if not k.endswith("__dims")}
    tmp = path + ".part.npz"
    np.savez(tmp, meta=json.dumps(meta), **arrs)
    os.replace(tmp, path)


def load_subset(path):
    np = _np()
    with np.load(path, allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    d["meta"] = json.loads(str(d["meta"]))
    return d


def looks_info(sub):
    """每视向：中心格时刻（unix）、是否覆盖（中心格 sss 非缺）。"""
    np = _np()
    out = []
    for lk in (0, 1):
        t = sub["time"][1, 1, lk]
        cov = bool(np.isfinite(sub["sss_smap_40km"][1, 1, lk]))
        out.append({"t": float(t) + EPOCH2000 if np.isfinite(t) else NAN, "covered": cov and bool(np.isfinite(t))})
    return out


def qc_valid(sss, iqc, wind2d, exclude):
    """V5：逐格有效掩码（3×3）。"""
    np = _np()
    q = np.where(np.isfinite(iqc), iqc, -1).astype(np.int64)
    ok = np.isfinite(sss) & (sss >= SSS_RANGE[0]) & (sss <= SSS_RANGE[1]) & (q >= 0)
    for b in tuple(QC_INVALID) + tuple(exclude):
        ok &= ((q >> b) & 1) == 0
    ok &= ~(np.isfinite(wind2d) & (wind2d > WIND_MAX))
    return ok


# ======================================================================== CMORPH 窗口（V20）
class WinExtractor:
    def __init__(self, stations, cache_root):
        self.st = {n: (la, lo % 360.0) for n, (la, lo) in stations.items()}
        self.root = os.path.join(cache_root, "cmorph")
        self.idx = None
        self.shape = None

    def _index(self, ds):
        np = _np()
        lat = np.asarray(ds["lat"][:], float)
        lon = np.asarray(ds["lon"][:], float)
        self.shape = (len(lat), len(lon))
        self.idx = {}
        for name, (la, lo) in self.st.items():
            i = int(np.argmin(np.abs(lat - la)))
            j = int(np.argmin(np.abs(((lon - lo) + 180.0) % 360.0 - 180.0)))
            nx = int(math.ceil(WIN_LAT_PIX / max(math.cos(math.radians(la)), 0.2)))
            rows = list(range(max(i - WIN_LAT_PIX, 0), min(i + WIN_LAT_PIX + 1, len(lat))))
            cols = [(j + k) % len(lon) for k in range(-nx, nx + 1)]
            self.idx[name] = (rows, cols)
            d = os.path.join(self.root, safe(name))
            os.makedirs(d, exist_ok=True)
            g = {"lat": [float(lat[r]) for r in rows], "lon": [float(lon[c]) for c in cols],
                 "center": [rows.index(i), cols.index(j)], "station": [la, lo]}
            gp = os.path.join(d, "grid.json")
            if os.path.exists(gp):
                old = json.load(open(gp))
                if old["lat"] != g["lat"] or old["lon"] != g["lon"]:
                    raise ValueError(f"CMORPH 网格与已存 grid.json 不一致：{name}")
            else:
                jdump(g, gp)

    def extract(self, path, h, names):
        import netCDF4
        np = _np()
        with netCDF4.Dataset(path) as ds:
            t = [int(x) for x in ds["time"][:]]
            if t != [h * 3600, h * 3600 + 1800]:
                return f"time_mismatch:{t}"
            shp = (ds.dimensions["lat"].size, ds.dimensions["lon"].size)
            if self.idx is None or shp != self.shape:
                self._index(ds)
            v = ds["cmorph"]
            v.set_auto_maskandscale(True)
            for name in names:
                rows, cols = self.idx[name]
                blk = np.ma.filled(np.ma.asarray(v[:, rows[0]:rows[-1] + 1, cols]).astype(np.float32), np.nan)
                d = os.path.join(self.root, safe(name))
                tmp = os.path.join(d, f"{h}.part.npy")
                np.save(tmp, blk)
                os.replace(tmp, os.path.join(d, f"{h}.npy"))
        return None


def cmorph_missing(cache_root):
    p = os.path.join(cache_root, "cmorph", "missing.jsonl")
    out = {}
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            for ln in f:
                try:
                    r = json.loads(ln)
                    out[int(r["hour"])] = r["why"]
                except Exception:
                    continue
    return out


def run_cmorph(need, cache_root, stations, out_dir, log, max_conc, raw=None):
    """need：{站: set(hours)}。缺窗口的小时：raw 档有原件就直接抽，没有就限速下载进 raw 档再抽；原件一律保留（V25）。"""
    miss = cmorph_missing(cache_root)
    by_hour = {}
    for st, hs in need.items():
        for h in hs:
            if not os.path.exists(os.path.join(cache_root, "cmorph", safe(st), f"{h}.npy")):
                by_hour.setdefault(h, []).append(st)
    todo = sorted(h for h in by_hour if h not in miss)
    all_h = set().union(*need.values()) if need else set()
    os.makedirs(os.path.join(cache_root, "cmorph"), exist_ok=True)
    raw = raw or {}
    fx = rs.CmorphRaw(raw.get("root", rs.RAW_ROOT_DEFAULT), raw.get("rate_mbps", 2.0), owner="p5",
                      source_job=raw.get("tag", "p5-sss"), exc_cls=p2.DataSourceError, heartbeat=True)
    rs.heartbeat_write(fx.root, "active")
    ex_ = WinExtractor(stations, cache_root)
    t0 = time.monotonic()
    n_done, failed, consec = 0, [], 0
    mpath = os.path.join(cache_root, "cmorph", "missing.jsonl")
    msg = f"[CMORPH] 需 {len(all_h)} 个唯一小时（站×小时 {sum(len(v) for v in need.values())}），待下 {len(todo)}"
    print(msg, flush=True)
    log.log(msg)

    def add_missing(h, why):
        with open(mpath, "a", encoding="utf-8") as f:
            f.write(json.dumps({"hour": h, "why": why}) + "\n")

    def handle(h, res):
        nonlocal n_done
        kind, val = res
        if kind == "missing":
            add_missing(h, val)
            return
        try:
            why = ex_.extract(val, h, by_hour[h])
        except Exception as e:
            why = f"read_error:{type(e).__name__}"
        if why:
            add_missing(h, why)
            log.log(f"CMORPH {h} 记缺测：{why}")
        n_done += 1
        if n_done % 250 == 0:
            el = time.monotonic() - t0
            rate = fx.bytes / 1e6 / max(el, 1e-6)
            m = f"[CMORPH] {n_done}/{len(todo)}，{fx.bytes / 1e9:.2f} GB，{rate:.2f} MB/s，剩余约 " \
                f"{(len(todo) - n_done) * p2.CMORPH_MEAN_MB / max(rate, 1e-6) / 3600:.2f} h"
            print(m, flush=True)
            log.log(m)

    it = iter(todo)
    with cf.ThreadPoolExecutor(max_workers=max_conc) as ex:
        pend = {}
        for h in it:
            pend[ex.submit(fx.fetch, h)] = h
            if len(pend) >= max_conc * 2:
                break
        while pend:
            done, _ = cf.wait(pend, return_when=cf.FIRST_COMPLETED)
            for fut in done:
                h = pend.pop(fut)
                try:
                    res = fut.result()
                except p2.DataSourceError as e:
                    failed.append(h)
                    consec += 1
                    log.log(f"CMORPH 失败 {h}：{e}")
                    if consec >= p2.MAX_CONSEC_FAIL:
                        for f in pend:
                            f.cancel()
                        raise DataError(f"CMORPH 连续 {consec} 个文件失败（重跑即续传）")
                    continue
                consec = 0
                handle(h, res)
            for h in it:
                pend[ex.submit(fx.fetch, h)] = h
                if len(pend) >= max_conc * 2:
                    break
    still = []
    for h in failed:
        try:
            handle(h, fx.fetch(h))
        except p2.DataSourceError as e:
            still.append(h)
            log.log(f"CMORPH 末轮仍失败 {h}：{e}")
    rs.heartbeat_write(fx.root, "done")
    if still:
        raise DataError(f"CMORPH 末轮仍失败 {len(still)} 小时（重跑即续传）")
    return {"unique_hours_needed": len(all_h), "station_hours_needed": sum(len(v) for v in need.values()),
            "downloaded_this_run": n_done, "bytes_this_run": fx.bytes, "http_requests": fx.n_req,
            "raw_root": fx.root, "raw_downloaded": fx.downloaded, "raw_reused": fx.reused, "rate_mbps": raw.get("rate_mbps", 2.0),
            "seconds": round(time.monotonic() - t0, 1), "missing_total": len(cmorph_missing(cache_root))}


class CmCache:
    """分析段读 CMORPH 站点窗口（LRU）。"""

    def __init__(self, cache_root, maxn=600):
        self.root = os.path.join(cache_root, "cmorph")
        self.grid = {}
        self.lru = {}
        self.order = []
        self.maxn = maxn
        self.missing = cmorph_missing(cache_root)

    def g(self, st):
        if st not in self.grid:
            p = os.path.join(self.root, safe(st), "grid.json")
            self.grid[st] = json.load(open(p)) if os.path.exists(p) else None
        return self.grid[st]

    def hour(self, st, h):
        np = _np()
        k = (st, h)
        if k in self.lru:
            return self.lru[k]
        p = os.path.join(self.root, safe(st), f"{h}.npy")
        a = np.load(p).astype(float) if os.path.exists(p) else None
        self.lru[k] = a
        self.order.append(k)
        if len(self.order) > self.maxn:
            self.lru.pop(self.order.pop(0), None)
        return a

    def halfsteps(self, st, q0, q1):
        """半步 q0..q1 → (n, ny, nx)；任一小时缺 → None。"""
        np = _np()
        out = []
        for q in range(q0, q1 + 1):
            a = self.hour(st, q // 2)
            if a is None:
                return None
            out.append(a[q % 2])
        return np.stack(out)


# ======================================================================== 足印核与 RIM（V8–V11、V19）
_KCACHE = {}


def cell_kernel(grid, clat, clon):
    np = _np()
    clat, clon = round(float(clat), 4), round(float(clon) % 360.0, 4)
    key = (tuple(grid["lat"][:2]), tuple(grid["lon"][:2]), len(grid["lat"]), len(grid["lon"]), clat, clon)
    if key in _KCACHE:
        return _KCACHE[key]
    lat = np.asarray(grid["lat"], float)[:, None]
    lon = np.asarray(grid["lon"], float)[None, :]
    dy = (lat - clat) * KM_PER_DEG
    dx = (((lon - clon) + 180.0) % 360.0 - 180.0) * KM_PER_DEG * math.cos(math.radians(clat))
    d2 = dy ** 2 + dx ** 2
    w = np.exp(-4.0 * math.log(2.0) * d2 / FWHM_KM ** 2) * np.cos(np.radians(lat))
    w = np.where(d2 <= TRUNC_KM ** 2, w, 0.0)
    s = w.sum()
    w = w / s if s > 0 else w * np.nan
    _KCACHE[key] = w
    return w


def vec_rim0(P, U):
    """P：(49, N) 半步雨强（q−48..q）；U：(49,)。返回 (N,) 的 z=0 稀释因子 F(q)，与 p2.rim_factor(P[:,n], U, 0)[0] 同式。"""
    np = _np()
    P = np.asarray(P, float)
    U = np.asarray(U, float)
    tc = np.array(p2.TC, float)[:, None]
    ti = np.array(p2.TI, float)[:, None]
    Pw = P[:48]
    Uw = np.broadcast_to(U[:48][:, None], Pw.shape)
    irr = Pw / 1000.0 / 3600.0
    kz = p2.KZ_COEF * Uw ** 2
    d0 = p2.d0_interp(Uw, Pw)
    Pc = P[48]
    Uc = np.full(Pc.shape, U[48])
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        term = p2.C1 * irr * tc / np.sqrt(kz * ti)
        prior = d0 / (d0 + term)
        irrc = Pc / 1000.0 / 3600.0
        kzc = p2.KZ_COEF * Uc ** 2
        d0c = p2.d0_interp(Uc, Pc)
        cterm = p2.C2 * irrc * 1800.0 / np.sqrt(kzc * p2.T_CURRENT_S)
        cur = d0c / (d0c + cterm)
    return np.prod(prior, axis=0) * cur


def wind_hourly(stn, h_lo, h_hi):
    """V10：小时 h_lo..h_hi（含）的 10 m 风；有缺 → None。"""
    np = _np()
    ws = stn["ser"]["wind"]
    i0 = h_lo - stn["h0"]
    n = h_hi - h_lo + 1
    pad = p2.WIND_GAP_MAX_H
    lo, hi = i0 - pad, i0 + n + pad
    seg = np.full(hi - lo, np.nan)
    a, b = max(lo, 0), min(hi, len(ws))
    if b > a:
        seg[a - lo:b - lo] = ws[a:b]
    w = p2.fill_wind(seg, pad)[pad:pad + n] * _rp.wind_factor(stn["name"], p2.WIND_FACTOR, p2.WIND_Z0_M)  # [repro] 风高开关 [options] wind_height；缺省 uniform_4m 时 _rp.wind_factor 原样返回原系数
    if not np.all(np.isfinite(w)):
        return None
    return np.maximum(w, p2.WIND_FLOOR)


def wind_halfsteps(stn, q0, q1):
    np = _np()
    h_lo, h_hi = q0 // 2, q1 // 2
    w = wind_hourly(stn, h_lo, h_hi)
    if w is None:
        return None
    return np.array([w[q // 2 - h_lo] for q in range(q0, q1 + 1)], float)


def weighted_F(Fp, w):
    """V9 丢弃规则：像元 F NaN 的权重 >5% → NaN；否则重归一。"""
    np = _np()
    bad = ~np.isfinite(Fp)
    wb = float(w[bad].sum())
    if wb > DROP_MAX:
        return NAN, wb
    ws = float(w[~bad].sum())
    return (float((w[~bad] * Fp[~bad]).sum() / ws) if ws > 0 else NAN), wb


class PassEval:
    """逐过境、逐视向的卫星值与足印模型值（带缓存）。"""

    def __init__(self, stations, cm):
        self.stations = stations
        self.cm = cm
        self.fpix = {}
        self.memo = {}

    def pix_F(self, st, gid, lk, q):
        """该视向半步 q 的全窗口 z=0 F（只算超集核内像元），缓存。返回 (F2d, U) 或 (None, 原因)。"""
        np = _np()
        k = (st, gid, lk)
        if k in self.fpix:
            return self.fpix[k]
        g = self.cm.g(st)
        if g is None:
            self.fpix[k] = (None, "cmorph_grid_missing")
            return self.fpix[k]
        Pst = self.cm.halfsteps(st, q - 48, q)
        if Pst is None:
            self.fpix[k] = (None, "cmorph_missing")
            return self.fpix[k]
        U = wind_halfsteps(self.stations[st], q - 48, q)
        if U is None:
            self.fpix[k] = (None, "wind_gap")
            return self.fpix[k]
        clat, clon = station_cell(*self.stations[st]["_coord"])
        lat = np.asarray(g["lat"], float)[:, None]
        lon = np.asarray(g["lon"], float)[None, :]
        d2 = ((lat - clat) * KM_PER_DEG) ** 2 + ((((lon - clon) + 180.0) % 360.0 - 180.0) * KM_PER_DEG
                                                 * math.cos(math.radians(clat))) ** 2
        ci, cj = g["center"]
        mask = np.broadcast_to(d2 <= SUPERSET_KM ** 2, Pst.shape[1:]).copy()
        mask[ci, cj] = True
        F = np.full(Pst.shape[1:], np.nan)
        F[mask] = vec_rim0(Pst[:, mask], U)
        self.fpix[k] = ((F, Pst, U), None)
        return self.fpix[k]

    def look(self, st, gid, sub, lk, v):
        np = _np()
        sss = sub["sss_smap_40km"][:, :, lk]
        cla, clo = sub["cellat"][:, :, lk], sub["cellon"][:, :, lk]
        ok = qc_valid(sss, sub["iqc_flag"][:, :, lk], sub["winspd"], v["qc_ex"]) & np.isfinite(cla) & np.isfinite(clo)
        t = sub["time"][1, 1, lk]
        if not np.isfinite(t):
            return None, "no_time"
        if not ok[1, 1]:
            return None, "center_qc"
        if v["agg"] == "A70":
            if int(ok.sum()) < A70_MIN_VALID:
                return None, "a70_few_valid"
            cells = [(i, j) for i in range(3) for j in range(3) if ok[i, j]]
        else:
            cells = [(1, 1)]
        sat = float(np.mean([sss[i, j] for i, j in cells]))
        unc_c = [sub["sss_smap_40km_unc"][i, j, lk] for i, j in cells]
        unc_c = [u for u in unc_c if np.isfinite(u)]
        sig = float(np.mean(unc_c) / math.sqrt(len(cells))) if unc_c else NAN
        ts = float(t) + EPOCH2000
        q = int(ts // 1800)
        res, why = self.pix_F(st, gid, lk, q)
        if res is None:
            return None, why
        F2d, Pst, U = res
        g = self.cm.g(st)
        w = np.zeros(F2d.shape)
        for i, j in cells:
            w += cell_kernel(g, cla[i, j], clo[i, j])
        w /= len(cells)
        if v["drive"] == "pixel":
            Fbar, dropped = weighted_F(F2d, w)
        else:
            m = w > 0
            Pm = Pst[:, m]
            wm = w[m]
            bad = ~np.all(np.isfinite(Pm), axis=0)
            dropped = float(wm[bad].sum())
            if dropped > DROP_MAX:
                Fbar = NAN
            else:
                Pbar = (Pm[:, ~bad] * wm[~bad]).sum(axis=1) / wm[~bad].sum()
                Fbar = float(p2.rim_factor(Pbar, U, 0.0)[0])
        if not isnum(Fbar):
            return None, "model_nan"
        ci, cj = g["center"]
        mfoot = w > 0
        rain_foot = float((Pst[-1][mfoot] * w[mfoot]).sum()) if np.all(np.isfinite(Pst[-1][mfoot])) else NAN
        return {"t": ts, "q": q, "sat": sat, "F": Fbar, "Fc": float(F2d[ci, cj]), "sig": sig, "n_cells": len(cells),
                "dropped_w": dropped, "rain_foot_cmorph": rain_foot}, None

    def evaluate(self, st, gid, sub, v):
        """一次过境（两视向合并）。返回 dict 或 (None, 原因)。"""
        np = _np()
        key = (st, gid, v["agg"], tuple(v["qc_ex"]), v["drive"])
        if key in self.memo:
            return self.memo[key]
        looks, whys = [], []
        for lk in (0, 1):
            r, why = self.look(st, gid, sub, lk, v)
            if r is None:
                whys.append(why)
            else:
                looks.append(r)
        if not looks:
            out = (None, whys[0] if whys else "no_look")
        else:
            t = float(np.mean([r["t"] for r in looks]))
            lon = self.stations[st]["_coord"][1]
            lst = ((t / 3600.0) + lon180(lon) / 15.0) % 24.0
            sref = sub["sss_ref"][1, 1]
            out = ({"t": t, "sat": float(np.mean([r["sat"] for r in looks])), "F": float(np.mean([r["F"] for r in looks])),
                    "Fc": float(np.mean([r["Fc"] for r in looks])),
                    "sig": float(math.sqrt(np.mean([r["sig"] ** 2 for r in looks]) / len(looks)))
                    if all(isnum(r["sig"]) for r in looks) else NAN,
                    "n_looks": len(looks), "node": "am" if lst < 12.0 else "pm", "lst": lst,
                    "sss_ref": float(sref) if np.isfinite(sref) else NAN,
                    "imerg_c": float(sub["rain"][1, 1]) if np.isfinite(sub["rain"][1, 1]) else NAN,
                    "rain_foot_cmorph": float(np.nanmean([r["rain_foot_cmorph"] for r in looks]))
                    if any(isnum(r["rain_foot_cmorph"]) for r in looks) else NAN,
                    "q": looks[0]["q"], "gid": gid, "dropped_w": max(r["dropped_w"] for r in looks)}, None)
        self.memo[key] = out
        return out


# ======================================================================== 事件级（V7、V11、V12）
def mean_near(arr, idx, half):
    np = _np()
    a = arr[max(0, idx - half):idx + half + 1]
    a = a[np.isfinite(a)]
    return float(a.mean()) if len(a) else NAN


def mooring_third(stn, e, post, cm, st):
    """V12：同一雨后过境时刻的系泊观测与站点像元 RIM-3。"""
    np = _np()
    ser = stn["ser"]
    i = e["i"]
    H = e["hour"]
    hs = int(post["t"] // 3600)
    q = int(post["t"] // 1800)
    ix = hs - stn["h0"]
    out = {}
    for key, lab in (("s1", "1"), ("sss05", "05")):
        base = p2.pre_median(ser[key], i)
        v = mean_near(ser[key], ix, OBS_HALF_H[key]) if 0 <= ix < len(ser[key]) else NAN
        out[f"obs{lab}"] = v - base if isnum(v) and isnum(base) else NAN
    g = cm.g(st)
    q0 = 2 * H - 60
    Pst = cm.halfsteps(st, q0, q) if g is not None else None
    U = wind_halfsteps(stn, q0, q)
    z1 = p2.model_depths(e, st)[1.0]
    out["s1_depth"] = z1
    if Pst is None or U is None:
        out.update({"rim1": NAN, "rim05": NAN, "why": "cmorph_or_wind"})
        return out
    ci, cj = g["center"]
    Pc = Pst[:, ci, cj]
    for z, lab, key in ((z1, "1", "s1"), (0.5, "05", "sss05")):
        F = p2.rim_factor(Pc, U, z)                     # 半步 2H−12 .. q
        Fpre = F[:12].reshape(-1, 2).mean(axis=1)
        s0 = p2.pre_median(ser[key], i)
        d = float(F[-1] - np.median(Fpre)) if np.all(np.isfinite(Fpre)) and np.isfinite(F[-1]) else NAN
        out[f"rim{lab}"] = s0 * d if isnum(s0) and isnum(d) else NAN
    return out


def event_eval(e, stn, st, passes, pe, v, cm, third=False):
    """passes：该站 [(t_cov, gid, sub)]（按时刻）。返回 (行 dict, None) 或 (None, 原因)。"""
    np = _np()
    t0 = e["hour"] * 3600.0
    plo, phi = v["post"]
    rlo, rhi = v["ref"]
    post, refs, why_post = None, [], None
    for tc, gid, sub in passes:
        dt = (tc - t0) / 3600.0
        if dt < rlo - 1.0 or dt >= phi + 1.0:
            continue
        r, why = pe.evaluate(st, gid, sub, v)
        if r is None:
            if plo <= dt < phi:
                why_post = why_post or why
            continue
        dt = (r["t"] - t0) / 3600.0
        if plo <= dt < phi and (post is None or r["t"] < post["t"]):
            post = r
        elif rlo <= dt < rhi:
            refs.append(r)
    if post is None:
        return None, f"no_valid_post:{why_post or 'none'}"
    if v["same_node"]:
        refs = [r for r in refs if r["node"] == post["node"]]
    if not refs:
        return None, "no_valid_ref"
    sat_ref = float(np.mean([r["sat"] for r in refs]))
    F_ref = float(np.mean([r["F"] for r in refs]))
    Fc_ref = float(np.mean([r["Fc"] for r in refs]))
    s0 = post["sss_ref"] if isnum(post["sss_ref"]) else sat_ref
    sig = math.sqrt(post["sig"] ** 2 + np.mean([r["sig"] ** 2 for r in refs]) / len(refs)) \
        if isnum(post["sig"]) and all(isnum(r["sig"]) for r in refs) else NAN
    row = {"station": st, "season": e["season"], "onset_utc": e["onset_utc"], "hour": e["hour"],
           "post_gid": post["gid"], "post_utc": iso(post["t"]), "dt_post_h": round((post["t"] - t0) / 3600.0, 3),
           "n_ref": len(refs), "node": post["node"], "n_looks_post": post["n_looks"],
           "sat_post": post["sat"], "sat_ref": sat_ref, "y": post["sat"] - sat_ref,
           "F_post": post["F"], "F_ref": F_ref, "S0": s0, "x": s0 * (post["F"] - F_ref),
           "x_pix0": s0 * (post["Fc"] - Fc_ref), "sig": sig,
           "imerg_c_post": post["imerg_c"], "rain_foot_cmorph_post": post["rain_foot_cmorph"],
           "dropped_w_post": post["dropped_w"]}
    if third:
        row.update(mooring_third(stn, e, post, cm, st))
    return row, None


# ======================================================================== 统计（V13–V16）
def boot_block(rows, num, den, pairs_extra=None):
    """站×季与站两种簇的 bootstrap；返回 (点, 站×季数组 dict, 站数组 dict, 簇数)。"""
    np = _np()
    ks = [(r["station"], r["season"]) for r in rows]
    kt = [r["station"] for r in rows]
    a = np.array([r[num] for r in rows], float)
    b = np.array([r[den] for r in rows], float)
    pairs = {"R": (a, b), "den": (b, np.ones_like(b)), "slope": (a * b, b * b)}
    for name, (f, g) in (pairs_extra or {}).items():
        pairs[name] = (np.array([f(r) for r in rows], float), np.array([g(r) for r in rows], float))
    bs, Kss = p2.cluster_boot(ks, pairs, B, SEED)
    bt, Kst = p2.cluster_boot(kt, pairs, B, SEED)
    return bs, bt, Kss, Kst


def ratio_block(rows, num="y", den="x", label="", pairs_extra=None):
    np = _np()
    rows = [r for r in rows if isnum(r.get(num)) and isnum(r.get(den))]
    n = len(rows)
    blk = {"label": label, "n": n, "stations": len({r["station"] for r in rows})}
    if n < 3:
        blk.update({"evaluable": False, "reason": "事件不足"})
        return blk, None
    sn = sum(r[num] for r in rows)
    sd = sum(r[den] for r in rows)
    bs, bt, Kss, Kst = boot_block(rows, num, den, pairs_extra)
    ci_ss, _ = pm.ci95(bs["R"])
    ci_st, _ = pm.ci95(bt["R"])
    frac = float(np.mean(bs["den"] >= 0))
    jk = pm.jackknife_station(rows, lambda rr: (sum(r[num] for r in rr) / sum(r[den] for r in rr))
                              if rr and sum(r[den] for r in rr) != 0 else NAN)
    slope = sum(r[num] * r[den] for r in rows) / sum(r[den] ** 2 for r in rows)
    blk.update({"clusters_station_season": Kss, "clusters_station": Kst, "sum_num": rnd(sn), "sum_den": rnd(sd),
                "R": rnd(sn / sd) if sd else None, "ci95_station_season": ci_ss, "ci95_station": ci_st,
                "frac_boot_den_nonneg": rnd(frac, 4), "jackknife_station": jk,
                "slope_through_origin": rnd(slope), "slope_ci95_station_season": pm.ci95(bs["slope"])[0],
                "evaluable": bool(n >= MIN_N and Kss >= MIN_G and frac <= DEN_NONNEG_MAX)})
    if not blk["evaluable"]:
        blk["reason"] = f"n={n}（需≥{MIN_N}）、簇={Kss}（需≥{MIN_G}）、Σx≥0 抽样比例={frac:.4f}（需≤{DEN_NONNEG_MAX}）"
    return blk, (bs, bt)


def judge(blk):
    """V15：按序取第一条成立者；另加稳健标注。"""
    if not blk.get("evaluable"):
        return {"category": "不可评", "robust": None, "rule": "第六节第 1 条"}
    lo, hi = blk["ci95_station_season"]
    if lo is None or hi is None:
        return {"category": "不可评", "robust": None, "rule": "第六节第 1 条（bootstrap 区间不可算）"}
    jk = blk.get("jackknife_station") or {}
    loo = jk.get("loo_range") if jk.get("evaluable") else None
    st_lo, st_hi = blk["ci95_station"]
    if hi < THRESH:
        if hi < 0:
            cat, rule = "卫星未见淡化（异常，须排查 QC 与参照）：不支持 R2", "第六节第 2 条（上端 <0）"
        else:
            cat, rule = "支持 R1（海面也显著小于 RIM-3）", "第六节第 2 条"
        robust = bool(st_hi is not None and st_hi < THRESH and loo is not None and loo[1] < THRESH)
    elif lo > THRESH:
        cat, rule = "支持 R2（海面与 RIM-3 相符或更大）", "第六节第 3 条"
        robust = bool(st_lo is not None and st_lo > THRESH and loo is not None and loo[0] > THRESH)
    else:
        return {"category": "不能区分", "robust": None, "rule": "第六节第 4 条"}
    return {"category": cat, "robust": "稳健" if robust else "依赖簇定义", "rule": rule}


def power_block(rows):
    """V16：只用模型 x 与形式不确定度；在任何 ΣΔSSS_sat 之前写盘。"""
    np = _np()
    xs = np.array([r["x"] for r in rows], float)
    sg = np.array([r["sig"] for r in rows], float)
    ok = np.isfinite(xs) & np.isfinite(sg)
    out = {"n": int(len(rows)), "n_with_formal_sigma": int(ok.sum()),
           "clusters_station_season": len({(r["station"], r["season"]) for r in rows}),
           "stations": sorted({r["station"] for r in rows}),
           "sum_x": rnd(float(np.nansum(xs))), "mean_x": rnd(float(np.nanmean(xs))) if len(xs) else None,
           "rms_x": rnd(float(np.sqrt(np.nanmean(xs ** 2)))) if len(xs) else None,
           "median_formal_sigma_event": rnd(float(np.nanmedian(sg))) if ok.any() else None,
           "note": "形式不确定度按独立近似，未含地球物理代表性误差与参照期真实变化，所以偏乐观；只作信息，不改判读"}
    if ok.sum() >= 3 and abs(xs[ok].sum()) > 0:
        se = float(math.sqrt((sg[ok] ** 2).sum()) / abs(xs[ok].sum()))
        hw = 1.959964 * se
        out.update({"expected_se_R_formal": rnd(se), "expected_ci_halfwidth_formal": rnd(hw),
                    "can_exclude_065_if_R_true_0.3": bool(0.3 + hw < THRESH),
                    "can_exclude_065_if_R_true_1.0": bool(1.0 - hw > THRESH),
                    "min_resolvable_abs_diff_from_065": rnd(hw)})
    return out


# ======================================================================== 主流程
def build(args, out_dir, log):
    stations, ev, val, _prim, _extra, _plan = p2.build_all(args, out_dir, log)
    sb = {}
    for s in stations:
        s["_coord"] = (float(s["lat"]), float(s["lon"]))
        sb[s["name"]] = s
    events = [e for e in ev[p2.PRIMARY_THRESHOLD] if e["hour"] >= ERA_START_H]
    for e in events:
        e["_st"] = e["_stn"]
    info = {"events_646": len(ev[p2.PRIMARY_THRESHOLD]), "events_smap_era": len(events),
            "rebuild_d4": val.get("d4_count") if isinstance(val, dict) else None,
            "stations_coords": {n: list(s["_coord"]) for n, s in sb.items()}}
    return sb, events, info


def events_by_station(events):
    out = {}
    for e in events:
        out.setdefault(e["_st"], []).append(e)
    return out


def cmr_for_station(st, coord, evs, cache_root, log, use_cache=True):
    path = os.path.join(cache_root, "cmr", f"{safe(st)}.json") if cache_root else None
    t0 = min(e["hour"] for e in evs) * 3600.0 - 4 * 86400.0
    t1 = max(e["hour"] for e in evs) * 3600.0 + 2 * 86400.0
    if use_cache and path and os.path.exists(path):
        d = json.load(open(path))
        if d.get("span") == [t0, t1]:
            return [tuple(g) for g in d["granules"]]
    gs = cmr_granules(COLL, coord[0], coord[1], t0, t1)
    if path:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        jdump({"span": [t0, t1], "coord": list(coord), "granules": [list(g) for g in gs]}, path)
    log.log(f"CMR {st}：{len(gs)} 轨")
    return gs


def candidates(sb, evs_by_st, cache_root, log):
    """V2：每站候选轨 {站: {gid: t_approx}}；另给估计过境时刻（plan 计数用）。"""
    cand, approx = {}, {}
    for st, evs in evs_by_st.items():
        coord = sb[st]["_coord"]
        gs = cmr_for_station(st, coord, evs, cache_root, log)
        pts = [(approx_pass(g[1], g[2], coord[1]), g[0]) for g in gs]
        approx[st] = sorted(pts)
        keep = {}
        for e in evs:
            t0 = e["hour"] * 3600.0
            for tp, gid in pts:
                if CAND_H[0] <= (tp - t0) / 3600.0 < CAND_H[1]:
                    keep[gid] = tp
        cand[st] = keep
    return cand, approx


def count_windows(evs_by_st, approx):
    out = {}
    for wname, (plo, phi) in (("post6", (0, 6)), ("post12", (0, 12)), ("post24", (0, 24))):
        for rname, (rlo, rhi) in (("ref72", (-72, 0)), ("ref48", (-48, 0))):
            n, cl, sts = 0, set(), {}
            for st, evs in evs_by_st.items():
                ts = [t for t, _ in approx.get(st, [])]
                for e in evs:
                    t0 = e["hour"] * 3600.0
                    d = [(t - t0) / 3600.0 for t in ts]
                    if any(plo <= x < phi for x in d) and any(rlo <= x < rhi for x in d):
                        n += 1
                        cl.add((st, e["season"]))
                        sts[st] = sts.get(st, 0) + 1
            out[f"{wname}&{rname}"] = {"events": n, "clusters": len(cl), "stations": sts}
    return out


def hours_needed(evs_by_st, passes_by_st):
    """V21：passes_by_st：{站: [(t_min, t_max, gid)]}（覆盖视向时刻）。"""
    need = {}
    for st, evs in evs_by_st.items():
        ps = passes_by_st.get(st, [])
        hs = set()
        for e in evs:
            H = e["hour"]
            t0 = H * 3600.0
            post = [p for p in ps if FETCH_POST_H[0] <= (p[0] - t0) / 3600.0 < FETCH_POST_H[1]]
            ref = [p for p in ps if FETCH_REF_H[0] <= (p[0] - t0) / 3600.0 < FETCH_REF_H[1]]
            if not post or not ref:
                continue
            for p in post:
                hs.update(range(H - 30, int(p[1] // 1800) // 2 + 1))
            for p in ref:
                qa, qb = int(p[0] // 1800), int(p[1] // 1800)
                hs.update(range((qa - 48) // 2, qb // 2 + 1))
        if hs:
            need[st] = hs
    return need


def run_plan(args, out_dir, log):
    sb, events, info = build(args, out_dir, log)
    ebs = events_by_station(events)
    cand, approx = candidates(sb, ebs, os.path.join(out_dir, "plan_cache"), log)
    counts = count_windows(ebs, approx)
    n_req = sum(len(v) for v in cand.values())
    pseudo = {st: [(t, t, g) for t, g in approx[st]] for st in approx}
    need = hours_needed(ebs, pseudo)
    uh = set().union(*need.values()) if need else set()
    aq_note = "Aquarius：2011-08–2015-06 有过境且有参照的事件 7 个／5 簇，不可评；本任务不重算"
    # 认证探测：.dmr 维序＋站点格坐标与时刻（不取盐度）
    token = read_token(args.token_file)
    probe = {"token_file": args.token_file, "token_read": True}
    st0 = sorted(cand, key=lambda s: -len(cand[s]))[0]
    gid0 = sorted(cand[st0])[0]
    xml = get_retry(OPENDAP + urllib.parse.quote(gid0) + ".dmr", token, 120)
    probe["dmr_granule"] = gid0
    probe["dmr_dims"] = {v: dmr_dims(xml, v) for v in VARS3 + VARS2}
    probe["dmr_dims_ok"] = all(tuple(probe["dmr_dims"][v] or ()) == DIMS3 for v in VARS3) and \
        all(tuple(probe["dmr_dims"][v] or ()) == DIMS3[:2] for v in VARS2)
    coords = {}
    for st in sorted(cand):
        if not cand[st]:
            continue
        gid = sorted(cand[st], key=lambda g: cand[st][g])[len(cand[st]) // 2]
        y, _ = cell_yx(*sb[st]["_coord"])
        clat, clon = station_cell(*sb[st]["_coord"])
        c = {"granule": gid, "y": y, "station_cell": [clat, clon], "approx_pass_utc": iso(cand[st][gid])}
        try:
            x, info = locate_x(gid, y, clat, clon, token)
            c["locate"] = dict(info, x=x)
            if x is not None:
                c["h1_at_located_x"] = probe_coords(gid, y, x, token, 1, True)
                sub = parse_nc(get_retry(coord_url(gid, y, x), token, 90, attempts=2), ("cellat", "cellon", "time"))
                c["center_dist_deg"] = rnd(check_center(sub, clat, clon), 4)
        except (DataError, NotFound, ValueError) as e:
            c["error"] = redact(f"{type(e).__name__}: {e}")[:300]
        coords[st] = c
        jdump({"partial": True, "station_cells": coords}, os.path.join(out_dir, "p5_plan_probe_partial.json"))
        log.log(f"probe {st}：{json.dumps(c, ensure_ascii=False)[:600]}")
    probe["station_cells"] = coords
    plan = {"script": os.path.basename(__file__), "version": VERSION, "mode": "plan",
            "run_utc": iso(time.time()), "events": info, "cmr_window_counts_approx": counts,
            "candidate_granules": {st: len(v) for st, v in cand.items()}, "opendap_requests_est": n_req,
            "cmorph_hours_est": {"unique_hours": len(uh), "station_hours": sum(len(v) for v in need.values()),
                                 "est_GB": round(len(uh) * p2.CMORPH_MEAN_MB / 1000, 1),
                                 "est_hours_at_pca_rate": round(len(uh) / PCA_FILES_PER_S / 3600, 2),
                                 "basis": "按 CMR 估计过境时刻、假定候选轨全部覆盖（上界）；V21"},
            "aquarius": aq_note, "auth_probe": probe}
    jdump(plan, os.path.join(out_dir, "p5_plan.json"))
    print(json.dumps({k: plan[k] for k in ("cmr_window_counts_approx", "opendap_requests_est", "cmorph_hours_est")},
                     ensure_ascii=False)[:2000], flush=True)
    ok_loc = all(c.get("center_dist_deg") is not None for c in coords.values())
    return 0 if probe["dmr_dims_ok"] and ok_loc else 3


def load_passes(cache_root, st):
    d = os.path.join(cache_root, "smap", safe(st))
    out = []
    if not os.path.isdir(d):
        return out
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".npz") or fn.endswith(".part.npz"):
            continue
        sub = load_subset(os.path.join(d, fn))
        li = looks_info(sub)
        ts = [l["t"] for l in li if l["covered"]]
        if ts:
            out.append((min(ts), max(ts), sub["meta"]["gid"], sub))
    out.sort(key=lambda p: p[0])
    return out


def run_fetch(args, out_dir, log, sb, events):
    t_start = time.monotonic()
    token = read_token(args.token_file)
    cache = args.cache
    os.makedirs(cache, exist_ok=True)
    ebs = events_by_station(events)
    cand, _approx = candidates(sb, ebs, cache, log)
    # 维序核对（V3）
    st0 = sorted(cand, key=lambda s: -len(cand[s]))[0]
    xml = get_retry(OPENDAP + urllib.parse.quote(sorted(cand[st0])[0]) + ".dmr", token, 120)
    for v in VARS3:
        if tuple(dmr_dims(xml, v) or ()) != DIMS3:
            raise ValueError(f"V3 .dmr 维序不符：{v}")
    jobs = []
    for st, gids in cand.items():
        y, _ = cell_yx(*sb[st]["_coord"])
        clat, clon = station_cell(*sb[st]["_coord"])
        d = os.path.join(cache, "smap", safe(st))
        os.makedirs(d, exist_ok=True)
        for gid in sorted(gids):
            p = os.path.join(d, f"{gid}.npz")
            if not any(os.path.exists(p + ext) for ext in ("", ".missing", ".nocover")):
                jobs.append((st, gid, y, (clat, clon), p))
    msg = f"[SMAP] 候选 {sum(len(v) for v in cand.values())} 轨×站，待取 {len(jobs)}"
    print(msg, flush=True)
    log.log(msg)
    n_ok = n_miss = n_nocov = 0

    raw_root = rs.ensure_root(args.raw_root)
    raw_man = rs.Manifest(raw_root)

    def keep(st, gid):
        def _k(buf, what, url):
            rel = f"smap_opendap_subset/{safe(st)}/{gid}.{what}.nc4"
            dst = os.path.join(raw_root, rel)
            import repro_io as _rio  # [repro] never overwrite a kept original: save as .dup<time> (same rule as p7a keep())
            dst, rel = _rio.dup_if_exists(raw_root, rel)  # [repro] (os.replace below is reached only for a new name)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            tmp = dst + f".{os.getpid()}.part"
            with open(tmp, "wb") as f:
                f.write(buf)
            os.replace(tmp, dst)
            raw_man.add(dataset="smap_opendap_subset", url=url, path=rel, bytes=len(buf),
                        sha256=hashlib.sha256(buf).hexdigest(), status="ok", source_job=args.raw_tag)
        return _k

    def one(j):
        st, gid, y, (clat, clon), p = j
        try:
            x, info = locate_x(gid, y, clat, clon, token, keep=keep(st, gid))
            if x is None:
                return "nocover", j, info
            buf = get_retry(subset_url(gid, y, x), token, 180)
            keep(st, gid)(buf, f"sub_y{y}_x{x}", subset_url(gid, y, x))
        except NotFound:
            return "missing", j, None
        sub = parse_nc(buf, VARS3 + VARS2)
        check_dims(sub)
        info["center_dist_deg"] = rnd(check_center(sub, clat, clon), 4)
        save_subset(p, sub, dict({"station": st, "gid": gid, "y": y, "x": x}, **info))
        return "ok", j, info

    with cf.ThreadPoolExecutor(max_workers=args.max_conc) as ex:
        futs = [ex.submit(one, j) for j in jobs]
        for k, f in enumerate(cf.as_completed(futs), 1):
            kind, j, info = f.result()
            if kind == "missing":
                open(j[4] + ".missing", "w").write("404\n")
                n_miss += 1
            elif kind == "nocover":
                open(j[4] + ".nocover", "w").write(json.dumps(info) + "\n")
                n_nocov += 1
            else:
                n_ok += 1
            if k % 100 == 0:
                print(f"[SMAP] {k}/{len(jobs)}", flush=True)
    smap_s = round(time.monotonic() - t_start, 1)
    passes = {st: [(p[0], p[1], p[2]) for p in load_passes(cache, st)] for st in cand}
    need = hours_needed(ebs, passes)
    cm = run_cmorph(need, cache, {st: sb[st]["_coord"] for st in need}, out_dir, log, args.max_conc,
                    raw={"root": raw_root, "rate_mbps": args.rate_mbps, "tag": args.raw_tag})
    man = {"version": VERSION, "run_utc": iso(time.time()), "candidates": {st: len(v) for st, v in cand.items()},
           "smap_fetched_this_run": n_ok, "smap_missing_this_run": n_miss,
           "smap_nocover_this_run": n_nocov, "smap_seconds": smap_s,
           "passes_covered": {st: len(v) for st, v in passes.items()}, "cmorph": cm,
           "seconds_total": round(time.monotonic() - t_start, 1)}
    jdump(man, os.path.join(cache, "fetch_manifest.json"))
    jdump(man, os.path.join(out_dir, "p5_fetch.json"))
    return man


def variant(**kw):
    v = dict(PRIMARY)
    v.update(kw)
    return v


def compute_rows(sb, events, passes_by_st, pe, v, cm, third=False):
    rows, why = [], {}
    for e in events:
        st = e["_st"]
        if st in v["excl"]:
            continue
        ps = [(p[0], p[2], p[3]) for p in passes_by_st.get(st, [])]
        r, w = event_eval(e, sb[st], st, ps, pe, v, cm, third)
        if r is None:
            k = w.split(":")[0] if w.startswith("no_valid_ref") else w
            why[k] = why.get(k, 0) + 1
        else:
            rows.append(r)
    return rows, why


def analyze_core(sb, events, passes_by_st, cm, out_dir, log):
    np = _np()
    pe = PassEval(sb, cm)
    rows, why = compute_rows(sb, events, passes_by_st, pe, PRIMARY, cm, third=True)
    power = power_block(rows)
    power["excluded_reasons"] = why
    if out_dir:
        jdump(power, os.path.join(out_dir, "p5_power.json"))       # V16：先写
    blk, boots = ratio_block(rows, "y", "x", "R_sfc 主（A70、[0,12) h、参照 [−72,0) h）")
    verdict = judge(blk)
    # 第三方（V12）：同一事件子集同一抽样
    sub = [r for r in rows if all(isnum(r.get(k)) for k in ("obs1", "rim1"))]
    third = {}
    if len(sub) >= 3:
        tb, (bs_sub, _bt) = ratio_block(sub, "y", "x", "R_sfc（第三方子集）", pairs_extra={
            "R1m": (lambda r: r["obs1"], lambda r: r["rim1"])})
        with np.errstate(divide="ignore", invalid="ignore"):
            qarr = bs_sub["R"] / bs_sub["R1m"]
        r1 = sum(r["obs1"] for r in sub) / sum(r["rim1"] for r in sub)
        third = {"n": len(sub), "R_sfc_subset": tb.get("R"), "R_sfc_subset_ci95": tb.get("ci95_station_season"),
                 "R_1m_matched": rnd(r1), "R_1m_matched_ci95": pm.ci95(bs_sub["R1m"])[0],
                 "Q_Rsfc_over_R1m": rnd(tb["R"] / r1) if tb.get("R") is not None and r1 else None,
                 "Q_ci95": pm.ci95(qarr)[0],
                 "reading": "R1 下 Q≈1；R2 下 Q≈1/R_1m（描述，不进判读）"}
        sub5 = [r for r in rows if all(isnum(r.get(k)) for k in ("obs05", "rim05"))]
        if len(sub5) >= 3:
            third["R_05m_matched"] = rnd(sum(r["obs05"] for r in sub5) / sum(r["rim05"] for r in sub5))
            third["n_05m"] = len(sub5)
        third["sat_over_obs1m_direct"] = rnd(sum(r["y"] for r in sub) / sum(r["obs1"] for r in sub)) \
            if sum(r["obs1"] for r in sub) else None
    sc = [r for r in rows if isnum(r.get("x_pix0")) and r["x_pix0"] != 0]
    desc = {"scale_foot_over_pix0": rnd(sum(r["x"] for r in sc) / sum(r["x_pix0"] for r in sc)) if sc else None,
            "per_station": {st: {"n": len(rr), "R": rnd(sum(r["y"] for r in rr) / sum(r["x"] for r in rr))
                                 if sum(r["x"] for r in rr) else None}
                            for st, rr in _group(rows).items()},
            "dt_post_h_median": rnd(float(np.median([r["dt_post_h"] for r in rows]))) if rows else None,
            "n_ref_median": rnd(float(np.median([r["n_ref"] for r in rows]))) if rows else None,
            "imerg_center_post_median": rnd(float(np.nanmedian([r["imerg_c_post"] for r in rows]))) if rows else None,
            "cmorph_foot_rain_post_median": rnd(float(np.nanmedian([r["rain_foot_cmorph_post"] for r in rows])))
            if rows else None}
    sens = {}
    for code, label, kw in SENS:
        rr, ww = compute_rows(sb, events, passes_by_st, pe, variant(**kw), cm)
        b, _ = ratio_block(rr, "y", "x", label)
        j = judge(b)
        sens[code] = {"label": label, "n": b["n"], "clusters": b.get("clusters_station_season"), "R": b.get("R"),
                      "ci95": b.get("ci95_station_season"), "evaluable": b.get("evaluable"),
                      "category": j["category"], "same_as_primary": j["category"] == verdict["category"],
                      "excluded_reasons": ww}
    return rows, power, blk, verdict, third, desc, sens


def _group(rows):
    g = {}
    for r in rows:
        g.setdefault(r["station"], []).append(r)
    return g


CSV_FIELDS = ["station", "season", "onset_utc", "post_gid", "post_utc", "dt_post_h", "n_ref", "node", "n_looks_post",
              "sat_post", "sat_ref", "y", "F_post", "F_ref", "S0", "x", "x_pix0", "sig", "imerg_c_post",
              "rain_foot_cmorph_post", "dropped_w_post", "obs1", "obs05", "rim1", "rim05", "s1_depth"]


def run_analyze(args, out_dir, log, sb, events, info):
    cm = CmCache(args.cache)
    passes = {st: load_passes(args.cache, st) for st in {e["_st"] for e in events}}
    rows, power, blk, verdict, third, desc, sens = analyze_core(sb, events, passes, cm, out_dir, log)
    with open(os.path.join(out_dir, "p5_events.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (rnd(r[k], 6) if isinstance(r.get(k), float) else r.get(k)) for k in CSV_FIELDS})
    here = os.path.dirname(os.path.abspath(__file__))
    shas = {fn: sha256_file(os.path.join(here, fn)) for fn in
            ("p5_sss_sat.py", "p1_events.py", "p1b_extend.py", "p2_rim_test.py", "p4_mech.py")}
    S = {"script": "p5_sss_sat.py", "version": VERSION, "run_utc": iso(time.time()),
         "code_sha256": shas, "events": info, "passes_loaded": {st: len(v) for st, v in passes.items()},
         "power_file": "p5_power.json（先于本文件写盘）", "primary": blk, "verdict": verdict,
         "primary_excluded_reasons": power.get("excluded_reasons"),
         "third_party_mooring": third, "descriptive": desc, "sensitivity": sens,
         "nature": "探索性；对卫星 SSS 盲"}
    jdump(S, os.path.join(out_dir, "p5_summary.json"))
    m = f"[P5] R_sfc={blk.get('R')} CI={blk.get('ci95_station_season')} n={blk['n']} 判读={verdict['category']}"
    print(m, flush=True)
    log.log(m)
    return 0


# ======================================================================== 自测（合成，无网络）
def selftest(out_dir=None):
    np = _np()
    res = []

    def chk(name, ok, detail=""):
        res.append({"name": name, "ok": bool(ok), "detail": str(detail)[:300]})

    rng = np.random.default_rng(7)
    # 1 vec_rim0 与 p2.rim_factor(z=0) 逐值一致
    P = np.where(rng.random((49, 40)) < 0.6, 0.0, rng.gamma(0.8, 6.0, (49, 40)))
    P[:, 0] = 0.0
    P[10, 1] = 250.0                                    # d0 表外 → NaN
    U = rng.uniform(0.3, 14.0, 49)
    Fv = vec_rim0(P, U)
    Fr = np.array([p2.rim_factor(P[:, n], U, 0.0)[0] for n in range(P.shape[1])])
    fin = np.isfinite(Fr)
    rel = np.max(np.abs(Fv[fin] - Fr[fin]) / np.abs(Fr[fin]))
    chk("V9 vec_rim0 == p2.rim_factor(z=0)", rel <= 1e-12 and np.array_equal(np.isfinite(Fv), fin) and Fv[0] == 1.0,
        f"max rel {rel:.2e}; NaN 模式一致={np.array_equal(np.isfinite(Fv), fin)}")
    # 2 核：半功率宽 40 km、截断、归一
    grid = {"lat": [0.125 + 0.009 * (k - 90) for k in range(181)], "lon": [165.125 + 0.009 * (k - 90) for k in range(181)]}
    w = cell_kernel(grid, 0.125, 165.125)
    lat = np.asarray(grid["lat"])[:, None]
    lon = np.asarray(grid["lon"])[None, :]
    d = np.sqrt(((lat - 0.125) * KM_PER_DEG) ** 2 + ((lon - 165.125) * KM_PER_DEG * math.cos(math.radians(0.125))) ** 2)
    c = w[90, 90]
    near20 = w[np.abs(d - 20.0) < 0.5]
    ratio = float(np.median(near20 / c)) if near20.size else NAN
    chk("V8 核归一＋FWHM＋截断", abs(w.sum() - 1) < 1e-12 and abs(ratio - 0.5) < 0.06 and np.all(w[d > TRUNC_KM + 1e-9] == 0),
        f"sum={w.sum():.12f} w(20km)/w(0)={ratio:.3f}")
    # 3 格号
    chk("V3 格号", cell_yx(0.0, 165.0) == (360, 660) and cell_yx(-8.0, 165.0) == (328, 660)
        and cell_yx(50.1, -144.9) == (560, 860) and cell_center(360, 660) == (0.125, 165.125),
        f"{cell_yx(0.0, 165.0)} {cell_yx(50.1, -144.9)}")
    # 3b 行定位（x 为逐轨网格，lon 随 x 每格 −0.25°）
    la_row = np.full((1560, 2), np.nan)
    lo_row = np.full((1560, 2), np.nan)
    for i in range(400, 441):
        la_row[i, :] = 0.12
        lo_row[i, :] = 170.125 - 0.25 * (i - 400)
    x1, d1, n1 = locate_in_row(la_row, lo_row, *station_cell(0.0, 165.0))
    x2, d2, n2 = locate_in_row(la_row, lo_row, *station_cell(0.0, 150.0))
    chk("V3 行定位", x1 == 420 and d1 is not None and d1 < 0.01 and n1 == 6 and x2 is None and station_cell(0.0, 165.0) == (0.125, 165.125),
        f"x1={x1} d1={d1} n1={n1} x2={x2}")
    # 4 QC：位 15 保留、位 10 剔除、位 0 无效、风>20 无效、范围
    sss = np.full((3, 3), 34.0)
    iqc = np.zeros((3, 3))
    iqc[0, 0] = 1 << 15
    iqc[0, 1] = 1 << 10
    iqc[0, 2] = 1
    iqc[1, 0] = (1 << 12) | (1 << 15)
    wind = np.full((3, 3), 8.0)
    wind[2, 2] = 21.0
    sss[2, 0] = 1.0
    ok = qc_valid(sss, iqc, wind, QC_EXCLUDE)
    ok4 = qc_valid(sss, iqc, wind, (5, 6, 7, 8, 9, 14))
    chk("V5 QC 位", ok[0, 0] and not ok[0, 1] and not ok[0, 2] and ok[1, 0] and not ok[2, 2] and not ok[2, 0] and ok4[0, 1],
        ok.astype(int).tolist())
    # 5 认证主机与 redact
    chk("V22 仅对 Earthdata 主机附 token", auth_host("opendap.earthdata.nasa.gov") and auth_host("x.earthdatacloud.nasa.gov")
        and not auth_host("s3.us-west-2.amazonaws.com") and not auth_host("cmr.earthdata.nasa.gov")
        and "Authorization" not in build_headers("www.ncei.noaa.gov", "tok") and
        build_headers("opendap.earthdata.nasa.gov", "tok")["Authorization"] == "Bearer tok")
    _SECRETS.append("SECRETTOKEN123")
    chk("V22 redact", "SECRETTOKEN123" not in redact("err SECRETTOKEN123 x"))
    _SECRETS.remove("SECRETTOKEN123")
    # 6 nc4 内存解析与维序核对
    try:
        import netCDF4
        nc = netCDF4.Dataset("t.nc", "w", memory=4096)
        for dname, n in (("ydim_grid", 3), ("xdim_grid", 3), ("look", 2)):
            nc.createDimension(dname, n)
        for vname in VARS3:
            vv = nc.createVariable(vname, "f8", DIMS3, fill_value=-9999.0)
            vv[:] = np.arange(18, dtype=float).reshape(3, 3, 2)
        for vname in VARS2:
            vv = nc.createVariable(vname, "f4", DIMS3[:2], fill_value=-9999.0)
            vv[:] = np.arange(9, dtype=float).reshape(3, 3)
        nc["sss_smap_40km"][0, 0, 0] = -9999.0
        buf = nc.close()
        sub = parse_nc(buf, VARS3 + VARS2)
        check_dims(sub)
        chk("V3 nc4 内存解析", np.isnan(sub["sss_smap_40km"][0, 0, 0]) and sub["time"][2, 2, 1] == 17.0
            and sub["rain"].shape == (3, 3))
    except Exception as e:
        chk("V3 nc4 内存解析", False, f"{type(e).__name__}: {e}")
    # 7 风：与 P2 event_forcing 同一换算（半步重复小时值）
    stn = {"name": "X", "h0": 1000, "ser": {"wind": np.r_[np.full(50, 5.0), np.nan, np.nan, np.full(60, 7.0)]}}
    Uh = wind_halfsteps(stn, 2 * 1040, 2 * 1060 + 1)
    ok_w = Uh is not None and abs(Uh[0] - 5.0 * p2.WIND_FACTOR) < 1e-12 and Uh[1] == Uh[0] and len(Uh) == 42
    stn2 = {"name": "X", "h0": 1000, "ser": {"wind": np.r_[np.full(50, 5.0), np.full(8, np.nan), np.full(60, 7.0)]}}
    chk("V10 风半步", ok_w and wind_halfsteps(stn2, 2 * 1045, 2 * 1060) is None)
    # 8 判读映射
    def fake(lo, hi, st=(None, None), loo=None, ev=True):
        return {"evaluable": ev, "ci95_station_season": [lo, hi], "ci95_station": list(st),
                "jackknife_station": {"evaluable": loo is not None, "loo_range": loo}}
    j1 = judge(fake(0.1, 0.5, (0.15, 0.45), [0.2, 0.4]))
    j2 = judge(fake(0.7, 1.4, (0.8, 1.2), [0.9, 1.1]))
    j3 = judge(fake(0.4, 0.9))
    j4 = judge(fake(-0.8, -0.1, (-0.7, -0.2), [-0.5, -0.3]))
    j5 = judge(fake(0.1, 0.5, ev=False))
    j6 = judge(fake(0.1, 0.5, (0.1, 0.7), [0.2, 0.4]))
    chk("V15 判读映射", j1["category"].startswith("支持 R1") and j1["robust"] == "稳健" and j2["category"].startswith("支持 R2")
        and j3["category"] == "不能区分" and j4["category"].startswith("卫星未见淡化") and j5["category"] == "不可评"
        and j6["robust"] == "依赖簇定义", [j1, j2, j3, j4, j5, j6])
    # 9 端到端：合成缓存 → 已知 R_true 精确复原；窗选择
    try:
        e2e = _e2e(np, rng)
        chk("E2E R_sfc 复原（无噪声）", abs(e2e["R"] - 0.3) < 1e-9, e2e)
        chk("E2E 窗选择（最早雨后、参照 [−72,0)）", e2e["post_dt"] == [3.0] * e2e["n"] and e2e["n_ref"] == [2] * e2e["n"], e2e)
        chk("E2E 可评/判读与功效先写", e2e["evaluable"] and e2e["power_first"], e2e)
        chk("E2E 系泊第三方可算", e2e["third_ok"], e2e)
    except Exception as e:
        chk("E2E", False, f"{type(e).__name__}: {e} {traceback.format_exc()[-600:]}")
    n_ok = sum(r["ok"] for r in res)
    out = {"version": VERSION, "passed": n_ok, "total": len(res), "all_ok": n_ok == len(res), "checks": res}
    if out_dir:
        jdump(out, os.path.join(out_dir, "p5_selftest.json"))
    for r in res:
        print(("  ok  " if r["ok"] else "  FAIL") + f" {r['name']}" + ("" if r["ok"] else f" :: {r['detail']}"), flush=True)
    print(f"[selftest] {n_ok}/{len(res)}", flush=True)
    return out


def _e2e(np, rng):
    """合成 40 个事件（2 站×4 季×5）：每事件参照过境 −50、−10 h，雨后 +3、+8 h，另有 −80、+30 h（窗外）。
    卫星值＝34＋R_true·S0·(F̄−1)，全 9 格同值 ⇒ R_sfc 精确等于 R_true。"""
    td = tempfile.mkdtemp(prefix="p5e2e_")
    try:
        R_true = 0.3
        sts = {"A": (0.0, 165.0), "B": (15.0, 90.0)}
        sb, events = {}, []
        h_base = 400000
        n_h = 24 * 400
        for k, (st, (la, lo)) in enumerate(sts.items()):
            ser = {kk: np.full(n_h, 35.0) for kk in ("s1", "sss05", "s5", "rain", "pco2", "sst")}
            ser["wind"] = np.full(n_h, 6.0)
            sb[st] = {"name": st, "h0": h_base, "n": n_h, "ser": ser, "_coord": (la, lo), "lat": la, "lon": lo}
        seasons = [("DJF", 0), ("MAM", 90), ("JJA", 180), ("SON", 270)]
        for st in sts:
            for s, day0 in seasons:
                for m in range(5):
                    H = h_base + (day0 + 12 * m + 5) * 24 + 7
                    events.append({"_st": st, "hour": H, "i": H - h_base, "season": s, "onset_utc": p1.hour_to_iso(H),
                                   "s1_depth_m": 1.0})
        # CMORPH 窗口：用真 WinExtractor 的索引逻辑造网格
        cache = os.path.join(td, "cache")
        clat = np.arange(-59.963, 60.0, 0.072771)
        clon = np.arange(0.0364, 360.0, 0.072756)

        class _DS:
            def __init__(self):
                self.v = {"lat": clat, "lon": clon}

            def __getitem__(self, k):
                return self.v[k]
        wx = WinExtractor(sts, cache)
        wx._index(_DS())
        hours = set()
        for e in events:
            H = e["hour"]
            hours.update(range(H - 80 - 30, H + 40))
        for st in sts:
            rows, cols = wx.idx[st]
            ny, nx = len(rows), len(cols)
            d = os.path.join(cache, "cmorph", st)
            ev_h = {e["hour"] for e in events if e["_st"] == st}
            for h in sorted(hours):
                a = np.zeros((2, ny, nx), np.float32)
                if any(0 <= h - H < 4 for H in ev_h):
                    a[:] = rng.gamma(1.0, 8.0, (2, ny, nx)).astype(np.float32)
                np.save(os.path.join(d, f"{h}.npy"), a)
        cm = CmCache(cache)
        pe = PassEval(sb, cm)
        passes = {}
        for st, (la, lo) in sts.items():
            y, x = cell_yx(la, lo)
            lst = []
            for e in [e for e in events if e["_st"] == st]:
                t0 = e["hour"] * 3600.0
                for dh in (-80.0, -50.0, -10.0, 3.0, 8.0, 30.0):
                    t = t0 + dh * 3600.0
                    clat, clon = station_cell(la, lo)
                    cla = np.array([[[clat + 0.25 * (i - 1)] * 2 for _j in range(3)] for i in range(3)])
                    clo = np.array([[[(clon - 0.25 * (j - 1)) % 360.0] * 2 for j in range(3)] for _i in range(3)])
                    sub = {"sss_smap_40km": np.full((3, 3, 2), 34.0), "sss_smap_40km_unc": np.full((3, 3, 2), 0.7),
                           "time": np.full((3, 3, 2), t - EPOCH2000), "iqc_flag": np.zeros((3, 3, 2)),
                           "cellat": cla, "cellon": clo,
                           "rain": np.zeros((3, 3)), "sss_ref": np.full((3, 3), 34.0), "winspd": np.full((3, 3), 6.0),
                           "surtep": np.full((3, 3), 300.0)}
                    gid = f"G{st}{int(t)}"
                    r, why = pe.evaluate(st, gid, sub, PRIMARY)
                    if r is None:
                        raise RuntimeError(f"合成过境无效：{why}")
                    sub["sss_smap_40km"][:] = 34.0 + R_true * 34.0 * (r["F"] - 1.0)
                    lst.append((t, t, gid, sub))
            lst.sort(key=lambda p: p[0])
            passes[st] = lst
        pe2 = PassEval(sb, cm)
        rows, why = compute_rows(sb, events, passes, pe2, PRIMARY, cm, third=True)
        blk, _ = ratio_block(rows, "y", "x", "e2e")
        power = power_block(rows)
        third_ok = all(isnum(r.get("rim1")) and isnum(r.get("obs1")) for r in rows)
        return {"R": blk.get("R"), "n": len(rows), "evaluable": blk.get("evaluable"), "why": why,
                "post_dt": [r["dt_post_h"] for r in rows], "n_ref": [r["n_ref"] for r in rows],
                "power_first": "expected_se_R_formal" in power, "third_ok": third_ok}
    finally:
        shutil.rmtree(td, ignore_errors=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="P5 卫星海表盐度共址检验")
    g = ap.add_mutually_exclusive_group(required=True)
    for m in ("selftest", "plan", "fetch", "analyze", "full"):
        g.add_argument(f"--{m}", action="store_true")
    ap.add_argument("--out")
    ap.add_argument("--cache", default=CACHE_DEFAULT)
    ap.add_argument("--p1-events", default=p2.P1_EVENTS_DEFAULT)
    ap.add_argument("--p1b-dir", default=p2.P1B_DIR_DEFAULT)
    ap.add_argument("--token-file", default=TOKEN_FILE_DEFAULT)
    ap.add_argument("--max-conc", type=int, default=MAX_CONC)
    ap.add_argument("--raw-root", default=rs.RAW_ROOT_DEFAULT)
    ap.add_argument("--rate-mbps", type=float, default=2.0)
    ap.add_argument("--raw-tag", default="p5-sss")
    args = ap.parse_args(argv)
    args.max_conc = max(1, min(MAX_CONC, args.max_conc))
    out_dir = args.out or os.environ.get("REPRO_OUTPUT_DIR")
    if args.selftest:
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        return 0 if selftest(out_dir)["all_ok"] else 4
    if not out_dir:
        print("需要 --out 或环境变量 REPRO_OUTPUT_DIR", file=sys.stderr)
        return 3
    os.makedirs(out_dir, exist_ok=True)
    log = p1.Log(os.path.join(out_dir, "p5_log.txt"))
    mode = next(m for m in ("plan", "fetch", "analyze", "full") if getattr(args, m))
    log.log(f"=== start {VERSION} mode={mode} out={out_dir} cache={args.cache}", echo=True)
    try:
        if not selftest(out_dir)["all_ok"]:
            log.log("自测不过，退出 4", echo=True)
            return 4
        if mode == "plan":
            rc = run_plan(args, out_dir, log)
        else:
            sb, events, info = build(args, out_dir, log)
            if mode in ("fetch", "full"):
                run_fetch(args, out_dir, log, sb, events)
            rc = run_analyze(args, out_dir, log, sb, events, info) if mode in ("analyze", "full") else 0
    except AuthError as e:
        log.log(f"FATAL 认证失败：{redact(e)}", echo=True)
        return 5
    except (DataError, p2.DataSourceError, p1.FetchError) as e:
        log.log(f"FATAL 数据源故障（重跑即续传）：{redact(e)}", echo=True)
        return 2
    except Exception:
        log.log("FATAL 未预期异常：\n" + redact(traceback.format_exc()), echo=True)
        return 3
    log.log(f"=== done rc={rc}", echo=True)
    log.close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
