#!/usr/bin/env python3
"""p7a_jpl_smap.py — P7a：换独立处理链的 SMAP 产品（JPL SMAP L2B CAP V5.0）重做 P5 的卫星海表淡化共址检验，
检验 P5（RSS L2C V6）「支持 R2」是否来自 RSS 反演在降雨下的伪淡化（矛盾 A 判据①）。

性质：探索性；方法与判读规则在运行前写定；P5 的 RSS 结果已看过；对 JPL 盐度值盲。判读在 judge_overall()。

复用（只 import，不改）：p5_sss_sat（事件重建 build、候选轨 candidates、逐事件匹配 event_eval／compute_rows、比值与判读
  ratio_block／judge、功效 power_block、RIM 向量化 vec_rim0、风 wind_halfsteps、丢弃规则 weighted_F、HTTP 与认证、
  CMORPH 窗口 WinExtractor／CmCache）；raw_store（原件落档、令牌桶、manifest）；raw_fetch_list（心跳接替 Relay／wait_turn）。
  与 P5 事件集逐条一致：同一 p5.build → 317 个 SMAP 期事件；自测 V0 用 P5 缓存重算 P5 主结果，逐事件核对 P5 的
  p5_events.csv（y、x 相对差 ≤1e−5）与 p5_summary.json 的 R，不一致即退出 3。

用法（产物写 --out 或 $REPRO_OUTPUT_DIR）：
  python p7a_jpl_smap.py --selftest     合成数据自测（无网络）
  python p7a_jpl_smap.py --plan         自测 → 事件与 P5 一致性 → 候选轨（P5 同一套）→ JPL 颗粒按轨道号对应 →
                                        认证与格式探测（.dmr、一轨 lat/lon/row_time、各站 3×3 的 lat/lon/row_time；不取盐度）
  python p7a_jpl_smap.py --full         自测 → 下载（JPL 子集原样落 bulk raw；缺的 CMORPH 小时）→ V0 → 分析
  python p7a_jpl_smap.py --fetch / --analyze   只下载／只读缓存分析
  --cache DIR（默认 <fast_root>/p7a-cache）  --p5-cache DIR（只读）
  --p5-results DIR（P5 正式产物，只读；V0 用）  --raw-root DIR  --rate-mbps R（≤2）  --raw-tag S  --max-conc N（≤2）
依赖：numpy、scipy、netCDF4＋同目录 p1_events／p1b_extend／p2_rim_test／p4_mech／p5_sss_sat／raw_store／raw_fetch_list。
产物：p7a_selftest.json、p7a_plan.json、p7a_fetch.json、p7a_v0.json、p7a_power.json（先写）、p7a_summary.json、
  p7a_events_jpl.csv、p7a_pairs.csv、p7a_log.txt。
退出码：0 跑完（无论判读）；2 数据源故障（重跑即续传）；3 其他异常（含 V0 不符、维序不符）；4 自测不过；5 认证失败；6 等心跳超时。

实现选择（Y 条）：
  Y1 产品：PO.DAAC `SMAP_JPL_L2B_SSS_CAP_V5`（CMR C2208420167-POCLOUD），每轨一个 HDF5；轨道号（REVNO）与 RSS L2C 的 r 号一一对应，
     候选轨＝P5 的候选轨（p5.candidates，读 P5 的 CMR 缓存副本）按轨道号换成 JPL 颗粒；JPL 颗粒清单由 CMR 时间查询（公开）得到。
  Y2 取数：OPeNDAP DAP4（Earthdata token，同 P5 V22），每轨先取整轨 lat、lon、row_time（定位用），再按站取 3×3 SWC 子集
     （smap_sss、smap_sss_uncertainty、quality_flag、lat、lon、anc_sss、anc_spd、smap_spd、anc_sst、row_time[j−1..j+1]）。
     数组维序 (ncti=76, nati=1624)（用户手册 §6.2、UMM-Var；plan 核 .dmr），不符即退出 3。返回字节原样落 raw/smap_jpl/。
  Y3 定位：站点坐标到各 SWC（lat/lon）的等距近似距离最小者，≤20 km 且 i∈[1,74]、j∈[1,1622] 为覆盖，否则 .nocover；
     子集取回后核中心格 lat/lon 与整轨数组同值。
  Y4 过境时刻＝中心行 row_time＋1420070400 s；JPL 前后视向已合并成单一反演（无 look 维），一次过境一个值。
  Y5 QC（主）：quality_flag 位 0（SSS usable；手册 §5.1 推荐）、5（辅助风 >20 m/s）、7（陆地）、8（海冰）任一置位或为填充值即无效；
     smap_sss ∉[2,42] 无效。JPL 无降雨标志位（手册 §6.2.24 只有位 0–9），所以没有「保留雨标志」这一项。
  Y6 卫星值：A70 类比＝中心 SWC 有效且 3×3 有效 ≥5 时取有效格均值，否则无效；A40 类比＝只中心 SWC（敏感性）。
  Y7 窗与差分：与 P5 V7／V11 相同（直接调用 p5.event_eval）：雨后 [0,12) h 最早有效过境，参照 [−72,0) h 全部有效过境，
     ΔSSS＝雨后−参照均值；ΔS_RIM,foot＝S0·[F̄(雨后)−mean F̄(参照)]，S0＝雨后中心 anc_sss（HYCOM），缺则参照卫星均值。
  Y8 足印核：每个有效 SWC 以其 lat/lon 为中心的圆高斯×cos φ，主 FWHM 40 km、截断 60 km（与 P5 V8 同核，便于与 RSS 同尺度比较）；
     A70 核＝有效格核均值。敏感性 J4 与分辨率对照：FWHM 60 km（JPL 标称分辨率）、截断 90 km（A70）。模型在 CMORPH 站点窗口全窗计算
     （P5 为站点 100 km 超集；两者在核内逐值相同）。
  Y9 RIM、风、丢弃规则：p5.vec_rim0、p5.wind_halfsteps、p5.weighted_F（像元 NaN 权重 >5% → 无效），q＝⌊t/1800⌋。
  Y10 IMERG：JPL 文件没有雨量；雨中／不雨分组用同一轨道号的 RSS L2C 站点格 `rain`（IMERG 按 40 km 足印平均，P5 缓存，只读），
     >0.1 mm/h 为「正在下雨」，缺值不进分组（计数）。次分组：JPL 足印核加权的 CMORPH 过境半步雨量 >0.1 mm/h。
  Y11 CMORPH：先读 P5 窗口缓存（只读），缺的小时从 bulk raw 原件抽窗口（没有原件才按 2 MB/s 下载进 raw，原件不删）到本缓存；
     P5 记为缺测的小时一律按缺测处理（不重取）。网格核对：本缓存 grid.json 须与 P5 的逐值相同。
  Y12 限速：OPeNDAP 与 CMORPH 共用一个令牌桶，2 MB/s；OPeNDAP 响应到手后按字节数分块扣令牌（长期平均＝限速）。
     下载段先 raw_fetch_list.wait_turn（等别人的 active 心跳结束），再 Relay 接替写 `raw/.p5_fetch_heartbeat`（补档任务据此降到
     2 MB/s），结束写 done。没有要下的就不写心跳。
  Y13 V0：用 P5 缓存与 p5.compute_rows(PRIMARY) 重算 RSS 主结果，逐事件核 p5_events.csv 的 (站, onset) 集合与 y、x
     （|差|≤1e−5·max(1,|值|)），并核 R 与 p5_summary.json（5 位小数）；RSS 端一律只用 P5 缓存（不看本缓存新抽的小时）。
  Y14 逐对（对比 d）：同一事件、同一轨道号：雨后取两产品都有效且各自 Δt∈[0,12) h 的最早轨道；参照取两产品都有效且 Δt∈[−72,0) h
     的全部轨道（≥1）；各自按 P5 公式算 y、x；d＝y_RSS−y_JPL；D＝Σd／Σx_RSS（RSS 比值中未被 JPL 复现的部分，以 P5 分母为单位）；
     分组（IMERG，Y10）：雨中／不雨组的平均 d（psu）之差 Δd̄＝d̄_rain−d̄_dry（主；不依赖分组分母），D_rain、D_dry 只描述；
     分辨率对照：同一对子集上 x_J60（JPL 60 km 核）→ ΔR60＝Σy_RSS/Σx_RSS − Σy_JPL/Σx_J60；
     站×季整簇 bootstrap 同一抽样（B=10000，seed 20260926）。
  Y15 判读：JPL 主结果按 P5 同一判读（p5.judge，阈值 0.65、可评条件同）；总判读 judge_overall()。
  Y16 功效：p7a_power.json 在任何 Σy 之前写盘（模型侧 x、JPL 形式不确定度、n、簇；逐对子集与分组计数）。
     JPL 的 smap_sss_uncertainty 是似然 FWHM 估计（手册 §3.4.1），与 RSS 形式不确定度定义不同，只作信息。
  Y17 敏感性 J1–J7（只报 n、簇、点估计、CI、类别与是否同类，不改判读）。

Change Log：
  2026-09-27 初版。
"""

import argparse
import concurrent.futures as cf
import csv
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
import threading
import time
import traceback
import urllib.parse

import p1_events as p1
import p2_rim_test as p2
import p4_mech as pm
import p5_sss_sat as p5
import raw_store as rs
import raw_fetch_list as rfl

VERSION = "p7a-2026-09-27a"
NAN = float("nan")
isnum, rnd, iso, jdump, safe = p5.isnum, p5.rnd, p5.iso, p5.jdump, p5.safe

# ---- 数据与格式（Y1–Y4） ----
JCOLL = "C2208420167-POCLOUD"                 # SMAP_JPL_L2B_SSS_CAP_V5
JOPENDAP = f"https://opendap.earthdata.nasa.gov/collections/{JCOLL}/granules/"
CMR_CSV = "https://cmr.earthdata.nasa.gov/search/granules.csv"
EPOCH2015 = 1420070400
NCT, NAT = 76, 1624
V2D = ("smap_sss", "smap_sss_uncertainty", "quality_flag", "lat", "lon", "anc_sss", "anc_spd", "smap_spd", "anc_sst")
LOCATE_MAX_KM = 20.0
QF_FILL = 65535
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
CACHE_DEFAULT = _rp.path("p7a_cache")  # [repro] 路径来自集中配置
P5_CACHE_DEFAULT = p5.CACHE_DEFAULT
P5_RESULTS_DEFAULT = _rp.upstream("p5_dir")  # [repro] 读 p5-sss 阶段产物
RAW_SUB = "smap_jpl"
MAX_CONC = 2
RATE_MAX = 2.0

# ---- 事先写定的参数 ----
RAIN_THR = 0.1
JQC_BITS = (0, 5, 7, 8)
PRIMARY_J = dict(p5.PRIMARY, name="primary", agg="A70", fwhm=40.0, trunc=60.0, qc_bits=JQC_BITS)
SENS_J = [
    ("J1", "A40 类比（只中心 SWC）", {"agg": "A40"}),
    ("J2", "雨后窗 [0,6) h", {"post": (0.0, 6.0)}),
    ("J3", "雨后窗 [0,24) h", {"post": (0.0, 24.0)}),
    ("J4", "核 FWHM 60 km、截断 90 km（JPL 标称分辨率）", {"fwhm": 60.0, "trunc": 90.0}),
    ("J5", "剔除 KEO、Papa", {"excl": ("KEO", "Papa")}),
    ("J6", "参照只用同交点", {"same_node": True}),
    ("J7", "QC 不用位 0（只剔位 5、7、8）", {"qc_bits": (5, 7, 8)}),
]
PAIR_MIN_N, PAIR_MIN_G = 30, 8
J60 = dict(PRIMARY_J, fwhm=60.0, trunc=90.0)          # 分辨率对照（Y14）
GROUP_MIN_N, GROUP_MIN_G = 10, 5
GROUP_DEN_BAD_MAX = 0.05


def _np():
    import numpy
    return numpy


def rev_of(gid):
    """RSS『RSS_SMAP_SSS_L2C_r21988_…』或 JPL『SMAP_L2B_SSS_21988_…』→ 21988。"""
    if gid.startswith("RSS_"):
        return int(gid.split("_r", 1)[1].split("_", 1)[0])
    return int(gid.split("_")[3])


# ======================================================================== 限速与原件（Y2、Y12）
def throttle(bucket, n):
    """按字节数分块扣令牌（单次不超过桶容量，避免 n>burst 时死等）。"""
    while n > 0:
        k = min(n, rs.CHUNK)
        bucket.consume(k)
        n -= k


class Net:
    """下载段：心跳接替＋共用令牌桶＋原件落档。"""

    def __init__(self, raw_root, rate_mbps, tag, log, wait_max_s=6 * 3600, settle_s=20.0):
        self.root = rs.ensure_root(raw_root)
        self.man = rs.Manifest(self.root)
        self.tag = tag
        self.log = log
        self.rate = rate_mbps * rs.MB
        self.bucket = rs.TokenBucket(self.rate)
        self.bucket.set_rate(0)
        self.relay = None
        self.wait_max_s, self.settle_s = wait_max_s, settle_s
        self.bytes = 0
        self.n_req = 0
        self.reused = 0
        self.lock = threading.Lock()

    def start(self):
        ok, waited = rfl.wait_turn(self.root, self.wait_max_s, self.log)
        if not ok:
            raise WaitTimeout(f"等心跳超过 {self.wait_max_s:.0f} s")
        self.relay = rfl.Relay(self.root, self.bucket, self.rate, self.log, settle=self.settle_s)
        self.relay.take()
        self.relay.start()
        return waited

    def stop(self):
        if self.relay is not None:
            self.relay.release()
            self.relay = None

    def get(self, url, token, rel, dataset="smap_jpl_opendap", timeout=180, reuse=False):
        """reuse=True：raw 档已有同名原件就直接读它（不重下）。"""
        dst = os.path.join(self.root, rel)
        if reuse and os.path.exists(dst):
            with open(dst, "rb") as f:
                buf = f.read()
            with self.lock:
                self.reused += 1
            return buf
        buf = p5.get_retry(url, token, timeout)
        throttle(self.bucket, len(buf))
        with self.lock:
            self.bytes += len(buf)
            self.n_req += 1
        self.keep(buf, url, rel, dataset)
        return buf

    def keep(self, buf, url, rel, dataset):
        dst = os.path.join(self.root, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.exists(dst):
            rel2 = rel + f".dup{int(time.time())}"       # 已有同名原件不覆盖：另存
            dst = os.path.join(self.root, rel2)
            rel = rel2
        tmp = dst + f".{os.getpid()}.{threading.get_ident()}.part"
        with open(tmp, "wb") as f:
            f.write(buf)
        os.replace(tmp, dst)
        self.man.add(dataset=dataset, url=url, path=rel, bytes=len(buf), sha256=hashlib.sha256(buf).hexdigest(),
                     status="ok", source_job=self.tag)


class WaitTimeout(Exception):
    pass


# ======================================================================== CMR（Y1）
def jpl_granules(t0, t1, cache_path, log):
    """JPL L2B 颗粒 {rev: [gid, start, end]}（CMR 公开元数据，时间查询，翻页）。"""
    if os.path.exists(cache_path):
        d = json.load(open(cache_path))
        if d.get("span") == [t0, t1]:
            return {int(k): v for k, v in d["by_rev"].items()}
    out, sa = {}, None
    while True:
        q = urllib.parse.urlencode({"collection_concept_id": JCOLL, "temporal": f"{iso(t0)},{iso(t1)}",
                                    "page_size": 2000, "sort_key": "start_date"})
        body, hdr = p5.get_retry(f"{CMR_CSV}?{q}", None, 120, want_headers=True,
                                 extra={"CMR-Search-After": sa} if sa else None)
        rows = list(csv.reader(body.decode("utf-8").splitlines()))[1:]
        for x in rows:
            if len(x) >= 4 and x[0].startswith("SMAP_L2B_SSS_"):
                out[rev_of(x[0])] = [x[0], p5.parse_iso(x[2]), p5.parse_iso(x[3])]
        sa = hdr.get("CMR-Search-After") or hdr.get("Cmr-Search-After")
        if len(rows) < 2000 or not sa:
            break
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    jdump({"span": [t0, t1], "by_rev": {str(k): v for k, v in out.items()}}, cache_path)
    log.log(f"CMR JPL：{len(out)} 轨（{iso(t0)}–{iso(t1)}）")
    return out


def p5_candidates(sb, ebs, cache, p5_cache, log):
    """Y1：复制 P5 的 CMR 缓存（只读源）到本缓存，调用 p5.candidates 得到与 P5 相同的候选轨；再核与 P5 子集文件清单一致。"""
    src = os.path.join(p5_cache, "cmr")
    dst = os.path.join(cache, "cmr")
    os.makedirs(dst, exist_ok=True)
    for fn in os.listdir(src):
        if fn.endswith(".json") and not os.path.exists(os.path.join(dst, fn)):
            shutil.copy2(os.path.join(src, fn), os.path.join(dst, fn))
    cand, _approx = p5.candidates(sb, ebs, cache, log)
    listed = {}
    for st in cand:
        d = os.path.join(p5_cache, "smap", safe(st))
        names = set()
        if os.path.isdir(d):
            for fn in os.listdir(d):
                for ext in (".npz", ".npz.missing", ".npz.nocover"):
                    if fn.endswith(ext) and not fn.endswith(".part.npz"):
                        names.add(fn[: -len(ext)])
        listed[st] = names
    same = all(set(cand[st]) == listed[st] for st in cand)
    return cand, same, {st: [len(cand[st]), len(listed[st])] for st in cand}


# ======================================================================== OPeNDAP（Y2、Y3）
def geo_url(gid):
    return JOPENDAP + urllib.parse.quote(gid) + ".dap.nc4?dap4.ce=" + urllib.parse.quote("/lat;/lon;/row_time", safe="")


def sub_url(gid, i, j, names=V2D):
    s = f"[{i - 1}:1:{i + 1}][{j - 1}:1:{j + 1}]"
    ce = ";".join([f"/{v}{s}" for v in names] + [f"/row_time[{j - 1}:1:{j + 1}]"])
    return JOPENDAP + urllib.parse.quote(gid) + ".dap.nc4?dap4.ce=" + urllib.parse.quote(ce, safe="")


def dmr_url(gid):
    return JOPENDAP + urllib.parse.quote(gid) + ".dmr"


def dmr_dim_sizes(xml_bytes, var):
    """.dmr → 变量 var 的维名与大小（Dimension 声明表查大小）。"""
    import xml.etree.ElementTree as ET
    root = ET.fromstring(xml_bytes)
    sizes = {}
    for el in root.iter():
        if el.tag.endswith("Dimension") and "name" in el.attrib:
            sizes[el.attrib["name"].split("/")[-1]] = int(el.attrib.get("size", "-1"))
    for el in root.iter():
        if el.attrib.get("name") == var and not el.tag.endswith("Dimension"):
            out = []
            for c in el:
                if c.tag.endswith("Dim"):
                    d = c.attrib.get("name", "").split("/")[-1]
                    out.append((d, int(c.attrib["size"]) if "size" in c.attrib else sizes.get(d)))
            return out
    return None


def parse_geo(buf):
    g = p5.parse_nc(buf, ("lat", "lon", "row_time"))
    if g["lat"].shape != (NCT, NAT) or g["lon"].shape != (NCT, NAT) or g["row_time"].shape != (NAT,):
        raise ValueError(f"Y2 维序／形状不符：lat {g['lat'].shape} row_time {g['row_time'].shape}")
    return g


def locate_swath(lat, lon, row_time, la, lo):
    """Y3：返回 (i, j, 距离 km, 20 km 内相距 >50 行的另一组数) 或 (None, None, 最近距离, 0)。"""
    np = _np()
    ok = np.isfinite(lat) & np.isfinite(lon)
    if not ok.any():
        return None, None, None, 0
    dy = (lat - la) * p5.KM_PER_DEG
    dx = (((lon - lo) + 180.0) % 360.0 - 180.0) * p5.KM_PER_DEG * math.cos(math.radians(la))
    d = np.where(ok, np.sqrt(dy ** 2 + dx ** 2), np.inf)
    k = int(np.argmin(d))
    i, j = divmod(k, lat.shape[1])
    dmin = float(d[i, j])
    if dmin > LOCATE_MAX_KM:
        return None, None, dmin, 0
    near = np.argwhere(d <= LOCATE_MAX_KM)
    n_alt = int(np.sum(np.abs(near[:, 1] - j) > 50))
    if not (1 <= i <= NCT - 2 and 1 <= j <= NAT - 2) or not np.isfinite(row_time[j]):
        return None, None, dmin, n_alt
    return i, j, dmin, n_alt


def parse_sub(buf, names=V2D):
    sub = p5.parse_nc(buf, tuple(names) + ("row_time",))
    for v in names:
        if sub[v].shape != (3, 3):
            raise ValueError(f"Y2 子集形状不符：{v} {sub[v].shape}")
    if sub["row_time"].shape != (3,):
        raise ValueError(f"Y2 row_time 形状不符：{sub['row_time'].shape}")
    return {k: v for k, v in sub.items() if not k.endswith("__dims")}


def save_npz(path, sub, meta):
    np = _np()
    tmp = path + ".part.npz"
    np.savez(tmp, meta=json.dumps(meta), **sub)
    os.replace(tmp, path)


def load_npz(path):
    np = _np()
    with np.load(path, allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    d["meta"] = json.loads(str(d["meta"]))
    return d


def load_jpl_passes(cache, st):
    """[(t, t, gid, sub)]（按时刻）；只收中心行时刻有限者。"""
    np = _np()
    d = os.path.join(cache, "smap_jpl", safe(st))
    out = []
    if not os.path.isdir(d):
        return out
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".npz") or fn.endswith(".part.npz"):
            continue
        sub = load_npz(os.path.join(d, fn))
        rt = sub["row_time"][1]
        if np.isfinite(rt):
            t = float(rt) + EPOCH2015
            out.append((t, t, sub["meta"]["gid"], sub))
    out.sort(key=lambda p: p[0])
    return out


# ======================================================================== QC、核、逐过境（Y5–Y10）
def jqc_valid(sss, qf, bits):
    np = _np()
    q = np.where(np.isfinite(qf), qf, -1).astype(np.int64)
    ok = np.isfinite(sss) & (sss >= p5.SSS_RANGE[0]) & (sss <= p5.SSS_RANGE[1]) & (q >= 0) & (q != QF_FILL)
    for b in bits:
        ok &= ((q >> b) & 1) == 0
    return ok


_JK = {}


def jkernel(grid, clat, clon, fwhm, trunc):
    np = _np()
    clat, clon = round(float(clat), 4), round(float(clon) % 360.0, 4)
    key = (tuple(grid["lat"][:2]), tuple(grid["lon"][:2]), len(grid["lat"]), len(grid["lon"]), clat, clon, fwhm, trunc)
    if key in _JK:
        return _JK[key]
    lat = np.asarray(grid["lat"], float)[:, None]
    lon = np.asarray(grid["lon"], float)[None, :]
    dy = (lat - clat) * p5.KM_PER_DEG
    dx = (((lon - clon) + 180.0) % 360.0 - 180.0) * p5.KM_PER_DEG * math.cos(math.radians(clat))
    d2 = dy ** 2 + dx ** 2
    w = np.exp(-4.0 * math.log(2.0) * d2 / fwhm ** 2) * np.cos(np.radians(lat))
    w = np.where(d2 <= trunc ** 2, w, 0.0)
    s = w.sum()
    w = w / s if s > 0 else w * np.nan
    _JK[key] = w
    return w


class MultiCm(p5.CmCache):
    """CMORPH 站点窗口：按 roots 顺序查找（Y11）；缺测＝各 root missing.jsonl 的并集。"""

    def __init__(self, roots, maxn=600):
        super().__init__(roots[0], maxn)
        self.roots = [os.path.join(r, "cmorph") for r in roots]
        self.missing = {}
        for r in roots:
            self.missing.update(p5.cmorph_missing(r))

    def g(self, st):
        if st not in self.grid:
            self.grid[st] = None
            for r in self.roots:
                p = os.path.join(r, safe(st), "grid.json")
                if os.path.exists(p):
                    self.grid[st] = json.load(open(p))
                    break
        return self.grid[st]

    def has(self, st, h):
        return any(os.path.exists(os.path.join(r, safe(st), f"{h}.npy")) for r in self.roots)

    def hour(self, st, h):
        np = _np()
        k = (st, h)
        if k in self.lru:
            return self.lru[k]
        a = None
        for r in self.roots:
            p = os.path.join(r, safe(st), f"{h}.npy")
            if os.path.exists(p):
                a = np.load(p).astype(float)
                break
        self.lru[k] = a
        self.order.append(k)
        if len(self.order) > self.maxn:
            self.lru.pop(self.order.pop(0), None)
        return a


class JPassEval:
    """JPL 逐过境：卫星值（Y5、Y6）与足印模型值（Y8、Y9），带缓存。接口与 p5.PassEval.evaluate 相同，可直接喂 p5.event_eval。"""

    def __init__(self, stations, cm, rss_rain):
        self.stations = stations
        self.cm = cm
        self.rss_rain = rss_rain
        self.fq = {}
        self.memo = {}

    def pix(self, st, q):
        np = _np()
        k = (st, q)
        if k in self.fq:
            return self.fq[k]
        if self.cm.g(st) is None:
            self.fq[k] = (None, "cmorph_grid_missing")
            return self.fq[k]
        Pst = self.cm.halfsteps(st, q - 48, q)
        if Pst is None:
            self.fq[k] = (None, "cmorph_missing")
            return self.fq[k]
        U = p5.wind_halfsteps(self.stations[st], q - 48, q)
        if U is None:
            self.fq[k] = (None, "wind_gap")
            return self.fq[k]
        ny, nx = Pst.shape[1:]
        F = p5.vec_rim0(Pst.reshape(49, -1), U).reshape(ny, nx)
        self.fq[k] = ((F, Pst[-1].copy()), None)
        return self.fq[k]

    def evaluate(self, st, gid, sub, v):
        np = _np()
        key = (st, gid, v["agg"], tuple(v["qc_bits"]), v["fwhm"], v["trunc"])
        if key in self.memo:
            return self.memo[key]
        out = self._eval(st, gid, sub, v)
        self.memo[key] = out
        return out

    def _eval(self, st, gid, sub, v):
        np = _np()
        sss, la, lo = sub["smap_sss"], sub["lat"], sub["lon"]
        ok = jqc_valid(sss, sub["quality_flag"], v["qc_bits"]) & np.isfinite(la) & np.isfinite(lo)
        rt = sub["row_time"][1]
        if not np.isfinite(rt):
            return None, "no_time"
        if not ok[1, 1]:
            return None, "center_qc"
        if v["agg"] == "A70":
            if int(ok.sum()) < p5.A70_MIN_VALID:
                return None, "a70_few_valid"
            cells = [(i, j) for i in range(3) for j in range(3) if ok[i, j]]
        else:
            cells = [(1, 1)]
        sat = float(np.mean([sss[i, j] for i, j in cells]))
        unc = [sub["smap_sss_uncertainty"][i, j] for i, j in cells]
        unc = [u for u in unc if np.isfinite(u)]
        sig = float(np.mean(unc) / math.sqrt(len(cells))) if unc else NAN
        t = float(rt) + EPOCH2015
        q = int(t // 1800)
        res, why = self.pix(st, q)
        if res is None:
            return None, why
        F2d, Plast = res
        g = self.cm.g(st)
        w = np.zeros(F2d.shape)
        for i, j in cells:
            w += jkernel(g, la[i, j], lo[i, j], v["fwhm"], v["trunc"])
        w /= len(cells)
        Fbar, dropped = p5.weighted_F(F2d, w)
        if not isnum(Fbar):
            return None, "model_nan"
        m = w > 0
        rain_foot = float((Plast[m] * w[m]).sum() / w[m].sum()) if np.all(np.isfinite(Plast[m])) else NAN
        ci, cj = g["center"]
        lon_st = self.stations[st]["_coord"][1]
        lst = ((t / 3600.0) + p5.lon180(lon_st) / 15.0) % 24.0
        sref = sub["anc_sss"][1, 1]
        rev = rev_of(gid)
        return {"t": t, "sat": sat, "F": Fbar, "Fc": float(F2d[ci, cj]), "sig": sig, "n_looks": 1,
                "node": "am" if lst < 12.0 else "pm", "lst": lst,
                "sss_ref": float(sref) if np.isfinite(sref) else NAN,
                "imerg_c": self.rss_rain.get((st, rev), NAN), "rain_foot_cmorph": rain_foot, "q": q, "gid": gid,
                "rev": rev, "dropped_w": dropped, "n_cells": len(cells)}, None


def rss_rain_table(p5_cache, stations):
    """Y10：{(站, rev): RSS 站点格 IMERG 雨量}（P5 缓存 3×3 子集的中心 rain；只读）。"""
    np = _np()
    out = {}
    for st in stations:
        d = os.path.join(p5_cache, "smap", safe(st))
        if not os.path.isdir(d):
            continue
        for fn in os.listdir(d):
            if fn.endswith(".npz") and not fn.endswith(".part.npz"):
                sub = p5.load_subset(os.path.join(d, fn))
                r = sub["rain"][1, 1]
                out[(st, rev_of(sub["meta"]["gid"]))] = float(r) if np.isfinite(r) else NAN
    return out


# ======================================================================== 逐对（Y14）
def _side(post, refs):
    np = _np()
    sat_ref = float(np.mean([r["sat"] for r in refs]))
    F_ref = float(np.mean([r["F"] for r in refs]))
    s0 = post["sss_ref"] if isnum(post["sss_ref"]) else sat_ref
    return post["sat"] - sat_ref, s0 * (post["F"] - F_ref)


def pair_rows(events, rss_passes, jpl_passes, pe_r, pe_j, vr=None, vj=None):
    vr = vr or p5.PRIMARY
    vj = vj or PRIMARY_J
    plo, phi = vr["post"]
    rlo, rhi = vr["ref"]
    rows, why = [], {}
    for e in events:
        st = e["_st"]
        t0 = e["hour"] * 3600.0
        R, J, jsub = {}, {}, {}
        for tc, _t1, gid, sub in rss_passes.get(st, []):
            if rlo - 1.0 <= (tc - t0) / 3600.0 < phi + 1.0:
                r, _w = pe_r.evaluate(st, gid, sub, vr)
                if r is not None:
                    R[rev_of(gid)] = r
        for tc, _t1, gid, sub in jpl_passes.get(st, []):
            if rlo - 1.0 <= (tc - t0) / 3600.0 < phi + 1.0:
                r, _w = pe_j.evaluate(st, gid, sub, vj)
                if r is not None:
                    J[rev_of(gid)] = r
                    jsub[gid] = sub
        common = sorted(set(R) & set(J))

        def dt(x):
            return (x["t"] - t0) / 3600.0
        posts = [k for k in common if plo <= dt(R[k]) < phi and plo <= dt(J[k]) < phi]
        if not posts:
            why["no_common_post"] = why.get("no_common_post", 0) + 1
            continue
        kp = min(posts, key=lambda k: R[k]["t"])
        refs = [k for k in common if rlo <= dt(R[k]) < rhi and rlo <= dt(J[k]) < rhi]
        if not refs:
            why["no_common_ref"] = why.get("no_common_ref", 0) + 1
            continue
        yR, xR = _side(R[kp], [R[k] for k in refs])
        yJ, xJ = _side(J[kp], [J[k] for k in refs])
        xJ60 = NAN
        g60 = [pe_j.evaluate(st, J[k]["gid"], jsub[J[k]["gid"]], J60)[0] for k in [kp] + refs]
        if all(x is not None for x in g60):
            xJ60 = _side(g60[0], g60[1:])[1]
        rows.append({"station": st, "season": e["season"], "onset_utc": e["onset_utc"], "post_rev": kp,
                     "dt_post_h": round(dt(R[kp]), 3), "n_ref": len(refs), "node": R[kp]["node"],
                     "y_R": yR, "x_R": xR, "y_J": yJ, "x_J": xJ, "x_J60": xJ60, "d": yR - yJ,
                     "imerg_post": R[kp]["imerg_c"], "cmorph_foot_post_J": J[kp]["rain_foot_cmorph"]})
    return rows, why


def pair_stats(rows):
    """Y14：D、D_rain、D_dry、ΔD 与配套比值，站×季整簇 bootstrap 同一抽样。"""
    np = _np()
    rows = [r for r in rows if all(isnum(r[k]) for k in ("y_R", "x_R", "y_J", "x_J"))]
    n = len(rows)
    out = {"n": n, "stations": len({r["station"] for r in rows}),
           "clusters_station_season": len({(r["station"], r["season"]) for r in rows})}
    if n < 3:
        out.update({"evaluable": False, "reason": "事件不足"})
        return out
    keys = [(r["station"], r["season"]) for r in rows]
    a = {k: np.array([r[k] for r in rows], float) for k in ("y_R", "x_R", "y_J", "x_J", "d")}
    im = np.array([r["imerg_post"] if isnum(r["imerg_post"]) else np.nan for r in rows], float)
    Ir = (np.isfinite(im) & (im > RAIN_THR)).astype(float)
    Id = (np.isfinite(im) & (im <= RAIN_THR)).astype(float)
    one = np.ones(n)
    pairs = {"D": (a["d"], a["x_R"]), "R_RSS": (a["y_R"], a["x_R"]), "R_JPL": (a["y_J"], a["x_J"]),
             "R_JPL_xR": (a["y_J"], a["x_R"]), "den": (a["x_R"], one),
             "D_rain": (a["d"] * Ir, a["x_R"] * Ir), "D_dry": (a["d"] * Id, a["x_R"] * Id),
             "den_rain": (a["x_R"] * Ir, one), "den_dry": (a["x_R"] * Id, one),
             "md_rain": (a["d"] * Ir, Ir), "md_dry": (a["d"] * Id, Id)}
    bs, K = p2.cluster_boot(keys, pairs, p5.B, p5.SEED)
    dD = bs["D_rain"] - bs["D_dry"]
    dmd = bs["md_rain"] - bs["md_dry"]
    bad_g = float(np.mean(~np.isfinite(dmd)))
    frac = float(np.mean(bs["den"] >= 0))

    def pt(num, den):
        s = float(den.sum())
        return rnd(float(num.sum()) / s) if s else None
    nr, nd = int(Ir.sum()), int(Id.sum())
    Gr = len({k for k, w in zip(keys, Ir) if w})
    Gd = len({k for k, w in zip(keys, Id) if w})
    Dr, Dd = pt(a["d"] * Ir, a["x_R"] * Ir), pt(a["d"] * Id, a["x_R"] * Id)
    jk = pm.jackknife_station(rows, lambda rr: (sum(r["d"] for r in rr) / sum(r["x_R"] for r in rr))
                              if rr and sum(r["x_R"] for r in rr) != 0 else NAN)
    out.update({
        "sum_d": rnd(float(a["d"].sum())), "sum_x_R": rnd(float(a["x_R"].sum())), "sum_x_J": rnd(float(a["x_J"].sum())),
        "R_RSS_paired": pt(a["y_R"], a["x_R"]), "R_RSS_paired_ci95": pm.ci95(bs["R_RSS"])[0],
        "R_JPL_paired": pt(a["y_J"], a["x_J"]), "R_JPL_paired_ci95": pm.ci95(bs["R_JPL"])[0],
        "R_JPL_over_xRSS": pt(a["y_J"], a["x_R"]), "R_JPL_over_xRSS_ci95": pm.ci95(bs["R_JPL_xR"])[0],
        "D": pt(a["d"], a["x_R"]), "D_ci95": pm.ci95(bs["D"])[0], "D_jackknife_station": jk,
        "mean_d_psu": rnd(float(a["d"].mean())),
        "frac_boot_den_nonneg": rnd(frac, 4),
        "evaluable": bool(n >= PAIR_MIN_N and K >= PAIR_MIN_G and frac <= p5.DEN_NONNEG_MAX),
        "groups": {"n_rain": nr, "n_dry": nd, "n_imerg_missing": int(n - nr - nd), "clusters_rain": Gr,
                   "clusters_dry": Gd, "D_rain": Dr, "D_rain_ci95": pm.ci95(bs["D_rain"])[0],
                   "D_dry": Dd, "D_dry_ci95": pm.ci95(bs["D_dry"])[0],
                   "dD": rnd(Dr - Dd) if Dr is not None and Dd is not None else None, "dD_ci95": pm.ci95(dD)[0],
                   "mean_d_rain_psu": rnd(float((a["d"] * Ir).sum() / nr)) if nr else None,
                   "mean_d_dry_psu": rnd(float((a["d"] * Id).sum() / nd)) if nd else None,
                   "mean_d_rain_ci95": pm.ci95(bs["md_rain"])[0], "mean_d_dry_ci95": pm.ci95(bs["md_dry"])[0],
                   "dmd_psu": rnd(float((a["d"] * Ir).sum() / nr) - float((a["d"] * Id).sum() / nd)) if nr and nd else None,
                   "dmd_ci95": pm.ci95(dmd)[0],
                   "frac_boot_den_rain_nonneg": rnd(float(np.mean(bs["den_rain"] >= 0)), 4),
                   "frac_boot_den_dry_nonneg": rnd(float(np.mean(bs["den_dry"] >= 0)), 4),
                   "frac_boot_group_bad": rnd(bad_g, 4),
                   "evaluable": bool(nr >= GROUP_MIN_N and nd >= GROUP_MIN_N and Gr >= GROUP_MIN_G and Gd >= GROUP_MIN_G
                                     and bad_g <= GROUP_DEN_BAD_MAX)}})
    sub60 = [r for r in rows if isnum(r.get("x_J60"))]
    res60 = {"n": len(sub60)}
    if len(sub60) >= 3:
        k60 = [(r["station"], r["season"]) for r in sub60]
        b60, _K = p2.cluster_boot(k60, {"RR": (np.array([r["y_R"] for r in sub60]), np.array([r["x_R"] for r in sub60])),
                                        "RJ60": (np.array([r["y_J"] for r in sub60]), np.array([r["x_J60"] for r in sub60]))},
                                  p5.B, p5.SEED)
        rr = sum(r["y_R"] for r in sub60) / sum(r["x_R"] for r in sub60)
        rj = sum(r["y_J"] for r in sub60) / sum(r["x_J60"] for r in sub60)
        res60.update({"R_RSS": rnd(rr), "R_JPL_60km": rnd(rj), "dR60": rnd(rr - rj),
                      "dR60_ci95": pm.ci95(b60["RR"] - b60["RJ60"])[0], "R_JPL_60km_ci95": pm.ci95(b60["RJ60"])[0]})
    out["resolution_control"] = res60
    if not out["evaluable"]:
        out["reason"] = f"n={n}（需≥{PAIR_MIN_N}）、簇={K}（需≥{PAIR_MIN_G}）、Σx_RSS≥0 抽样比例={frac:.4f}（需≤{p5.DEN_NONNEG_MAX}）"
    return out


# ======================================================================== 判读（Y15）
def group_label(blk):
    """(c) 分组读法：n≥10、簇≥5 才读；CI 下端 >1 →「超出 1」；上端 <1 →「低于 1」；否则「不能区分」。"""
    if not blk or blk.get("n", 0) < GROUP_MIN_N or (blk.get("clusters_station_season") or 0) < GROUP_MIN_G:
        return "不可评"
    lo, hi = blk.get("ci95_station_season") or [None, None]
    if lo is None or hi is None:
        return "不可评"
    if lo > 1.0:
        return "超出 1（CI 下端 >1）"
    if hi < 1.0:
        return "低于 1（CI 上端 <1）"
    return "不能区分"


def judge_overall(cat_j, pair):
    """Y15 总判读：对「P5 的 R2 来自 RSS 反演（链特有的）降雨伪淡化」。"""
    if cat_j == "不可评" or not pair.get("evaluable"):
        return {"answer": "不可评", "code": "U", "A1": None, "A2": None, "A3": None}
    D_lo = (pair.get("D_ci95") or [None, None])[0]
    g = pair.get("groups") or {}
    A1 = bool(D_lo is not None and D_lo > 0)
    A2 = None
    if g.get("evaluable"):
        dhi = (g.get("dmd_ci95") or [None, None])[1]
        A2 = bool(dhi is not None and dhi < 0)
    A3 = not cat_j.startswith("支持 R2")
    if A1 and A2 and A3:
        ans, code = "支持（P5 的 R2 来自 RSS 链特有的降雨伪淡化）", "S"
        if cat_j.startswith("支持 R1"):
            ans += "；JPL 本身支持 R1（强）"
    elif (not A3) and (not A1):
        ans, code = "反对（JPL 独立复现 R2，且 RSS 并不显著更淡）", "O"
    elif A1 and A2 and not A3:
        ans, code = "不能区分：链间差异存在且集中在雨中，但 JPL 仍支持 R2（伪淡化只解释一部分）", "N1"
    elif A3 and not (A1 and A2):
        ans, code = "不能区分：JPL 不再支持 R2，但与 RSS 的逐对差异不显著或不集中在雨中", "N2"
    elif (not A3) and A1 and not A2:
        ans, code = "不能区分：RSS 显著更淡但不集中在雨中（或雨中分组不可评），JPL 仍支持 R2", "N3"
    else:
        ans, code = "不能区分（其他组合）", "N4"
    res = None
    if code in ("S", "N1"):
        lo60 = ((pair.get("resolution_control") or {}).get("dR60_ci95") or [None, None])[0]
        res = "分辨率对照后仍成立" if (lo60 is not None and lo60 > 0) else "可能来自分辨率差异（分辨率对照 CI 下端 ≤0 或不可算）"
    return {"answer": ans, "code": code, "A1": A1, "A2": A2, "A3": A3, "resolution_note": res}


# ======================================================================== 下载（Y1–Y3、Y11、Y12）
def run_fetch(args, out_dir, log, sb, events, plan_only=False):
    np = _np()
    t_start = time.monotonic()
    cache = args.cache
    os.makedirs(cache, exist_ok=True)
    token = p5.read_token(args.token_file)
    ebs = p5.events_by_station(events)
    cand, same_as_p5, cnt = p5_candidates(sb, ebs, cache, args.p5_cache, log)
    if not same_as_p5:
        raise ValueError(f"Y1 候选轨与 P5 子集清单不一致：{cnt}")
    revs_need = {st: sorted({rev_of(g) for g in gs}) for st, gs in cand.items()}
    all_t = [t for gs in cand.values() for t in gs.values()]
    jg = jpl_granules(min(all_t) - 86400.0, max(all_t) + 86400.0, os.path.join(cache, "cmr_jpl.json"), log)
    unmatched = {st: [r for r in rv if r not in jg] for st, rv in revs_need.items()}
    by_gid = {}
    for st, rv in revs_need.items():
        for r in rv:
            if r in jg:
                by_gid.setdefault(jg[r][0], []).append(st)
    info = {"candidates_rss": {st: len(v) for st, v in cand.items()}, "candidates_same_as_p5_files": same_as_p5,
            "jpl_granules_needed": len(by_gid), "station_granules_needed": sum(len(v) for v in by_gid.values()),
            "rev_unmatched": {st: len(v) for st, v in unmatched.items()}}
    log.log(f"[JPL] {json.dumps(info, ensure_ascii=False)}", echo=True)
    net = Net(args.raw_root, args.rate_mbps, args.raw_tag, log)
    if plan_only:
        return plan_probe(args, out_dir, log, sb, by_gid, token, net, info, t_start)
    todo = []
    for gid, sts in sorted(by_gid.items()):
        need = [st for st in sts if not any(os.path.exists(os.path.join(cache, "smap_jpl", safe(st), f"{gid}.npz" + x))
                                              for x in ("", ".missing", ".nocover"))]
        if need:
            todo.append((gid, need))
    counts = {"ok": 0, "nocover": 0, "missing": 0, "multi_alt": 0}
    lock = threading.Lock()

    def one(item):
        gid, sts = item
        try:
            geo_b = net.get(geo_url(gid), token, f"{RAW_SUB}/geo/{gid}.geo.nc4", reuse=True)
        except p5.NotFound:
            return [(st, gid, "missing", None) for st in sts]
        geo = parse_geo(geo_b)
        res = []
        for st in sts:
            la, lo = sb[st]["_coord"]
            i, j, dist, n_alt = locate_swath(geo["lat"], geo["lon"], geo["row_time"], la, lo)
            if i is None:
                res.append((st, gid, "nocover", {"dist_km": rnd(dist, 3) if dist is not None else None}))
                continue
            b = net.get(sub_url(gid, i, j), token, f"{RAW_SUB}/{safe(st)}/{gid}.sub_i{i}_j{j}.nc4")
            sub = parse_sub(b)
            if abs(sub["lat"][1, 1] - geo["lat"][i, j]) > 1e-4 or abs(sub["lon"][1, 1] - geo["lon"][i, j]) > 1e-4 \
                    or abs(sub["row_time"][1] - geo["row_time"][j]) > 1e-3:
                raise ValueError(f"Y3 子集中心与整轨数组不符：{gid} {st}")
            meta = {"station": st, "gid": gid, "rev": rev_of(gid), "i": i, "j": j, "dist_km": rnd(dist, 3), "n_alt": n_alt}
            res.append((st, gid, "ok", (sub, meta)))
        return res

    n_done = 0
    if todo:
        net.start()
    try:
        with cf.ThreadPoolExecutor(max_workers=args.max_conc) as ex:
            for fut in cf.as_completed([ex.submit(one, it) for it in todo]):
                for st, gid, kind, val in fut.result():
                    d = os.path.join(cache, "smap_jpl", safe(st))
                    os.makedirs(d, exist_ok=True)
                    p = os.path.join(d, f"{gid}.npz")
                    with lock:
                        if kind == "ok":
                            save_npz(p, *val)
                            counts["ok"] += 1
                            counts["multi_alt"] += int(val[1]["n_alt"] > 0)
                        else:
                            open(p + "." + kind, "w").write(json.dumps(val or {}) + "\n")
                            counts[kind] += 1
                n_done += 1
                if n_done % 50 == 0:
                    el = time.monotonic() - t_start
                    log.log(f"[JPL] {n_done}/{len(todo)} 轨，{net.bytes / 1e6:.1f} MB，{net.bytes / 1e6 / max(el, 1):.2f} MB/s",
                            echo=True)
        smap_s = round(time.monotonic() - t_start, 1)
        # CMORPH（Y11）
        passes = {st: [(p[0], p[1], p[2]) for p in load_jpl_passes(cache, st)
                       if np.isfinite(p[3]["smap_sss"][1, 1])] for st in cand}
        need = p5.hours_needed(ebs, passes)
        cm = MultiCm([args.p5_cache, cache])
        by_hour = {}
        for st, hs in need.items():
            for h in hs:
                if h in cm.missing or cm.has(st, h):
                    continue
                by_hour.setdefault(h, []).append(st)
        cmi = {"station_hours_needed": sum(len(v) for v in need.values()), "hours_to_extract": len(by_hour)}
        log.log(f"[CMORPH] {json.dumps(cmi)}", echo=True)
        if by_hour:
            if net.relay is None:
                net.start()
            fx = rs.CmorphRaw(net.root, args.rate_mbps, owner="p7a", source_job=args.raw_tag, exc_cls=p2.DataSourceError,
                              heartbeat=False, bucket=net.bucket, manifest=net.man)
            wx = p5.WinExtractor({st: sb[st]["_coord"] for st in need}, cache)
            mpath = os.path.join(cache, "cmorph", "missing.jsonl")
            os.makedirs(os.path.dirname(mpath), exist_ok=True)
            for h in sorted(by_hour):
                kind, val = fx.fetch(h)
                why = None
                if kind == "ok":
                    try:
                        why = wx.extract(val, h, by_hour[h])
                    except Exception as e:
                        why = f"read_error:{type(e).__name__}"
                else:
                    why = val
                if why:
                    with open(mpath, "a", encoding="utf-8") as f:
                        f.write(json.dumps({"hour": h, "why": why}) + "\n")
            cmi.update({"raw_reused": fx.reused, "raw_downloaded": fx.downloaded, "bytes": fx.bytes})
            for st in need:
                a = os.path.join(args.p5_cache, "cmorph", safe(st), "grid.json")
                b = os.path.join(cache, "cmorph", safe(st), "grid.json")
                if os.path.exists(a) and os.path.exists(b) and json.load(open(a)) != json.load(open(b)):
                    raise ValueError(f"Y11 CMORPH 网格与 P5 不一致：{st}")
    finally:
        net.stop()
    man = {"version": VERSION, "run_utc": iso(time.time()), **info, "smap_todo": len(todo), "smap_counts": counts,
           "smap_seconds": smap_s, "opendap_bytes": net.bytes, "opendap_requests": net.n_req,
           "raw_geo_reused": net.reused, "cmorph": cmi,
           "seconds_total": round(time.monotonic() - t_start, 1)}
    jdump(man, os.path.join(cache, "fetch_manifest.json"))
    jdump(man, os.path.join(out_dir, "p7a_fetch.json"))
    return man


def plan_probe(args, out_dir, log, sb, by_gid, token, net, info, t_start):
    """plan：只取 .dmr 与坐标／时刻（lat、lon、row_time），不取任何盐度或质量位。原样返回也落 raw。"""
    np = _np()
    probe = {"token_read": True}
    gids = sorted(by_gid)
    net.start()
    try:
        gid0 = gids[len(gids) // 2]
        xml = p5.get_retry(dmr_url(gid0), token, 120)
        throttle(net.bucket, len(xml))
        net.keep(xml, dmr_url(gid0), f"{RAW_SUB}/plan/{gid0}.dmr.xml", "smap_jpl_opendap")
        probe["dmr_granule"] = gid0
        probe["dmr_dims"] = {v: dmr_dim_sizes(xml, v) for v in V2D + ("row_time",)}
        probe["dmr_dims_ok"] = all([s for _n, s in (probe["dmr_dims"][v] or [])] == [NCT, NAT] for v in V2D) and \
            [s for _n, s in (probe["dmr_dims"]["row_time"] or [])] == [NAT]
        per_st = {}
        for st in sorted({s for v in by_gid.values() for s in v}):
            gs = [g for g in gids if st in by_gid[g]]
            c = {}
            for gid in gs[len(gs) // 2: len(gs) // 2 + 3]:
                try:
                    geo = parse_geo(net.get(geo_url(gid), token, f"{RAW_SUB}/geo/{gid}.geo.nc4", reuse=True))
                    la, lo = sb[st]["_coord"]
                    i, j, dist, n_alt = locate_swath(geo["lat"], geo["lon"], geo["row_time"], la, lo)
                    r = {"i": i, "j": j, "dist_km": rnd(dist, 3), "n_alt": n_alt,
                         "geo_lon_range": [rnd(float(np.nanmin(geo["lon"])), 3), rnd(float(np.nanmax(geo["lon"])), 3)],
                         "row_time_utc": iso(float(geo["row_time"][j]) + EPOCH2015) if j is not None else None}
                    if i is not None:
                        b = net.get(sub_url(gid, i, j, ("lat", "lon")), token,
                                    f"{RAW_SUB}/plan/{safe(st)}/{gid}.coord_i{i}_j{j}.nc4")
                        s = parse_sub(b, ("lat", "lon"))
                        r["sub_center_matches_geo"] = bool(abs(s["lat"][1, 1] - geo["lat"][i, j]) < 1e-4 and
                                                           abs(s["lon"][1, 1] - geo["lon"][i, j]) < 1e-4)
                        r["sub_row_time_utc"] = [iso(float(x) + EPOCH2015) if np.isfinite(x) else None for x in s["row_time"]]
                    c[gid] = r
                except Exception as e:
                    c[gid] = {"error": p5.redact(f"{type(e).__name__}: {e}")[:300]}
            per_st[st] = c
            jdump({"partial": True, "probe": probe, "per_station": per_st}, os.path.join(out_dir, "p7a_plan_partial.json"))
    finally:
        net.stop()
    probe["per_station"] = per_st
    geo_bytes = net.bytes / max(1, sum(len(v) for v in per_st.values()))
    plan = {"script": os.path.basename(__file__), "version": VERSION, "mode": "plan",
            "run_utc": iso(time.time()), **info, "probe": probe, "opendap_bytes_plan": net.bytes,
            "est_download_MB": round((info["jpl_granules_needed"] * geo_bytes + info["station_granules_needed"] * 40e3) / 1e6, 1),
            "seconds": round(time.monotonic() - t_start, 1)}
    jdump(plan, os.path.join(out_dir, "p7a_plan.json"))
    ok = probe["dmr_dims_ok"] and all(any(r.get("sub_center_matches_geo") for r in c.values()) for c in per_st.values())
    return 0 if ok else 3


# ======================================================================== 事件一致性与 V0（Y13）
def read_p5_events(p5_results):
    with open(os.path.join(p5_results, "p5_events.csv"), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def event_set_check(events, p5_rows):
    keys = {(e["_st"], e["onset_utc"]) for e in events}
    miss = [(r["station"], r["onset_utc"]) for r in p5_rows if (r["station"], r["onset_utc"]) not in keys]
    return {"n_events": len(events), "n_p5_rows": len(p5_rows), "p5_rows_not_in_events": miss,
            "ok": len(events) == 317 and not miss}


def v0_check(rows_r, p5_rows, p5_R):
    got = {(r["station"], r["onset_utc"]): r for r in rows_r}
    ref = {(r["station"], r["onset_utc"]): r for r in p5_rows}
    bad = []
    for k, r in ref.items():
        g = got.get(k)
        if g is None:
            bad.append([*k, "missing"])
            continue
        for f in ("y", "x"):
            a, b = float(r[f]), g[f]
            if not (isnum(b) and abs(a - b) <= 1e-5 * max(1.0, abs(a))):
                bad.append([*k, f, a, b])
    extra = [list(k) for k in got if k not in ref]
    R = sum(r["y"] for r in rows_r) / sum(r["x"] for r in rows_r) if rows_r else NAN
    ok = not bad and not extra and isnum(R) and round(R, 5) == round(p5_R, 5)
    return {"ok": ok, "n": len(rows_r), "R": rnd(R), "R_p5_summary": p5_R, "bad": bad[:20], "extra": extra[:20]}


# ======================================================================== 分析
def analyze(args, out_dir, log, sb, events, info):
    np = _np()
    p5_rows = read_p5_events(args.p5_results)
    p5_sum = json.load(open(os.path.join(args.p5_results, "p5_summary.json")))
    es = event_set_check(events, p5_rows)
    sts = sorted({e["_st"] for e in events})
    cm_r = p5.CmCache(args.p5_cache)
    rss = {st: p5.load_passes(args.p5_cache, st) for st in sts}
    pe_r = p5.PassEval(sb, cm_r)
    rows_r, why_r = p5.compute_rows(sb, events, rss, pe_r, p5.PRIMARY, cm_r)
    v0 = v0_check(rows_r, p5_rows, p5_sum["primary"]["R"])
    v0.update({"event_set": es})
    jdump(v0, os.path.join(out_dir, "p7a_v0.json"))
    if not (v0["ok"] and es["ok"]):
        raise ValueError(f"Y13 V0 不符：{json.dumps({k: v0[k] for k in ('R', 'R_p5_summary', 'bad', 'extra')})[:600]}")
    log.log(f"[V0] RSS 主结果复算 R={v0['R']}（P5 {v0['R_p5_summary']}），87 事件逐条一致；事件集 317 一致", echo=True)
    rain = rss_rain_table(args.p5_cache, sts)
    cm = MultiCm([args.p5_cache, args.cache])
    jpl = {st: load_jpl_passes(args.cache, st) for st in sts}
    pe_j = JPassEval(sb, cm, rain)
    rows, why = p5.compute_rows(sb, events, jpl, pe_j, PRIMARY_J, cm)
    prs, why_p = pair_rows(events, rss, jpl, pe_r, pe_j)
    # 功效（Y16）：先写，只含模型 x、形式不确定度与计数
    power = p5.power_block(rows)
    power["excluded_reasons"] = why
    power["sigma_note"] = "JPL smap_sss_uncertainty＝似然 FWHM 估计（手册 §3.4.1），定义与 RSS 形式不确定度不同；只作信息"
    im = [r["imerg_c_post"] for r in rows]
    power["groups_imerg"] = {"n_rain": sum(1 for v in im if isnum(v) and v > RAIN_THR),
                             "n_dry": sum(1 for v in im if isnum(v) and v <= RAIN_THR),
                             "n_missing": sum(1 for v in im if not isnum(v))}
    power["pairs"] = {"n": len(prs), "clusters": len({(r["station"], r["season"]) for r in prs}),
                      "sum_x_R": rnd(sum(r["x_R"] for r in prs)), "sum_x_J": rnd(sum(r["x_J"] for r in prs)),
                      "n_rain": sum(1 for r in prs if isnum(r["imerg_post"]) and r["imerg_post"] > RAIN_THR),
                      "n_dry": sum(1 for r in prs if isnum(r["imerg_post"]) and r["imerg_post"] <= RAIN_THR),
                      "excluded": why_p}
    jdump(power, os.path.join(out_dir, "p7a_power.json"))
    # (a) 主结果
    blk, _ = p5.ratio_block(rows, "y", "x", "R_sfc,JPL 主（A70 类比、[0,12) h、参照 [−72,0) h、QC 位 0/5/7/8）")
    verdict = p5.judge(blk)
    # (c) 分组
    def grp(rr, f):
        a = [r for r in rr if isnum(r.get(f)) and r[f] > RAIN_THR]
        b = [r for r in rr if isnum(r.get(f)) and r[f] <= RAIN_THR]
        out = {}
        for name, sub in (("rain", a), ("dry", b)):
            bl, _ = p5.ratio_block(sub, "y", "x", f"{f} {name}")
            out[name] = {k: bl.get(k) for k in ("n", "clusters_station_season", "stations", "R", "ci95_station_season",
                                                 "ci95_station", "sum_num", "sum_den", "evaluable")}
            out[name]["category_p5_rule6"] = p5.judge(bl)["category"]
            out[name]["label_vs_1"] = group_label(bl)
        return out
    groups = {"jpl_by_imerg": grp(rows, "imerg_c_post"), "jpl_by_cmorph_foot": grp(rows, "rain_foot_cmorph_post"),
              "rss_by_imerg_V0": grp(rows_r, "imerg_c_post")}
    # (d) 逐对
    pair = pair_stats(prs)
    pair["excluded"] = why_p
    overall = judge_overall(verdict["category"], pair)
    # 敏感性
    sens = {}
    for code, label, kw in SENS_J:
        v = dict(PRIMARY_J, **kw)
        rr, ww = p5.compute_rows(sb, events, jpl, pe_j, v, cm)
        b, _ = p5.ratio_block(rr, "y", "x", label)
        j = p5.judge(b)
        sens[code] = {"label": label, "n": b["n"], "clusters": b.get("clusters_station_season"), "R": b.get("R"),
                      "ci95": b.get("ci95_station_season"), "evaluable": b.get("evaluable"), "category": j["category"],
                      "same_as_primary": j["category"] == verdict["category"], "excluded_reasons": ww}
    # 描述
    both = {(r["station"], r["onset_utc"]): r for r in rows_r}
    same_rev = [rev_of(r["post_gid"]) == rev_of(both[(r["station"], r["onset_utc"])]["post_gid"])
                for r in rows if (r["station"], r["onset_utc"]) in both]
    desc = {"per_station": {st: {"n": len(rr), "R": rnd(sum(r["y"] for r in rr) / sum(r["x"] for r in rr))
                                 if sum(r["x"] for r in rr) else None} for st, rr in p5._group(rows).items()},
            "dt_post_h_median": rnd(float(np.median([r["dt_post_h"] for r in rows]))) if rows else None,
            "n_ref_median": rnd(float(np.median([r["n_ref"] for r in rows]))) if rows else None,
            "events_in_both_primary": len(same_rev), "same_post_rev_as_rss": int(sum(same_rev)),
            "rss_excluded_reasons_V0": why_r,
            "passes_loaded_jpl": {st: len(v) for st, v in jpl.items()}}
    with open(os.path.join(out_dir, "p7a_events_jpl.csv"), "w", newline="", encoding="utf-8") as f:
        flds = ["station", "season", "onset_utc", "post_gid", "post_utc", "dt_post_h", "n_ref", "node", "sat_post", "sat_ref",
                "y", "F_post", "F_ref", "S0", "x", "x_pix0", "sig", "imerg_c_post", "rain_foot_cmorph_post", "dropped_w_post"]
        w = csv.DictWriter(f, fieldnames=flds, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (rnd(r[k], 6) if isinstance(r.get(k), float) else r.get(k)) for k in flds})
    with open(os.path.join(out_dir, "p7a_pairs.csv"), "w", newline="", encoding="utf-8") as f:
        flds = ["station", "season", "onset_utc", "post_rev", "dt_post_h", "n_ref", "node", "y_R", "x_R", "y_J", "x_J", "x_J60", "d",
                "imerg_post", "cmorph_foot_post_J"]
        w = csv.DictWriter(f, fieldnames=flds, extrasaction="ignore")
        w.writeheader()
        for r in prs:
            w.writerow({k: (rnd(r[k], 6) if isinstance(r.get(k), float) else r.get(k)) for k in flds})
    here = os.path.dirname(os.path.abspath(__file__))
    shas = {fn: p5.sha256_file(os.path.join(here, fn)) for fn in
            ("p7a_jpl_smap.py", "p5_sss_sat.py", "p1_events.py", "p1b_extend.py", "p2_rim_test.py", "p4_mech.py",
             "raw_store.py", "raw_fetch_list.py")}
    S = {"script": "p7a_jpl_smap.py", "version": VERSION, "run_utc": iso(time.time()),
         "code_sha256": shas, "events": info, "v0": {k: v0[k] for k in ("ok", "n", "R", "R_p5_summary")},
         "power_file": "p7a_power.json（先于本文件写盘）",
         "a_primary": blk, "a_verdict_rule6": verdict, "b_rain_flag": "不可执行：JPL L2B V5 无降雨标志位（手册 §6.2.24）",
         "c_groups": groups, "d_pairs": pair, "overall": overall, "sensitivity": sens, "descriptive": desc,
         "nature": "探索性；P5 RSS 结果已看过；对 JPL 盐度盲"}
    jdump(S, os.path.join(out_dir, "p7a_summary.json"))
    log.log(f"[P7a] R_JPL={blk.get('R')} CI={blk.get('ci95_station_season')} n={blk['n']} 类={verdict['category']}；"
            f"D={pair.get('D')} {pair.get('D_ci95')}；总判读={overall['answer']}", echo=True)
    return 0


def build(args, out_dir, log):
    return p5.build(args, out_dir, log)


# ======================================================================== 自测（合成，无网络）
def selftest(out_dir=None):
    np = _np()
    res = []

    def chk(name, ok, detail=""):
        res.append({"name": name, "ok": bool(ok), "detail": str(detail)[:300]})

    chk("Y1 轨道号解析", rev_of("RSS_SMAP_SSS_L2C_r21988_20190314T235118_2019073_FNL_V06.0") == 21988 and
        rev_of("SMAP_L2B_SSS_03258_20150911T072905_R17000_V5.0") == 3258)
    # Y5 QC
    sss = np.full((3, 3), 34.0)
    qf = np.zeros((3, 3))
    qf[0, 0] = 1
    qf[0, 1] = 1 << 6
    qf[0, 2] = 1 << 7
    qf[1, 0] = QF_FILL
    qf[2, 0] = np.nan
    qf[2, 1] = 1 << 5
    sss[2, 2] = 50.0
    ok = jqc_valid(sss, qf, JQC_BITS)
    ok7 = jqc_valid(sss, qf, (5, 7, 8))
    chk("Y5 QC 位", (not ok[0, 0]) and ok[0, 1] and not ok[0, 2] and not ok[1, 0] and not ok[2, 0] and not ok[2, 1]
        and not ok[2, 2] and ok[1, 1] and ok7[0, 0], ok.astype(int).tolist())
    # Y3 定位
    lat = np.full((NCT, NAT), np.nan)
    lon = np.full((NCT, NAT), np.nan)
    rt = np.full(NAT, np.nan)
    for jj in range(300, 700):
        for ii in range(10, 66):
            lat[ii, jj] = -10.0 + 0.225 * (jj - 300)
            lon[ii, jj] = 150.0 + 0.23 * (ii - 10)
        rt[jj] = 1e8 + 3.6 * jj
    i1, j1, d1, _ = locate_swath(lat, lon, rt, 0.1, 155.1)
    i2, _j2, _d2, _ = locate_swath(lat, lon, rt, 0.1, 175.0)
    chk("Y3 定位", i1 is not None and abs(lat[i1, j1] - 0.1) < 0.12 and abs(lon[i1, j1] - 155.1) < 0.12 and d1 < 20
        and i2 is None, f"i={i1} j={j1} d={d1}")
    # Y8 核
    grid = {"lat": [0.125 + 0.009 * (k - 120) for k in range(241)], "lon": [165.125 + 0.009 * (k - 120) for k in range(241)]}
    for fw, tr in ((40.0, 60.0), (60.0, 90.0)):
        w = jkernel(grid, 0.125, 165.125, fw, tr)
        la_ = np.asarray(grid["lat"])[:, None]
        lo_ = np.asarray(grid["lon"])[None, :]
        d = np.sqrt(((la_ - 0.125) * p5.KM_PER_DEG) ** 2 + ((lo_ - 165.125) * p5.KM_PER_DEG * math.cos(math.radians(0.125))) ** 2)
        ratio = float(np.median(w[np.abs(d - fw / 2) < 0.5] / w[120, 120]))
        chk(f"Y8 核 FWHM {fw:.0f}", abs(w.sum() - 1) < 1e-12 and abs(ratio - 0.5) < 0.06 and np.all(w[d > tr + 1e-9] == 0),
            f"ratio={ratio:.3f}")
    w40 = jkernel(grid, 0.2, 165.2, 40.0, 60.0)
    wp5 = p5.cell_kernel(grid, 0.2, 165.2)
    chk("Y8 主核与 P5 cell_kernel 逐值相同", np.allclose(w40, wp5, rtol=0, atol=1e-15))
    # Y12 限速：n>桶容量不死等
    b = rs.TokenBucket(50 * rs.MB)
    t0 = time.monotonic()
    throttle(b, 3 * 1024 * 1024)
    el = time.monotonic() - t0
    chk("Y12 分块扣令牌（3 MiB＠50 MB/s）", 0.02 < el < 1.5, f"{el:.3f}s")
    # Y15 判读映射
    def P(Dlo, dlo, ev=True, gev=True):
        return {"evaluable": ev, "D_ci95": [Dlo, Dlo + 1], "groups": {"evaluable": gev, "dmd_ci95": [-dlo - 1, -dlo]},
                "resolution_control": {"dR60_ci95": [0.1, 0.5]}}
    cases = [(judge_overall("支持 R1（海面也显著小于 RIM-3）", P(0.3, 0.2)), "S"),
             (judge_overall("不能区分", P(0.3, 0.2)), "S"),
             (judge_overall("支持 R2（海面与 RIM-3 相符或更大）", P(-0.2, 0.2)), "O"),
             (judge_overall("支持 R2（海面与 RIM-3 相符或更大）", P(0.3, 0.2)), "N1"),
             (judge_overall("不能区分", P(-0.1, -0.5)), "N2"),
             (judge_overall("支持 R2（海面与 RIM-3 相符或更大）", P(0.3, -0.5)), "N3"),
             (judge_overall("支持 R2（海面与 RIM-3 相符或更大）", P(0.3, 0.2, gev=False)), "N3"),
             (judge_overall("不能区分", P(0.3, 0.2, gev=False)), "N2"),
             (judge_overall("不可评", P(0.3, 0.2)), "U"),
             (judge_overall("支持 R2（海面与 RIM-3 相符或更大）", P(0.3, 0.2, ev=False)), "U")]
    chk("Y15 总判读映射", all(c["code"] == want for c, want in cases), [(c["code"], w_) for c, w_ in cases])
    chk("Y15 分组读法", group_label({"n": 12, "clusters_station_season": 6, "ci95_station_season": [1.1, 2]}).startswith("超出")
        and group_label({"n": 12, "clusters_station_season": 6, "ci95_station_season": [0.1, 0.9]}).startswith("低于")
        and group_label({"n": 9, "clusters_station_season": 6, "ci95_station_season": [1.1, 2]}) == "不可评")
    # 端到端
    try:
        e = _e2e(np)
        chk("E2E JPL R 复原（无噪声）", abs(e["R_J"] - 0.3) < 1e-9, e)
        chk("E2E 走 p5.event_eval 的窗选择（最早雨后 +3 h、参照 2 个）", e["post_dt"] == {3.0} and e["n_ref"] == {2}, e)
        chk("E2E 逐对 D、D_rain、D_dry 复原", abs(e["D"] - e["D_expect"]) < 1e-5 and abs(e["D_rain"] - 0.5) < 1e-5
            and abs(e["D_dry"]) < 1e-5 and abs(e["R_RSS_paired"] - e["R_RSS_expect"]) < 1e-5
            and e["dmd"] is not None and e["dmd"] < 0 and e["n60"] == e["n_pairs"], e)
        chk("E2E 事件集核对与 V0 核对函数", e["es_ok"] and e["v0_ok"] and not e["v0_bad_ok"], e)
        chk("E2E 缓存存取（npz）", e["npz_ok"], e)
    except Exception as ex:
        chk("E2E", False, f"{type(ex).__name__}: {ex} {traceback.format_exc()[-800:]}")
    n_ok = sum(r["ok"] for r in res)
    out = {"version": VERSION, "passed": n_ok, "total": len(res), "all_ok": n_ok == len(res), "checks": res}
    if out_dir:
        jdump(out, os.path.join(out_dir, "p7a_selftest.json"))
    for r in res:
        print(("  ok  " if r["ok"] else "  FAIL") + f" {r['name']}" + ("" if r["ok"] else f" :: {r['detail']}"), flush=True)
    print(f"[selftest] {n_ok}/{len(res)}", flush=True)
    return out


def _e2e(np):
    """合成：2 站×4 季×5 事件；过境 −50、−10、+3、+8 h（另 −80、+30 h 窗外）。JPL 与 RSS 的格位置、S0 相同 ⇒ x_R＝x_J。
    JPL 卫星＝34＋0.3·S0·(F̄−1)；RSS＝JPL＋0.5·S0·(F̄−1)·[该过境 IMERG>0.1]；奇数事件雨后过境 IMERG＝2 mm/h，其余 0。"""
    td = tempfile.mkdtemp(prefix="p7ae2e_")
    try:
        rng = np.random.default_rng(11)
        sts = {"A": (0.0, 165.0), "B": (15.0, 90.0)}
        sb, events = {}, []
        h_base, n_h = 400000, 24 * 400
        for st, (la, lo) in sts.items():
            ser = {k: np.full(n_h, 35.0) for k in ("s1", "sss05", "s5", "rain", "pco2", "sst")}
            ser["wind"] = np.full(n_h, 6.0)
            sb[st] = {"name": st, "h0": h_base, "n": n_h, "ser": ser, "_coord": (la, lo), "lat": la, "lon": lo}
        for st in sts:
            for s, day0 in (("DJF", 0), ("MAM", 90), ("JJA", 180), ("SON", 270)):
                for m in range(5):
                    H = h_base + (day0 + 12 * m + 5) * 24 + 7
                    events.append({"_st": st, "hour": H, "i": H - h_base, "season": s, "onset_utc": p1.hour_to_iso(H)})
        cache = os.path.join(td, "cache")
        clat = np.arange(-59.963, 60.0, 0.072771)
        clon = np.arange(0.0364, 360.0, 0.072756)

        class _DS:
            def __getitem__(self, k):
                return {"lat": clat, "lon": clon}[k]
        wx = p5.WinExtractor(sts, cache)
        wx._index(_DS())
        hours = set()
        for e in events:
            hours.update(range(e["hour"] - 110, e["hour"] + 40))
        for st in sts:
            rows_, cols_ = wx.idx[st]
            ev_h = {e["hour"] for e in events if e["_st"] == st}
            for h in sorted(hours):
                a = np.zeros((2, len(rows_), len(cols_)), np.float32)
                if any(0 <= h - H < 4 for H in ev_h):
                    a[:] = rng.gamma(1.0, 8.0, a.shape).astype(np.float32)
                np.save(os.path.join(cache, "cmorph", st, f"{h}.npy"), a)
        cm = MultiCm([cache])
        cmr = p5.CmCache(cache)
        rain = {}
        rss, jpl = {}, {}
        pe0 = p5.PassEval(sb, cmr)
        rev = 10000
        wet_post = {}
        for k, e in enumerate(events):
            wet_post[e["onset_utc"] + e["_st"]] = (k % 2 == 1)
        for st, (la, lo) in sts.items():
            lr, lj = [], []
            for e in [e for e in events if e["_st"] == st]:
                t0 = e["hour"] * 3600.0
                for dh in (-80.0, -50.0, -10.0, 3.0, 8.0, 30.0):
                    rev += 1
                    t = t0 + dh * 3600.0
                    cla_, clo_ = p5.station_cell(la, lo)
                    cla = np.array([[[cla_ + 0.25 * (i - 1)] * 2 for _j in range(3)] for i in range(3)])
                    clo = np.array([[[(clo_ - 0.25 * (j - 1)) % 360.0] * 2 for j in range(3)] for _i in range(3)])
                    wet = (dh == 3.0 and wet_post[e["onset_utc"] + st])
                    rain[(st, rev)] = 2.0 if wet else 0.0
                    subr = {"sss_smap_40km": np.full((3, 3, 2), 34.0), "sss_smap_40km_unc": np.full((3, 3, 2), 0.7),
                            "time": np.full((3, 3, 2), t - p5.EPOCH2000), "iqc_flag": np.zeros((3, 3, 2)),
                            "cellat": cla, "cellon": clo, "rain": np.full((3, 3), rain[(st, rev)]),
                            "sss_ref": np.full((3, 3), 34.0), "winspd": np.full((3, 3), 6.0), "surtep": np.full((3, 3), 300.0)}
                    gr = f"RSS_SMAP_SSS_L2C_r{rev:05d}_X"
                    gj = f"SMAP_L2B_SSS_{rev:05d}_X"
                    r0, why = pe0.evaluate(st, gr, subr, p5.PRIMARY)
                    if r0 is None:
                        raise RuntimeError(f"合成 RSS 过境无效：{why}")
                    base = 34.0 + 0.3 * 34.0 * (r0["F"] - 1.0)
                    subr["sss_smap_40km"][:] = base + (0.5 * 34.0 * (r0["F"] - 1.0) if wet else 0.0)
                    subj = {"smap_sss": np.full((3, 3), base), "smap_sss_uncertainty": np.full((3, 3), 0.8),
                            "quality_flag": np.zeros((3, 3)), "lat": cla[:, :, 0].copy(), "lon": clo[:, :, 0].copy(),
                            "anc_sss": np.full((3, 3), 34.0), "anc_spd": np.full((3, 3), 6.0), "smap_spd": np.full((3, 3), 6.0),
                            "anc_sst": np.full((3, 3), 300.0), "row_time": np.full(3, t - EPOCH2015)}
                    lr.append((t, t, gr, subr))
                    lj.append((t, t, gj, subj))
            rss[st] = sorted(lr, key=lambda p: p[0])
            jpl[st] = sorted(lj, key=lambda p: p[0])
        # npz 往返
        p = os.path.join(td, "x.npz")
        save_npz(p, jpl["A"][0][3], {"gid": jpl["A"][0][2]})
        back = load_npz(p)
        npz_ok = back["meta"]["gid"] == jpl["A"][0][2] and np.array_equal(back["smap_sss"], jpl["A"][0][3]["smap_sss"])
        pe_r = p5.PassEval(sb, cmr)
        pe_j = JPassEval(sb, cm, rain)
        rows_j, why = p5.compute_rows(sb, events, jpl, pe_j, PRIMARY_J, cm)
        blk, _ = p5.ratio_block(rows_j, "y", "x", "e2e")
        prs, _w = pair_rows(events, rss, jpl, pe_r, pe_j)
        ps = pair_stats(prs)
        xr = np.array([r["x_R"] for r in prs])
        wet = np.array([r["imerg_post"] > RAIN_THR for r in prs])
        D_expect = 0.5 * xr[wet].sum() / xr.sum()
        rows_r, _ = p5.compute_rows(sb, events, rss, pe_r, p5.PRIMARY, cmr)
        fake_csv = [{"station": r["station"], "onset_utc": r["onset_utc"], "y": f"{r['y']:.6f}", "x": f"{r['x']:.6f}"}
                    for r in rows_r]
        Rr = sum(r["y"] for r in rows_r) / sum(r["x"] for r in rows_r)
        v0 = v0_check(rows_r, fake_csv, round(Rr, 5))
        bad_csv = [dict(fake_csv[0], y="9.0")] + fake_csv[1:]
        v0b = v0_check(rows_r, bad_csv, round(Rr, 5))
        ev317 = [dict(e) for e in events] + [{"_st": "A", "onset_utc": f"pad{k}"} for k in range(317 - len(events))]
        es = event_set_check(ev317, fake_csv)
        return {"R_J": blk.get("R"), "n": len(rows_j), "why": why, "post_dt": {r["dt_post_h"] for r in rows_j},
                "n_ref": {r["n_ref"] for r in rows_j}, "D": ps.get("D"), "D_expect": float(D_expect),
                "D_rain": ps["groups"]["D_rain"], "D_dry": ps["groups"]["D_dry"],
                "R_RSS_paired": ps["R_RSS_paired"], "R_RSS_expect": 0.3 + float(D_expect),
                "dmd": ps["groups"]["dmd_psu"], "n60": ps["resolution_control"]["n"], "n_pairs": ps["n"],
                "es_ok": es["ok"], "v0_ok": v0["ok"], "v0_bad_ok": v0b["ok"], "npz_ok": bool(npz_ok)}
    finally:
        shutil.rmtree(td, ignore_errors=True)


# ======================================================================== main
def main(argv=None):
    ap = argparse.ArgumentParser(description="P7a：JPL SMAP L2B 独立链共址检验")
    g = ap.add_mutually_exclusive_group(required=True)
    for m in ("selftest", "plan", "fetch", "analyze", "full"):
        g.add_argument(f"--{m}", action="store_true")
    ap.add_argument("--out")
    ap.add_argument("--cache", default=CACHE_DEFAULT)
    ap.add_argument("--p5-cache", default=P5_CACHE_DEFAULT)
    ap.add_argument("--p5-results", default=P5_RESULTS_DEFAULT)
    ap.add_argument("--p1-events", default=p2.P1_EVENTS_DEFAULT)
    ap.add_argument("--p1b-dir", default=p2.P1B_DIR_DEFAULT)
    ap.add_argument("--token-file", default=p5.TOKEN_FILE_DEFAULT)
    ap.add_argument("--max-conc", type=int, default=MAX_CONC)
    ap.add_argument("--raw-root", default=rs.RAW_ROOT_DEFAULT)
    ap.add_argument("--rate-mbps", type=float, default=2.0)
    ap.add_argument("--raw-tag", default="p7a-jpl")
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
    if not (0 < args.rate_mbps <= RATE_MAX):
        print("限速越界：本任务只用 2 MB/s 名额", file=sys.stderr)
        return 3
    os.makedirs(out_dir, exist_ok=True)
    log = p1.Log(os.path.join(out_dir, "p7a_log.txt"))
    mode = next(m for m in ("plan", "fetch", "analyze", "full") if getattr(args, m))
    log.log(f"=== start {VERSION} mode={mode} out={out_dir} cache={args.cache}", echo=True)
    try:
        if not selftest(out_dir)["all_ok"]:
            log.log("自测不过，退出 4", echo=True)
            return 4
        sb, events, info = build(args, out_dir, log)
        es = event_set_check(events, read_p5_events(args.p5_results))
        log.log(f"[事件] {json.dumps({k: es[k] for k in ('n_events', 'n_p5_rows', 'ok')})}", echo=True)
        if not es["ok"]:
            log.log(f"事件集与 P5 不一致：{es['p5_rows_not_in_events'][:10]}，退出 3", echo=True)
            return 3
        if mode == "plan":
            rc = run_fetch(args, out_dir, log, sb, events, plan_only=True)
        else:
            if mode in ("fetch", "full"):
                run_fetch(args, out_dir, log, sb, events)
            rc = analyze(args, out_dir, log, sb, events, info) if mode in ("analyze", "full") else 0
    except p5.AuthError as e:
        log.log(f"FATAL 认证失败：{p5.redact(e)}", echo=True)
        return 5
    except WaitTimeout as e:
        log.log(f"FATAL {e}", echo=True)
        return 6
    except (p5.DataError, p2.DataSourceError, p1.FetchError) as e:
        log.log(f"FATAL 数据源故障（重跑即续传）：{p5.redact(e)}", echo=True)
        return 2
    except Exception:
        log.log("FATAL 未预期异常：\n" + p5.redact(traceback.format_exc()), echo=True)
        return 3
    log.log(f"=== done rc={rc}", echo=True)
    log.close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
