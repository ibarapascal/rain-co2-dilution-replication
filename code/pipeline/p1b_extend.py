#!/usr/bin/env python3
"""P1b：一次性扩样（R 组新时段＋O 组新站），在「P1 原 204 事件＋新增」合并集上重判 D1/D2/D4/D5。

设计（运行前写定），逐条：
  - 「新增数据」1（R 组）：8 站 NDBC OceanSITES 实时（R 模式）文件中，不与 P1 已用 D/M 文件时段重叠的小时；重叠小时丢弃 R。
  - 「新增数据」2（O 组）：PMEL ERDDAP `pmel_co2_moorings_*` 全部数据集，去掉 P1 的 8 站与 Chuuk，同址（≤25 km）有匿名可取
    的 ≤1 h 雨量与 ≤1.5 m 盐度者全部纳入；站单机械生成写入 p1b_stations.json。缺 1 m 盐度时以最浅 ≤1.5 m 盐度代「1 m」并标注。
  - 「处理与判定」：事件定义/质控/对照窗/ΔS/统计与 P1 相同（复用 `p1_events.py` 的 I1–I17 实现，I1 对 R 组按上条例外）；
    判定集＝P1 原 204 个合格事件（从 P1 产物 p1_events.csv 读入，不重新判定）＋新增合格事件；D1/D2/D4/D5 阈值同 P1；
    D4 过＝合并集中 ΔS(0.5 m) ≤ −0.2 psu 事件 ≥30。新增部分单独统计只作信息。
  - 本脚本不改任何事先写定的阈值；阈值常量全部取自 p1_events.py。

用法：
  python3 p1b_extend.py --r-rule union|per-var  # 全量（输出到 $REPRO_OUTPUT_DIR）
  python3 p1b_extend.py --r-rule union|per-var --out DIR   # 本地覆盖输出目录（--r-rule 必填，见 J2/J2b）
  python3 p1b_extend.py --plan --out DIR        # 只做元数据：O 组站单发现＋R 组覆盖小时估计（不下数据、不数事件）
  python3 p1b_extend.py --selftest              # 合成数据逻辑自测（无网络）
  --p1-events PATH   P1 产物 p1_events.csv（默认 p1 阶段产物；同目录 p1_summary.json 用于一致性校验）
  --p1-cache DIR     P1 缓存目录（默认 p1-events 同目录下 cache/；只读，用于复原 P1 文件清单与 D/M 文件时段、免重下）

依赖：仅 Python 3.12 标准库＋同目录 p1_events.py。数据全部匿名 HTTP（不注册任何账号）：
  - PMEL ERDDAP：allDatasets（站单）、tabledap `pmel_co2_moorings_*`（MAPCO2 3 小时 SSS/pCO₂）
  - NDBC OceanSITES GDAC：索引 oceansites_index.txt（O 组候选文件的位置/时段/参数/数据模式），THREDDS OPeNDAP（DDS/DAS/ASCII）
缓存：<out>/cache/（gzip，已下载跳过，可续跑）。产物：p1b_summary.json、p1b_events_new.csv、p1b_stations.json、p1b_log.txt。
退出码：0＝跑完（无论门过不过）；2＝数据源故障（重跑可续）；3＝其他异常（含 P1 csv 校验不一致）。

实现选择清单（设计未规定或有歧义之处；均未放宽任何门）：
  J1 P1 事件来源：只读 P1 产物 p1_events.csv（列已含 D2/D4/D5 所需：station/season/dS1_0_6h/n_ctrl/ctrl_dS1/pre_S1_median/
     pre_S05_median/dS05_0_6h），行数须＝204；同目录有 p1_summary.json 时，用 csv 复算 P1 的 D1/D2/D4/D5 并与之比对
     （计数须完全一致，D2/D5 数值容差见 P1_TOL），不一致即退出 3。csv 数值为 4 位小数，合并判定按此精度使用。缺列即退出 3。
  J2 R 组「P1 已用 D/M 文件时段」：用 P1 缓存里的 catalog 快照按 P1 同一函数（list_gtmba_files，include_rt=False）复原 P1 文件
     清单；每个文件时段＝其 TIME 最小值到最大值所在的 UTC 小时（闭区间）；8 站各自取所有 D/M 文件（RAIN/WIND/SALT 不分类）
     时段之并集为「P1 小时」。R 文件＝同一函数 include_rt=True 时被选中的 R 模式文件（即该部署×变量无 D/M 者）；
     R 值只写入非「P1 小时」；P1 小时对 R 组全为缺测。
  J2b（替代口径，`--r-rule per-var`）：逐变量判定重叠——某变量（雨/风/盐）在某小时落在 P1 该变量 D/M 文件时段内，则保留
     P1 D/M 值（由 P1 同一函数从 P1 缓存重载，I4 平均）、丢弃该变量 R 值；否则用 R 值补。新增事件＝合并序列上合格、窗口
     （雨 [t0−24, t0+24)，盐/风 [t0−6, t0+24)）内至少用到一个 R 值、且起点不在 P1 204 之列者；对照按 I9 在合并序列上选。
     字面口径（J2，`--r-rule union`）下 R 组新增小时约为 0；两口径差异巨大，全量运行必须显式指定（参考运行用 per-var）。
  J3 R 组对照（仅 union）：除 P1 的 I9 规则外，对照窗 [c−24 h, c+24 h) 不得触及任何「P1 小时」（P1 小时内雨量对 R 组不可见，无法按 I9 排除
     与 P1 雨事件重叠的对照）；事件窗 [t0−24 h, t0+24 h) 触及 P1 小时者剔除（由雨量有效性要求已自动满足，另作防御检查）。
  J4 O 组候选文件：GDAC 索引中与 MAPCO2 站点（ERDDAP 标称经纬度）距离 ≤25 km（索引经纬度范围中点，大圆距离）、数据模式 D/M
     （沿用 I1：不用 R/P 模式）、时段与该站 MAPCO2 覆盖期相交、时段 ≥2 天（排除 CTD 单次剖面）、参数含雨量
     （rainfall_rate / lwe_precipitation_rate）或风（wind_speed / eastward_wind+northward_wind）或（盐度且文件最小深度 ≤1.5 m）。
  J5 O 组逐文件判定（DDS/DAS 元数据）：时间分辨率＝(索引末−索引首)/(nt−1) ≤1 h（容差 2%）；nt > MAX_NT 的文件跳过（体量）；
     雨量单位须可换算为 mm/h（不可换算即弃该变量）；盐度层深取自深度坐标变量（positive=up 取负），每个文件取 1.0 m（±0.05）层，
     无则取该文件最浅的 0–1.5 m 层并标注「代 1 m」；5.0 m 层只记 has_5m（信息量）。无深度坐标而有逐时刻传感器深度变量
     <V>_H（如 IMOS SOFS 的 PSAL_H，标注 height/positive=down 而值为负）时，取其约 50 点抽样中位数的绝对值作层深。
  J6 「开阔洋」：O 组的雨量与盐度一律来自 OceanSITES（按其计划定义＝开阔洋长期参照站），不再另设水深/离岸距离判据；
     因此 CCE1/CCE2 若满足同址条件即纳入（需确认）。
  J7 O 组质控：沿用 I2（OceanSITES QC ∈ {1,2}）；文件无 QC 变量时（如 WHOI 的 WHOTS/Stratus 文件），按 I3 先例把 D/M 模式文件中
     非填充、在物理范围内的值视为已质控（需确认）。
  J8 O 组多源重叠：同一站同一变量按优先级逐小时「先到先得」（不做 P1 I4 的跨部署平均）：模式 D 先于 M → 原生时间步长短者先
     → 文件时段短者先 → 文件名；低优先级文件只补高优先级缺测的小时。
  J9 O 组风：有 wind_speed 用之；否则同一文件的 eastward/northward 风合成风速；只用于事件/对照的数据可用性要求（同 I7）。
  J10 O 组型态（regime）＝「O:<站名>」；D1 的「≥10 事件站」对合并集逐站计数（R 组事件计入其 P1 站）。
  J11 P1 的「无 0.5 m 分支」（I11 信息量）不在 P1b 重算（P1 csv 不含该分支事件）；D5 若在合并集不过，只报数，不在脚本里裁定（同 I11）。
  J12 网络失败（重试用尽）＝退出 2 可续跑；单个 O 组文件格式不符（缺时间坐标/深度等）＝跳过该文件并在 p1b_stations.json 记录原因。
  J13 站单由脚本运行时机械生成（--plan 与全量用同一函数）；GDAC 索引每日更新，全量运行时的站单以产物为准。

Change Log:
  2026-09-26 初版（复用 p1_events.py，仅给其 build_events 增可选 ctrl_ok 参数）。
"""

import argparse
import csv
import gzip
import json
import math
import os
import re
import shutil
import statistics
import sys
import time
import traceback
import urllib.parse
from array import array
from datetime import datetime, timedelta, timezone

import p1_events as p1

VERSION = "p1b-2026-09-26a"

import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
P1_EVENTS_DEFAULT = _rp.upstream("p1_events")  # [repro] 读 p1 阶段产物
P1_EXPECTED_N = 204
P1_TOL = {"d2_diff": 1e-3, "d2_ci": 2e-3, "d5_median": 1e-3, "d5_corr": 5e-3}
P1_CSV_NEED = ["station", "onset_utc", "season", "dS1_0_6h", "n_ctrl", "ctrl_dS1",
               "pre_S1_median", "pre_S05_median", "dS05_0_6h"]

CHUUK_ID = "pmel_co2_moorings_3e17_4ac9_c1f7"
OS_INDEX_URL = "https://dods.ndbc.noaa.gov/oceansites/oceansites_index.txt"
DODS_OS = p1.THREDDS + "/dodsC/oceansites/"

# ---- 事先写定的 P1b 参数 ----
CO_LOC_KM = 25.0
MAX_S1_DEPTH = 1.5
MAX_DT_S = 3600.0

# ---- 实现层参数 ----
DT_TOL = 1.02            # J5
MIN_SPAN_DAYS = 2.0      # J4
MAX_NT = 3_000_000       # J5
S1_NOMINAL = 1.0
DEPTH_EQ_TOL = 0.05
RAIN_SN = {"rainfall_rate", "lwe_precipitation_rate"}
SAL_SN = {"sea_water_practical_salinity", "sea_water_salinity"}
WSPD_SN = "wind_speed"
U_SN, V_SN = "eastward_wind", "northward_wind"
RANGES = {"rain": (-50.0, 500.0), "wind": (0.0, 60.0), "s1": (0.0, 45.0), "s5": (0.0, 45.0)}
MODE_RANK = {"D": 0, "M": 1}
NAN = float("nan")


class FormatError(Exception):
    """单个文件格式不符：跳过该文件（J12）。"""


# ======================================================================== HTTP（可从 P1 缓存播种）
class SeededHttp(p1.Http):
    """dest 不存在时，若 P1 缓存同相对路径下有文件，先复制过来（P1 缓存只读）。"""

    def __init__(self, log, cache, seed_cache=None):
        super().__init__(log)
        self.cache = os.path.abspath(cache)
        self.seed = os.path.abspath(seed_cache) if seed_cache and os.path.isdir(seed_cache) else None
        self.n_seeded = 0

    def get(self, url, dest, empty_ok=False):
        if self.seed and not os.path.exists(dest) and not os.path.exists(dest + ".empty"):
            rel = os.path.relpath(os.path.abspath(dest), self.cache)
            if not rel.startswith(".."):
                for suf in ("", ".empty"):
                    src = os.path.join(self.seed, rel) + suf
                    if os.path.exists(src):
                        os.makedirs(os.path.dirname(dest), exist_ok=True)
                        shutil.copyfile(src, dest + suf)
                        self.n_seeded += 1
                        break
        return super().get(url, dest, empty_ok)


# ======================================================================== 小工具
def haversine_km(la1, lo1, la2, lo2):
    p = math.radians
    x = (math.sin(p(la2 - la1) / 2) ** 2
         + math.cos(p(la1)) * math.cos(p(la2)) * math.sin(p(lo2 - lo1) / 2) ** 2)
    return 2 * 6371.0 * math.asin(math.sqrt(min(1.0, x)))


_TS_RE = re.compile(r"(\d{4})-(\d{2})-(\d{2})(?:[T ](\d{1,2}):(\d{2})(?::(\d{2}))?)?")


def parse_ts(s):
    """宽松解析 ISO 时间（容忍 T24:00:00）；返回 UTC 小时序号，失败 None。"""
    m = _TS_RE.search(s or "")
    if not m:
        return None
    y, mo, d = int(m[1]), int(m[2]), int(m[3])
    hh, mi, ss = int(m[4] or 0), int(m[5] or 0), int(m[6] or 0)
    try:
        dt = datetime(y, mo, d, tzinfo=timezone.utc) + timedelta(hours=hh, minutes=mi, seconds=ss)
    except ValueError:
        return None
    return int(dt.timestamp() // 3600)


_TU_RE = re.compile(r"^\s*(days?|hours?|minutes?|mins?|seconds?|secs?|s)\s+since\s+(\d{4})-(\d{1,2})-(\d{1,2})"
                    r"(?:[T ](\d{1,2}):(\d{1,2})(?::(\d{1,2})(?:\.\d*)?)?)?", re.I)
_TU_FACTOR = {"d": 1.0, "h": 1.0 / 24, "m": 1.0 / 1440, "s": 1.0 / 86400}


def time_units_to_days1950(units):
    """返回 (factor, offset)：days_since_1950 = v*factor + offset。"""
    m = _TU_RE.match(units or "")
    if not m:
        raise FormatError(f"无法解析时间单位 {units!r}")
    u = m[1].lower()
    f = _TU_FACTOR["m" if u.startswith("min") else u[0]]
    ref = datetime(int(m[2]), int(m[3]), int(m[4]), int(m[5] or 0), int(m[6] or 0), int(m[7] or 0))
    return f, (ref - datetime(1950, 1, 1)).total_seconds() / 86400.0


def rain_factor(units):
    """雨量单位 → mm/h 的乘子；不认识返回 None（J5）。"""
    u = (units or "").strip().lower().replace(" ", "").replace("**", "").replace("^", "").replace(".", "")
    u = u.replace("millimeters", "mm").replace("millimetres", "mm").replace("millimeter", "mm").replace("millimetre", "mm")
    table = {"mm/hr": 1.0, "mm/hour": 1.0, "mm/h": 1.0, "mmh-1": 1.0, "mmhr-1": 1.0, "mmhour-1": 1.0,
             "millimeters/hour": 1.0, "millimeter/hour": 1.0,
             "mm/min": 60.0, "mmmin-1": 60.0, "mm/s": 3600.0, "mms-1": 3600.0,
             "m/s": 3.6e6, "ms-1": 3.6e6, "kgm-2s-1": 3600.0, "kg/m2/s": 3600.0, "kgm-2/s": 3600.0}
    return table.get(u)


def hour_to_days1950(h):
    return (h + p1.EPOCH1950_H) / 24.0


# ======================================================================== DDS / DAS / ASCII
_DDS_DECL = re.compile(r"\b(Float32|Float64|Int8|UInt8|Int16|UInt16|Int32|UInt32|Int64|UInt64|Byte|String|Url)"
                       r"\s+([\w.\-]+)((?:\s*\[\s*[\w.\-]+\s*=\s*\d+\s*\])*)\s*;")
_DDS_DIM = re.compile(r"\[\s*([\w.\-]+)\s*=\s*(\d+)\s*\]")
_DDS_GRID = re.compile(r"Grid\s*\{\s*ARRAY:\s*\w+\s+([\w.\-]+)\s*\[")


def dds_parse(text):
    """返回 ({变量: [(维名, 长度)]}, grid 变量名集合)。首次声明为准（Grid 的 ARRAY 先于 MAPS）。"""
    vars_ = {}
    for m in _DDS_DECL.finditer(text):
        name = m.group(2)
        if name not in vars_:
            vars_[name] = [(d, int(n)) for d, n in _DDS_DIM.findall(m.group(3))]
    return vars_, set(_DDS_GRID.findall(text))


_DAS_OPEN = re.compile(r"^([\w.\-]+)\s*\{\s*$")
_DAS_ATTR = re.compile(r"^(\w+)\s+([\w.\-]+)\s+(.*);\s*$")


def das_parse(text):
    out, stack = {}, []
    for line in text.splitlines():
        s = line.strip()
        m = _DAS_OPEN.match(s)
        if m:
            stack.append(m.group(1))
            if len(stack) == 2:
                out.setdefault(m.group(1), {})
            continue
        if s == "}":
            if stack:
                stack.pop()
            continue
        if len(stack) == 2:
            m = _DAS_ATTR.match(s)
            if m:
                v = m.group(3).strip()
                if len(v) >= 2 and v[0] == '"' and v[-1] == '"':
                    v = v[1:-1]
                out[stack[1]].setdefault(m.group(2), v)
    return out


def parse_ascii_values(path, name):
    """读 .ascii 响应里某个变量（数组或标量）的全部数值。"""
    got = p1.parse_dap_ascii(path, {name}).get(name)
    if got:
        return list(got)
    in_data = False
    with gzip.open(path, "rt", encoding="latin-1") as f:
        for line in f:
            if not in_data:
                in_data = line.startswith("-----")
                continue
            m = re.match(r"^\s*([\w.\-]+)\s*,\s*(\S+)\s*$", line)
            if m and m.group(1).split(".")[-1] == name:
                try:
                    return [float(m.group(2))]
                except ValueError:
                    return []
    return []


def q(s):
    return urllib.parse.quote(s, safe=":,._-")


# ======================================================================== P1 事件（J1）
def _f(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return NAN
    return v


def read_p1_events(path):
    with open(path, encoding="utf-8", newline="") as f:
        rd = csv.DictReader(f)
        missing = [c for c in P1_CSV_NEED if c not in (rd.fieldnames or [])]
        if missing:
            raise RuntimeError(f"P1 csv 缺列 {missing}：需按 P1 代码重算（本版未实现此路径，停止）")
        evs = []
        for r in rd:
            n_ctrl = int(r["n_ctrl"])
            cds = [(_f(x) if x != "" else NAN) for x in r["ctrl_dS1"].split(";")] if r["ctrl_dS1"] else []
            if len(cds) != n_ctrl:
                raise RuntimeError(f"P1 csv 行 {r['station']} {r['onset_utc']}：ctrl_dS1 个数 {len(cds)} ≠ n_ctrl {n_ctrl}")
            evs.append({"group": "P1", "station": r["station"], "onset_utc": r["onset_utc"], "season": r["season"],
                        "ds1": _f(r["dS1_0_6h"]), "pre_s1": _f(r["pre_S1_median"]),
                        "ds05": _f(r["dS05_0_6h"]), "pre_s05": _f(r["pre_S05_median"]),
                        "n_ctrl": n_ctrl, "ctrl_ds1": cds})
    return evs


def gates(evs):
    """在给定事件集上判 D1/D2/D4/D5（阈值取自 p1_events）。evs: 规范化事件 dict 列表。"""
    by_st = {}
    for e in evs:
        by_st[e["station"]] = by_st.get(e["station"], 0) + 1
    n = len(evs)
    n10 = sum(1 for v in by_st.values() if v >= p1.D1_MIN_EVENTS_PER_STATION)
    d1 = {"n_events_total": n, "n_stations_ge10": n10, "events_by_station": by_st,
          "threshold": f"总数 ≥{p1.D1_MIN_TOTAL} 且 ≥{p1.D1_MIN_EVENTS_PER_STATION} 事件的站 ≥{p1.D1_MIN_STATIONS}",
          "pass": n >= p1.D1_MIN_TOTAL and n10 >= p1.D1_MIN_STATIONS}
    full = [(e["station"], e["season"], e["ds1"], e["ctrl_ds1"]) for e in evs if e["n_ctrl"] == p1.CTRL_PER_EVENT]
    ge1 = [(e["station"], e["season"], e["ds1"], e["ctrl_ds1"]) for e in evs if e["n_ctrl"] >= 1]
    d2 = p1.d2_test(full, "主判：配满 5 个对照的事件")
    d2["threshold"] = f"|事件均值−对照均值| ≥{p1.D2_MIN_ABS_DIFF} psu 且 95% CI 不含 0"
    d2["sensitivity_ge1ctrl_info"] = p1.d2_test(ge1, "敏感性：≥1 个对照（不参与判定）")
    d5m = p1.d5_metrics(evs)
    if d5m["median_pre_diff_psu"] is None or d5m["corr_event_anom"] is None:
        d5 = dict(d5m, evaluable=False, reason="配对不足（I16）")
        d5["pass"] = False
    else:
        fail = (abs(d5m["median_pre_diff_psu"]) > p1.D5_MAX_ABS_MEDDIFF or d5m["corr_event_anom"] < p1.D5_MIN_CORR)
        d5 = dict(d5m, evaluable=True)
        d5["pass"] = not fail
    d5["threshold"] = f"|雨前中位差| ≤{p1.D5_MAX_ABS_MEDDIFF} psu 且 事件异常相关 ≥{p1.D5_MIN_CORR}"
    n_d4 = sum(1 for e in evs if e["ds05"] == e["ds05"] and e["ds05"] <= p1.D4_DS_THRESH)
    d4 = {"n_events_dS05_le_-0.2": n_d4, "threshold": f"≥{p1.D4_MIN_EVENTS}", "pass": n_d4 >= p1.D4_MIN_EVENTS}
    return {"D1": d1, "D2": d2, "D4": d4, "D5": d5}


def validate_p1(evs, summary_path, log):
    """J1：csv 复算 P1 门并与 p1_summary.json 比对。返回校验记录；不一致抛 RuntimeError。"""
    rec = {"n_rows": len(evs), "expected": P1_EXPECTED_N}
    if len(evs) != P1_EXPECTED_N:
        raise RuntimeError(f"P1 csv 行数 {len(evs)} ≠ {P1_EXPECTED_N}")
    keys = [(e["station"], e["onset_utc"]) for e in evs]
    if len(set(keys)) != len(keys):
        raise RuntimeError("P1 csv 存在重复 (station, onset_utc)")
    g = gates(evs)
    rec["recomputed"] = {"D1_n": g["D1"]["n_events_total"], "D1_st10": g["D1"]["n_stations_ge10"],
                         "D2_diff": g["D2"].get("diff_psu"), "D2_ci": g["D2"].get("ci95_cluster_bootstrap"),
                         "D4_n": g["D4"]["n_events_dS05_le_-0.2"],
                         "D5_median": g["D5"].get("median_pre_diff_psu"), "D5_corr": g["D5"].get("corr_event_anom")}
    if not os.path.exists(summary_path):
        rec["p1_summary"] = "缺失：未比对（仅校验行数与唯一性）"
        log.log(f"WARN P1 summary 不存在：{summary_path}")
        return rec
    with open(summary_path, encoding="utf-8") as f:
        s = json.load(f)["gates"]
    ref = {"D1_n": s["D1"]["n_events_total"], "D1_st10": s["D1"]["n_stations_ge10"],
           "D2_diff": s["D2"].get("diff_psu"), "D2_ci": s["D2"].get("ci95_cluster_bootstrap"),
           "D4_n": s["D4_info"]["n_events_dS05_le_-0.2"],
           "D5_median": s["D5"].get("median_pre_diff_psu"), "D5_corr": s["D5"].get("corr_event_anom")}
    rec["p1_summary"] = ref
    r = rec["recomputed"]
    bad = []
    for k in ("D1_n", "D1_st10", "D4_n"):
        if r[k] != ref[k]:
            bad.append(k)
    for k, tol in (("D2_diff", P1_TOL["d2_diff"]), ("D5_median", P1_TOL["d5_median"]), ("D5_corr", P1_TOL["d5_corr"])):
        if r[k] is None or ref[k] is None or abs(r[k] - ref[k]) > tol:
            bad.append(k)
    if not r["D2_ci"] or not ref["D2_ci"] or max(abs(a - b) for a, b in zip(r["D2_ci"], ref["D2_ci"])) > P1_TOL["d2_ci"]:
        bad.append("D2_ci")
    rec["mismatch"] = bad
    if bad:
        raise RuntimeError(f"P1 csv 复算与 p1_summary 不一致：{bad}；recomputed={r} ref={ref}")
    return rec


# ======================================================================== 事件（沿用 P1 逻辑）
def prefix(mask):
    ps = array("i", [0]) * (len(mask) + 1)
    acc = 0
    for i, v in enumerate(mask):
        acc += v
        ps[i + 1] = acc
    return ps


def masked_any(ps, a, b):
    a, b = max(0, a), min(len(ps) - 1, b)
    return b > a and ps[b] - ps[a] > 0


def run_events(series, n, h0, st_name, regime, lon, group, log, ps_mask=None, s1_depth=None):
    """对一个站的小时序列跑 P1 的事件/对照（I5–I9），返回规范化事件 dict 列表与计数。"""
    dry_run, onsets, ostats = p1.find_events(series, n)
    ctrl_ok = None
    if ps_mask is not None:
        def ctrl_ok(c):  # J3
            return not masked_any(ps_mask, c - p1.DRY_BEFORE_H, c + p1.POST_H)
    raw, n_drop = p1.build_events(series, n, dry_run, onsets, ["sss05", "pco2", "s1", "wind"], h0, ctrl_ok)
    evs, n_guard = [], 0
    s5, pco2 = series["s5"], series["pco2"]
    for e in raw:
        i = e["i"]
        if ps_mask is not None and masked_any(ps_mask, i - p1.DRY_BEFORE_H, i + p1.POST_H):
            n_guard += 1
            continue
        dep = S1_NOMINAL
        if s1_depth is not None:
            ds = [s1_depth[k] for k in range(max(0, i - p1.PRE_H), min(n, i + p1.DS_WIN_H)) if s1_depth[k] == s1_depth[k]]
            dep = statistics.median(ds) if ds else NAN
        evs.append({
            "group": group, "station": st_name, "regime": regime, "hour": e["hour"],
            "onset_utc": p1.hour_to_iso(e["hour"]), "season": p1.season_of(e["hour"]),
            "local_solar_hour": round((e["hour"] % 24 + lon / 15.0) % 24, 2),
            "acc24": e["acc24"], "first_mm": e["first_mm"], "n_wet": e["n_wet"],
            "ds1": e["ds1"], "pre_s1": e["pre_s1"], "ds05": e.get("ds05", NAN), "pre_s05": e.get("pre_s05", NAN),
            "n_ctrl": len(e["ctrl"]), "ctrl_ds1": e["ctrl_ds1"], "ctrl_hours": [h0 + c for c in e["ctrl"]],
            "has_5m": p1.seg_has(s5, i - p1.PRE_H, i) and p1.seg_has(s5, i, i + p1.DS_WIN_H),
            "n_pco2": sum(1 for x in pco2[max(0, i - p1.PRE_H):i + p1.POST_H] if x == x),
            "s1_depth_m": dep,
            "s1_substitute": not (dep == dep and abs(dep - S1_NOMINAL) <= DEPTH_EQ_TOL),
        })
    stats = {"rain_onset_stats": ostats, "n_rain_events_10mm": len(onsets), "n_events": len(evs),
             "n_dropped_insufficient_data": n_drop, "n_dropped_touch_p1_hours": n_guard,
             "n_events_5ctrl": sum(1 for e in evs if e["n_ctrl"] == p1.CTRL_PER_EVENT)}
    log.log(f"[{group}:{st_name}] ≥10mm 雨事件 {len(onsets)}；新增合格 {len(evs)}（数据不足剔 {n_drop}，触及P1小时剔 {n_guard}）")
    return evs, stats


def new_axis(hs, he):
    h0 = hs - p1.PAD_H
    return h0, he + p1.PAD_H - h0 + 1


def mapco2_series(mrows, h0, n):
    ser = {k: p1.Hourly(h0, n) for k in ("sss05", "pco2")}
    for h, sss, pc in mrows:
        if sss is not None:
            ser["sss05"].add(h, sss)
        if pc is not None:
            ser["pco2"].add(h, pc)
    return ser


# ======================================================================== R 组（J2/J3）
def file_mode(base):
    m = re.search(r"_([DMR])_", base.rsplit("/", 1)[-1])
    return m.group(1) if m else "?"


def dm_file_span(http, cache, base, log):
    """P1 D/M 文件的 TIME 首末小时。优先读（P1 播种的）完整 ascii 缓存，否则只取 TIME。"""
    fname = base.rsplit("/", 1)[-1]
    sub = os.path.join(cache, "thredds", fname.split("_")[1])
    full = os.path.join(sub, fname + ".ascii.gz")
    path = full
    if not os.path.exists(full):
        seeded = os.path.join(http.seed, os.path.relpath(full, http.cache)) if http.seed else None
        if seeded and os.path.exists(seeded):
            path = seeded  # 直接读 P1 缓存（只读）
        else:
            path = os.path.join(sub, fname + ".TIME.ascii.gz")
            http.get(base + ".ascii?TIME", path)
    t = p1.parse_dap_ascii(path, {"TIME"}).get("TIME")
    if not t:
        raise p1.FetchError(f"{fname} 取不到 TIME")
    return p1.days1950_to_hour(min(t)), p1.days1950_to_hour(max(t))


def process_r_station(http, st, cache, log, rule, p1_keys):
    """R 组。rule="union"：J2 字面口径（任何 D/M 文件时段内的小时对 R 全缺测）；
    rule="per-var"：J2b 逐变量口径（某变量有 D/M 文件时段的小时保留 P1 D/M 值、丢弃该变量 R 值，其余小时用 R 补；
    P1 D/M 数据由 P1 同一函数从 P1 缓存重载；新增事件＝窗口内用到 R 值、且不在 P1 204 之列的合格事件）。"""
    hs, he, mrows = p1.fetch_mapco2(http, st, cache, log)
    h0, n = new_axis(hs, he)
    ser = mapco2_series(mrows, h0, n)
    del mrows
    dm = p1.list_gtmba_files(http, st, cache, hs, he, False, include_rt=False)
    allf = p1.list_gtmba_files(http, st, cache, hs, he, False, include_rt=True)
    rfiles = [f for f in allf if file_mode(f[0]) == "R"]
    keys = ("rain", "wind", "s1", "s5")
    for k in keys:
        ser[k] = p1.Hourly(h0, n)
    mask = bytearray(n)                          # union：任一 D/M 文件时段
    vmask = {k: bytearray(n) for k in keys}      # per-var：该变量的 D/M 文件时段
    spans, qc_drop = [], {}
    for base, kind, per_hour, _sh in dm:
        a, b = dm_file_span(http, cache, base, log)
        spans.append([base.rsplit("/", 1)[-1], p1.hour_to_iso(a), p1.hour_to_iso(b)])
        lo_i, hi_i = max(0, a - h0), min(n, b - h0 + 1)
        for i in range(lo_i, hi_i):
            mask[i] = 1
        if rule == "per-var":
            got = p1.load_gtmba_file(http, cache, base, kind, per_hour, log, {})  # 与 P1 相同（I4 平均）
            for key, pairs in got.items():
                for i in range(lo_i, hi_i):
                    vmask[key][i] = 1
                for h, v in pairs:
                    ser[key].add(h, v)
            del got
    rflag = {k: bytearray(n) for k in keys}
    added, dropped = {}, {}
    for base, kind, per_hour, _sh in rfiles:
        got = p1.load_gtmba_file(http, cache, base, kind, per_hour, log, qc_drop)
        for key, pairs in got.items():
            m = mask if rule == "union" else vmask[key]
            for h, v in pairs:
                i = h - h0
                if not (0 <= i < n):
                    continue
                if m[i]:
                    dropped[key] = dropped.get(key, 0) + 1
                    continue
                ser[key].add(h, v)
                rflag[key][i] = 1
                added[key] = added.get(key, 0) + 1
    series = {k: v.values() for k, v in ser.items()}
    del ser
    need = ("rain", "s1", "wind", "sss05", "pco2")
    joint = sum(1 for i in range(n) if (rflag["rain"][i] or rflag["s1"][i] or rflag["wind"][i])
                and all(series[k][i] == series[k][i] for k in need))
    cov = {"rule": rule, "axis_hours": n, "mapco2_start": p1.hour_to_iso(hs), "mapco2_end": p1.hour_to_iso(he),
           "p1_dm_files": len(dm), "p1_hours_in_axis_union": sum(mask),
           "p1_hours_in_axis_by_var": {k: sum(v) for k, v in vmask.items()} if rule == "per-var" else None,
           "r_files": [f[0].rsplit("/", 1)[-1] for f in rfiles],
           "r_hour_values_added_by_var": added, "r_hour_values_dropped_overlap_p1": dropped,
           "new_hours_all5_valid": joint,
           "qc_dropped_by_code": {k: {str(c): x for c, x in v.items()} for k, v in qc_drop.items()},
           "p1_dm_file_spans": spans}
    if rule == "union":
        evs, stats = run_events(series, n, h0, st["name"], st["regime"], st["lon"], "R", log, ps_mask=prefix(mask))
        return evs, dict(stats, coverage=cov)
    allev, stats = run_events(series, n, h0, st["name"], st["regime"], st["lon"], "R", log)
    pr = {k: prefix(v) for k, v in rflag.items()}
    evs, recon, n_dup = [], [], 0
    for e in allev:
        i = e["hour"] - h0
        uses_r = (masked_any(pr["rain"], i - p1.DRY_BEFORE_H, i + p1.POST_H)
                  or any(masked_any(pr[k], i - p1.PRE_H, i + p1.POST_H) for k in ("s1", "wind")))
        if (st["name"], e["onset_utc"]) in p1_keys:
            n_dup += 1
            if uses_r:
                log.log(f"  WARN {st['name']} {e['onset_utc']} 与 P1 事件同起点且窗内含 R 值（保留 P1 原判，不计新增）")
            continue
        if uses_r:
            evs.append(e)
        else:
            recon.append(e["onset_utc"])
    p1_here = {o for s_, o in p1_keys if s_ == st["name"]}
    cov["p1_reconstruction_check"] = {
        "combined_events_matching_p1": n_dup, "p1_events_this_station": len(p1_here),
        "combined_events_without_R_not_in_p1": recon[:20], "n_combined_events_without_R_not_in_p1": len(recon),
        "note": "无 R 值的合并事件应与 P1 事件集一致；差异来自 R 雨改变对照/前 24 h 干燥判定或重建偏差，只作诊断"}
    stats = dict(stats, n_events=len(evs), n_events_5ctrl=sum(1 for e in evs if e["n_ctrl"] == p1.CTRL_PER_EVENT),
                 n_combined_qualified_all=len(allev))
    log.log(f"[R:{st['name']}] per-var：合并序列合格 {len(allev)}，其中新增（用到 R）{len(evs)}，与 P1 同起点 {n_dup}，"
            f"无 R 且非 P1 {len(recon)}")
    return evs, dict(stats, coverage=cov)


# ======================================================================== O 组发现（J4–J7）
def fetch_all_co2(http, cache):
    path = os.path.join(cache, "erddap", "allDatasets_pmel_co2_moorings.csv.gz")
    qs = ("datasetID,title,minLongitude,maxLongitude,minLatitude,maxLatitude,minTime,maxTime"
          '&datasetID=~"pmel_co2_moorings_.*"')
    http.get(f"{p1.ERDDAP}/tabledap/allDatasets.csv?" + urllib.parse.quote(qs, safe="=&,~"), path)
    rows = list(csv.reader(p1.read_gz_text(path).splitlines()))
    hdr = rows[0]
    out = []
    for r in rows[2:]:
        if len(r) != len(hdr):
            continue
        d = dict(zip(hdr, r))
        out.append({"erddap": d["datasetID"], "title": d["title"], "name": d["title"].split(" NOAA")[0].strip(),
                    "lat": (float(d["minLatitude"]) + float(d["maxLatitude"])) / 2,
                    "lon": (float(d["minLongitude"]) + float(d["maxLongitude"])) / 2,
                    "t0": parse_ts(d["minTime"]), "t1": parse_ts(d["maxTime"])})
    return out


def scan_os_index(http, cache, sites, log):
    """流式扫 GDAC 索引，按 J4 返回 {站名: [条目]}，另返回 {文件名: (首小时, 末小时)} 供 R 组估计。"""
    path = os.path.join(cache, "oceansites", "oceansites_index.txt.gz")
    http.get(OS_INDEX_URL, path)
    hits = {s["name"]: {} for s in sites}
    spans = {}
    n_lines = 0
    with gzip.open(path, "rt", encoding="latin-1") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fld = line.rstrip("\r\n").split(",", 15)
            if len(fld) < 16:
                continue
            n_lines += 1
            relp = fld[0].strip()
            base = relp.rsplit("/", 1)[-1]
            a, b = parse_ts(fld[2]), parse_ts(fld[3])
            if a is None or b is None:
                continue
            spans.setdefault(base, (a, b))
            try:
                la = (float(fld[4]) + float(fld[5])) / 2
                lo = (float(fld[6]) + float(fld[7])) / 2
            except ValueError:
                continue
            try:
                dmin = float(fld[8])
            except ValueError:
                dmin = None
            mode = fld[14].strip()
            params = set(fld[15].split())
            for s in sites:
                if haversine_km(s["lat"], s["lon"], la, lo) > CO_LOC_KM:
                    continue
                if mode not in MODE_RANK or b < s["t0"] or a > s["t1"] or (b - a) < MIN_SPAN_DAYS * 24:
                    continue
                has_rain = bool(params & RAIN_SN)
                has_wind = WSPD_SN in params or {U_SN, V_SN} <= params
                has_sal = bool(params & SAL_SN) and dmin is not None and dmin <= MAX_S1_DEPTH
                if not (has_rain or has_wind or has_sal):
                    continue
                cur = hits[s["name"]].get(base)
                if cur is None or (cur["path"].startswith("DATA_GRIDDED") and relp.startswith("DATA/")):
                    hits[s["name"]][base] = {"path": relp, "file": base, "h0": a, "h1": b, "mode": mode,
                                             "dist_km": round(haversine_km(s["lat"], s["lon"], la, lo), 2),
                                             "index_min_depth": dmin, "index_params": sorted(params)}
    log.log(f"GDAC 索引 {n_lines} 行")
    return {k: sorted(v.values(), key=lambda e: e["file"]) for k, v in hits.items()}, spans


def _depthish(name, attrs):
    a = attrs.get(name, {})
    return (a.get("standard_name") == "depth" or a.get("axis") == "Z" or name.upper().startswith("DEPTH"))


def _depth_values(http, sub, url, name, attrs):
    path = os.path.join(sub, f"{os.path.basename(url)}.{name}.ascii.gz")
    http.get(url + ".ascii?" + q(name), path)
    vals = parse_ascii_values(path, name)
    if (attrs.get(name, {}).get("positive") or "").lower() == "up":
        vals = [-v for v in vals]
    return vals


def _qc_name(v, vars_, attrs):
    for c in (attrs.get(v, {}).get("ancillary_variables") or "").split():
        if c.endswith("_QC") and c in vars_ and [d for d in vars_[c]] == vars_[v]:
            return c
    c = v + "_QC"
    return c if c in vars_ and vars_[c] == vars_[v] else None


def analyze_file(http, cache, ent, log):
    """J5：由 DDS/DAS（及深度坐标值）判定一个候选文件能提供哪些源。返回 (源列表, 文件记录)。"""
    url = DODS_OS + ent["path"]
    sub = os.path.join(cache, "oceansites", ent["path"].rsplit("/", 1)[0])
    rec = {"file": ent["file"], "path": ent["path"], "mode": ent["mode"], "dist_km": ent["dist_km"],
           "start": p1.hour_to_iso(ent["h0"]), "end": p1.hour_to_iso(ent["h1"]), "sources": [], "rejected": []}
    dds_p = os.path.join(sub, ent["file"] + ".dds.gz")
    das_p = os.path.join(sub, ent["file"] + ".das.gz")
    http.get(url + ".dds", dds_p)
    http.get(url + ".das", das_p)
    vars_, grids = dds_parse(p1.read_gz_text(dds_p))
    attrs = das_parse(p1.read_gz_text(das_p))
    tname = None
    for v, dims in vars_.items():
        if len(dims) == 1 and dims[0][0] == v and (v.upper() == "TIME" or attrs.get(v, {}).get("standard_name") == "time"):
            tname = v
            break
    if tname is None:
        raise FormatError("无 TIME 坐标")
    nt = vars_[tname][0][1]
    rec["nt"] = nt
    if nt < 2:
        raise FormatError("nt<2")
    dt = (ent["h1"] - ent["h0"]) * 3600.0 / (nt - 1)
    rec["dt_approx_s"] = round(dt, 1)
    if dt > MAX_DT_S * DT_TOL:
        raise FormatError(f"时间分辨率约 {dt:.0f} s > 1 h")
    if nt > MAX_NT:
        raise FormatError(f"nt={nt} > {MAX_NT}（体量）")
    tf = time_units_to_days1950(attrs.get(tname, {}).get("units"))
    base = {"file": ent["file"], "url": url, "mode": ent["mode"], "dt_s": dt, "span_h": ent["h1"] - ent["h0"],
            "tname": tname, "tfac": tf, "nt": nt}

    def slicer(dims, level_dim=None, level_idx=0):
        spec = []
        for d, size in dims:
            if d == tname:
                spec.append(None)
            elif d == level_dim:
                spec.append((level_idx, level_idx))
            elif size == 1:
                spec.append((0, 0))
            else:
                return None
        return spec

    uv = {}
    for v, dims in vars_.items():
        a = attrs.get(v, {})
        sn = a.get("standard_name", "")
        dnames = [d for d, _ in dims]
        if tname not in dnames or v == tname:
            continue
        if sn in RAIN_SN or sn == WSPD_SN:
            spec = slicer(dims)
            kind = "rain" if sn in RAIN_SN else "wind"
            if spec is None:
                rec["rejected"].append(f"{v}: 非时间维长度>1")
                continue
            fac = rain_factor(a.get("units")) if kind == "rain" else 1.0
            if fac is None:
                rec["rejected"].append(f"{v}: 雨量单位 {a.get('units')!r} 不可换算")
                continue
            qc = _qc_name(v, vars_, attrs)
            base_src = dict(base, kind=kind, vars=[v], qcs=[qc], grid=[v in grids], specs=[spec],
                            factor=fac, depth=None)
            rec["sources"].append(base_src)
        elif sn in (U_SN, V_SN):
            spec = slicer(dims)
            if spec is not None:
                uv[sn] = (v, spec, _qc_name(v, vars_, attrs))
        elif sn in SAL_SN:
            others = [(d, s) for d, s in dims if d != tname]
            big = [(d, s) for d, s in others if s > 1]
            if len(big) > 1:
                rec["rejected"].append(f"{v}: 多个非时间维")
                continue
            level_dim, depths = None, None
            if big:
                level_dim = big[0][0]
                if level_dim not in vars_:
                    rec["rejected"].append(f"{v}: 无深度坐标变量 {level_dim}")
                    continue
                depths = _depth_values(http, sub, url, level_dim, attrs)
            else:
                for d, _s in others:
                    if d in vars_ and _depthish(d, attrs):
                        level_dim, depths = d, _depth_values(http, sub, url, d, attrs)
                        break
                if depths is None:
                    for c in (a.get("coordinates") or "").split():
                        if c in vars_ and c != tname and _depthish(c, attrs) and len(vars_[c]) <= 1:
                            depths = _depth_values(http, sub, url, c, attrs)[:1]
                            break
                if depths is None:  # 如 IMOS SOFS：逐时刻传感器深度变量 <V>_H（positive=down），取抽样中位数
                    hv = v + "_H"
                    ha = attrs.get(hv, {})
                    if (hv in vars_ and vars_[hv] == dims and ha.get("standard_name") in ("height", "depth")
                            and (ha.get("positive") or "").lower() in ("down", "up")):
                        step = max(1, nt // 50)
                        hp = os.path.join(sub, f"{ent['file']}.{hv}.sample.ascii.gz")
                        http.get(url + ".ascii?" + q(f"{hv}[0:{step}:{nt - 1}]"), hp)
                        hs_ = [x for x in parse_ascii_values(hp, hv) if -5.0 <= x <= 6000.0]
                        if hs_:
                            dmed = statistics.median(hs_)
                            depths = [abs(dmed)]  # 盐度传感器必在水下；该类变量 height/positive 标注常自相矛盾（J5）
                            rec.setdefault("notes", []).append(f"{v} 深度取 {hv} 抽样中位数 {depths[0]:.3f} m（n={len(hs_)}）")
                if depths is None:
                    for k in ("sensor_depth", "nominal_depth", "depth"):
                        try:
                            depths = [float(a[k])]
                            break
                        except (KeyError, ValueError):
                            pass
            if not depths:
                rec["rejected"].append(f"{v}: 深度不明")
                continue
            idx1 = [j for j, d in enumerate(depths) if abs(d - S1_NOMINAL) <= DEPTH_EQ_TOL]
            if idx1:
                j1 = idx1[0]
            else:
                shallow = [(d, j) for j, d in enumerate(depths) if -DEPTH_EQ_TOL <= d <= MAX_S1_DEPTH]
                j1 = min(shallow)[1] if shallow else None
            qc = _qc_name(v, vars_, attrs)
            for kind, j in (("s1", j1), ("s5", next((j for j, d in enumerate(depths)
                                                        if abs(d - 5.0) <= DEPTH_EQ_TOL), None))):
                if j is None:
                    continue
                spec = slicer(dims, level_dim, j)
                if spec is None:
                    rec["rejected"].append(f"{v}: 切片失败")
                    continue
                rec["sources"].append(dict(base, kind=kind, vars=[v], qcs=[qc], grid=[v in grids], specs=[spec],
                                           factor=1.0, depth=round(depths[j], 3)))
            if j1 is None:
                rec["rejected"].append(f"{v}: 无 ≤{MAX_S1_DEPTH} m 层（层深 {[round(d, 2) for d in depths][:6]}…）")
    has_wspd = any(s["kind"] == "wind" for s in rec["sources"])
    if not has_wspd and U_SN in uv and V_SN in uv:
        (vu, su, qu), (vv, sv, qv) = uv[U_SN], uv[V_SN]
        rec["sources"].append(dict(base, kind="wind", vars=[vu, vv], qcs=[qu, qv], grid=[vu in grids, vv in grids],
                                   specs=[su, sv], factor=1.0, depth=None))
    for s in rec["sources"]:
        s["no_qc"] = any(x is None for x in s["qcs"])
    return rec


def discover_o(http, cache, log):
    """O 组站单（J4–J7、J13）。返回 (站列表, 发现记录, 全 CO2 数据集, 索引时段表)。"""
    allds = fetch_all_co2(http, cache)
    p1_ids = {s["erddap"] for s in p1.STATIONS}
    cands, excluded = [], []
    for d in allds:
        if d["erddap"] in p1_ids:
            excluded.append({"name": d["name"], "erddap": d["erddap"], "reason": "P1 已用 8 站之一"})
        elif d["erddap"] == CHUUK_ID or d["name"].lower().startswith("chuuk"):
            excluded.append({"name": d["name"], "erddap": d["erddap"], "reason": "Chuuk（礁湖，按设计排除）"})
        elif d["t0"] is None or d["t1"] is None:
            excluded.append({"name": d["name"], "erddap": d["erddap"], "reason": "ERDDAP 无时间覆盖"})
        else:
            cands.append(d)
    hits, spans = scan_os_index(http, cache, cands, log)
    stations, records = [], []
    for d in cands:
        ents = hits[d["name"]]
        rec = {"name": d["name"], "erddap": d["erddap"], "lat": d["lat"], "lon": d["lon"],
               "mapco2_start": p1.hour_to_iso(d["t0"]), "mapco2_end": p1.hour_to_iso(d["t1"]),
               "n_index_candidates": len(ents), "files": []}
        srcs = []
        for ent in ents:
            try:
                fr = analyze_file(http, cache, ent, log)
            except FormatError as e:
                fr = {"file": ent["file"], "path": ent["path"], "mode": ent["mode"], "skipped": str(e)}
            rec["files"].append({k: v for k, v in fr.items() if k != "sources"}
                                | {"provides": sorted({(s["kind"], s["depth"]) for s in fr.get("sources", [])},
                                                      key=lambda x: (x[0], x[1] if x[1] is not None else -1))})
            srcs += fr.get("sources", [])
        kinds = {s["kind"] for s in srcs}
        s1_depths = sorted({s["depth"] for s in srcs if s["kind"] == "s1"})
        rec["has"] = {"rain": "rain" in kinds, "s1": "s1" in kinds, "wind": "wind" in kinds, "s5": "s5" in kinds}
        rec["s1_depths_m"] = s1_depths
        rec["s1_substitute"] = any(abs(x - S1_NOMINAL) > DEPTH_EQ_TOL for x in s1_depths)
        rec["no_qc_sources"] = sorted({s["file"] for s in srcs if s["no_qc"]})
        mapco2_ok = None
        if "rain" in kinds and "s1" in kinds:
            ipath = os.path.join(cache, "erddap", d["erddap"] + "_info.csv.gz")
            http.get(f"{p1.ERDDAP}/info/{d['erddap']}/index.csv", ipath)
            vnames = {r[1] for r in csv.reader(p1.read_gz_text(ipath).splitlines()) if len(r) > 1 and r[0] == "variable"}
            mapco2_ok = {"SSS", "pCO2_sw", "time"} <= vnames
            rec["mapco2_has_SSS_pCO2_sw"] = mapco2_ok
        if "rain" in kinds and "s1" in kinds and mapco2_ok:
            rec["decision"] = "纳入"
            rec["reason"] = ("同址 ≤25 km 有 ≤1 h 雨量与 ≤1.5 m 盐度（OceanSITES，D/M）"
                             + ("；缺 1 m，以最浅 ≤1.5 m 层代「1 m」" if rec["s1_substitute"] else "")
                             + ("" if "wind" in kinds else "；无风（事件数据要求 I7 将无法满足）"))
            stations.append({"name": d["name"], "regime": "O:" + d["name"], "erddap": d["erddap"],
                             "lon": d["lon"], "lat": d["lat"], "sources": srcs})
        else:
            rec["decision"] = "排除"
            miss = [k for k, lab in (("rain", "≤1 h 雨量"), ("s1", "≤1.5 m 盐度")) if k not in kinds]
            rec["reason"] = ("同址 ≤25 km 无 OceanSITES D/M 文件" if not ents else
                             "MAPCO2 数据集缺 SSS/pCO2_sw 列" if not miss else
                             "缺 " + "、".join(lab for k, lab in (("rain", "≤1 h 雨量"), ("s1", "≤1.5 m 盐度")) if k in miss))
        records.append(rec)
    return stations, {"excluded_upfront": excluded, "candidates": records}, allds, spans


# ======================================================================== O 组加载（J8/J9）
def _constraint(src, k, i0, i1):
    v = src["vars"][k]
    name = f"{v}.{v}" if src["grid"][k] else v
    sl = "".join(f"[{i0}:1:{i1}]" if s is None else f"[{s[0]}:1:{s[1]}]" for s in src["specs"][k])
    parts = [name + sl]
    qc = src["qcs"][k]
    if qc:
        parts.append((f"{qc}.{qc}" if src["grid"][k] else qc) + sl)
    return parts


def load_source(http, cache, src, h_lo, h_hi, log, qc_drop):
    """下载一个源在 [h_lo, h_hi] 内的数据并小时化（P1 file_to_hourly）。返回 [(hour, value)]。"""
    sub = os.path.join(cache, "oceansites", "data", src["file"].rsplit(".", 1)[0])
    tpath = os.path.join(sub, "TIME.ascii.gz")
    http.get(src["url"] + ".ascii?" + q(src["tname"]), tpath)
    traw = p1.parse_dap_ascii(tpath, {src["tname"]}).get(src["tname"])
    if not traw or len(traw) != src["nt"]:
        raise FormatError(f"{src['file']} TIME 长度不符")
    fac, off = src["tfac"]
    d_lo, d_hi = hour_to_days1950(h_lo), hour_to_days1950(h_hi + 1)
    i0 = i1 = None
    for i, t in enumerate(traw):
        if d_lo <= t * fac + off < d_hi:
            if i0 is None:
                i0 = i
            i1 = i
    if i0 is None:
        return [], None
    times = array("d", (traw[i] * fac + off for i in range(i0, i1 + 1)))
    del traw
    if any(times[i] < times[i - 1] for i in range(1, len(times))):
        log.log(f"  WARN {src['file']} TIME 非单调")
    diffs = [times[i] - times[i - 1] for i in range(1, min(len(times), 20001)) if times[i] > times[i - 1]]
    dt_s = statistics.median(diffs) * 86400.0 if diffs else 3600.0
    per_hour = max(1, round(3600.0 / dt_s))
    cols = []
    for k in range(len(src["vars"])):
        parts = _constraint(src, k, i0, i1)
        tag = "_".join(src["vars"][k:k + 1]) + f".{src['kind']}.{i0}-{i1}"
        apath = os.path.join(sub, tag + ".ascii.gz")
        http.get(src["url"] + ".ascii?" + q(",".join(parts)), apath)
        want = {src["vars"][k]} | ({src["qcs"][k]} if src["qcs"][k] else set())
        arr = p1.parse_dap_ascii(apath, want)
        vals = arr.get(src["vars"][k])
        nexp = i1 - i0 + 1
        if vals is None or len(vals) != nexp:
            raise FormatError(f"{src['file']} {src['vars'][k]} 长度 {None if vals is None else len(vals)} ≠ {nexp}")
        qcs = arr.get(src["qcs"][k]) if src["qcs"][k] else None
        if src["qcs"][k] and (qcs is None or len(qcs) != nexp):
            raise FormatError(f"{src['file']} {src['qcs'][k]} 长度不符")
        cols.append((vals, qcs))
    if len(cols) == 2:  # J9：u/v 合成风速
        (u, qu), (v, qv) = cols
        vals = array("d", (math.hypot(a, b) if (a == a and b == b) else NAN for a, b in zip(u, v)))
        if qu is not None and qv is not None:
            qcs = array("d", (a if (a == a and b == b and int(a) in p1.QC_OK and int(b) in p1.QC_OK) else 4.0
                              for a, b in zip(qu, qv)))
        else:
            qcs = None
    else:
        vals, qcs = cols[0]
        if src["factor"] != 1.0:
            vals = array("d", (x * src["factor"] for x in vals))
    if qcs is None:  # J7
        qcs = array("d", [1.0]) * len(vals)
    lo, hi = RANGES[src["kind"]]
    pairs = p1.file_to_hourly(times, vals, qcs, per_hour, qc_drop.setdefault(src["kind"], {}), lo, hi)
    return pairs, dt_s


def _prio(s):
    return (MODE_RANK.get(s["mode"], 9), round(s["dt_s"]), s["span_h"], s["file"])


def process_o_station(http, st, cache, log):
    hs, he, mrows = p1.fetch_mapco2(http, st, cache, log)
    h0, n = new_axis(hs, he)
    ser = mapco2_series(mrows, h0, n)
    del mrows
    for k in ("rain", "wind", "s1", "s5"):
        ser[k] = p1.Hourly(h0, n)
    s1_depth = array("d", [NAN]) * n
    used, skipped, qc_drop = [], [], {}
    for kind in ("rain", "wind", "s1", "s5"):
        for src in sorted((s for s in st["sources"] if s["kind"] == kind), key=_prio):
            try:
                pairs, dt_s = load_source(http, cache, src, h0, h0 + n - 1, log, qc_drop)
            except FormatError as e:
                skipped.append({"file": src["file"], "kind": kind, "reason": str(e)})
                log.log(f"  跳过 {src['file']} {kind}：{e}")
                continue
            hr = ser[kind]
            added = 0
            for h, v in pairs:
                i = h - h0
                if 0 <= i < n and hr.c[i] == 0:
                    hr.s[i] = v
                    hr.c[i] = 1
                    added += 1
                    if kind == "s1":
                        s1_depth[i] = src["depth"]
            used.append({"file": src["file"], "kind": kind, "depth": src["depth"], "mode": src["mode"],
                         "native_dt_s": None if dt_s is None else round(dt_s, 1), "hours_added": added,
                         "no_qc": src["no_qc"]})
            log.log(f"  [{st['name']}] {kind} ← {src['file']}（depth={src['depth']}）+{added} h")
    series = {k: v.values() for k, v in ser.items()}
    joint = sum(1 for i in range(n) if all(series[k][i] == series[k][i] for k in ("rain", "s1", "wind", "sss05", "pco2")))
    dep_hours = {}
    for x in s1_depth:
        if x == x:
            dep_hours[str(x)] = dep_hours.get(str(x), 0) + 1
    cov = {"axis_hours": n, "mapco2_start": p1.hour_to_iso(hs), "mapco2_end": p1.hour_to_iso(he),
           "valid_hours": {k: ser[k].n_valid() for k in ser}, "hours_all5_valid": joint,
           "s1_depth_hours": dep_hours, "sources_used": used, "sources_skipped": skipped,
           "qc_dropped_by_code": {k: {str(c): x for c, x in v.items()} for k, v in qc_drop.items()}}
    del ser
    evs, stats = run_events(series, n, h0, st["name"], st["regime"], st["lon"], "O", log, s1_depth=s1_depth)
    return evs, dict(stats, coverage=cov)


# ======================================================================== R 组覆盖估计（--plan，只用索引时段）
def plan_r_coverage(http, cache, allds, spans, log):
    by_id = {d["erddap"]: d for d in allds}
    out = {}
    for st in p1.STATIONS:
        d = by_id.get(st["erddap"])
        if not d:
            out[st["name"]] = {"error": "ERDDAP 无此数据集"}
            continue
        hs, he = d["t0"], d["t1"]
        h0, n = new_axis(hs, he)
        dm = p1.list_gtmba_files(http, st, cache, hs, he, False, include_rt=False)
        allf = p1.list_gtmba_files(http, st, cache, hs, he, False, include_rt=True)
        rf = [f for f in allf if file_mode(f[0]) == "R"]
        mask = bytearray(n)
        dmk = {}
        miss = []
        for base, kind, *_ in dm:
            sp = spans.get(base.rsplit("/", 1)[-1])
            if not sp:
                miss.append(base.rsplit("/", 1)[-1])
                continue
            mk = dmk.setdefault(kind, bytearray(n))
            for i in range(max(0, sp[0] - h0), min(n, sp[1] - h0 + 1)):
                mask[i] = 1
                mk[i] = 1
        cover = {}
        for base, kind, *_ in rf:
            sp = spans.get(base.rsplit("/", 1)[-1])
            if not sp:
                miss.append(base.rsplit("/", 1)[-1])
                continue
            m = cover.setdefault(kind, bytearray(n))
            for i in range(max(0, sp[0] - h0), min(n, sp[1] - h0 + 1)):
                m[i] = 1
        lo, hi = hs - h0, he - h0 + 1
        per_kind = {k: sum(1 for i in range(lo, hi) if m[i] and not mask[i]) for k, m in cover.items()}
        kinds = ["RAIN", "SALT", "WIND"] if st["src"] == "tao" else ["TM", "SALT"]
        allk = [cover.get(k) for k in kinds]
        joint = 0 if any(m is None for m in allk) else sum(
            1 for i in range(lo, hi) if not mask[i] and all(m[i] for m in allk))
        z = bytearray(n)
        rk = kinds[0]
        pv = sum(1 for i in range(lo, hi) if cover.get(rk, z)[i] and not dmk.get(rk, z)[i]
                 and all(dmk.get(k, z)[i] or cover.get(k, z)[i] for k in kinds[1:]))
        out[st["name"]] = {"mapco2": [p1.hour_to_iso(hs), p1.hour_to_iso(he)], "p1_dm_files": len(dm),
                           "per_var_new_r_rain_hours_with_salt_wind_span_est": pv,
                           "r_files": len(rf), "p1_hours_in_mapco2": sum(mask[lo:hi]),
                           "new_r_hours_by_kind_span_est": per_kind,
                           "new_r_hours_all_kinds_span_est": joint, "files_not_in_index": miss}
    return out


# ======================================================================== 输出
def event_row(e):
    return {"group": e["group"], "station": e["station"], "regime": e["regime"], "onset_utc": e["onset_utc"],
            "local_solar_hour": e["local_solar_hour"], "season": e["season"],
            "rain_first_hour_mm": round(e["first_mm"], 3), "rain_24h_mm": round(e["acc24"], 3),
            "n_rain_hours_24h": e["n_wet"], "pre_S1_median": p1._r(e["pre_s1"], 4), "dS1_0_6h": p1._r(e["ds1"], 4),
            "pre_S05_median": p1._r(e["pre_s05"], 4), "dS05_0_6h": p1._r(e["ds05"], 4),
            "has_5m": int(e["has_5m"]), "n_pco2_in_window": e["n_pco2"], "n_ctrl": e["n_ctrl"],
            "ctrl_dS1_mean": p1._r(p1.nanmean(e["ctrl_ds1"]), 4),
            "ctrl_onsets_utc": ";".join(p1.hour_to_iso(h) for h in e["ctrl_hours"]),
            "ctrl_dS1": ";".join("" if x != x else f"{x:.4f}" for x in e["ctrl_ds1"]),
            "s1_depth_m": p1._r(e["s1_depth_m"], 3), "s1_is_substitute": int(e["s1_substitute"])}


EVENT_FIELDS = ["group", "station", "regime", "onset_utc", "local_solar_hour", "season", "rain_first_hour_mm",
                "rain_24h_mm", "n_rain_hours_24h", "pre_S1_median", "dS1_0_6h", "pre_S05_median", "dS05_0_6h",
                "has_5m", "n_pco2_in_window", "n_ctrl", "ctrl_dS1_mean", "ctrl_onsets_utc", "ctrl_dS1",
                "s1_depth_m", "s1_is_substitute"]


def jdump(obj, path):
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, allow_nan=False, default=str)
    os.replace(tmp, path)


def clean(o):
    """JSON 前把 NaN 换成 None。"""
    if isinstance(o, float):
        return None if o != o else o
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    return o


def strip_sources(stations):
    return [{k: v for k, v in s.items() if k != "sources"} | {"n_sources": len(s["sources"])} for s in stations]


def verdict_lines(g):
    d2 = g["D2"]
    return [
        f"D1 N={g['D1']['n_events_total']}，≥10事件站={g['D1']['n_stations_ge10']} → {'过' if g['D1']['pass'] else '不过'}",
        f"D2 diff={d2.get('diff_psu')} psu，CI={d2.get('ci95_cluster_bootstrap')} → {'过' if d2['pass'] else '不过'}"
        + ("" if d2.get("evaluable", False) else f"（{d2.get('reason')}）"),
        f"D4 ΔS0.5≤−0.2 事件数={g['D4']['n_events_dS05_le_-0.2']}（阈值≥{p1.D4_MIN_EVENTS}）→ {'过' if g['D4']['pass'] else '不过'}",
        f"D5 中位差={g['D5'].get('median_pre_diff_psu')}，相关={g['D5'].get('corr_event_anom')} → {'过' if g['D5']['pass'] else '不过'}",
    ]


# ======================================================================== 主流程
def run_full(args, out_dir, cache, log, http, t0):
    p1_csv = args.p1_events
    p1_evs = read_p1_events(p1_csv)
    val = validate_p1(p1_evs, os.path.join(os.path.dirname(p1_csv), "p1_summary.json"), log)
    print(f"P1 事件读入 {len(p1_evs)}，与 p1_summary 校验通过", flush=True)
    log.log("P1 校验：" + json.dumps(clean(val), ensure_ascii=False))

    p1_keys = {(e["station"], e["onset_utc"]) for e in p1_evs}
    stations_o, disc, allds, _spans = discover_o(http, cache, log)
    print(f"O 组站单：{[s['name'] for s in stations_o]}", flush=True)
    st_json = {"script": "p1b_extend.py", "version": VERSION, "mode": "full",
               "o_group": {"rule": "PMEL pmel_co2_moorings_* 去 P1 8 站与 Chuuk；同址 ≤25 km 有 ≤1 h 雨量与 ≤1.5 m 盐度"
                                   "（OceanSITES GDAC，D/M 模式）", "included": strip_sources(stations_o), **disc},
               "r_group": {}}
    jdump(clean(st_json), os.path.join(out_dir, "p1b_stations.json"))

    new_evs, per_station = [], {}
    for st in p1.STATIONS:
        print(f"[R:{st['name']}] 下载与处理中…", flush=True)
        evs, stats = process_r_station(http, st, cache, log, args.r_rule, p1_keys)
        new_evs += evs
        per_station["R:" + st["name"]] = stats
        st_json["r_group"][st["name"]] = stats["coverage"]
        print(f"[R:{st['name']}] 新增小时（五项全有效）{stats['coverage']['new_hours_all5_valid']}，"
              f"新增合格事件 {len(evs)}", flush=True)
    for st in stations_o:
        print(f"[O:{st['name']}] 下载与处理中…", flush=True)
        evs, stats = process_o_station(http, st, cache, log)
        new_evs += evs
        per_station["O:" + st["name"]] = stats
        print(f"[O:{st['name']}] 小时（五项全有效）{stats['coverage']['hours_all5_valid']}，"
              f"新增合格事件 {len(evs)}", flush=True)
    jdump(clean(st_json), os.path.join(out_dir, "p1b_stations.json"))

    dup = [(e["station"], e["onset_utc"]) for e in new_evs if (e["station"], e["onset_utc"]) in p1_keys]
    if dup:
        raise RuntimeError(f"新增事件与 P1 事件重复：{dup[:5]}")
    merged = p1_evs + new_evs
    g = gates(merged)
    g_new = gates(new_evs) if new_evs else None
    overall = all(g[k]["pass"] for k in ("D1", "D2", "D4", "D5"))
    summary = {
        "script": "p1b_extend.py", "version": VERSION, "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "runtime_s": round(time.monotonic() - t0, 1),
        "design_ref": "门同 P1；实现 I1–I17（p1_events.py）＋J1–J13（本脚本）",
        "p1_events_csv": p1_csv, "p1_validation": val, "r_rule": args.r_rule,
        "n_events": {"p1": len(p1_evs), "new_R": sum(1 for e in new_evs if e["group"] == "R"),
                     "new_O": sum(1 for e in new_evs if e["group"] == "O"), "merged": len(merged)},
        "gates_merged": g,
        "gates_verdict": {k: g[k]["pass"] for k in ("D1", "D2", "D4", "D5")},
        "p1b_pass_all": overall,
        "verdict_note": ("P1b 过：D1/D2/D4/D5 在合并集全过；进入 P2" if overall else
                         "P1b 不过：D4 不过；不再扩样"
                         if not g["D4"]["pass"] else "D4 过但 D1/D2/D5 有门不过：见 gates_merged"),
        "new_only_info": {"note": "仅信息量，不参与判定", "gates": g_new},
        "new_by_station": per_station,
        "o_group_stations": [s["name"] for s in stations_o],
        "downloads": {"http_requests": http.n_requests, "cache_hits": http.n_cached,
                      "seeded_from_p1_cache": http.n_seeded, "bytes_downloaded_gz": http.bytes},
    }
    jdump(clean(summary), os.path.join(out_dir, "p1b_summary.json"))
    with open(os.path.join(out_dir, "p1b_events_new.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=EVENT_FIELDS)
        w.writeheader()
        w.writerows(event_row(e) for e in new_evs)
    log.log(json.dumps(summary["gates_verdict"], ensure_ascii=False))
    print(f"---- P1b 合并集门判定（P1 {len(p1_evs)} ＋ 新增 {len(new_evs)}）----")
    for line in verdict_lines(g):
        print(line)
    print("P1b 总判：" + ("过" if overall else "不过") + f"；{summary['verdict_note']}")
    print(f"下载 {http.n_requests} 次请求 / {http.bytes / 1e6:.1f} MB(gz)，缓存命中 {http.n_cached}，"
          f"P1 缓存播种 {http.n_seeded}；耗时 {time.monotonic() - t0:.0f}s；产物在 {out_dir}")


def run_plan(out_dir, cache, log, http):
    stations_o, disc, allds, spans = discover_o(http, cache, log)
    r_est = plan_r_coverage(http, cache, allds, spans, log)
    st_json = {"script": "p1b_extend.py", "version": VERSION,
               "mode": "plan（仅元数据：未下数据、未数事件）",
               "o_group": {"included": strip_sources(stations_o), **disc},
               "r_group_estimate_from_index_spans": r_est}
    jdump(clean(st_json), os.path.join(out_dir, "p1b_stations.json"))
    print("---- O 组站单（plan）----")
    for c in disc["candidates"]:
        if c["n_index_candidates"] or c["decision"] == "纳入":
            print(f"{c['decision']} {c['name']}：{c['reason']}；s1 层深 {c['s1_depths_m']}；无QC文件 {len(c['no_qc_sources'])}")
    print(f"其余 {sum(1 for c in disc['candidates'] if not c['n_index_candidates'])} 站：同址无 OceanSITES D/M 文件 → 排除")
    print("---- R 组覆盖估计（按索引文件时段；小时数，不含事件）----")
    for k, v in r_est.items():
        print(f"{k}: R 文件 {v.get('r_files')}；union 口径新增小时≈{v.get('new_r_hours_all_kinds_span_est')}；"
              f"per-var 口径 R 雨新小时（盐/风有 D/M 或 R）≈{v.get('per_var_new_r_rain_hours_with_salt_wind_span_est')}")


def main():
    ap = argparse.ArgumentParser(description="P1b 一次性扩样")
    ap.add_argument("--out", help="输出目录（默认 $REPRO_OUTPUT_DIR）")
    ap.add_argument("--p1-events", default=P1_EVENTS_DEFAULT)
    ap.add_argument("--p1-cache", help="P1 缓存目录（默认 p1-events 同目录 cache/）")
    ap.add_argument("--plan", action="store_true", help="只做元数据发现与 R 组覆盖估计")
    ap.add_argument("--r-rule", choices=["union", "per-var"],
                    help="R 组「与 P1 D/M 文件时段重叠」的口径（J2/J2b，全量运行必填；参考运行用 per-var）")
    ap.add_argument("--selftest", action="store_true", help="合成数据自测")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    out_dir = args.out or os.environ.get("REPRO_OUTPUT_DIR")
    if not out_dir:
        print("需要 --out 或环境变量 REPRO_OUTPUT_DIR", file=sys.stderr)
        return 3
    if not args.plan and not args.r_rule:
        print("全量运行须指定 --r-rule union|per-var（J2/J2b；参考运行用 per-var）", file=sys.stderr)
        return 3
    os.makedirs(out_dir, exist_ok=True)
    cache = os.path.join(out_dir, "cache")
    os.makedirs(cache, exist_ok=True)
    log = p1.Log(os.path.join(out_dir, "p1b_log.txt"))
    seed = args.p1_cache or os.path.join(os.path.dirname(args.p1_events), "cache")
    http = SeededHttp(log, cache, seed)
    t0 = time.monotonic()
    log.log(f"=== start {VERSION} plan={args.plan} r_rule={args.r_rule} out={out_dir} p1_events={args.p1_events} "
            f"seed={'有' if http.seed else '无'}:{seed}", echo=True)
    try:
        if args.plan:
            run_plan(out_dir, cache, log, http)
        else:
            run_full(args, out_dir, cache, log, http, t0)
    except p1.FetchError as e:
        log.log(f"FATAL 数据源故障：{e}", echo=True)
        return 2
    except Exception:
        log.log("FATAL 未预期异常：\n" + traceback.format_exc(), echo=True)
        return 3
    log.log("=== done")
    log.close()
    return 0


# ======================================================================== 自测（合成数据，无网络）
def selftest():
    import random
    import tempfile
    ok = []

    def check(cond, msg):
        if not cond:
            raise AssertionError(msg)
        ok.append(msg)

    # 1) 解析器
    dds = """Dataset {
    Float64 TIME[TIME = 4];
    Float32 DEPTH_PSAL[DEPTH_PSAL = 3];
    Float64 DEPTH;
    Grid {
     ARRAY:
        Float32 PSAL[DEPTH_PSAL = 3][TIME = 4];
     MAPS:
        Float32 DEPTH_PSAL[DEPTH_PSAL = 3];
        Float64 TIME[TIME = 4];
    } PSAL;
    Float64 RAIN[TIME = 4];
} x;"""
    v, g = dds_parse(dds)
    check(v["PSAL"] == [("DEPTH_PSAL", 3), ("TIME", 4)] and v["DEPTH"] == [] and g == {"PSAL"}, "DDS 解析（Grid/标量/维序）")
    das = 'Attributes {\n    RAIN {\n        String units "mm/hour";\n        Float64 FillValue 1.0E35;\n    }\n' \
          '    NC_GLOBAL {\n        String summary "a { b";\n    }\n}\n'
    a = das_parse(das)
    check(a["RAIN"]["units"] == "mm/hour" and a["NC_GLOBAL"]["summary"] == "a { b", "DAS 解析")
    check(rain_factor("mm/hour") == 1.0 and rain_factor("kg m-2 s-1") == 3600.0 and rain_factor("furlong") is None
          and rain_factor("millimeters hour-1") == 1.0,
          "雨量单位换算")
    f_, o_ = time_units_to_days1950("seconds since 1970-01-01T00:00:00Z")
    check(abs(0 * f_ + o_ - 7305.0) < 1e-9 and abs(86400 * f_ - 1.0) < 1e-12, "时间单位换算")
    check(parse_ts("2011-06-13T24:00:00Z") == parse_ts("2011-06-14T00:00:00Z"), "T24 时间解析")
    td = tempfile.mkdtemp(prefix="p1b_selftest_")
    pth = os.path.join(td, "a.ascii.gz")
    with gzip.open(pth, "wt") as f:
        f.write("Dataset {\n} x;\n---------------------------------------------\nPSAL.PSAL[1][4]\n[0], 35.1, 35.2, NaN, 35.0\n\n"
                "DEPTH, 0.85\n")
    check(parse_ascii_values(pth, "PSAL")[1] == 35.2 and parse_ascii_values(pth, "DEPTH") == [0.85], "ASCII 数组/标量")
    src = {"vars": ["PSAL"], "grid": [True], "specs": [[(0, 0), None]], "qcs": ["PSAL_QC"]}
    check(_constraint(src, 0, 5, 9) == ["PSAL.PSAL[0:1:0][5:1:9]", "PSAL_QC.PSAL_QC[0:1:0][5:1:9]"], "OPeNDAP 约束")

    # 2) 合成站：R 组掩膜逻辑（J2/J3）
    rng = random.Random(1)
    n = 24 * 400
    h0 = 400000
    ser = {k: array("d", [NAN]) * n for k in ("rain", "wind", "s1", "s5", "sss05", "pco2")}
    mask = bytearray(n)
    for i in range(0, 24 * 150):
        mask[i] = 1                                   # 前 150 天＝P1 小时
    for i in range(n):
        if mask[i]:
            continue
        ser["rain"][i] = 0.0
        ser["wind"][i] = 5.0
        ser["s1"][i] = 35.0 + rng.gauss(0, 0.005)
        if i % 3 == 0:
            ser["sss05"][i] = ser["s1"][i] + 0.001
            ser["pco2"][i] = 400.0
    onsets_true = []
    for d in list(range(152, 390, 9)):
        i = d * 24 + 7
        for k in range(5):
            ser["rain"][i + k] = 4.0                  # 20 mm
        for k in range(0, 8):
            ser["s1"][i + k] -= 0.3
            if (i + k) % 3 == 0:
                ser["sss05"][i + k] -= 0.3
        onsets_true.append(i)
    # 一个紧贴 P1 边界的雨事件：前 24 h 落在 P1 小时 → 不应成事件
    i_edge = 24 * 150 + 5
    ser["rain"][i_edge] = 12.0
    class _L:
        def log(self, *a, **k):
            pass
    evs, stt = run_events(ser, n, h0, "SYN", "syn", 0.0, "R", _L(), ps_mask=prefix(mask))
    got = sorted(e["hour"] - h0 for e in evs)
    check(got == onsets_true, f"R 组事件起点＝合成真值（{len(got)} 个），边界事件被排除")
    check(all(not masked_any(prefix(mask), c - h0 - 24, c - h0 + 24) for e in evs for c in e["ctrl_hours"]),
          "对照窗不触及 P1 小时（J3）")
    check(all(e["n_ctrl"] == 5 for e in evs) and all(e["ds1"] < -0.25 for e in evs), "对照配满 5 且 ΔS1 符号正确")
    check(all(e["ds05"] <= -0.2 for e in evs), "ΔS0.5 进入 D4 计数口径")
    # 无掩膜时，边界事件应被识别（验证掩膜确实起作用）
    ser2 = {k: array("d", v) for k, v in ser.items()}
    for i in range(0, 24 * 150):
        ser2["rain"][i] = 0.0
    evs2, _ = run_events(ser2, n, h0, "SYN", "syn", 0.0, "R", _L())
    check(i_edge in [e["hour"] - h0 for e in evs2], "无掩膜时边界事件出现（掩膜有效）")

    # 3) 合并判定与 P1 csv 读取/校验
    fields = ["station", "regime", "onset_utc", "local_solar_hour", "season", "rain_first_hour_mm", "rain_24h_mm",
              "n_rain_hours_24h", "pre_S1_median", "dS1_0_6h", "pre_S05_median", "dS05_0_6h", "has_5m",
              "n_pco2_in_window", "n_ctrl", "ctrl_dS1_mean", "ctrl_onsets_utc", "ctrl_dS1"]
    cpath = os.path.join(td, "p1_events.csv")
    rows = []
    for k in range(P1_EXPECTED_N):
        stn = ["A", "B", "C", "D"][k % 4]
        ds = -0.25 if k < 14 else -0.05 + rng.gauss(0, 0.02)
        rows.append({"station": stn, "regime": "x", "onset_utc": p1.hour_to_iso(h0 + 100 * k), "local_solar_hour": 0,
                     "season": ["DJF", "MAM", "JJA", "SON"][(k // 4) % 4], "rain_first_hour_mm": 1,
                     "rain_24h_mm": 12, "n_rain_hours_24h": 3, "pre_S1_median": 35.0, "dS1_0_6h": f"{ds:.4f}",
                     "pre_S05_median": 35.001, "dS05_0_6h": f"{ds + rng.gauss(0, 0.01):.4f}" if k >= 14 else "-0.2500",
                     "has_5m": 1, "n_pco2_in_window": 8, "n_ctrl": 5 if k % 7 else 3, "ctrl_dS1_mean": 0,
                     "ctrl_onsets_utc": "", "ctrl_dS1": ";".join(f"{rng.gauss(0, 0.02):.4f}"
                                                                   for _ in range(5 if k % 7 else 3))})
    with open(cpath, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    pe = read_p1_events(cpath)
    gp = gates(pe)
    spath = os.path.join(td, "p1_summary.json")
    with open(spath, "w", encoding="utf-8") as f:
        json.dump({"gates": {"D1": gp["D1"], "D2": gp["D2"], "D5": gp["D5"],
                             "D4_info": {"n_events_dS05_le_-0.2": gp["D4"]["n_events_dS05_le_-0.2"]}}}, f)
    val = validate_p1(pe, spath, _L())
    check(val["mismatch"] == [] and gp["D4"]["n_events_dS05_le_-0.2"] == 14, "P1 csv 读取、复算与 summary 比对")
    with open(spath, "w", encoding="utf-8") as f:
        bad = dict(gp["D1"], n_events_total=203)
        json.dump({"gates": {"D1": bad, "D2": gp["D2"], "D5": gp["D5"],
                             "D4_info": {"n_events_dS05_le_-0.2": 14}}}, f)
    try:
        validate_p1(pe, spath, _L())
        check(False, "不一致应报错")
    except RuntimeError:
        ok.append("P1 summary 不一致时报错")
    gm = gates(pe + evs)
    check(gm["D4"]["n_events_dS05_le_-0.2"] == 14 + len(evs) and gm["D4"]["pass"] == (14 + len(evs) >= 30),
          f"合并 D4 计数＝P1 14＋新增 {len(evs)}")
    check(gm["D1"]["n_events_total"] == P1_EXPECTED_N + len(evs) and "SYN" in gm["D1"]["events_by_station"],
          "合并 D1 逐站计数")
    check(gm["D2"]["n_events"] == sum(1 for e in pe + evs if e["n_ctrl"] == 5), "合并 D2 只用配满 5 对照事件")
    # 4) R 组两种口径（J2 union / J2b per-var），打桩替换数据获取函数
    H0 = 24 * 15000
    ND = 400
    ev_days = {"p1": 100, "gap": 200, "new": 350}
    rain_on = {d * 24 + 7 for d in ev_days.values()}

    def fake_mapco2(http_, st_, cache_, log_, smoke_month=None):
        rows = [(h, 35.0 + (-0.3 if any(0 <= h - r < 8 for r in [H0 + x for x in rain_on]) else 0.0), 400.0)
                for h in range(H0, H0 + ND * 24, 3)]
        return H0, H0 + ND * 24 - 1, rows

    files = {"D": [("RAIN", 0, 150), ("SALT", 0, 300), ("WIND", 0, 300)],
             "R": [("RAIN", 150, 400), ("SALT", 300, 400), ("WIND", 300, 400)]}

    def fname(mode, kind, d0):
        return f"x/OS_T0N0E_DM{d0:03d}A-20100101_{mode}_{kind}_10min.nc"

    def fake_list(http_, st_, cache_, hs_, he_, smoke_, include_rt=False):
        out = [(fname("D", k, a), k, 6, H0 + a * 24) for k, a, b_ in files["D"]]
        if include_rt:
            out += [(fname("R", k, a), k, 6, H0 + a * 24) for k, a, b_ in files["R"]]
        return out

    def span_of(base):
        for mode in ("D", "R"):
            for k, a, b_ in files[mode]:
                if base == fname(mode, k, a):
                    return k, H0 + a * 24, H0 + b_ * 24 - 1
        raise KeyError(base)

    def fake_span(http_, cache_, base, log_):
        return span_of(base)[1:]

    def fake_load(http_, cache_, base, kind, per_hour, log_, qc_drop):
        k, a, b_ = span_of(base)
        if k == "RAIN":
            return {"rain": [(h, 4.0 if any(0 <= h - (H0 + r) < 5 for r in rain_on) else 0.0) for h in range(a, b_ + 1)]}
        if k == "WIND":
            return {"wind": [(h, 5.0) for h in range(a, b_ + 1)]}
        return {"s1": [(h, 35.0 - (0.3 if any(0 <= h - (H0 + r) < 8 for r in rain_on) else 0.0))
                       for h in range(a, b_ + 1)]}

    saved = (p1.fetch_mapco2, p1.list_gtmba_files, p1.load_gtmba_file, globals()["dm_file_span"])
    p1.fetch_mapco2, p1.list_gtmba_files, p1.load_gtmba_file = fake_mapco2, fake_list, fake_load
    globals()["dm_file_span"] = fake_span
    try:
        st0 = {"name": "SYNR", "regime": "syn", "lon": 0.0}
        keys = {("SYNR", p1.hour_to_iso(H0 + ev_days["p1"] * 24 + 7))}
        eu, _ = process_r_station(None, st0, td, _L(), "union", keys)
        ep, sp = process_r_station(None, st0, td, _L(), "per-var", keys)
    finally:
        p1.fetch_mapco2, p1.list_gtmba_files, p1.load_gtmba_file = saved[:3]
        globals()["dm_file_span"] = saved[3]
    du = sorted((e["hour"] - H0) // 24 for e in eu)
    dp = sorted((e["hour"] - H0) // 24 for e in ep)
    check(du == [ev_days["new"]], f"union 口径只得 D/M 全覆盖外的新事件（{du}）")
    check(dp == [ev_days["gap"], ev_days["new"]], f"per-var 口径另得 D 盐/风＋R 雨的新事件（{dp}）")
    check(sp["coverage"]["p1_reconstruction_check"]["combined_events_matching_p1"] == 1
          and sp["coverage"]["p1_reconstruction_check"]["n_combined_events_without_R_not_in_p1"] == 0,
          "per-var：P1 事件被识别为 P1、不计新增")
    shutil.rmtree(td, ignore_errors=True)
    for m in ok:
        print("PASS", m)
    print(f"selftest：{len(ok)} 项全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
