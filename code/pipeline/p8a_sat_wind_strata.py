#!/usr/bin/env python3
"""p8a_sat_wind_strata.py — P8a：按过境时刻风速分层的卫星共址检验。把 P5（RSS SMAP）与 P7a（JPL SMAP）已有的
逐事件卫星海表淡化 y 与足印 RIM-3 海面预测 x，按过境时刻系泊 U10 分档，逐档给与 P5 主结果同口径的比值和 R＝Σy/Σx 与三种
95% 区间（站×季簇 bootstrap、站簇 bootstrap、站级 jackknife-t），判强风档卫星比值落在 0.65 哪一侧：强风下 0–1 m 应已混匀、
不形成近表淡水透镜，若卫星比值仍 >0.65，卫星超出就不能用真实近表透镜解释；若 <0.65，则强风下卫星与「海面也被高估」一致。

性质：探索性；方法与判读在运行前写定；P5／P7a 全体、时间分档、雨中分组结果已看过；对「按风速分层的卫星比值」盲。
只读现有产物，不下载、不重跑 P5／P7a：P5 `p5_events.csv`、P7a `p7a_events_jpl.csv`（y、x、过境时刻）、P7c `p7c_events.csv`
（事件级 u10，敏感性用）；过境时刻风速从 P1／P1b 系泊缓存重建（p5.build → p2.build_all，只读缓存，网络调用一律拦截）。

用法（产物写 --out，未给则 $REPRO_OUTPUT_DIR）：
  python p8a_sat_wind_strata.py --selftest-synthetic   只跑合成自测（不读任何输入）
  python p8a_sat_wind_strata.py --check                核 sha256 → 合成自测 → 复现 P5／P7a 主结果 → 重建系泊 → 过境风速与 u10 连接
                                                       → 只写各档 n／簇／站计数（不算任何分档 Σy、比值）→ p8a_check.json
  python p8a_sat_wind_strata.py                        正式：--check 的全部步骤 → 各档计算 → p8a_summary.json、p8a_bins.csv、
                                                       p8a_events_wind.csv、p8a_log.txt
  --p5-dir／--p7a-dir／--p7c-events／--p1-events／--p1b-dir   输入位置（默认上游阶段产物）
  --build-dir DIR   重建用目录（其下 cache/ 即系泊缓存；默认 <out>/_build，对 P1／P1b 缓存建符号链接）。Windows 上（不便建
                    符号链接）用预先合并好的 P1＋P1b 缓存副本：先放 P1 cache、再以 P1b cache 覆盖（与 p2.LinkSeedHttp 种子
                    顺序 P1b 优先相同），重建时不再需要链接
依赖：numpy、scipy＋同目录 p5_sss_sat.py、p1_events.py、p1b_extend.py、p2_rim_test.py、p4_mech.py、raw_store.py
  （只 import，不改）。
退出码：0 跑完（无论判读）；3 输入不符（sha256、行数、复现、连接、风速缺）；4 自测不过；5 其他异常。

实现选择（X 条）：
  X1 行：P5 CSV 87 行（RSS）、P7a CSV 96 行（JPL），每行一个事件的主口径雨后过境；数值列空串＝NaN；x、y 非数的行由 ratio_block 丢弃。
  X2 过境时刻 t_s＝p5.parse_iso(post_utc)（CSV 写到整秒、截断，不跨小时）；过境小时 h＝⌊t_s/3600⌋。
  X3 过境风速 U_tr＝p5.wind_hourly(stn, h, h)：站点小时风，P2 K9（≤6 h 缺口线性插补、×p2.WIND_FACTOR、下限 0.1 m/s）——
     与 P5／P7a 在过境半步 q＝⌊t_s/1800⌋ 驱动 RIM-3 的风（小时 ⌊q/2⌋＝h）是同一个值。缺 → 该行无 U_tr；任一缺即退出 3
     （P5／P7a 的 x 需要该风有效，缺说明重建或连接有误）。
  X4 站点序列：p5.build（→ p2.build_all，K1 逐条校验 646 事件）；缓存目录＝<build-dir>/cache（默认 <out>/_build/cache，对 P1／P1b
     缓存建符号链接；Windows 上用预合并副本，见 --build-dir）。
     网络拦截：urllib.request.urlopen、http.client.HTTP(S)Connection 替换为抛 NoNetwork（BaseException 子类，不被重试捕获）。
  X5 连接：CSV 的 (station, onset_utc) 须在重建的 SMAP 期事件里（RSS 87/87、JPL 96/96），否则退出 3。
     p7c u10 按 (station, onset_utc) 连接；缺或非数的行不进 E 档（计数报出，不退出）。
  X6 分档：主＝过境风速三档，切点 CUT1＝6.031179、CUT2＝9.325534（P4f／P7c 在 646 事件事件级 u10 上的三分位）：
     TL U≤CUT1；TM CUT1<U≤CUT2；TH U>CUT2。二分（敏感性）：B8L U≤8、B8H U>8。事件级 u10（敏感性）：EL／EM／EH 同切点、
     同「≤ 归下档」（CSV 6 位小数 u10 对 6 位小数切点，同 P4f）。
  X7 区间：p5.ratio_block(rows, "y", "x")——站×季与站两种簇 p2.cluster_boot(B=10000, seed=20260926) 百分位 95%，
     pm.jackknife_station 逐站剔除 t 区间；每档在该档自己的行上重抽（同 P5 S1–S9、P7d）；过原点斜率及其站×季 CI 同函数给出。
  X8 判读 classify()：0 不可评＝Σx≥0 抽样比例 >1% 或区间不可算或 n<3；4 上端 <0（卫星未见淡化，异常）；
     1 上端 <0.65；2 下端 >0.65；3 其余。P5 可评门槛（n≥30、站×季簇≥8）只作标注「低于门槛，探索性」，不改类别（同 P7d）。
     同侧标注：站簇 CI、jackknife-t 各自按 1–3 分类，与主类别比较。
  X9 总映射 overall()：对 TH 的 (RSS, JPL) 类别组合写死 M-A…M-G；B8H、EH 只报与 TH 是否同类（敏感性）。
  X10 描述项（不进判读）：各档过境时下雨比例（IMERG 中心格 >0.1 mm/h，缺测单列）；过原点斜率；各档 U_tr 中位数、
     dt_post_h 中位数、x 均值、逐站 n；TH／B8H 剔 KEO、Papa；雨量计等效缩放 R×c_k（c＝P7c D2 的 CMORPH/雨量计 0–6 h 雨量比
     0.63493／0.7919／1.14921，区间端点同乘，只作描述）；RSS 同刻系泊 1 m 比值 Σobs1/Σrim1（ratio_block，描述）。
  X11 自测对照：P5 主块 R 1.65273（容差 0.001）、三种区间端点（容差 0.005）、簇 21／7、站 7；P7a 主块 R 1.12056、[0.4839, 1.98601]、
     [0.73029, 1.66143]、jk [0.56418, 1.67694]、簇 22／7、站 7（同容差）。切点溯源（p7c CSV 中 u10、obs1、rim1 皆有限者的
     np.quantile 1/3、2/3 与 CUT1／CUT2 差 ≤1e−5）只记录，不退出。
  X12 Windows：`raw_store` 在模块顶层 import fcntl（Unix 专有，只用于原件归档写锁）；P8a 不调用 raw_store 的任何函数，
     在 import p5_sss_sat 之前向 sys.modules 放一个空 fcntl 模块，只为让 import 通过；被 import 的文件一行不改。

Change Log：
  2026-09-27 初版。正式运行前两处修改（均在任何分档计算之前）：加 --build-dir（Windows 用预合并缓存、
     不建符号链接）与 X12（Windows 上 import raw_store 时因无 fcntl 退出 1）；--out 改为优先于 $REPRO_OUTPUT_DIR。
"""

import argparse
import csv
import hashlib
import http.client
import json
import os
import statistics
import sys
import time
import traceback
import urllib.request

if os.name == "nt" and "fcntl" not in sys.modules:      # X12：raw_store 只在原件归档 I/O 用 fcntl（Unix 专有）；P8a 不调用
    import types
    sys.modules["fcntl"] = types.ModuleType("fcntl")

import p1_events as p1
import p4_mech as pm
import p5_sss_sat as p5

VERSION = "p8a-2026-09-27a"
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
P5_DIR_DEFAULT = _rp.upstream("p5_dir")  # [repro] 以下三个默认值读上游阶段产物
P7A_DIR_DEFAULT = _rp.upstream("p7a_dir")
P7C_EVENTS_DEFAULT = _rp.upstream("p7c_events")
# [repro] 输入 sha 闸门（只在 [options] check_upstream_sha = true 时核对；缺省不核对，只记录实际 sha）：期望值缺省是参考运行的产物 sha，可在配置 [upstream_sha] 换（同 P7d）
SHA = {"p5_events.csv": _rp.expect_sha("p8a.p5_events_csv", "7a6d70758d575d180d9f0100579a9b7f6bc014350c9080ca9b5ab84ffd6985fa"),
       "p5_summary.json": _rp.expect_sha("p8a.p5_summary_json", "e5f66e65fb13d812683eee58c0f9182906f51f1cd3356247bd284bd00385ab4f"),
       "p7a_events_jpl.csv": _rp.expect_sha("p8a.p7a_events_jpl_csv", "1b529937eb13dbd02e79d7b66e22ed08810471a6f2087e30c3c559903c05f6a8"),
       "p7a_summary.json": _rp.expect_sha("p8a.p7a_summary_json", "8e030de8599390dd65600657c67ac4190346264e62d5f28bd55dbc2698453801"),
       "p7c_events.csv": _rp.expect_sha("p8a.p7c_events_csv", "247f05856bdeb73ad8261dba8009929298f5f77ea9fc66cf6901739d1b0ab4db")}
N_ROWS = {"RSS": 87, "JPL": 96}
MAIN_REF = {"RSS": {"R": 1.65273, "ss": [1.05832, 2.3908], "st": [0.65053, 2.53267], "jk": [0.32114, 2.98431],
                    "Kss": 21, "Kst": 7, "stations": 7},
            "JPL": {"R": 1.12056, "ss": [0.4839, 1.98601], "st": [0.73029, 1.66143], "jk": [0.56418, 1.67694],
                    "Kss": 22, "Kst": 7, "stations": 7}}
TOL_R, TOL_CI = 0.001, 0.005
THRESH = p5.THRESH
R1_REF, R2_REF = 0.3, 1.0
CUT1, CUT2 = 6.031179, 9.325534
CUT8 = 8.0
C_GAUGE = {"L": 0.63493, "M": 0.7919, "H": 1.14921}
MIDLAT = ("KEO", "Papa")
IMERG_RAIN = 0.1

LOG = []


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    LOG.append(line)
    print(line, flush=True)


class NoNetwork(BaseException):
    pass


def _deny(*a, **k):
    raise NoNetwork("P8a 只读缓存：拦截到网络调用")


def block_network():
    urllib.request.urlopen = _deny
    http.client.HTTPSConnection = _deny
    http.client.HTTPConnection = _deny


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def load_rows(path, product):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            d = {"product": product, "station": r["station"], "season": r["season"], "onset_utc": r["onset_utc"],
                 "post_utc": r["post_utc"], "dt_post_h": fnum(r["dt_post_h"]), "x": fnum(r["x"]), "y": fnum(r["y"]),
                 "imerg_c_post": fnum(r.get("imerg_c_post"))}
            if product == "RSS":
                d["obs1"], d["rim1"] = fnum(r.get("obs1")), fnum(r.get("rim1"))
            rows.append(d)
    return rows


def load_p7c_u10(path):
    out, allrows = {}, []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            u = fnum(r["u10"])
            out[(r["station"], r["onset_utc"])] = u
            allrows.append((u, fnum(r.get("obs1")), fnum(r.get("rim1"))))
    return out, allrows


# ======================================================================== 分档（X6）
def tbin(u):
    if not p5.isnum(u):
        return None
    return "L" if u <= CUT1 else ("M" if u <= CUT2 else "H")


def ebin(u):
    if not p5.isnum(u):
        return None
    return "L" if round(u, 6) <= CUT1 else ("M" if round(u, 6) <= CUT2 else "H")


def bin8(u):
    if not p5.isnum(u):
        return None
    return "H" if u > CUT8 else "L"


def is_w1(r):
    return p5.isnum(r["imerg_c_post"]) and r["imerg_c_post"] > IMERG_RAIN


def is_w0(r):
    return p5.isnum(r["imerg_c_post"]) and r["imerg_c_post"] <= IMERG_RAIN


# ======================================================================== 判读（X8、X9）
CAT_TEXT = {0: "不可评（分母不稳、区间不可算或事件 <3）", 1: "上端 <0.65：与 R1 相容而与 R2 不容",
            2: "下端 >0.65：与 R2 相容而与 R1 不容", 3: "不能区分", 4: "上端 <0：卫星未见淡化（异常，不作 R1 读法）"}


def cat_of(lo, hi):
    if lo is None or hi is None or not (p5.isnum(lo) and p5.isnum(hi)):
        return 0
    if hi < 0:
        return 4
    if hi < THRESH:
        return 1
    if lo > THRESH:
        return 2
    return 3


def classify(blk):
    ci = blk.get("ci95_station_season") or [None, None]
    frac = blk.get("frac_boot_den_nonneg")
    if blk.get("n", 0) < 3 or frac is None or frac > p5.DEN_NONNEG_MAX or ci[0] is None or ci[1] is None:
        main = 0
    else:
        main = cat_of(*ci)
    st = blk.get("ci95_station") or [None, None]
    jk = blk.get("jackknife_station") or {}
    jt = jk.get("t_ci95") if jk.get("evaluable") else None
    c_st = cat_of(*st)
    c_jk = cat_of(*jt) if jt else 0

    def contains(v, c):
        return None if (not c or c[0] is None or c[1] is None) else bool(c[0] <= v <= c[1])
    return {"category_code": main, "category": CAT_TEXT[main],
            "station_cluster_category": CAT_TEXT[c_st], "station_same_side": (c_st == main) if main in (1, 2, 3) else None,
            "jackknife_category": CAT_TEXT[c_jk] if jt else "jackknife 不可算",
            "jackknife_same_side": (c_jk == main) if (jt and main in (1, 2, 3)) else None,
            "ci_ss_contains_0.3": contains(R1_REF, ci), "ci_ss_contains_1": contains(R2_REF, ci),
            "ci_st_contains_0.3": contains(R1_REF, st), "ci_st_contains_1": contains(R2_REF, st),
            "jk_contains_0.3": contains(R1_REF, jt), "jk_contains_1": contains(R2_REF, jt)}


READ_TEXT = {1: "强风下卫星与幅度偏大读法（R1）一致", 2: "强风下卫星仍超出，卫星超出不能用近表透镜解释",
             3: "强风档不能区分", 0: "强风档不可评", 4: "强风档卫星未见淡化（异常）"}


def overall(c_rss, c_jpl):
    """X9 总映射：TH 档 (RSS, JPL) 类别组合。"""
    s = {c_rss, c_jpl}
    if s & {0, 4}:
        return "M-G", "至少一条链不可评或异常：逐链按实际写，不下合并结论"
    if c_rss == c_jpl == 2:
        return "M-A", "两条链强风档都仍超出（下端 >0.65）→ 卫星超出不能用近表淡水透镜解释"
    if c_rss == c_jpl == 1:
        return "M-B", "两条链强风档都 <0.65 → 强风下卫星与「海面也被高估」（R1）读法一致"
    if s == {2, 3}:
        return "M-C", "一条链强风档仍超出、另一条不能区分 → 倾向 M-A，不足以下结论"
    if s == {1, 3}:
        return "M-D", "一条链强风档 <0.65、另一条不能区分 → 倾向 M-B，不足以下结论"
    if s == {3}:
        return "M-E", "两条链强风档都不能区分"
    return "M-F", "两条链强风档类别相反（一条 <0.65、一条 >0.65）→ 矛盾，不下结论"


def block(rows, label, kind, num="y", den="x"):
    blk, _ = p5.ratio_block(rows, num, den, label=label)
    blk.pop("evaluable", None)
    blk.pop("reason", None)
    blk["p5_gate_met"] = bool(blk.get("n", 0) >= p5.MIN_N and (blk.get("clusters_station_season") or 0) >= p5.MIN_G)
    blk["kind"] = kind
    blk["judgement"] = classify(blk) if "ci95_station_season" in blk else {"category_code": 0, "category": CAT_TEXT[0]}
    return blk


def describe(rows):
    """X10 各档描述：下雨比例、中位数、逐站 n。"""
    w1, w0 = sum(is_w1(r) for r in rows), sum(is_w0(r) for r in rows)
    us = [r["u_tr"] for r in rows if p5.isnum(r.get("u_tr"))]
    dts = [r["dt_post_h"] for r in rows if p5.isnum(r["dt_post_h"])]
    xs = [r["x"] for r in rows if p5.isnum(r["x"])]
    per = {}
    for r in rows:
        per[r["station"]] = per.get(r["station"], 0) + 1
    return {"n": len(rows), "raining_at_transit_W1": w1, "not_raining_W0": w0,
            "imerg_missing": len(rows) - w1 - w0,
            "frac_raining_of_nonmissing": p5.rnd(w1 / (w1 + w0), 4) if (w1 + w0) else None,
            "median_u_tr": p5.rnd(statistics.median(us), 3) if us else None,
            "median_dt_post_h": p5.rnd(statistics.median(dts), 3) if dts else None,
            "mean_x": p5.rnd(sum(xs) / len(xs), 4) if xs else None,
            "n_by_station": dict(sorted(per.items())),
            "frac_KEO_Papa": p5.rnd(sum(per.get(s, 0) for s in MIDLAT) / len(rows), 4) if rows else None}


def gauge_scaled(blk, c):
    """X10：雨量计等效缩放（描述）。"""
    if blk.get("R") is None:
        return None
    ss = blk.get("ci95_station_season") or [None, None]
    lo, hi = (None, None) if ss[0] is None else (p5.rnd(ss[0] * c), p5.rnd(ss[1] * c))
    return {"c": c, "R_scaled": p5.rnd(blk["R"] * c), "ci95_station_season_scaled": [lo, hi],
            "category_scaled": CAT_TEXT[cat_of(lo, hi)] if lo is not None else None}


# ======================================================================== 自测
def selftest_synthetic():
    np = p5._np()
    res = []

    def chk(name, ok, detail=""):
        res.append({"name": name, "ok": bool(ok), "detail": str(detail)})
    rng = np.random.Generator(np.random.PCG64(1))
    for r_true in (0.2, 1.3):
        rows = []
        for k in range(12):
            for j in range(4):
                x = -float(rng.uniform(0.1, 1.0))
                rows.append({"station": f"S{k % 6}", "season": "DJF" if k < 6 else "JJA", "x": x, "y": r_true * x})
        b = block(rows, f"syn{r_true}", "synthetic")
        chk(f"合成比值 {r_true} 点估计", abs(b["R"] - r_true) <= 1e-9, b["R"])
        chk(f"合成比值 {r_true} 类别", b["judgement"]["category_code"] == (1 if r_true < THRESH else 2),
            b["judgement"]["category"])
    fake = lambda lo, hi, frac=0.0, n=40: {"n": n, "ci95_station_season": [lo, hi], "frac_boot_den_nonneg": frac,
                                           "ci95_station": [lo, hi], "jackknife_station": {"evaluable": False}}
    chk("判读 上端<0.65→1", classify(fake(0.1, 0.6))["category_code"] == 1)
    chk("判读 下端>0.65→2", classify(fake(0.7, 1.5))["category_code"] == 2)
    chk("判读 跨→3", classify(fake(0.2, 1.2))["category_code"] == 3)
    chk("判读 恰 0.65 上端→3", classify(fake(0.1, 0.65))["category_code"] == 3)
    chk("判读 分母不稳→0", classify(fake(0.1, 0.6, frac=0.02))["category_code"] == 0)
    chk("判读 上端<0→4", classify(fake(-2.0, -0.1))["category_code"] == 4)
    chk("判读 n<3→0", classify(fake(0.1, 0.6, n=2))["category_code"] == 0)
    chk("分档 切点归下档", tbin(CUT1) == "L" and tbin(CUT1 + 1e-9) == "M" and tbin(CUT2) == "M" and tbin(CUT2 + 1e-9) == "H")
    chk("分档 8 m/s 归下档", bin8(8.0) == "L" and bin8(8.0 + 1e-9) == "H" and tbin(float("nan")) is None)
    chk("分档 事件 u10 6 位小数", ebin(6.0311794) == "L" and ebin(6.0311796) == "M")
    exp = {(2, 2): "M-A", (1, 1): "M-B", (2, 3): "M-C", (3, 2): "M-C", (1, 3): "M-D", (3, 1): "M-D", (3, 3): "M-E",
           (1, 2): "M-F", (2, 1): "M-F", (0, 2): "M-G", (4, 3): "M-G"}
    chk("总映射", all(overall(a, b)[0] == v for (a, b), v in exp.items()))
    return res


def reproduce(rows, product):
    res = []

    def chk(name, ok, detail=""):
        res.append({"name": f"{product} {name}", "ok": bool(ok), "detail": str(detail)})
    ref = MAIN_REF[product]
    chk(f"CSV 行数 {N_ROWS[product]}", len(rows) == N_ROWS[product], len(rows))
    b = block(rows, f"{product} 全体（复现主结果）", "selftest")
    chk("复现 R", abs(b["R"] - ref["R"]) <= TOL_R, f"{b['R']} vs {ref['R']}")
    chk("复现 站×季 CI", all(abs(a - c) <= TOL_CI for a, c in zip(b["ci95_station_season"], ref["ss"])),
        f"{b['ci95_station_season']} vs {ref['ss']}")
    chk("复现 站簇 CI", all(abs(a - c) <= TOL_CI for a, c in zip(b["ci95_station"], ref["st"])),
        f"{b['ci95_station']} vs {ref['st']}")
    jt = (b.get("jackknife_station") or {}).get("t_ci95")
    chk("复现 jackknife-t", bool(jt) and all(abs(a - c) <= TOL_CI for a, c in zip(jt, ref["jk"])), f"{jt} vs {ref['jk']}")
    chk("复现 簇／站", b["clusters_station_season"] == ref["Kss"] and b["clusters_station"] == ref["Kst"]
        and b["stations"] == ref["stations"], f"{b['clusters_station_season']}/{b['clusters_station']}/{b['stations']}")
    return res, b


# ======================================================================== 输入与风速
def load_inputs(args):
    paths = {"p5_events.csv": os.path.join(args.p5_dir, "p5_events.csv"),
             "p5_summary.json": os.path.join(args.p5_dir, "p5_summary.json"),
             "p7a_events_jpl.csv": os.path.join(args.p7a_dir, "p7a_events_jpl.csv"),
             "p7a_summary.json": os.path.join(args.p7a_dir, "p7a_summary.json"),
             "p7c_events.csv": args.p7c_events}
    shas = {k: sha256_file(p) for k, p in paths.items()}
    for k, v in shas.items():
        log(f"输入 {k} sha256 {v}{'' if v == SHA[k] else '  ≠ 预期值'}")
    ok = not any(_rp.upstream_sha_mismatch(shas[k], SHA[k]) for k in SHA)
    rows = {"RSS": load_rows(paths["p5_events.csv"], "RSS"), "JPL": load_rows(paths["p7a_events_jpl.csv"], "JPL")}
    u10, p7c_all = load_p7c_u10(paths["p7c_events.csv"])
    return ok, shas, paths, rows, u10, p7c_all


def attach_wind(args, out_dir, rows, u10):
    """X2–X5：重建系泊站点序列（只读缓存），给每行过境风速与事件 u10。返回问题列表。"""
    block_network()
    bdir = args.build_dir or os.path.join(out_dir, "_build")
    os.makedirs(bdir, exist_ok=True)
    blog = p1.Log(os.path.join(out_dir, "p8a_build_log.txt"))
    ns = argparse.Namespace(p1_events=args.p1_events, p1b_dir=args.p1b_dir)
    sb, events, info = p5.build(ns, bdir, blog)
    blog.close()
    log(f"重建：646 事件 {info['events_646']}、SMAP 期 {info['events_smap_era']}、D4 {info['rebuild_d4']}")
    keys = {(e["_st"], e["onset_utc"]) for e in events}
    problems = []
    if info["events_646"] != 646 or info["events_smap_era"] != 317:
        problems.append(f"重建事件数 {info['events_646']}/{info['events_smap_era']}（应 646/317）")
    for prod, rr in rows.items():
        miss_k, miss_u, miss_e = 0, 0, 0
        for r in rr:
            if (r["station"], r["onset_utc"]) not in keys or r["station"] not in sb:
                miss_k += 1
                r["u_tr"] = float("nan")
            else:
                h = int(p5.parse_iso(r["post_utc"]) // 3600)
                w = p5.wind_hourly(sb[r["station"]], h, h)
                r["u_tr"] = float(w[0]) if w is not None else float("nan")
                miss_u += w is None
            r["u10_event"] = u10.get((r["station"], r["onset_utc"]), float("nan"))
            miss_e += not p5.isnum(r["u10_event"])
            r["bin_t"], r["bin_8"], r["bin_e"] = tbin(r["u_tr"]), bin8(r["u_tr"]), ebin(r["u10_event"])
        log(f"{prod}：连接缺 {miss_k}、过境风速缺 {miss_u}、事件 u10 缺 {miss_e}（共 {len(rr)} 行）")
        if miss_k:
            problems.append(f"{prod} 有 {miss_k} 行不在重建事件里")
        if miss_u:
            problems.append(f"{prod} 有 {miss_u} 行过境风速缺")
    return problems


def cut_provenance(p7c_all):
    np = p5._np()
    v = np.array([u for u, o, m in p7c_all if p5.isnum(u) and p5.isnum(o) and p5.isnum(m)], float)
    q1, q2 = np.quantile(v, [1 / 3, 2 / 3])
    return {"n": int(len(v)), "q1": float(q1), "q2": float(q2),
            "ok": bool(abs(q1 - CUT1) <= 1e-5 and abs(q2 - CUT2) <= 1e-5)}


def counts(rr):
    out = {}
    for key, lab in (("bin_t", "T"), ("bin_8", "B8"), ("bin_e", "E")):
        for b in ("L", "M", "H"):
            sub = [r for r in rr if r[key] == b]
            if not sub and b == "M" and key == "bin_8":
                continue
            out[f"{lab}{b}"] = {"n": len(sub), "clusters_station_season": len({(r["station"], r["season"]) for r in sub}),
                                "stations": len({r["station"] for r in sub})}
    return out


# ======================================================================== 主体
def prepare(args, out_dir):
    ok, shas, paths, rows, u10, p7c_all = load_inputs(args)
    if not ok:
        log("输入 sha256 与预期不符，退出 3")
        return 3, None
    st = selftest_synthetic()
    reps = {}
    for prod in ("RSS", "JPL"):
        r, b = reproduce(rows[prod], prod)
        st += r
        reps[prod] = b
    for s in st:
        log(f"自测 {'OK ' if s['ok'] else 'FAIL'} {s['name']} {s['detail']}")
    if not all(s["ok"] for s in st):
        rc = 4 if any((not s["ok"]) and s["name"].startswith(("合成", "判读", "分档", "总映射")) for s in st) else 3
        with open(os.path.join(out_dir, "p8a_selftest.json"), "w") as f:
            json.dump(st, f, ensure_ascii=False, indent=1)
        log(f"自测不过，退出 {rc}")
        return rc, None
    problems = attach_wind(args, out_dir, rows, u10)
    prov = cut_provenance(p7c_all)
    log(f"切点溯源（信息）：p7c n={prov['n']} q1={prov['q1']:.6f} q2={prov['q2']:.6f} ok={prov['ok']}")
    if problems:
        for p in problems:
            log(f"输入问题：{p}")
        log("连接或风速不符，退出 3")
        return 3, None
    return 0, {"shas": shas, "paths": paths, "rows": rows, "selftest": st, "reproduced": reps, "cut_provenance": prov}


def run_check(args, out_dir):
    rc, ctx = prepare(args, out_dir)
    if rc:
        return rc
    out = {"script": os.path.basename(__file__), "version": VERSION, "mode": "check", "inputs": ctx["shas"],
           "selftest": ctx["selftest"], "cut_provenance": ctx["cut_provenance"],
           "counts_only": {p: counts(rr) for p, rr in ctx["rows"].items()},
           "note": "只含各档 n／簇／站计数；未算任何分档 Σy 或比值"}
    with open(os.path.join(out_dir, "p8a_check.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    for p, c in out["counts_only"].items():
        log(f"{p} 计数：" + "；".join(f"{k} n={v['n']} K={v['clusters_station_season']} 站={v['stations']}" for k, v in c.items()))
    return 0


def run_full(args, out_dir):
    rc, ctx = prepare(args, out_dir)
    if rc:
        return rc
    here = os.path.dirname(os.path.abspath(__file__))
    out = {"script": os.path.basename(__file__), "version": VERSION,
           "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "code_sha256": {n: sha256_file(os.path.join(here, n)) for n in
                           ("p8a_sat_wind_strata.py", "p5_sss_sat.py", "raw_store.py", "p1_events.py", "p1b_extend.py",
                            "p2_rim_test.py", "p4_mech.py")},
           "inputs": {"paths": ctx["paths"], "sha256": ctx["shas"]},
           "boot": {"B": p5.B, "seed": p5.SEED, "threshold": THRESH, "den_nonneg_max": p5.DEN_NONNEG_MAX},
           "cuts": {"CUT1": CUT1, "CUT2": CUT2, "CUT8": CUT8}, "cut_provenance": ctx["cut_provenance"],
           "selftest": ctx["selftest"], "reproduced": ctx["reproduced"],
           "main": {}, "sensitivity": {}, "descriptive": {}, "reading_TH": {}, "overall": {}}
    for prod in ("RSS", "JPL"):
        rr = ctx["rows"][prod]
        m, s, d = {}, {}, {}
        for b in ("L", "M", "H"):
            sub = [r for r in rr if r["bin_t"] == b]
            m[f"T{b}"] = block(sub, f"{prod} 过境风速 {b} 档", "main" if b == "H" else "reported")
            d[f"T{b}"] = describe(sub)
            d[f"T{b}_gauge_scaled"] = gauge_scaled(m[f"T{b}"], C_GAUGE[b])
            esub = [r for r in rr if r["bin_e"] == b]
            s[f"E{b}"] = block(esub, f"{prod} 事件级 u10 {b} 档", "sensitivity")
            d[f"E{b}"] = describe(esub)
        for b in ("L", "H"):
            sub = [r for r in rr if r["bin_8"] == b]
            s[f"B8{b}"] = block(sub, f"{prod} 过境风速 {'>8' if b == 'H' else '≤8'} m/s", "sensitivity")
            d[f"B8{b}"] = describe(sub)
        d["TH_noKEOPapa"] = block([r for r in rr if r["bin_t"] == "H" and r["station"] not in MIDLAT],
                                  f"{prod} TH 剔 KEO/Papa", "descriptive")
        d["B8H_noKEOPapa"] = block([r for r in rr if r["bin_8"] == "H" and r["station"] not in MIDLAT],
                                   f"{prod} >8 m/s 剔 KEO/Papa", "descriptive")
        if prod == "RSS":
            for key, lab in (("bin_t", "T"), ("bin_8", "B8")):
                for b in (("L", "M", "H") if key == "bin_t" else ("L", "H")):
                    sub = [r for r in rr if r[key] == b]
                    d[f"{lab}{b}_R1m_matched"] = block(sub, f"RSS {lab}{b} 同刻系泊 1 m Σobs1/Σrim1", "descriptive",
                                                        num="obs1", den="rim1")
        cth = m["TH"]["judgement"]["category_code"]
        s["same_side_as_TH"] = {k: s[k]["judgement"]["category_code"] == cth for k in ("B8H", "EH")}
        out["main"][prod], out["sensitivity"][prod], out["descriptive"][prod] = m, s, d
        out["reading_TH"][prod] = READ_TEXT[cth]
    code, text = overall(out["main"]["RSS"]["TH"]["judgement"]["category_code"],
                         out["main"]["JPL"]["TH"]["judgement"]["category_code"])
    out["overall"] = {"code": code, "text": text,
                      "sensitivity_B8H": {p: READ_TEXT[out["sensitivity"][p]["B8H"]["judgement"]["category_code"]]
                                          for p in ("RSS", "JPL")},
                      "sensitivity_EH": {p: READ_TEXT[out["sensitivity"][p]["EH"]["judgement"]["category_code"]]
                                         for p in ("RSS", "JPL")}}

    with open(os.path.join(out_dir, "p8a_summary.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    cols = ["product", "key", "kind", "label", "n", "clusters_station_season", "clusters_station", "stations", "sum_num",
            "sum_den", "R", "ss_lo", "ss_hi", "st_lo", "st_hi", "jk_lo", "jk_hi", "loo_lo", "loo_hi", "frac_den_nonneg",
            "slope", "slope_lo", "slope_hi", "p5_gate_met", "category", "station_same_side", "jackknife_same_side"]
    with open(os.path.join(out_dir, "p8a_bins.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for prod in ("RSS", "JPL"):
            items = (list(out["main"][prod].items()) +
                     [(k, v) for k, v in out["sensitivity"][prod].items() if isinstance(v, dict) and "judgement" in v] +
                     [(k, v) for k, v in out["descriptive"][prod].items() if isinstance(v, dict) and "judgement" in v])
            for k, b in items:
                jk = b.get("jackknife_station") or {}
                jt = jk.get("t_ci95") or [None, None]
                loo = jk.get("loo_range") or [None, None]
                ss = b.get("ci95_station_season") or [None, None]
                st = b.get("ci95_station") or [None, None]
                sl = b.get("slope_ci95_station_season") or [None, None]
                j = b["judgement"]
                w.writerow([prod, k, b["kind"], b["label"], b.get("n"), b.get("clusters_station_season"),
                            b.get("clusters_station"), b.get("stations"), b.get("sum_num"), b.get("sum_den"), b.get("R"),
                            ss[0], ss[1], st[0], st[1], jt[0], jt[1], loo[0], loo[1], b.get("frac_boot_den_nonneg"),
                            b.get("slope_through_origin"), sl[0], sl[1], b["p5_gate_met"], j["category"],
                            j.get("station_same_side"), j.get("jackknife_same_side")])
                log(f"{prod} {k:<18} n={b.get('n')} K={b.get('clusters_station_season')}/{b.get('clusters_station')} "
                    f"R={b.get('R')} ss={ss} st={st} jk={jt} frac={b.get('frac_boot_den_nonneg')} → {j['category']}")
    with open(os.path.join(out_dir, "p8a_events_wind.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["product", "station", "season", "onset_utc", "post_utc", "u_tr", "u10_event", "bin_t", "bin_8", "bin_e"])
        for prod in ("RSS", "JPL"):
            for r in ctx["rows"][prod]:
                w.writerow([prod, r["station"], r["season"], r["onset_utc"], r["post_utc"], p5.rnd(r["u_tr"], 6),
                            p5.rnd(r["u10_event"], 6), r["bin_t"], r["bin_8"], r["bin_e"]])
    log(f"总映射：{code} {text}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="P8a 按风速分层的卫星共址检验")
    ap.add_argument("--selftest-synthetic", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--p5-dir", default=P5_DIR_DEFAULT)
    ap.add_argument("--p7a-dir", default=P7A_DIR_DEFAULT)
    ap.add_argument("--p7c-events", default=P7C_EVENTS_DEFAULT)
    ap.add_argument("--p1-events", default=p5.p2.P1_EVENTS_DEFAULT)
    ap.add_argument("--p1b-dir", default=p5.p2.P1B_DIR_DEFAULT)
    ap.add_argument("--build-dir", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    if args.selftest_synthetic:
        res = selftest_synthetic()
        for s in res:
            print(("OK  " if s["ok"] else "FAIL"), s["name"], s["detail"])
        return 0 if all(s["ok"] for s in res) else 4
    out_dir = args.out or os.environ.get("REPRO_OUTPUT_DIR")        # --out 优先（同 P5）
    if not out_dir:
        print("需要 $REPRO_OUTPUT_DIR 或 --out", file=sys.stderr)
        return 5
    os.makedirs(out_dir, exist_ok=True)
    try:
        log(f"{VERSION} start mode={'check' if args.check else 'full'}")
        rc = run_check(args, out_dir) if args.check else run_full(args, out_dir)
    except NoNetwork as e:
        log(f"缓存不全、拦截到网络调用（退出 3）：{e}")
        rc = 3
    except Exception:
        log("异常：" + traceback.format_exc())
        rc = 5
    log(f"done rc={rc}")
    with open(os.path.join(out_dir, "p8a_log.txt" if not args.check else "p8a_check_log.txt"), "w") as f:
        f.write("\n".join(LOG) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
