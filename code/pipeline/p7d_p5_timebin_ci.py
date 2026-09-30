#!/usr/bin/env python3
"""p7d_p5_timebin_ci.py — P7d：给 P5 卫星共址结果的雨后过境时间分档（T1–T4）与「过境时下雨」分组（W1／W0）
补上与 P5 主结果同口径的三种 95% 区间（站×季簇 bootstrap、站簇 bootstrap、站级 jackknife-t），判断各档能否区分 R1（约 0.3）
与 R2（约 1）。探索性、非盲（分档点估计在 P5 结果里已看过）；方法与判读在运行前写定。
只读 P5 已有产物 `p5_events.csv`（不下载、不重跑 P5）；区间直接调用 `p5_sss_sat.ratio_block`（与 P5 主结果同一函数、同一 seed）。

用法：
  python p7d_p5_timebin_ci.py --selftest-synthetic     只跑合成自测（不读 CSV）
  python p7d_p5_timebin_ci.py [--p5-dir DIR] [--out DIR]
      正式：核 CSV sha256 → 自测 1–4 → 各档计算 → p7d_summary.json、p7d_bins.csv、p7d_log.txt
      产物写 $REPRO_OUTPUT_DIR（未设则 --out）。
  --p5-dir  P5 产物目录（默认 p5-sss 阶段输出目录）
依赖：numpy、scipy＋同目录 p5_sss_sat.py（及其 import 的 raw_store／p1_events／p1b_extend／p2_rim_test／p4_mech），只 import 不改。
退出码：0 跑完；3 输入不符（sha256、行数、分档复现）；4 自测不过；5 其他异常。

实现选择（W 条）：
  W1 行＝CSV 一行（一个事件的主口径雨后过境）；数值列空串＝NaN；x、y 任一非数的行由 ratio_block 丢弃（87 行实际全有）。
  W2 分档：T1–T4 按 dt_post_h 左闭右开 [0,3)/[3,6)/[6,9)/[9,12)；W1＝imerg_c_post>0.1，W0＝imerg_c_post≤0.1，缺测不进 W。
  W3 区间：p5.ratio_block(rows, "y", "x")——站×季与站两种簇 p2.cluster_boot(B=10000, seed=20260926) 百分位 95%，
     pm.jackknife_station 逐站剔除 t 区间；每档在该档自己的行上重抽（与 P5 S1–S9 同）。
  W4 判读 classify()：类别 0–3；阈值 0.65＝p5.THRESH；分母不稳＝frac_boot_den_nonneg>p5.DEN_NONNEG_MAX（0.01）。
     同侧标注：站簇 CI、jackknife-t 区间各自按同一规则分类，与主类别比较。
  W5 D1–D4 描述项；D2 剔除集＝{KEO, Papa}（同 P5 S8）。
  W6 自测对照：p5_summary.json 主块（R 容差 0.001，区间端点 0.005）；P5 结果的分档 n 与点估计（容差 0.0051）。

Change Log：
  2026-09-27 初版。
"""

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
import traceback

import p5_sss_sat as p5

VERSION = "p7d-2026-09-27a"
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
P5_DIR_DEFAULT = _rp.upstream("p5_dir")  # [repro] 读 p5-sss 阶段输出目录
# [repro] 输入 sha 闸门（只在 [options] check_upstream_sha = true 时核对；缺省不核对，只记录实际 sha）：期望值缺省是参考运行的产物 sha，可在配置 [upstream_sha] 换
EVENTS_SHA = _rp.expect_sha("p7d.p5_events_csv", "7a6d70758d575d180d9f0100579a9b7f6bc014350c9080ca9b5ab84ffd6985fa")
SUMMARY_SHA = _rp.expect_sha("p7d.p5_summary_json", "e5f66e65fb13d812683eee58c0f9182906f51f1cd3356247bd284bd00385ab4f")
N_EVENTS = 87
THRESH = p5.THRESH
R1_REF, R2_REF = 0.3, 1.0
MIDLAT = ("KEO", "Papa")
IMERG_RAIN = 0.1

TBINS = [("T1", 0.0, 3.0), ("T2", 3.0, 6.0), ("T3", 6.0, 9.0), ("T4", 9.0, 12.0)]
# P5 结果第六节已报（n, R）；W0 报的是 (n, Σx, Σy)
P5_REPORTED = {"T1": (16, 3.52), "T2": (31, 2.51), "T3": (23, 0.80), "T4": (17, 0.87), "W1": (58, 1.98)}
P5_REPORTED_W0 = (25, -1.71, 6.45)
P5_PRIMARY = {"R": 1.65273, "ci95_station_season": [1.05832, 2.3908], "ci95_station": [0.65053, 2.53267],
              "t_ci95": [0.32114, 2.98431], "Kss": 21, "Kst": 7, "stations": 7}
TOL_R, TOL_CI, TOL_BIN = 0.001, 0.005, 0.0051

LOG = []


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    LOG.append(line)
    print(line, flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fnum(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return float("nan")
    return v


def load_rows(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append({"station": r["station"], "season": r["season"], "onset_utc": r["onset_utc"],
                         "dt_post_h": fnum(r["dt_post_h"]), "x": fnum(r["x"]), "y": fnum(r["y"]),
                         "imerg_c_post": fnum(r["imerg_c_post"])})
    return rows


def in_tbin(r, lo, hi):
    return p5.isnum(r["dt_post_h"]) and lo <= r["dt_post_h"] < hi


def is_w1(r):
    return p5.isnum(r["imerg_c_post"]) and r["imerg_c_post"] > IMERG_RAIN


def is_w0(r):
    return p5.isnum(r["imerg_c_post"]) and r["imerg_c_post"] <= IMERG_RAIN


def cat_of(lo, hi):
    """W4 类别 1–3（不含分母条件）。"""
    if lo is None or hi is None or not (p5.isnum(lo) and p5.isnum(hi)):
        return 0
    if hi < THRESH:
        return 1
    if lo > THRESH:
        return 2
    return 3


CAT_TEXT = {0: "不可评（分母不稳或区间不可算）", 1: "与 R1 相容而与 R2 不容", 2: "与 R2 相容而与 R1 不容", 3: "该档不能区分"}


def classify(blk):
    """W4：0 分母不稳 → 1 上端<0.65 → 2 下端>0.65 → 3 其余；另加同侧标注。"""
    ci = blk.get("ci95_station_season") or [None, None]
    frac = blk.get("frac_boot_den_nonneg")
    if frac is None or frac > p5.DEN_NONNEG_MAX or ci[0] is None or ci[1] is None:
        main = 0
    else:
        main = cat_of(*ci)
    st = blk.get("ci95_station") or [None, None]
    jk = blk.get("jackknife_station") or {}
    jt = jk.get("t_ci95") if jk.get("evaluable") else None
    c_st = cat_of(*st)
    c_jk = cat_of(*jt) if jt else 0
    def contains(v, c):
        return None if c[0] is None or c[1] is None else bool(c[0] <= v <= c[1])
    return {"category_code": main, "category": CAT_TEXT[main],
            "station_cluster_category": CAT_TEXT[c_st], "station_same_side": (c_st == main) if main in (1, 2, 3) else None,
            "jackknife_category": CAT_TEXT[c_jk] if jt else "jackknife 不可算",
            "jackknife_same_side": (c_jk == main) if (jt and main in (1, 2, 3)) else None,
            "ci_ss_contains_0.3": contains(R1_REF, ci), "ci_ss_contains_1": contains(R2_REF, ci),
            "ci_st_contains_0.3": contains(R1_REF, st), "ci_st_contains_1": contains(R2_REF, st),
            "jk_contains_0.3": contains(R1_REF, jt) if jt else None, "jk_contains_1": contains(R2_REF, jt) if jt else None}


def block(rows, label, kind):
    blk, _ = p5.ratio_block(rows, "y", "x", label=label)
    blk.pop("evaluable", None)
    blk.pop("reason", None)
    blk["p5_gate_met"] = bool(blk.get("n", 0) >= p5.MIN_N and (blk.get("clusters_station_season") or 0) >= p5.MIN_G)
    blk["kind"] = kind
    if blk.get("n", 0) >= 3 and "ci95_station_season" in blk:
        blk["judgement"] = classify(blk)
    else:
        blk["judgement"] = {"category_code": 0, "category": CAT_TEXT[0] + "（事件 <3）"}
    return blk


def fmt_ci(c):
    return "—" if not c or c[0] is None else f"[{c[0]:.3f}, {c[1]:.3f}]"


# ======================================================================== 自测
def selftest_synthetic():
    np = p5._np()
    res = []
    def chk(name, ok, detail=""):
        res.append({"name": name, "ok": bool(ok), "detail": str(detail)})
    rng = np.random.Generator(np.random.PCG64(1))
    for r_true in (0.3, 1.0):
        rows = []
        for k in range(12):
            for j in range(4):
                x = -float(rng.uniform(0.1, 1.0))
                rows.append({"station": f"S{k % 6}", "season": "DJF" if k < 6 else "JJA", "x": x, "y": r_true * x})
        b = block(rows, f"syn{r_true}", "synthetic")
        chk(f"合成比值 {r_true} 点估计", abs(b["R"] - r_true) <= 1e-9, b["R"])
        exp = 1 if r_true < THRESH else 2
        chk(f"合成比值 {r_true} 类别", b["judgement"]["category_code"] == exp, b["judgement"]["category"])
    fake = lambda lo, hi, frac=0.0: {"ci95_station_season": [lo, hi], "frac_boot_den_nonneg": frac,
                                    "ci95_station": [lo, hi], "jackknife_station": {"evaluable": False}}
    chk("判读 上端<0.65→1", classify(fake(0.1, 0.6))["category_code"] == 1)
    chk("判读 下端>0.65→2", classify(fake(0.7, 1.5))["category_code"] == 2)
    chk("判读 跨→3", classify(fake(0.2, 1.2))["category_code"] == 3)
    chk("判读 分母不稳→0", classify(fake(0.1, 0.6, frac=0.02))["category_code"] == 0)
    chk("判读 恰 0.65 上端→3", classify(fake(0.1, 0.65))["category_code"] == 3)
    return res


def selftest_real(rows, summ):
    res = []
    def chk(name, ok, detail=""):
        res.append({"name": name, "ok": bool(ok), "detail": str(detail)})
    chk("CSV 行数 87", len(rows) == N_EVENTS, len(rows))
    b = block(rows, "全 87（复现 P5 主结果）", "selftest")
    pp = summ["primary"]
    chk("R 对 p5_summary", abs(b["R"] - pp["R"]) <= TOL_R, f"{b['R']} vs {pp['R']}")
    chk("R 对参照值 1.65273", abs(b["R"] - P5_PRIMARY["R"]) <= TOL_R, b["R"])
    for key in ("ci95_station_season", "ci95_station"):
        ok = all(abs(a - c) <= TOL_CI for a, c in zip(b[key], pp[key]))
        chk(f"{key} 对 p5_summary", ok, f"{b[key]} vs {pp[key]}")
        ok2 = all(abs(a - c) <= TOL_CI for a, c in zip(b[key], P5_PRIMARY[key]))
        chk(f"{key} 对参照值", ok2, f"{b[key]} vs {P5_PRIMARY[key]}")
    jt = b["jackknife_station"].get("t_ci95")
    chk("jackknife-t 对 p5_summary", jt and all(abs(a - c) <= TOL_CI for a, c in
                                                zip(jt, pp["jackknife_station"]["t_ci95"])), f"{jt}")
    chk("簇数／站数", b["clusters_station_season"] == P5_PRIMARY["Kss"] and b["clusters_station"] == P5_PRIMARY["Kst"]
        and b["stations"] == P5_PRIMARY["stations"], f"{b['clusters_station_season']}/{b['clusters_station']}/{b['stations']}")
    for name, lo, hi in TBINS:
        sub = [r for r in rows if in_tbin(r, lo, hi)]
        R = sum(r["y"] for r in sub) / sum(r["x"] for r in sub)
        n0, R0 = P5_REPORTED[name]
        chk(f"{name} n 与点估计复现", len(sub) == n0 and abs(R - R0) <= TOL_BIN, f"n {len(sub)} R {R:.4f}")
    w1 = [r for r in rows if is_w1(r)]
    R = sum(r["y"] for r in w1) / sum(r["x"] for r in w1)
    chk("W1 n 与点估计复现", len(w1) == P5_REPORTED["W1"][0] and abs(R - P5_REPORTED["W1"][1]) <= TOL_BIN, f"n {len(w1)} R {R:.4f}")
    w0 = [r for r in rows if is_w0(r)]
    sx, sy = sum(r["x"] for r in w0), sum(r["y"] for r in w0)
    n0, sx0, sy0 = P5_REPORTED_W0
    chk("W0 n、Σx、Σy 复现", len(w0) == n0 and abs(sx - sx0) <= TOL_BIN and abs(sy - sy0) <= TOL_BIN,
        f"n {len(w0)} Σx {sx:.4f} Σy {sy:.4f}")
    return res, b


# ======================================================================== 主体
def run(args, out_dir):
    ev_path = os.path.join(args.p5_dir, "p5_events.csv")
    sm_path = os.path.join(args.p5_dir, "p5_summary.json")
    sha_ev, sha_sm = sha256_file(ev_path), sha256_file(sm_path)
    log(f"p5_events.csv sha256 {sha_ev}；p5_summary.json sha256 {sha_sm}")
    if _rp.upstream_sha_mismatch(sha_ev, EVENTS_SHA) or _rp.upstream_sha_mismatch(sha_sm, SUMMARY_SHA):
        log("输入 sha256 与预期不符，退出 3")
        return 3
    rows = load_rows(ev_path)
    with open(sm_path) as f:
        summ = json.load(f)
    syn = selftest_synthetic()
    real, prim = selftest_real(rows, summ)
    st_all = syn + real
    for s in st_all:
        log(f"自测 {'OK ' if s['ok'] else 'FAIL'} {s['name']} {s['detail']}")
    if not all(s["ok"] for s in st_all):
        input_like = ("复现", "行数")
        rc = 3 if all(s["ok"] or any(t in s["name"] for t in input_like) for s in st_all) else 4
        with open(os.path.join(out_dir, "p7d_selftest.json"), "w") as f:
            json.dump(st_all, f, ensure_ascii=False, indent=1)
        log(f"自测不过，退出 {rc}")
        return rc

    out = {"script": os.path.basename(__file__), "version": VERSION,
           "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "code_sha256": {n: sha256_file(os.path.join(os.path.dirname(os.path.abspath(__file__)), n)) for n in
                           ("p7d_p5_timebin_ci.py", "p5_sss_sat.py", "raw_store.py", "p1_events.py", "p1b_extend.py",
                            "p2_rim_test.py", "p4_mech.py")},
           "inputs": {"p5_dir": args.p5_dir, "p5_events_sha256": sha_ev, "p5_summary_sha256": sha_sm},
           "boot": {"B": p5.B, "seed": p5.SEED, "threshold": THRESH, "den_nonneg_max": p5.DEN_NONNEG_MAX},
           "selftest": st_all, "primary_reproduced": prim, "main": {}, "descriptive": {}}

    # 主判：T1–T4、W1、W0
    for name, lo, hi in TBINS:
        out["main"][name] = block([r for r in rows if in_tbin(r, lo, hi)], f"{name} 雨后过境 [{lo:g},{hi:g}) h", "main")
    out["main"]["W1"] = block([r for r in rows if is_w1(r)], "W1 过境中心格 IMERG>0.1 mm/h", "main")
    out["main"]["W0"] = block([r for r in rows if is_w0(r)], "W0 过境中心格 IMERG≤0.1 mm/h", "main")

    # 描述 D1–D4
    d = out["descriptive"]
    t34 = [r for r in rows if in_tbin(r, 6.0, 12.0)]
    d["D1_6to12"] = block(t34, "D1 [6,12) h 合并", "descriptive")
    for name, lo, hi in TBINS + [("D1", 6.0, 12.0)]:
        sub = [r for r in rows if in_tbin(r, lo, hi) and r["station"] not in MIDLAT]
        d[f"D2_{name}_noKEOPapa"] = block(sub, f"D2 {name} [{lo:g},{hi:g}) h 剔 KEO/Papa", "descriptive")
    xt = {}
    for name, lo, hi in TBINS:
        sub = [r for r in rows if in_tbin(r, lo, hi)]
        xt[name] = {"W1": sum(is_w1(r) for r in sub), "W0": sum(is_w0(r) for r in sub),
                    "IMERG_missing": sum(not p5.isnum(r["imerg_c_post"]) for r in sub)}
    d["D3_crosstab"] = xt
    d["D3_6to12_W0"] = block([r for r in t34 if is_w0(r)], "D3 [6,12) h 且 IMERG≤0.1", "descriptive")
    d["D3_6to12_W1"] = block([r for r in t34 if is_w1(r)], "D3 [6,12) h 且 IMERG>0.1（对照）", "descriptive")
    per = {}
    for st in sorted({r["station"] for r in t34}):
        sub = [r for r in t34 if r["station"] == st]
        sx, sy = sum(r["x"] for r in sub), sum(r["y"] for r in sub)
        per[st] = {"n": len(sub), "sum_x": p5.rnd(sx), "sum_y": p5.rnd(sy), "R": p5.rnd(sy / sx) if sx else None}
    d["D4_6to12_by_station"] = per

    # 矛盾 B 映射
    c3, c4 = out["main"]["T3"]["judgement"]["category_code"], out["main"]["T4"]["judgement"]["category_code"]
    if c3 == 3 and c4 == 3:
        mp = "T3、T4 都不能区分 → 6–12 h 的 0.8–0.9 不能与 R1 的约 0.3 区分；矛盾 B 在统计上不成立为独立反证，只是点估计"
    elif 1 in (c3, c4):
        mp = "至少一档与 R1 相容（按实际写）"
    elif 2 in (c3, c4):
        mp = "6–12 h 至少一档与 R1 不容，矛盾 B 保留为未解释"
    else:
        mp = "至少一档不可评（按实际写）"
    out["contradiction_B_mapping"] = {"T3": c3, "T4": c4, "text": mp}

    with open(os.path.join(out_dir, "p7d_summary.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    cols = ["key", "kind", "label", "n", "clusters_station_season", "clusters_station", "stations", "sum_num", "sum_den",
            "R", "ss_lo", "ss_hi", "st_lo", "st_hi", "jk_lo", "jk_hi", "loo_lo", "loo_hi", "frac_den_nonneg",
            "p5_gate_met", "category", "station_same_side", "jackknife_same_side"]
    with open(os.path.join(out_dir, "p7d_bins.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        items = list(out["main"].items()) + [(k, v) for k, v in d.items() if isinstance(v, dict) and "judgement" in v]
        for k, b in items:
            jk = b.get("jackknife_station") or {}
            jt = jk.get("t_ci95") or [None, None]
            loo = jk.get("loo_range") or [None, None]
            ss = b.get("ci95_station_season") or [None, None]
            st = b.get("ci95_station") or [None, None]
            j = b["judgement"]
            w.writerow([k, b["kind"], b["label"], b.get("n"), b.get("clusters_station_season"), b.get("clusters_station"),
                        b.get("stations"), b.get("sum_num"), b.get("sum_den"), b.get("R"), ss[0], ss[1], st[0], st[1],
                        jt[0], jt[1], loo[0], loo[1], b.get("frac_boot_den_nonneg"), b["p5_gate_met"], j["category"],
                        j.get("station_same_side"), j.get("jackknife_same_side")])
            log(f"{k:<22} n={b.get('n')} K={b.get('clusters_station_season')}/{b.get('clusters_station')} "
                f"R={b.get('R')} ss={fmt_ci(ss)} st={fmt_ci(st)} jk={fmt_ci(jt if jt[0] is not None else None)} "
                f"frac={b.get('frac_boot_den_nonneg')} → {j['category']}")
    log(f"矛盾 B 映射：{mp}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest-synthetic", action="store_true")
    ap.add_argument("--p5-dir", default=P5_DIR_DEFAULT)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    if args.selftest_synthetic:
        res = selftest_synthetic()
        for s in res:
            print(("OK  " if s["ok"] else "FAIL"), s["name"], s["detail"])
        return 0 if all(s["ok"] for s in res) else 4
    out_dir = os.environ.get("REPRO_OUTPUT_DIR") or args.out
    if not out_dir:
        print("需要 $REPRO_OUTPUT_DIR 或 --out", file=sys.stderr)
        return 5
    os.makedirs(out_dir, exist_ok=True)
    try:
        log(f"{VERSION} start")
        rc = run(args, out_dir)
    except Exception:
        log("异常：" + traceback.format_exc())
        rc = 5
    log(f"done rc={rc}")
    with open(os.path.join(out_dir, "p7d_log.txt"), "w") as f:
        f.write("\n".join(LOG) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
