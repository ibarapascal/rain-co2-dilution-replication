#!/usr/bin/env python3
"""C1＋C2 诊断（阶段 cal-surface-diag；秒级，只读已有事件表）。

C1：热带子集（剔 KEO、Papa）卫星比值 ÷ SPURS-2 现场海面比值 Q 与 ln Q 的
    两样本独立簇 bootstrap 区间（卫星与现场各用已报口径的 p2.cluster_boot、种子 20260926 以逐位复现已有区间；
    现场重抽值用 PCG64(20260929) 置换后与卫星配对）；保守变体（卫星以站为簇、现场以事件为簇）与严格热带变体。
C2：各检验合并比值的分母份额诊断 share_max、share_top2、n_opposite、R_drop_max、n_eff。

用法：python3 cal_surface_diag.py --sp <输入目录根> --out <输出目录>
      （读 <sp>/p5sss-out/p5_events.csv、<sp>/p7a-out/p7a_events_jpl.csv、<sp>/p6-out/p6_events.csv、
        <sp>/p8a-out/p8a_events_wind.csv、<sp>/p7e-out/p7e_spurs1_events.csv、
        <sp>/p8b-out/p8b_events.csv；先核 sha256，不符退出 3；自测不过退出 3）
依赖：numpy；import 同目录 p2_rim_test（cluster_boot，未改）。
输出：cal_c1_c2.json。
Change Log:
  2026-09-29：初版。
"""
import argparse
import csv
import hashlib
import json
import math
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_rim_test as p2  # noqa: E402

SEED, B = p2.BOOT_SEED, p2.BOOT_B
PAIR_SEED = 20260929
MIDLAT = {"KEO", "Papa"}
STRICT_OUT = {"KEO", "Papa", "MOSEAN/WHOTS"}           # |lat|>20°
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
# [repro] --sp 下的相对路径由 run.py 摆放；sha 闸门（只在 [options] check_upstream_sha = true 时核对；缺省不核对，只记录实际 sha）期望值缺省是参考运行的产物前缀，
# [repro] 可在配置 [upstream_sha] 换（同 P7d）；p7a 不比对
INPUTS = {"p5": ("p5sss-out/p5_events.csv", _rp.expect_sha("cal12.p5_events_csv", "7a6d70758d575d180d9f0100579a9b7f6bc014350c9080ca9b5ab84ffd6985fa")),
          "p7a": ("p7a-out/p7a_events_jpl.csv", _rp.expect_sha("cal12.p7a_events_jpl_csv", None)),  # [repro] 原 None：只记录不比对
          "p6": ("p6-out/p6_events.csv", _rp.expect_sha("cal12.p6_events_csv", "2de4313d")),
          "p8a": ("p8a-out/p8a_events_wind.csv", _rp.expect_sha("cal12.p8a_events_wind_csv", "b8952fc6")),
          "p7e": ("p7e-out/p7e_spurs1_events.csv", _rp.expect_sha("cal12.p7e_spurs1_events_csv", "7c2e0a56")),
          "p8b": ("p8b-out/p8b_events.csv", _rp.expect_sha("cal12.p8b_events_csv", "ddba60b7"))}  # [repro]
EXPECT = {  # 自测：产物全精度值（RSS＝P5SS sensitivity.S8；JPL＝P7AS J5；P6＝P6S primary）；容差 0.002
    "RSS_TROP": (37, 13, 0.79986, (0.14645, 1.29223), 0.002), "JPL_TROP": (41, 14, 1.075, (0.539, 2.025), 0.002),
    "P6": (29, 8, 0.364, (0.150, 0.607), 0.002),
}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fnum(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def pct(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]


def boot(keys, num, den):
    out, K = p2.cluster_boot(keys, {"R": (np.asarray(num, float), np.asarray(den, float))}, B, SEED)
    return out["R"], K


def sat_rows(rows, drop):
    return [r for r in rows if r["station"] not in drop and fnum(r["y"]) is not None and fnum(r["x"]) is not None]


def sat_block(rows, by="ss"):
    keys = [(r["station"], r["season"]) for r in rows] if by == "ss" else [r["station"] for r in rows]
    y = [fnum(r["y"]) for r in rows]
    x = [fnum(r["x"]) for r in rows]
    bR, K = boot(keys, y, x)
    return {"n": len(rows), "clusters": K, "R": sum(y) / sum(x), "ci95": pct(bR)}, bR


def p6_block(rows, by="cluster"):
    keys = [r["cluster"] for r in rows] if by == "cluster" else list(range(len(rows)))
    o = [fnum(r["snake_obs"]) for r in rows]
    m = [fnum(r["snake_rim"]) for r in rows]
    bR, K = boot(keys, o, m)
    return {"n": len(rows), "clusters": K, "R": sum(o) / sum(m), "ci95": pct(bR)}, bR


def q_block(bs, bp, rs, rp):
    perm = np.random.Generator(np.random.PCG64(PAIR_SEED)).permutation(len(bp))
    q = bs / bp[perm]
    lo, hi = pct(q)
    point = rs / rp
    read = ("热带子集内卫星比值显著大于现场海面比值" if lo > 1 else
            "热带子集内卫星比值显著小于现场海面比值" if hi < 1 else
            "热带子集内两者之差不超出两样本簇抽样不确定度")
    return {"Q": point, "Q_ci95": [lo, hi], "lnQ": math.log(point) if point > 0 else None,
            "lnQ_ci95": [math.log(lo) if lo > 0 else None, math.log(hi) if hi > 0 else None],
            "lnQ_note": "下端 ≤0，ln 不定义" if lo <= 0 else "", "P_Q_le_1": float(np.mean(q[np.isfinite(q)] <= 1)),
            "frac_q_nonpos": float(np.mean(q <= 0)), "reading": read}


def share_block(obs, mod, label, expect=None):
    obs = np.asarray(obs, float)
    mod = np.asarray(mod, float)
    R = obs.sum() / mod.sum()
    if expect is not None and abs(R - expect) > 0.002:
        raise SystemExit(f"C2 自测不过：{label} R={R:.5f} 预期 {expect}")
    w = mod / mod.sum()
    order = np.argsort(-w)
    i = order[0]
    keep = np.ones(len(w), bool)
    keep[i] = False
    return {"label": label, "n": int(len(w)), "R": float(R), "share_max": float(w[i]), "share_top2": float(w[order[:2]].sum()),
            "n_opposite": int((w < 0).sum()), "R_drop_max": float(obs[keep].sum() / mod[keep].sum()),
            "n_eff": float(1.0 / np.sum(w ** 2)), "max_event_index": int(i)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sp", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ins = {}
    for k, (rel, want) in INPUTS.items():
        p = os.path.join(a.sp, rel)
        h = sha(p)
        if _rp.upstream_sha_mismatch(h, want):
            print(f"输入 sha256 不符：{rel} {h[:12]} ≠ {want}", file=sys.stderr)
            sys.exit(3)
        ins[k] = {"path": rel, "sha256": h}
    rd = lambda k: list(csv.DictReader(open(os.path.join(a.sp, INPUTS[k][0]), encoding="utf-8")))
    p5, p7a, p6 = rd("p5"), rd("p7a"), rd("p6")
    p6 = [r for r in p6 if fnum(r["snake_obs"]) is not None and fnum(r["snake_rim"]) is not None]

    # ---- C1 ----
    res = {"inputs": ins, "seed": SEED, "B": B, "pair_seed": PAIR_SEED, "c1": {}, "c2": {}}
    rss_t, jpl_t = sat_rows(p5, MIDLAT), sat_rows(p7a, MIDLAT)
    brss, b_rss = sat_block(rss_t)
    bjpl, b_jpl = sat_block(jpl_t)
    bp6, b_p6 = p6_block(p6)
    self_ok = {}
    for nm, blk in (("RSS_TROP", brss), ("JPL_TROP", bjpl), ("P6", bp6)):
        n, K, R, ci, tol = EXPECT[nm]
        ok = (blk["n"] == n and blk["clusters"] == K and abs(blk["R"] - R) <= tol
              and all(abs(u - v) <= 0.002 for u, v in zip(blk["ci95"], ci)))
        self_ok[nm] = {"ok": bool(ok), "this": blk, "registered": {"n": n, "clusters": K, "R": R, "ci95": ci}}
    res["selftest"] = self_ok
    if not all(v["ok"] for v in self_ok.values()):
        json.dump(res, open(os.path.join(a.out, "cal_c1_c2_selftest_fail.json"), "w"), ensure_ascii=False, indent=1)
        print("C1 自测不过", json.dumps(self_ok, ensure_ascii=False), file=sys.stderr)
        sys.exit(3)
    res["c1"]["inputs_blocks"] = {"RSS_TROP": brss, "JPL_TROP": bjpl, "P6": bp6}
    res["c1"]["primary"] = {"RSS": q_block(b_rss, b_p6, brss["R"], bp6["R"]), "JPL": q_block(b_jpl, b_p6, bjpl["R"], bp6["R"])}
    crss, cb_rss = sat_block(rss_t, "st")
    cjpl, cb_jpl = sat_block(jpl_t, "st")
    cp6, cb_p6 = p6_block(p6, "event")
    res["c1"]["conservative"] = {"blocks": {"RSS_TROP_st": crss, "JPL_TROP_st": cjpl, "P6_event": cp6},
                                 "RSS": q_block(cb_rss, cb_p6, crss["R"], cp6["R"]), "JPL": q_block(cb_jpl, cb_p6, cjpl["R"], cp6["R"])}
    srss, sb_rss = sat_block(sat_rows(p5, STRICT_OUT))
    sjpl, sb_jpl = sat_block(sat_rows(p7a, STRICT_OUT))
    res["c1"]["strict_tropics"] = {"blocks": {"RSS": srss, "JPL": sjpl},
                                   "RSS": q_block(sb_rss, b_p6, srss["R"], bp6["R"]), "JPL": q_block(sb_jpl, b_p6, sjpl["R"], bp6["R"])}

    # ---- C2 ----
    c2 = res["c2"]
    yx = lambda rows: ([fnum(r["y"]) for r in rows], [fnum(r["x"]) for r in rows])
    p5ok = sat_rows(p5, set())
    p7ok = sat_rows(p7a, set())
    c2["RSS"] = share_block(*yx(p5ok), "RSS 主结果", 1.65273)
    c2["JPL"] = share_block(*yx(p7ok), "JPL 主结果", 1.12056)
    c2["P6"] = share_block([fnum(r["snake_obs"]) for r in p6], [fnum(r["snake_rim"]) for r in p6], "SPURS-2 海面", 0.36386)
    w = rd("p8a")
    i5 = {(r["station"], r["onset_utc"]): r for r in p5ok}
    i7 = {(r["station"], r["onset_utc"]): r for r in p7ok}
    th_r = [i5[(r["station"], r["onset_utc"])] for r in w if r["product"] == "RSS" and r["bin_t"] == "H"]
    th_j = [i7[(r["station"], r["onset_utc"])] for r in w if r["product"] == "JPL" and r["bin_t"] == "H"]
    c2["RSSTH"] = share_block(*yx(th_r), "RSS 强风档（P8a TH）", 5.7437)
    c2["JPLTH"] = share_block(*yx(th_j), "JPL 强风档（P8a TH）", 3.2071)
    c2["RSSTROP"] = share_block(*yx(rss_t), "RSS 剔 KEO、Papa", brss["R"])
    c2["JPLTROP"] = share_block(*yx(jpl_t), "JPL 剔 KEO、Papa", bjpl["R"])
    s1 = [r for r in rd("p7e") if fnum(r["s086_obs"]) is not None and fnum(r["s086_rim"]) is not None]
    c2["S1"] = share_block([fnum(r["s086_obs"]) for r in s1], [fnum(r["s086_rim"]) for r in s1], "SPURS-1 0.865 m", 0.254)
    la = [r for r in rd("p8b") if fnum(r["la1m_obs"]) is not None and fnum(r["la1m_rim"]) is not None]
    c2["LA1M"] = share_block([fnum(r["la1m_obs"]) for r in la], [fnum(r["la1m_rim"]) for r in la], "Lady Amber 1 m", 0.4265)
    c2["_not_computed"] = {"SPURS2_buoy_1m": "无逐事件表（p6_events.csv 只有船载层），未算"}
    out = os.path.join(a.out, "cal_c1_c2.json")
    json.dump(res, open(out, "w"), ensure_ascii=False, indent=1)
    print(out, sha(out))


if __name__ == "__main__":
    main()
