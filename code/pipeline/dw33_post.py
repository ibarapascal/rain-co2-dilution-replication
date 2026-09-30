#!/usr/bin/env python3
"""C3（阶段 dw33）后处理（秒级）：dw33_bins.json → 汇总量。

U(s_T1,s_T2,s_T3)＝dep_pct＋(1−c)·pct(skin＋Σ_b[int0_b＋s_T(b)·(int1_b−int0_b)])；主估计量取系泊三档点估计；区间＝P4 第一轮
646 事件三档同一次站×季簇 bootstrap（p4_mech.terciles／boot_pairs，B＝10,000、种子 20260926，先复现 P4S D_U.cells）逐次代入；
界＝三档同取各自 95% 区间下端／上端；占比、系泊风况覆盖、线性组合对逐像元检查情形的误差、对 P3B 的核对；敏感性分档同报。

用法：python3 dw33_post.py --bins dw33_bins.json --p4e p4_mech_events.csv --p4s p4_mech_summary.json \
          --p3b p3b_curves.json --p3c p3c_decision.json --out DIR
      （bins 的 mode＝bench 时只测通流程，不做对 P3B 的全年核对，产物标 bench）
输出：dw33_post.json。
依赖：numpy、scipy（p4_mech 的 jackknife import）；import 同目录 p4_mech、p2_rim_test（未改）。
Change Log:
  2026-09-29：初版。
  2026-10-01：去掉 dw33_keys.md 输出（只写 dw33_post.json）；输入 sha 门按 [options] check_upstream_sha（缺省关）。
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
import p4_mech as pm  # noqa: E402

C = 0.30
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
# [repro] 输入 sha 闸门（只在 [options] check_upstream_sha = true 时核对；缺省不核对，只记录实际 sha）：期望值缺省是参考运行的产物前缀，可在配置 [upstream_sha] 换（同 P7d）
SHA = {k: _rp.expect_sha("dw33post." + k, v) for k, v in  # [repro]
       {"p4e": "c5feedaa", "p4s": "e3e2bf7b", "p3b": "2bdbe143", "p3c": "8a88e222"}.items()}  # [repro] 参考运行的产物前缀
COVER = {"minmax": (1.641173, 26.633166), "q95": (2.581536, 15.071225)}
GR = 2.5


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


def mooring_draws(p4e_path, p4s):
    """复现 P4 D-U 三档的站×季簇 bootstrap，返回 (点估计[3], 各档重抽[3,B], 复现核对)。"""
    rows = []
    for r in csv.DictReader(open(p4e_path, encoding="utf-8")):
        if r["model_ok"] != "True":
            continue
        d = {"station": r["station"], "season": r["season"], "obs1": fnum(r["obs1"]), "rim1": fnum(r["rim1"]),
             "u10": fnum(r["u10"]), "low1": fnum(r["low1"])}
        if None in (d["obs1"], d["rim1"], d["u10"], d["low1"]):
            continue
        rows.append(d)
    lab, cuts, _ = pm.terciles(rows, lambda r: r["u10"])
    pairs = {f"t{k}": ([r["obs1"] if lab[j] == k else 0.0 for j, r in enumerate(rows)],
                       [r["rim1"] if lab[j] == k else 0.0 for j, r in enumerate(rows)]) for k in range(3)}
    bb, K = pm.boot_pairs([(r["station"], r["season"]) for r in rows], pairs)
    pts = [pm.ratio_of([r for j, r in enumerate(rows) if lab[j] == k], "obs1", "rim1") for k in range(3)]
    ref = p4s["D_U"]["cells"]
    chk = []
    for k in range(3):
        ci = pm.ci95(bb[f"t{k}"])[0]
        chk.append({"k": k, "n": int((lab == k).sum()), "point": pts[k], "ci95": ci, "p4s_n": ref[k]["n"], "p4s_point": ref[k]["point"],
                    "p4s_ci95": ref[k]["ci95"],
                    "ok": bool(int((lab == k).sum()) == ref[k]["n"] and abs(pts[k] - ref[k]["point"]) <= 1e-4
                               and all(abs(a - b) <= 1e-4 for a, b in zip(ci, ref[k]["ci95"])))})
    cut_ok = all(abs(a - b) <= 1e-6 for a, b in zip(cuts, p4s["D_U"]["cuts"]))
    return np.array(pts), np.vstack([bb[f"t{k}"] for k in range(3)]), {"cells": chk, "cuts": cuts, "cuts_ok": cut_ok, "K": K,
                                                                        "n": len(rows)}


def analyze_variant(doc, v, pts, draws, lows, highs, u_m1, bench):
    b = {q: np.array(x, float) for q, x in doc["bacc"][v].items()}
    tert = np.array(doc["tercile_of_bin"])
    ue = np.array(doc["upper_edges"], float)
    nb = doc["n_bins"]
    lo_edge = np.r_[0.0, ue, np.inf][:nb]            # 细档下界（NaN 档无意义）
    hi_edge = np.r_[ue, np.inf, np.nan][:nb]
    valid = tert >= 0
    W = b["W_all"].sum()
    pct = lambda x: 100.0 * x / W
    dep = -100.0 * b["dep"].sum() / W
    skin = b["skin"].sum()
    dI = b["int1"] - b["int0"]
    dD = b["dil1"] - b["dil0"]
    U0 = dep + (1 - C) * pct(skin + b["int0"].sum())
    U1 = dep + (1 - C) * pct(skin + b["int1"].sum())
    a = np.array([(1 - C) * pct(dI[valid & (tert == k)].sum()) for k in range(3)])
    Ufun = lambda s: U0 + a @ np.asarray(s, float)
    Up = Ufun(pts)
    Ud = U0 + a @ draws
    ci = [float(np.percentile(Ud, 2.5)), float(np.percentile(Ud, 97.5))]
    dd = Ud - U1
    out = {"W_all_gC": W, "dep_pct": dep, "skin_pct": pct(skin), "U0_bins": U0, "U1_bins": U1, "slope_by_tercile_pp_per_s": a.tolist(),
           "U_wind": Up, "U_wind_ci95": ci, "Delta_vs_U1": Up - U1, "Delta_vs_U1_ci95": [float(np.percentile(dd, 2.5)), float(np.percentile(dd, 97.5))],
           "U_wind_minus_U_M1": Up - u_m1, "U_wind_minus_U_M1_ci95": [float(np.percentile(Ud - u_m1, 2.5)), float(np.percentile(Ud - u_m1, 97.5))],
           "bounds_all_low_high": [Ufun(lows), Ufun(highs)],
           "gr_position": ("Δ 区间整体低于 −2.5 pp（在『改写』判据带之外）" if np.percentile(dd, 97.5) < -GR else
                           "Δ 区间整体高于 +2.5 pp" if np.percentile(dd, 2.5) > GR else "Δ 区间与 ±2.5 pp 判据带相交"),
           "share_dI_by_tercile": [float(dI[valid & (tert == k)].sum() / dI[valid].sum()) for k in range(3)],
           "share_dil_by_tercile": [float(dD[valid & (tert == k)].sum() / dD[valid].sum()) for k in range(3)],
           "share_Wall_by_tercile": [float(b["W_all"][valid & (tert == k)].sum() / W) for k in range(3)],
           "nan_bin": {q: float(b[q][nb - 1]) for q in ("n_W", "W_all", "int1", "dil1")}}
    for nm, (lo, hi) in COVER.items():
        m = valid & (lo_edge >= lo - 1e-9) & (hi_edge <= hi + 1e-9)
        out[f"cover_{nm}_share_dI"] = float(dI[m].sum() / dI[valid].sum())
        out[f"cover_{nm}_share_dil"] = float(dD[m].sum() / dD[valid].sum())
        out[f"cover_{nm}_bins"] = [float(lo_edge[m].min()), float(hi_edge[m].max())] if m.any() else None
    out["fine_bins"] = [{"lo": float(lo_edge[k]), "hi": float(hi_edge[k]) if np.isfinite(hi_edge[k]) else None, "tercile": int(tert[k]),
                         "share_dI": float(dI[k] / dI[valid].sum()), "share_W": float(b["W_all"][k] / W)} for k in range(nb - 1)]
    if "intP" in doc["bacc"][v] and np.any(b["intP"] != 0):
        UP = dep + (1 - C) * pct(skin + b["intP"].sum())
        out["linearity_check"] = {"U_exact_pixelwise": UP, "U_linear": Up, "diff_pp": UP - Up}
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ("bins", "p4e", "p4s", "p3b", "p3c", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ins = {}
    for k in ("bins", "p4e", "p4s", "p3b", "p3c"):
        h = sha(getattr(a, k))
        if k in SHA and _rp.upstream_sha_mismatch(h, SHA[k]):
            print(f"输入 sha256 不符：{k} {h[:12]}", file=sys.stderr)
            sys.exit(3)
        ins[k] = {"path": getattr(a, k), "sha256": h}
    doc = json.load(open(a.bins, encoding="utf-8"))
    bench = doc.get("mode") == "bench"
    p4s = json.load(open(a.p4s, encoding="utf-8"))
    p3b = json.load(open(a.p3b, encoding="utf-8"))
    p3c = json.load(open(a.p3c, encoding="utf-8"))
    pts, draws, rep = mooring_draws(a.p4e, p4s)
    if not (all(c["ok"] for c in rep["cells"]) and rep["cuts_ok"]):
        print("系泊三档复现不过", json.dumps(rep, ensure_ascii=False), file=sys.stderr)
        sys.exit(3)
    if any(abs(x - y) > 1e-4 for x, y in zip(pts, doc["s_check_by_tercile"])):
        print("逐像元检查情形所用三档点估计与系泊复现不一致", file=sys.stderr)
        sys.exit(3)
    lows = np.array([c["ci95"][0] for c in rep["cells"]])
    highs = np.array([c["ci95"][1] for c in rep["cells"]])
    u_m1 = p3c["M1"]["U_curve"]
    res = {"inputs": ins, "mode": doc.get("mode"), "mooring": rep, "variants": {}}
    for v in doc["variants"]:
        res["variants"][v] = analyze_variant(doc, v, pts, draws, lows, highs, u_m1, bench)
    main_v = res["variants"]["cur"]
    g = p3b["global"]["whole"]
    ref = {"U1": g["U_pp"][g["s"].index(1.0)], "U0": g["U_pp"][g["s"].index(0.0)], "dep": p3b["deposition_pct"],
           "skin": p3b["skin_term_nonD_pct"]}
    res["p3b_check"] = {"applicable": not bench, "tol_pp": 1e-6,
                        **{k: {"bins": main_v[kk], "p3b": ref[k], "diff": main_v[kk] - ref[k]}
                           for k, kk in (("U1", "U1_bins"), ("U0", "U0_bins"), ("dep", "dep_pct"), ("skin", "skin_pct"))}}
    res["p3b_check"]["pass"] = (None if bench else all(abs(res["p3b_check"][k]["diff"]) <= 1e-6 for k in ("U1", "U0", "dep", "skin")))
    op = os.path.join(a.out, "dw33_post.json")
    json.dump(res, open(op, "w"), ensure_ascii=False, indent=1)
    h = sha(op)
    print(op, h)
    print("P3B 核对：", json.dumps(res["p3b_check"], ensure_ascii=False))


if __name__ == "__main__":
    main()
