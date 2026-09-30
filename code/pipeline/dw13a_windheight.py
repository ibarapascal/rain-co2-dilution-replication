#!/usr/bin/env python3
"""dw13a 风速计高度敏感性：P2 主集只改风速计高度的重算。

唯一改动：B 群（O 组 5 站）逐小时按「提供该小时风值的（文件, 变量）」取风速计实际高度 z，风速按 K9 同一中性对数律
×ln(10/z0)/ln(z/z0) 换到 10 m；A 群维持 4 m。其余全部复用 p2_rim_test.py（p2-2026-09-26c，不改一字）：
import 后只包装四处——p1b_extend.analyze_file（记下风变量的 DAS coordinates 属性）、p1b_extend.load_source 与
process_o_station（逐小时记下风值来源，条件与 J8 先到先得相同）、p2.rebuild_stations（重建后把 B 群小时风乘
r(z)＝f(z)/f(4)，原代码再 ×f(4)）、p2.event_record（只旁记全精度模型 ΔS，不改返回值）。
高度来源：dw13_meta.json（D13M，cal-dw13-meta 元数据核查产物）的 wind[文件].heights_m；取值规则 a／b／c 见 resolve_height()。

用法（cwd＝本目录；--meta 给出 dw13_meta.json）：
  dw13a_windheight.py --selftest
  dw13a_windheight.py --stage nominal  --resume-from ORIG --out O/nominal     # 全 4 m（r≡1.0），正确性门用
  dw13a_windheight.py --gate --orig ORIG --nominal O/nominal --out O          # 与原 P2 产物逐值比对；不过退出 4
  dw13a_windheight.py --stage measured --resume-from ORIG --out O/measured    # 实测高度
  dw13a_windheight.py --compare --orig ORIG --nominal O/nominal --measured O/measured --out O
产物：各 stage 目录＝p2_rim_test.py --part primary 的全部产物＋dw13a_heights.json、dw13a_events_full.csv；
  O/dw13a_gate.json、O/dw13a_compare.json。
退出码：0 完成；3 异常（含来源追踪核对不符、D13M sha 不符）；4 正确性门不过；其余透传 p2_rim_test.main。
Change Log:
  2026-09-29：初版。
"""
import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

VERSION = "dw13a-2026-09-29a"
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
# [repro] 核 import 的 p2_rim_test.py 等于 pipeline/ 里当前文件的 sha256（保证包装的是同目录这份 p2_rim_test.py）
P2_SHA = _rp.module_sha("p2_rim_test.py")
META_SHA = _rp.expect_sha("dw13a.dw13_meta", "fd5a036b")  # [repro] 输入 sha 闸门（只在 [options] check_upstream_sha = true 时核对；缺省不核对，只记录实际 sha）；期望值缺省＝参考运行的 D13M 前缀，可在 [upstream_sha] 换
META_DEFAULT = _rp.upstream("dw13_meta")  # [repro] 读 cal-dw13-meta 阶段产物
ORIG_DEFAULT = _rp.upstream("p2_dir")  # [repro] 读 p2 阶段输出目录
Z_NOM = 4.0
STRICT = ["primary.R_RIM_1m.point", "primary.R_RIM_1m.ci95", "primary.R_RIM_1m.ci90",
          "primary.R_S20_05m.point", "primary.R_S20_05m.ci95", "D3.verdict", "exit"]
ALLOW_SUMMARY = ["run_utc", "runtime_s", "inputs.resume_from", "cmorph.downloaded_this_run", "cmorph.bytes_this_run",
                 "cmorph.http_requests", "cmorph.seconds", "cmorph.pixels", "D3.written_utc"]
ALLOW_D3 = ["written_utc"]
NAN = float("nan")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def f10(z, z0=1e-4):
    """K9 同式：z 高度风速换到 10 m 的系数（与 p2.WIND_FACTOR 同一表达式，z=4.0 时逐位相同）。"""
    return math.log(10.0 / z0) / math.log(z / z0)


# ======================================================================== 高度查表
def load_meta(path):
    got = sha256_file(path)
    if _rp.upstream_sha_mismatch(got, META_SHA):
        raise RuntimeError(f"dw13_meta.json sha256 {got[:8]}… ≠ D13M {META_SHA}…")
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    by_file, prim = {}, {}
    for w in d["wind"]:
        base = os.path.basename(w["file"])
        by_file[base] = dict(w["heights_m"])
        prim.setdefault(w["station"], []).append(next(iter(w["heights_m"].values())))
    med = {st: statistics.median(v) for st, v in prim.items()}
    return {"sha256": got, "by_file": by_file, "station_median": med}


def resolve_height(meta, station, file, coords):
    """返回 (高度, 规则)：a＝coordinates 指到的高度坐标；b＝文件只登记一个高度；c＝该站中位数回退。"""
    hm = meta["by_file"].get(file)
    if hm:
        hit = [t for t in (coords or "").split() if t in hm]
        if len(hit) == 1:
            return float(hm[hit[0]]), "a:" + hit[0]
        if len(hm) == 1:
            k = next(iter(hm))
            return float(hm[k]), "b:" + k
    if station not in meta["station_median"]:
        raise RuntimeError(f"{station} 在 D13M 中无任何风文件，无法回退")
    return float(meta["station_median"][station]), "c:station_median"


# ======================================================================== 包装（不改原代码）
class Tracker:
    def __init__(self):
        self.cur = None
        self.hour_src = {}      # 站 -> {小时: (文件, 变量)}
        self.coords = {}        # (文件, 变量) -> DAS coordinates
        self.zarr = {}          # 站 -> 小时高度数组（无风值处 NaN）
        self.recs = []
        self.heights = {}


def install(tr, mode, meta):
    import numpy as np
    import p1_events as p1
    import p1b_extend as p1b
    import p2_rim_test as p2

    if sha256_file(p2.__file__) != P2_SHA:
        raise RuntimeError("p2_rim_test.py sha256 与预期不符")
    orig_af, orig_ls, orig_pos = p1b.analyze_file, p1b.load_source, p1b.process_o_station
    orig_rb, orig_er = p2.rebuild_stations, p2.event_record

    def analyze_file(http, cache, ent, log):
        rec = orig_af(http, cache, ent, log)
        winds = [s for s in rec.get("sources", []) if s["kind"] == "wind"]
        if winds:
            das_p = os.path.join(cache, "oceansites", ent["path"].rsplit("/", 1)[0], ent["file"] + ".das.gz")
            attrs = p1b.das_parse(p1.read_gz_text(das_p))
            for s in winds:
                tr.coords[(s["file"], s["vars"][0])] = attrs.get(s["vars"][0], {}).get("coordinates")
        return rec

    def load_source(http, cache, src, h_lo, h_hi, log, qc_drop):
        pairs, dt_s = orig_ls(http, cache, src, h_lo, h_hi, log, qc_drop)
        if src["kind"] == "wind" and tr.cur is not None:
            m = tr.hour_src.setdefault(tr.cur, {})
            key = (src["file"], src["vars"][0])
            for h, _v in pairs:                    # 与 process_o_station「0<=i<n 且 c[i]==0」同条件：先到先得
                if h_lo <= h <= h_hi and h not in m:
                    m[h] = key
        return pairs, dt_s

    def process_o_station(http, st, cache, log):
        tr.cur = st["name"]
        try:
            return orig_pos(http, st, cache, log)
        finally:
            tr.cur = None

    def rebuild_stations(http, cache, log, names_p1=None):
        stations = orig_rb(http, cache, log, names_p1)
        apply_heights(stations, tr, mode, meta, p2, np)
        return stations

    def event_record(stn, e, pix, variants=("cmorph",), windows=(p2.PRIMARY_WIN,), k27=False):
        rec = orig_er(stn, e, pix, variants, windows, k27)
        tr.recs.append(rec)
        return rec

    p1b.analyze_file, p1b.load_source, p1b.process_o_station = analyze_file, load_source, process_o_station
    p2.rebuild_stations, p2.event_record = rebuild_stations, event_record
    return p2


def apply_heights(stations, tr, mode, meta, p2, np):
    base = p2.WIND_FACTOR
    if f10(Z_NOM, p2.WIND_Z0_M) != base:
        raise RuntimeError("f10(4.0) 与 p2.WIND_FACTOR 不逐位相同")
    info = {"mode": mode, "stations": {}}
    for stn in stations:
        name = stn["name"]
        w = stn["ser"]["wind"]
        n = len(w)
        if stn["group"] != p2.GROUP_B:
            tr.zarr[name] = np.where(np.isfinite(w), Z_NOM, np.nan)
            continue
        hs = tr.hour_src.get(name, {})
        fin = set(np.nonzero(np.isfinite(w))[0].tolist())
        trk = {h - stn["h0"] for h in hs}
        if fin != trk:
            raise RuntimeError(f"{name}：来源追踪小时 {len(trk)} 与风序列有值小时 {len(fin)} 不一致"
                               f"（差 {len(fin ^ trk)}）")
        z = np.full(n, np.nan)
        r = np.ones(n)
        src_tab = {}
        for h, key in hs.items():
            if key not in src_tab:
                zz, rule = resolve_height(meta, name, key[0], tr.coords.get(key))
                src_tab[key] = {"file": key[0], "var": key[1], "coordinates": tr.coords.get(key),
                                "height_m": zz, "rule": rule, "hours": 0}
            src_tab[key]["hours"] += 1
            zz = Z_NOM if mode == "nominal" else src_tab[key]["height_m"]
            i = h - stn["h0"]
            z[i] = zz
            r[i] = f10(zz, p2.WIND_Z0_M) / base
        if mode == "nominal" and not np.all(r == 1.0):
            raise RuntimeError("nominal：r 不恒为 1.0")
        stn["ser"]["wind"] = w * r               # 新数组，不触及 series_raw
        tr.zarr[name] = z
        cnt = {}
        for v in z[np.isfinite(z)]:
            k = f"{v:.3f}"
            cnt[k] = cnt.get(k, 0) + 1
        info["stations"][name] = {"sources": sorted(src_tab.values(), key=lambda d: (d["file"], d["var"])),
                                  "hours_by_height_applied": dict(sorted(cnt.items())),
                                  "fallback_c": [d for d in src_tab.values() if d["rule"].startswith("c")]}
    tr.heights = info


def event_heights(tr, rec, p2):
    z = tr.zarr[rec["station"]]
    i = rec["i"]
    zo = float(z[i]) if 0 <= i < len(z) else NAN
    lo, hi = i + p2.WIN_LO_H, i + p2.WIN_HI_H
    seg = [float(x) for x in z[max(lo, 0):min(hi, len(z))] if x == x]
    kinds = sorted({round(x, 3) for x in seg})
    return zo, kinds


def write_stage_extras(tr, out, p2):
    W = p2.PRIMARY_WIN
    rows, by_st = [], {}
    for r in tr.recs:
        zo, kinds = event_heights(tr, r, p2)
        m = r["model"]
        row = {"station": r["station"], "group": r["group"], "onset_utc": r["onset_utc"],
               "z_onset_m": zo, "z_window_kinds": ";".join(f"{k:g}" for k in kinds), "n_z_kinds": len(kinds),
               "any_z_ne_4": int(any(k != Z_NOM for k in kinds)),
               "dS_rim_05": m.get(("cmorph", "rim", 0.5, W)), "dS_rim_1": m.get(("cmorph", "rim", 1.0, W)),
               "dS_rim_5": m.get(("cmorph", "rim", 5.0, W)), "dS_s20": m.get(("cmorph", "s20", 0.0, W)),
               "dS1_obs": r.get("obs", {}).get((1.0, W)), "dS05_obs": r.get("obs", {}).get((0.5, W))}
        rows.append(row)
        s = by_st.setdefault(r["station"], {"n": 0, "any_z_ne_4": 0, "straddle": 0, "z_onset": {}})
        s["n"] += 1
        s["any_z_ne_4"] += row["any_z_ne_4"]
        s["straddle"] += int(len(kinds) > 1)
        k = f"{zo:.3f}"
        s["z_onset"][k] = s["z_onset"].get(k, 0) + 1
    with open(os.path.join(out, "dw13a_events_full.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ["station"])
        w.writeheader()
        for row in rows:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in row.items()})
    info = dict(tr.heights, version=VERSION, events_by_station=by_st,
                n_events=len(rows), n_events_any_z_ne_4=sum(x["any_z_ne_4"] for x in rows),
                n_events_straddle=sum(1 for x in rows if x["n_z_kinds"] > 1))
    with open(os.path.join(out, "dw13a_heights.json"), "w", encoding="utf-8") as f:
        json.dump(info, f, ensure_ascii=False, indent=1)


def run_stage(args):
    meta = load_meta(args.meta)
    tr = Tracker()
    p2 = install(tr, args.stage, meta)
    os.makedirs(args.out, exist_ok=True)
    rc = p2.main(["--part", "primary", "--resume-from", args.resume_from, "--out", args.out])
    if rc != 0:
        print(f"p2_rim_test.main 返回 {rc}", flush=True)
        return rc
    write_stage_extras(tr, args.out, p2)
    print(f"[{args.stage}] 完成：{len(tr.recs)} 事件；heights → dw13a_heights.json", flush=True)
    return 0


# ======================================================================== 比对
def _decimals(x):
    s = repr(float(x))
    if "e" in s or "E" in s:
        return 12
    return len(s.split(".")[1]) if "." in s else 0


def _isnum(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _allowed(path, allow):
    return any(path == a or path.startswith(a + ".") or path.startswith(a + "[") for a in allow)


def deep_cmp(a, b, path, st, allow):
    if _allowed(path, allow):
        st["skipped"].append(path)
        return
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            st["fail"].append({"path": path, "why": "keys", "only_orig": sorted(set(a) - set(b))[:10],
                               "only_new": sorted(set(b) - set(a))[:10]})
        for k in sorted(set(a) & set(b)):
            deep_cmp(a[k], b[k], f"{path}.{k}" if path else k, st, allow)
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            st["fail"].append({"path": path, "why": "len", "orig": len(a), "new": len(b)})
            return
        for i, (x, y) in enumerate(zip(a, b)):
            deep_cmp(x, y, f"{path}[{i}]", st, allow)
        return
    st["n"] += 1
    if _isnum(a) and _isnum(b):
        fa, fb = float(a), float(b)
        if (fa != fa) or (fb != fb):
            ok = (fa != fa) and (fb != fb)
            st["exact" if ok else "fail_n"] += 1
            if not ok:
                st["fail"].append({"path": path, "orig": a, "new": b})
            return
        if fa == fb:
            st["exact"] += 1
            return
        d = abs(fa - fb)
        dec = max(_decimals(fa), _decimals(fb))
        st["max_abs_diff"] = max(st["max_abs_diff"], d)
        if dec >= 3 and d <= 10.0 ** (-dec) * (1 + 1e-9) and not _allowed(path, STRICT_CUR):
            st["tol"] += 1
            st["within_tol"].append({"path": path, "orig": a, "new": b})
        else:
            st["fail"].append({"path": path, "orig": a, "new": b})
        return
    if a == b:
        st["exact"] += 1
    else:
        st["fail"].append({"path": path, "orig": a, "new": b})


STRICT_CUR = []


def cmp_json(pa, pb, allow, strict):
    global STRICT_CUR
    STRICT_CUR = strict
    with open(pa, encoding="utf-8") as f:
        a = json.load(f)
    with open(pb, encoding="utf-8") as f:
        b = json.load(f)
    st = {"n": 0, "exact": 0, "tol": 0, "fail_n": 0, "max_abs_diff": 0.0, "within_tol": [], "fail": [], "skipped": []}
    deep_cmp(a, b, "", st, allow)
    return st


def cmp_csv(pa, pb):
    def rd(p):
        with open(p, encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))
        return {(r["station"], r["onset_utc"]): r for r in rows}, (list(rows[0]) if rows else [])
    A, ha = rd(pa)
    B, hb = rd(pb)
    st = {"n": 0, "exact": 0, "tol": 0, "fail_n": 0, "max_abs_diff": 0.0, "within_tol": [], "fail": [], "skipped": []}
    if ha != hb or set(A) != set(B):
        st["fail"].append({"why": "header_or_keys", "n_orig": len(A), "n_new": len(B)})
    for k in sorted(set(A) & set(B)):
        for c in ha:
            x, y = A[k][c], B[k][c]
            st["n"] += 1
            if x == y:
                st["exact"] += 1
                continue
            try:
                fx, fy = float(x), float(y)
            except ValueError:
                st["fail"].append({"path": f"{k}.{c}", "orig": x, "new": y})
                continue
            d = abs(fx - fy)
            dec = max(len(x.split(".")[1]) if "." in x else 0, len(y.split(".")[1]) if "." in y else 0)
            st["max_abs_diff"] = max(st["max_abs_diff"], d)
            if dec >= 3 and d <= 10.0 ** (-dec) * (1 + 1e-9):
                st["tol"] += 1
                st["within_tol"].append({"path": f"{k[0]}|{k[1]}.{c}", "orig": x, "new": y})
            else:
                st["fail"].append({"path": f"{k[0]}|{k[1]}.{c}", "orig": x, "new": y})
    return st


def _trim(st):
    return dict(st, within_tol=st["within_tol"][:60], fail=st["fail"][:60], n_within_tol_listed=len(st["within_tol"]),
                n_fail=len(st["fail"]), skipped=sorted(set(st["skipped"])))


def run_gate(args):
    res = {"version": VERSION, "rule": "运行记录类字段除外逐值比对；数值容差≤原值末位 1 单位（≥3 位小数时）；"
                                       "STRICT 列表须完全相等", "strict": STRICT, "orig": args.orig, "nominal": args.nominal}
    res["summary"] = _trim(cmp_json(os.path.join(args.orig, "p2_summary.json"),
                                    os.path.join(args.nominal, "p2_summary.json"), ALLOW_SUMMARY, STRICT))
    res["d3"] = _trim(cmp_json(os.path.join(args.orig, "p2_d3.json"), os.path.join(args.nominal, "p2_d3.json"),
                               ALLOW_D3, []))
    res["events_csv"] = _trim(cmp_csv(os.path.join(args.orig, "p2_events.csv"),
                                      os.path.join(args.nominal, "p2_events.csv")))
    for nm in ("p2_summary.json", "p2_d3.json", "p2_events.csv"):
        res.setdefault("sha256_orig", {})[nm] = sha256_file(os.path.join(args.orig, nm))
        res.setdefault("sha256_nominal", {})[nm] = sha256_file(os.path.join(args.nominal, nm))
    res["pass"] = all(res[k]["n_fail"] == 0 for k in ("summary", "d3", "events_csv"))
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "dw13a_gate.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
    print(f"[gate] pass={res['pass']} summary 精确 {res['summary']['exact']}/{res['summary']['n']} 容差内 "
          f"{res['summary']['tol']} 失败 {res['summary']['n_fail']}；d3 失败 {res['d3']['n_fail']}；events 精确 "
          f"{res['events_csv']['exact']}/{res['events_csv']['n']} 容差内 {res['events_csv']['tol']} 失败 "
          f"{res['events_csv']['n_fail']}", flush=True)
    return 0 if res["pass"] else 4


# ======================================================================== 新旧对比
def _read_full(p):
    with open(p, encoding="utf-8", newline="") as f:
        return {(r["station"], r["onset_utc"]): r for r in csv.DictReader(f)}


def _fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NAN


def _delta(o, n):
    if _isnum(o) and _isnum(n):
        return round(float(n) - float(o), 6)
    if isinstance(o, list) and isinstance(n, list) and len(o) == len(n):
        return [_delta(x, y) for x, y in zip(o, n)]
    return None


def _pair(o, n, keys):
    return {k: {"orig": o.get(k), "new": n.get(k), "delta": _delta(o.get(k), n.get(k))} for k in keys if k in o or k in n}


def run_compare(args):
    with open(os.path.join(args.orig, "p2_summary.json"), encoding="utf-8") as f:
        O = json.load(f)
    with open(os.path.join(args.measured, "p2_summary.json"), encoding="utf-8") as f:
        N = json.load(f)
    with open(os.path.join(args.measured, "dw13a_heights.json"), encoding="utf-8") as f:
        H = json.load(f)
    keys = ["point", "ci95", "ci90", "se_boot", "mde_ratio_units", "sum_num", "sum_den", "n_events", "n_clusters",
            "tost", "beta_model_pCO2_pooled", "beta_model_fCO2_pooled"]
    out = {"version": VERSION, "status": "敏感性：只改 B 群风速计高度；不改变任何既有判定",
           "orig": args.orig, "measured": args.measured, "nominal": args.nominal,
           "sha256_measured": {nm: sha256_file(os.path.join(args.measured, nm))
                               for nm in ("p2_summary.json", "p2_d3.json", "p2_events.csv", "dw13a_heights.json",
                                          "dw13a_events_full.csv")},
           "primary": {}, "strata_by_group": {}, "per_station_points": {}}
    for k in ("R_RIM_1m", "R_S20_05m", "P5", "beta_D4subset", "beta_all_info"):
        out["primary"][k] = _pair(O["primary"][k], N["primary"][k], keys)
        if k == "P5":
            out["primary"][k]["P5_RIM3"] = _pair(O["primary"][k]["P5_RIM3"], N["primary"][k]["P5_RIM3"],
                                                 ["point", "ci95", "diff_main_minus_this_ci95"])
    for g in O["strata_by_group"]:
        out["strata_by_group"][g] = {k: _pair(O["strata_by_group"][g][k], N["strata_by_group"][g][k], keys)
                                     for k in ("R_RIM_1m", "R_S20_05m", "P5", "beta_D4subset")}
    for s in O["per_station_points"]:
        o, n = O["per_station_points"][s], N["per_station_points"][s]
        out["per_station_points"][s] = {k: {"orig": o[k].get("point"), "new": n[k].get("point"),
                                            "delta": _delta(o[k].get("point"), n[k].get("point")), "n": n[k].get("n")}
                                        for k in ("R_RIM_1m", "R_S20_05m", "P5_obs")}
    out["k28_R_RIM_1m_by_layer"] = {k: _pair(O["k28_R_RIM_1m_by_layer"][k], N["k28_R_RIM_1m_by_layer"][k], ["point", "n_events"])
                                    for k in ("true_1m", "substitute_layer")}
    out["D3"] = {z: _pair(O["D3"][z], N["D3"][z], ["mean_dS_rim", "mean_dS_s20", "diff", "threshold_2se", "pass"])
                 for z in ("z=0.5m", "z=1.0m")}
    out["D3"]["verdict"] = {"orig": O["D3"]["verdict"], "new": N["D3"]["verdict"]}
    out["exit"] = {"orig": O["exit"]["label"], "new": N["exit"]["label"], "same": O["exit"] == N["exit"]}
    out["model_availability"] = {"orig": O["model_availability_k29"], "new": N["model_availability_k29"]}
    out["wind_hours_floored"] = {"orig": O["wind"]["hours_floored"], "new": N["wind"]["hours_floored"]}
    out["exploratory_k27_info"] = {k: _pair(O["exploratory_k27"][k], N["exploratory_k27"][k], ["point", "ci95"])
                                   for k in ("R_RIM_0.5m", "R_RIM_1.0m")}
    # 逐事件（全精度，nominal vs measured）
    A = _read_full(os.path.join(args.nominal, "dw13a_events_full.csv"))
    B = _read_full(os.path.join(args.measured, "dw13a_events_full.csv"))
    ev = {}
    for key in sorted(B):
        b, a = B[key], A.get(key)
        s = ev.setdefault(b["station"], {"n": 0, "any_z_ne_4": 0, "changed_dS_rim_1": 0, "rel_change_dS_rim_1": [],
                                         "rel_change_dS_s20": [], "z_onset": {}})
        s["n"] += 1
        s["any_z_ne_4"] += int(b["any_z_ne_4"])
        s["z_onset"][b["z_onset_m"]] = s["z_onset"].get(b["z_onset_m"], 0) + 1
        for col, lab in (("dS_rim_1", "rel_change_dS_rim_1"), ("dS_s20", "rel_change_dS_s20")):
            x, y = _fl(a[col]), _fl(b[col])
            if col == "dS_rim_1" and x == x and y == y and x != y:
                s["changed_dS_rim_1"] += 1
            if x == x and y == y and abs(x) > 1e-3:
                s[lab].append(y / x - 1.0)
    for s in ev.values():
        for lab in ("rel_change_dS_rim_1", "rel_change_dS_s20"):
            v = sorted(s.pop(lab))
            s[lab + "_summary"] = ({"n": len(v), "min": round(v[0], 5), "median": round(statistics.median(v), 5),
                                    "max": round(v[-1], 5)} if v else {"n": 0})
    out["events_by_station"] = ev
    out["n_events_any_z_ne_4"] = H["n_events_any_z_ne_4"]
    out["n_events_straddle"] = H["n_events_straddle"]
    out["fallback_c"] = {st: v["fallback_c"] for st, v in H["stations"].items() if v["fallback_c"]}
    with open(os.path.join(args.out, "dw13a_compare.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    r = out["primary"]["R_RIM_1m"]
    print(f"[compare] R_RIM(1 m) {r['point']['orig']}→{r['point']['new']}；CI95 {r['ci95']['orig']}→{r['ci95']['new']}；"
          f"z≠4 事件 {out['n_events_any_z_ne_4']}；出口相同 {out['exit']['same']}", flush=True)
    return 0


# ======================================================================== 自测（无网络，秒级）
def selftest():
    ok = True

    def check(c, msg):
        nonlocal ok
        print(("PASS " if c else "FAIL ") + msg)
        ok = ok and bool(c)

    base = math.log(10.0 / 1e-4) / math.log(4.0 / 1e-4)
    check(f10(4.0) == base and f10(4.0) / base == 1.0, "f10(4.0) 与 K9 表达式逐位相同、r=1.0")
    check(abs(f10(3.3) / base - 1 - 0.0185) < 5e-4, f"3.3 m 的系数比 4 m 大约 1.85%（{f10(3.3) / base - 1:.4f}）")
    meta = {"by_file": {"F_MET.nc": {"HEIGHT_WIND": 3.84, "HEIGHT_WIND2": 4.2}, "F_M.nc": {"HEIGHT_WND": 3.3},
                        "F_TV.nc": {"HEIGHTWIND": 4.0}, "F_S.nc": {"WSPD_H": 3.305}},
            "station_median": {"X": 3.5}}
    check(resolve_height(meta, "X", "F_MET.nc", "TIME HEIGHT_WIND2 LATITUDE LONGITUDE") == (4.2, "a:HEIGHT_WIND2"),
          "规则 a：coordinates 指到第二传感器")
    check(resolve_height(meta, "X", "F_MET.nc", "TIME HEIGHT_WIND LATITUDE LONGITUDE") == (3.84, "a:HEIGHT_WIND"),
          "规则 a：主传感器")
    check(resolve_height(meta, "X", "F_S.nc", "TIME NOMINAL_DEPTH LATITUDE LONGITUDE") == (3.305, "b:WSPD_H"),
          "规则 b：单高度文件")
    check(resolve_height(meta, "X", "F_TV.nc", None) == (4.0, "b:HEIGHTWIND"), "规则 b：无 coordinates")
    check(resolve_height(meta, "X", "UNKNOWN.nc", None) == (3.5, "c:station_median"), "规则 c：站中位数回退")
    check(resolve_height(meta, "X", "F_MET.nc", "TIME LATITUDE") == (3.5, "c:station_median"),
          "规则 c：多高度文件无法对应")
    st = {"n": 0, "exact": 0, "tol": 0, "fail_n": 0, "max_abs_diff": 0.0, "within_tol": [], "fail": [], "skipped": []}
    global STRICT_CUR
    STRICT_CUR = ["p.x"]
    deep_cmp({"p": {"x": 0.29616, "y": 0.12345, "z": 1.0, "w": 0.5, "t": "a"}, "run_utc": "u1"},
             {"p": {"x": 0.29616, "y": 0.12346, "z": 1.0, "w": 0.5, "t": "a"}, "run_utc": "u2"}, "", st, ["run_utc"])
    check(st["exact"] == 4 and st["tol"] == 1 and not st["fail"], "比对：末位 1 单位入容差，运行记录字段跳过")
    st = {"n": 0, "exact": 0, "tol": 0, "fail_n": 0, "max_abs_diff": 0.0, "within_tol": [], "fail": [], "skipped": []}
    deep_cmp({"p": {"x": 0.29616, "z": 1.0}}, {"p": {"x": 0.29617, "z": 1.1}}, "", st, [])
    check(len(st["fail"]) == 2, "比对：STRICT 键差 1 单位即失败；少于 3 位小数须完全相等")
    # 先到先得追踪与 process_o_station 同条件（合成：两个源重叠）
    m = {}
    for key, pairs in ((("A", "WSPD"), [(10, 1.0), (11, 2.0)]), (("A", "WSPD2"), [(11, 9.0), (12, 3.0), (99, 1.0)])):
        for h, _v in pairs:
            if 10 <= h <= 20 and h not in m:
                m[h] = key
    check(m == {10: ("A", "WSPD"), 11: ("A", "WSPD"), 12: ("A", "WSPD2")}, "来源追踪：先到先得、越界小时不记")
    print("selftest", "OK" if ok else "FAILED")
    return 0 if ok else 3


def main(argv=None):
    ap = argparse.ArgumentParser(description="dw13a 风速计高度敏感性")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--stage", choices=["nominal", "measured"])
    ap.add_argument("--gate", action="store_true")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--resume-from", default=ORIG_DEFAULT)
    ap.add_argument("--orig", default=ORIG_DEFAULT)
    ap.add_argument("--nominal")
    ap.add_argument("--measured")
    ap.add_argument("--out")
    ap.add_argument("--meta", default=META_DEFAULT)
    a = ap.parse_args(argv)
    try:
        if a.selftest:
            return selftest()
        if a.stage:
            return run_stage(a)
        if a.gate:
            return run_gate(a)
        if a.compare:
            return run_compare(a)
    except Exception:
        import traceback
        traceback.print_exc()
        return 3
    ap.print_help()
    return 3


if __name__ == "__main__":
    sys.exit(main())
