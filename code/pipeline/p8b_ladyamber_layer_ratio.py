#!/usr/bin/env python3
"""p8b_ladyamber_layer_ratio.py — P8b：SPURS-2 Lady Amber 同平台「海面（~1 cm）对 1 m」观测淡化层比，与 RIM-3 同事件层比对照，
直接检验 R1（RIM-3 幅度整体偏大：观测层比≈RIM 层比，Q_L≈1）与 R2（海面对、剖面太深：观测层比约为 RIM 层比的 1/R(1 m)≈3.3 倍）。

性质：探索性（设计晚于 P2、P5、P6、P7a–e；P6 结果已知）；方法与判读在读 Lady Amber 盐度数值之前写定。

用法：
  python p8b_ladyamber_layer_ratio.py --selftest [--out DIR]          合成数据自测（不读原件）
  python p8b_ladyamber_layer_ratio.py --load-check --root LA --out DIR
        装载核对：自测＋读原件＋小时化＋事件；只输出各小时序列有值小时数、每事件各层前 6 h／后 6 h／后 12 h 有值小时数、
        位移与雨量（不输出任何盐度数值或淡化量）→ p8b_loadcheck.json
  python p8b_ladyamber_layer_ratio.py --full --root LA --out DIR
        自测 → 读原件（只读）→ 小时化 → 事件 → RIM-3 → 先写 p8b_power.json → 再写 p8b_summary.json、p8b_events.csv
依赖：numpy、scipy、netCDF4＋同目录 p6_spurs2_profile.py（小时化、事件、RIM-3 调用、ΔS 配对、簇 bootstrap、判读工具；只 import 不改）
  → p2_rim_test.py（rim_factor、fill_wind、cluster_boot、model_hourly）→ p1_events.py、p1b_extend.py。路径全由命令行给。
产物（<out>/）：p8b_selftest.json、p8b_power.json（先写）、p8b_summary.json、p8b_events.csv、p8b_log.txt；--load-check 只写 p8b_loadcheck.json。
退出码：0 跑完（无论判读）；3 异常；4 自测不过。

实现选择（Y 条）：
  Y1 数据：Lady Amber（LAPS，帆船）第 2–7 航次各一个 `SPURS2_LadyAmber{2..7}.nc`（1 s，time＝days since 1950-01-01）；航次键 LA2…LA7。
     借用雨源（中心浮标、Revelle）不用：plan 实测 Lady Amber 离浮标 ≤50 km 的小时共 167 h、离 Revelle ≤50 km 41 h。
  Y2 雨：Lady Amber 自带 RM Young 雨量计。`rainfall_rate`（mm/hr）无雨时存为缺测（NaN），`rain_accumulation` 几乎全程有值 →
     「雨量计在线」＝该秒 rain_accumulation 有限；在线秒里 rainfall_rate 为 NaN 记 0、负值置 0；离线秒为缺测。
     事件用小时雨（≥50% 期望秒数在线，P1 I4）；RIM 用 30 min 半步（≥50% 在线，否则置 0 并计数，P6 W3／P2 K21）。
     S11 敏感性改用严格口径（rainfall_rate 本身有值的秒才算，≥50%）。
  Y3 风：√(zonal_wind_speed²＋meridional_wind_speed²)（Airmar PB200，~18 m，standard_name eastward/northward_wind＝对地），小时均值
     （≥50%），×ln(10/z0)/ln(18/z0)（z0＝p2.WIND_Z0_M，P2 K9 同式）到 10 m；之后 p6.forcing 原样（≤6 h 插补、下限 0.1 m/s）。
  Y4 层：`salinity_at_1cm`（SBE45，~1 cm）→ 键 "snake"、模型 z＝0（与 P6 海面、Witte 全球同；S8 改 0.01 m）；
     `salinity_at_1m`（SBE49，龙骨）→ "la1m"、z＝1.0；`salinity_at_2m` → 键 "tsg5"（沿用 p6.rel5_block 的键名；本检验里一律指 2 m）、z＝2.0。
     小时值＝小时内有效样本均值，需 ≥25% 期望样本（P6 W2）。船位 latitude／longitude 小时中位点（≥50%）。
  Y5 事件：p6.find_events（雨小时 ≥0.4 mm/h；onset 前 6 h 全有效无雨；[t0, t0+24 h) 全有效、累积 ≥10 mm）；
     [t0−30 h, t0+12 h) 须在雨量记录内且风可补齐（p6.forcing）。簇＝航次×ISO 周。
  Y6 RIM-3：p6.event_block → p6.model_hourly_F → p2.rim_factor（K2–K4：d0 7×7 双线性、48 历史半步＋当前半步、当前项 t=1 s），
     一行未改。自测：各深度与 p2.rim_factor 半步平均逐值相同（差 0），0.5／1／5 m 与 p2.model_hourly 逐值相同，z=0 与 p2._rim_literal <1e−12。
  Y7 ΔS：p6.delta_pair（雨后窗有值小时均值 − 雨前 6 h 中位数；模型按同一批小时配对；S0＝该层观测雨前中位数）；主窗 [0, 6) h。
  Y8 主估计量（配对子集＝同一事件 snake 与 la1m 都有效）：L_obs＝ΣΔS_obs(1 cm)/ΣΔS_obs(1 m)，L_RIM＝ΣΔS_RIM(0)/ΣΔS_RIM(1.0)，
     Q_L＝L_obs/L_RIM（≡ 同子集 R_sfc/R_1m）；同一簇抽样（p2.cluster_boot：B=10000、PCG64 seed 20260926、百分位 95%）；另报事件级 CI 与逐航次留一。
  Y9 可评（Q_L）：n≥15、簇≥4、≥2 个航次各 ≥3 事件；bootstrap 中 ΣΔS_RIM(0)≥0、ΣΔS_RIM(1.0)≥0 的比例各 ≤1%；ΣΔS_obs(1 m)≥0 的比例 ≤2.5%
     （1 m 观测淡化须可分辨，否则层比无定义）。
  Y10 判读 judge_q()：不可评 → CI 上端 <1.8「支持 R1」→ 下端 >1.8「支持 R2」→ 其余「不能区分」；稳健标注＝事件级 CI 与全部逐航次留一值同侧。
     1.8＝两预言 Q_L＝1 与 Q_L＝1/0.296≈3.38 的对数中点 1.84 向下取整到一位小数（1 与 3 的对数中点 1.73，两者之间取 1.8）。
  Y11 另报：R_sfc_obs（snake 有效的全部事件）按 P6 第六节阈值 0.65 分类（描述）；R_1m_obs（la1m 有效的全部事件）按 P6 浮标读法
     （CI 上端 <0.65 →「与 P2 同侧」）并标是否含 0.30；R_2m 描述；可评条件同 Y9 前四条。
  Y12 形状旁证 Q_12（la1m 与 2 m 都有效）：(ΣO_1/ΣO_2)/(ΣM_1/ΣM_2)；n≥8、簇≥4 且 ΣΔS_obs(2 m)≥0 的比例 ≤2.5% 才读：CI 下端 >1 →
     「观测淡化比 RIM-3 更集中在上层（R2 型形状方向）」；上端 <1 →「比 RIM-3 更深」；其余「不能区分」。不进主判读。
  Y13 雨中／雨停：[t0, t0+12 h) 内有值小时按该小时雨 ≥0.4 mm/h 分组（p6.rain_split 原样，分别对 snake 与 la1m）；
     Q_L 分组（两层同组都有效的事件，同一抽样）；各组 n<8 →「不可评」。
  Y14 R_rel2（1 cm − 2 m 同小时差分，p6.rel5_block，平流订正）：描述，标注与 Q_L 无关。
  Y15 敏感性（只报，不改判读；Q_L 与 R_1m 各报 n、点估计、CI、类别）：S1 干期 24 h；S2 干期 12 h；S3 阈值 5 mm；S4 窗 [0,3)；
     S5 窗 [0,12)；S6 位移 ≤30 km；S8 海面 z＝0.01 m；S9 当前项 t＝1800 s；S10 前后窗各 ≥4 小时；S11 雨量严格口径；
     S12 剔除 24 h 累积 >150 mm 的事件（雨量计尖峰防护：装载核对见 LA7 一个事件 352.65 mm、秒级雨量率文件级最大 2409 mm/h）。
  Y16 功效（p8b_power.json，先于任何观测 ΔS 合计写盘）：配对子集与 1 m 全集的 n、簇数、模型侧 Σx（0 m、1 m）；无雨对照窗
     （[c−12, c+6) 全有效无雨、c 每 6 h）上 snake、la1m 观测 ΔS 的标准差与推得的 SE；只作信息。
Change Log：
  2026-09-27 初版（P8b）。
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
from datetime import datetime, timezone

import numpy as np

import p6_spurs2_profile as p6
import p2_rim_test as p2

VERSION = "p8b-2026-09-27a"
NAN = float("nan")
Q_THRESH = 1.8
R_THRESH = p6.THRESH           # 0.65
R_REF = 0.30
OBS_DEN_NONNEG_MAX = 0.025
PRE_H, POST_MAX_H, RAIN_HOUR_MM = p6.PRE_H, p6.POST_MAX_H, p6.RAIN_HOUR_MM
WIND_Z_LA = 18.0
WIND_FACTOR_LA = math.log(10.0 / p2.WIND_Z0_M) / math.log(WIND_Z_LA / p2.WIND_Z0_M)
CRUISES = {f"LA{i}": f"SPURS2_LadyAmber{i}.nc" for i in range(2, 8)}
LAYERS = {"snake": ("salinity_at_1cm", 0.0), "la1m": ("salinity_at_1m", 1.0), "tsg5": ("salinity_at_2m", 2.0)}
PRIMARY = {"name": "primary", "dry_h": 6, "thr": 10.0, "win": (0, 6), "sfc": "snake", "z_sfc": 0.0,
           "t_cur": p2.T_CURRENT_S, "disp_max": None, "min_hours": 1, "rain": "on"}
SENS = [dict(PRIMARY, name="S1_dry24", dry_h=24), dict(PRIMARY, name="S2_dry12", dry_h=12),
        dict(PRIMARY, name="S3_thr5", thr=5.0), dict(PRIMARY, name="S4_win0_3", win=(0, 3)),
        dict(PRIMARY, name="S5_win0_12", win=(0, 12)), dict(PRIMARY, name="S6_disp30", disp_max=30.0),
        dict(PRIMARY, name="S8_z001", z_sfc=0.01), dict(PRIMARY, name="S9_tcur1800", t_cur=1800.0),
        dict(PRIMARY, name="S10_min4h", min_hours=4), dict(PRIMARY, name="S11_rain_strict", rain="strict"),
        dict(PRIMARY, name="S12_acc_le150", acc_max=150.0)]
LOG = []


def log(msg):
    line = f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} {msg}"
    LOG.append(line)
    print(line, flush=True)


def isnum(x):
    return p6.isnum(x)


def rnd(x, nd=5):
    return p6.rnd(x, nd)


def jdump(obj, path):
    p6.jdump(obj, path)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ======================================================================== 读原件与小时化（Y1–Y4）
def rain_online(rate, acc):
    """Y2：在线（acc 有限）秒里 rate 的 NaN 记 0、负值置 0；离线秒 NaN。"""
    r = np.where(np.isfinite(rate), np.maximum(rate, 0.0), 0.0)
    return np.where(np.isfinite(acc), r, np.nan)


def load_la(root, S, info):
    import netCDF4
    for cr, fn in CRUISES.items():
        ds = netCDF4.Dataset(os.path.join(root, fn))
        t = p6.nc_time(ds)
        rate, acc = p6.nc_var(ds, "rainfall_rate"), p6.nc_var(ds, "rain_accumulation")
        on = rain_online(rate, acc)
        strict = np.maximum(rate, 0.0)
        for tag, r in (("", on), ("_strict", strict)):
            k0, hr = p6.bin_mean(t, r, 3600.0, p6.MIN_FRAC_RAIN)
            S[f"{cr}|rain_h{tag}"] = p6.Series(k0, hr)
            k0, hs = p6.bin_mean(t, r, 1800.0, p6.MIN_FRAC_RAIN)
            S[f"{cr}|rain_hs{tag}"] = p6.Series(k0, hs)
        del rate, acc, on, strict
        spd = np.hypot(p6.nc_var(ds, "zonal_wind_speed"), p6.nc_var(ds, "meridional_wind_speed"))
        k0, w = p6.bin_mean(t, spd, 3600.0, p6.MIN_FRAC_RAIN)
        S[f"{cr}|wind_h"] = p6.Series(k0, w * WIND_FACTOR_LA if len(w) else w)
        del spd
        for q, vn in (("lat_h", "latitude"), ("lon_h", "longitude")):
            k0, m = p6.bin_mean(t, p6.nc_var(ds, vn), 3600.0, p6.MIN_FRAC_RAIN)
            S[f"{cr}|{q}"] = p6.Series(k0, m)
        present = {}
        for ln, (vn, _z) in LAYERS.items():
            if vn in ds.variables:
                k0, m = p6.bin_mean(t, p6.nc_var(ds, vn), 3600.0, p6.MIN_FRAC_SAL)
                S[f"{cr}|{ln}"] = p6.Series(k0, m)
                present[ln] = int(np.isfinite(m).sum()) if len(m) else 0
            else:
                S[f"{cr}|{ln}"] = p6.Series(0, np.array([]))
                present[ln] = None
        info[cr] = {"file": fn, "wind_factor": WIND_FACTOR_LA, "layer_valid_hours": present,
                    "rain_valid_hours": int(np.isfinite(S[f"{cr}|rain_h"].a).sum()),
                    "rain_strict_valid_hours": int(np.isfinite(S[f"{cr}|rain_h_strict"].a).sum()),
                    "wind_valid_hours": int(np.isfinite(S[f"{cr}|wind_h"].a).sum())}
        ds.close()
        del t
    return S, info


def rain_view(S, variant):
    """S11：把严格口径换进 rain_h／rain_hs（浅拷贝）。"""
    if variant == "on":
        return S
    T = dict(S)
    for cr in CRUISES:
        T[f"{cr}|rain_h"] = S[f"{cr}|rain_h_strict"]
        T[f"{cr}|rain_hs"] = S[f"{cr}|rain_hs_strict"]
    return T


# ======================================================================== 事件集（Y5–Y7）
def build(S0, cfg):
    S = rain_view(S0, cfg.get("rain", "on"))
    recs, drop = [], {}
    uid = 0
    for cr in CRUISES:
        for H, acc in p6.find_events(S[f"{cr}|rain_h"], cfg["dry_h"], cfg["thr"]):
            if cfg.get("acc_max") is not None and acc > cfg["acc_max"]:
                drop["acc_max"] = drop.get("acc_max", 0) + 1
                continue
            layers = {ln: (f"{cr}|{ln}", z) for ln, (_, z) in LAYERS.items()}
            b = p6.event_block(S, cr, H, layers, cfg)
            if b["diag"].get("reason"):
                drop[b["diag"]["reason"]] = drop.get(b["diag"]["reason"], 0) + 1
                continue
            la1, lo1 = np.nanmedian(S[f"{cr}|lat_h"].get(H - 6, H)), np.nanmedian(S[f"{cr}|lon_h"].get(H - 6, H))
            la2, lo2 = np.nanmedian(S[f"{cr}|lat_h"].get(H, H + 6)), np.nanmedian(S[f"{cr}|lon_h"].get(H, H + 6))
            disp = p6.haversine(la1, lo1, la2, lo2) if all(np.isfinite([la1, lo1, la2, lo2])) else NAN
            if cfg["disp_max"] is not None and not (np.isfinite(disp) and disp <= cfg["disp_max"]):
                drop["disp"] = drop.get("disp", 0) + 1
                continue
            _, iw, _ = datetime.fromtimestamp(H * 3600, timezone.utc).isocalendar()
            uid += 1
            recs.append({"uid": uid, "cruise": cr, "H": H, "onset": p6.iso(H), "acc24": acc, "cluster": f"{cr}-W{iw:02d}",
                         "disp_km": disp, "diag": b["diag"], "L": b["layers"],
                         "z": {ln: (cfg["z_sfc"] if ln == "snake" else z) for ln, (_, z) in LAYERS.items()},
                         "lat": la1, "rain12": S[f"{cr}|rain_h"].get(H, H + POST_MAX_H)})
    return recs, drop


# ======================================================================== 统计（Y8–Y14）
def per_cruise(rows):
    per = {}
    for r in rows:
        per[r["cruise"]] = per.get(r["cruise"], 0) + 1
    return per


def base_conds(n, K, per, n_min=15, k_min=4):
    return {f"n>={n_min}": n >= n_min, f"clusters>={k_min}": (K or 0) >= k_min,
            ">=2_cruises_with>=3": sum(1 for v in per.values() if v >= 3) >= 2}


def boot_sum_frac(rows, f, clus_f):
    """同 p2.cluster_boot 的抽样（同 seed、同排序簇）下，Σf ≥0 的比例。"""
    uk = sorted({clus_f(r) for r in rows})
    kid = {k: j for j, k in enumerate(uk)}
    cid = np.array([kid[clus_f(r)] for r in rows])
    s = np.bincount(cid, weights=np.array([f(r) for r in rows], float), minlength=len(uk))
    rng = np.random.Generator(np.random.PCG64(p6.SEED))
    draw = rng.integers(0, len(uk), size=(p6.B, len(uk)))
    return float((s[draw].sum(axis=1) >= 0).mean())


def q_stats(rows, a, b, clus_f):
    """同一抽样下 R_a、R_b、Q＝R_a/R_b、观测与模型层比。"""
    bc, K = p6.boot(rows, {"Ra": (p6.f_obs(a), p6.f_rim(a)), "Rb": (p6.f_obs(b), p6.f_rim(b)),
                           "LRo": (p6.f_obs(a), p6.f_obs(b)), "LRm": (p6.f_rim(a), p6.f_rim(b))}, clus_f)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = bc["Ra"] / bc["Rb"]
    return q, bc, K


def q_block(rows, a, b):
    out = {"n": len(rows)}
    if len(rows) < 3:
        return dict(out, note="事件 <3，不算", clusters=len({r["cluster"] for r in rows}), per_cruise=per_cruise(rows),
                    ci95_cluster=[None, None], ci95_event=[None, None], point=None)
    so, sm = sum(r["L"][a]["obs"] for r in rows), sum(r["L"][a]["rim"] for r in rows)
    oo, om = sum(r["L"][b]["obs"] for r in rows), sum(r["L"][b]["rim"] for r in rows)
    q, bc, K = q_stats(rows, a, b, p6.clus)
    qe, _, _ = q_stats(rows, a, b, lambda r: r["uid"])
    loo = {}
    for cr in sorted({r["cruise"] for r in rows}):
        rs = [r for r in rows if r["cruise"] != cr]
        if rs:
            A = sum(r["L"][a]["obs"] for r in rs) / sum(r["L"][a]["rim"] for r in rs)
            Bv = sum(r["L"][b]["obs"] for r in rs) / sum(r["L"][b]["rim"] for r in rs)
            loo[f"drop_{cr}"] = rnd(A / Bv) if Bv else None
        else:
            loo[f"drop_{cr}"] = None
    out.update({"clusters": K, "per_cruise": per_cruise(rows),
                "sum_obs_a": rnd(so), "sum_rim_a": rnd(sm), "sum_obs_b": rnd(oo), "sum_rim_b": rnd(om),
                f"R_{a}": rnd(so / sm) if sm else None, f"R_{a}_ci95": p6.pct(bc["Ra"]),
                f"R_{b}": rnd(oo / om) if om else None, f"R_{b}_ci95": p6.pct(bc["Rb"]),
                "point": rnd((so / sm) / (oo / om)) if (sm and om and oo) else None,
                "ci95_cluster": p6.pct(q), "ci95_event": p6.pct(qe), "loo_cruise": loo,
                "layer_ratio_obs": rnd(so / oo) if oo else None, "layer_ratio_obs_ci95": p6.pct(bc["LRo"]),
                "layer_ratio_rim": rnd(sm / om) if om else None, "layer_ratio_rim_ci95": p6.pct(bc["LRm"]),
                "frac_boot_rim_a_nonneg": rnd(boot_sum_frac(rows, p6.f_rim(a), p6.clus), 4),
                "frac_boot_rim_b_nonneg": rnd(boot_sum_frac(rows, p6.f_rim(b), p6.clus), 4),
                "frac_boot_obs_b_nonneg": rnd(boot_sum_frac(rows, p6.f_obs(b), p6.clus), 4)})
    return out


def q_evaluable(blk):
    per = blk.get("per_cruise") or {}
    conds = base_conds(blk["n"], blk.get("clusters"), per)
    fa, fb, fo = blk.get("frac_boot_rim_a_nonneg"), blk.get("frac_boot_rim_b_nonneg"), blk.get("frac_boot_obs_b_nonneg")
    conds["frac_rim_sfc_nonneg<=0.01"] = fa is not None and fa <= p6.DEN_NONNEG_MAX
    conds["frac_rim_1m_nonneg<=0.01"] = fb is not None and fb <= p6.DEN_NONNEG_MAX
    conds["frac_obs_1m_nonneg<=0.025"] = fo is not None and fo <= OBS_DEN_NONNEG_MAX
    return all(conds.values()), conds


def q_category(lo, hi):
    if lo is None or hi is None:
        return "不可评"
    if hi < Q_THRESH:
        return "支持 R1（观测层比不高于 RIM-3 层比的 1.8 倍）"
    if lo > Q_THRESH:
        return "支持 R2（观测层比高于 RIM-3 层比的 1.8 倍）"
    return "不能区分"


def judge_q(blk, ok):
    """Y10：层比 Q_L 的判读。"""
    if not ok:
        return {"category": "不可评", "rule": "6.1 第 1 条", "robust": None}
    lo, hi = blk["ci95_cluster"]
    cat = q_category(lo, hi)
    rule = {"不可评": "6.1 第 1 条（CI 不可算）", "不能区分": "6.1 第 4 条"}.get(cat, "6.1 第 2 条" if cat.startswith("支持 R1") else "6.1 第 3 条")
    robust = None
    if cat.startswith("支持"):
        elo, ehi = blk["ci95_event"]
        lv = [v for v in blk["loo_cruise"].values() if v is not None]
        if cat.startswith("支持 R1"):
            same = ehi is not None and ehi < Q_THRESH and len(lv) == len(blk["loo_cruise"]) and max(lv) < Q_THRESH
        else:
            same = elo is not None and elo > Q_THRESH and len(lv) == len(blk["loo_cruise"]) and min(lv) > Q_THRESH
        robust = "稳健" if same else "依赖簇定义"
    return {"category": cat, "rule": rule, "robust": robust}


def layer_R(recs, ln):
    rows = [r for r in recs if p6.lay_ok(r, ln)]
    if len(rows) < 3:
        return {"n": len(rows), "note": "事件 <3，不算", "per_cruise": per_cruise(rows), "evaluable": False}
    blk = p6.strip(p6.ratio_block(rows, p6.f_obs(ln), p6.f_rim(ln), p6.clus))
    per = per_cruise(rows)
    conds = base_conds(blk["n"], blk.get("clusters"), per)
    conds["frac_den_nonneg<=0.01"] = blk.get("frac_boot_den_nonneg") is not None and blk["frac_boot_den_nonneg"] <= p6.DEN_NONNEG_MAX
    blk.update({"per_cruise": per, "evaluable": all(conds.values()), "evaluable_conditions": conds,
                "loo_cruise": p6.loo_cruise(rows, p6.f_obs(ln), p6.f_rim(ln))})
    return blk


def reading_1m(blk):
    lo, hi = blk.get("ci95_cluster") or [None, None]
    if not blk.get("evaluable") or lo is None:
        return "不可评"
    s = "与 P2 同侧（CI 上端 <0.65）" if hi < R_THRESH else "未落在 P2 同侧（CI 上端 ≥0.65）"
    return s + ("；CI 含 0.30" if lo <= R_REF <= hi else "；CI 不含 0.30")


def reading_12(blk):
    if blk.get("n", 0) < 8 or (blk.get("clusters") or 0) < 4 or blk.get("ci95_cluster", [None])[0] is None:
        return "不可评（n<8 或簇<4 或 CI 不可算）"
    fo = blk.get("frac_boot_obs_b_nonneg")
    if fo is None or fo > OBS_DEN_NONNEG_MAX:
        return "不可评（2 m 观测淡化不可分辨）"
    lo, hi = blk["ci95_cluster"]
    if lo > 1.0:
        return "观测淡化比 RIM-3 更集中在上层（R2 型形状方向）"
    if hi < 1.0:
        return "观测淡化比 RIM-3 更深"
    return "不能区分（CI 含 1）"


def q_rain_split(recs, cfg):
    """Y13：Q_L 分组（两层同组都有效）。"""
    out = {}
    for grp in ("rain", "dry"):
        rows = []
        for r in recs:
            wet = np.isfinite(r["rain12"]) & (r["rain12"] >= RAIN_HOUR_MM)
            mask = wet if grp == "rain" else (np.isfinite(r["rain12"]) & ~wet)
            L = {}
            for ln in ("snake", "la1m"):
                X = r["L"].get(ln)
                if X is None:
                    break
                dso, dsm, _, _, npost = p6.delta_pair(X["obs18"], X["F18"], (0, POST_MAX_H), 1, mask_post=mask)
                if not (isnum(dso) and isnum(dsm)):
                    break
                L[ln] = {"obs": dso, "rim": dsm}
            if len(L) == 2:
                rows.append(dict(r, L=L))
        blk = q_block(rows, "snake", "la1m")
        blk["reading"] = ("不可评（n<8）" if blk["n"] < 8 or blk.get("ci95_cluster", [None])[0] is None
                          else q_category(*blk["ci95_cluster"]) + "（描述，不进主判读）")
        out[grp] = blk
    return out


def power_block(recs, S, cfg):
    """Y16：只用模型 x 与无雨对照窗观测噪声；先于任何观测 ΔS 合计写盘。"""
    pair = [r for r in recs if p6.lay_ok(r, "snake") and p6.lay_ok(r, "la1m")]
    one = [r for r in recs if p6.lay_ok(r, "la1m")]
    sfc = [r for r in recs if p6.lay_ok(r, "snake")]
    out = {"pair": {"n": len(pair), "clusters": len({r["cluster"] for r in pair}), "per_cruise": per_cruise(pair),
                    "sum_rim_sfc": rnd(sum(r["L"]["snake"]["rim"] for r in pair)) if pair else None,
                    "sum_rim_1m": rnd(sum(r["L"]["la1m"]["rim"] for r in pair)) if pair else None},
           "sfc_all": {"n": len(sfc), "clusters": len({r["cluster"] for r in sfc}), "per_cruise": per_cruise(sfc),
                       "sum_rim": rnd(sum(r["L"]["snake"]["rim"] for r in sfc)) if sfc else None},
           "one_m_all": {"n": len(one), "clusters": len({r["cluster"] for r in one}), "per_cruise": per_cruise(one),
                         "sum_rim": rnd(sum(r["L"]["la1m"]["rim"] for r in one)) if one else None}}
    for ln in ("snake", "la1m"):
        ctrl = []
        for cr in CRUISES:
            rh = S[f"{cr}|rain_h"]
            a, b = rh.span()
            for c in range(a + 12, b - 6, 6):
                seg = rh.get(c - 12, c + 6)
                if not (np.isfinite(seg).all() and (seg < RAIN_HOUR_MM).all()):
                    continue
                dso = p6.delta_pair(S[f"{cr}|{ln}"].get(c - PRE_H, c + POST_MAX_H), np.zeros(18), cfg["win"], 1)[0]
                if isnum(dso):
                    ctrl.append(dso)
        ctrl = np.array(ctrl, float)
        sd = rnd(ctrl.std(ddof=1)) if len(ctrl) > 2 else None
        rows = one if ln == "la1m" else sfc
        sx = sum(r["L"][ln]["rim"] for r in rows) if rows else 0.0
        out[f"ctrl_{ln}"] = {"n_windows": int(len(ctrl)), "sd_ctrl_dS": sd,
                             "se_R_ctrl": rnd(sd * math.sqrt(len(rows)) / abs(sx)) if (sd and rows and sx) else None}
    out["note"] = "只含事件数、模型侧 Σx 与无雨对照窗观测 ΔS 的离散度；只作信息，不改判读"
    return out


def run_cfg(S, cfg):
    recs, drop = build(S, cfg)
    pair = [r for r in recs if p6.lay_ok(r, "snake") and p6.lay_ok(r, "la1m")]
    qb = q_block(pair, "snake", "la1m")
    ok, conds = q_evaluable(qb) if pair else (False, {})
    j = judge_q(qb, ok)
    return {"cfg": dict(cfg), "n_events_all": len(recs), "dropped": drop,
            "Q_L": dict(qb, evaluable=ok, evaluable_conditions=conds), "judgement": j}, recs


def events_csv(recs, path):
    lays = list(LAYERS)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cruise", "onset", "acc24", "cluster", "disp_km", "lat"] +
                   [f"{ln}_{k}" for ln in lays for k in ("z", "obs", "rim", "s0", "npre", "npost")] + ["rain_hours_0_12"])
        for r in recs:
            row = [r["cruise"], r["onset"], rnd(r["acc24"], 2), r["cluster"], rnd(r["disp_km"], 1), rnd(r["lat"], 3)]
            for ln in lays:
                L = r["L"].get(ln, {})
                row += [r["z"].get(ln), rnd(L.get("obs")), rnd(L.get("rim")), rnd(L.get("s0"), 4), L.get("n_pre"), L.get("n_post")]
            row.append(int(np.nansum(r["rain12"] >= RAIN_HOUR_MM)))
            w.writerow(row)


def code_shas():
    here = os.path.dirname(os.path.abspath(__file__))
    return {fn: sha256_file(os.path.join(here, fn)) for fn in
            ("p8b_ladyamber_layer_ratio.py", "p6_spurs2_profile.py", "p2_rim_test.py", "p1_events.py", "p1b_extend.py")}


def avail(S, key, H):
    a = S[key].get(H - PRE_H, H + POST_MAX_H)
    f = np.isfinite(a)
    return [int(f[:PRE_H].sum()), int(f[PRE_H:PRE_H + 6].sum()), int(f[PRE_H:].sum())]


def run_load_check(a):
    os.makedirs(a.out, exist_ok=True)
    st = selftest()
    out = {"version": VERSION, "selftest": st, "code_sha256": code_shas(),
           "note": "只输出有值小时计数、事件时刻、雨量、位移；不含任何盐度数值或淡化量"}
    S, info = load_la(a.root, {}, {})
    out["series"] = info
    cnt = {}
    for cfg in [PRIMARY] + SENS:
        Sv = rain_view(S, cfg.get("rain", "on"))
        evs = []
        for cr in CRUISES:
            for H, acc in p6.find_events(Sv[f"{cr}|rain_h"], cfg["dry_h"], cfg["thr"]):
                P, U, dg = p6.forcing(Sv, cr, H)
                evs.append({"cruise": cr, "onset": p6.iso(H), "acc24": rnd(acc, 2), "model_ok": P is not None,
                            "reason": dg.get("reason"), "avail_pre6_post6_post12": {ln: avail(Sv, f"{cr}|{ln}", H) for ln in LAYERS}})
        ok = [e for e in evs if e["model_ok"]]

        def both(e, ln):
            v = e["avail_pre6_post6_post12"][ln]
            return v[0] >= 1 and v[1] >= 1
        cnt[cfg["name"]] = {"n_events": len(evs), "n_model_ok": len(ok),
                            "n_snake": sum(both(e, "snake") for e in ok), "n_1m": sum(both(e, "la1m") for e in ok),
                            "n_2m": sum(both(e, "tsg5") for e in ok),
                            "n_pair_snake_1m": sum(both(e, "snake") and both(e, "la1m") for e in ok),
                            "n_pair_1m_2m": sum(both(e, "la1m") and both(e, "tsg5") for e in ok),
                            "events": evs if cfg["name"] == "primary" else None}
    out["counts"] = cnt
    jdump(out, os.path.join(a.out, "p8b_loadcheck.json"))
    log("load-check done")
    return 0


def run_full(a):
    t_start = time.time()
    os.makedirs(a.out, exist_ok=True)
    st = selftest()
    jdump(st, os.path.join(a.out, "p8b_selftest.json"))
    if not st["all_ok"]:
        log("自测不过，退出 4")
        return 4
    code_sha = code_shas()
    log(f"code sha {json.dumps(code_sha)}")
    S, info = load_la(a.root, {}, {})
    log(f"series {len(S)}")
    recs, drop = build(S, PRIMARY)
    log(f"primary events {len(recs)} dropped {drop}")
    pw = power_block(recs, S, PRIMARY)
    pw.update({"version": VERSION, "code_sha256": code_sha, "written": datetime.now(timezone.utc).isoformat()})
    jdump(pw, os.path.join(a.out, "p8b_power.json"))
    log("power written")
    time.sleep(1.0)
    main_res, recs = run_cfg(S, PRIMARY)
    r_sfc = layer_R(recs, "snake")
    r_sfc["category_P6_threshold"] = (p6.category(*r_sfc["ci95_cluster"]) if r_sfc.get("evaluable") and r_sfc.get("ci95_cluster")
                                      else "不可评")
    r_1m = layer_R(recs, "la1m")
    r_1m["reading"] = reading_1m(r_1m)
    r_2m = layer_R(recs, "tsg5")
    pair12 = [r for r in recs if p6.lay_ok(r, "la1m") and p6.lay_ok(r, "tsg5")]
    q12 = q_block(pair12, "la1m", "tsg5")
    q12["reading"] = reading_12(q12)
    rain = {"R_sfc": p6.rain_split(recs, PRIMARY), "R_1m": p6.rain_split(recs, dict(PRIMARY, sfc="la1m")),
            "Q_L": q_rain_split(recs, PRIMARY)}
    rel2 = p6.strip(p6.rel5_block(recs, PRIMARY))
    sens = {}
    for c in SENS:
        res, rr = run_cfg(S, c)
        b1 = layer_R(rr, "la1m")
        sens[c["name"]] = {"Q_L": {"n": res["Q_L"]["n"], "clusters": res["Q_L"].get("clusters"), "point": res["Q_L"].get("point"),
                                   "ci95_cluster": res["Q_L"].get("ci95_cluster"), "evaluable": res["Q_L"]["evaluable"],
                                   "category": res["judgement"]["category"],
                                   "same_as_primary": res["judgement"]["category"] == main_res["judgement"]["category"]},
                           "R_1m": {"n": b1.get("n"), "point": b1.get("point"), "ci95_cluster": b1.get("ci95_cluster"),
                                    "evaluable": b1.get("evaluable"), "reading": reading_1m(b1)}}
    summ = {"version": VERSION, "code_sha256": code_sha, "data_info": info,
            "primary": main_res, "verdict": main_res["judgement"], "R_sfc_obs": r_sfc, "R_1m_obs": r_1m, "R_2m_obs": r_2m,
            "Q_12_shape": q12, "rain_vs_dry": rain, "R_rel2_sfc_minus_2m": rel2, "sensitivity": sens,
            "disp_km_primary": {"median": rnd(np.nanmedian([r["disp_km"] for r in recs])) if recs else None,
                                "max": rnd(np.nanmax([r["disp_km"] for r in recs])) if recs else None},
            "elapsed_s": round(time.time() - t_start, 1)}
    jdump(summ, os.path.join(a.out, "p8b_summary.json"))
    events_csv(recs, os.path.join(a.out, "p8b_events.csv"))
    log(f"verdict {json.dumps(main_res['judgement'], ensure_ascii=False)}")
    log("done rc=0")
    return 0


# ======================================================================== 自测
def _synth(S, cr, n, H0, rain_h, rho, zs):
    S[f"{cr}|rain_h"] = p6.Series(H0, rain_h)
    S[f"{cr}|rain_hs"] = p6.Series(2 * H0, np.repeat(rain_h, 2))
    S[f"{cr}|rain_h_strict"] = S[f"{cr}|rain_h"]
    S[f"{cr}|rain_hs_strict"] = S[f"{cr}|rain_hs"]
    S[f"{cr}|wind_h"] = p6.Series(H0, np.full(n, 5.0))
    S[f"{cr}|lat_h"] = p6.Series(H0, np.full(n, 12.0))
    S[f"{cr}|lon_h"] = p6.Series(H0, np.full(n, -125.0))
    U = np.full(2 * n, 5.0)
    P = np.repeat(rain_h, 2)
    for ln, z in zs.items():
        F = np.full(n, np.nan)
        F[24:] = p2.rim_factor(P, U, z).reshape(-1, 2).mean(axis=1)
        S[f"{cr}|{ln}"] = p6.Series(H0, 34.0 * (1.0 + rho[ln] * (F - 1.0)))


def selftest():
    res = {}
    res["T0_p6_selftest_all_ok"] = bool(p6.selftest()["all_ok"])
    rng = np.random.default_rng(13)
    P = np.where(rng.random(84) < 0.3, rng.gamma(1.5, 6.0, 84), 0.0)
    U = np.repeat(rng.uniform(1, 11, 42), 2)
    mx = max(float(np.nanmax(np.abs(p6.model_hourly_F(P, U, z) - p2.rim_factor(P, U, z).reshape(-1, 2).mean(axis=1))))
             for z in (0.0, 0.01, 1.0, 2.0))
    res["T1_model_equals_p2_rim_factor"] = mx == 0.0
    ref, _ = p2.model_hourly(P, U, 34.0, {0.5: 0.5, 1.0: 1.0, 5.0: 5.0})
    res["T1b_model_equals_p2_model_hourly"] = max(float(np.nanmax(np.abs(p6.model_hourly_F(P, U, z) - ref[z]))) for z in (0.5, 1.0, 5.0)) == 0.0
    lit = p2._rim_literal(P, U, np.full(len(P), 34.0))
    res["T2_z0_equals_literal"] = bool(np.nanmax(np.abs(np.asarray(lit) / 34.0 - p2.rim_factor(P, U, 0.0))) < 1e-12)
    r = rain_online(np.array([np.nan, 2.0, -1.0, np.nan, 3.0]), np.array([1.0, 1.0, 1.0, np.nan, np.nan]))
    res["T3_rain_online"] = (r[0] == 0.0 and r[1] == 2.0 and r[2] == 0.0 and np.isnan(r[3]) and np.isnan(r[4]))
    res["T4_wind_factor"] = abs(WIND_FACTOR_LA - math.log(1e5) / math.log(18.0 / 1e-4)) < 1e-12 and WIND_FACTOR_LA < 1.0
    n = 24 * 60
    H0 = 400000
    rain_h = np.zeros(n)
    for k in range(60, n - 40, 48):
        rain_h[k:k + 3] = [8.0, 6.0, 4.0]
    zs = {"snake": 0.0, "la1m": 1.0, "tsg5": 2.0}
    for tag, rho, want_q, want in (("R1", {"snake": 0.3, "la1m": 0.3, "tsg5": 0.3}, 1.0, "支持 R1"),
                                   ("R2", {"snake": 1.0, "la1m": 0.3, "tsg5": 0.3}, 1.0 / 0.3, "支持 R2")):
        S = {}
        for cr in ("LA2", "LA5", "LA7"):
            _synth(S, cr, n, H0, rain_h, rho, zs)
        for cr in CRUISES:
            if cr not in ("LA2", "LA5", "LA7"):
                for k in ("rain_h", "rain_hs", "rain_h_strict", "rain_hs_strict", "wind_h", "lat_h", "lon_h", *zs):
                    S[f"{cr}|{k}"] = p6.Series(0, np.array([]))
        out, recs = run_cfg(S, PRIMARY)
        q = out["Q_L"]
        ok = (len(recs) >= 20 and q["evaluable"] and abs(q["point"] - want_q) < 0.03 and out["judgement"]["category"].startswith(want)
              and abs(layer_R(recs, "la1m")["point"] - 0.3) < 0.02)
        res[f"T5_e2e_{tag}"] = bool(ok)
        if tag == "R1":
            q12 = q_block([r for r in recs if p6.lay_ok(r, "la1m") and p6.lay_ok(r, "tsg5")], "la1m", "tsg5")
            res["T5c_Q12_R1_is_1"] = abs(q12["point"] - 1.0) < 0.03 and reading_12(q12).startswith("不能区分")
            qr = q_rain_split(recs, PRIMARY)
            res["T5d_rain_split_runs"] = qr["rain"]["n"] >= 8 and abs(qr["rain"]["point"] - 1.0) < 0.05
    res["T6_q_category"] = (q_category(0.5, 1.7).startswith("支持 R1") and q_category(1.9, 5).startswith("支持 R2")
                            and q_category(1.2, 2.5) == "不能区分" and q_category(None, 1) == "不可评")
    blk = {"n": 20, "clusters": 6, "per_cruise": {"LA5": 12, "LA7": 8}, "frac_boot_rim_a_nonneg": 0.0,
           "frac_boot_rim_b_nonneg": 0.0, "frac_boot_obs_b_nonneg": 0.10}
    res["T7_obs_den_gate"] = q_evaluable(blk)[0] is False and q_evaluable(dict(blk, frac_boot_obs_b_nonneg=0.01))[0] is True
    res["T7b_cruise_gate"] = q_evaluable(dict(blk, per_cruise={"LA5": 18, "LA7": 2}, frac_boot_obs_b_nonneg=0.0))[0] is False
    res["T8_reading_1m"] = (reading_1m({"evaluable": True, "ci95_cluster": [0.1, 0.5]}).startswith("与 P2 同侧")
                            and reading_1m({"evaluable": False, "ci95_cluster": [0.1, 0.5]}) == "不可评")
    res = {k: (bool(v) if isinstance(v, (bool, np.bool_)) else v) for k, v in res.items()}
    res["all_ok"] = all(v is True for v in res.values())
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--load-check", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--root")
    ap.add_argument("--out", default=os.environ.get("REPRO_OUTPUT_DIR"))
    a = ap.parse_args(argv)
    rc = 0
    try:
        if a.selftest:
            st = selftest()
            print(json.dumps(st, ensure_ascii=False, indent=1))
            if a.out:
                os.makedirs(a.out, exist_ok=True)
                jdump(st, os.path.join(a.out, "p8b_selftest.json"))
            rc = 0 if st["all_ok"] else 4
        elif a.load_check or a.full:
            if not (a.root and a.out):
                ap.error("--root、--out 必须给")
            rc = run_load_check(a) if a.load_check else run_full(a)
        else:
            ap.error("need --selftest / --load-check / --full")
    except Exception:
        log("异常：" + traceback.format_exc())
        rc = 3
    if a.out and (a.full or a.load_check):
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, "p8b_log.txt"), "a", encoding="utf-8") as f:
            f.write("\n".join(LOG) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
