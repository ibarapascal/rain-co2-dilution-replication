"""P3c 判门（纯查表，laptop 上跑）：读 p3a_repro.json（须已存在，L23）、p3b_curves.json、P2 p2_summary.json／p2_events.csv，
按 L15/L16/L18/L19/L42 算 M1、S1–S7、E1、型态，判 G-X/G-R/G-N/G-0，写 p3c_decision.json。
用法：python3 p3c_decide.py IN_DIR（含 p3a_repro.json、p3b_curves.json、p2/p2_summary.json、p2/p2_events.csv；输出写回 IN_DIR），省略时用本脚本目录。
P2 bootstrap 复用 p2_rim_test.cluster_boot 同算法（B=10000、seed=20260926、PCG64、站×季簇）。"""
import csv, json, math, os, sys
import numpy as np

D = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
assert os.path.exists(os.path.join(D, "p3a_repro.json")), "L23：p3a_repro.json 不在"
A3 = json.load(open(os.path.join(D, "p3a_repro.json")))
C = json.load(open(os.path.join(D, "p3b_curves.json")))
S = json.load(open(os.path.join(D, "p2", "p2_summary.json")))
EV = list(csv.DictReader(open(os.path.join(D, "p2", "p2_events.csv"))))

c = C["c"]; dep = C["deposition_pct"]; turb = C["turbulence_pct"]; W = C["wind_PgC"]
sg = np.array(C["global"]["whole"]["s"])
ANCH = (10.65, 10.0)


def interp(s, ys, flag=None):
    ys = np.asarray(ys, float)
    if s < 0:
        if flag is not None: flag.append(("clip_s<0", s))
        s = 0.0
    if s > sg[-1]:
        if flag is not None: flag.append(("extrap_s>3", s))
        return float(ys[-1] + (s - sg[-1]) * (ys[-1] - ys[-2]) / (sg[-1] - sg[-2]))
    return float(np.interp(s, sg, ys))


def curve(kind, key):
    return C["global"][kind][key]


U1 = interp(1.0, curve("whole", "U_pp"))


def Dl(s, kind="whole"):
    return interp(s, curve(kind, "Delta_pp"))


def fnum(x):
    try:
        v = float(x); return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def boot_ratio(rows, B=10000, seed=20260926, rng=None):
    keys = [(r["station"], r["season"]) for r in rows]
    uk = sorted(set(keys)); kid = {k: j for j, k in enumerate(uk)}; K = len(uk)
    idx = np.array([kid[k] for k in keys], int)
    sn = np.bincount(idx, weights=np.array([r["_n"] for r in rows]), minlength=K)
    sd = np.bincount(idx, weights=np.array([r["_d"] for r in rows]), minlength=K)
    rng = rng or np.random.Generator(np.random.PCG64(seed))
    draw = rng.integers(0, K, size=(B, K))
    with np.errstate(divide="ignore", invalid="ignore"):
        b = sn[draw].sum(1) / sd[draw].sum(1)
    point = sum(r["_n"] for r in rows) / sum(r["_d"] for r in rows)
    bf = b[np.isfinite(b)]
    return point, [float(np.percentile(bf, 2.5)), float(np.percentile(bf, 97.5))], b, K


def rows_for(numc, denc, pred=lambda r: True):
    out = []
    for r in EV:
        a, d = fnum(r[numc]), fnum(r[denc])
        if a is not None and d is not None and pred(r):
            rr = dict(r); rr["_n"] = a; rr["_d"] = d; out.append(rr)
    return out


out = {"reads": {"p3a_all_pass": A3["all_pass"], "p2_summary_exit": S["exit"]["exit"]}, "U1_curve": U1}

# ---- M1 + L16
m1 = S["primary"]["R_RIM_1m"]; s_pt, s_lo, s_hi = m1["point"], m1["ci95"][0], m1["ci95"][1]
r1 = rows_for("dS1_obs", "dS_rim_1")
p_re, ci_re, _, K1 = boot_ratio(r1)
l16 = {"n": len(r1), "K": K1, "point_re": p_re, "ci_re": ci_re, "ci_p2": m1["ci95"],
       "max_abs_diff": max(abs(ci_re[0] - s_lo), abs(ci_re[1] - s_hi)),
       "alarm_gt_0.01": max(abs(ci_re[0] - s_lo), abs(ci_re[1] - s_hi)) > 0.01, "used": "P2 ci95（L16：以 P2 为准）"}
mono = C["global"]["whole"]["monotone"]
d_pt, d_lo, d_hi = Dl(s_pt), Dl(s_lo), Dl(s_hi)
dci = sorted([d_lo, d_hi])
M1 = {"s": s_pt, "s_ci": [s_lo, s_hi], "monotone": mono, "Delta": d_pt, "Delta_ci": dci,
      "U_curve": U1 + d_pt, "U_curve_ci": [U1 + dci[0], U1 + dci[1]],
      "U_new": {str(a): [a + d_pt, [a + dci[0], a + dci[1]]] for a in ANCH},
      "PgC_uptake_U": (U1 + d_pt) / 100 * abs(W)}
out["M1"] = M1; out["L16"] = l16

# P3-0 速查
p30 = lambda s: 3.978 + 6.671 * s
out["P3_0"] = {"U": p30(s_pt), "U_ci": [p30(s_lo), p30(s_hi)], "Delta": p30(s_pt) - p30(1), "note": "未作判门依据（P3a 过）"}

# ---- 判门
def in_out(lo, hi):
    return hi < 5 or lo > 15


GX = all(in_out(a + dci[0], a + dci[1]) for a in ANCH)
GR = dci[1] < -2.5 or dci[0] > 2.5
GN = all(5 <= a + dci[0] and a + dci[1] <= 15 for a in ANCH) and (dci[1] - dci[0]) <= 5
gate = "G-X" if GX else "G-R" if GR else "G-N" if GN else "G-0"
out["gates"] = {"G-X": GX, "G-R": GR, "G-N_info": GN, "hit": gate}

# ---- 稳健性集合
# S1
r05 = rows_for("dS05_obs", "dS_rim_05")
s1p, s1ci, _, _ = boot_ratio(r05)
# S2：A/B 各自簇，同一 PCG64(seed) 生成器先 A 后 B（独立流）；另报「两群各用 seed」变体
isA = lambda r: r["group"].startswith("A")
rA = [r for r in r1 if isA(r)]; rB = [r for r in r1 if not isA(r)]
g = np.random.Generator(np.random.PCG64(20260926))
sA, sAci, bA, KA = boot_ratio(rA, rng=g)
sB, sBci, bB, KB = boot_ratio(rB, rng=g)
_, sAci_same, bA2, _ = boot_ratio(rA)
_, sBci_same, bB2, _ = boot_ratio(rB)
BD = C["bands"]; lo = np.array(BD["lat_lo"]); hi = np.array(BD["lat_hi"])
Iband = np.array(BD["whole_I_int_PgC"])  # [band][s]


def U_split(sa, sb, cut):
    trop = (lo >= -cut) & (hi <= cut)
    tot = 0.0
    for j in range(len(lo)):
        tot += interp(sa if trop[j] else sb, Iband[j])
    I_pct = -tot / abs(W) * 100
    return dep + I_pct * (1 - c)


chk_I1 = -Iband.sum(0) / abs(W) * 100
band_check = {"I_int_pct_from_bands_s1": float(interp(1.0, chk_I1)), "global": interp(1.0, curve("whole", "I_int_pct"))}
S2 = {}
for cut in (20, 15, 25):
    up = U_split(sA, sB, cut)
    ub = np.array([U_split(a, b, cut) for a, b in zip(bA[:10000], bB[:10000])]) if cut == 20 else None
    ent = {"U": up, "Delta": up - U1}
    if ub is not None:
        ub = ub[np.isfinite(ub)]
        ent["U_ci"] = [float(np.percentile(ub, 2.5)), float(np.percentile(ub, 97.5))]
        ub2 = np.array([U_split(a, b, cut) for a, b in zip(bA2, bB2)])
        ent["U_ci_same_seed_variant"] = [float(np.percentile(ub2, 2.5)), float(np.percentile(ub2, 97.5))]
    S2[str(cut)] = ent
k27 = S["exploratory_k27"]["R_RIM_1.0m"]
rob = {
    "S1": {"s": s1p, "s_ci": s1ci, "n": len(r05), "Delta": Dl(s1p), "Delta_ci": sorted([Dl(s1ci[0]), Dl(s1ci[1])])},
    "S2": {"s_A": sA, "s_A_ci": sAci, "s_A_ci_same_seed": sAci_same, "K_A": KA, "s_B": sB, "s_B_ci": sBci, "s_B_ci_same_seed": sBci_same, "K_B": KB,
           "Delta": S2["20"]["Delta"], "by_cut": S2, "band_check": band_check},
    "S3": {"s": s_pt, "Delta": Dl(s_pt, "hist"), "Delta_ci": sorted([Dl(s_lo, "hist"), Dl(s_hi, "hist")])},
    "E1": {"s": k27["point"], "s_ci": k27["ci95"], "Delta": Dl(k27["point"]), "Delta_ci": sorted([Dl(k27["ci95"][0]), Dl(k27["ci95"][1])])},
}
sign = np.sign(d_pt)
n_ok = sum(1 for v in rob.values() if abs(v["Delta"]) >= 2.5 and np.sign(v["Delta"]) == sign)
n_opp = sum(1 for v in rob.values() if np.sign(v["Delta"]) == -sign and v["Delta"] != 0)
robust = n_ok >= 3 and n_opp == 0
out["robustness"] = {"set": rob, "n_meet_GR_point": n_ok, "n_opposite": n_opp, "robust": robust}

# ---- 型态
REG = {"暖池": ["TAO165E"], "SPCZ": ["TAO8S165E"], "孟加拉湾": ["BOBOA"],
       "中东太平洋": ["TAO170W", "TAO155W", "TAO140W", "TAO125W", "TAO110W"], "中纬／副热带": None}
Rp = s_pt
reg = {}
for k, st in REG.items():
    rr = [r for r in r1 if (not isA(r) if st is None else r["station"] in st)]
    n = len(rr)
    R = sum(r["_n"] for r in rr) / sum(r["_d"] for r in rr) if n else None
    counted = n >= 10
    reg[k] = {"n_valid": n, "R": R, "counted": counted,
              "GXR_ok": bool(counted and (R - 1) * (Rp - 1) > 0 and abs(R - 1) >= 0.2),
              "GN_ok": bool(counted and 0.65 <= R <= 1.35)}
n_shape = sum(v["GXR_ok"] for v in reg.values())
shape_ok = n_shape >= 3
out["regimes"] = {"by_regime": reg, "n_counted": sum(v["counted"] for v in reg.values()), "n_GXR_ok": n_shape, "shape_ok_for_GXR": shape_ok}

# ---- 信息项 S4–S7
s20 = C["global"]["s20_1"]
I0 = interp(0.0, curve("whole", "I_int_pct"))
Rs20 = S["primary"]["R_S20_05m"]
def U_S5(R, base): return dep + (base + R * (s20["I_int_pct"] - base)) * (1 - c)
beta = S["primary"]["beta_D4subset"]; rb = beta["point"] / beta["beta_model_pCO2_pooled"]
def U_S6(s, r):
    I = interp(s, curve("whole", "I_int_pct")) + (r - 1) * interp(s, curve("whole", "G_int_pct"))
    return dep + I * (1 - c)
P5 = S["primary"]["P5"]; q = P5["point"] / P5["P5_RIM3"]["point"]; cq = min(max(0.30 * q, 0), 1)
def U_S7(s, cc): return dep + interp(s, curve("whole", "I_int_pct")) * (1 - cc)
out["info"] = {
    "S4": {"U": s20["U_pp"], "Delta": s20["U_pp"] - U1},
    "S5_L42": {"R": Rs20["point"], "R_ci": Rs20["ci95"], "I0_not_scaled": I0, "U": U_S5(Rs20["point"], I0),
               "U_ci": [U_S5(Rs20["ci95"][0], I0), U_S5(Rs20["ci95"][1], I0)]},
    "S5_scale_with_skin_explor": {"base_not_scaled": turb, "U": U_S5(Rs20["point"], turb),
                                  "U_ci": [U_S5(Rs20["ci95"][0], turb), U_S5(Rs20["ci95"][1], turb)]},
    "S6": {"r_beta": rb, "beta_obs": beta["point"], "beta_model_pCO2": beta["beta_model_pCO2_pooled"],
           "U": U_S6(s_pt, rb), "U_ci_sCI": sorted([U_S6(s_lo, rb), U_S6(s_hi, rb)]),
           "explor_U_rbeta_ci": sorted([U_S6(s_pt, beta["ci95"][0] / beta["beta_model_pCO2_pooled"]), U_S6(s_pt, beta["ci95"][1] / beta["beta_model_pCO2_pooled"])])},
    "S7": {"q": q, "c_prime": cq, "U": U_S7(s_pt, cq), "U_ci": sorted([U_S7(s_lo, cq), U_S7(s_hi, cq)])},
}
# ---- L30（探索性）
skin_nonD = C["skin_term_nonD_pct"]; dil0 = interp(0.0, curve("whole", "I_dil_pct")); dil1 = interp(1.0, curve("whole", "I_dil_pct"))
out["L30_explor"] = {"skin_nonD_pct": skin_nonD, "dilution_at_s0_pct": dil0, "interfacial_at_s0_pct": I0,
                     "share_of_dilution_s1_global_const": dil0 / dil1, "share_of_dilution_s1_nonD": skin_nonD / dil1,
                     "U1_minus_skin": U1 - dil0 * (1 - c), "U_M1_minus_skin": U1 + d_pt - dil0 * (1 - c),
                     "curve_slope_dU_ds_0_1": interp(1, curve("whole", "U_pp")) - interp(0, curve("whole", "U_pp")),
                     "P3_0_slope": 6.671}
out["non_gate"] = {k: C["global"][kind][k2] for k, kind, k2 in
                   [("whole_monotone", "whole", "monotone"), ("hist_monotone", "hist", "monotone"),
                    ("whole_dev_P3_0", "whole", "max_abs_dev_from_P3_0_line_pp"), ("hist_dev_P3_0", "hist", "max_abs_dev_from_P3_0_line_pp"),
                    ("whole_dev_own", "whole", "max_abs_dev_from_own_endpoints_line_pp"), ("hist_dev_own", "hist", "max_abs_dev_from_own_endpoints_line_pp")]}
json.dump(out, open(os.path.join(D, "p3c_decision.json"), "w"), ensure_ascii=False, indent=1, default=float)
print(json.dumps(out, ensure_ascii=False, default=float)[:200])
