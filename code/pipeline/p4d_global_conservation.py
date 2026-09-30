"""P4d 模型推导（秒级、只读 P3 小产物）：把 RIM-3 按淡水守恒缩放后，全球 U(s) 的均匀因子近似。

背景：P4b 描述项／P4c 解析式 M＝c1·√π·(tc/1800)/d0 表明 Witte 式 (4) 在系泊事件上多造约 4.46 倍淡水。
本脚本在 P3b 实算曲线（p3b_curves.json，s 网格 0–3，分段线性插值同 L15）上查：
  - s_c＝1/M_RIM（M_RIM＝4.46 [4.26, 4.68]，P4b `p4b_summary.json`）与 p4c 风速档范围 s∈[1/10, 1/1.7]；
  - 与 U(1)、U(M1＝0.29616 [0.22631, 0.36672]) 的差；稀释分量 U(s)−U(0) 的相对变化；
  - whole 曲线（S_dil＝S0＋s·ΔS_RIM，当前项一起缩放）为主，hist 曲线（只缩放历史项）为信息；
  - d0 表逐格的解析 M（w＝1、0.95）→ s＝1/M → U(s)，只作量级参照。
口径：曲线在任何 s 下都不缩放 L30 皮层常数项（它在 s=0 截距里），与 L42（S5 不缩皮层常数项）一致。
限定：均匀 s 近似风速依赖的守恒因子（真值随格点 U、R、雨龄 w 变）；精确版需按格点缩放全球重跑（约 7 h，未实现）。
用法：python3 p4d_global_conservation.py P3B_JSON [OUT_JSON]；仅标准库＋numpy。
"""
import json
import sys

import numpy as np

C = json.load(open(sys.argv[1]))
OUT = sys.argv[2] if len(sys.argv) > 2 else None
sg = np.array(C["s_grid"], float)

M_RIM, M_RIM_CI = 4.46, (4.26, 4.68)          # P4b 细网格淡水/雨量（580 事件，CMORPH 驱动）
M1, M1_CI = 0.29616, (0.22631, 0.36672)       # P2 R_RIM(1 m)，K12 bootstrap
M_RANGE = (10.0, 1.7)                          # p4c 风速档（d0 表 1.0–5.8 m，w=1）
C1, W_EFF = 5.7, 0.95


def U(s, kind="whole"):
    assert 0.0 <= s <= sg[-1]
    return float(np.interp(s, sg, C["global"][kind]["U_pp"]))


def row(s, kind="whole"):
    u, u0, u1, um1 = U(s, kind), U(0.0, kind), U(1.0, kind), U(M1, kind)
    return {"s": s, "U": u, "minus_U1": u - u1, "minus_UM1": u - um1,
            "dil_component": u - u0, "dil_component_rel_change": (u - u0) / (u1 - u0) - 1}


out = {"inputs": {"M_RIM": M_RIM, "M_RIM_ci": M_RIM_CI, "M1": M1, "M1_ci": M1_CI, "M_range": M_RANGE},
       "U0_whole": U(0.0), "U1_whole": U(1.0), "UM1_whole": U(M1), "UM1_ci_whole": [U(M1_CI[0]), U(M1_CI[1])]}
for kind in ("whole", "hist"):
    sc = 1 / M_RIM
    out[kind] = {
        "s_c": row(sc, kind),
        "s_c_ci_from_M_RIM": [row(1 / M_RIM_CI[1], kind), row(1 / M_RIM_CI[0], kind)],
        "wind_bin_bounds": [row(1 / M_RANGE[0], kind), row(1 / M_RANGE[1], kind)],
        "M1_row": row(M1, kind),
    }
# 系泊量级闭合：若守恒缩放后的模型是对的，1 m 观测/守恒模型＝M1×M_RIM
out["closure_M1_times_MRIM"] = M1 * M_RIM
# d0 表逐格解析 M（胶囊第 41–50 行实质格点 U∈{2..10}×R∈{2..50}）
DL = {2: [1.2, 1.0, 1.0, 1.0, 1.1], 4: [1.6, 1.9, 2.0, 2.0, 2.2], 6: [1.5, 2.5, 2.9, 3.1, 3.5],
      8: [1.9, 2.7, 3.6, 4.2, 4.7], 10: [2.4, 3.0, 4.0, 5.0, 5.8]}
RR = [2, 5, 10, 20, 50]
tab = []
for u10, d0s in DL.items():
    for r, d0 in zip(RR, d0s):
        m = C1 * np.sqrt(np.pi) * W_EFF / d0
        s = min(1 / m, sg[-1])
        tab.append({"U10": u10, "R": r, "d0": d0, "M_w095": m, "s": s, "U_whole": U(s)})
out["d0_table_analytic"] = tab

for kind in ("whole", "hist"):
    o = out[kind]
    print(f"[{kind}] U(0)={U(0.0, kind):.3f} U(1)={U(1.0, kind):.3f} U(M1)={U(M1, kind):.3f}")
    for lab, r in [("s_c", o["s_c"]), ("s_c(M_RIM hi)", o["s_c_ci_from_M_RIM"][0]), ("s_c(M_RIM lo)", o["s_c_ci_from_M_RIM"][1]),
                   ("s=1/10", o["wind_bin_bounds"][0]), ("s=1/1.7", o["wind_bin_bounds"][1])]:
        print(f"  {lab:14s} s={r['s']:.4f} U={r['U']:.3f} dU1={r['minus_U1']:+.3f} dUM1={r['minus_UM1']:+.3f} "
              f"dil={r['dil_component']:.3f} ({r['dil_component_rel_change']*100:+.1f}%)")
print(f"M1*M_RIM={out['closure_M1_times_MRIM']:.3f}")
for t in tab:
    print(f"  U10={t['U10']:2d} R={t['R']:2d} d0={t['d0']:.1f} M={t['M_w095']:.2f} s={t['s']:.3f} U={t['U_whole']:.3f}")
if OUT:
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
