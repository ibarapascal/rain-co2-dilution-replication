#!/usr/bin/env python3
"""fill_upstream_sha.py — 可选工具（只在配置 [options] check_upstream_sha = true 时需要）：列出脚本里写死的「上游输入 sha 闸门」在本次 out_root 上的实际值，供操作者核过上游产物后
粘贴进配置 [upstream_sha]（只读文件、只打印，不改配置、不做任何计算）。

闸门（缺省＝参考运行的上游产物 sha，由 repro_paths.expect_sha 读取）：P7d 读 P5；P8a 读 P5／P7a／P7c；C1/C2（cal-surface-diag）读
P5／P6／P8a／P7e／P8b；dw33-post 读 P4／P3B／P3C；dw13a 读 D13M。闸门缺省关闭（check_upstream_sha = false：只记录实际 sha，不停）；打开时用来防止上游产物漂移：
上游一变这些阶段就退出 3——这是「人工确认上游新产物后再放行」的检查点，不要自动放行。
每行输出：新值与缺省值相同（前缀匹配）则注释掉（不必覆盖）；不同则给出可粘贴的 TOML 行，并标出所依赖的上游阶段。
用法：python3 tools/fill_upstream_sha.py OUT_ROOT [--only p8a,cal12]
Change Log：2026-09-29 初版；2026-09-30 在全量重跑的每个闸门上用过；2026-10-01 闸门改为缺省关闭，本工具变为可选。
"""

import hashlib
import os
import sys

# 键 → (上游阶段, 文件, 缺省值或前缀)
GATES = {
    "p7d.p5_events_csv": ("p5-sss", "p5_events.csv", "7a6d70758d575d180d9f0100579a9b7f6bc014350c9080ca9b5ab84ffd6985fa"),
    "p7d.p5_summary_json": ("p5-sss", "p5_summary.json", "e5f66e65fb13d812683eee58c0f9182906f51f1cd3356247bd284bd00385ab4f"),
    "p8a.p5_events_csv": ("p5-sss", "p5_events.csv", "7a6d70758d575d180d9f0100579a9b7f6bc014350c9080ca9b5ab84ffd6985fa"),
    "p8a.p5_summary_json": ("p5-sss", "p5_summary.json", "e5f66e65fb13d812683eee58c0f9182906f51f1cd3356247bd284bd00385ab4f"),
    "p8a.p7a_events_jpl_csv": ("p7a", "p7a_events_jpl.csv", "1b529937eb13dbd02e79d7b66e22ed08810471a6f2087e30c3c559903c05f6a8"),
    "p8a.p7a_summary_json": ("p7a", "p7a_summary.json", "8e030de8599390dd65600657c67ac4190346264e62d5f28bd55dbc2698453801"),
    "p8a.p7c_events_csv": ("p7c", "p7c_events.csv", "247f05856bdeb73ad8261dba8009929298f5f77ea9fc66cf6901739d1b0ab4db"),
    "cal12.p5_events_csv": ("p5-sss", "p5_events.csv", "7a6d70758d575d180d9f0100579a9b7f6bc014350c9080ca9b5ab84ffd6985fa"),
    "cal12.p6_events_csv": ("p6", "p6_events.csv", "2de4313d"),
    "cal12.p8a_events_wind_csv": ("p8a", "p8a_events_wind.csv", "b8952fc6"),
    "cal12.p7e_spurs1_events_csv": ("p7e", "p7e_spurs1_events.csv", "7c2e0a56"),
    "cal12.p8b_events_csv": ("p8b", "p8b_events.csv", "ddba60b7"),
    "dw33post.p4e": ("p4-mech", "p4_mech_events.csv", "c5feedaa"),
    "dw33post.p4s": ("p4-mech", "p4_mech_summary.json", "e3e2bf7b"),
    "dw33post.p3b": ("p3-global", "p3b_curves.json", "2bdbe143"),
    "dw33post.p3c": ("p3c", "p3c_decision.json", "8a88e222"),
    "dw13a.dw13_meta": ("cal-dw13-meta", "dw13_meta.json", "fd5a036b"),   # 2026-09-29c
}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def report(root, only=None):
    lines = ["[upstream_sha]  # 由 tools/fill_upstream_sha.py 列出；逐条核过上游新产物再粘贴"]
    for key, (stage, fn, frozen) in GATES.items():
        if only and key.split(".", 1)[0] not in only:
            continue
        p = os.path.join(root, stage, fn)
        if not os.path.exists(p):
            lines.append(f"# {key}：上游 {stage}/{fn} 尚不存在")
            continue
        s = sha256(p)
        if s.startswith(frozen):
            lines.append(f"# \"{key}\" 与缺省值相同（{frozen[:12]}…），不必覆盖")
        else:
            lines.append(f"\"{key}\" = \"{s}\"  # 上游 {stage}/{fn}；缺省 {frozen[:12]}…")
    return lines


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    only = None
    if "--only" in sys.argv:
        only = set(sys.argv[sys.argv.index("--only") + 1].split(","))
    print("\n".join(report(sys.argv[1], only)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
