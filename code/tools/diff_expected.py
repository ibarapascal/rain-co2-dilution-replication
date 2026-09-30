#!/usr/bin/env python3
"""diff_expected.py — 重跑之后用：把 <out_root>/<阶段>/ 的新产物与本 repo 发布的数据（data/products/）逐项对照
（只读新产物，不做任何计算；只用标准库）。

对 expected/products.tsv 的每行：
  新文件是否存在、sha256 是否与参考运行的产物 sha（前缀，或「前缀…后缀」）一致 → same；
  也接受 data/products/MANIFEST.tsv 里发布副本的 sha → same(released-copy)（用于对发布副本自比）；
  compare=rows 的 CSV（行序随并发而变）改比「全部行排序后的 sha」（rows_sha256）→ same(rows-multiset)；
  compare=json 的 JSON：sha 不同时与基准副本（缺省 data/products/，--baseline 可改）逐叶比较——数值 |新−旧| ≤ 1e-12＋1e-9·|旧|，
    字符串与键按 docs/translation-table.tsv 中英互认；跳过 ignore 列的叶路径（逗号分隔，fnmatch，如 .inputs.*.path；运行时间、
    代码 sha、输入路径、下载计数等运行记录，以及参考运行输出里键名或措辞与现行脚本不同的文字字段；不忽略任何数值字段）→ same(json-leaves)；
    其余任何叶不同 → differs（计入失败）；没有基准副本时按下一条处理；
  sha 不同且 expect_identical 以 no 开头 → differs(expected)，否则 differs（计入失败）；
  阶段不属于 run.py --all（stages.py 的 in_all 为假，如 dw33-bench）且阶段目录不存在 → not-run（不计失败）；
  对 data/products 自比时，未发布的文件（日志、自测等不在 MANIFEST）→ not-released。
键名本身含点（R_RIM_1.0m、z=1.0m）时按「最长匹配」识别；中文键（代码输出）与发布副本的英文键按翻译表双向互认。
用法：python3 tools/diff_expected.py OUT_ROOT [--report report.tsv] [--baseline DIR] [--present-only]
      --present-only：只比 OUT_ROOT 里已有阶段目录的行（重跑中途的早期比对用）；
      自比：python3 tools/diff_expected.py ../data/products（发布副本对自身，应 0 失败）。
退出码：0 无失败项（只有 same／not-run／not-released／预期不同）；1 有文件缺失或非预期的差异。
Change Log：2026-09-27 初版。
            2026-09-30 带点键名（最长匹配）、* 通配、中英键互认；compare=rows（行多重集）；接受发布副本 sha；--present-only。
            2026-10-01 compare=json（与发布副本逐叶比、ignore 列跳过运行记录字段）；非 in_all 阶段缺失记 not-run；
                       去掉按文中数字表的数值比对，JSON 产物一律逐叶比。
"""

import csv
import fnmatch
import hashlib
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.normpath(os.path.join(HERE, "..", "expected"))
PRODUCTS = os.path.normpath(os.path.join(HERE, "..", "..", "data", "products"))
TRANSLATION = os.path.normpath(os.path.join(HERE, "..", "..", "docs", "translation-table.tsv"))
BAD = {"missing", "differs"}


def rows(path):
    with open(path, encoding="utf-8") as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rows_sha256(p):
    """全部行（含表头，按字节）排序后的 sha256：只看行多重集，不看行序。"""
    with open(p, "rb") as f:
        lines = f.read().splitlines()
    h = hashlib.sha256()
    for ln in sorted(lines):
        h.update(ln + b"\n")
    return h.hexdigest()


def sha_match(full, want):
    if not want:
        return False
    if "…" in want:
        a, b = want.split("…", 1)
        return full.startswith(a) and full.endswith(b)
    return full.startswith(want)


def load_translation():
    alias = {}                                # 一个中文串可有多个译法（如 JSON 键与正文），双向登记
    if os.path.exists(TRANSLATION):
        for r in rows(TRANSLATION):
            alias.setdefault(r["original_zh"], []).append(r["released_en"])
            alias.setdefault(r["released_en"], []).append(r["original_zh"])
    return alias


ALIAS = load_translation()


def same_text(a, b):
    return a == b or b in ALIAS.get(a, ()) or any(b in ALIAS.get(z, ()) for z in ALIAS.get(a, ()))


def lookup(d, key):
    if key in d:
        return key
    for k in d:
        if isinstance(k, str) and same_text(key, k):
            return k
    return None


def leaf_diffs(a, b, path="", out=None, rtol=1e-9, atol=1e-12, ignore=()):
    """逐叶比较（新 a 对基准 b）：键与字符串按翻译表互认，数值 |a−b| ≤ atol＋rtol·|b|；
    ignore＝叶路径的 fnmatch 模式（如 .inputs.*.path），命中的子树两侧都不比（含只在一侧出现）。"""
    out = [] if out is None else out
    skip = lambda p: any(fnmatch.fnmatchcase(p, g) for g in ignore)  # noqa: E731
    if isinstance(a, dict) and isinstance(b, dict):
        seen = set()
        for k in a:
            kb = lookup(b, k)
            if kb is not None:
                seen.add(kb)
            if skip(f"{path}.{k}"):
                continue
            if kb is None:
                out.append(f"{path}.{k}: only in new")
            else:
                leaf_diffs(a[k], b[kb], f"{path}.{k}", out, rtol, atol, ignore)
        out += [f"{path}.{k}: only in baseline" for k in b if k not in seen and not skip(f"{path}.{k}")]
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: length {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            leaf_diffs(x, y, f"{path}[{i}]", out, rtol, atol, ignore)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
        if not (a == b or (math.isnan(a) and math.isnan(b)) or abs(a - b) <= atol + rtol * abs(b)):
            out.append(f"{path}: {a} vs {b}")
    elif not (a == b or (isinstance(a, str) and isinstance(b, str) and same_text(a, b))):
        out.append(f"{path}: {str(a)[:40]!r} vs {str(b)[:40]!r}")
    return out


def load_json(p, cache={}):
    if p not in cache:
        with open(p, encoding="utf-8") as f:
            cache[p] = json.load(f)
    return cache[p]


def released_shas():
    p = os.path.join(PRODUCTS, "MANIFEST.tsv")
    return {r["file"]: r["released_sha256"] for r in rows(p)} if os.path.exists(p) else {}


def not_in_all():
    """stages.py 里 in_all 为假的阶段（run.py --all 不跑）；读不到阶段表时返回空集（全部按 in_all 处理）。"""
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..")))
    try:
        import stages
    except ImportError:
        return set()
    finally:
        sys.path.pop(0)
    return {s["id"] for s in stages.STAGES if not s.get("in_all")}


def compare_files(root, present_only, baseline=PRODUCTS):
    out, rel, skip_all = [], released_shas(), not_in_all()
    self_cmp = os.path.isdir(root) and os.path.samefile(root, PRODUCTS)
    for r in rows(os.path.join(EXP, "products.tsv")):
        if present_only and not os.path.isdir(os.path.join(root, r["stage"])):
            continue
        name = f"{r['stage']}/{r['file']}"
        p = os.path.join(root, r["stage"], r["file"])
        if not os.path.exists(p):
            if r["stage"] in skip_all and not os.path.isdir(os.path.join(root, r["stage"])):
                state = "not-run"                 # 不属于 run.py --all 的阶段，没跑不算缺
            else:
                state = "not-released" if self_cmp and name not in rel else "missing"
            out.append(("file", r["stage"], r["file"], state, ""))
            continue
        s, why = sha256(p), ""
        if sha_match(s, r["sha256"]):
            state = "same"
        elif rel.get(name) == s:
            state = "same(released-copy)"
        elif r.get("compare") == "rows":
            want = r.get("rows_sha256", "")
            state = "manual(rows)" if not want else ("same(rows-multiset)" if rows_sha256(p) == want else "differs")
        elif r.get("compare") == "json" and os.path.exists(os.path.join(baseline, r["stage"], r["file"])):
            ign = [g.strip() for g in (r.get("ignore") or "").split(",") if g.strip()]
            d = leaf_diffs(load_json(p), load_json(os.path.join(baseline, r["stage"], r["file"])), ignore=ign)
            state = "differs" if d else "same(json-leaves)"
            why = "; ".join(d[:3]) or ("ignored: " + ",".join(ign) if ign else "")
        else:
            state = "differs(expected)" if r["expect_identical"].startswith("no") else "differs"
        out.append(("file", r["stage"], r["file"], state, s[:16] + (f" {why}" if why else "")))
    return out


def main():
    if len(sys.argv) < 2 or sys.argv[1].startswith("--"):
        print(__doc__)
        return 2
    root = sys.argv[1]
    report = sys.argv[sys.argv.index("--report") + 1] if "--report" in sys.argv else None
    baseline = sys.argv[sys.argv.index("--baseline") + 1] if "--baseline" in sys.argv else PRODUCTS
    present_only = "--present-only" in sys.argv
    out = compare_files(root, present_only, baseline)
    for o in out:
        print("\t".join(str(x) for x in o))
    if report:
        with open(report, "w", encoding="utf-8") as f:
            f.write("type\tkey\tfile\tstate\tdetail\n")
            for o in out:
                f.write("\t".join(str(x) for x in o) + "\n")
    cnt = {}
    for o in out:
        cnt[f"{o[0]}:{o[3]}"] = cnt.get(f"{o[0]}:{o[3]}", 0) + 1
    print(f"# 汇总 {dict(sorted(cnt.items()))}", file=sys.stderr)
    return 1 if any(o[3] in BAD for o in out) else 0


if __name__ == "__main__":
    sys.exit(main())
