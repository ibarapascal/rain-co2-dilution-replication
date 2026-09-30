"""test_public.py — static self-tests of the pipeline (standard library only; no network; seconds; reads only the small
released products in data/products/).

Covers: no hard-coded data paths in code; stage table complete (scripts exist, dependencies acyclic, every mode expands without
placeholders under a temporary config, no descriptive-only fields other than ref_run); repro_paths placeholders, derivation and
overrides; minimal TOML parser equals tomllib; switches (anemometer height, input-SHA gates off by default); module_sha (SHA-256 of the pipeline file an importing
stage checks); all .py files compile; every stage's completion markers are written by its script; the released products compare
as identical to themselves; comparison-tool rules on modified copies of released files; the raw archive is never overwritten.
Not covered: the stages' own --selftest modes (synthetic data, but need numpy etc.; see docs/reproduction.md).
Usage: cd code && python3 -m unittest tests/test_public.py -v
Change Log: 2026-09-27 first version; 2026-09-30 done-marker and expected-table self-compare tests; 2026-10-01 raw archive
            never overwritten, comparison-tool rules (json leaves, not-run); module_sha replaces the assembly manifest check;
            2026-10-01b input-SHA gates: off by default ([options] check_upstream_sha), on = stop as before.
"""

import ast
import glob
import hashlib
import json
import os
import py_compile
import re
import subprocess
import sys
import tempfile
import unittest

REPRO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # = code/
PIPE = os.path.join(REPRO, "pipeline")
sys.path.insert(0, REPRO)
sys.path.insert(0, PIPE)

import stages as S  # noqa: E402
import repro_paths as rp  # noqa: E402

CFG_TEXT = """[paths]
raw_root = "/tmp/repro-t/raw"
fast_root = "/tmp/repro-t/fast"
out_root = "{fast_root}/runs/T1"
token_file = "/tmp/repro-t/token-not-read"
p3_cache = "{fast_root}/p3-cache"
p3_work = "{fast_root}/p3-work"
p5_cache = "{fast_root}/p5-cache"
p6_cache = "{fast_root}/p6-cache"
p7a_cache = "{fast_root}/p7a-cache"
raw_selftest_root = "/tmp/repro-t/raw-selftest"
dw33_work = "{fast_root}/dw33-work"

[options]
raw_rate_mbps = 2

[upstream]
p2_dir = "/hist/p2"

[upstream_sha]
"p7d.p5_events_csv" = "abc"
"""


def read(p):
    with open(p, encoding="utf-8") as f:
        return f.read()


def write_cfg(d):
    p = os.path.join(d, "cfg.toml")
    with open(p, "w", encoding="utf-8") as f:
        f.write(CFG_TEXT)
    return p


def docstring_nodes(tree):
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(getattr(first, "value", None), ast.Constant) \
                    and isinstance(first.value.value, str):
                out.add(id(first.value))
    return out


def load_json(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


class TestPipelineStatic(unittest.TestCase):
    def pyfiles(self):
        return [p for p in glob.glob(os.path.join(PIPE, "*.py"))]

    def test_no_hardcoded_data_paths(self):
        bad = tuple(os.sep + x for x in ("Volumes" + os.sep + "bulk", "Volumes" + os.sep + "fast", "Users" + os.sep, "private" + os.sep + "tmp"))
        for p in self.pyfiles():
            tree = ast.parse(read(p))
            docs = docstring_nodes(tree)
            for node in ast.walk(tree):
                if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docs:
                    for b in bad:
                        self.assertNotIn(b, node.value, f"{os.path.basename(p)}:{node.lineno} still hard-codes {b}")

    def test_raw_delete_paths_removed(self):
        for p in self.pyfiles():
            tree = ast.parse(read(p))
            docs = docstring_nodes(tree)
            for node in ast.walk(tree):
                if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docs:
                    self.assertNotEqual(node.value, "deleted_after_derive", os.path.basename(p))
        with open(os.path.join(PIPE, "p2_rim_test.py"), encoding="utf-8") as f:
            p2 = f.read()
        code_only = "\n".join(ln.split("#", 1)[0] for ln in p2.splitlines())
        self.assertNotIn("os.remove(val)", code_only)
        self.assertIn("_rio.cmorph_fetcher(DataSourceError, log)", p2)
        self.assertNotIn("urllib.request.urlopen(req, timeout=300)", code_only)

    def test_upstream_keys_known(self):
        for p in self.pyfiles():
            with open(p, encoding="utf-8") as f:
                txt = f.read()
            for k in re.findall(r'_rp\.upstream\("([a-z0-9_]+)"\)', txt):
                self.assertIn(k, rp.UPSTREAM_DEFAULTS, f"{os.path.basename(p)} uses an undefined upstream key {k}")

    def test_py_compile_all(self):
        files = glob.glob(os.path.join(REPRO, "**", "*.py"), recursive=True)
        self.assertGreater(len(files), 30)
        with tempfile.TemporaryDirectory() as td:
            for i, f in enumerate(files):
                py_compile.compile(f, cfile=os.path.join(td, f"{i}.pyc"), doraise=True)


class TestPathsAndStages(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.cfg = write_cfg(self.td.name)
        rp.load(self.cfg, force=True)

    def tearDown(self):
        rp._CFG = None
        self.td.cleanup()

    def test_placeholder_when_unconfigured(self):
        rp._CFG = {}
        self.assertTrue(rp.path("raw_root").startswith(rp.PLACEHOLDER))
        self.assertTrue(rp.upstream("p2_dir").startswith(rp.PLACEHOLDER))

    def test_derivation_and_override(self):
        self.assertEqual(rp.path("out_root"), "/tmp/repro-t/fast/runs/T1")
        self.assertEqual(rp.upstream("p1_events"), "/tmp/repro-t/fast/runs/T1/p1/p1_events.csv")
        self.assertEqual(rp.upstream("p2_dir"), "/hist/p2")
        self.assertEqual(rp.upstream("p4_events"), "/tmp/repro-t/fast/runs/T1/p4-mech/p4_mech_events.csv")
        self.assertEqual(rp.raw_sub("spurs2"), "/tmp/repro-t/raw/spurs2")
        self.assertEqual(rp.expect_sha("p7d.p5_events_csv", "default"), "abc")
        self.assertEqual(rp.expect_sha("p7d.p5_summary_json", "default"), "default")
        vol = os.path.join(os.sep, "Volumes", "bulk")
        self.assertEqual(rp.volume_of(os.path.join(vol, "x", "y")), vol)
        with self.assertRaises(KeyError):
            rp.upstream("no_such_key")

    def test_min_toml_matches_tomllib(self):
        import tomllib
        for text in (CFG_TEXT, read(os.path.join(REPRO, "config.example.toml"))):
            self.assertEqual(rp._parse_toml_min(text), tomllib.loads(text))

    def test_stage_table(self):
        ids = [s["id"] for s in S.STAGES]
        self.assertEqual(len(ids), len(set(ids)))
        order = S.topo()
        self.assertEqual(sorted(order), sorted(ids))
        import run
        for st in S.STAGES:
            self.assertTrue(os.path.exists(os.path.join(PIPE, st["script"])), st["script"])
            for d in st["deps"]:
                self.assertIn(d, S.BY_ID)
                self.assertLess(order.index(d), order.index(st["id"]))
            self.assertIn(st["default"], st["modes"])
            for mode in st["modes"]:
                _st, _m, argvs, inputs, ensure = run.build(st["id"], mode)
                self.assertTrue(argvs and all(isinstance(v, list) for v in argvs), f"{st['id']}/{mode}")
                for a in [x for v in argvs for x in v] + list(inputs.values()) + ensure:
                    self.assertNotIn(rp.PLACEHOLDER, a, f"{st['id']}/{mode}")
                    self.assertNotIn("{", a, f"{st['id']}/{mode} unexpanded placeholder: {a}")

    def test_done_markers_written_by_script(self):
        """每个 done 标记都能在阶段脚本（及其 import 的 pipeline 模块）或该阶段 argv 里找到写出点（raw-p7e 曾因标记与实际输出不符出错）。"""
        def local_imports(fn, seen):
            if fn in seen or not os.path.exists(os.path.join(PIPE, fn)):
                return
            seen.add(fn)
            for node in ast.walk(ast.parse(read(os.path.join(PIPE, fn)))):
                names = [a.name for a in node.names] if isinstance(node, ast.Import) else \
                    [node.module] if isinstance(node, ast.ImportFrom) and node.module else []
                for n in names:
                    local_imports(n.split(".")[0] + ".py", seen)
        for st in S.STAGES:
            seen = set()
            local_imports(st["script"], seen)
            argvs = [a for v in st["modes"].values() for a in (v["steps"] if isinstance(v, dict) else [v])]
            argv = [x for a in argvs for x in a]
            text = "\n".join([read(os.path.join(PIPE, f)) for f in sorted(seen)] + [str(a) for a in argv])
            for d in st["done"]:
                *dirs, base = d.split("/")
                stem, ext = os.path.splitext(base)
                ok = base in text or (f'"{stem}"' in text and ext in text)
                self.assertTrue(ok and all(x in text for x in dirs), f"{st['id']}: done marker {d} not written by {sorted(seen)}")

    @unittest.skipUnless(os.path.isdir(os.path.join(REPRO, "..", "data", "products")), "data/products not present")
    def test_expected_tables_self_compare(self):
        """发布产物对自身比对无失败项；每个发布文件都按发布副本 sha 认出。"""
        import importlib.util
        spec = importlib.util.spec_from_file_location("diff_expected", os.path.join(REPRO, "tools", "diff_expected.py"))
        de = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(de)
        out = de.compare_files(de.PRODUCTS, False)
        self.assertEqual([o for o in out if o[3] in de.BAD], [])
        states = {o[3] for o in out}
        self.assertTrue(states <= {"same", "same(released-copy)", "same(rows-multiset)", "not-released", "not-run"}, states)
        self.assertIn("same(released-copy)", states)

    def test_stage_fields(self):
        """阶段表只含约定字段；每个阶段都有 ref_run 说明；共 39 个阶段，其中 33 个属于 --all。"""
        ids = [s["id"] for s in S.STAGES]
        allowed = {"id", "script", "title", "modes", "default", "deps", "done", "inputs", "stage_dir", "ensure_dirs", "net",
                   "token", "machine", "res", "ref_run", "in_all", "wind"}
        self.assertEqual(len(ids), 39)
        for st in S.STAGES:
            self.assertLessEqual(set(st), allowed, st["id"])
            self.assertIsInstance(st.get("ref_run"), str, st["id"])
            for d in st["deps"]:
                self.assertIn(d, ids)
        self.assertEqual(sum(1 for s in S.STAGES if s["in_all"]), 33)

    def test_raw_archive_constants_follow_config(self):
        env = dict(os.environ, REPRO_CONFIG=self.cfg)
        code = ("import raw_archive as ra, raw_store as rs; "
                "print(ra.P3_CACHE, ra.P1_CACHE, ra.P2_PIX, rs.RAW_ROOT_DEFAULT, ra.TOKEN_FILE == '/tmp/repro-t/token-not-read')")
        out = subprocess.run([sys.executable, "-c", code], cwd=PIPE, env=env, capture_output=True, text=True, check=True).stdout
        self.assertEqual(out.split(), ["/tmp/repro-t/fast/p3-cache", "/tmp/repro-t/fast/runs/T1/p1/cache",
                                       "/hist/p2/cmorph_pixels.csv", "/tmp/repro-t/raw", "True"])

    def test_raw_archive_never_overwrites(self):
        """raw 档只增不改：Copier.copy 遇同名不同大小抛错且原文件不变；p5 子集落档改为另存 .dup。"""
        env = dict(os.environ, REPRO_CONFIG=self.cfg)
        root, src = os.path.join(self.td.name, "raw"), os.path.join(self.td.name, "src")
        os.makedirs(src)
        code = f"""
import os, raw_archive as ra, repro_io as rio
class M:
    def add(self, **k): pass
class L:
    def log(self, m): pass
root, src = {root!r}, {src!r}
os.makedirs(os.path.join(root, ".tmp"), exist_ok=True)
open(os.path.join(src, "a.bin"), "wb").write(b"old")
cp = ra.Copier(root, M(), L(), "t")
print(cp.copy(os.path.join(src, "a.bin"), "d/a.bin", "x")[0])
print(cp.copy(os.path.join(src, "a.bin"), "d/a.bin", "x")[0])
open(os.path.join(src, "a.bin"), "wb").write(b"newer")
try:
    cp.copy(os.path.join(src, "a.bin"), "d/a.bin", "x")
    print("overwritten")
except FileExistsError:
    print("refused")
print(open(os.path.join(root, "d", "a.bin"), "rb").read().decode())
print(rio.dup_if_exists(root, "d/b.bin")[1])
d1 = rio.dup_if_exists(root, "d/a.bin")[1]
open(os.path.join(root, d1), "w").close()
d2 = rio.dup_if_exists(root, "d/a.bin")[1]
print(d1.startswith("d/a.bin.dup"), d2 == d1 + "_1")
"""
        out = subprocess.run([sys.executable, "-c", code], cwd=PIPE, env=env, capture_output=True, text=True, check=True).stdout
        self.assertEqual(out.split(), ["copied", "exists", "refused", "old", "d/b.bin", "True", "True"])
        txt = read(os.path.join(PIPE, "p5_sss_sat.py"))
        self.assertEqual(txt.count("dst, rel = _rio.dup_if_exists(raw_root, rel)"), 1)


def load_tool(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, os.path.join(REPRO, "tools", name + ".py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


DATA = os.path.normpath(os.path.join(REPRO, "..", "data"))


@unittest.skipUnless(os.path.isdir(os.path.join(DATA, "products")), "data/ not present")
class TestCompareTools(unittest.TestCase):
    """比对工具的规则：只差路径字段、非 in_all 阶段、译文文字格不算失败；数值或数字变了仍报 differs。"""

    def setUp(self):
        self.td = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.td.cleanup()

    def put_json(self, root, rel, obj):
        p = os.path.join(root, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w", encoding="utf-8") as f:
            json.dump(obj, f, ensure_ascii=False)

    def test_diff_expected_json_leaves_and_not_run(self):
        de = load_tool("diff_expected")
        root = self.td.name
        px = load_json(os.path.join(DATA, "products", "p4-explore", "p4_explore.json"))
        c12 = load_json(os.path.join(DATA, "products", "cal-surface-diag", "cal_c1_c2.json"))
        px["inputs_dir"] = "/somewhere/else"
        for v in c12["inputs"].values():
            v["path"] = "elsewhere/" + os.path.basename(v["path"])
        self.put_json(root, "p4-explore/p4_explore.json", px)
        self.put_json(root, "cal-surface-diag/cal_c1_c2.json", c12)
        st = {(o[1], o[3]) for o in de.compare_files(root, True)}
        self.assertEqual(st, {("p4-explore", "same(json-leaves)"), ("cal-surface-diag", "same(json-leaves)")})
        full = {o[1]: o[3] for o in de.compare_files(root, False)}
        self.assertEqual(full["dw33-bench"], "not-run")             # 不属于 --all，没跑不算缺
        self.assertEqual(full["p1"], "missing")                      # in_all 阶段缺文件仍是失败

        def bump(o):                                                  # 改第一个数值叶
            for k, v in o.items():
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    o[k] = v + 1
                    return True
                if isinstance(v, dict) and bump(v):
                    return True
            return False
        self.assertTrue(bump(c12))
        self.put_json(root, "cal-surface-diag/cal_c1_c2.json", c12)
        de = load_tool("diff_expected")                               # 新实例：load_json 按路径缓存
        st = {o[1]: o[3] for o in de.compare_files(root, True)}
        self.assertEqual(st["cal-surface-diag"], "differs")
        self.assertIn("differs", de.BAD)


class TestSwitches(unittest.TestCase):
    """风高开关与 module_sha：缺省逐位不变，非缺省按表生效。"""

    def setUp(self):
        self.td = tempfile.TemporaryDirectory()

    def tearDown(self):
        rp._CFG = None
        self.td.cleanup()

    def cfg(self, extra):
        p = os.path.join(self.td.name, "c.toml")
        with open(p, "w", encoding="utf-8") as f:
            f.write(CFG_TEXT.replace("[options]\n", "[options]\n" + extra))
        rp.load(p, force=True)

    def test_wind_default_is_identity(self):
        import math
        f4 = math.log(10.0 / 1e-4) / math.log(4.0 / 1e-4)
        rp._CFG = {}
        self.assertIs(rp.wind_factor("MOSEAN/WHOTS", f4, 1e-4), f4)
        self.cfg("")
        self.assertIs(rp.wind_factor("MOSEAN/WHOTS", f4, 1e-4), f4)
        self.cfg('wind_height = "uniform_4m"\n')
        self.assertIs(rp.wind_factor("SOFS", f4, 1e-4), f4)

    def test_wind_per_station(self):
        import math
        f4 = math.log(10.0 / 1e-4) / math.log(4.0 / 1e-4)
        self.cfg('wind_height = "per_station"\n')
        w = rp.wind_factor("MOSEAN/WHOTS", f4, 1e-4)
        self.assertAlmostEqual(w, math.log(1e5) / math.log(3.355 / 1e-4), places=14)
        self.assertTrue(1.015 < w / f4 < 1.02)                     # 约 +1.7%（WIND_H_BIAS 1.4–1.9% 之间）
        for st in ("KEO", "Papa", "TAO165E", "SYNB", "X"):
            self.assertIs(rp.wind_factor(st, f4, 1e-4), f4)
        with open(os.path.join(self.td.name, "c.toml"), "a", encoding="utf-8") as f:
            f.write("\n[wind_height_m]\nWHOTS = 3.5\n")
        rp.load(os.path.join(self.td.name, "c.toml"), force=True)
        self.assertAlmostEqual(rp.wind_factor("MOSEAN/WHOTS", f4, 1e-4), math.log(1e5) / math.log(3.5 / 1e-4), places=14)
        self.cfg('wind_height = "per_deployment"\n')
        with self.assertRaises(ValueError):
            rp.wind_factor("SOFS", f4, 1e-4)

    def test_wind_per_source(self):
        import math
        import run
        f4 = math.log(10.0 / 1e-4) / math.log(4.0 / 1e-4)
        self.cfg('wind_height = "per_source"\n')
        for st in ("MOSEAN/WHOTS", "SOFS", "KEO", "TAO165E"):
            self.assertIs(rp.wind_factor(st, f4, 1e-4), f4)       # 逐源高度由 dw13a 包装层施加，通用换算保持 4 m
        self.assertIsNotNone(run.wind_conflict(S.BY_ID["p2"], "per_source"))
        self.assertIsNone(run.wind_conflict(S.BY_ID["dw13a"], "per_source"))
        self.assertIsNotNone(run.wind_conflict(S.BY_ID["dw13a"], "per_station"))
        self.assertIsNone(run.wind_conflict(S.BY_ID["p2"], "uniform_4m"))
        _st, mode, argvs, _i, _e = run.build("dw13a")
        self.assertEqual(mode, "full")
        self.assertEqual([v[2] for v in argvs], ["--selftest", "--stage", "--gate", "--stage", "--compare"])

    def test_dw13a_patch_sites(self):
        txt = read(os.path.join(PIPE, "dw13a_windheight.py"))
        for s in ('P2_SHA = _rp.module_sha("p2_rim_test.py")', 'META_SHA = _rp.expect_sha("dw13a.dw13_meta"',
                  'META_DEFAULT = _rp.upstream("dw13_meta")', 'ORIG_DEFAULT = _rp.upstream("p2_dir")'):
            self.assertEqual(txt.count(s), 1, s)

    def test_wind_patch_sites(self):
        for fn, n in (("p2_rim_test.py", 1), ("p4b_posthoc.py", 1), ("p5_sss_sat.py", 1)):
            txt = read(os.path.join(PIPE, fn))
            self.assertEqual(txt.count("_rp.wind_factor(stn[\"name\"]"), n, fn)

    def test_module_sha(self):
        """module_sha 取 pipeline/ 当前文件的 sha256；dw33 与 dw13a 都用它核 import 的模块，与静态场 meta 的写法一致。"""
        for rel in ("p3_global.py", "p2_rim_test.py"):
            with open(os.path.join(PIPE, rel), "rb") as f:
                self.assertEqual(rp.module_sha(rel), hashlib.sha256(f.read()).hexdigest())
        with self.assertRaises(FileNotFoundError):
            rp.module_sha("no_such_module.py")
        self.assertEqual(read(os.path.join(PIPE, "dw33_wind_bins.py")).count('P3_CODE_SHA = _rp.module_sha("p3_global.py")'), 1)
        self.assertEqual(read(os.path.join(PIPE, "p3_global.py")).count("code_sha = sha256_file(os.path.abspath(__file__))"), 2)
        self.assertFalse(os.path.exists(os.path.join(PIPE, "ASSEMBLY.json")))
        for p in glob.glob(os.path.join(PIPE, "*.py")) + [os.path.join(REPRO, "run.py")]:
            self.assertNotIn("assembled_sha", read(p), p)
            self.assertNotIn("ASSEMBLY.json", read(p), p)



class TestUpstreamShaGate(unittest.TestCase):
    """输入 sha 门：[options] check_upstream_sha 缺省（或 false）时不停、只记录；为 true 时不符即停（与原行为相同）。"""
    SITES = {"p7d_p5_timebin_ci.py": 2, "p8a_sat_wind_strata.py": 1, "dw13a_windheight.py": 1, "cal_surface_diag.py": 1,
             "dw33_post.py": 1}

    def setUp(self):
        self.td = tempfile.TemporaryDirectory()

    def tearDown(self):
        rp._CFG = None
        self.td.cleanup()

    def cfg(self, extra, sha_line=""):
        p = os.path.join(self.td.name, "cfg.toml")
        with open(p, "w", encoding="utf-8") as f:
            f.write(CFG_TEXT.replace("[options]\n", "[options]\n" + extra) + sha_line)
        rp.load(p, force=True)

    def test_switch(self):
        for extra, stop in (("", False), ("check_upstream_sha = false\n", False), ("check_upstream_sha = true\n", True)):
            self.cfg(extra)
            self.assertEqual(rp.upstream_sha_mismatch("abcdef", "abd"), stop, extra)
            self.assertFalse(rp.upstream_sha_mismatch("abcdef", "abc"))
            self.assertFalse(rp.upstream_sha_mismatch("abcdef", None))

    def test_gate_sites_use_switch(self):
        for fn, n in self.SITES.items():
            self.assertEqual(read(os.path.join(PIPE, fn)).count("_rp.upstream_sha_mismatch("), n, fn)

    def test_dw13a_meta_gate(self):
        """dw13a 读 D13M 的门：与参考运行 sha 不同的 dw13_meta.json 在缺省下照常读入；打开开关时停下，
        [upstream_sha] 填入该文件的 sha 后放行（与原行为相同）。"""
        import importlib
        import dw13a_windheight as dw
        other = os.path.join(self.td.name, "dw13_meta.json")
        with open(other, "w", encoding="utf-8") as f:
            json.dump({"wind": [{"file": "/x/a.nc", "station": "S", "heights_m": {"z": 3.5}}]}, f)
        self.cfg("")
        dw = importlib.reload(dw)
        self.assertEqual(dw.load_meta(other)["station_median"], {"S": 3.5})
        self.cfg("check_upstream_sha = true\n")
        dw = importlib.reload(dw)
        with self.assertRaises(RuntimeError):
            dw.load_meta(other)
        with open(other, "rb") as f:
            h = hashlib.sha256(f.read()).hexdigest()
        self.cfg("check_upstream_sha = true\n", f'"dw13a.dw13_meta" = "{h}"\n')
        dw = importlib.reload(dw)
        self.assertEqual(dw.load_meta(other)["sha256"], h)

if __name__ == "__main__":
    unittest.main()
