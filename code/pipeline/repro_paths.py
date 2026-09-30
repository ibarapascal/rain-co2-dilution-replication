"""repro_paths.py — 集中路径配置的读取层（只解析路径，不含任何分析逻辑，不读凭据内容）。

各阶段脚本的路径默认值一律通过本模块取得：
  path(key)            [paths] 表里的路径（raw_root、fast_root、out_root、token_file、p3_cache …），支持 ~、$VAR、{其他键}
  upstream(key)        上游阶段产物位置：[upstream] 表显式给出的优先，否则按 out_root/<阶段>/… 推出（见 UPSTREAM_DEFAULTS）
  raw_sub(name)        raw_root 下的子目录（如 spurs2）
  expect_sha(key, v)   脚本里写死的「上游输入 sha256」闸门：[upstream_sha] 表可换成新值（重跑后由操作者显式填写）
  upstream_sha_mismatch(got, want)  输入 sha 门是否应当停下：只有 [options] check_upstream_sha = true 时才核对
                       （缺省 false＝各门只记录实际 sha，不因上游 sha 不同而 exit 3）；want 为空时不核对；按前缀比
  get(key, default)    [options] 表的标量（如 raw_rate_mbps）
  volume_of(p)         <mount-root>/<卷> 前缀（挂载哨兵用）
  path_or(key, d)      [paths] 表里有该键就用它，否则返回 d（可选路径）
  module_sha(rel)      pipeline/ 里 rel 文件当前内容的 sha256；给「核 import 模块自身 sha」的脚本用
                       （dw13a 核 p2_rim_test.py，dw33 核 p3_global.py）
  wind_factor(name, f4, z0)   B 群风高开关（[options] wind_height）：缺省 "uniform_4m" 原样返回 f4（＝p2.WIND_FACTOR，逐位不变）；
                       "per_station" 按站实测风速计高度 z 返回 ln(10/z0)/ln(z/z0)（表见 WIND_HEIGHT_PER_STATION，可由
                       [wind_height_m] 覆盖）；表里没有的站（A 群 TAO／BOBOA、合成自测站）仍返回 f4；
                       "per_source"＝逐源高度（B 群逐小时按提供该小时风值的（文件, 变量）取高度），只由 dw13a
                       阶段的 dw13a_windheight.py 包装层实现（它在站点小时风上先乘 r(z)＝f(z)/f(4)），故本函数在此模式下返回 f4、
                       通用换算保持 4 m；其余直接受影响的阶段没有逐源追踪，run.py 在 per_source 下拒绝运行它们

配置文件路径由环境变量 REPRO_CONFIG 给出（../run.py 负责设置）。未设置时 path()/upstream() 返回
`/__REPRO_UNCONFIGURED__/<键>` 占位串：模块照常可 import（合成自测不受影响），但任何真实读写都会因路径不存在而立即失败，
不会静默落到某个默认盘上。token_file 只返回路径字符串，本模块从不打开它。

Change Log：
  2026-09-27 初版。
  2026-09-29 加 path_or、wind_factor（风高开关，缺省＝4 m）；upstream 键补 P6–P8、dw33。
  2026-09-29c 风高开关加 "per_source"（dw13a 阶段实现）；upstream 键加 dw13_meta、dw13a_dir。
  2026-10-01 module_sha（核 import 模块自身 sha，取 pipeline/ 当前文件）。
  2026-10-01b upstream_sha_mismatch：输入 sha 门改为默认关闭（[options] check_upstream_sha）。
"""

import os
import re

ENV = "REPRO_CONFIG"
PLACEHOLDER = "/__REPRO_UNCONFIGURED__"

# upstream 键 → (相对 out_root 的默认位置)；{p1_events} 这类引用在展开时按同表递归解析
UPSTREAM_DEFAULTS = {
    "p1_dir": "{out_root}/p1",
    "p1_events": "{p1_dir}/p1_events.csv",
    "p1_cache": "{p1_dir}/cache",
    "p1b_dir": "{out_root}/p1b",
    "p1b_cache": "{p1b_dir}/cache",
    "p2_dir": "{out_root}/p2",
    "p2_pixels": "{p2_dir}/cmorph_pixels.csv",
    "p2sens_dir": "{out_root}/p2-sens",
    "p2sens_pixels": "{p2sens_dir}/cmorph_pixels.csv",
    "p3_dir": "{out_root}/p3-global",
    "p3b_curves": "{p3_dir}/p3b_curves.json",
    "p3c": "{out_root}/p3c/p3c_decision.json",
    "p4_dir": "{out_root}/p4-mech",
    "p4_events": "{p4_dir}/p4_mech_events.csv",
    "p4b_dir": "{out_root}/p4b",
    "p5_dir": "{out_root}/p5-sss",
    "p6_dir": "{out_root}/p6",
    "p7a_dir": "{out_root}/p7a",
    "p7c_dir": "{out_root}/p7c",
    "p7c_events": "{p7c_dir}/p7c_events.csv",
    "p7e_dir": "{out_root}/p7e",
    "p8a_dir": "{out_root}/p8a",
    "p8b_dir": "{out_root}/p8b",
    "dw33_dir": "{out_root}/dw33",
    "dw13_meta": "{out_root}/cal-dw13-meta/dw13_meta.json",
    "dw13a_dir": "{out_root}/dw13a",
}

# 风高开关：B 群逐站主风速计高度（m）＝D13M `dw13_meta.json`（sha256 fd5a036b…，cal_dw13_meta.py）
# 里该站逐部署文件（D_M／D_MET／ASIMET；不含跨多年的合并文件）主传感器高度的中位数。KEO、Papa 的中位数就是 4.0，与缺省 4 m 相同；
# Stratus 无合格事件，列出只为完整。A 群（TAO、BOBOA）未核，按 4 m（返回原系数）。
WIND_HEIGHT_PER_STATION = {"MOSEAN/WHOTS": 3.355, "SOFS": 3.3125, "Stratus": 3.32, "KEO": 4.0, "Papa": 4.0}
WIND_HEIGHT_ALIASES = {"WHOTS": "MOSEAN/WHOTS"}          # [wind_height_m] 表的键不能含「/」（最小 TOML 解析）
WIND_HEIGHT_MODES = ("uniform_4m", "per_station", "per_source")

_CFG = None


def _parse_toml_min(text):
    """tomllib 不可用时（Python <3.11）的最小解析：只认 [表]、key = "字符串" / 数字 / true|false、# 注释。"""
    out = {}
    cur = out
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip() if '"' not in raw else raw.strip()
        if not line or line.startswith("#"):
            continue
        m = re.fullmatch(r"\[([A-Za-z0-9_.\-]+)\]", line)
        if m:
            cur = out.setdefault(m.group(1), {})
            continue
        m = re.fullmatch(r'"?([A-Za-z0-9_.\-]+)"?\s*=\s*"([^"]*)"\s*(#.*)?', line)
        if m:
            cur[m.group(1)] = m.group(2)
            continue
        m = re.fullmatch(r'"?([A-Za-z0-9_.\-]+)"?\s*=\s*([^#\s]+)\s*(#.*)?', line)
        if m:
            v = m.group(2)
            if v in ("true", "false"):
                cur[m.group(1)] = v == "true"
            else:
                try:
                    cur[m.group(1)] = int(v)
                except ValueError:
                    cur[m.group(1)] = float(v)
            continue
        raise ValueError(f"配置行无法解析：{raw!r}")
    return out


def load(path=None, force=False):
    """读配置（缓存）。path 缺省取环境变量 REPRO_CONFIG；都没有则返回空配置（占位模式）。"""
    global _CFG
    if _CFG is not None and not force and path is None:
        return _CFG
    p = path or os.environ.get(ENV)
    if not p:
        _CFG = {}
        return _CFG
    with open(p, "rb") as f:
        data = f.read()
    try:
        import tomllib
        cfg = tomllib.loads(data.decode("utf-8"))
    except ImportError:
        cfg = _parse_toml_min(data.decode("utf-8"))
    cfg["_file"] = os.path.abspath(p)
    _CFG = cfg
    return cfg


def configured():
    return bool(load())


def _expand(s, table, seen):
    def rep(m):
        k = m.group(1)
        if k in seen:
            raise ValueError(f"配置循环引用：{k}")
        return _resolve(k, table, seen | {k})
    s = re.sub(r"\{([A-Za-z0-9_]+)\}", rep, s)
    return os.path.expanduser(os.path.expandvars(s))


def _resolve(key, table, seen=frozenset()):
    cfg = load()
    paths = cfg.get("paths", {})
    ups = cfg.get("upstream", {})
    if key in ups and ups[key]:
        return _expand(str(ups[key]), table, seen)
    if key in paths and paths[key]:
        return _expand(str(paths[key]), table, seen)
    if key in UPSTREAM_DEFAULTS:
        return _expand(UPSTREAM_DEFAULTS[key], table, seen)
    return f"{PLACEHOLDER}/{key}"


def path(key):
    """[paths] 表中的路径；未配置返回占位串。"""
    if not configured():
        return f"{PLACEHOLDER}/{key}"
    return _resolve(key, "paths")


def upstream(key):
    """上游产物位置：[upstream] 显式值 > out_root 推出的默认值；未配置返回占位串。"""
    if key not in UPSTREAM_DEFAULTS and key not in load().get("upstream", {}):
        raise KeyError(f"未知 upstream 键：{key}")
    if not configured():
        return f"{PLACEHOLDER}/{key}"
    return _resolve(key, "upstream")


def raw_sub(name):
    return os.path.join(path("raw_root"), name)


def expect_sha(key, frozen):
    """脚本内写死的输入 sha256 闸门：默认沿用参考运行的产物 sha；[upstream_sha] 里给了该键就用它（重跑后由操作者显式填写并记录）。"""
    v = load().get("upstream_sha", {}).get(key)
    return v if v else frozen


def upstream_sha_mismatch(got, want):
    """输入 sha 门：[options] check_upstream_sha 为 true 且 want 非空、got 不以 want 开头时返回 True（调用方据此停下）。"""
    if get("check_upstream_sha", False) is not True or not want:
        return False
    return not got.startswith(want)


def get(key, default=None):
    return load().get("options", {}).get(key, default)


def volume_of(p):
    """<mount-root>/<卷>/… → <mount-root>/<卷>；其他路径返回其根（/）。"""
    parts = os.path.abspath(p).split(os.sep)
    if len(parts) > 2 and parts[1] == "Volumes":
        return os.sep + os.path.join(parts[1], parts[2])
    return os.sep


def path_or(key, default):
    """[paths] 里显式给出的路径；没给（或未配置）时返回 default（不产生占位串）。"""
    v = load().get("paths", {}).get(key)
    return _resolve(key, "paths") if v else default


def module_sha(rel):
    """pipeline/ 里 rel 文件当前内容的 sha256（导入它的脚本据此确认 import 的正是同目录这份文件）。"""
    import hashlib
    here = os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha256()
    with open(os.path.join(here, rel), "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def wind_height_mode():
    m = get("wind_height", "uniform_4m")
    if m not in WIND_HEIGHT_MODES:
        raise ValueError(f"[options] wind_height 只能是 {WIND_HEIGHT_MODES}，得到 {m!r}")
    return m


def wind_heights():
    """per_station 模式下实际使用的逐站高度表（内置表＋[wind_height_m] 覆盖）。"""
    tab = dict(WIND_HEIGHT_PER_STATION)
    for k, v in load().get("wind_height_m", {}).items():
        tab[WIND_HEIGHT_ALIASES.get(k, k)] = float(v)
    return tab


def wind_factor(station, factor_4m, z0):
    """站点风速换算到 10 m 的系数。uniform_4m（缺省）原样返回 factor_4m（同一个 float，逐位不变）；per_source 同样原样返回
    （逐源高度由 dw13a 包装层在小时风上先乘 r(z)，此处再乘 4 m 系数；见模块头）。"""
    if wind_height_mode() in ("uniform_4m", "per_source"):
        return factor_4m
    z = wind_heights().get(station)
    if z is None or z == 4.0:
        return factor_4m
    import math
    return math.log(10.0 / z0) / math.log(z / z0)


def run_tag(stage):
    """写进 raw/manifest.jsonl 的 source_job：<阶段>:<$REPRO_ATTEMPT_ID 或 local>。"""
    return f"repro-{stage}:{os.environ.get('REPRO_ATTEMPT_ID', 'local')}"
