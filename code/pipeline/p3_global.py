#!/usr/bin/env python3
"""P3 全球重算：Witte 等 2026 胶囊逐式移植＋稀疏 RIM-3＋s 网格，P3a（复现）与 P3b（盲算曲线）同一遍。

设计：X1–X10 与 L1–L24 在运行前写定；L25 起为本实现的选择（见下；L26 留给环境决定，不在本文件）。
公式与常数以 Code Ocean 胶囊为准（code/CO2_Rain_Flux_Toolbox.py sha256 cbce94db…、code/main.py sha256 53017ba3…；
注释里的「TB Lnn」「main Lnn」是这两个文件的行号）；原文描述与代码不一致时用代码（L11）。
**本脚本不读任何 P2 产物（L23）**：只 import 同目录 p2_rim_test 的 rim_factor 作等价性自测的参考实现（L4）。P3c（读 P2、
判门定档）是另一步骤。

用法（产物写 $REPRO_OUTPUT_DIR，未设则 --out）：
  p3_global.py --mode selftest [--cache-root DIR]     代码等价性自测（也是门）→ selftest.json；不过则退出 4。
                                                      合成部分无需数据；--cache-root 下 200001 已就绪时加跑真实块。
  p3_global.py --mode bench --cache-root DIR          自测（真实块必跑）＋2000-01-01 一天（含 1999-12-31 起转）逐步计时、
                                                      峰值 RSS、L25 对精确 gsw 逐值比较、整年墙钟外推 → bench.json
  p3_global.py --mode full --cache-root DIR --workers 4
                                                      自测 → 按月并行（每月前补 24 h 起转）→ p3a_repro.json（先写）、
                                                      p3b_curves.json、p3_bands.csv、p3_monthly.csv、p3_run.json
  （内部）--mode capsule-ref --ref-in X --ref-out Y --capsule-dir D：import 胶囊原函数（xarray）算参考值，供自测比对；
         解释器可用 --capsule-python 另指（L26 若把 xarray 放进单独 venv）。
其他参数：--work-dir D（静态场与月结果；默认 <out>/_work，跨次续跑须给固定路径）--months 1,2（只跑这些月，调试用；
  full 缺月不出 p3a）--chunk N（有雨像元分块，控内存）--wait-hours H（等下载，
  默认 30）--capsule-dir D（胶囊两文件所在目录；缺省依次找 src/witte_capsule、<cache>/witte_capsule，都没有则按
  Code Ocean 公开 API 取到 <out>/_capsule 并核 sha256）。

数据布局（--cache-root，下载侧 p3_fetch.py 约定；所有路径与变量名假设集中在 Layout 类）：
  cmorph/YYYY/MM/CMORPH_V1.0_ADJ_8km-30min_YYYYMMDDHH.nc   era5_ws10/YYYYMM.nc（ws10 int16×0.01，60N→60S，0–359.75）
  hycom_sss/YYYYMMDD.nc（salinity）  oisst/oisst-avhrr-v02r01.YYYYMMDD.nc  watson/*.nc（tar 已解）  glodap/GLODAPv2.2016b.*.nc
  month_ready/YYYYMM.json（该月齐；200001 含起转日与静态源）  可选 woa09/*basin*.csv（WOA09 分区；缺则只按纬带）
依赖：numpy、scipy、netCDF4、gsw、PyCO2SYS（自测另需 xarray＋pandas，仅 capsule-ref 子进程用）。
内存：主进程与每个 worker 只常驻「当日/当月」场与 48 槽稀疏缓冲，逐时读 CMORPH/ERA5，不整年入内存（目标 ≈2 GB/worker）。
退出码：0 完成；2 数据缺失/超时；3 其他异常（含 L8/L46：直接 fCO2 与公开 pCO2 都缺，或 pCO2 单位不是 µatm）；4 等价性自测不过（不做后续计算）。

实现选择（设计未写死处；阈值、判门、估计量一律照设计）：
  L25 s≠1 的溶解度按 ΔS 二阶展开（中心 S0，每像元每日一次精确值与一、二阶导；S0/SST 逐日常数，故与逐步重算逐值相同）；
      s=1 复现路径逐值精确 gsw；基准时逐值比较相对误差 <1e-6。
  L27 展开在 ln α 上做：Weiss 部分 ln 对 S 严格线性（精确），只有 ln ρ 二阶截断；ρ 导数用 gsw.rho_first/second_derivatives
      （照胶囊把 SST 当 CT 传入 gsw.rho）；SA 对 SP 逐像元仿射，斜率 r 用 gsw 取；|ΔS|>3 psu 或 S<15 psu 的元素改用精确 gsw。
  L28 整体缩放 s>1 使 S<0 时截到 0（照胶囊对 S20 的 where(sal>0, other=0)），逐情形计数。
  L29 稀疏 RIM：条目＝本步「非 (P==0 且 0<U≤200)」的像元（有雨、NaN、P>200、U==0），因子逐式照胶囊算术（0/0→NaN、
      U=0 有雨→0）；对数域累加；某像元 24 h 内任一条目 NaN → RIM 全族（整体/历史/z5，含 s=0）同记缺测。
  L30 胶囊稀释情形大气侧溶解度用 S_dil 而不加 0.1（TB L498 对 L501/L504），此不对称对**所有**海洋像元（含无雨）成立，
      按代码保留：非 D（24 h 无雨）像元贡献一个与 s 无关的常数项「skin」，单独累计并报出；非 D 像元 kH≡kw、Fdep=0。
  L31 空间插值：ERA5/OISST/HYCOM 双线性、NaN 传播（与 xarray interp 默认一致）、经度周期；OISST 冰插值后 fillna(0)，
      值域 ≤1 视作分数再 ×100（main L73）；SST、S0 在 UTC 日内常数（与 L9 一致）。
  L32 风时间（L6）：HH:30＝前后整点插值场的平均；2000-12-31 23:30 缺 2001-01-01 00Z（不在下载范围）时沿用 23:00 并计数；
      每月起转日用前月文件。
  L33 Watson：变量 sfco2、fco2atm_skinT（L8 替补名单见 WATSON_AIR_FALLBACK；缺则走 L46）；时间按 units 解码取 2000 年；
      Watson 与 GLODAP 斜率均照 main L52–60 做「经度 %360、排序、补 0/360、最近邻」；GLODAP 斜率向量化调用 PyCO2SYS，
      斜率存 float32（main L39）；海陆掩码＝斜率 NaN（L10）。
  L34 累计：面积＝胶囊 grid_cell_areas（TB L628）于 CMORPH 坐标；逐步×面积求和到 gC，全年 17,568 步；百分比分母＝全体有效
      Fwind（非配对），各效应按 L12 与同像元同时次 Fwind 配对相减；另报配对分母与「Witte 式非配对」两种诊断；分区＝5° 纬带
      （×WOA09 分区，若有掩码文件）。
  L35 CMORPH：文件缺或时间轴不符 → 该半步整场 NaN（L12）；S20 在 P 为 NaN 时记缺测（胶囊会按 NaN>0 为假置 0）。
  L36 全程 float64（胶囊在 float32 文件输入上会以 float32 算 IRR、Kz、SSTk 等中间量）；等价性自测以 float64 调胶囊，float32 只作信息。
  L37 自测判据：通量逐像元 |a−b|≤1e-6·max(|b|, 1e-9 gC m⁻² 每半步)，盐度同式（floor 1e-9 psu），NaN 模式须一致；RIM 与 P2 逐值
      相对差 ≤1e-12；合成块＋真实块（1999-12-31 起转、2000-01-01 00Z 小时，按起转日雨像元数选 50×50 块）；另核 GLODAP 斜率抽样与面积。
  L38 z=5 m（只复现 62%，非门）用 L25 展开＋兜底；S20（门：RIM/S20 比）与 s=1 用精确 gsw。
  L39 每情形另累计稀释通量中对 β 线性的分量 G＝Tr·β·ΔS·αw·(1−ice)，P3c 可对任意 r_β 精确得 S6（X10）；S20 只算 s=1。
  L40 按月并行（spawn，ProcessPoolExecutor，worker 被杀不挂起），每月前补前一日 48 个半步起转；worker 等 month_ready(m)、
      month_ready(m−1)、era5_ws10/(m+1)；月结果写 <work-dir>/month_YYYYMM.*，参数哈希与天数一致可续跑（跨次续跑须 --work-dir
      固定路径）；静态场按 VERSION＋代码 sha256＋cache 根判沿用；12 个月齐才汇总。
  L41 bench 外推：整年墙钟＝Σ月(天数×(48×单步＋日开销)＋48×起转步) 按 LPT 排到 W 个 worker，再加静态场与自测；
      单步耗时扣除 bench 专有的 L25 抽检时间（full 不抽检）。
  L46 公开 RECCAP2 包只有 pCO2（spco2、pco2atm，µatm），无 sfco2/fco2atm_skinT：
      优先路径不变（sfco2＋大气 fCO2 齐→直接用，"direct"）；否则两侧都用 pCO2×FugFac(T)（Weiss 1974＝CO2SYS 同式，
      P＝1.01325 bar，R＝83.14462618），"pco2_fugacity_L46"；pCO2 也缺或单位非 µatm→退出 3。顺序：月、1° pCO2 照 L33 最近邻
      落到 CMORPH 像元（不改），再逐日按像元 SST 换算——海侧 T＝海侧溶解度所用 SST，大气侧 T＝SST−SKIN_DT（与溶解度同一常数）。
      p3a diagnostics 记 watson_fco2_source 与全场 FugFac min/mean/max；自测新增 FugFac 对 PyCO2SYS gas.fugacity_factor ≤1e-9。
  L47 RIM 历史含 NaN 型条目（P/U 为 NaN、P>200、P=U=0、整场缺测等，即
      RimRing.history 的 logH 为 NaN）→ 此后 24 h RIM 全族缺测，L3/L12/L29 原样，引擎不动；胶囊 np.prod 作用于 DataArray 实为
      xarray nanprod（NaN 槽按 1 跳过），此缺测副作用不采用（同 L35）。胶囊逐像元自测对这类像元×半步只不比 RIM 三场（salR、
      Fdil_RIMv3、Fint_RIMv3），并正向断言本实现三场全 NaN、胶囊 f64 三场全有限；check_accum 的 dil_RIM、int_RIM 同掩码扣除；
      selftest 各块报 rim_excluded_nan_history。p3a coverage 另记 u10_nan_frac、precip_u10_both_zero_count（只加计数）。

Change Log:
  2026-09-26 a：初版（未在任何数据上运行过，只做 py_compile）。
  2026-09-26 a（独立审查修订，VERSION 不变）：--work-dir 与续跑判据（天数、静态场代码 sha）；mp.Pool→ProcessPoolExecutor；
      自测 KEEP 的无雨 salS20 照胶囊 where(>0, 0)（S0 为 NaN 时 0）；bench 单步扣 L25 抽检；WOA09 只读 *basin*.csv；
      main 兜底异常退出 3；bands 符号说明。
  2026-09-27 a（VERSION p3-2026-09-27a）：L46——Watson 公开包只有 pCO2，两侧改为 pCO2×FugFac(T) 逐日换算（直接 fCO2 仍为优先
      路径）；p3a/bench 记 watson_fco2_source 与 FugFac 统计；自测新增 fugacity_pyco2sys 项（计入自测门）；真实块在 L46 下
      引擎喂 pCO2、胶囊喂测试侧换算的 fCO2。
  2026-09-27 b（VERSION p3-2026-09-27b）：L47——引擎数值路径不变；胶囊自测对「24 h 历史含 NaN 型条目」的像元×半步不比 RIM
      三场（掩码＝RimRing.history 的 np.isnan(logH)，与 Engine.step 同一调用）＋正向断言＋累计同掩码扣除，selftest 报
      rim_excluded_nan_history；BASE_Q 加 n_Unan、n_PU0，p3a coverage 加 u10_nan_frac、precip_u10_both_zero_count。
"""
import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import argparse
import concurrent.futures as cf
import glob
import hashlib
import json
import math
import multiprocessing as mp
import re
import resource
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import traceback
import urllib.request
from datetime import date, datetime, timedelta, timezone

import numpy as np

VERSION = "p3-2026-09-27b"
SRC_DIR = os.path.dirname(os.path.abspath(__file__))
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
DEFAULT_CACHE = _rp.path("p3_cache")  # [repro] 路径来自集中配置（本常量在主程序里未被使用）
YEAR = 2000
UTC = timezone.utc

# ======================================================================== 胶囊常数（TB＝CO2_Rain_Flux_Toolbox.py）
C1 = 5.7                                                   # TB L37
C2 = 1 / 11                                                # TB L38
TC = np.arange(600, 1800, 25).astype(float)               # TB L39（48 个）
TI = (-np.arange(0, 24 * 3600, 1800) + 24 * 3600).astype(float)   # TB L40
KZ_COEF = 2.5 * (10 ** (-5))                               # TB L65、L70
U_GRID = np.array([0, 2, 4, 6, 8, 10, 200], float)        # TB L51
R_GRID = np.array([0, 2, 5, 10, 20, 50, 200], float)
D0_TABLE = np.array([[1.2, 1.2, 1.0, 1.0, 1.0, 1.1, 1.1],  # TB L43–L50（行＝U，列＝R）
                     [1.2, 1.2, 1.0, 1.0, 1.0, 1.1, 1.1],
                     [1.6, 1.6, 1.9, 2.0, 2.0, 2.2, 2.2],
                     [1.5, 1.5, 2.5, 2.9, 3.1, 3.5, 3.5],
                     [1.9, 1.9, 2.7, 3.6, 4.2, 4.7, 4.7],
                     [2.4, 2.4, 3.0, 4.0, 5.0, 5.8, 5.8],
                     [2.4, 2.4, 3.0, 4.0, 5.0, 5.8, 5.8]])
S20_A, S20_B = -0.35, 0.77                                 # TB L97–L98
SC_A, SC_B, SC_C, SC_D = 2073.1, 125.62, 3.6276, 0.043219  # TB L120–L123
NSA1, NSA2, NSA3 = -60.2409, 93.4517, 23.3585              # TB L144–L146
NSB1, NSB2, NSB3 = 0.023517, -0.023656, 0.0047036          # TB L147–L149
GAMMA_K = 0.266                                            # TB L181
TCONV_K = 10 ** -2 * 0.5                                   # TB L212（'30min'）
TCONV_R = 10 ** -3 * 0.5                                   # TB L396（'30min'）
SKIN_DT, SKIN_DS = 0.17, 0.1                               # TB L498
DIC_RAIN, K_CARBONIC = 25, 17                              # TB L414
EARTH_RADIUS = 6371000.0                                   # TB L554
Z5 = 5.0                                                   # 原文「z=5」（L2）
T_CURRENT_S = 1.0                                          # K3/L2：当前项深度因子的 t
AGE_TC = TC[48 - np.arange(1, 49)]                         # 历史年龄 a=1..48 ↔ 胶囊窗口下标 48−a
AGE_TI = TI[48 - np.arange(1, 49)]
AGE_W = AGE_TC / np.sqrt(AGE_TI)

# ======================================================================== 事先写定的设计常数（不得改）
S_GRID = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0)   # 第三节 P3b、L15
CAPTURE_C = 0.30                                           # L13
P3_0_LINE = (3.978, 6.671)                                 # 1.4 节 U(s)≈3.978+6.671·s
BAND_DEG, LAT_LO, LAT_HI = 5.0, -60.0, 60.0                # 8.1 节 5° 纬带
NBANDS = int((LAT_HI - LAT_LO) / BAND_DEG)
REPRO_GATES = [  # 第四节复现门：(键, 目标, 类型, 容差)
    ("wind_PgC", -1.468, "rel", 0.05),
    ("dilution_RIM_pct", 7.97, "rel", 0.10),
    ("interfacial_RIM_pct", 9.67, "rel", 0.10),
    ("deposition_pct", 3.88, "rel", 0.05),
    ("dilution_ratio_RIM_S20", 2.05, "rel", 0.10),
    ("turbulence_pct", 0.14, "abs", 0.10),
    ("U1_pp", 10.65, "abs", 1.0),
]
WITTE_Z5_CAPTURE = 0.62

# ======================================================================== 实现层参数（L 条）
L25_DMAX, L25_SMIN, L25_TOL = 3.0, 15.0, 1e-6              # L27 兜底阈值；L25 容差
TOL_CAPSULE, FLOOR_FLUX, FLOOR_SAL, TOL_P2 = 1e-6, 1e-9, 1e-9, 1e-12   # L37
CHUNK = 262144
POLL_S = 300
WATSON_SEA = ("sfco2",)
WATSON_AIR = "fco2atm_skinT"
WATSON_AIR_FALLBACK = ("fco2atm_skin", "fco2atm", "fco2_atm", "fco2atm_subskinT")   # L8/L33（只接受 fCO2 类大气字段）
WATSON_SEA_PCO2, WATSON_AIR_PCO2 = "spco2", "pco2atm"      # L46：公开包的海侧／大气 pCO2
PCO2_UNITS_OK = ("uatm", "µatm", "μatm", "microatm", "microatmosphere", "microatmospheres", "micro-atmospheres")
FUG_P_BAR, FUG_R, FUG_T0 = 1.01325, 83.14462618, 273.15    # L46：Weiss 1974／CO2SYS（R＝CODATA 2018，cm³·bar/(mol·K)）
FUG_TOL, FUG_TEST_T = 1e-9, (-2.0, 0.0, 10.0, 20.0, 30.0)  # L46 自测
SRC_DIRECT, SRC_L46 = "direct", "pco2_fugacity_L46"
GLODAP_VARS = ("TAlk", "TCO2", "salinity", "temperature")
CAPSULE_FILES = {"CO2_Rain_Flux_Toolbox.py": "cbce94db03fc4a3bc988146ceea71509bd3c85980e7c6b1827126eb2d863df7b",
                 "main.py": "53017ba313bc094c423ece48e22e3a0ef94ac0db912c2947e1f8c69d804002b1"}
CAPSULE_API = ("https://codeocean.com/api/capsules/6b8891ef-0ca8-43d0-ba62-edfdad8e65d5/blob"
               "?owner_id=verified&path=code/{name}&commit=HEAD&version=1")
UA = "rain-co2-dilution-replication-p3/1.0 (research; python)"

CASES = [("whole", s) for s in S_GRID] + [("hist", s) for s in S_GRID] + [("s20", 1.0), ("z5", 1.0)]   # 24 情形
LABELS = [f"{k}_{s:g}" for k, s in CASES]
BASE_Q = ["n_ocean", "n_W", "W_all", "skin", "u_nonD_dil", "u_nonD_W", "n_D", "n_C", "n_Pnan", "n_Pgt200", "n_U0",
          "n_Unan", "n_PU0", "n_poisonR", "dep", "Wm_dep", "nm_dep", "u_dep", "turb", "Wm_turb", "nm_turb", "u_turb_D"]
CASE_Q = ["dil", "int", "gd", "gi", "Wm_dil", "nm_dil", "Wm_int", "nm_int", "nfb", "nclip"]
UNPAIRED_CASES = ("whole_1", "s20_1")
KEEP_FIELDS = ("Fwind", "Fdep", "Fturb", "Fdil_RIMv3", "Fdil_S20", "Fint_RIMv3", "Fint_S20", "salR", "salS20")


def qnames():
    q = list(BASE_Q)
    for lab in LABELS:
        q += [f"{n}|{lab}" for n in CASE_Q]
    for lab in UNPAIRED_CASES:
        q += [f"u_dil_D|{lab}", f"u_int_D|{lab}"]
    return q


class DataError(Exception):
    """输入缺失或等待超时：退出 2。"""


class StopAsk(Exception):
    """设计要求停下、需人工处理的情形（如 L8/L46 直接 fCO2 与公开 pCO2 都缺）：退出 3。"""


_LOG_PATH = None


def log(msg):
    line = f"{datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')} [{os.getpid()}] {msg}"
    print(line, flush=True)
    if _LOG_PATH:
        try:
            with open(_LOG_PATH, "a", encoding="utf-8") as f:
                f.write(line + "\n")
        except OSError:
            pass


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jdump(obj, path):
    def clean(x):
        if isinstance(x, dict):
            return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        if isinstance(x, (np.floating, float)):
            x = float(x)
            return x if math.isfinite(x) else None
        if isinstance(x, (np.integer,)):
            return int(x)
        if isinstance(x, np.bool_):
            return bool(x)
        if isinstance(x, np.ndarray):
            return clean(x.tolist())
        return x
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(clean(obj), f, ensure_ascii=False, indent=1, allow_nan=False, default=str)
    os.replace(tmp, path)


def peak_rss_bytes():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(r if sys.platform == "darwin" else r * 1024)


def epoch_hour(d):
    if d.tzinfo is None:
        d = d.replace(tzinfo=UTC)
    return int(d.timestamp()) // 3600


def hour_dt(h):
    return datetime.fromtimestamp(h * 3600, tz=UTC)


def days_of_month(y, m):
    d0 = date(y, m, 1)
    d1 = date(y + (m == 12), m % 12 + 1, 1)
    return [d0 + timedelta(i) for i in range((d1 - d0).days)]


def _gsw():
    import gsw
    return gsw


# ======================================================================== 胶囊函数的逐式移植（数值核）
def d0_bilinear(u, r):
    """TB L51–L52、L62、L71：d0.interp(U=wind, R=precip)，xarray 默认 linear＝双线性，坐标外 NaN（K2/L3）。"""
    u = np.asarray(u, float)
    r = np.asarray(r, float)
    out = np.full(u.shape, np.nan)
    with np.errstate(invalid="ignore"):
        ok = (u >= U_GRID[0]) & (u <= U_GRID[-1]) & (r >= R_GRID[0]) & (r <= R_GRID[-1])
    if ok.any():
        uu, rr = u[ok], r[ok]
        i = np.clip(np.searchsorted(U_GRID, uu, side="right") - 1, 0, len(U_GRID) - 2)
        j = np.clip(np.searchsorted(R_GRID, rr, side="right") - 1, 0, len(R_GRID) - 2)
        tu = (uu - U_GRID[i]) / (U_GRID[i + 1] - U_GRID[i])
        tr = (rr - R_GRID[j]) / (R_GRID[j + 1] - R_GRID[j])
        out[ok] = ((1 - tu) * (1 - tr) * D0_TABLE[i, j] + (1 - tu) * tr * D0_TABLE[i, j + 1]
                   + tu * (1 - tr) * D0_TABLE[i + 1, j] + tu * tr * D0_TABLE[i + 1, j + 1])
    return out


def sc_co2(sst):
    """TB L106–L126 calculate_Sc_CO2。"""
    return SC_A - SC_B * (sst) + SC_C * (sst ** 2) - SC_D * (sst ** 3)


def weiss_ns(sst, S):
    """TB L151–L155：ln α_Weiss（mol/kg/atm）。"""
    SSTk = sst + 273.15
    ns = NSB1 + NSB2 * (SSTk / 100) + NSB3 * (SSTk / 100) ** 2
    ns = ns * S
    ns = ns + NSA1 + NSA2 * (100 / SSTk) + NSA3 * np.log(SSTk / 100)
    return ns


def weiss_slope(sst):
    """∂ ln α_Weiss / ∂S（TB L153 的系数），L27。"""
    SSTk = sst + 273.15
    return NSB1 + NSB2 * (SSTk / 100) + NSB3 * (SSTk / 100) ** 2


def fug_factor(t_c):
    """L46：fCO2/pCO2＝exp[(B＋2δ)·P/(R·T_K)]（Weiss 1974，与 CO2SYS／PyCO2SYS gas.fugacity_factor 同式；B、δ 单位 cm³/mol）。"""
    TK = np.asarray(t_c, float) + FUG_T0
    B = -1636.75 + 12.0408 * TK - 0.0327957 * TK ** 2 + 3.16528e-5 * TK ** 3
    delta = 57.7 - 0.118 * TK
    return np.exp((B + 2 * delta) * FUG_P_BAR / (FUG_R * TK))


def sol_exact(sst, S, lon, lat):
    """TB L130–L165 calculate_CO2_solubility_per_kg（mol/m³/µatm）；照胶囊把 SST 作 CT 传 gsw.rho。"""
    gsw = _gsw()
    a = np.exp(weiss_ns(sst, S))
    SA = gsw.SA_from_SP(S, 0, lon, lat)
    rho = gsw.rho(SA, sst, 0)
    return a * rho * (10 ** -6)


def sol_pair(sst, S, lon, lat):
    """稀释情形两侧溶解度（TB L500–L501）：αw＝sol(sst,S)，αa＝sol(sst−0.17,S)（不加 0.1，L30）；SA 只算一次，与分别调用逐值相同。"""
    gsw = _gsw()
    SA = gsw.SA_from_SP(S, 0, lon, lat)
    ta = sst - SKIN_DT
    aw = np.exp(weiss_ns(sst, S)) * gsw.rho(SA, sst, 0) * (10 ** -6)
    aa = np.exp(weiss_ns(ta, S)) * gsw.rho(SA, ta, 0) * (10 ** -6)
    return aw, aa


def k_harrison(sck, u10, R):
    """TB L343–L372 calculate_ktotal_Harrison2012（cm/h）；sck＝(Sc/600)^(−1/2)。"""
    kwind = 0.266 * (u10 ** 2)
    KEF_rain = 0.0112 * R
    krain = 63.02 * (KEF_rain ** (0.6242))
    ustar2 = ((1.03e-3 + 0.04e-3 * (u10 ** 1.48)) / (u10 ** 0.21)) * (u10 ** 2)
    ustar3 = ustar2 * (ustar2 ** (1 / 2))
    KEF_wind = ustar3 * 1.22
    beta = KEF_rain / KEF_wind
    return (kwind + (krain * (1 - np.exp(-0.3677 * beta)))) * sck


def day_coefs(T, S0, lon, lat, r_sa):
    """逐日缓存（SST、S0 日内常数，L31）：两侧溶解度、Schmidt、淡水溶解度与 L25/L27 展开系数。"""
    gsw = _gsw()
    c = {}
    with np.errstate(all="ignore"):
        Sc = sc_co2(T)
        c["sck"] = (Sc / 600) ** (-1 / 2)
        c["g_kw"] = GAMMA_K * c["sck"]                            # TB L183：Γ*(Sc/600)**(-1/2) * u10**2
        Ta = T - SKIN_DT
        SA0 = gsw.SA_from_SP(S0, 0, lon, lat)
        SA0p = gsw.SA_from_SP(S0 + SKIN_DS, 0, lon, lat)
        rho_w = gsw.rho(SA0, T, 0)
        rho_a = gsw.rho(SA0, Ta, 0)
        c["aw0"] = np.exp(weiss_ns(T, S0)) * rho_w * (10 ** -6)                                  # TB L497
        c["aa0"] = np.exp(weiss_ns(Ta, S0 + SKIN_DS)) * gsw.rho(SA0p, Ta, 0) * (10 ** -6)        # TB L498
        c["aa0d"] = np.exp(weiss_ns(Ta, S0)) * rho_a * (10 ** -6)                                # TB L501 在 S=S0
        c["af"] = sol_exact(T, 0, 0, 0)                                                          # TB L521
        for side, tt, rho in (("w", T, rho_w), ("a", Ta, rho_a)):
            d1 = gsw.rho_first_derivatives(SA0, tt, 0)[0]
            d2 = gsw.rho_second_derivatives(SA0, tt, 0)[0]
            q = d1 / rho
            c["K" + side + "1"] = weiss_slope(tt) + q * r_sa
            c["K" + side + "2"] = (d2 / rho - q * q) * r_sa * r_sa
    return c


def alpha_expand(S, S0, co, T, lon, lat):
    """L25/L27：α(S)＝α(S0)·exp(K1·ΔS＋½K2·ΔS²)，两侧各一；|ΔS|>3 psu 或 min(S,S0)<15 psu 改精确 gsw。返回 (αw, αa, 兜底掩码)。"""
    with np.errstate(all="ignore"):
        d = S - S0
        d2 = d * d
        aw = co["aw0"] * np.exp(co["Kw1"] * d + 0.5 * co["Kw2"] * d2)
        aa = co["aa0d"] * np.exp(co["Ka1"] * d + 0.5 * co["Ka2"] * d2)
        fb = np.isfinite(S) & ((np.abs(d) > L25_DMAX) | (np.minimum(S, S0) < L25_SMIN))
    if fb.any():
        aw[fb], aa[fb] = sol_pair(T[fb], S[fb], lon[fb], lat[fb])
    return aw, aa, fb


def _guess_bounds(points, bound_position=0.5):
    """TB L557–L580（dennissergeev gist）。"""
    diffs = np.diff(points)
    diffs = np.insert(diffs, 0, diffs[0])
    diffs = np.append(diffs, diffs[-1])
    min_bounds = points - diffs[:-1] * bound_position
    max_bounds = points + diffs[1:] * (1 - bound_position)
    return np.array([min_bounds, max_bounds]).transpose()


def grid_cell_areas(lon1d, lat1d, radius=EARTH_RADIUS):
    """TB L583–L648：球面格点面积 (nlat, nlon)，m²。"""
    lon_b = np.deg2rad(_guess_bounds(np.asarray(lon1d, float)))
    lat_b = np.deg2rad(_guess_bounds(np.asarray(lat1d, float)))
    ylen = np.sin(lat_b[:, 1].astype(np.float64)) - np.sin(lat_b[:, 0].astype(np.float64))
    xlen = lon_b[:, 1].astype(np.float64) - lon_b[:, 0].astype(np.float64)
    return np.abs((radius ** 2) * np.outer(ylen, xlen))


def glodap_slopes(ta, tco2, sss, sst, chunk=4000):
    """TB L414–L465 calculate_fco2_dilution_per_psu 的向量化（main L41–L51 逐点循环）：返回 float32 斜率（main L39），任一输入 NaN → NaN。"""
    import PyCO2SYS as pyco2
    shape = np.shape(ta)
    ta, tco2, sss, sst = (np.asarray(x, float).ravel() for x in (ta, tco2, sss, sst))
    out = np.full(ta.shape, np.nan, dtype=np.float32)
    ok = np.isfinite(ta) & np.isfinite(tco2) & np.isfinite(sss) & np.isfinite(sst)
    ii = np.flatnonzero(ok)
    res = np.full(len(ii), np.nan)

    def fit(TA, DIC, S, T, sals):
        SR = sals / S[:, None]                                     # TB L439
        talk_in = TA[:, None] * SR                                 # TB L440
        dic_in = DIC[:, None] * SR + DIC_RAIN * (1 - SR)           # TB L441
        r = pyco2.sys(par1=talk_in, par2=dic_in, par1_type=1, par2_type=2, salinity=sals,
                      temperature=np.broadcast_to(T[:, None], sals.shape).copy(), opt_k_carbonic=K_CARBONIC)
        y = np.asarray(r["fCO2"], float)
        xm = sals - sals.mean(axis=1, keepdims=True)               # scipy.stats.linregress：cov(bias=1) 的 ssxym/ssxm
        ym = y - y.mean(axis=1, keepdims=True)
        return (xm * ym).mean(axis=1) / (xm * xm).mean(axis=1)

    S = sss[ii]
    hi = np.flatnonzero(S > 9)
    for a in range(0, len(hi), chunk):
        k = hi[a:a + chunk]
        sals = S[k][:, None] - np.arange(10)[None, :].astype(float)   # TB L434
        res[k] = fit(ta[ii][k], tco2[ii][k], S[k], sst[ii][k], sals)
    for k in np.flatnonzero(~(S > 9)):                               # TB L436（低盐）
        sals = np.arange(S[k], 0, -1)[None, :]
        if sals.shape[1] >= 2:
            res[k] = fit(ta[ii][k:k + 1], tco2[ii][k:k + 1], S[k:k + 1], sst[ii][k:k + 1], sals)[0]
    out[ii] = res.astype(np.float32)
    return out.reshape(shape)


# ======================================================================== 网格工具（插值、分区）
def _lin_idx(x, xn):
    """x 升序；返回 (i, w, valid)：xn 落在 [x_i, x_{i+1}]，w＝相对位置；越界 invalid（→NaN，同 scipy interp1d bounds_error=False）。"""
    xn = np.asarray(xn, float)
    valid = (xn >= x[0]) & (xn <= x[-1])
    i = np.clip(np.searchsorted(x, xn, side="right") - 1, 0, len(x) - 2)
    w = (xn - x[i]) / (x[i + 1] - x[i])
    w = np.where(valid, w, np.nan)
    return i, w, valid


def _nearest_idx(x, xn):
    """scipy interp1d(kind='nearest')：中点处取低侧；越界 invalid。"""
    xn = np.asarray(xn, float)
    mid = (x[:-1] + x[1:]) / 2.0
    i = np.clip(np.searchsorted(mid, xn, side="left"), 0, len(x) - 1)
    valid = (xn >= x[0]) & (xn <= x[-1])
    return i, valid


class Bilin:
    """源规则网格（纬度可降序、可不等距；经度周期）→ 目标（域的唯一行纬度 × 唯一列经度）的可分离双线性，NaN 传播（L31）。"""

    def __init__(self, src_lat, src_lon, tgt_lat, tgt_lon):
        lat = np.asarray(src_lat, float)
        self.flip = bool(lat[0] > lat[-1])
        if self.flip:
            lat = lat[::-1]
        lon = np.asarray(src_lon, float) % 360.0
        ulon, first = np.unique(lon, return_index=True)          # 去掉 −180/180 这类重复列
        lon_ext = np.r_[ulon, ulon[0] + 360.0]
        col_ext = np.r_[first, first[0]]
        tl = np.asarray(tgt_lon, float) % 360.0
        tl = np.where(tl < lon_ext[0], tl + 360.0, tl)
        self.i0, self.wy, self.vy = _lin_idx(lat, tgt_lat)
        j, self.wx, self.vx = _lin_idx(lon_ext, tl)
        self.c0, self.c1 = col_ext[j], col_ext[j + 1]
        self.src_shape = (len(lat), len(lon))

    def __call__(self, F2d):
        F = np.asarray(F2d, float)
        if F.shape != self.src_shape:
            raise ValueError(f"插值源形状 {F.shape} ≠ {self.src_shape}")
        if self.flip:
            F = F[::-1]
        wy = self.wy[:, None]
        Fy = F[self.i0] * (1 - wy) + F[self.i0 + 1] * wy
        G = Fy[:, self.c0] * (1 - self.wx) + Fy[:, self.c1] * self.wx
        return G


class Nearest:
    """main L52–L60：经度 %360、排序、补 0 与 360 两列（分别复制首末列），再做最近邻（interp_like method='nearest'）。"""

    def __init__(self, src_lat, src_lon, tgt_lat, tgt_lon):
        lat = np.asarray(src_lat, float)
        self.lat_order = np.argsort(lat, kind="stable")
        lat_s = lat[self.lat_order]
        lon = np.asarray(src_lon, float) % 360.0
        order = np.argsort(lon, kind="stable")
        lon_s = lon[order]
        lon_aug = np.r_[0.0, lon_s, 360.0]
        col_aug = np.r_[order[0], order, order[-1]]
        ri, self.vr = _nearest_idx(lat_s, tgt_lat)
        ci, self.vc = _nearest_idx(lon_aug, np.asarray(tgt_lon, float) % 360.0)
        self.rows = self.lat_order[ri]
        self.cols = col_aug[ci]

    def __call__(self, F2d):
        G = np.asarray(F2d)[self.rows][:, self.cols].astype(np.float64)
        G[~self.vr, :] = np.nan
        G[:, ~self.vc] = np.nan
        return G


class Seg:
    """已按分区升序的一段像元：分区连续段的起点与分区号（供 np.add.reduceat）。"""

    def __init__(self, regs):
        regs = np.asarray(regs)
        self.n = len(regs)
        if self.n:
            self.starts = np.flatnonzero(np.r_[True, regs[1:] != regs[:-1]])
            self.ids = regs[self.starts]
        else:
            self.starts = self.ids = np.zeros(0, np.int64)


class Domain:
    """一组海洋像元（全场或自测小块）：CMORPH 网格行列、面积、β、分区。像元顺序＝分区升序（分区内按平面索引）。"""

    def __init__(self, lat1d, lon1d, row, col, area, beta, reg, nreg, r_sa=None):
        self.lat1d = np.asarray(lat1d, float)
        self.lon1d = np.asarray(lon1d, float)
        self.row = np.asarray(row, np.int32)
        self.col = np.asarray(col, np.int32)
        self.N = len(self.row)
        self.flat = self.row.astype(np.int64) * len(self.lon1d) + self.col
        self.area = np.asarray(area, float)
        self.beta = np.asarray(beta, np.float32)
        self.reg = np.asarray(reg, np.int32)
        self.nreg = int(nreg)
        if self.N and np.any(np.diff(self.reg) < 0):
            raise ValueError("Domain 像元须按分区升序排列")
        self.seg_full = Seg(self.reg)
        self.reg_count = np.bincount(self.reg, minlength=self.nreg).astype(float)
        self.urow, ri = np.unique(self.row, return_inverse=True)
        self.ucol, ci = np.unique(self.col, return_inverse=True)
        self.ri = ri.astype(np.int32)
        self.ci = ci.astype(np.int32)
        if r_sa is None:
            gsw = _gsw()
            lon, lat = self.lon(), self.lat()
            r_sa = gsw.SA_from_SP(36.0, 0, lon, lat) - gsw.SA_from_SP(35.0, 0, lon, lat)   # SA 对 SP 的仿射斜率（L27）
        self.r_sa = np.asarray(r_sa, float)

    def lat(self):
        return self.lat1d[self.row]

    def lon(self):
        return self.lon1d[self.col]

    def gather(self, G):
        return G[self.ri, self.ci]

    def from_native(self, F2d):
        return np.asarray(F2d).ravel()[self.flat]

    def tgt(self):
        return self.lat1d[self.urow], self.lon1d[self.ucol]


# ======================================================================== RIM-3 稀疏历史（L29）
class RimRing:
    """48 槽环形缓冲：每槽只存「非恒等」条目（像元位置、v＝c1·IRR/(d0·√Kz)、g＝−z²/(4Kz)）；
    历史因子 1/(1+v·tc/√ti)（＝胶囊 d0/(d0+c1·IRR·tc/√(Kz·ti))），z=5 m 另乘 exp(g/ti)。"""

    def __init__(self, N):
        self.N = N
        self.slots = [None] * 48
        self.n = 0
        self.entry_counts = []

    @staticmethod
    def entries(P, U):
        P = np.asarray(P, float)
        U = np.asarray(U, float)
        allnan = bool(len(P) > 0 and not np.isfinite(P).any())
        with np.errstate(all="ignore"):
            nz = ~((P == 0) & (U > 0) & (U <= U_GRID[-1]))
            pos = np.flatnonzero(nz).astype(np.int32)
            p, u = P[pos], U[pos]
            irr = p / 1000 / 3600                                  # TB L59、L67
            kz = KZ_COEF * (u ** 2)                                # TB L61、L69
            d0 = d0_bilinear(u, p)                                 # TB L62、L71
            v = (C1 * irr) / (d0 * np.sqrt(kz))                    # TB L65（除以 d0 后）
            g = -(Z5 ** 2) / (4.0 * kz)
            cterm = C2 * irr * 1800 / (np.sqrt(kz))                # TB L73
            cur0 = d0 / (d0 + cterm)
            cur5 = d0 / (d0 + cterm * np.exp(-(Z5 ** 2) / (4.0 * kz * T_CURRENT_S)))   # K3：当前项 t=1 s
        return {"pos": pos, "v": v, "g": g, "cur0": cur0, "cur5": cur5, "allnan": allnan}

    def history(self):
        N = self.N
        pos_l, l0, l5 = [], [], []
        allnan = False
        with np.errstate(all="ignore"):
            for a in range(1, 49):
                sl = self.slots[(self.n - a) % 48]
                if sl is None:
                    raise RuntimeError("RIM 历史不足 48 个半步（起转未完成）")
                if sl["allnan"]:
                    allnan = True
                    continue
                if len(sl["pos"]) == 0:
                    continue
                x = sl["v"] * AGE_W[a - 1]
                pos_l.append(sl["pos"])
                l0.append(-np.log1p(x))
                l5.append(-np.log1p(x * np.exp(sl["g"] / AGE_TI[a - 1])))
        hmask = np.zeros(N, bool)
        if pos_l:
            pos = np.concatenate(pos_l)
            logH = np.bincount(pos, weights=np.concatenate(l0), minlength=N)
            logH5 = np.bincount(pos, weights=np.concatenate(l5), minlength=N)
            hmask[pos] = True
        else:
            logH = np.zeros(N)
            logH5 = np.zeros(N)
        if allnan:
            logH[:] = np.nan
            logH5[:] = np.nan
            hmask[:] = True
        return logH, logH5, hmask

    def push(self, ent):
        if ent["allnan"]:
            self.slots[self.n % 48] = {"allnan": True}
        else:
            self.slots[self.n % 48] = {"pos": ent["pos"], "v": ent["v"], "g": ent["g"], "allnan": False}
        self.entry_counts.append(len(ent["pos"]))
        self.n += 1


# ======================================================================== 通量引擎
class Engine:
    """逐半步：稠密风驱通量＋有雨像元（D）上的沉降、湍流、24 情形稀释/界面通量，按分区累计（不落全球时间序列）。"""

    def __init__(self, dom, chunk=CHUNK, keep=False):
        self.dom = dom
        self.N = dom.N
        self.ring = RimRing(dom.N)
        self.chunk = int(chunk)
        self.keep = keep
        self.l25_check = False
        self.qn = qnames()
        self.qi = {q: i for i, q in enumerate(self.qn)}
        self.acc = np.zeros((len(self.qn), dom.nreg))
        self.diag = {"l25_max_rel": {}, "l25_n_checked": 0, "D_frac": [], "C_entries": [], "steps": 0}
        self.tim = {}
        self.sf = self.fa = None                                    # 引擎实际用的 fCO2（L46 下为逐日换算值）
        self.sf_in = self.fa_in = None                              # 月场输入：fCO2（direct）或 pCO2（L46）
        self.w_pco2 = False
        self.fug = {"sea": [np.inf, -np.inf, 0.0, 0], "air": [np.inf, -np.inf, 0.0, 0]}   # min, max, 和, 个数（L46）
        self.co = None
        self.kept = []

    def _t(self, key, t0):
        t1 = time.perf_counter()
        self.tim[key] = self.tim.get(key, 0.0) + t1 - t0
        return t1

    def _add(self, q, seg, x):
        if seg.n:
            self.acc[self.qi[q], seg.ids] += np.add.reduceat(np.asarray(x, float), seg.starts)

    def _add_rows(self, qs, seg, rows):
        if seg.n:
            M = np.vstack([np.asarray(r, float) for r in rows])
            sums = np.add.reduceat(M, seg.starts, axis=1)
            self.acc[np.ix_([self.qi[q] for q in qs], seg.ids)] += sums

    # ---- 场 ----
    def set_month(self, sf, fa, pco2=False):
        """sf/fa：像元上的 Watson 月值；pco2=True 时是 pCO2（L46），在 _mix 里逐日按 SST 换成 fCO2。"""
        self.sf_in = np.asarray(sf, float)
        self.fa_in = np.asarray(fa, float)
        self.w_pco2 = bool(pco2)
        if self.co is not None:
            self._mix()

    def set_day(self, sst, ice_pct, s0):
        self.co = self.dC0 = self.dC0d = None                       # 先释放前一日缓存，避免两日并存
        self.T = np.asarray(sst, float)
        self.S0 = np.asarray(s0, float)
        self.icef = 1 - (np.asarray(ice_pct, float) / 100)          # TB L223
        self.co = day_coefs(self.T, self.S0, self.dom.lon(), self.dom.lat(), self.dom.r_sa)
        if self.sf_in is not None:
            self._mix()

    def _fug_acc(self, side, ff, base):
        m = np.isfinite(ff) & np.isfinite(base)
        if m.any():
            v = ff[m]
            a = self.fug[side]
            a[0], a[1] = min(a[0], float(v.min())), max(a[1], float(v.max()))
            a[2] += float(v.sum())
            a[3] += int(m.sum())

    def fug_summary(self):
        if not self.w_pco2:
            return None
        return {k: {"min": v[0], "max": v[1], "sum": v[2], "n": v[3], "mean": v[2] / v[3] if v[3] else None}
                for k, v in self.fug.items()}

    def _mix(self):
        c = self.co
        if self.w_pco2:                                              # L46：海侧 T＝溶解度用的 SST，大气侧 T＝SST−SKIN_DT（TB L498）
            with np.errstate(all="ignore"):
                fw = fug_factor(self.T)
                fa_ = fug_factor(self.T - SKIN_DT)
            self._fug_acc("sea", fw, self.sf_in)
            self._fug_acc("air", fa_, self.fa_in)
            self.sf = self.sf_in * fw
            self.fa = self.fa_in * fa_
            del fw, fa_
        else:
            self.sf, self.fa = self.sf_in, self.fa_in
        with np.errstate(all="ignore"):
            self.dC0 = (self.sf * c["aw0"]) - (self.fa * c["aa0"])     # TB L223（风驱／湍流）
            self.dC0d = (self.sf * c["aw0"]) - (self.fa * c["aa0d"])   # 稀释情形在 S=S0（L30）

    # ---- 步 ----
    def push_only(self, P, U):
        self.ring.push(self.ring.entries(P, U))

    def step(self, P, U):
        dom, N = self.dom, self.N
        P = np.asarray(P, float)
        U = np.asarray(U, float)
        t0 = time.perf_counter()
        ent = self.ring.entries(P, U)
        t0 = self._t("entries", t0)
        logH, logH5, Dmask = self.ring.history()
        t0 = self._t("history", t0)
        Dmask[ent["pos"]] = True
        c = self.co
        seg = dom.seg_full
        A = dom.area
        with np.errstate(all="ignore"):
            kw = c["g_kw"] * (U ** 2)                                 # TB L183
            Trw = 12 * (kw * TCONV_K)                                 # TB L218
            Fwind = Trw * self.dC0 * self.icef                        # TB L223
            okW = np.isfinite(Fwind)
            Fdil0 = Trw * self.dC0d * self.icef
            nonD = ~Dmask
            okS = nonD & okW & np.isfinite(Fdil0)
            self.acc[self.qi["n_ocean"]] += dom.reg_count
            self._add("n_W", seg, okW)
            self._add("W_all", seg, np.where(okW, A * Fwind, 0.0))
            self._add("skin", seg, np.where(okS, A * (Fdil0 - Fwind), 0.0))
            self._add("u_nonD_dil", seg, np.where(nonD & np.isfinite(Fdil0), A * Fdil0, 0.0))
            self._add("u_nonD_W", seg, np.where(nonD & okW, A * Fwind, 0.0))
            self._add("n_Pnan", seg, ~np.isfinite(P))
            self._add("n_Pgt200", seg, P > R_GRID[-1])
            self._add("n_U0", seg, U == 0)
            self._add("n_Unan", seg, ~np.isfinite(U))                 # L47：只计数（RIM 历史 NaN 型条目的来源之一）
            self._add("n_PU0", seg, (P == 0) & (U == 0))              # L47：只计数（0/0）
        if self.keep:
            K = {k: np.full(N, np.nan) for k in KEEP_FIELDS}
            with np.errstate(all="ignore"):
                K["Fwind"][:] = Fwind
                K["Fdep"][:] = 12 * (P * TCONV_R) * c["af"] * self.fa
                K["Fturb"][:] = Fwind
                for k in ("Fdil_RIMv3", "Fint_RIMv3", "Fdil_S20", "Fint_S20"):
                    K[k][:] = Fdil0
                K["salR"][:] = self.S0
                K["salS20"][:] = np.where(self.S0 > 0, self.S0, 0.0)          # 无雨像元胶囊 S20＝S0−0.0 再 where(>0, 0)：S0 为 NaN 时得 0（仅自测比对用）
            self._K = K
        t0 = self._t("dense", t0)
        D = np.flatnonzero(Dmask)
        cz0 = np.ones(len(D))
        cz5 = np.ones(len(D))
        k = np.searchsorted(D, ent["pos"])
        cz0[k] = ent["cur0"]
        cz5[k] = ent["cur5"]
        for s in range(0, len(D), self.chunk):
            self._chunk(D[s:s + self.chunk], cz0[s:s + self.chunk], cz5[s:s + self.chunk], P, U, logH, logH5, Trw, Fwind)
        t0 = time.perf_counter()
        self.ring.push(ent)
        self.diag["D_frac"].append(len(D) / max(N, 1))
        self.diag["C_entries"].append(len(ent["pos"]))
        self.diag["steps"] += 1
        if self.keep:
            self.kept.append(self._K)
        self._t("push", t0)

    def _chunk(self, idx, cz0, cz5, P, U, logH, logH5, Trw_f, Fwind_f):
        dom, c = self.dom, self.co
        t0 = time.perf_counter()
        seg = Seg(dom.reg[idx])
        A = dom.area[idx]
        Fw = Fwind_f[idx]
        okW = np.isfinite(Fw)
        S0, T, icef = self.S0[idx], self.T[idx], self.icef[idx]
        sf, fa, beta = self.sf[idx], self.fa[idx], dom.beta[idx]
        aw0, aa0d, dC0 = c["aw0"][idx], c["aa0d"][idx], self.dC0[idx]
        p, u, Trw = P[idx], U[idx], Trw_f[idx]
        lon, lat = dom.lon1d[dom.col[idx]], dom.lat1d[dom.row[idx]]
        with np.errstate(all="ignore"):
            # ---- 湍流、沉降（s 无关） ----
            needH = ~((p == 0) & (u > 0))            # 其余 P==0 像元 Harrison 式退化为 kw（L30）
            TrH = Trw.copy()
            if needH.any():
                TrH[needH] = 12 * (k_harrison(c["sck"][idx][needH], u[needH], p[needH]) * TCONV_K)
            Fturb = TrH * dC0 * icef                                    # TB L532–L534
            Fdep = 12 * (p * TCONV_R) * c["af"][idx] * fa               # TB L402、L408
            okT = okW & np.isfinite(Fturb)
            okP = okW & np.isfinite(Fdep)
            self._add_rows(["n_D", "n_C", "turb", "Wm_turb", "nm_turb", "u_turb_D", "dep", "Wm_dep", "nm_dep", "u_dep"], seg,
                           [okW, okW & (p != 0),
                            np.where(okT, A * (Fturb - Fw), 0.0), np.where(okW & ~okT, A * Fw, 0.0), okW & ~okT,
                            np.where(np.isfinite(Fturb), A * Fturb, 0.0),
                            np.where(okP, A * Fdep, 0.0), np.where(okW & ~okP, A * Fw, 0.0), okW & ~okP,
                            np.where(np.isfinite(Fdep), A * Fdep, 0.0)])
            # ---- RIM s=1：精确 gsw（L25） ----
            logHd = logH[idx]
            S_R = S0 * (np.exp(logHd) * cz0)                            # TB L76：sal*(prod(prior)*current)
            okR = np.isfinite(S_R)
            missR = ~okR                                                # L29：RIM 全族共同缺测集
            self._add("n_poisonR", seg, okW & missR)
        t0 = self._t("d_base", t0)
        aw_R = np.full(len(idx), np.nan)
        aa_R = np.full(len(idx), np.nan)
        if okR.any():
            aw_R[okR], aa_R[okR] = sol_pair(T[okR], S_R[okR], lon[okR], lat[okR])
        t0 = self._t("d_exact", t0)
        co = None
        logH5d = None
        for ci_, (kind, s) in enumerate(CASES):
            lab = LABELS[ci_]
            nclip = None
            fb = None
            miss = None
            with np.errstate(all="ignore"):
                if kind in ("whole", "hist") and s == 1.0:
                    S, aw, aa = S_R, aw_R, aa_R
                    miss = missR
                elif kind == "whole" and s == 0.0:
                    S, aw, aa = S0, aw0, aa0d
                    miss = missR
                elif kind == "s20":
                    S = S0 + S20_A * p * (u ** -S20_B)                  # TB L99–L101
                    S = np.where(S > 0, S, 0.0)                         # main L79：where(sal_S20>0, other=0)
                    miss = ~np.isfinite(p)                              # L35/L12
                    ch = np.isfinite(S) & (S != S0)
                    aw, aa = aw0.copy(), aa0d.copy()
                    if ch.any():
                        aw[ch], aa[ch] = sol_pair(T[ch], S[ch], lon[ch], lat[ch])
                else:
                    if co is None:
                        co = {k: c[k][idx] for k in ("aw0", "aa0d", "Kw1", "Kw2", "Ka1", "Ka2")}
                    if kind == "whole":
                        S = S0 + s * (S_R - S0)                         # S_dil＝S0＋s·(S_RIM3−S0)
                        neg = S < 0
                        nclip = okW & neg
                        S = np.where(neg, 0.0, S)                       # L28
                    elif kind == "hist":
                        E = np.exp(s * logHd) if s != 0.0 else 1.0      # S3：ln Π(历史)×s，当前项不缩放
                        S = S0 * (E * cz0)
                    else:  # z5
                        if logH5d is None:
                            logH5d = logH5[idx]
                        S = S0 * (np.exp(logH5d) * cz5)
                    miss = missR
                    aw, aa, fb = alpha_expand(S, S0, co, T, lon, lat)
                    if self.l25_check:
                        tl = time.perf_counter()
                        self._l25_compare(lab, S, aw, aa, fb, T, lon, lat)
                        self._t("l25_check", tl)                        # bench 外推时从单步耗时里扣除
                dS = S - S0                                             # main L80–L81
                Xc = ((sf + beta * dS) * aw) - (fa * aa)                # main L84–L85、TB L223
                Fdil = Trw * Xc * icef                                  # TB L528–L529
                Fint = TrH * Xc * icef                                  # TB L537–L538
                Y = beta * dS * aw * icef                               # L39
                okd = okW & np.isfinite(Fdil) & ~miss
                oki = okW & np.isfinite(Fint) & ~miss
                rows = [np.where(okd, A * (Fdil - Fw), 0.0), np.where(oki, A * (Fint - Fw), 0.0),
                        np.where(okd, A * Trw * Y, 0.0), np.where(oki, A * TrH * Y, 0.0),
                        np.where(okW & ~okd, A * Fw, 0.0), okW & ~okd, np.where(okW & ~oki, A * Fw, 0.0), okW & ~oki,
                        (okW & fb) if fb is not None else np.zeros(len(idx)),
                        nclip if nclip is not None else np.zeros(len(idx))]
                qs = [f"{n}|{lab}" for n in CASE_Q]
                if lab in UNPAIRED_CASES:
                    rows += [np.where(np.isfinite(Fdil) & ~miss, A * Fdil, 0.0), np.where(np.isfinite(Fint) & ~miss, A * Fint, 0.0)]
                    qs += [f"u_dil_D|{lab}", f"u_int_D|{lab}"]
                self._add_rows(qs, seg, rows)
                if self.keep:
                    K = self._K
                    if lab == "whole_1":
                        K["Fdil_RIMv3"][idx], K["Fint_RIMv3"][idx], K["salR"][idx] = Fdil, Fint, S
                        K["Fdep"][idx], K["Fturb"][idx] = Fdep, Fturb
                    elif lab == "s20_1":
                        K["Fdil_S20"][idx], K["Fint_S20"][idx], K["salS20"][idx] = Fdil, Fint, S
        self._t("d_cases", t0)

    def _l25_compare(self, lab, S, aw, aa, fb, T, lon, lat):
        ok = np.isfinite(S) & ~fb
        if not ok.any():
            return
        ew, ea = sol_pair(T[ok], S[ok], lon[ok], lat[ok])
        with np.errstate(all="ignore"):
            r = np.nanmax(np.r_[np.abs(aw[ok] / ew - 1), np.abs(aa[ok] / ea - 1)])
        d = self.diag["l25_max_rel"]
        d[lab] = max(d.get(lab, 0.0), float(r) if np.isfinite(r) else 0.0)
        self.diag["l25_n_checked"] += int(ok.sum())


# ======================================================================== 读取层（数据布局假设全在这里）
def _nc_var(ds, names):
    low = {k.lower(): k for k in ds.variables}
    for n in names:
        if n in ds.variables:
            return ds.variables[n]
        if n.lower() in low:
            return ds.variables[low[n.lower()]]
    return None


def _nc_coord(ds, kind):
    names = ("lat", "latitude", "y", "Latitude") if kind == "lat" else ("lon", "longitude", "x", "Longitude")
    v = _nc_var(ds, names)
    if v is None:
        raise KeyError(f"{ds.filepath()} 无 {kind} 坐标")
    return np.asarray(v[:], float)


def _nc_take2d(ds, v, fixed=None):
    """取 (lat, lon) 二维切片：lat/lon 全取；fixed 指定的维取给定下标；depth/lev 维取坐标为 0 的层（无坐标取 0）；其余取 0。"""
    idx, order = [], []
    for d in v.dimensions:
        dl = d.lower()
        if fixed and d in fixed:
            idx.append(fixed[d])
        elif "lat" in dl or dl == "y":
            idx.append(slice(None))
            order.append("lat")
        elif "lon" in dl or dl == "x":
            idx.append(slice(None))
            order.append("lon")
        elif "depth" in dl or "lev" in dl or dl == "z":
            k = 0
            if d in ds.variables:
                cz = np.asarray(ds.variables[d][:], float)
                hit = np.flatnonzero(cz == 0)
                k = int(hit[0]) if len(hit) else 0
            idx.append(k)
        else:
            idx.append(0)
    a = np.ma.filled(np.ma.asarray(v[tuple(idx)]).astype(np.float64), np.nan)
    if order == ["lon", "lat"]:
        a = a.T
    return a


def _decode_ym(ds, v):
    """返回变量时间维每个下标的 (年, 月)；units 为 'months since' 时手算（cftime 不支持标准历法的月单位）。"""
    tdim = [d for d in v.dimensions if "time" in d.lower() or d.lower() in ("t", "month", "mtime")]
    if not tdim:
        return None, None
    tname = tdim[0]
    tv = ds.variables.get(tname)
    n = ds.dimensions[tname].size
    if tv is None:
        return tname, None
    vals = np.asarray(tv[:], float)
    units = getattr(tv, "units", "")
    try:
        import netCDF4
        dts = netCDF4.num2date(vals, units, calendar=getattr(tv, "calendar", "standard"))
        return tname, [(int(d.year), int(d.month)) for d in dts]
    except Exception:
        m = re.match(r"\s*months since\s+(\d{4})-(\d{1,2})", units)
        if m:
            y0, m0 = int(m.group(1)), int(m.group(2))
            out = []
            for x in vals:
                k = y0 * 12 + (m0 - 1) + int(round(x))
                out.append((k // 12, k % 12 + 1))
            return tname, out
    return tname, [None] * n


class Layout:
    """输入布局（根＝--cache-root，约定见 p3_fetch.py 文件头）。所有「文件在哪、叫什么、变量名」的假设集中在本类；
    数据落地后若实际布局不同，只改这里。"""

    def __init__(self, root):
        self.root = os.path.abspath(root)
        self._cm_grid = None

    def path(self, *a):
        return os.path.join(self.root, *a)

    def _glob(self, pattern):
        return sorted(glob.glob(self.path(pattern), recursive=True))

    # ---- 就绪 ----
    def month_ready(self, ym):
        p = self.path("month_ready", f"{ym}.json")
        if not os.path.exists(p):
            return None
        with open(p, encoding="utf-8") as f:
            return json.load(f)

    def ready(self, ym, need_next=True):
        """L40：本月 month_ready；m>1 另需前月 month_ready（起转日与前月 ERA5）；m<12 另需下月 ws10 文件（月末 HH:30；
        bench 只跑首日时不需要）。"""
        y, m = divmod(ym, 100)
        if self.month_ready(ym) is None:
            return False, f"month_ready/{ym}.json 未写"
        if m > 1 and self.month_ready(ym - 1) is None:
            return False, f"起转需要 month_ready/{ym - 1}.json"
        if need_next and m < 12 and not os.path.exists(self.era5_path(ym + 1)):
            return False, f"月末 HH:30 需要 {os.path.relpath(self.era5_path(ym + 1), self.root)}"
        return True, "ok"

    # ---- CMORPH ----
    def cmorph_path(self, dt):
        p = self.path("cmorph", f"{dt:%Y}", f"{dt:%m}", f"CMORPH_V1.0_ADJ_8km-30min_{dt:%Y%m%d%H}.nc")
        if os.path.exists(p):
            return p
        hits = self._glob(os.path.join("cmorph", "**", f"*_{dt:%Y%m%d%H}.nc"))
        return hits[0] if hits else None

    def cmorph_grid(self, dt=None):
        if self._cm_grid is None:
            import netCDF4
            p = self.cmorph_path(dt) if dt else None
            if p is None:
                hits = self._glob(os.path.join("cmorph", "**", "CMORPH_*.nc"))
                if not hits:
                    raise DataError("cmorph/ 下没有文件")
                p = hits[0]
            with netCDF4.Dataset(p) as ds:
                self._cm_grid = (np.asarray(ds["lat"][:], float), np.asarray(ds["lon"][:], float))
        return self._cm_grid

    def read_cmorph(self, dt):
        """返回 ([P(HH:00), P(HH:30)] 两个 float32 二维场，NaN＝缺测, 原因)；缺文件／时间轴不符 → (None, 原因)（L35、K8）。"""
        import netCDF4
        p = self.cmorph_path(dt)
        if p is None:
            return None, "absent"
        h = epoch_hour(dt)
        try:
            with netCDF4.Dataset(p) as ds:
                t = [int(x) for x in ds["time"][:]]
                if t != [h * 3600, h * 3600 + 1800]:
                    return None, f"time_mismatch:{t}"
                v = ds["cmorph"]
                v.set_auto_maskandscale(True)
                out = [np.ma.filled(np.ma.asarray(v[k]).astype(np.float32), np.nan) for k in (0, 1)]
        except (OSError, RuntimeError, KeyError, IndexError) as e:
            return None, f"read_error:{type(e).__name__}"
        return out, None

    # ---- ERA5 ws10 ----
    def era5_path(self, ym):
        return self.path("era5_ws10", f"{ym}.nc")

    # ---- OISST ----
    def oisst_path(self, d):
        p = self.path("oisst", f"oisst-avhrr-v02r01.{d:%Y%m%d}.nc")
        if os.path.exists(p):
            return p
        hits = self._glob(os.path.join("oisst", "**", f"*{d:%Y%m%d}*.nc"))
        return hits[0] if hits else None

    def read_oisst(self, d):
        import netCDF4
        p = self.oisst_path(d)
        if p is None:
            raise DataError(f"OISST {d} 缺")
        with netCDF4.Dataset(p) as ds:
            lat, lon = _nc_coord(ds, "lat"), _nc_coord(ds, "lon")
            sst = _nc_take2d(ds, _nc_var(ds, ("sst",)))
            iv = _nc_var(ds, ("ice", "icec"))
            ice = _nc_take2d(ds, iv) if iv is not None else np.full(sst.shape, np.nan)
        return lat, lon, sst, ice

    # ---- HYCOM ----
    def hycom_path(self, d):
        p = self.path("hycom_sss", f"{d:%Y%m%d}.nc")
        return p if os.path.exists(p) else None

    def read_hycom(self, d):
        import netCDF4
        p = self.hycom_path(d)
        if p is None:
            raise DataError(f"HYCOM {d} 缺")
        with netCDF4.Dataset(p) as ds:
            v = _nc_var(ds, ("salinity", "sss", "so", "salt"))
            if v is None:
                raise DataError(f"{p} 无盐度变量")
            return _nc_coord(ds, "lat"), _nc_coord(ds, "lon"), _nc_take2d(ds, v)

    # ---- tar 兜底（下载侧应已解出成员） ----
    def _members(self, sub, want, tmpdir):
        """want：{键: 文件名谓词}。先找已解出的 .nc，缺的再从 sub/ 下的 tar 流式解到 tmpdir（每键取第一个匹配）。"""
        found = {}
        for p in self._glob(os.path.join(sub, "**", "*.nc")):
            for k, pred in want.items():
                if k not in found and pred(os.path.basename(p)):
                    found[k] = p
        if len(found) < len(want):
            for tp in self._glob(os.path.join(sub, "*.tar*")):
                with tarfile.open(tp, "r:*") as tf:
                    for m in tf:
                        base = os.path.basename(m.name)
                        for k, pred in want.items():
                            if k not in found and m.isfile() and pred(base):
                                found[k] = self._extract(tf, m, tmpdir)
        return found

    @staticmethod
    def _extract(tf, m, tmpdir):
        dst = os.path.join(tmpdir, os.path.basename(m.name))
        with tf.extractfile(m) as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo, 8 << 20)
        return dst

    def _all_nc(self, sub, tmpdir):
        """sub/ 下全部 .nc；没有则把 sub/ 下 tar 里的全部 .nc 成员解到 tmpdir。"""
        hits = self._glob(os.path.join(sub, "**", "*.nc"))
        if hits:
            return hits
        out = []
        for tp in self._glob(os.path.join(sub, "*.tar*")):
            with tarfile.open(tp, "r:*") as tf:
                for m in tf:
                    if m.isfile() and m.name.endswith(".nc") and not os.path.basename(m.name).startswith("."):
                        out.append(self._extract(tf, m, tmpdir))
        return sorted(out)

    # ---- Watson 2020（RECCAP2） ----
    def load_watson(self, year, tmpdir):
        """返回 {lat, lon, sf(12,nla,nlo), fa(12,…), source, sea_var, air_var, air_substituted, files}（L8/L33/L46）。
        source＝"direct"：sf/fa 是 fCO2（sfco2＋fco2atm_skinT 或 L8 替补），照旧直接用；
        source＝"pco2_fugacity_L46"：直接字段不齐，sf/fa 是公开 pCO2（spco2、pco2atm，µatm），由 Engine 逐日×FugFac(T)。"""
        import netCDF4
        cands = self._all_nc("watson", tmpdir)

        def find(name):
            for p in cands:
                with netCDF4.Dataset(p) as ds:
                    names = {k.lower(): k for k in ds.variables}
                    if name.lower() in names:
                        return p, names[name.lower()]
            return None

        sea = find(WATSON_SEA[0])
        air = find(WATSON_AIR)
        air_sub = False
        if air is None:
            for alt in WATSON_AIR_FALLBACK:
                air = find(alt)
                if air:
                    air_sub = True
                    break
        source = SRC_DIRECT
        if sea is None or air is None:                                  # L46：海侧缺（L8 未覆盖）或大气 fCO2 全缺 → 两侧都走 pCO2
            psea, pair = find(WATSON_SEA_PCO2), find(WATSON_AIR_PCO2)
            if psea is None or pair is None:
                raise StopAsk(f"L8/L46：Watson 直接 fCO2 不齐（sfco2{'缺' if sea is None else '有'}，"
                              f"大气 fCO2{'缺' if air is None else '有'}），公开 pCO2 也不齐（{WATSON_SEA_PCO2}"
                              f"{'缺' if psea is None else '有'}，{WATSON_AIR_PCO2}{'缺' if pair is None else '有'}）；"
                              f"按设计停下，需人工处理（候选文件 {len(cands)} 个）")
            sea, air, air_sub, source = psea, pair, False, SRC_L46
        out = {"source": source, "sea_var": sea[1], "air_var": air[1], "air_substituted": air_sub,
               "files": {"sea": sea[0], "air": air[0]}, "units": {}}
        for key, (p, name) in (("sf", sea), ("fa", air)):
            with netCDF4.Dataset(p) as ds:
                v = ds.variables[name]
                u = str(getattr(v, "units", "")).strip()
                out["units"][key] = u
                if source == SRC_L46 and u.lower() not in PCO2_UNITS_OK:
                    raise StopAsk(f"L46：{os.path.basename(p)} 的 {name} 单位为 '{u}'，不是 µatm；停下，需人工处理")
                tname, yms = _decode_ym(ds, v)
                if yms is None or any(x is None for x in yms):
                    n = ds.dimensions[tname].size if tname else 0
                    if n == 420:                                        # 1985-01…2019-12 月值（文件名所示）
                        yms = [(1985 + k // 12, k % 12 + 1) for k in range(n)]
                        out["time_assumed"] = "1985-01 起逐月（units 无法解码）"
                        tv = ds.variables.get(tname)
                        if tv is not None:                              # 旁证：公开包 time 无 units，值像「days since 1980-01-01」月中
                            vals = np.asarray(tv[:], float)
                            d80 = [date(1980, 1, 1) + timedelta(days=float(x)) for x in vals]
                            out.setdefault("time_check_days_since_1980", {})[key] = bool(
                                all((d.year, d.month) == ym for d, ym in zip(d80, yms)))
                    else:
                        raise DataError(f"{p} 时间轴无法解码（{n} 个）")
                arr = []
                for m in range(1, 13):
                    hit = [i for i, x in enumerate(yms) if x == (year, m)]
                    if not hit:
                        raise DataError(f"{p} 无 {year}-{m:02d}")
                    arr.append(_nc_take2d(ds, v, fixed={tname: hit[0]}))
                out[key] = np.stack(arr)
                if "lat" not in out:
                    out["lat"], out["lon"] = _nc_coord(ds, "lat"), _nc_coord(ds, "lon")
                else:
                    la, lo = _nc_coord(ds, "lat"), _nc_coord(ds, "lon")
                    if la.shape != out["lat"].shape or lo.shape != out["lon"].shape or not (
                            np.array_equal(la, out["lat"]) and np.array_equal(lo, out["lon"])):
                        raise DataError(f"Watson 海侧与大气侧网格不一致（{os.path.basename(p)}）")
        return out

    # ---- GLODAPv2.2016b mapped ----
    def load_glodap(self, tmpdir):
        import netCDF4
        want = {v: (lambda b, v=v: b == f"GLODAPv2.2016b.{v}.nc") for v in GLODAP_VARS}
        files = self._members("glodap", want, tmpdir)
        if len(files) < len(GLODAP_VARS):
            raise DataError(f"GLODAP 成员不全：{sorted(files)}")
        out = {"files": files}
        for var in GLODAP_VARS:
            with netCDF4.Dataset(files[var]) as ds:
                out[var] = _nc_take2d(ds, ds.variables[var])           # main L31–L34：.sel(depth_surface=0)
                if "lat" not in out:
                    out["lat"], out["lon"] = _nc_coord(ds, "lat"), _nc_coord(ds, "lon")
        return out

    # ---- WOA09 分区（可选） ----
    def load_basins(self):
        """woa09/ 下任一 *basin*.csv（逗号或空白分隔，前三列＝纬度、经度、表层分区码，1° 格）；没有则 None。"""
        # 只认 CSV（p3_fetch_extra 写 WOA09_basin_surface_1deg.csv）；原始 basin.msk 是 10F8.0 定宽文本，按行拆三列会得到垃圾分区
        hits = [p for p in self._glob(os.path.join("woa09", "**", "*basin*.csv")) if os.path.isfile(p)]
        if not hits:
            return None
        rows = []
        with open(hits[0], encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = re.split(r"[,\s]+", line.strip())
                try:
                    la, lo, code = float(parts[0]), float(parts[1]), float(parts[2])
                except (ValueError, IndexError):
                    continue
                rows.append((la, lo, code))
        if not rows:
            return None
        a = np.array(rows)
        lats = np.unique(a[:, 0])
        lons = np.unique(a[:, 1] % 360.0)
        grid = np.full((len(lats), len(lons)), np.nan)
        grid[np.searchsorted(lats, a[:, 0]), np.searchsorted(lons, a[:, 1] % 360.0)] = a[:, 2]
        return {"lat": lats, "lon": lons, "code": grid, "file": hits[0]}

    # ---- 清单 ----
    def manifest_records(self, rels):
        p = self.path("manifest.jsonl")
        want, out = set(rels), {}
        if not os.path.exists(p):
            return out
        with open(p, encoding="utf-8") as f:
            for line in f:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                if r.get("path") in want:
                    out[r["path"]] = {k: r.get(k) for k in ("bytes", "sha256", "status")}
        return out


class Fields:
    """把各源读成域像元向量（插值器按源网格缓存）。"""

    def __init__(self, lay, dom):
        self.lay, self.dom = lay, dom
        self.bil = {}
        self.diag = {"ice_as_fraction": None, "wind_persist_last": 0, "cmorph_missing": {}}
        self.wfiles = {}
        self.wcache = {}
        self.last_ice_frac = None

    def _bilin(self, key, lat, lon):
        k = (key, len(lat), len(lon), float(lat[0]), float(lat[-1]), float(lon[0]), float(lon[-1]))
        if k not in self.bil:
            tla, tlo = self.dom.tgt()
            self.bil[k] = Bilin(lat, lon, tla, tlo)
        return self.bil[k]

    def cmorph(self, dt):
        out, why = self.lay.read_cmorph(dt)
        if out is not None and out[0].shape != (len(self.dom.lat1d), len(self.dom.lon1d)):
            raise DataError(f"CMORPH {dt:%Y%m%d%H} 网格 {out[0].shape} 与静态网格不符")
        if out is None:
            self.diag["cmorph_missing"][dt.strftime("%Y%m%d%H")] = why
            nan = np.full(self.dom.N, np.nan)
            return nan, nan.copy()
        return self.dom.from_native(out[0]).astype(np.float64), self.dom.from_native(out[1]).astype(np.float64)

    def day(self, d):
        """(SST ℃, 冰 %, S0 psu) 于 UTC 日 d（L31、L9）。"""
        lat, lon, sst2d, ice2d = self.lay.read_oisst(d)
        b = self._bilin("oisst", lat, lon)
        sst = self.dom.gather(b(sst2d))
        icei = self.dom.gather(b(ice2d))
        mx = np.nanmax(ice2d) if np.isfinite(ice2d).any() else 0.0
        frac = bool(mx <= 1.0 + 1e-6)
        self.diag["ice_as_fraction"] = frac
        self.last_ice_frac = icei if frac else icei / 100.0             # 自测把同一分数场交给胶囊（main L73 自己 fillna×100）
        ice = np.where(np.isfinite(icei), icei, 0.0)                   # main L73：fillna(0)
        ice = ice * 100 if frac else ice
        hla, hlo, sss2d = self.lay.read_hycom(d)
        s0 = self.dom.gather(self._bilin("hycom", hla, hlo)(sss2d))
        return sst, ice, s0

    def _wind_hour(self, h):
        if h in self.wcache:
            return self.wcache[h]
        d = hour_dt(h)
        ym = d.year * 100 + d.month
        if ym not in self.wfiles:
            import netCDF4
            p = self.lay.era5_path(ym)
            if not os.path.exists(p):
                self.wfiles[ym] = None
            else:
                ds = netCDF4.Dataset(p)
                t = ds.variables["time"]
                dts = netCDF4.num2date(t[:], t.units, calendar=getattr(t, "calendar", "standard"),
                                       only_use_cftime_datetimes=False, only_use_python_datetimes=True)
                self.wfiles[ym] = (ds, {epoch_hour(x): i for i, x in enumerate(dts)})
        f = self.wfiles[ym]
        if f is None or h not in f[1]:
            return None
        ds, hmap = f
        w = np.ma.filled(np.ma.asarray(ds.variables["ws10"][hmap[h]]).astype(np.float64), np.nan)
        U = self.dom.gather(self._bilin("era5", np.asarray(ds.variables["latitude"][:], float),
                                        np.asarray(ds.variables["longitude"][:], float))(w))
        self.wcache[h] = U
        for k in [k for k in self.wcache if k < h - 1]:
            del self.wcache[k]
        return U

    def wind_pair(self, h):
        """(U(HH:00), U(HH:30))（L6/L32）。"""
        u0 = self._wind_hour(h)
        if u0 is None:
            raise DataError(f"ERA5 ws10 缺 {hour_dt(h):%Y-%m-%dT%H}Z")
        u1 = self._wind_hour(h + 1)
        if u1 is None:
            if h + 1 >= epoch_hour(datetime(YEAR + 1, 1, 1, tzinfo=UTC)):
                self.diag["wind_persist_last"] += 1
                return u0, u0
            raise DataError(f"ERA5 ws10 缺 {hour_dt(h + 1):%Y-%m-%dT%H}Z")
        return u0, (u0 + u1) / 2

    def close(self):
        for f in self.wfiles.values():
            if f:
                f[0].close()
        self.wfiles = {}


# ======================================================================== 静态场（主进程建一次，worker 以 mmap 共享）
STATIC_KEYS = ("lat1d", "lon1d", "row", "col", "area", "beta", "reg", "r_sa", "w_rows", "w_cols", "w_ok", "sf12", "fa12",
               "w_lat", "w_lon")


def build_static(lay, work_dir, tmpdir):
    sd = os.path.join(work_dir, "static")
    mp_ = os.path.join(sd, "meta.json")
    code_sha = sha256_file(os.path.abspath(__file__))
    if os.path.exists(mp_):
        with open(mp_, encoding="utf-8") as f:
            meta = json.load(f)
        if meta.get("version") == VERSION and meta.get("code_sha256") == code_sha and meta.get("cache_root") == lay.root:
            log(f"静态场沿用 {sd}")
            return sd
    os.makedirs(sd, exist_ok=True)
    t0 = time.perf_counter()
    lat1d, lon1d = lay.cmorph_grid(datetime(1999, 12, 31, 0, tzinfo=UTC))
    nla, nlo = len(lat1d), len(lon1d)
    gl = lay.load_glodap(tmpdir)
    slope_g = glodap_slopes(gl["TAlk"], gl["TCO2"], gl["salinity"], gl["temperature"])
    beta2d = Nearest(gl["lat"], gl["lon"], lat1d, lon1d)(slope_g)       # main L52–L60（float32 值原样保留）
    ocean = np.isfinite(beta2d)                                          # main L67（L10）
    rows, cols = np.nonzero(ocean)
    area2d = grid_cell_areas(lon1d, lat1d)
    band = np.clip(np.floor((lat1d[rows] - LAT_LO) / BAND_DEG).astype(int), 0, NBANDS - 1)
    bas = lay.load_basins()
    basin_codes = [0]
    if bas is not None:
        codes2d = Nearest(bas["lat"], bas["lon"], lat1d, lon1d)(bas["code"])
        cp = codes2d[rows, cols]
        cp = np.where(np.isfinite(cp), cp, -1).astype(int)
        basin_codes = sorted(set(cp.tolist()))
        bidx = np.searchsorted(np.array(basin_codes), cp)
    else:
        bidx = np.zeros(len(rows), int)
    nb = len(basin_codes)
    reg = band * nb + bidx
    flat = rows.astype(np.int64) * nlo + cols
    order = np.lexsort((flat, reg))
    rows, cols, reg = rows[order], cols[order], reg[order]
    gsw = _gsw()
    lonp, latp = lon1d[cols], lat1d[rows]
    r_sa = gsw.SA_from_SP(36.0, 0, lonp, latp) - gsw.SA_from_SP(35.0, 0, lonp, latp)
    wat = lay.load_watson(YEAR, tmpdir)
    wn = Nearest(wat["lat"], wat["lon"], lat1d, lon1d)
    arrs = {"lat1d": lat1d, "lon1d": lon1d, "row": rows.astype(np.int32), "col": cols.astype(np.int32),
            "area": area2d[rows, cols], "beta": beta2d[rows, cols].astype(np.float32), "reg": reg.astype(np.int32),
            "r_sa": np.asarray(r_sa, float), "w_rows": wn.rows[rows].astype(np.int32), "w_cols": wn.cols[cols].astype(np.int32),
            "w_ok": (wn.vr[rows] & wn.vc[cols]), "sf12": wat["sf"], "fa12": wat["fa"],
            "w_lat": np.asarray(wat["lat"], float), "w_lon": np.asarray(wat["lon"], float)}
    for k in STATIC_KEYS:
        np.save(os.path.join(sd, k + ".npy"), arrs[k])
    np.save(os.path.join(sd, "glodap_slope.npy"), slope_g)
    rels = [os.path.relpath(p, lay.root) for p in list(gl["files"].values()) + list(wat["files"].values())
            if p.startswith(lay.root)]
    meta = {"version": VERSION, "code_sha256": code_sha, "cache_root": lay.root, "nlat": nla, "nlon": nlo, "N": int(len(rows)), "nreg": int(NBANDS * nb), "nbasin": nb,
            "basin_codes": basin_codes, "basin_file": bas["file"] if bas else None, "band_deg": BAND_DEG,
            "cmorph_lat_range": [float(lat1d.min()), float(lat1d.max())],
            "ocean_frac_of_grid": float(ocean.mean()), "area_ocean_m2": float(arrs["area"].sum()),
            "glodap_files": {k: os.path.relpath(v, lay.root) if v.startswith(lay.root) else v for k, v in gl["files"].items()},
            "glodap_slope_cells": int(np.isfinite(slope_g).sum()),
            "watson": {"fco2_source": wat["source"], "sea_var": wat["sea_var"], "air_var": wat["air_var"],
                       "air_substituted": wat["air_substituted"], "units": wat["units"],
                       "files": {k: os.path.relpath(v, lay.root) if v.startswith(lay.root) else v for k, v in wat["files"].items()},
                       "time_assumed": wat.get("time_assumed"),
                       "time_check_days_since_1980": wat.get("time_check_days_since_1980")},
            "manifest": lay.manifest_records(rels), "seconds": round(time.perf_counter() - t0, 1)}
    jdump(meta, mp_)
    log(f"静态场：海洋像元 N={meta['N']:,}，分区 {meta['nreg']}，Watson {wat['source']}：海侧 {wat['sea_var']}、大气 {wat['air_var']}"
        f"{'（L8 替补）' if wat['air_substituted'] else ''}{'（L46 逐日 pCO2→fCO2）' if wat['source'] == SRC_L46 else ''}，"
        f"用时 {meta['seconds']} s")
    return sd


def load_static(sd):
    with open(os.path.join(sd, "meta.json"), encoding="utf-8") as f:
        meta = json.load(f)
    st = {k: np.load(os.path.join(sd, k + ".npy"), mmap_mode="r") for k in STATIC_KEYS}
    return st, meta


def domain_from_static(st, meta):
    return Domain(np.asarray(st["lat1d"]), np.asarray(st["lon1d"]), st["row"], st["col"], st["area"], st["beta"], st["reg"],
                  meta["nreg"], r_sa=st["r_sa"])


def watson_is_pco2(meta):
    """L46：静态场里的 sf12/fa12 是否为 pCO2（需 Engine 逐日换算）。"""
    return (meta.get("watson") or {}).get("fco2_source") == SRC_L46


def watson_month(st, m):
    """月 m 的 Watson 两侧值落到域像元（L33 最近邻，不改）；L46 下是 pCO2，换算在 Engine._mix。"""
    sf = np.asarray(st["sf12"][m - 1])[st["w_rows"], st["w_cols"]]
    fa = np.asarray(st["fa12"][m - 1])[st["w_rows"], st["w_cols"]]
    ok = np.asarray(st["w_ok"])
    return np.where(ok, sf, np.nan), np.where(ok, fa, np.nan)


# ======================================================================== 按月计算（worker）
def config_hash(cfg_core):
    return hashlib.sha256(json.dumps(cfg_core, sort_keys=True).encode()).hexdigest()[:16]


def wait_ready(lay, ym, wait_s, need_next=True):
    t_end = time.time() + wait_s
    last = None
    while True:
        ok, why = lay.ready(ym, need_next)
        if ok:
            return
        if time.time() > t_end:
            raise DataError(f"{ym} 输入未齐（等待超时）：{why}")
        if why != last:
            log(f"[{ym}] 等待输入：{why}")
            last = why
        time.sleep(POLL_S)


def process_month(cfg):
    """子进程入口：一个月（含前一日 48 个半步起转）；写 _work/month_YYYYMM.npz/json；返回摘要（异常也返回，不抛）。"""
    global _LOG_PATH
    ym = cfg["ym"]
    y, m = divmod(ym, 100)
    _LOG_PATH = os.path.join(cfg["work_dir"], f"log_{ym}.txt")
    base = os.path.join(cfg["work_dir"], f"month_{ym}")
    try:
        if cfg.get("resume") and os.path.exists(base + ".json"):
            with open(base + ".json", encoding="utf-8") as f:
                old = json.load(f)
            if (old.get("config_hash") == cfg["config_hash"] and old.get("days") == len(days_of_month(y, m))
                    and os.path.exists(base + ".npz")):                   # 不续用 bench 的 1 天结果
                log(f"[{ym}] 已有月结果，续用")
                return {"ym": ym, "ok": True, "resumed": True}
        lay = Layout(cfg["cache_root"])
        wait_ready(lay, ym, cfg["wait_s"], need_next=not cfg.get("max_days"))
        t_start = time.perf_counter()
        st, meta = load_static(cfg["static_dir"])
        dom = domain_from_static(st, meta)
        eng = Engine(dom, chunk=cfg["chunk"])
        fld = Fields(lay, dom)
        tim = {}

        def tk(key, t0):
            t1 = time.perf_counter()
            tim[key] = tim.get(key, 0.0) + t1 - t0
            return t1

        first = date(y, m, 1)
        prev = first - timedelta(days=1)
        t0 = time.perf_counter()
        for hh in range(24):                                            # 起转（L40）：前一日 48 个半步，只入缓冲
            dt = datetime(prev.year, prev.month, prev.day, hh, tzinfo=UTC)
            P0, P1 = fld.cmorph(dt)
            U0, U30 = fld.wind_pair(epoch_hour(dt))
            eng.push_only(P0, U0)
            eng.push_only(P1, U30)
        t0 = tk("spinup", t0)
        eng.set_month(*watson_month(st, m), pco2=watson_is_pco2(meta))
        days = days_of_month(y, m)
        if cfg.get("max_days"):
            days = days[:cfg["max_days"]]
        nstep = 0
        for d in days:
            t0 = time.perf_counter()
            sst, ice, s0 = fld.day(d)
            t0 = tk("day_read", t0)
            eng.set_day(sst, ice, s0)
            del sst, ice, s0
            t0 = tk("day_coefs", t0)
            for hh in range(24):
                dt = datetime(d.year, d.month, d.day, hh, tzinfo=UTC)
                P0, P1 = fld.cmorph(dt)
                t0 = tk("read_cmorph", t0)
                U0, U30 = fld.wind_pair(epoch_hour(dt))
                t0 = tk("wind", t0)
                for P, U in ((P0, U0), (P1, U30)):
                    le = cfg.get("l25_every") or 0
                    eng.l25_check = bool(le and nstep % le == 0)
                    eng.step(P, U)
                    nstep += 1
                t0 = tk("steps", t0)
            log(f"[{ym}] {d} 完成（累计 {nstep} 步，峰值 RSS {peak_rss_bytes() / 2**30:.2f} GiB）")
        fld.close()
        ec = np.asarray(eng.ring.entry_counts, float)
        info = {"ym": ym, "version": VERSION, "config_hash": cfg["config_hash"], "qnames": eng.qn, "nreg": dom.nreg,
                "days": len(days), "steps": nstep, "N": dom.N, "seconds": round(time.perf_counter() - t_start, 1),
                "timing_s": {**{k: round(v, 2) for k, v in tim.items()}, **{"eng_" + k: round(v, 2) for k, v in eng.tim.items()}},
                "peak_rss_bytes": peak_rss_bytes(),
                "entries_per_step": {"mean": float(ec.mean()) if len(ec) else 0, "max": float(ec.max()) if len(ec) else 0},
                "D_frac": {"mean": float(np.mean(eng.diag["D_frac"])) if eng.diag["D_frac"] else 0,
                           "max": float(np.max(eng.diag["D_frac"])) if eng.diag["D_frac"] else 0},
                "l25_max_rel": eng.diag["l25_max_rel"], "l25_n_checked": eng.diag["l25_n_checked"],
                "cmorph_missing": fld.diag["cmorph_missing"], "wind_persist_last": fld.diag["wind_persist_last"],
                "ice_as_fraction": fld.diag["ice_as_fraction"], "month_ready": lay.month_ready(ym),
                "watson_fco2_source": SRC_L46 if eng.w_pco2 else SRC_DIRECT, "watson_fugfac": eng.fug_summary()}
        if os.path.exists(base + ".json"):
            os.remove(base + ".json")                                    # json 最后写＝续跑判据；先删旧的，防 npz 写一半时被续用
        np.savez(base + ".npz", acc=eng.acc)
        jdump(info, base + ".json")
        log(f"[{ym}] 月结果已写（{info['seconds']} s）")
        return {"ym": ym, "ok": True, "resumed": False, "info_path": base + ".json"}
    except DataError as e:
        log(f"[{ym}] 数据错误：{e}")
        return {"ym": ym, "ok": False, "kind": "data", "error": str(e)}
    except StopAsk as e:
        log(f"[{ym}] 停下，需人工处理：{e}")
        return {"ym": ym, "ok": False, "kind": "stop", "error": str(e)}
    except Exception as e:
        log(f"[{ym}] 异常：{e}\n{traceback.format_exc()}")
        return {"ym": ym, "ok": False, "kind": "error", "error": f"{type(e).__name__}: {e}"}


# ======================================================================== 汇总 → p3a_repro.json、p3b_curves.json
def aggregate(work_dir, months, meta, out_dir, extra):
    qn = qnames()
    qi = {q: i for i, q in enumerate(qn)}
    nreg, nb = meta["nreg"], meta["nbasin"]
    Atot = np.zeros((len(qn), nreg))
    monthly = []
    infos = {}
    for ym in months:
        with open(os.path.join(work_dir, f"month_{ym}.json"), encoding="utf-8") as f:
            info = json.load(f)
        if info["qnames"] != qn:
            raise RuntimeError(f"{ym} 量名与当前版本不符")
        a = np.load(os.path.join(work_dir, f"month_{ym}.npz"))["acc"]
        Atot += a
        monthly.append((ym, a.sum(axis=1)))
        infos[ym] = info
    tot = Atot.sum(axis=1)
    bands = Atot.reshape(len(qn), NBANDS, nb).sum(axis=2)
    fug = None                                                          # L46：全年全场（像元×日）FugFac 统计
    if watson_is_pco2(meta):
        fug = {}
        for side in ("sea", "air"):
            ss = [infos[y_]["watson_fugfac"][side] for y_ in infos if infos[y_].get("watson_fugfac")]
            n = sum(x["n"] for x in ss)
            fug[side] = {"min": min(x["min"] for x in ss) if n else None, "max": max(x["max"] for x in ss) if n else None,
                         "mean": sum(x["sum"] for x in ss) / n if n else None, "n_pixel_days": n}
        fug["temperature"] = {"sea": "海侧溶解度所用 SST（CMORPH 像元、逐日）", "air": f"SST−{SKIN_DT} K（SKIN_DT，TB L498）"}
    T = lambda q: float(tot[qi[q]])
    W = T("W_all")
    pg = lambda x: x / 1e15
    pct = lambda x: 100.0 * x / W                      # 效应（ΣF_x−F_wind，负＝多吸收）÷ 风驱（负）→ 正＝多吸收（L14）
    skin = T("skin")
    I = lambda kind, lab: skin + T(f"{kind}|{lab}")      # 非 D 常数项＋D 上配对和（L30）
    dep_pct = -100.0 * T("dep") / W
    turb_pct = pct(T("turb"))
    dil_R, int_R = pct(I("dil", "whole_1")), pct(I("int", "whole_1"))
    dil_S, int_S = pct(I("dil", "s20_1")), pct(I("int", "s20_1"))
    U1 = dep_pct + int_R * (1 - CAPTURE_C)
    vals = {"wind_PgC": pg(W), "dilution_RIM_pct": dil_R, "interfacial_RIM_pct": int_R, "deposition_pct": dep_pct,
            "dilution_ratio_RIM_S20": dil_R / dil_S if dil_S else float("nan"), "turbulence_pct": turb_pct, "U1_pp": U1}
    gates = []
    for key, target, kind, tol in REPRO_GATES:
        v = vals[key]
        err = abs(v - target) / abs(target) if kind == "rel" else abs(v - target)
        gates.append({"item": key, "target": target, "value": v, "error_type": kind, "error": err, "tolerance": tol,
                      "pass": bool(np.isfinite(err) and err <= tol)})
    nW = T("n_W")

    def paired(kind, lab, base_pct):
        wp = W - T(f"Wm_{kind}|{lab}")
        return {"denominator_paired_PgC": pg(wp), "pct_paired_denominator": 100.0 * (base_pct / 100.0 * W) / wp if wp else None,
                "missing_frac": T(f"nm_{kind}|{lab}") / nW if nW else None}

    u_W = W
    unp = {
        "dilution_RIM_pct": 100.0 * (T("u_nonD_dil") + T("u_dil_D|whole_1") - u_W) / u_W,
        "interfacial_RIM_pct": 100.0 * (T("u_nonD_dil") + T("u_int_D|whole_1") - u_W) / u_W,
        "dilution_S20_pct": 100.0 * (T("u_nonD_dil") + T("u_dil_D|s20_1") - u_W) / u_W,
        "interfacial_S20_pct": 100.0 * (T("u_nonD_dil") + T("u_int_D|s20_1") - u_W) / u_W,
        "turbulence_pct": 100.0 * (T("u_nonD_W") + T("u_turb_D") - u_W) / u_W,
        "deposition_pct": -100.0 * T("u_dep") / u_W}
    int_z5 = pct(I("int", "z5_1"))
    int_0 = pct(I("int", "whole_0"))
    repro = {
        "version": VERSION, "code_sha256": extra["code_sha256"], "written_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "reads_p2": False, "note_L23": "本文件写盘前、写盘时均未读取任何 P2 产物；P3c 须在本文件存在后才运行。",
        "months": months, "all_pass": all(g["pass"] for g in gates), "gates": gates,
        "table_pct_of_wind": {"wind_PgC": pg(W), "interfacial_RIM": int_R, "dilution_RIM": dil_R, "turbulence": turb_pct,
                              "deposition": dep_pct, "total_RIM": dep_pct + int_R, "interfacial_S20": int_S,
                              "dilution_S20": dil_S, "total_S20": dep_pct + int_S, "U1_reconstructed": U1},
        "table_PgC_uptake_positive": {"interfacial_RIM": -pg(I("int", "whole_1")), "dilution_RIM": -pg(I("dil", "whole_1")),
                                      "turbulence": -pg(T("turb")), "deposition": pg(T("dep")),
                                      "interfacial_S20": -pg(I("int", "s20_1")), "dilution_S20": -pg(I("dil", "s20_1"))},
        "witte_fig1": {"wind_PgC": -1.468, "interfacial_RIM": 9.67, "dilution_RIM": 7.97, "turbulence": 0.14, "deposition": 3.88,
                       "total_RIM": 13.56, "interfacial_S20": 4.50, "dilution_S20": 3.88, "total_S20": 8.38, "U1": 10.65},
        "non_gate": {
            "z5_capture_of_interfacial": int_z5 / int_R if int_R else None, "witte_z5_capture": WITTE_Z5_CAPTURE,
            "z5_capture_excluding_s0_term": (int_z5 - int_0) / (int_R - int_0) if int_R != int_0 else None,
            "skin_term_nonD_pct": pct(skin),
            "interfacial_at_s0_whole_pct": int_0,
            "dilution_at_s0_whole_pct": pct(I("dil", "whole_0")),
            "note_L30": "胶囊稀释情形大气侧不加 0.1 psu 皮层盐度（TB L498 vs L501/L504），无雨像元也产生与 s 无关的界面项；"
                        "dilution_at_s0_whole_pct 即该项全局量。"},
        "diagnostics": {
            "paired_denominator": {"dilution_RIM": paired("dil", "whole_1", dil_R), "interfacial_RIM": paired("int", "whole_1", int_R),
                                   "dilution_S20": paired("dil", "s20_1", dil_S), "interfacial_S20": paired("int", "s20_1", int_S)},
            "witte_style_unpaired_pct": unp,
            "coverage": {"pixel_steps_ocean": T("n_ocean"), "pixel_steps_wind_valid": nW,
                         "wind_valid_frac": nW / T("n_ocean") if T("n_ocean") else None,
                         "rain_history_D_frac": T("n_D") / nW if nW else None, "current_rain_frac": T("n_C") / nW if nW else None,
                         "precip_nan_frac": T("n_Pnan") / T("n_ocean") if T("n_ocean") else None,
                         "precip_gt200_count": T("n_Pgt200"), "u10_zero_count": T("n_U0"),
                         "u10_nan_frac": T("n_Unan") / T("n_ocean") if T("n_ocean") else None,
                         "precip_u10_both_zero_count": T("n_PU0"),
                         "rim_missing_frac_of_wind_valid": T("n_poisonR") / nW if nW else None,
                         "dep_missing_frac": T("nm_dep") / nW if nW else None, "turb_missing_frac": T("nm_turb") / nW if nW else None},
            "l25_fallback_frac": {lab: (T(f"nfb|{lab}") / T("n_D") if T("n_D") else None) for lab in LABELS},
            "negative_salinity_clipped": {lab: T(f"nclip|{lab}") for lab in LABELS if lab.startswith("whole_")},
            "watson_fco2_source": (meta.get("watson") or {}).get("fco2_source", SRC_DIRECT),
            "watson_fugfac": fug,
            "cmorph_missing_halfhours": sum(len(infos[y_]["cmorph_missing"]) * 2 for y_ in infos),
            "wind_persist_last": sum(infos[y_]["wind_persist_last"] for y_ in infos),
            "l25_max_rel_bench_checks": {y_: infos[y_]["l25_max_rel"] for y_ in infos if infos[y_]["l25_max_rel"]}},
        "inputs": extra.get("inputs", {}), "static": {k: meta[k] for k in ("N", "nreg", "nbasin", "basin_file", "watson",
                                                                           "glodap_files", "area_ocean_m2")},
    }
    jdump(repro, os.path.join(out_dir, "p3a_repro.json"))                   # 先写（L23 顺序）
    log(f"p3a_repro.json 已写：复现门 {'全过' if repro['all_pass'] else '未全过'}")

    curves = {"version": VERSION, "reads_p2": False, "s_grid": list(S_GRID), "c": CAPTURE_C, "wind_PgC": pg(W),
              "deposition_pct": dep_pct, "turbulence_pct": turb_pct, "skin_term_nonD_pct": pct(skin),
              "definitions": {"I_int_pct": "界面项（稀释＋湍流＋交互＋L30 常数项）占 |风驱| 百分比",
                              "U_pp": "Dep＋I_int·(1−c)", "Delta_pp": "U(s)−U(1)",
                              "G_int_pct": "I_int 对 β 线性的分量；S6：I(s; r_β)＝I(s)＋(r_β−1)·G(s)（L39，精确）",
                              "bands": "5° 纬带 PgC（Fwind、turb、skin、I、G 正＝出海；dep_PgC 为 ΣFdep，正＝入海；I、G 为配对和，"
                                       "含该带 L30 常数项）；P3c 自行按全局风驱换算百分比"},
              "global": {}, "bands": {"lat_lo": [LAT_LO + BAND_DEG * b for b in range(NBANDS)],
                                      "lat_hi": [LAT_LO + BAND_DEG * (b + 1) for b in range(NBANDS)],
                                      "Fwind_PgC": pg(bands[qi["W_all"]]), "dep_PgC": pg(bands[qi["dep"]]),
                                      "turb_PgC": pg(bands[qi["turb"]]), "skin_PgC": pg(bands[qi["skin"]])}}
    for kind in ("whole", "hist"):
        g = {"s": list(S_GRID), "I_int_pct": [], "I_dil_pct": [], "U_pp": [], "Delta_pp": [], "G_int_pct": [], "G_dil_pct": []}
        bI, bG = [], []
        for s in S_GRID:
            lab = f"{kind}_{s:g}"
            ii = pct(I("int", lab))
            g["I_int_pct"].append(ii)
            g["I_dil_pct"].append(pct(I("dil", lab)))
            g["U_pp"].append(dep_pct + ii * (1 - CAPTURE_C))
            g["G_int_pct"].append(pct(T(f"gi|{lab}")))
            g["G_dil_pct"].append(pct(T(f"gd|{lab}")))
            bI.append(pg(bands[qi["skin"]] + bands[qi[f"int|{lab}"]]))
            bG.append(pg(bands[qi[f"gi|{lab}"]]))
        u1 = g["U_pp"][S_GRID.index(1.0)]
        g["Delta_pp"] = [u - u1 for u in g["U_pp"]]
        du = np.diff(g["U_pp"])
        g["monotone"] = bool(np.all(du >= 0) or np.all(du <= 0))
        line = [P3_0_LINE[0] + P3_0_LINE[1] * s for s in S_GRID]
        g["max_abs_dev_from_P3_0_line_pp"] = float(np.max(np.abs(np.array(g["U_pp"]) - line)))
        u0 = g["U_pp"][0]
        sg = np.array(S_GRID)
        dev_own = np.abs(np.array(g["U_pp"]) - (u0 + (u1 - u0) * sg))
        g["max_abs_dev_from_own_endpoints_line_pp"] = float(dev_own[sg <= 1.0].max())
        curves["global"][kind] = g
        curves["bands"][f"{kind}_I_int_PgC"] = [list(x) for x in np.array(bI).T]
        curves["bands"][f"{kind}_G_int_PgC"] = [list(x) for x in np.array(bG).T]
    for lab in ("s20_1", "z5_1"):
        ii = pct(I("int", lab))
        curves["global"][lab] = {"I_int_pct": ii, "I_dil_pct": pct(I("dil", lab)), "U_pp": dep_pct + ii * (1 - CAPTURE_C),
                                 "G_int_pct": pct(T(f"gi|{lab}"))}
        curves["bands"][f"{lab}_I_int_PgC"] = pg(bands[qi["skin"]] + bands[qi[f"int|{lab}"]])
        curves["bands"][f"{lab}_G_int_PgC"] = pg(bands[qi[f"gi|{lab}"]])
    curves["checks_non_gate"] = {"note": "非门检查：P3b 与 P3-0 直线最大偏差（>0.5 pp 表示非线性）、U(s) 单调性；"
                                         "own_endpoints 只取 s≤1 段，扣除 L30 常数项与复现偏差后看非线性"}
    jdump(curves, os.path.join(out_dir, "p3b_curves.json"))

    with open(os.path.join(out_dir, "p3_bands.csv"), "w", encoding="utf-8") as f:
        cols = ["W_all", "dep", "turb", "skin"] + [f"int|{lab}" for lab in LABELS] + [f"gi|{lab}" for lab in LABELS]
        f.write("lat_lo,lat_hi," + ",".join(c + "_PgC" for c in cols) + "\n")
        for b in range(NBANDS):
            f.write(f"{LAT_LO + BAND_DEG * b:g},{LAT_LO + BAND_DEG * (b + 1):g}," +
                    ",".join(f"{pg(bands[qi[c], b]):.9g}" for c in cols) + "\n")
    with open(os.path.join(out_dir, "p3_monthly.csv"), "w", encoding="utf-8") as f:
        cols = ["W_all", "dep", "turb", "skin", "int|whole_1", "dil|whole_1", "int|s20_1", "dil|s20_1", "n_W", "n_D"]
        f.write("month," + ",".join(cols) + "\n")
        for ym, v in monthly:
            f.write(f"{ym}," + ",".join(f"{(pg(v[qi[c]]) if not c.startswith('n_') else v[qi[c]]):.9g}" for c in cols) + "\n")
    if nb > 1:
        with open(os.path.join(out_dir, "p3_regions.csv"), "w", encoding="utf-8") as f:
            cols = ["W_all", "dep", "turb", "skin", "int|whole_1", "int|s20_1"]
            f.write("lat_lo,basin_code," + ",".join(c + "_PgC" for c in cols) + "\n")
            for r in range(nreg):
                if Atot[qi["n_ocean"], r] > 0:
                    f.write(f"{LAT_LO + BAND_DEG * (r // nb):g},{meta['basin_codes'][r % nb]}," +
                            ",".join(f"{pg(Atot[qi[c], r]):.9g}" for c in cols) + "\n")
    return repro, curves, infos


# ======================================================================== 代码等价性自测（L37）
def cmp_fields(ours, ref, floor, tol, exclude=None):
    o = np.asarray(ours, float)
    r = np.asarray(ref, float)
    m = np.ones(o.shape, bool) if exclude is None else ~exclude
    mism = (np.isnan(o) ^ np.isnan(r)) & m
    fin = np.isfinite(o) & np.isfinite(r) & m
    rel = np.abs(o[fin] - r[fin]) / np.maximum(np.abs(r[fin]), floor)
    mx = float(rel.max()) if rel.size else 0.0
    return {"n_compared": int(fin.sum()), "nan_mismatch": int(mism.sum()), "max_rel": mx, "tol": tol,
            "excluded": int((~m).sum()), "pass": bool(mism.sum() == 0 and mx <= tol and fin.sum() > 0)}


def synth_series(rng, M=400, L=110):
    """RIM 合成序列：普通雨、强雨、以及胶囊算术的边界（P=U=0→NaN、U=0 有雨→0、P>200→NaN、P 为 NaN、U 极小）。"""
    U = np.clip(0.3 + 14 * rng.random((1, M)) + rng.normal(0, 0.8, (L, M)), 0.3, None)
    P = np.where(rng.random((L, M)) < 0.15, rng.gamma(0.8, 6.0, (L, M)), 0.0)
    P[:, rng.random(M) < 0.3] = 0.0
    P[40:52, 5] = 120.0
    P[60, 0], U[60, 0] = 0.0, 0.0
    P[62, 1], U[62, 1] = 5.0, 0.0
    P[70, 2] = 250.0
    P[75, 3] = np.nan
    U[:, 4] = 1e-3
    P[50:60, 4] = 3.0
    U[80, 6] = 250.0
    P[80, 6] = 1.0
    return P, U


def test_p2_equiv(P, U, tag):
    """z=0 与 z=5 m：RimRing（本实现）对 p2_rim_test.rim_factor（P2 参考实现，L4）逐值比较。"""
    sys.path.insert(0, SRC_DIR)
    import p2_rim_test as p2
    L, M = P.shape
    ring = RimRing(M)
    o0, o5 = [], []
    for t in range(L):
        ent = ring.entries(P[t], U[t])
        if t >= 48:
            logH, logH5, _ = ring.history()
            c0 = np.ones(M)
            c5 = np.ones(M)
            c0[ent["pos"]] = ent["cur0"]
            c5[ent["pos"]] = ent["cur5"]
            with np.errstate(all="ignore"):
                o0.append(np.exp(logH) * c0)
                o5.append(np.exp(logH5) * c5)
        ring.push(ent)
    o0, o5 = np.array(o0), np.array(o5)
    with np.errstate(all="ignore"):
        r0 = np.stack([p2.rim_factor(P[:, j], U[:, j], 0.0) for j in range(M)], axis=1)
        r5 = np.stack([p2.rim_factor(P[:, j], U[:, j], Z5) for j in range(M)], axis=1)
    res = {"ran": True, "tag": tag, "pixels": M, "outputs_per_pixel": L - 48,
           "z0": cmp_fields(o0, r0, 1e-300, TOL_P2), "z5": cmp_fields(o5, r5, 1e-300, TOL_P2),
           "n_nan_ref_z0": int(np.isnan(r0).sum()), "n_zero_ref_z0": int((r0 == 0).sum())}
    res["pass"] = res["z0"]["pass"] and res["z5"]["pass"]
    return res


def test_l25_grid():
    """L25：展开（含 L27 兜底）对精确 gsw，覆盖 s 网格全部取值、两种缩放方式、z5 型因子，SST −1.8–32 ℃、S0 5–40。"""
    Ts = np.array([-1.8, 0, 5, 10, 15, 20, 25, 28, 32], float)
    S0s = np.array([5, 10, 15, 20, 25, 30, 33, 35, 37, 40], float)
    Fs = np.array([1 - 1e-7, 0.9999, 0.999, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.6, 0.3, 0.05], float)
    locs = [(0.0, 0.0), (150.0, 5.0), (330.0, 40.0), (20.0, -55.0), (20.0, 58.0), (270.0, 25.0)]
    TT, SS, FF, LL = np.meshgrid(Ts, S0s, Fs, np.arange(len(locs)), indexing="ij")
    T, S0, F = TT.ravel(), SS.ravel(), FF.ravel()
    lon = np.array([locs[k][0] for k in LL.ravel()])
    lat = np.array([locs[k][1] for k in LL.ravel()])
    gsw = _gsw()
    r = gsw.SA_from_SP(36.0, 0, lon, lat) - gsw.SA_from_SP(35.0, 0, lon, lat)
    co = day_coefs(T, S0, lon, lat, r)
    logF = np.log(F)
    S_R = S0 * F
    worst, worst_nofb, nfb, n = 0.0, 0.0, 0, 0
    for kind in ("whole", "hist"):
        for s in S_GRID:
            if s == 1.0 or (kind == "whole" and s == 0.0):
                continue
            if kind == "whole":
                S = np.maximum(S0 + s * (S_R - S0), 0.0)
            else:
                S = S0 * (np.exp(s * 0.7 * logF) * np.exp(0.3 * logF))
            aw, aa, fb = alpha_expand(S, S0, co, T, lon, lat)
            ew, ea = sol_pair(T, S, lon, lat)
            with np.errstate(all="ignore"):
                rel = np.maximum(np.abs(aw / ew - 1), np.abs(aa / ea - 1))
            ok = np.isfinite(rel)
            worst = max(worst, float(rel[ok].max()))
            if (ok & ~fb).any():
                worst_nofb = max(worst_nofb, float(rel[ok & ~fb].max()))
            nfb += int(fb.sum())
            n += int(ok.sum())
    return {"ran": True, "n": n, "fallback_n": nfb, "max_rel_all": worst, "max_rel_expansion_only": worst_nofb, "tol": L25_TOL,
            "thresholds": {"dmax_psu": L25_DMAX, "smin_psu": L25_SMIN}, "pass": bool(worst < L25_TOL)}


def test_fugacity():
    """L46：fug_factor 对 PyCO2SYS gas.fugacity_factor（1.8.x；opt_k_carbonic=10 即非 GEOSECS/Peng，pressure_bar=0 即 1 atm，
    R＝constants.RGasConstant_CODATA2018）在 T＝−2、0、10、20、30 ℃ 相对差 ≤1e-9。"""
    T = np.array(FUG_TEST_T, float)
    ours = fug_factor(T)
    res = {"ran": True, "T_C": list(FUG_TEST_T), "ours": ours.tolist(), "tol": FUG_TOL}
    try:
        import PyCO2SYS as pyco2
        from PyCO2SYS import constants as pc, gas as pg_
        R = float(pc.RGasConstant_CODATA2018)
        ref = np.asarray(pg_.fugacity_factor(T, 10, R, 0.0), float)
        res.update({"reference": f"PyCO2SYS {pyco2.__version__} gas.fugacity_factor(T, 10, RGasConstant_CODATA2018, 0.0)",
                    "pyco2sys_R": R, "R_equal": bool(R == FUG_R), "pyco2sys_Tzero": float(getattr(pc, "Tzero", float("nan"))),
                    "ref": ref.tolist()})
    except Exception as e:
        res.update({"pass": False, "error": f"PyCO2SYS 逸度系数函数不可用：{type(e).__name__}: {e}"})
        return res
    rel = np.abs(ours / ref - 1)
    res["max_rel"] = float(rel.max())
    res["pass"] = bool(np.all(np.isfinite(rel)) and rel.max() <= FUG_TOL)
    return res


def synth_block(rng, n=50, L=50):
    lat = 5.0 + 0.072756 * np.arange(n)
    lon = 150.0 + 0.072756 * np.arange(n)
    times = (np.datetime64("1999-12-31T00:00:00") + np.arange(L) * np.timedelta64(1800, "s")).astype("datetime64[s]").astype(np.int64)
    base = 1.0 + 12.0 * rng.random((n, n, 1))
    U = np.clip(base + rng.normal(0, 0.7, (n, n, L)), 0.4, None)
    P = np.where(rng.random((n, n, L)) < 0.2, rng.gamma(0.8, 5.0, (n, n, L)), 0.0)
    P[rng.random((n, n)) < 0.5] = 0.0
    P[10:15, 10:15, 38:50] = 60.0
    U[10:13, 10:13, :] = 0.8
    P[20, 20, 30] = 250.0
    icefrac = np.full((n, n), np.nan)
    icefrac[0:3, 0:3] = 0.4
    slope = (8.0 + 10.0 * rng.random((n, n))).astype(np.float32)
    slope[45:, 45:] = np.nan
    return {"lat": lat, "lon": lon, "times": times, "P": P, "U": U, "S0": 33.0 + 3.0 * rng.random((n, n)),
            "sst": 2.0 + 28.0 * rng.random((n, n)), "icefrac": icefrac, "sf": 350.0 + 100.0 * rng.random((n, n)),
            "fa": 360.0 + 10.0 * rng.random((n, n)), "slope": slope}


def ours_on_block(blk):
    """本实现（Engine，与全场同一代码路径）在 50×50 块上跑 48 个起转半步＋2 个通量半步，返回逐像元场（la, lo, 2）。"""
    slope = blk["slope"]
    rows, cols = np.nonzero(np.isfinite(slope))
    dom = Domain(blk["lat"], blk["lon"], rows, cols, np.ones(len(rows)), slope[rows, cols], np.zeros(len(rows)), 1)
    eng = Engine(dom, keep=True)
    eng.l25_check = True
    for t in range(48):
        eng.push_only(blk["P"][rows, cols, t], blk["U"][rows, cols, t])
    if "sf_p" in blk:                                                    # L46：引擎吃 pCO2 自己换算；胶囊吃测试侧换算好的 blk["sf"/"fa"]
        eng.set_month(blk["sf_p"][rows, cols], blk["fa_p"][rows, cols], pco2=True)
    else:
        eng.set_month(blk["sf"][rows, cols], blk["fa"][rows, cols])
    fr = blk["icefrac"][rows, cols]
    eng.set_day(blk["sst"][rows, cols], np.where(np.isfinite(fr), fr, 0.0) * 100, blk["S0"][rows, cols])
    hnan = np.zeros(slope.shape + (2,), bool)
    for j, t in enumerate((48, 49)):
        # L47：24 h 历史含 NaN 型条目＝RimRing.history() 的 logH 为 NaN——与 Engine.step 随后所用同一调用（只读，不改环）
        hnan[rows, cols, j] = np.isnan(eng.ring.history()[0])
        eng.step(blk["P"][rows, cols, t], blk["U"][rows, cols, t])
    out = {}
    for k in KEEP_FIELDS:
        a = np.full(slope.shape + (2,), np.nan)
        for j in (0, 1):
            a[rows, cols, j] = eng.kept[j][k]
        out[k] = a
    return out, eng.diag["l25_max_rel"], eng, hnan


RIM_FIELDS_L47 = ("salR", "Fdil_RIMv3", "Fint_RIMv3")


def rim_excluded_report(ours, ref, hnan, nmax=20):
    """L47：被排除像元×半步（历史含 NaN 型条目）的清单与正向断言——本实现 RIM 三场全 NaN、胶囊 f64 三场全有限。"""
    ours_nan = {k: bool(np.all(np.isnan(ours[k][hnan]))) for k in RIM_FIELDS_L47}
    cap_fin = {k: bool(np.all(np.isfinite(ref[f"f64_{k}"][hnan]))) for k in RIM_FIELDS_L47}
    idx = np.argwhere(hnan)
    return {"count": int(hnan.sum()), "fields": list(RIM_FIELDS_L47),
            "criterion": "np.isnan(logH)，logH＝RimRing.history()[0]（Engine.step 同一调用）",
            "coords_row_col_halfstep": [[int(r), int(c), 48 + int(j)] for r, c, j in idx[:nmax]], "coords_listed_max": nmax,
            "positive_assertion": {"ours_all_nan": ours_nan, "capsule_f64_all_finite": cap_fin,
                                   "pass": bool(all(ours_nan.values()) and all(cap_fin.values()))}}


def check_accum(eng, ref, pnan, rim_excl):
    """累计管线对胶囊逐像元场之和（面积≡1）：风驱、沉降、湍流、RIM/S20 稀释与界面（非 D 常数项＋D 配对和）。
    逐像元自测不覆盖 reduceat 行序、skin 拼接、配对掩码，这里补上（审查修订）。rim_excl：L47 掩码，从胶囊 RIM 求和里扣除。"""
    tot = eng.acc.sum(axis=1)
    T = lambda q: float(tot[eng.qi[q]])
    g = lambda k: ref[f"f64_{k}"]
    Fw = g("Fwind")
    okW = np.isfinite(Fw)
    res, ok_all = {}, True

    def one(name, ours, terms):
        nonlocal ok_all
        refv = float(terms.sum())
        scale = float(np.abs(terms).sum()) + 1e-30
        err = abs(ours - refv) / scale
        passed = bool(np.isfinite(ours) and err <= 1e-9)
        ok_all &= passed
        res[name] = {"ours": ours, "capsule_sum": refv, "err_rel_to_abs_sum": err, "pass": passed}

    one("W_all", T("W_all"), Fw[okW])
    dep = g("Fdep")
    m = okW & np.isfinite(dep)
    one("dep", T("dep"), dep[m])
    for name, key, q in (("turb", "Fturb", "turb"), ("dil_RIM", "Fdil_RIMv3", "dil|whole_1"), ("int_RIM", "Fint_RIMv3", "int|whole_1"),
                         ("dil_S20", "Fdil_S20", "dil|s20_1"), ("int_S20", "Fint_S20", "int|s20_1")):
        F = g(key)
        m = okW & np.isfinite(F)
        if key.endswith("S20"):
            m &= ~pnan                                                   # L35：本实现 P 缺测记缺测，胶囊置 0
        elif key in RIM_FIELDS_L47:
            m &= ~rim_excl                                               # L47：历史含 NaN 型条目，本实现缺测，胶囊 nanprod 计入
        ours = T(q) if q == "turb" else T("skin") + T(q)
        one(name, ours, (F - Fw)[m])
    res["pass"] = ok_all
    return res


def run_capsule_ref(args):
    """（子进程）用胶囊原函数算参考值：main.py L23–L24、L67–L87 的流程，输入来自 npz。"""
    import xarray as xr
    sys.path.insert(0, args.capsule_dir)
    import CO2_Rain_Flux_Toolbox as tb
    d = np.load(args.ref_in)
    out = {}
    lat, lon = d["lat"], d["lon"]
    tt = d["times"].astype("datetime64[s]").astype("datetime64[ns]")
    dims3 = ("latitude", "longitude", "time")
    for tag, dt in (("f64", np.float64), ("f32", np.float32)):
        if tag == "f32" and not int(d["with_f32"]):
            continue
        co3 = {"latitude": lat, "longitude": lon, "time": tt}
        co2 = {"latitude": lat, "longitude": lon, "time": tt[48:50]}
        rep2 = lambda a: np.repeat(np.asarray(a)[:, :, None], 2, axis=2).astype(dt)
        precip = xr.DataArray(d["P"].astype(dt), coords=co3, dims=dims3)
        wind = xr.DataArray(d["U"].astype(dt), coords=co3, dims=dims3)
        sal = xr.DataArray(np.repeat(d["S0"][:, :, None], len(tt), axis=2).astype(dt), coords=co3, dims=dims3)
        sal_RIMv3 = tb.RIMv3(precip, wind, sal)                                                       # main L23
        sal_S20 = tb.S20(precip.isel(time=slice(48, 50)), wind.isel(time=slice(48, 50)), sal.isel(time=slice(48, 50)))  # L24
        slope = xr.DataArray(d["slope"], coords={"latitude": lat, "longitude": lon}, dims=("latitude", "longitude"))
        landmask = np.isnan(slope)                                                                    # L67
        sst = xr.DataArray(rep2(d["sst"]), coords=co2, dims=dims3).where(~landmask)                  # L71–L72
        icec = xr.DataArray(rep2(d["icefrac"]), coords=co2, dims=dims3).where(~landmask)
        ice = icec.fillna(0) * 100                                                                    # L73
        precip_h = precip.isel(time=slice(48, 50)).where(~landmask)                                   # L74
        u10 = wind.isel(time=slice(48, 50)).where(~landmask)                                          # L75
        salinity = xr.Dataset({"so": sal.isel(time=slice(48, 50))}).transpose("latitude", "longitude", "time").where(~landmask)  # L77
        salinity["sal_diluted_RIMv3"] = sal_RIMv3.where(~landmask)                                    # L78
        salinity["sal_diluted_S20"] = sal_S20.where(sal_S20 > 0, other=0).where(~landmask)            # L79
        sal_change_RIMv3 = salinity.sal_diluted_RIMv3 - salinity.so                                   # L80
        sal_change_S20 = salinity.sal_diluted_S20 - salinity.so                                       # L81
        fCO2 = xr.Dataset({"sfco2": xr.DataArray(rep2(d["sf"]), coords=co2, dims=dims3),
                           "fco2atm_skinT": xr.DataArray(rep2(d["fa"]), coords=co2, dims=dims3)}).where(~landmask)  # L83
        fCO2["sfco2_diluted_RIMv3"] = fCO2.sfco2 + slope * sal_change_RIMv3                           # L84
        fCO2["sfco2_diluted_S20"] = fCO2.sfco2 + slope * sal_change_S20                               # L85
        fl = tb.calculate_allFluxes_RIMv3_S20(fCO2, sst, salinity, u10, ice, precip_h)                # L87
        tr = lambda x: np.asarray(x.transpose("latitude", "longitude", "time").values, np.float64)
        for name in ("Fwind", "Fdep", "Fturb", "Fdil_RIMv3", "Fdil_S20", "Fint_RIMv3", "Fint_S20"):
            out[f"{tag}_{name}"] = tr(fl[name])
        out[f"{tag}_salR"] = tr(salinity.sal_diluted_RIMv3)
        out[f"{tag}_salS20"] = tr(salinity.sal_diluted_S20)
    if "g_ta" in d.files:
        out["slopes"] = np.array([tb.calculate_fco2_dilution_per_psu(ta=float(a), tco2=float(b), sst=float(c), sss=float(e))
                                  for a, b, c, e in zip(d["g_ta"], d["g_tco2"], d["g_sst"], d["g_sss"])])
    if "area_lon" in d.files:
        out["area"] = tb.grid_cell_areas(d["area_lon"], d["area_lat"])
    np.savez(args.ref_out, **out)
    return 0


def find_capsule(args, lay, out_dir):
    """胶囊两文件：--capsule-dir → src/witte_capsule → <cache>/witte_capsule → 公开 API 取到 <out>/_capsule；sha256 必须一致（L37）。"""
    cands = [args.capsule_dir] if args.capsule_dir else []
    cands += [os.path.join(SRC_DIR, "witte_capsule")] + ([lay.path("witte_capsule")] if lay else [])
    for d in cands:
        if d and all(os.path.exists(os.path.join(d, n)) for n in CAPSULE_FILES):
            bad = {n: sha256_file(os.path.join(d, n)) for n in CAPSULE_FILES if sha256_file(os.path.join(d, n)) != CAPSULE_FILES[n]}
            if bad:
                return None, f"{d} sha256 不符：{bad}"
            return d, "local"
    d = os.path.join(out_dir, "_capsule")
    os.makedirs(d, exist_ok=True)
    for n, sha in CAPSULE_FILES.items():
        req = urllib.request.Request(CAPSULE_API.format(name=n), headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
        except Exception as e:
            return None, f"取胶囊 {n} 失败：{e}"
        if hashlib.sha256(body).hexdigest() != sha:
            return None, f"胶囊 {n} sha256 不符（API 返回 {len(body)} B）"
        with open(os.path.join(d, n), "wb") as f:
            f.write(body)
    return d, "fetched"


def capsule_call(args, cap_dir, blk, extra, out_dir, tag):
    inp = os.path.join(out_dir, "_selftest", f"{tag}_in.npz")
    outp = os.path.join(out_dir, "_selftest", f"{tag}_ref.npz")
    os.makedirs(os.path.dirname(inp), exist_ok=True)
    np.savez(inp, with_f32=1, **{k: v for k, v in blk.items()}, **extra)
    py = args.capsule_python or sys.executable
    cmd = [py, os.path.abspath(__file__), "--mode", "capsule-ref", "--ref-in", inp, "--ref-out", outp, "--capsule-dir", cap_dir]
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0 or not os.path.exists(outp):
        return None, {"returncode": r.returncode, "stderr_tail": r.stderr[-3000:], "seconds": round(time.perf_counter() - t0, 1)}
    return dict(np.load(outp)), {"returncode": 0, "seconds": round(time.perf_counter() - t0, 1), "python": py}


def test_capsule_block(args, cap_dir, blk, out_dir, tag, extra=None):
    res = {"ran": True, "tag": tag, "block_shape": list(blk["slope"].shape)}
    ours, l25, eng, hnan = ours_on_block(blk)
    ref, meta = capsule_call(args, cap_dir, blk, extra or {}, out_dir, tag)
    res["capsule_run"] = meta
    res["l25_max_rel_on_block"] = l25
    res["l25_pass"] = bool(all(v < L25_TOL for v in l25.values()))
    if ref is None:
        res["rim_excluded_nan_history"] = {"count": int(hnan.sum()), "positive_assertion": {"pass": False, "why": "胶囊未跑成"}}
        res["pass"] = False
        return res, None
    pnan = np.isnan(blk["P"][:, :, 48:50])                      # L35：S20 在 P 缺测时本实现记缺测，胶囊置 0 → 不比
    fields = {}
    for k in KEEP_FIELDS:
        floor = FLOOR_SAL if k.startswith("sal") else FLOOR_FLUX
        excl = pnan if k in ("Fdil_S20", "Fint_S20", "salS20") else (hnan if k in RIM_FIELDS_L47 else None)   # L47：只 RIM 三场
        fields[k] = cmp_fields(ours[k], ref[f"f64_{k}"], floor, TOL_CAPSULE, excl)
        if f"f32_{k}" in ref:
            fields[k]["info_max_rel_vs_capsule_float32_inputs"] = cmp_fields(ours[k], ref[f"f32_{k}"], floor, 1.0, excl)["max_rel"]
    res["fields"] = fields
    res["rim_excluded_nan_history"] = rim_excluded_report(ours, ref, hnan)
    res["accumulation"] = check_accum(eng, ref, pnan, hnan)
    res["pass"] = bool(all(f["pass"] for f in fields.values()) and res["l25_pass"] and res["accumulation"]["pass"]
                       and res["rim_excluded_nan_history"]["positive_assertion"]["pass"])
    return res, ref


def pick_block(lay, st, meta, n=50):
    """真实块：按 1999-12-31 整日（48 个半步）有雨像元×半步数，在海洋占比 ≥95% 的 50×50 窗里取最大者（L37）。"""
    nla, nlo = meta["nlat"], meta["nlon"]
    ocean = np.zeros((nla, nlo), bool)
    ocean[np.asarray(st["row"]), np.asarray(st["col"])] = True
    cnt = np.zeros((nla, nlo), np.int32)
    for hh in range(24):
        out, _ = lay.read_cmorph(datetime(1999, 12, 31, hh, tzinfo=UTC))
        if out is not None:
            for a in out:
                cnt += (a > 0)
    cnt[~ocean] = 0

    def csum(a):
        c = np.zeros((a.shape[0] + 1, a.shape[1] + 1), np.int64)
        c[1:, 1:] = a.astype(np.int64).cumsum(0).cumsum(1)
        return c

    cc, co = csum(cnt), csum(ocean)
    best = None
    for r0 in range(0, nla - n, 10):
        for c0 in range(0, nlo - n, 10):
            box = lambda c: c[r0 + n, c0 + n] - c[r0, c0 + n] - c[r0 + n, c0] + c[r0, c0]
            if box(co) < 0.95 * n * n:
                continue
            v = box(cc)
            if best is None or v > best[0]:
                best = (int(v), r0, c0)
    if best is None:
        raise DataError("找不到海洋占比 ≥95% 的 50×50 块")
    return best


def real_block(lay, st, meta, n=50):
    """读真实块输入：1999-12-31T00 起 25 个小时文件（50 个半步；输出＝2000-01-01 00:00、00:30，同 main L16–L24）。"""
    rain, r0, c0 = pick_block(lay, st, meta, n)
    lat1d, lon1d = np.asarray(st["lat1d"]), np.asarray(st["lon1d"])
    rr, cc = np.meshgrid(np.arange(r0, r0 + n), np.arange(c0, c0 + n), indexing="ij")
    full = Domain(lat1d, lon1d, rr.ravel(), cc.ravel(), np.ones(n * n), np.zeros(n * n, np.float32), np.zeros(n * n), 1)
    fld = Fields(lay, full)
    P = np.full((n, n, 50), np.nan)
    U = np.full((n, n, 50), np.nan)
    for k in range(25):
        dt = datetime(1999, 12, 31, 0, tzinfo=UTC) + timedelta(hours=k)
        P0, P1 = fld.cmorph(dt)
        U0, U30 = fld.wind_pair(epoch_hour(dt))
        P[:, :, 2 * k], P[:, :, 2 * k + 1] = P0.reshape(n, n), P1.reshape(n, n)
        U[:, :, 2 * k], U[:, :, 2 * k + 1] = U0.reshape(n, n), U30.reshape(n, n)
    sst, ice, s0 = fld.day(date(2000, 1, 1))
    icefrac = fld.last_ice_frac.reshape(n, n)
    fld.close()
    # 斜率与 Watson：取自静态场（与全场同一来源）
    slope = np.full((n, n), np.nan, np.float32)
    rows, cols = np.asarray(st["row"]), np.asarray(st["col"])
    sel = (rows >= r0) & (rows < r0 + n) & (cols >= c0) & (cols < c0 + n)
    slope[rows[sel] - r0, cols[sel] - c0] = np.asarray(st["beta"])[sel]
    wn = Nearest(np.asarray(st["w_lat"]), np.asarray(st["w_lon"]), lat1d[r0:r0 + n], lon1d[c0:c0 + n])
    sf = wn(np.asarray(st["sf12"][0]))
    fa = wn(np.asarray(st["fa12"][0]))
    times = (np.datetime64("1999-12-31T00:00:00") + np.arange(50) * np.timedelta64(1800, "s")).astype("datetime64[s]").astype(np.int64)
    sst2 = sst.reshape(n, n)
    blk = {"lat": lat1d[r0:r0 + n], "lon": lon1d[c0:c0 + n], "times": times, "P": P, "U": U, "S0": s0.reshape(n, n),
           "sst": sst2, "icefrac": icefrac, "sf": sf, "fa": fa, "slope": slope}
    l46 = watson_is_pco2(meta)
    if l46:                                                              # L46：胶囊照 main L83 吃 fCO2，这里按同一 SST 换算
        with np.errstate(all="ignore"):
            blk["sf_p"], blk["fa_p"] = sf, fa
            blk["sf"], blk["fa"] = sf * fug_factor(sst2), fa * fug_factor(sst2 - SKIN_DT)
    return blk, {"row0": int(r0), "col0": int(c0), "lat_range": [float(blk["lat"][0]), float(blk["lat"][-1])],
                 "lon_range": [float(blk["lon"][0]), float(blk["lon"][-1])], "rain_pixel_halfsteps": rain,
                 "ocean_pixels": int(np.isfinite(slope).sum()), "watson_fco2_source": SRC_L46 if l46 else SRC_DIRECT}


def run_selftest(args, out_dir, lay, real_required):
    t0 = time.perf_counter()
    res = {"version": VERSION, "code_sha256": sha256_file(os.path.abspath(__file__)), "started_utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "tolerances": {"capsule_rel": TOL_CAPSULE, "flux_floor": FLOOR_FLUX, "sal_floor": FLOOR_SAL, "p2_rel": TOL_P2, "l25_rel": L25_TOL},
           "reads_p2_products": False, "parts": {}}
    parts = res["parts"]
    rng = np.random.default_rng(20260926)
    try:
        P, U = synth_series(rng)
        parts["p2_rim_synthetic"] = test_p2_equiv(P, U, "synthetic")
    except Exception as e:
        parts["p2_rim_synthetic"] = {"ran": True, "pass": False, "error": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()[-2000:]}
    try:
        parts["fugacity_pyco2sys"] = test_fugacity()
    except Exception as e:
        parts["fugacity_pyco2sys"] = {"ran": True, "pass": False, "error": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()[-2000:]}
    try:
        parts["l25_grid"] = test_l25_grid()
    except Exception as e:
        parts["l25_grid"] = {"ran": True, "pass": False, "error": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()[-2000:]}
    cap_dir, how = find_capsule(args, lay, out_dir)
    res["capsule"] = {"dir": cap_dir, "how": how, "sha256": CAPSULE_FILES}
    if cap_dir is None:
        parts["capsule_synthetic"] = {"ran": True, "pass": False, "error": how}
    else:
        try:
            parts["capsule_synthetic"], _ = test_capsule_block(args, cap_dir, synth_block(rng), out_dir, "synthetic")
        except Exception as e:
            parts["capsule_synthetic"] = {"ran": True, "pass": False, "error": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()[-2000:]}
    real_ok = lay is not None and lay.month_ready(200001) is not None and os.path.exists(
        os.path.join(args.work_dir, "static", "meta.json"))
    if real_ok and cap_dir is not None:
        try:
            st, meta = load_static(os.path.join(args.work_dir, "static"))
            blk, binfo = real_block(lay, st, meta)
            gl = lay.load_glodap(tempfile.mkdtemp(dir=args.work_dir))
            ok = np.flatnonzero(np.isfinite(gl["TAlk"]) & np.isfinite(gl["TCO2"]) & np.isfinite(gl["salinity"])
                                & np.isfinite(gl["temperature"]))
            pick = np.random.default_rng(7).choice(ok, size=min(300, len(ok)), replace=False)
            low = ok[gl["salinity"].ravel()[ok] <= 9][:20]
            pick = np.r_[pick, low]
            g = {k: gl[v].ravel()[pick] for k, v in (("g_ta", "TAlk"), ("g_tco2", "TCO2"), ("g_sst", "temperature"), ("g_sss", "salinity"))}
            extra = {**g, "area_lon": np.asarray(st["lon1d"]), "area_lat": np.asarray(st["lat1d"])}
            part, ref = test_capsule_block(args, cap_dir, blk, out_dir, "real", extra)
            part["block"] = binfo
            L_ = blk["P"].shape
            Pt = blk["P"].reshape(-1, L_[2]).T
            Ut = blk["U"].reshape(-1, L_[2]).T
            oc = np.isfinite(blk["slope"]).ravel()
            part["p2_rim_real"] = test_p2_equiv(Pt[:, oc], Ut[:, oc], "real_block")
            if ref is not None:
                ours_s = glodap_slopes(g["g_ta"], g["g_tco2"], g["g_sss"], g["g_sst"]).astype(np.float64)
                part["glodap_slopes"] = cmp_fields(ours_s, ref["slopes"].astype(np.float32).astype(np.float64), 1e-12, TOL_CAPSULE)
                part["glodap_slopes"]["float32_equal_frac"] = float(np.mean(ours_s == ref["slopes"].astype(np.float32).astype(np.float64)))
                part["area"] = cmp_fields(grid_cell_areas(np.asarray(st["lon1d"]), np.asarray(st["lat1d"])), ref["area"], 1e-12, 1e-12)
                part["pass"] = bool(part["pass"] and part["p2_rim_real"]["pass"] and part["glodap_slopes"]["pass"] and part["area"]["pass"])
            else:
                part["pass"] = False
            parts["real_block"] = part
        except Exception as e:
            parts["real_block"] = {"ran": True, "pass": False, "error": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()[-3000:]}
    else:
        parts["real_block"] = {"ran": False, "pass": not real_required,
                               "why": "无 --cache-root 或 200001 未就绪/静态场未建" if cap_dir else "胶囊不可用"}
    res["pass"] = bool(all(p.get("pass") for p in parts.values()))
    res["seconds"] = round(time.perf_counter() - t0, 1)
    jdump(res, os.path.join(out_dir, "selftest.json"))
    log(f"自测 {'通过' if res['pass'] else '不通过'}（{res['seconds']} s）："
        + "，".join(f"{k}={'过' if v.get('pass') else ('未跑' if not v.get('ran', True) else '不过')}" for k, v in parts.items()))
    return res["pass"]


# ======================================================================== bench 外推（L41）
def lpt_makespan(jobs, W):
    loads = [0.0] * W
    for j in sorted(jobs, reverse=True):
        k = loads.index(min(loads))
        loads[k] += j
    return max(loads)


def bench_report(info, t_static, t_selftest, out_dir):
    tm = info["timing_s"]
    steps = max(info["steps"], 1)
    t_step = (tm.get("read_cmorph", 0) + tm.get("wind", 0) + tm.get("steps", 0) - tm.get("eng_l25_check", 0)) / steps  # full 不做 L25 抽检
    t_day = (tm.get("day_read", 0) + tm.get("day_coefs", 0)) / max(info["days"], 1)
    t_spin = tm.get("spinup", 0) / 48.0
    months = [(m, len(days_of_month(YEAR, m))) for m in range(1, 13)]
    per_month = [nd * (48 * t_step + t_day) + 48 * t_spin for _, nd in months]
    serial = sum(per_month)
    ext = {str(W): round((lpt_makespan(per_month, W) + t_static + t_selftest) / 3600, 2) for W in (1, 2, 3, 4, 6, 8)}
    rss = info["peak_rss_bytes"]
    phases = {k: round(v / steps, 4) for k, v in tm.items() if k.startswith("eng_")}
    rep = {"version": VERSION, "day": "2000-01-01（起转 1999-12-31）", "steps": info["steps"], "N": info["N"],
           "seconds_per_step": round(t_step, 3), "seconds_per_day_overhead": round(t_day, 2), "seconds_per_spinup_step": round(t_spin, 3),
           "engine_phase_seconds_per_step": phases, "static_build_s": round(t_static, 1), "selftest_s": round(t_selftest, 1),
           "serial_year_hours": round(serial / 3600, 2), "wall_hours_by_workers": ext,
           "peak_rss_GiB_single_process": round(rss / 2**30, 3),
           "entries_per_step": info["entries_per_step"], "D_frac": info["D_frac"],
           "l25_max_rel_vs_exact_gsw": info["l25_max_rel"], "l25_n_checked": info["l25_n_checked"],
           "l25_pass": bool(all(v < L25_TOL for v in info["l25_max_rel"].values())) if info["l25_max_rel"] else None,
           "watson_fco2_source": info.get("watson_fco2_source"), "watson_fugfac_day1": info.get("watson_fugfac"),
           "suggest_jobspec_p3_global": {
               "workers": 4, "ram_bytes": int(math.ceil((4 * rss * 1.3 + 2**30) / 2**30)) * 2**30,
               "estimate_seconds": int(ext["4"] * 3600), "hard_walltime_seconds": int(3 * ext["4"] * 3600),
               "note": "ram＝4×单进程峰值×1.3＋1 GiB（单进程峰值含静态场与自测，偏保守）；walltime＝估计×3（执行计划 §3 E）；未含等下载"}}
    jdump(rep, os.path.join(out_dir, "bench.json"))
    log(f"基准：每步 {t_step:.2f} s，日开销 {t_day:.1f} s，峰值 RSS {rss / 2**30:.2f} GiB；整年墙钟（4 worker）约 {ext['4']} h")
    return rep


# ======================================================================== 主程序
def main():
    global _LOG_PATH
    ap = argparse.ArgumentParser(description="P3 全球重算（P3a＋P3b）")
    ap.add_argument("--mode", choices=("selftest", "bench", "full", "capsule-ref"), required=True)
    ap.add_argument("--cache-root", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--work-dir", default=None,
                    help="静态场与月结果目录（默认 <out>/_work）；跨次续跑须给固定路径（REPRO_OUTPUT_DIR 可能每次不同）")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--months", default=None)
    ap.add_argument("--chunk", type=int, default=CHUNK)
    ap.add_argument("--wait-hours", type=float, default=30.0)
    ap.add_argument("--l25-every", type=int, default=12)
    ap.add_argument("--capsule-dir", default=None)
    ap.add_argument("--capsule-python", default=None)
    ap.add_argument("--ref-in", default=None)
    ap.add_argument("--ref-out", default=None)
    args = ap.parse_args()
    if args.mode == "capsule-ref":
        return run_capsule_ref(args)

    out_dir = os.environ.get("REPRO_OUTPUT_DIR") or args.out
    if not out_dir:
        print("需要 $REPRO_OUTPUT_DIR 或 --out", file=sys.stderr)
        return 3
    os.makedirs(out_dir, exist_ok=True)
    args.work_dir = os.path.abspath(args.work_dir) if args.work_dir else os.path.join(out_dir, "_work")
    os.makedirs(args.work_dir, exist_ok=True)
    _LOG_PATH = os.path.join(out_dir, "p3_log.txt")
    code_sha = sha256_file(os.path.abspath(__file__))
    log(f"{VERSION}（sha256 {code_sha[:12]}…）mode={args.mode} out={out_dir}")
    lay = Layout(args.cache_root) if args.cache_root else None
    tmpdir = tempfile.mkdtemp(prefix="p3tmp_", dir=args.work_dir)
    try:
        t_static = 0.0
        if args.mode in ("bench", "full"):
            if lay is None:
                log("bench/full 需要 --cache-root")
                return 3
            wait_ready_static(lay, args.wait_hours * 3600)
            t0 = time.perf_counter()
            build_static(lay, args.work_dir, tmpdir)
            t_static = time.perf_counter() - t0
        elif lay is not None and lay.month_ready(200001) is not None:
            build_static(lay, args.work_dir, tmpdir)
        t0 = time.perf_counter()
        ok = run_selftest(args, out_dir, lay, real_required=(args.mode == "bench"))
        t_self = time.perf_counter() - t0
        if args.mode == "selftest":
            return 0 if ok else 4
        sd = os.path.join(args.work_dir, "static")
        _, meta = load_static(sd)
        core = {"version": VERSION, "code_sha256": code_sha, "chunk": args.chunk, "cache_root": lay.root,
                "static_N": meta["N"], "nreg": meta["nreg"]}
        chash = config_hash(core)
        if args.mode == "bench":
            cfg = {"ym": 200001, "cache_root": lay.root, "static_dir": sd, "work_dir": args.work_dir, "chunk": args.chunk,
                   "wait_s": args.wait_hours * 3600, "max_days": 1, "l25_every": args.l25_every, "config_hash": chash, "resume": False}
            r = process_month(cfg)
            _LOG_PATH = os.path.join(out_dir, "p3_log.txt")
            if not r["ok"]:
                log(f"bench 失败：{r}")
                return 2 if r.get("kind") == "data" else 3
            with open(r["info_path"], encoding="utf-8") as f:
                info = json.load(f)
            os.replace(r["info_path"], os.path.join(args.work_dir, "bench_month_200001_day1.json"))
            os.replace(r["info_path"][:-5] + ".npz", os.path.join(args.work_dir, "bench_month_200001_day1.npz"))
            rep = bench_report(info, t_static, t_self, out_dir)
            if not ok:
                return 4
            return 0 if rep["l25_pass"] in (True, None) else 4
        # ---- full ----
        if not ok:
            log("等价性自测不过：不做全年计算")
            return 4
        months = [YEAR * 100 + int(x) for x in args.months.split(",")] if args.months else [YEAR * 100 + m for m in range(1, 13)]
        cfgs = [{"ym": ym, "cache_root": lay.root, "static_dir": sd, "work_dir": args.work_dir, "chunk": args.chunk,
                 "wait_s": args.wait_hours * 3600, "max_days": None, "l25_every": 0, "config_hash": chash, "resume": True}
                for ym in months]
        t0 = time.perf_counter()
        # ProcessPoolExecutor：worker 被杀（如 OOM）时抛 BrokenProcessPool 而不是像 mp.Pool 那样永久挂起
        ctx = mp.get_context("spawn")
        results = []
        with cf.ProcessPoolExecutor(max_workers=max(1, min(args.workers, len(cfgs))), mp_context=ctx,
                                    max_tasks_per_child=1) as ex:
            futs = {ex.submit(process_month, c): c["ym"] for c in cfgs}
            for fu in cf.as_completed(futs):
                try:
                    r = fu.result()
                except Exception as e:
                    r = {"ym": futs[fu], "ok": False, "kind": "error", "error": f"worker 异常退出：{type(e).__name__}: {e}"}
                results.append(r)
                log(f"月 {r['ym']}：{'完成' if r['ok'] else '失败 ' + r.get('error', '')}")
        bad = [r for r in results if not r["ok"]]
        run = {"version": VERSION, "code_sha256": code_sha, "config_hash": chash, "months": months, "results": results,
               "wall_seconds": round(time.perf_counter() - t0, 1), "static_build_s": round(t_static, 1), "selftest_s": round(t_self, 1)}
        if bad:
            jdump(run, os.path.join(out_dir, "p3_run.json"))
            log(f"有 {len(bad)} 个月失败，不汇总（可原样重提续跑）")
            return 2 if all(r.get("kind") == "data" for r in bad) else 3
        if sorted(months) != [YEAR * 100 + m for m in range(1, 13)]:
            jdump(run, os.path.join(out_dir, "p3_run.json"))
            log("--months 非全年：只留月结果，不写 p3a/p3b")
            return 0
        inputs = {"month_ready": {ym: (lay.month_ready(ym) or {}).get("cmorph_missing") for ym in months}}
        repro, curves, infos = aggregate(args.work_dir, months, meta, out_dir, {"code_sha256": code_sha, "inputs": inputs})
        run["month_info"] = {ym: {k: infos[ym][k] for k in ("seconds", "peak_rss_bytes", "timing_s", "entries_per_step", "D_frac")}
                             for ym in infos}
        jdump(run, os.path.join(out_dir, "p3_run.json"))
        log(f"完成：复现门 {'全过' if repro['all_pass'] else '未全过'}；U(1)={repro['table_pct_of_wind']['U1_reconstructed']:.3f} pp")
        return 0
    except DataError as e:
        log(f"数据错误：{e}")
        return 2
    except StopAsk as e:
        log(f"停下，需人工处理：{e}")
        return 3
    except Exception as e:
        log(f"异常：{type(e).__name__}: {e}\n{traceback.format_exc()}")
        return 3
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def wait_ready_static(lay, wait_s):
    t_end = time.time() + wait_s
    while lay.month_ready(200001) is None:
        if time.time() > t_end:
            raise DataError("month_ready/200001.json 未写（静态源与起转日未齐）")
        log("等待 month_ready/200001.json（静态源与起转日）")
        time.sleep(POLL_S)


if __name__ == "__main__":
    sys.exit(main())
