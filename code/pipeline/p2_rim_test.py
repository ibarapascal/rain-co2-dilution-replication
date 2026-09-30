#!/usr/bin/env python3
"""P2：CMORPH 驱动 RIM-3／S20 的系泊检验——D3 判别力、R_RIM(1 m)、R_S20(0.5 m)、P5、β，按事先写定的规则判出口。

设计（运行前写定，不改阈值与统计方法）；门 D1/D2/D4/D5 沿用 P1b 合并集判定；实现选择 K1–K30 见下。
RIM-3 与 S20、β_model 的实现以 Witte 等 2026 的 Code Ocean 代码为准（`code/CO2_Rain_Flux_Toolbox.py` 的 RIMv3／S20／
calculate_fco2_dilution_per_psu）；深度因子取 Witte et al. 式 (4)。

用法：
  python p2_rim_test.py                      # 全量（输出到 $REPRO_OUTPUT_DIR）；等价 --part all
  python p2_rim_test.py --part primary       # 只做主集：只下主集 CMORPH 小时，出 p2_summary.json／p2_d3.json／事件 CSV（含 K27），不跑敏感性
  python p2_rim_test.py --part sensitivity --resume-from <主集输出目录>
                                             # 只跑四项敏感性（K21）→ p2_sensitivity.json；不判出口、不写 summary；
                                             # 须 --resume-from 指向含 p2_d3.json 的主集输出目录（或本段先前输出目录，含 p2_d3_from_primary.json）
  python p2_rim_test.py --plan               # 只重建事件（读 P1/P1b 缓存）并算唯一小时数、预计下载量（仅 HEAD 请求，不下大文件）；
                                             # p2_plan.json 的 parts 分别给 primary／sensitivity 两段的唯一小时数与字节估计
  python p2_rim_test.py --selftest           # 合成数据自测（无网络）
  python p2_rim_test.py --smoke              # 冒烟：只取主集前 3 个事件的 CMORPH（约 120 个文件），结论无效
  --resume-from DIR   继承上次 P2 输出目录里的 cmorph_pixels.csv 与 p2_glodap.json（断点续跑；DIR 不存在即退出 3）
  --p1-events PATH    P1 产物 p1_events.csv（默认 p1 阶段产物）；同目录 cache/ 作只读种子缓存
  --p1b-dir DIR       P1b 产物目录（p1b_summary.json、p1b_events_new.csv、cache/）
  --max-conc N        CMORPH 并发下载数（≤4）

依赖：numpy、scipy、netCDF4、PyCO2SYS、matplotlib＋同目录 p1_events.py、p1b_extend.py（事件定义复用，不改）。
数据（匿名 HTTP，无需登录）：
  - CMORPH CDR 8 km/30 min：NCEI access/30min/8km/YYYY/MM/DD/CMORPH_V1.0_ADJ_8km-30min_YYYYMMDDHH.nc（每小时约 1.64 MB）
  - GLODAPv2.2016b 气候态：glodap.info/glodap_files/v2.2023/GLODAPv2.2016b.MappedProduct.tar.gz（211 MB，流式只取 4 个成员）
  - P1/P1b 已缓存的 MAPCO2／GTMBA／OceanSITES（以符号链接引用，只读；缺失时才联网补取）
产物（<out>/）：all／primary：p2_summary.json、p2_d3.json（先于观测比值写盘）、p2_events.csv、p2_plan.json、p2_glodap.json、
  cmorph_pixels.csv（小时级像元缓存）、p2_fig_ratio.png、p2_log.txt；K29 中止时只有 p2_abort.json（无 d3/summary）。
  sensitivity：p2_sensitivity.json、p2_plan_sensitivity.json、p2_d3_from_primary.json（主集 D3 的只读副本）、cmorph_pixels.csv、p2_log.txt。
退出码：0＝跑完（无论出口）；2＝数据源故障（可 --resume-from 续跑，不产出 summary）；3＝其他异常（含事件重建校验不一致、
  K28 缺 s1_depth_m、sensitivity 段缺主集 D3）或 K29 可用率保护中止（以 p2_abort.json 区分）。

实现选择清单（设计未写死之处；均未放宽任何门）：
  K1 事件重建：P1 八站复用 p1_events.process_station（D/M 文件），O 组五站复用 p1b_extend.discover_o＋process_o_station，
     读 P1b／P1 缓存（符号链接播种）；R 组（P1b 新增 0）不重建。10 mm 事件逐条与 p1_events.csv、p1b_events_new.csv 比对
     （(站,起点) 键集合、ΔS1／ΔS0.5 四位小数、对照数），且总数＝646、D4 计数＝34，否则退出 3。
  K2 d0：Witte 代码 RIMv3 里的 DL 是 7×7（U=[0,2,4,6,8,10,200]、R=[0,2,5,10,20,50,200]，首末行列为复制的边界），
     不是 6×6；按 xarray DataArray.interp 默认（linear）做双线性插值，坐标外（>200）为 NaN，不外推。
  K3 RIM-3 公式以代码为准：48 个历史半步＋当前半步；tc=600:25:1775（48 个）、ti=86400…1800；当前项 c2·IRR·1800/√Kz。
     深度因子取 Witte et al. 式 (4) 的 exp(−z²/(4·Kz·t))：历史项 t=ti；当前项沿用代码隐含的 t=1 s（代码 √Kz 即 √(Kz·1 s)），
     故 z≥0.5 m 时当前项≈1。z=0 时与代码逐值一致（自测对照逐行移植版）。
  K4 原文描述与代码的出入（原文当前项写 √(Kz·t)、称 tc 为 49 元素向量；代码为 48＋1，当前项无 t）：采用代码。
  K5 时间对齐：CMORPH time＝半小时起点；模型「小时值」＝该 UTC 小时两个半步（HH:00、HH:30）输出的平均，对应观测 I4 小时箱。
  K6 模型 ΔS 与观测同一定义（I8）：ΔS_RIM(z)＝S0(z)·[mean F(z, 0–6 h) − median F(z, 前 6 h)]，F＝RIM 稀释因子，
     S0(z)＝该深度观测雨前 6 h 中位数；ΔS_S20＝小时平均 ΔS_S20 的同样差值；S20 按 Witte 下限 S≥0。
  K7 CMORPH 像元＝站点坐标的最近格点（中心像元）；另存 3×3 均值只作记录、不进任何分析。站点坐标：P1 八站取 GTMBA
     名义站位（TAO110W 取 0°,110°W＝雨量／盐度所在），O 组取 P1b 发现记录的 MAPCO2 ERDDAP 名义中点。
  K8 CMORPH 缺测（掩码、HTTP 404、时间轴不符）→ 该半步 NaN → 该事件模型预测 NaN，从比值剔除并计数；不插补。
  K9 风：用浮标小时风，按中性对数律由 4 m 换算到 10 m（×ln(10/z0)/ln(4/z0)，z0=1e-4 m，≈1.0865），全站统一按 4 m；
     Witte 用 ERA5 U10（需 CDS 登录，未用）。缺测只在 ≤6 h 且两端有值时线性插补，否则该事件模型 NaN；U<0.1 m/s 置 0.1
     （避免 Kz=0 除零，计数报告）。小时风重复到两个半步。
  K10 主事件集＝10 mm 合并 646 事件中「该深度观测 ΔS 与模型 ΔS 都有效」者；不要求配满 5 对照（对照只用于 D3 的 SE）。
  K11 D3：z∈{0.5, 1} m 各自比较 |mean ΔS_RIM(z) − mean ΔS_S20| 与 2·SE_obs(z)，SE_obs(z)＝该深度全部对照窗 ΔS_obs(z)
     的样本 SD／√N_events(z)。两个深度都 ≥2·SE 才过；任一不足即「不确定」（从严，逐深度结果都写出）。D3 在 P2 新算的
     事件观测量（ΔS(5 m)、ΔpCO₂、各窗口观测 ΔS）与任何比值之前写盘（p2_d3.json）。事件 ΔS(1 m)／ΔS(0.5 m) 已在 P1/P1b
     算出并公开，K1 校验会重算它们、但 D3 不读；模型 ΔS 用到的 S0（雨前 6 h 观测中位数）是换算尺度，不是结果量。
  K12 CI：站×季整簇 bootstrap（同 I10：B=10000，seed=20260926，numpy PCG64），比值＝重抽后 Σ分子/Σ分母，百分位区间。
  K13 TOST：α=0.05、界 [0.65, 1.35]；p_lower＝P_boot(R≤0.65)，p_upper＝P_boot(R≥1.35)，两者 <0.05 即等价（⇔90% 区间在界内）；
     出口 2 另按设计要求 95% CI 在界内。MDE＝2.80×SE_boot(R_RIM)（双侧 α=0.05、80% 功效）。
  K14 「R_S20 被排除」＝R_S20(0.5 m) 的 95% CI 完全落在 [0.65, 1.35] 之外（与 R_RIM 同一 SESOI，从严）；CI 是否排除 1 只作信息。
  K15 出口顺序：出口 3 条件（D1/D2/D5 任一不过或 D3 不确定）优先；否则出口 1 条件成立即出口 1（另记 D4
     是否过）；否则出口 2 条件；都不成立（CI 跨 SESOI 边界）记「未落入事先写定的出口」，不自行归类。
     出口 1 与 2 条件同时成立时两者都报、主标出口 1。
  K16 D1/D2/D4/D5 不重判：读 p1b_summary.json 的合并集判定（须 merged=646 且四门全过，否则退出 3），P2 只加 D3。
  K17 P5：观测 ΣΔS_obs(5 m)/ΣΔS_obs(1 m) 与模型 ΣΔS_RIM(5)/ΣΔS_RIM(1) 用同一事件集（5 m 雨前与 0–6 h 均有值、模型有效）
     与同一 bootstrap 抽样；报两者及差值 CI。
  K18 β_obs：pCO₂ 温度归一 pCO₂·exp(0.0423·(T_pre−T))（Takahashi 1993），T＝MAPCO2 同小时 SST，T_pre＝雨前 6 h SST 中位数；
     ΔpCO₂＝0–6 h 均值 − 雨前 6 h 中位数；β_obs＝ΣΔpCO₂/ΣΔS(0.5 m)。主集＝D4 子集（ΔS0.5≤−0.2，设计定义的 pCO₂ 通道），
     全事件集只作信息。
  K19 β_model：照 Witte calculate_fco2_dilution_per_psu（S 到 S−9 共 10 点、TA·SR、DIC·SR＋25(1−SR)、opt_k_carbonic=17、
     其余 PyCO2SYS 默认，linregress 斜率），输入 GLODAPv2.2016b 表层（depth_surface 第 0 层）最近 1° 格点的 TA/DIC/S/T；
     与观测比较用 pCO₂ 输出的斜率（观测为 pCO₂），fCO₂ 斜率（Witte 原值）并报；多站合并＝Σβ_st·ΔS/ΣΔS；格点 NaN 不回填。
  K20 站群分层（P1b 披露限制 3）：A＝P1 热带 8 站（TAO/RAMA），B＝O 组 5 站；各群单独报 R_RIM、R_S20、P5、β 的点估计与
     簇 bootstrap CI；逐站只报点估计与 n；B 群「1 m 替代层」事件数另报。
  K21 敏感性只做设计列出的四项：阈值 5／20 mm（p1 同一事件定义，只改累积阈值）；窗口 0–3、6–12 h（观测与模型同窗）；
     雨量计驱动（小时雨量重复到两个半步，负值与缺测置 0 并计数）；剔除 BOBOA。每项只报 R_RIM(1 m)、R_S20(0.5 m) 的
     点估计与 95% CI，不参与出口判定。
  K22 CMORPH 小时窗：每事件 [t0−30 h, t0+12 h)（雨前 6 h 基线需 24 h 雨史；+12 h 覆盖 6–12 h 敏感窗），三个阈值事件集取并；
     主集（10 mm）的小时先下。
  K23 下载：唯一小时逐文件下载（≤4 并发，线程内持久连接）→ 主线程抽 13 站像元 → 立即删除；小时级缓存 cmorph_pixels.csv
     （每小时一块＋DONE 标记，重跑跳过）；单文件至多 5 次尝试（退避 5/10/20/40 s），404 连续 2 次记缺测；连续 40 个文件失败、
     或末轮补试仍失败 → 退出 2（不写 summary）。
  K24 GLODAP：NCEI 归档链接 2026-09-26 超时，改用 glodap.info 同版 tar（GLODAPv2.2016b.MappedProduct.tar.gz），流式只取
     TAlk/TCO2/salinity/temperature 四个成员，读完即删。
  K25 描述性信息：报告对照窗 ΔS_obs(1 m, 0–6 h) 均值（比值按设计不减对照；不另算对照校正比值）。
  K26 出口判定只用 CMORPH 驱动、10 mm、0–6 h、全站的主结果；P5、β、分层与敏感性都不进出口判定。
  K27 探索性变体（不参与出口、不进 D3）：当前项深度衰减因子 exp(−z²/(4·Kz·t)) 取
     t=1800 s（量级项 c2·IRR·1800/√Kz 不变），只对 CMORPH、10 mm、0–6 h 报 R_RIM(0.5 m)、R_RIM(1 m) 的点估计与同一
     bootstrap 的 95% CI、按站群分层点估计，写入 p2_summary.json 的 exploratory_k27；主结果数值路径不变。
  K28 替代层深度：凡模型计算「1 m 观测层」处（R_RIM(1 m)、D3 的 1 m、P5 分母、
     K27 的 1 m 变体、敏感性各项的 1 m），深度因子用该事件实际观测深度 s1_depth_m（P1 热带 8 站＝1.0；O 组＝该站 ≤1.5 m
     最浅层，事件窗内中位数，来自 p1b_extend.run_events）；0.5 m、5 m 层仍按名义深度。缺 s1_depth_m（或非有限、不在
     (0, 1.5] m）即报错退出 3，不回退 1.0（事件重建后、下载前对三阈值全部事件先校验，--plan 亦然，早失败）。
     事件 CSV 记 model_depth_s1_layer_m；summary 另报按「替代层／真 1 m」分层的
     R_RIM(1 m) 点估计（k28_R_RIM_1m_by_layer，仅信息，不进出口）。
  K29 可用率保护：主集（10 mm 全部事件，CMORPH、0–6 h）中模型可用（ΔS_RIM(1 m 层)
     与 ΔS_S20 均非 NaN）的比例 <50% → 写 p2_abort.json（总数、可用数、各缺测原因计数），不写 p2_d3.json 与
     p2_summary.json，不做出口判定，退出 3。检查放在 D3 之前（中止时 D3 也不写盘，不产生任何出口结论）。
     sensitivity 段不判出口，只把同一可用率写进 p2_sensitivity.json，不中止。
  K30 已知限制（不改）：(a) 观测 0.5 m（MAPCO2）多为 3 小时一次，模型 ΔS 用满 6 小时均值，期望无偏但增噪；
     (b) 模型小时值（HH:00、HH:30 两半步输出的平均，各为该半步末的状态）代表 HH:30–HH+1:00，观测在 HH:00–HH:17，
     模型相对观测超前约 15–45 分钟。两条都只作披露。
  分段（--part）：all＝原行为；primary＝只下主集小时、只算主结果＋K27（事件 CSV 的雨量计列留空），不跑敏感性；
     sensitivity＝只算四项敏感性（K21），不做 GLODAP／D3／出口，要求 D3 已由主集写盘（K11 顺序）。

Change Log:
  2026-09-26 初版（复用 p1_events.py／p1b_extend.py，不改它们）。
  2026-09-26 b：加 K27 探索性变体（rim_factor 增 t_cur_depth 参数，默认值即主路径）；K11 措辞按实际顺序改准；
     补计数（雨量计负值置 0、CMORPH 半步 >200 mm/h 致模型 NaN）；主结果数值路径不变。
  2026-09-26 c：K28（模型「1 m」层改用事件实际深度 s1_depth_m，主路径数值只因此变）、K29（主集模型可用率
     <50% 中止，exit 3＋p2_abort.json）、K30（已知限制记录）；加 --part {all,primary,sensitivity} 分段与 --plan 分段估计；
     事件 CSV 字段 wind_floor_halfsteps→wind_floor_hours、summary wind.halfsteps_floored→hours_floored（数的是小时）；
     --resume-from 目录不存在即退出 3；自测增 K28/K29/分段/端到端。
"""

import argparse
import concurrent.futures as cf
import csv
import glob
import gzip
import http.client
import json
import math
import os
import shutil
import statistics
import sys
import tarfile
import threading
import time
import traceback
import urllib.request
from datetime import datetime, timezone

import p1_events as p1
import p1b_extend as p1b

VERSION = "p2-2026-09-26c"

P1_EVENTS_DEFAULT = p1b.P1_EVENTS_DEFAULT
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
import repro_io as _rio  # [repro] 原件落 raw 档的取数小工具
P1B_DIR_DEFAULT = _rp.upstream("p1b_dir")  # [repro] 读 p1b 阶段输出目录
EXPECTED_MERGED = 646
EXPECTED_D4 = 34
EXPECTED_O_STATIONS = ["KEO", "MOSEAN/WHOTS", "Papa", "SOFS", "Stratus"]

# ---- 事先写定的参数（不得改） ----
SESOI_LO, SESOI_HI = 0.65, 1.35
D3_K_SE = 2.0
PRIMARY_WIN = (0, 6)
SENS_WINDOWS = [(0, 3), (6, 12)]
SENS_THRESHOLDS = [5.0, 20.0]
PRIMARY_THRESHOLD = 10.0
DEPTHS = (0.5, 1.0, 5.0)

# ---- Witte 代码常数（CO2_Rain_Flux_Toolbox.py RIMv3/S20/calculate_fco2_dilution_per_psu） ----
C1 = 5.7
C2 = 1 / 11
TC = [600 + 25 * k for k in range(48)]                 # np.arange(600,1800,25)
TI = [86400 - 1800 * k for k in range(48)]             # -np.arange(0,24*3600,1800)+24*3600
KZ_COEF = 2.5e-5
T_CURRENT_S = 1.0                                      # K3
K27_T_CUR_DEPTH_S = 1800.0                             # K27（探索性，只改当前项深度因子）
K27_DEPTHS = (0.5, 1.0)
K27_TAG = "cmorph_k27"
U_GRID = [0, 2, 4, 6, 8, 10, 200]
R_GRID = [0, 2, 5, 10, 20, 50, 200]
DL = [[1.2, 1.2, 1.0, 1.0, 1.0, 1.1, 1.1],
      [1.2, 1.2, 1.0, 1.0, 1.0, 1.1, 1.1],
      [1.6, 1.6, 1.9, 2.0, 2.0, 2.2, 2.2],
      [1.5, 1.5, 2.5, 2.9, 3.1, 3.5, 3.5],
      [1.9, 1.9, 2.7, 3.6, 4.2, 4.7, 4.7],
      [2.4, 2.4, 3.0, 4.0, 5.0, 5.8, 5.8],
      [2.4, 2.4, 3.0, 4.0, 5.0, 5.8, 5.8]]
S20_A, S20_B = -0.35, 0.77
DIC_RAIN = 25.0
K_CARBONIC = 17

# ---- 实现层参数（K 条） ----
BOOT_B = 10000
BOOT_SEED = 20260926
WIND_Z_M, WIND_Z0_M = 4.0, 1e-4                        # K9
WIND_FACTOR = math.log(10.0 / WIND_Z0_M) / math.log(WIND_Z_M / WIND_Z0_M)
WIND_FLOOR = 0.1
WIND_GAP_MAX_H = 6
HIST_H = 24
WIN_LO_H = -(p1.PRE_H + HIST_H)                        # −30
WIN_HI_H = 12                                          # K22
N_MODEL_H = p1.PRE_H + WIN_HI_H                        # 18 个模型小时（t0−6 … t0+11）
T_COEF = 0.0423                                        # K18
MDE_Z = 1.959964 + 0.841621                            # K13
K29_MIN_AVAIL = 0.5                                    # K29：主集模型可用率下限
K28_MAX_DEPTH = p1b.MAX_S1_DEPTH                       # K28：「1 m」层实际深度上限（1.5 m）
PCA_FILES_PER_S = 2.0                                  # 到 NCEI 4 并发实测约 2 文件/s（参考运行的 plan 模式）
PARTS = ("all", "primary", "sensitivity")
MAX_CONC = 4
MAX_ATTEMPTS = 5
MAX_CONSEC_FAIL = 40
CMORPH_BASE = "https://www.ncei.noaa.gov/data/cmorph-high-resolution-global-precipitation-estimates/access/30min/8km"
CMORPH_HOST = "www.ncei.noaa.gov"
CMORPH_MEAN_MB = 1.64                                  # 2026-09-26 对 12 个年份文件 HEAD 的均值
GLODAP_URL = "https://glodap.info/glodap_files/v2.2023/GLODAPv2.2016b.MappedProduct.tar.gz"
GLODAP_VARS = {"TAlk": "GLODAPv2.2016b.TAlk.nc", "TCO2": "GLODAPv2.2016b.TCO2.nc",
               "salinity": "GLODAPv2.2016b.salinity.nc", "temperature": "GLODAPv2.2016b.temperature.nc"}
UA = "rain-co2-dilution-replication-p2/1.0 (research; python)"

# P1 八站 CMORPH/GLODAP 取点坐标（GTMBA 名义站位，K7）
P1_COORDS = {"TAO165E": (0.0, 165.0), "TAO8S165E": (-8.0, 165.0), "BOBOA": (15.0, 90.0), "TAO170W": (0.0, -170.0),
             "TAO155W": (0.0, -155.0), "TAO140W": (0.0, -140.0), "TAO125W": (0.0, -125.0), "TAO110W": (0.0, -110.0)}
GROUP_A, GROUP_B = "A_P1热带8站", "B_O组5站"
NAN = float("nan")


class DataSourceError(Exception):
    """网络/数据源故障：退出 2，可续跑。"""


# ======================================================================== 小工具
def _np():
    import numpy as np
    return np


def isnum(x):
    return x is not None and x == x and not (isinstance(x, float) and math.isinf(x))


def rnd(x, nd=5):
    return round(float(x), nd) if isnum(x) else None


def nanmean(xs):
    v = [x for x in xs if x == x]
    return sum(v) / len(v) if v else NAN


def nanmedian(xs):
    v = [x for x in xs if x == x]
    return statistics.median(v) if v else NAN


def sha256_file(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def jdump(obj, path):
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(p1b.clean(obj), f, ensure_ascii=False, indent=2, allow_nan=False, default=str)
    os.replace(tmp, path)


# ======================================================================== 缓存（符号链接播种，P1/P1b 缓存只读）
class LinkSeedHttp(p1.Http):
    """dest 不存在时，若任一种子缓存同相对路径有文件，建符号链接（不复制、不改种子）。"""

    def __init__(self, log, cache, seeds):
        super().__init__(log)
        self.cache = os.path.abspath(cache)
        self.seeds = [os.path.abspath(s) for s in seeds if s and os.path.isdir(s)]
        self.n_linked = 0

    def get(self, url, dest, empty_ok=False):
        if not os.path.exists(dest) and not os.path.exists(dest + ".empty"):
            rel = os.path.relpath(os.path.abspath(dest), self.cache)
            if not rel.startswith(".."):
                for seed in self.seeds:
                    hit = False
                    for suf in ("", ".empty"):
                        src = os.path.join(seed, rel) + suf
                        if os.path.exists(src):
                            os.makedirs(os.path.dirname(dest), exist_ok=True)
                            if suf:
                                open(dest + suf, "w").close()
                            else:
                                os.symlink(os.path.realpath(src), dest)
                            self.n_linked += 1
                            hit = True
                            break
                    if hit:
                        break
        return super().get(url, dest, empty_ok)


# ======================================================================== 站点序列重建（K1）
def read_mapco2_sst(cache, erddap_id):
    """从已缓存的 ERDDAP 年块 CSV 读 SST，返回 {hour: 小时均值}（β 的温度归一，K18）。"""
    acc = {}
    for path in sorted(glob.glob(os.path.join(cache, "erddap", f"{erddap_id}_*.csv.gz"))):
        if path.endswith("_info.csv.gz"):
            continue
        with gzip.open(path, "rt", encoding="utf-8") as f:
            rd = csv.reader(f)
            hdr = next(rd)
            next(rd)
            it, isst = hdr.index("time"), hdr.index("SST")
            for r in rd:
                v = p1._num(r[isst])
                if v is None or not (-2.0 < v < 40.0):
                    continue
                h = p1.iso_to_hour(r[it])
                s = acc.setdefault(h, [0.0, 0])
                s[0] += v
                s[1] += 1
    return {h: s / c for h, (s, c) in acc.items()}


def to_station(name, group, regime, lat, lon, erddap, series, n, h0, s1_depth, sst, raw_events):
    np = _np()
    ser = {k: np.asarray(series[k], dtype=float) for k in ("rain", "wind", "s1", "s5", "sss05", "pco2")}
    t = np.full(n, np.nan)
    for h, v in sst.items():
        i = h - h0
        if 0 <= i < n:
            t[i] = v
    ser["sst"] = t
    return {"name": name, "group": group, "regime": regime, "lat": lat, "lon": lon, "erddap": erddap, "h0": h0, "n": n,
            "series_raw": series, "ser": ser, "s1_depth": s1_depth, "events_rebuild_raw": raw_events}


def rebuild_stations(http, cache, log, names_p1=None):
    """返回站列表（含原始 array 序列供 p1b.run_events、numpy 序列供分析）。"""
    stations = []
    for st in p1.STATIONS:
        if names_p1 is not None and st["name"] not in names_p1:
            continue
        cap = {}
        orig_fe, orig_fm = p1.find_events, p1.fetch_mapco2

        def fm(*a, **k):
            r = orig_fm(*a, **k)
            cap["hs"], cap["he"] = r[0], r[1]
            return r

        def fe(series, n):
            cap["series"], cap["n"] = series, n
            return orig_fe(series, n)

        p1.find_events, p1.fetch_mapco2 = fe, fm
        try:
            r = p1.process_station(http, st, cache, log, False, False)
        finally:
            p1.find_events, p1.fetch_mapco2 = orig_fe, orig_fm
        h0 = cap["hs"] - p1.PAD_H
        la, lo = P1_COORDS[st["name"]]
        stations.append(to_station(st["name"], GROUP_A, st["regime"], la, lo, st["erddap"], cap["series"], cap["n"], h0,
                                   None, read_mapco2_sst(cache, st["erddap"]), r["events"]))
        print(f"[重建] {st['name']}：P1 逻辑合格事件 {len(r['events'])}", flush=True)
    stations_o, _disc, _allds, _spans = p1b.discover_o(http, cache, log)
    for st in stations_o:
        cap = {}
        orig = p1b.run_events

        def re_(series, n, h0, *a, **k):
            cap.update(series=series, n=n, h0=h0, s1_depth=k.get("s1_depth"))
            return orig(series, n, h0, *a, **k)

        p1b.run_events = re_
        try:
            evs, _stats = p1b.process_o_station(http, st, cache, log)
        finally:
            p1b.run_events = orig
        stations.append(to_station(st["name"], GROUP_B, st["regime"], st["lat"], st["lon"], st["erddap"], cap["series"],
                                   cap["n"], cap["h0"], cap["s1_depth"], read_mapco2_sst(cache, st["erddap"]), evs))
        print(f"[重建] O:{st['name']}：合格事件 {len(evs)}", flush=True)
    return stations


class _QuietLog:
    def log(self, *a, **k):
        pass


def events_at(stn, thr, log=None):
    """用 p1 同一事件定义（只改累积阈值）在站序列上求事件（K21）。"""
    old = p1.RAIN_EVENT_MM
    p1.RAIN_EVENT_MM = thr
    try:
        grp = "P1" if stn["group"] == GROUP_A else "O"
        evs, _ = p1b.run_events(stn["series_raw"], stn["n"], stn["h0"], stn["name"], stn["regime"], stn["lon"], grp,
                                log or _QuietLog(), s1_depth=stn["s1_depth"])
    finally:
        p1.RAIN_EVENT_MM = old
    for e in evs:
        e["i"] = e["hour"] - stn["h0"]
    return evs


def read_event_csv(path, station_col="station"):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def validate_rebuild(events10, p1_csv, p1b_csv):
    """K1：重建的 10 mm 事件须与 P1/P1b 产物逐条一致。返回记录；不一致抛 RuntimeError。"""
    ref = {}
    for path in (p1_csv, p1b_csv):
        for r in read_event_csv(path):
            ref[(r["station"], r["onset_utc"])] = r
    got = {(e["station"], e["onset_utc"]): e for e in events10}
    miss = sorted(set(ref) - set(got))
    extra = sorted(set(got) - set(ref))
    bad = []
    for k in set(ref) & set(got):
        r, e = ref[k], got[k]
        for col, key in (("dS1_0_6h", "ds1"), ("dS05_0_6h", "ds05")):
            a = p1b._f(r[col]) if r[col] != "" else NAN
            b = e[key]
            if (a == a) != (b == b) or (a == a and abs(a - round(b, 4)) > 6e-5):
                bad.append((k, col, r[col], b))
        if int(r["n_ctrl"]) != e["n_ctrl"]:
            bad.append((k, "n_ctrl", r["n_ctrl"], e["n_ctrl"]))
    n_d4 = sum(1 for e in events10 if e["ds05"] == e["ds05"] and e["ds05"] <= p1.D4_DS_THRESH)
    rec = {"n_ref": len(ref), "n_rebuilt": len(got), "missing_in_rebuild": [list(x) for x in miss[:20]],
           "extra_in_rebuild": [list(x) for x in extra[:20]], "value_mismatch": [list(map(str, x)) for x in bad[:20]],
           "n_missing": len(miss), "n_extra": len(extra), "n_value_mismatch": len(bad), "d4_count": n_d4}
    if miss or extra or bad or len(got) != EXPECTED_MERGED or n_d4 != EXPECTED_D4:
        raise RuntimeError(f"事件重建与 P1/P1b 产物不一致（K1）：{json.dumps(rec, ensure_ascii=False)[:1500]}")
    return rec


# ======================================================================== 小时集合（K22）
def needed_hours(events):
    s = set()
    for e in events:
        s.update(range(e["hour"] + WIN_LO_H, e["hour"] + WIN_HI_H))
    return s


def hour_url(h):
    d = datetime.fromtimestamp(h * 3600, tz=timezone.utc)
    return (f"/data/cmorph-high-resolution-global-precipitation-estimates/access/30min/8km/"
            f"{d:%Y/%m/%d}/CMORPH_V1.0_ADJ_8km-30min_{d:%Y%m%d%H}.nc")


# ======================================================================== CMORPH 像元缓存与下载（K23）
class PixelCache:
    """cmorph_pixels.csv：每行 hour,station,v0,v1,m0,m1,n0,n1；每小时块以 hour,__DONE__ 或 hour,__MISSING__,原因 结束。"""

    def __init__(self, path):
        self.path = path
        self.data = {}        # hour -> {station: (v0, v1, m0, m1, n0, n1)}
        self.missing = {}     # hour -> reason

    @staticmethod
    def parse(path):
        data, missing, pend = {}, {}, {}
        if not os.path.exists(path):
            return data, missing
        with open(path, encoding="utf-8") as f:
            for line in f:
                p = line.rstrip("\n").split(",")
                if len(p) < 2:
                    continue
                try:
                    h = int(p[0])
                except ValueError:
                    continue
                if p[1] == "__DONE__":
                    data[h] = pend.pop(h, {})
                elif p[1] == "__MISSING__":
                    missing[h] = p[2] if len(p) > 2 else ""
                    pend.pop(h, None)
                elif len(p) == 8:
                    try:
                        vals = tuple(float(x) for x in p[2:6]) + (int(p[6]), int(p[7]))
                    except ValueError:
                        continue
                    pend.setdefault(h, {})[p[1]] = vals
        return data, missing

    def load(self, *paths):
        for pth in paths:
            if pth and os.path.exists(pth):
                d, m = self.parse(pth)
                self.data.update(d)
                for h, why in m.items():
                    if h not in self.data:
                        self.missing[h] = why
        tmp = self.path + ".part"
        with open(tmp, "w", encoding="utf-8") as f:
            for h in sorted(self.data):
                f.write(self._block(h, self.data[h]))
            for h in sorted(self.missing):
                f.write(f"{h},__MISSING__,{self.missing[h]}\n")
        os.replace(tmp, self.path)

    @staticmethod
    def _fmt(x):
        return "nan" if x != x else f"{x:.3f}"

    def _block(self, h, rows):
        out = [f"{h},{st},{self._fmt(v[0])},{self._fmt(v[1])},{self._fmt(v[2])},{self._fmt(v[3])},{v[4]},{v[5]}\n"
               for st, v in sorted(rows.items())]
        out.append(f"{h},__DONE__\n")
        return "".join(out)

    def _append(self, text):
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())

    def add(self, h, rows):
        self._append(self._block(h, rows))
        self.data[h] = rows

    def add_missing(self, h, why):
        self._append(f"{h},__MISSING__,{why}\n")
        self.missing[h] = why

    def get(self, h, st):
        r = self.data.get(h)
        if r is None:
            return None
        return r.get(st)


class Fetcher:
    """线程内持久 HTTPS 连接；下载到临时文件并校验长度。"""

    def __init__(self, tmpdir, log):
        self.tmpdir = tmpdir
        self.log = log
        self.local = threading.local()
        self.bytes = 0
        self.n_req = 0
        self.lock = threading.Lock()

    def _conn(self, fresh=False):
        c = getattr(self.local, "c", None)
        if fresh and c is not None:
            try:
                c.close()
            except Exception:
                pass
            c = None
        if c is None:
            c = http.client.HTTPSConnection(CMORPH_HOST, timeout=180)
            self.local.c = c
        return c

    def fetch(self, h):
        """返回 ("ok", path) 或 ("missing", 原因)；重试用尽抛 DataSourceError。"""
        path = hour_url(h)
        dest = os.path.join(self.tmpdir, f"cmorph_{h}.nc")
        last, n404 = None, 0
        for attempt in range(1, MAX_ATTEMPTS + 1):
            try:
                c = self._conn(fresh=attempt > 1)
                c.request("GET", path, headers={"User-Agent": UA, "Connection": "keep-alive"})
                r = c.getresponse()
                with self.lock:
                    self.n_req += 1
                if r.status == 200:
                    clen = int(r.headers.get("Content-Length") or -1)
                    got = 0
                    with open(dest, "wb") as f:
                        while True:
                            b = r.read(1 << 20)
                            if not b:
                                break
                            f.write(b)
                            got += len(b)
                    if clen >= 0 and got != clen:
                        raise IOError(f"长度 {got} ≠ {clen}")
                    with self.lock:
                        self.bytes += got
                    return "ok", dest
                r.read()
                if r.status == 404:
                    n404 += 1
                    if n404 >= 2:
                        return "missing", "HTTP404"
                    last = "HTTP 404"
                else:
                    last = f"HTTP {r.status}"
                    if r.status not in (408, 429, 500, 502, 503, 504):
                        raise DataSourceError(f"{path} → {last}")
            except DataSourceError:
                raise
            except Exception as e:  # 连接/读写错误：换新连接重试
                last = f"{type(e).__name__}: {e}"
                self._conn(fresh=True)
            if os.path.exists(dest):
                os.remove(dest)
            if attempt < MAX_ATTEMPTS:
                time.sleep(min(40, 5 * 2 ** (attempt - 1)))
        raise DataSourceError(f"{path} → {MAX_ATTEMPTS} 次失败：{last}")


class Extractor:
    """主线程：netCDF4 打开单个小时文件，抽 13 站中心像元与 3×3（K7/K8）。"""

    def __init__(self, stations):
        self.st = [(s["name"], s["lat"], s["lon"] % 360.0) for s in stations]
        self.idx = None
        self.shape = None

    def _index(self, ds):
        np = _np()
        lat = np.asarray(ds["lat"][:], dtype=float)
        lon = np.asarray(ds["lon"][:], dtype=float)
        self.shape = (len(lat), len(lon))
        self.idx = {}
        self.px = {}
        for name, la, lo in self.st:
            i = int(np.argmin(np.abs(lat - la)))
            dl = np.abs(((lon - lo) + 180.0) % 360.0 - 180.0)
            j = int(np.argmin(dl))
            i = min(max(i, 1), len(lat) - 2)
            js = [(j - 1) % len(lon), j, (j + 1) % len(lon)]
            self.idx[name] = (i, js)
            self.px[name] = {"lat": float(lat[i]), "lon": float(lon[j]),
                             "dist_km": round(p1b.haversine_km(la, lo if lo <= 180 else lo - 360, float(lat[i]),
                                                               float(lon[j]) if lon[j] <= 180 else float(lon[j]) - 360), 2)}

    def extract(self, path, h):
        import netCDF4
        np = _np()
        with netCDF4.Dataset(path) as ds:
            t = [int(x) for x in ds["time"][:]]
            if t != [h * 3600, h * 3600 + 1800]:
                return None, f"time_mismatch:{t}"
            shp = (ds.dimensions["lat"].size, ds.dimensions["lon"].size)
            if self.idx is None or shp != self.shape:
                self._index(ds)
            v = ds["cmorph"]
            v.set_auto_maskandscale(True)
            rows = {}
            for name, (i, js) in self.idx.items():
                blk = np.ma.filled(np.ma.asarray(v[:, i - 1:i + 2, js]).astype(float), np.nan)  # (2,3,3)
                c0, c1 = float(blk[0, 1, 1]), float(blk[1, 1, 1])
                n0 = int(np.isfinite(blk[0]).sum())
                n1 = int(np.isfinite(blk[1]).sum())
                m0 = float(np.nanmean(blk[0])) if n0 else NAN
                m1 = float(np.nanmean(blk[1])) if n1 else NAN
                rows[name] = (c0, c1, m0, m1, n0, n1)
        return rows, None


def run_downloads(hours_ordered, cache, stations, out_dir, log, max_conc):
    todo = [h for h in hours_ordered if h not in cache.data and h not in cache.missing]
    tmpdir = os.path.join(out_dir, "tmp_cmorph")
    os.makedirs(tmpdir, exist_ok=True)
    fx = _rio.cmorph_fetcher(DataSourceError, log)  # [repro] 原 Fetcher 下到 tmp_cmorph 后删；改经 raw_store.CmorphRaw 落 raw 档、先查档后下载
    ex_ = Extractor(stations)
    t0 = time.monotonic()
    n_done = 0
    failed = []
    consec = 0
    print(f"[CMORPH] 需 {len(hours_ordered)} 小时，已缓存 {len(hours_ordered) - len(todo)}，待下 {len(todo)}", flush=True)

    def handle(h, res):
        nonlocal n_done
        kind, val = res
        if kind == "missing":
            cache.add_missing(h, val)
            return
        try:
            rows, why = ex_.extract(val, h)
        except Exception as e:
            rows, why = None, f"read_error:{type(e).__name__}"
        finally:
            pass  # [repro] 原件留在 raw 档，不删；原为 os.remove(val)
        if rows is None:
            cache.add_missing(h, why)
            log.log(f"CMORPH {h} 记缺测：{why}")
        else:
            cache.add(h, rows)
        n_done += 1
        if n_done % 250 == 0:
            el = time.monotonic() - t0
            rate = fx.bytes / 1e6 / max(el, 1e-6)
            left = (len(todo) - n_done) * CMORPH_MEAN_MB / max(rate, 1e-6)
            msg = f"[CMORPH] {n_done}/{len(todo)}，{fx.bytes / 1e9:.2f} GB，{rate:.2f} MB/s，预计剩余 {left / 3600:.2f} h"
            print(msg, flush=True)
            log.log(msg)

    it = iter(todo)
    with cf.ThreadPoolExecutor(max_workers=max_conc) as ex:
        pend = {}
        for h in it:
            pend[ex.submit(fx.fetch, h)] = h
            if len(pend) >= max_conc * 2:
                break
        while pend:
            done, _ = cf.wait(pend, return_when=cf.FIRST_COMPLETED)
            for fut in done:
                h = pend.pop(fut)
                try:
                    res = fut.result()
                except DataSourceError as e:
                    failed.append(h)
                    consec += 1
                    log.log(f"CMORPH 失败 {h}：{e}")
                    if consec >= MAX_CONSEC_FAIL:
                        for f in pend:
                            f.cancel()
                        raise DataSourceError(f"CMORPH 连续 {consec} 个文件失败，停止（可 --resume-from 续跑）")
                    continue
                consec = 0
                handle(h, res)
            for h in it:
                pend[ex.submit(fx.fetch, h)] = h
                if len(pend) >= max_conc * 2:
                    break
    still = []
    for h in failed:  # 末轮串行补试
        try:
            handle(h, fx.fetch(h))
        except DataSourceError as e:
            still.append(h)
            log.log(f"CMORPH 末轮仍失败 {h}：{e}")
    shutil.rmtree(tmpdir, ignore_errors=True)
    if still:
        raise DataSourceError(f"CMORPH 末轮仍失败 {len(still)} 小时（例 {still[:5]}），可 --resume-from 续跑")
    return {"hours_needed": len(hours_ordered), "downloaded_this_run": n_done, "bytes_this_run": fx.bytes,
            "http_requests": fx.n_req, "seconds": round(time.monotonic() - t0, 1),
            "pixels": getattr(ex_, "px", None)}


# ======================================================================== RIM-3 / S20（K2–K6）
_D0_INTERP = None


def d0_interp(u, r):
    global _D0_INTERP
    np = _np()
    if _D0_INTERP is None:
        from scipy.interpolate import RegularGridInterpolator
        _D0_INTERP = RegularGridInterpolator((np.array(U_GRID, float), np.array(R_GRID, float)), np.array(DL, float),
                                             method="linear", bounds_error=False, fill_value=np.nan)
    u = np.asarray(u, float)
    r = np.asarray(r, float)
    pts = np.stack([u.ravel(), r.ravel()], axis=-1)
    out = _D0_INTERP(pts).reshape(u.shape)
    bad = ~(np.isfinite(u) & np.isfinite(r))
    out[bad] = np.nan
    return out


def rim_factor(P, U, z, t_cur_depth=T_CURRENT_S):
    """P、U：半步序列（mm/h、m/s，长 L≥49）。返回长 L−48 的稀释因子 F(z)，对应半步 48..L−1（K3）。
    t_cur_depth：只改当前项深度衰减因子 exp(−z²/(4·Kz·t)) 里的 t；量级项恒为 c2·IRR·1800/√(Kz·1 s)。
    默认＝T_CURRENT_S（主结果路径）；K27 探索性变体传 1800。"""
    np = _np()
    from numpy.lib.stride_tricks import sliding_window_view as swv
    P = np.asarray(P, float)
    U = np.asarray(U, float)
    L = len(P)
    tc = np.array(TC, float)
    ti = np.array(TI, float)
    Pw = swv(P, 48)[:L - 48]            # 行 k：半步 k..k+47（＝当前半步 48+k 的前 48 个）
    Uw = swv(U, 48)[:L - 48]
    irr = Pw / 1000.0 / 3600.0
    kz = KZ_COEF * Uw ** 2
    d0 = d0_interp(Uw, Pw)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        term = C1 * irr * tc / np.sqrt(kz * ti)
        if z:
            term = term * np.exp(-(z ** 2) / (4.0 * kz * ti))
        prior = d0 / (d0 + term)
        Pc, Uc = P[48:], U[48:]
        irrc = Pc / 1000.0 / 3600.0
        kzc = KZ_COEF * Uc ** 2
        d0c = d0_interp(Uc, Pc)
        cterm = C2 * irrc * 1800.0 / np.sqrt(kzc * T_CURRENT_S)
        if z:
            cterm = cterm * np.exp(-(z ** 2) / (4.0 * kzc * t_cur_depth))
        cur = d0c / (d0c + cterm)
    F = np.prod(prior, axis=-1) * cur
    return F


def s20_delta(P, U):
    np = _np()
    with np.errstate(divide="ignore", invalid="ignore"):
        return S20_A * np.asarray(P, float) * np.asarray(U, float) ** (-S20_B)


def fill_wind(seg, max_gap):
    """线性插补 ≤max_gap 小时的内部缺口（两端有值），返回新数组（K9）。"""
    np = _np()
    w = np.array(seg, float)
    ok = np.isfinite(w)
    idx = np.nonzero(ok)[0]
    if len(idx) < 2:
        return w
    for a, b in zip(idx[:-1], idx[1:]):
        gap = b - a - 1
        if 0 < gap <= max_gap:
            w[a + 1:b] = w[a] + (w[b] - w[a]) * (np.arange(1, gap + 1) / (gap + 1))
    return w


def event_forcing(stn, H, variant, pix):
    """返回 (P 半步, U 半步, 诊断)；[H−30, H+12) 共 42 小时＝84 半步。"""
    np = _np()
    i0 = H + WIN_LO_H - stn["h0"]
    nh = WIN_HI_H - WIN_LO_H
    diag = {"wind_floor": 0, "gauge_filled": 0}
    pad = WIND_GAP_MAX_H
    ws = stn["ser"]["wind"]
    lo, hi = i0 - pad, i0 + nh + pad
    seg = np.full(hi - lo, np.nan)
    a, b = max(lo, 0), min(hi, len(ws))
    if b > a:
        seg[a - lo:b - lo] = ws[a:b]
    w = fill_wind(seg, WIND_GAP_MAX_H)[pad:pad + nh] * _rp.wind_factor(stn["name"], WIND_FACTOR, WIND_Z0_M)  # [repro] 风高开关 [options] wind_height；缺省 uniform_4m 时 _rp.wind_factor 原样返回原系数
    if not np.all(np.isfinite(w)):
        return None, None, dict(diag, reason="wind_gap")
    diag["wind_floor"] = int((w < WIND_FLOOR).sum())
    w = np.maximum(w, WIND_FLOOR)
    U = np.repeat(w, 2)
    if variant == "cmorph":
        P = np.empty(2 * nh)
        for k in range(nh):
            r = pix.get(H + WIN_LO_H + k, stn["name"])
            if r is None:
                return None, None, dict(diag, reason="cmorph_missing")
            P[2 * k], P[2 * k + 1] = r[0], r[1]
        if not np.all(np.isfinite(P)):
            return None, None, dict(diag, reason="cmorph_missing")
        diag["p_gt_d0_range"] = int((P > R_GRID[-1]).sum())      # d0 表外 → 模型 NaN（K2），只计数
    else:
        rs = stn["ser"]["rain"]
        g = np.full(nh, np.nan)
        a, b = max(i0, 0), min(i0 + nh, len(rs))
        g[a - i0:b - i0] = rs[a:b]
        diag["gauge_filled"] = int((~np.isfinite(g)).sum())
        diag["gauge_neg"] = int((g < 0).sum())                   # K21：负值置 0 也计数
        g = np.where(np.isfinite(g), np.maximum(g, 0.0), 0.0)
        P = np.repeat(g, 2)
    return P, U, dict(diag, reason=None)


def model_hourly(P, U, s0_for_s20, zmodel):
    """返回 {层: F 小时值(18)} 与 S20 ΔS 小时值(18)，对应小时 t0−6 … t0+11（K5）。
    zmodel：{层标签: 模型深度}；「1 m」层标签 1.0 的模型深度＝事件实际观测深度 s1_depth_m（K28）。"""
    np = _np()
    out = {}
    for z in DEPTHS:
        F = rim_factor(P, U, zmodel[z])          # 半步 48..83 → 36 个
        out[z] = F.reshape(-1, 2).mean(axis=1)
    ds = s20_delta(P[48:], U[48:])
    if isnum(s0_for_s20):
        ds = np.maximum(ds, -s0_for_s20)         # Witte：S20 结果 <0 置 0
    return out, ds.reshape(-1, 2).mean(axis=1)


def win_delta(hourly, win):
    """hourly：t0−6 … 的小时序列；ΔX＝mean[t0+a, t0+b) − median[t0−6, t0)。"""
    np = _np()
    pre = hourly[:p1.PRE_H]
    post = hourly[p1.PRE_H + win[0]:p1.PRE_H + win[1]]
    if not (np.all(np.isfinite(pre)) and np.all(np.isfinite(post))):
        return NAN
    return float(np.mean(post) - np.median(pre))


# ======================================================================== 观测（I8 同法）
def obs_delta(arr, i, win):
    np = _np()
    pre = arr[max(0, i - p1.PRE_H):i]
    post = arr[i + win[0]:i + win[1]]
    pre = pre[np.isfinite(pre)]
    post = post[np.isfinite(post)]
    if not len(pre) or not len(post):
        return NAN
    return float(post.mean() - np.median(pre))


def pre_median(arr, i):
    np = _np()
    pre = arr[max(0, i - p1.PRE_H):i]
    pre = pre[np.isfinite(pre)]
    return float(np.median(pre)) if len(pre) else NAN


def dpco2_norm(stn, i, win=PRIMARY_WIN):
    """K18：温度归一后的 ΔpCO₂。"""
    np = _np()
    pc, t = stn["ser"]["pco2"], stn["ser"]["sst"]
    tref = pre_median(t, i)
    if not isnum(tref):
        return NAN

    def norm(sl):
        a, b = pc[sl], t[sl]
        ok = np.isfinite(a) & np.isfinite(b)
        return a[ok] * np.exp(T_COEF * (tref - b[ok]))

    pre = norm(slice(max(0, i - p1.PRE_H), i))
    post = norm(slice(i + win[0], i + win[1]))
    if not len(pre) or not len(post):
        return NAN
    return float(post.mean() - np.median(pre))


def model_depths(e, station=""):
    """K28：层标签 → 模型深度。「1 m」层＝事件实际观测深度 s1_depth_m；缺失或越界即报错，不回退 1.0。"""
    z1 = e.get("s1_depth_m")
    if not isnum(z1) or not (0.0 < float(z1) <= K28_MAX_DEPTH + 1e-9):
        raise ValueError(f"K28：事件 {station} {e.get('onset_utc')} 的 s1_depth_m 缺失或越界（{z1!r}），不回退 1.0")
    return {0.5: 0.5, 1.0: float(z1), 5.0: 5.0}


def event_record(stn, e, pix, variants=("cmorph",), windows=(PRIMARY_WIN,), k27=False):
    """一个事件的观测与模型 ΔS（各窗口、各驱动）。k27：另算探索性变体（K27，只 CMORPH、0–6 h、0.5/1 m）。
    模型「1 m」层用事件实际观测深度（K28）；键仍以层标签 1.0 记。"""
    ser = stn["ser"]
    i = e["i"]
    zm = model_depths(e, stn["name"])
    rec = {"station": stn["name"], "group": stn["group"], "season": e["season"], "hour": e["hour"],
           "onset_utc": e["onset_utc"], "acc24": e["acc24"], "i": i, "s1_substitute": bool(e.get("s1_substitute")),
           "s1_depth_m": zm[1.0], "zmodel": zm,
           "n_ctrl": e["n_ctrl"], "ctrl_hours": e["ctrl_hours"], "ctrl_ds1": list(e["ctrl_ds1"]),
           "s0": {0.5: pre_median(ser["sss05"], i), 1.0: pre_median(ser["s1"], i), 5.0: pre_median(ser["s5"], i)},
           "obs": {}, "model": {}, "reason": {}, "diag": {}, "windows": list(windows)}
    for v in variants:
        P, U, diag = event_forcing(stn, e["hour"], v, pix)
        rec["reason"][v] = diag["reason"]
        rec["diag"][v] = diag
        if P is None:
            for w in windows:
                for z in DEPTHS:
                    rec["model"][(v, "rim", z, w)] = NAN
                rec["model"][(v, "s20", 0.0, w)] = NAN
            if v == "cmorph" and k27:
                for z in K27_DEPTHS:
                    rec["model"][(K27_TAG, "rim", z, PRIMARY_WIN)] = NAN
            continue
        Fh, s20h = model_hourly(P, U, rec["s0"][0.5], zm)
        for w in windows:
            for z in DEPTHS:
                s0 = rec["s0"][z]
                d = win_delta(Fh[z], w)
                rec["model"][(v, "rim", z, w)] = s0 * d if isnum(s0) and isnum(d) else NAN
            rec["model"][(v, "s20", 0.0, w)] = win_delta(s20h, w)
        if v == "cmorph" and k27:                  # K27：独立键，D3/出口/主比值都不读它
            for z in K27_DEPTHS:
                s0 = rec["s0"][z]
                d = win_delta(rim_factor(P, U, zm[z], t_cur_depth=K27_T_CUR_DEPTH_S).reshape(-1, 2).mean(axis=1),
                              PRIMARY_WIN)
                rec["model"][(K27_TAG, "rim", z, PRIMARY_WIN)] = s0 * d if isnum(s0) and isnum(d) else NAN
    return rec


def fill_obs(rec, stn):
    """事件观测 ΔS(z) 与 ΔpCO₂——只在 D3 写盘之后调用（K11）。"""
    key = {0.5: "sss05", 1.0: "s1", 5.0: "s5"}
    for w in rec["windows"]:
        for z in DEPTHS:
            rec["obs"][(z, w)] = obs_delta(stn["ser"][key[z]], rec["i"], w)
    rec["dpco2"] = dpco2_norm(stn, rec["i"])


# ======================================================================== 统计（K12–K15）
def cluster_boot(keys, pairs, B=BOOT_B, seed=BOOT_SEED):
    """keys：每事件簇键；pairs：{名: (num 数组, den 数组)}（NaN 事件须已剔除）。返回 {名: 比值 bootstrap 数组}。"""
    np = _np()
    uk = sorted(set(keys))
    kid = {k: j for j, k in enumerate(uk)}
    K = len(uk)
    idx = np.array([kid[k] for k in keys], int)
    sums = {}
    for name, (num, den) in pairs.items():
        sums[name] = (np.bincount(idx, weights=np.asarray(num, float), minlength=K),
                      np.bincount(idx, weights=np.asarray(den, float), minlength=K))
    rng = np.random.Generator(np.random.PCG64(seed))
    draw = rng.integers(0, K, size=(B, K))
    out = {}
    for name, (sn, sd) in sums.items():
        with np.errstate(divide="ignore", invalid="ignore"):
            out[name] = sn[draw].sum(axis=1) / sd[draw].sum(axis=1)
    return out, K


def ratio_stats(recs, num_f, den_f, label, extra_pairs=None, tost=False):
    """recs 上的 Σnum/Σden，含簇 bootstrap CI；extra_pairs 共用抽样（P5 模型/差值）。"""
    np = _np()
    rows = []
    for r in recs:
        a, b = num_f(r), den_f(r)
        ext = [(f(r), g(r)) for _, f, g in (extra_pairs or [])]
        if isnum(a) and isnum(b) and all(isnum(x) and isnum(y) for x, y in ext):
            rows.append((r, a, b, ext))
    res = {"label": label, "n_events": len(rows)}
    if len(rows) < 2:
        res.update({"evaluable": False, "reason": "事件不足"})
        return res
    keys = [(r["station"], r["season"]) for r, *_ in rows]
    pairs = {"main": ([x[1] for x in rows], [x[2] for x in rows])}
    for j, (nm, _f, _g) in enumerate(extra_pairs or []):
        pairs[nm] = ([x[3][j][0] for x in rows], [x[3][j][1] for x in rows])
    boots, K = cluster_boot(keys, pairs)
    num, den = sum(x[1] for x in rows), sum(x[2] for x in rows)
    point = num / den if den else NAN
    bm = boots["main"]
    bm = bm[np.isfinite(bm)]
    res.update({"evaluable": K >= 2 and len(bm) > 100, "n_clusters": K, "point": rnd(point),
                "sum_num": rnd(num, 5), "sum_den": rnd(den, 5), "n_boot_finite": int(len(bm))})
    if not res["evaluable"]:
        res["reason"] = "簇数 <2 或 bootstrap 有效样本不足"
        return res
    res["ci95"] = [rnd(np.percentile(bm, 2.5)), rnd(np.percentile(bm, 97.5))]
    res["ci90"] = [rnd(np.percentile(bm, 5)), rnd(np.percentile(bm, 95))]
    res["se_boot"] = rnd(float(np.std(bm, ddof=1)))
    if tost:
        pl = float(np.mean(bm <= SESOI_LO))
        pu = float(np.mean(bm >= SESOI_HI))
        res["tost"] = {"alpha": 0.05, "bounds": [SESOI_LO, SESOI_HI], "p_lower": rnd(pl), "p_upper": rnd(pu),
                       "equivalent": pl < 0.05 and pu < 0.05}
        res["mde_ratio_units"] = rnd(MDE_Z * float(np.std(bm, ddof=1)))
    lo, hi = res["ci95"]
    res["ci95_outside_sesoi"] = bool(hi < SESOI_LO or lo > SESOI_HI)
    res["ci95_inside_sesoi"] = bool(lo >= SESOI_LO and hi <= SESOI_HI)
    res["ci95_excludes_1_info"] = bool(hi < 1.0 or lo > 1.0)
    for j, (nm, _f, _g) in enumerate(extra_pairs or []):
        b2 = boots[nm]
        n2 = sum(x[3][j][0] for x in rows)
        d2 = sum(x[3][j][1] for x in rows)
        ok = np.isfinite(boots["main"]) & np.isfinite(b2)
        diff = boots["main"] - b2
        diff = diff[np.isfinite(diff)]
        res[nm] = {"point": rnd(n2 / d2 if d2 else NAN),
                   "ci95": [rnd(np.percentile(b2[np.isfinite(b2)], 2.5)), rnd(np.percentile(b2[np.isfinite(b2)], 97.5))],
                   "diff_main_minus_this_ci95": [rnd(np.percentile(diff, 2.5)), rnd(np.percentile(diff, 97.5))]
                   if len(diff) else None, "n_boot_paired": int(ok.sum())}
    return res


def d3_gate(recs, stations_by_name):
    """K11：只用模型预测与对照窗观测（不用事件观测值）。"""
    np = _np()
    out = {"rule": f"|mean ΔS_RIM(z) − mean ΔS_S20| ≥ {D3_K_SE}·SE_obs(z)，z=0.5 与 1 m 都满足才过（K11）"}
    passes = []
    for z, key in ((0.5, "sss05"), (1.0, "s1")):
        use = [r for r in recs if isnum(r["model"][("cmorph", "rim", z, PRIMARY_WIN)])
               and isnum(r["model"][("cmorph", "s20", 0.0, PRIMARY_WIN)])
               and _obs_available(stations_by_name[r["station"]], key, r["i"])]
        ctrl = []
        for r in use:
            arr = stations_by_name[r["station"]]["ser"][key]
            h0 = stations_by_name[r["station"]]["h0"]
            for ch in r["ctrl_hours"]:
                v = obs_delta(arr, ch - h0, PRIMARY_WIN)
                if isnum(v):
                    ctrl.append(v)
        n = len(use)
        rim = [r["model"][("cmorph", "rim", z, PRIMARY_WIN)] for r in use]
        s20 = [r["model"][("cmorph", "s20", 0.0, PRIMARY_WIN)] for r in use]
        d = {"n_events": n, "n_control_windows": len(ctrl)}
        if n < 2 or len(ctrl) < 2:
            d.update({"evaluable": False, "pass": False})
            passes.append(False)
        else:
            se = float(np.std(ctrl, ddof=1)) / math.sqrt(n)
            diff = float(np.mean(rim) - np.mean(s20))
            d.update({"evaluable": True, "mean_dS_rim": rnd(np.mean(rim)), "mean_dS_s20": rnd(np.mean(s20)),
                      "diff": rnd(diff), "sd_control": rnd(np.std(ctrl, ddof=1)), "se_obs": rnd(se),
                      "threshold_2se": rnd(D3_K_SE * se), "pass": abs(diff) >= D3_K_SE * se})
            passes.append(d["pass"])
        out[f"z={z}m"] = d
    out["pass"] = all(passes)
    out["verdict"] = "过（可做模型裁决）" if out["pass"] else "不确定 → 出口 3"
    return out


def _obs_available(stn, key, i):
    """观测在雨前 6 h 与 0–6 h 都有有效值（不读取数值本身，只判可用性）。"""
    np = _np()
    a = stn["ser"][key]
    return bool(np.isfinite(a[max(0, i - p1.PRE_H):i]).any() and np.isfinite(a[i:i + PRIMARY_WIN[1]]).any())


def decide_exit(gates_p1b, d3, r_rim, r_s20):
    """K15。"""
    g = {k: bool(gates_p1b[k]) for k in ("D1", "D2", "D4", "D5")}
    ex3 = [k for k in ("D1", "D2", "D5") if not g[k]] + ([] if d3["pass"] else ["D3"])
    res = {"gates": dict(g, D3=bool(d3["pass"])), "exit3_triggers": ex3}
    if ex3:
        res.update({"exit": 3, "label": "出口 3：数据或信号死（判死并回填）", "also_satisfied": []})
        return res
    c1_rim = bool(r_rim.get("evaluable") and r_rim.get("ci95_outside_sesoi"))
    c1_s20 = bool(r_s20.get("evaluable") and r_s20.get("ci95_outside_sesoi"))
    c2 = bool(r_rim.get("evaluable") and r_rim.get("ci95_inside_sesoi") and r_rim.get("tost", {}).get("equivalent"))
    res["conditions"] = {"exit1_R_RIM_ci_outside_sesoi": c1_rim, "exit1_R_S20_excluded": c1_s20,
                         "exit2_R_RIM_ci_inside_and_tost": c2}
    sat = ([1] if (c1_rim or c1_s20) else []) + ([2] if c2 else [])
    res["also_satisfied"] = sat
    if 1 in sat:
        res.update({"exit": 1, "label": "出口 1：信号显著 → P3 全球重算"})
        if 2 in sat:
            res["note"] = "出口 1 与出口 2 条件同时成立（K15），主标出口 1"
    elif 2 in sat:
        res.update({"exit": 2, "label": "出口 2：有信息的 null（RIM-3 在 0.5–5 m 获独立原位支持，不得称海面值已证实）"})
    else:
        res.update({"exit": None, "label": "未落入出口：R_RIM 95% CI 跨 SESOI 边界且 R_S20 未被排除（K15）"})
    return res


# ======================================================================== GLODAP 与 β_model（K19/K24）
def glodap_values(stations, out_dir, log, resume_dir=None):
    path = os.path.join(out_dir, "p2_glodap.json")
    for src in ([resume_dir] if resume_dir else []) + [out_dir]:
        p = os.path.join(src, "p2_glodap.json")
        if os.path.exists(p):
            with open(p, encoding="utf-8") as f:
                got = json.load(f)
            if set(got.get("stations", {})) >= {s["name"] for s in stations}:
                if p != path:
                    shutil.copyfile(p, path)
                return got
    tmp = os.path.join(out_dir, "tmp_glodap")
    os.makedirs(tmp, exist_ok=True)
    want = {v: k for k, v in GLODAP_VARS.items()}
    last = None
    for attempt in range(1, 4):
        try:
            found = {}
            _tar = _rio.raw_file(GLODAP_URL, "p3-cache/glodap/GLODAPv2.2016b.MappedProduct.tar.gz", "glodap_tar",
                                 DataSourceError)  # [repro] 原为网络流式解包、不落原件；改为 tar 原件先落 raw 档
            with open(_tar, "rb") as r:  # [repro] 从本地原件流式解包，成员选取与读法不变
                with tarfile.open(fileobj=r, mode="r|gz") as tf:
                    for m in tf:
                        base = os.path.basename(m.name)
                        if base in want and m.isfile():
                            dst = os.path.join(tmp, base)
                            with tf.extractfile(m) as src, open(dst, "wb") as f:
                                shutil.copyfileobj(src, f, 1 << 20)
                            found[want[base]] = dst
                            if len(found) == len(want):
                                break
            if len(found) != len(want):
                raise DataSourceError(f"GLODAP tar 缺成员：{sorted(set(want.values()) - set(found))}")
            break
        except DataSourceError:
            raise
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
            log.log(f"GLODAP 第 {attempt} 次失败：{last}")
            time.sleep(10 * attempt)
    else:
        raise DataSourceError(f"GLODAP 下载失败：{last}")
    vals = read_glodap_points(found, stations)
    shutil.rmtree(tmp, ignore_errors=True)
    got = {"source": GLODAP_URL, "rule": "depth_surface 第 0 层、最近 1° 格点（经度按 mod 360 环绕），NaN 不回填（K19）",
           "stations": vals}
    jdump(got, path)
    return got


def read_glodap_points(files, stations):
    import netCDF4
    np = _np()
    out = {s["name"]: {} for s in stations}
    for var, path in files.items():
        with netCDF4.Dataset(path) as ds:
            v = ds[var]
            dims = v.dimensions
            lat = np.asarray(ds["lat"][:], float)
            lon = np.asarray(ds["lon"][:], float)
            ddim = [d for d in dims if "depth" in d.lower()][0]
            dco = np.asarray(ds[ddim][:], float) if ddim in ds.variables else np.arange(ds.dimensions[ddim].size)
            k0 = int(np.nonzero(dco == 0)[0][0]) if np.any(dco == 0) else 0
            for s in stations:
                i = int(np.argmin(np.abs(lat - s["lat"])))
                j = int(np.argmin(np.abs(((lon - s["lon"]) + 180.0) % 360.0 - 180.0)))
                ix = []
                for d in dims:
                    if d == ddim:
                        ix.append(k0)
                    elif "lat" in d.lower():
                        ix.append(i)
                    elif "lon" in d.lower():
                        ix.append(j)
                    else:
                        ix.append(0)
                x = np.ma.filled(np.ma.asarray(v[tuple(ix)]).astype(float), np.nan)
                out[s["name"]][var] = rnd(float(x), 4)
                out[s["name"]]["cell"] = [rnd(lat[i], 3), rnd(lon[j], 3)]
    return out


def beta_model(ta, dic, sst, sss):
    """照 Witte calculate_fco2_dilution_per_psu；返回 (fCO2 斜率, pCO2 斜率) [µatm/psu]。"""
    import PyCO2SYS as pyco2
    from scipy import stats
    np = _np()
    sals = np.array([sss - k for k in range(10)]) if sss > 9 else np.arange(sss, 0, -1)
    sr = sals / sss
    res = pyco2.sys(par1=ta * sr, par2=dic * sr + DIC_RAIN * (1 - sr), par1_type=1, par2_type=2,
                    salinity=sals, temperature=sst, opt_k_carbonic=K_CARBONIC)
    return (float(stats.linregress(x=sals, y=res["fCO2"]).slope), float(stats.linregress(x=sals, y=res["pCO2"]).slope))


# ======================================================================== 分析汇总
def g_rim(z, v="cmorph", w=PRIMARY_WIN):
    return lambda r: r["model"].get((v, "rim", z, w), NAN)


def g_s20(v="cmorph", w=PRIMARY_WIN):
    return lambda r: r["model"].get((v, "s20", 0.0, w), NAN)


def g_obs(z, w=PRIMARY_WIN):
    return lambda r: r["obs"].get((z, w), NAN)


def block_metrics(recs, betas, full=True):
    """一组事件的主指标：R_RIM(1 m)、R_S20(0.5 m)，full 时加 P5、β。"""
    out = {"R_RIM_1m": ratio_stats(recs, g_obs(1.0), g_rim(1.0), "ΣΔS_obs(1 m)/ΣΔS_RIM3(1 m)", tost=True),
           "R_S20_05m": ratio_stats(recs, g_obs(0.5), g_s20(), "ΣΔS_obs(0.5 m)/ΣΔS_S20（深度不对等，上界对照）")}
    if not full:
        return out
    out["P5"] = ratio_stats(recs, g_obs(5.0), g_obs(1.0), "P5_obs＝ΣΔS_obs(5 m)/ΣΔS_obs(1 m)",
                            extra_pairs=[("P5_RIM3", g_rim(5.0), g_rim(1.0))])
    d4 = [r for r in recs if isnum(r["obs"][(0.5, PRIMARY_WIN)]) and r["obs"][(0.5, PRIMARY_WIN)] <= p1.D4_DS_THRESH]
    for nm, sub in (("beta_D4subset", d4), ("beta_all_info", recs)):
        st = ratio_stats(sub, lambda r: r.get("dpco2", NAN), g_obs(0.5), "β_obs＝ΣΔpCO₂(T 归一)/ΣΔS(0.5 m) [µatm/psu]")
        use = [r for r in sub if isnum(r.get("dpco2", NAN)) and isnum(r["obs"][(0.5, PRIMARY_WIN)])
               and betas.get(r["station"]) is not None]
        den = sum(r["obs"][(0.5, PRIMARY_WIN)] for r in use)
        st["beta_model_pCO2_pooled"] = rnd(sum(betas[r["station"]][1] * r["obs"][(0.5, PRIMARY_WIN)] for r in use) / den) \
            if use and den else None
        st["beta_model_fCO2_pooled"] = rnd(sum(betas[r["station"]][0] * r["obs"][(0.5, PRIMARY_WIN)] for r in use) / den) \
            if use and den else None
        st["n_events_with_beta_model"] = len(use)
        out[nm] = st
    return out


def per_station_points(recs):
    out = {}
    for name in sorted({r["station"] for r in recs}):
        sub = [r for r in recs if r["station"] == name]
        d = {"n_events": len(sub), "group": sub[0]["group"]}
        for key, nf, df in (("R_RIM_1m", g_obs(1.0), g_rim(1.0)), ("R_S20_05m", g_obs(0.5), g_s20()),
                            ("P5_obs", g_obs(5.0), g_obs(1.0))):
            rows = [(nf(r), df(r)) for r in sub if isnum(nf(r)) and isnum(df(r))]
            den = sum(b for _, b in rows)
            d[key] = {"n": len(rows), "point": rnd(sum(a for a, _ in rows) / den) if rows and den else None}
        out[name] = d
    return out


def k27_block(recs):
    """K27 探索性变体：当前项深度因子取 t=1800 s。不参与出口、不进 D3；只在出口判定之后调用。"""
    out = {"status": "探索性，不参与出口，不进 D3（K27）",
           "rule": (f"当前项深度衰减因子 exp(−z²/(4·Kz·t)) 取 t={K27_T_CUR_DEPTH_S:.0f} s（主结果 t={T_CURRENT_S:.0f} s）；"
                    "量级项 c2·IRR·1800/√Kz、历史项、d0、风、窗口、事件集、bootstrap（B、种子、站×季簇）都与主结果相同"),
           "n_events_model_valid": {f"{z}m": sum(1 for r in recs if isnum(g_rim(z, K27_TAG)(r))) for z in K27_DEPTHS}}
    drop = ("ci95_outside_sesoi", "ci95_inside_sesoi", "ci95_excludes_1_info")
    for z in K27_DEPTHS:
        st = ratio_stats(recs, g_obs(z), g_rim(z, K27_TAG), f"探索性 K27：ΣΔS_obs({z} m)/ΣΔS_RIM3,K27({z} m)")
        out[f"R_RIM_{z}m"] = {k: v for k, v in st.items() if k not in drop}
    ref = [(g_obs(0.5)(r), g_rim(0.5)(r)) for r in recs if isnum(g_obs(0.5)(r)) and isnum(g_rim(0.5)(r))]
    den = sum(b for _, b in ref)
    out["reference_main_path_R_RIM_0.5m_point"] = {
        "n": len(ref), "point": rnd(sum(a for a, _ in ref) / den) if ref and den else None,
        "note": "主路径（t=1 s）在 0.5 m 的比值点估计，只作 K27 对照"}
    strata = {}
    for grp in (GROUP_A, GROUP_B):
        sub = [r for r in recs if r["group"] == grp]
        d = {"n_events": len(sub)}
        for z in K27_DEPTHS:
            rows = [(g_obs(z)(r), g_rim(z, K27_TAG)(r)) for r in sub if isnum(g_obs(z)(r)) and isnum(g_rim(z, K27_TAG)(r))]
            dd = sum(b for _, b in rows)
            d[f"R_RIM_{z}m"] = {"n": len(rows), "point": rnd(sum(a for a, _ in rows) / dd) if rows and dd else None}
        strata[grp] = d
    out["strata_by_group_points"] = strata
    return out


def k28_layer_points(recs):
    """K28 信息块：R_RIM(1 m) 按「替代层 vs 真 1 m」分层的点估计（不进出口）。"""
    out = {"status": "仅信息，不进出口（K28）",
           "rule": "替代层＝|s1_depth_m − 1.0| > 0.05 m（p1b_extend 同一判据）；模型「1 m」层深度＝s1_depth_m"}
    for nm, sub in (("true_1m", [r for r in recs if not r["s1_substitute"]]),
                    ("substitute_layer", [r for r in recs if r["s1_substitute"]])):
        rows = [(g_obs(1.0)(r), g_rim(1.0)(r)) for r in sub if isnum(g_obs(1.0)(r)) and isnum(g_rim(1.0)(r))]
        den = sum(b for _, b in rows)
        depths = {}
        for r in sub:
            k = f"{r['s1_depth_m']:.3f}"
            depths[k] = depths.get(k, 0) + 1
        out[nm] = {"n_events": len(sub), "n_ratio": len(rows),
                   "point": rnd(sum(a for a, _ in rows) / den) if rows and den else None,
                   "model_depth_m_counts": dict(sorted(depths.items()))}
    return out


def model_availability(recs):
    """K29：主集（CMORPH、0–6 h）模型可用＝ΔS_RIM(1 m 层) 与 ΔS_S20 均非 NaN；返回计数与缺测原因。"""
    n_ok, why = 0, {}
    for r in recs:
        if isnum(g_rim(1.0)(r)) and isnum(g_s20()(r)):
            n_ok += 1
            continue
        d = r["diag"].get("cmorph", {})
        k = (r["reason"].get("cmorph") or ("d0_out_of_range_gt_200mmh" if d.get("p_gt_d0_range", 0) > 0 else
                                            "s0_1m_pre_missing" if not isnum(r["s0"][1.0]) else "other_nan"))
        why[k] = why.get(k, 0) + 1
    n = len(recs)
    return {"n_total": n, "n_available": n_ok, "fraction": rnd(n_ok / n, 4) if n else None,
            "threshold": K29_MIN_AVAIL, "pass": bool(n and n_ok / n >= K29_MIN_AVAIL),
            "missing_reason_counts": dict(sorted(why.items())),
            "rule": "主集（10 mm 全部事件，CMORPH、0–6 h）中 ΔS_RIM(1 m 层) 与 ΔS_S20 均非 NaN 的比例 ≥50% 才继续（K29）"}


def sensitivity_block(R, recs):
    """K21：四项敏感性，只报 R_RIM(1 m)、R_S20(0.5 m)，不进出口。"""
    sens = {}
    for w in SENS_WINDOWS:
        sens[f"window_{w[0]}-{w[1]}h"] = {
            "R_RIM_1m": ratio_stats(R, g_obs(1.0, w), g_rim(1.0, w=w), f"窗口 {w}", tost=False),
            "R_S20_05m": ratio_stats(R, g_obs(0.5, w), g_s20(w=w), f"窗口 {w}")}
    sens["gauge_driven"] = {"R_RIM_1m": ratio_stats(R, g_obs(1.0), g_rim(1.0, "gauge"), "雨量计驱动"),
                            "R_S20_05m": ratio_stats(R, g_obs(0.5), g_s20("gauge"), "雨量计驱动"),
                            "gauge_hours_filled_0": sum(r["diag"].get("gauge", {}).get("gauge_filled", 0) for r in R),
                            "gauge_hours_negative_set_0": sum(r["diag"].get("gauge", {}).get("gauge_neg", 0) for r in R)}
    sens["exclude_BOBOA"] = block_metrics([r for r in R if r["station"] != "BOBOA"], {}, full=False)
    for thr in SENS_THRESHOLDS:
        sens[f"threshold_{int(thr)}mm"] = dict(block_metrics(recs[thr], {}, full=False), n_events=len(recs[thr]))
    return sens


# ======================================================================== 出图
def plot(recs, path):
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {GROUP_A: "#2a78d6", GROUP_B: "#eb6834"}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), facecolor="#fcfcfb")
    for ax, (xf, yf, xl, yl, ttl) in zip(axes, [
            (g_rim(1.0), g_obs(1.0), "ΔS RIM-3 (1 m, 0–6 h) [psu]", "ΔS observed (1 m) [psu]", "RIM-3 at 1 m"),
            (g_s20(), g_obs(0.5), "ΔS S20 (surface) [psu]", "ΔS observed (0.5 m) [psu]", "S20 vs 0.5 m")]):
        ax.set_facecolor("#fcfcfb")
        lim = 0.05
        for grp in (GROUP_A, GROUP_B):
            pts = [(xf(r), yf(r)) for r in recs if r["group"] == grp and isnum(xf(r)) and isnum(yf(r))]
            if pts:
                xs, ys = zip(*pts)
                lim = max(lim, max(abs(v) for v in xs + ys))
                ax.scatter(xs, ys, s=16, c=col[grp], alpha=0.7, linewidths=0.4, edgecolors="#fcfcfb",
                           label=f"{'TAO/RAMA (8)' if grp == GROUP_A else 'O group (5)'}  n={len(pts)}")
        lim *= 1.05
        for k, ls in ((1.0, "-"), (SESOI_LO, ":"), (SESOI_HI, ":")):
            ax.plot([-lim, lim], [-k * lim, k * lim], ls, color="#52514e", lw=1)
        ax.axhline(0, color="#d6d5d0", lw=0.8)
        ax.axvline(0, color="#d6d5d0", lw=0.8)
        ax.set_xlim(-lim, lim * 0.2)
        ax.set_ylim(-lim, lim * 0.2)
        ax.set_xlabel(xl, color="#0b0b0b")
        ax.set_ylabel(yl, color="#0b0b0b")
        ax.set_title(ttl + "  (solid 1:1, dotted 0.65/1.35)", color="#0b0b0b", fontsize=10)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=fig.get_facecolor())
    plt.close(fig)


# ======================================================================== 主流程
EVENT_FIELDS = ["station", "group", "onset_utc", "season", "rain_24h_mm", "s1_substitute", "n_ctrl",
                "dS05_obs", "dS1_obs", "dS5_obs", "dS_rim_05", "dS_rim_1", "dS_rim_5", "dS_s20",
                "dS_rim_1_gauge", "dS_s20_gauge", "dpCO2_Tnorm", "model_reason_cmorph", "wind_floor_hours",
                "dS_rim_05_k27_explor", "dS_rim_1_k27_explor", "model_depth_s1_layer_m"]


def event_row(r):
    m = r["model"]
    W = PRIMARY_WIN
    return {"station": r["station"], "group": r["group"], "onset_utc": r["onset_utc"], "season": r["season"],
            "rain_24h_mm": rnd(r["acc24"], 3), "s1_substitute": int(r["s1_substitute"]), "n_ctrl": r["n_ctrl"],
            "dS05_obs": rnd(r["obs"][(0.5, W)], 4), "dS1_obs": rnd(r["obs"][(1.0, W)], 4),
            "dS5_obs": rnd(r["obs"][(5.0, W)], 4),
            "dS_rim_05": rnd(m.get(("cmorph", "rim", 0.5, W)), 4), "dS_rim_1": rnd(m.get(("cmorph", "rim", 1.0, W)), 4),
            "dS_rim_5": rnd(m.get(("cmorph", "rim", 5.0, W)), 4), "dS_s20": rnd(m.get(("cmorph", "s20", 0.0, W)), 4),
            "dS_rim_1_gauge": rnd(m.get(("gauge", "rim", 1.0, W)), 4),
            "dS_s20_gauge": rnd(m.get(("gauge", "s20", 0.0, W)), 4),
            "dpCO2_Tnorm": rnd(r.get("dpco2", NAN), 3), "model_reason_cmorph": r["reason"].get("cmorph") or "",
            "wind_floor_hours": r["diag"].get("cmorph", {}).get("wind_floor", ""),
            "dS_rim_05_k27_explor": rnd(m.get((K27_TAG, "rim", 0.5, W)), 4),
            "dS_rim_1_k27_explor": rnd(m.get((K27_TAG, "rim", 1.0, W)), 4),
            "model_depth_s1_layer_m": rnd(r["s1_depth_m"], 3)}


def read_p1b_gates(p1b_dir):
    with open(os.path.join(p1b_dir, "p1b_summary.json"), encoding="utf-8") as f:
        s = json.load(f)
    gv = s["gates_verdict"]
    if s["n_events"]["merged"] != EXPECTED_MERGED or not all(gv[k] for k in ("D1", "D2", "D4", "D5")):
        raise RuntimeError(f"P1b 合并集判定不符预期（K16）：n={s['n_events']}, verdict={gv}")
    if sorted(s.get("o_group_stations", [])) != sorted(EXPECTED_O_STATIONS):
        raise RuntimeError(f"P1b O 组站单不符：{s.get('o_group_stations')}")
    g = s["gates_merged"]
    return gv, {"D1": {k: g["D1"][k] for k in ("n_events_total", "n_stations_ge10", "pass")},
                "D2": {k: g["D2"].get(k) for k in ("diff_psu", "ci95_cluster_bootstrap", "pass")},
                "D4": {k: g["D4"][k] for k in ("n_events_dS05_le_-0.2", "pass")},
                "D5": {k: g["D5"].get(k) for k in ("median_pre_diff_psu", "corr_event_anom", "pass")}}


def build_all(args, out_dir, log):
    cache = os.path.join(out_dir, "cache")
    os.makedirs(cache, exist_ok=True)
    seeds = [os.path.join(args.p1b_dir, "cache"), os.path.join(os.path.dirname(args.p1_events), "cache")]
    http = LinkSeedHttp(log, cache, seeds)
    stations = rebuild_stations(http, cache, log)
    ev = {thr: [] for thr in [PRIMARY_THRESHOLD] + SENS_THRESHOLDS}
    for stn in stations:
        for thr in ev:
            evs = events_at(stn, thr, log if thr == PRIMARY_THRESHOLD else None)
            for e in evs:
                e["_stn"] = stn["name"]
            ev[thr] += evs
    val = validate_rebuild(ev[PRIMARY_THRESHOLD], args.p1_events, os.path.join(args.p1b_dir, "p1b_events_new.csv"))
    print(f"[重建] 10 mm 事件 {len(ev[PRIMARY_THRESHOLD])}（与 P1/P1b 产物逐条一致），D4={val['d4_count']}；"
          f"5 mm {len(ev[5.0])}，20 mm {len(ev[20.0])}", flush=True)
    k28 = {}
    for thr, evs in ev.items():                    # K28：下载前先校验全部事件的 s1_depth_m（缺即报错，早失败）
        for e in evs:
            z1 = model_depths(e, e["_stn"])[1.0]
            key = f"{thr:g}mm"
            k28.setdefault(key, {})
            k28[key][f"{e['_stn']}@{z1:.3f}"] = k28[key].get(f"{e['_stn']}@{z1:.3f}", 0) + 1
    hrs = {thr: needed_hours(v) for thr, v in ev.items()}
    prim = sorted(hrs[PRIMARY_THRESHOLD])
    extra = sorted((hrs[5.0] | hrs[20.0]) - hrs[PRIMARY_THRESHOLD])
    plan = {"k28_model_depth_counts": {k: dict(sorted(v.items())) for k, v in k28.items()},"window_h": [WIN_LO_H, WIN_HI_H], "events": {str(k): len(v) for k, v in ev.items()},
            "unique_hours": {"10mm_primary": len(prim), "extra_for_5mm_20mm": len(extra),
                             "total": len(prim) + len(extra)},
            "est_download_GB": {"primary": round(len(prim) * CMORPH_MEAN_MB / 1024, 1),
                                "total": round((len(prim) + len(extra)) * CMORPH_MEAN_MB / 1024, 1)},
            "est_hours_at_MBps": {str(r): round((len(prim) + len(extra)) * CMORPH_MEAN_MB / r / 3600, 2)
                                  for r in (2, 4, 8)},
            "cache_links": http.n_linked, "cache_http_requests": http.n_requests,
            "parts": plan_parts(len(prim), len(extra))}
    return stations, ev, val, prim, extra, plan


def plan_parts(n_prim, n_extra):
    """--part 两段的唯一小时数与字节估计：小时数×CMORPH_MEAN_MB（1e6 B）；时长按 workstation 实测 PCA_FILES_PER_S。"""
    def est(nh):
        return {"unique_hours": nh, "est_bytes": int(round(nh * CMORPH_MEAN_MB * 1e6)),
                "est_GB": round(nh * CMORPH_MEAN_MB / 1000, 1),
                "est_hours_at_pca_rate": round(nh / PCA_FILES_PER_S / 3600, 2)}
    return {"basis": (f"每小时文件 {CMORPH_MEAN_MB} MB（HEAD 均值，GB＝1e9 B）；时长按 workstation 4 并发约 {PCA_FILES_PER_S} 文件/s"
                      "（09-26 实测）；sensitivity 段需主集＋5/20 mm 全部小时，续主集像元缓存时只下主集之外的部分"),
            "primary": est(n_prim),
            "sensitivity_beyond_primary": dict(est(n_extra), note="--resume-from 主集输出目录时，本段只需下这些小时"),
            "sensitivity_if_no_resume": dict(est(n_prim + n_extra), note="不继承主集像元缓存时本段需下的全部小时"),
            "all": est(n_prim + n_extra)}


def run_plan(args, out_dir, log):
    stations, ev, val, prim, extra, plan = build_all(args, out_dir, log)
    heads = {}
    for url in [GLODAP_URL, "https://" + CMORPH_HOST + hour_url(prim[len(prim) // 2])]:
        try:
            req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                heads[url] = {"status": r.status, "content_length": r.headers.get("Content-Length")}
        except Exception as e:
            heads[url] = {"error": f"{type(e).__name__}: {e}"}
    plan.update({"mode": "plan（未下载 CMORPH/GLODAP 数据文件）", "rebuild_validation": val, "head_checks": heads,
                 "stations": [{k: s[k] for k in ("name", "group", "lat", "lon")} for s in stations]})
    jdump(plan, os.path.join(out_dir, "p2_plan.json"))
    print(json.dumps({k: plan[k] for k in ("events", "unique_hours", "est_download_GB", "est_hours_at_MBps", "parts")},
                     ensure_ascii=False))
    return 0


def primary_d3_ref(resume_from, out_dir):
    """sensitivity 段：主集 D3 须已写盘（K11 顺序）。在 resume 目录找 p2_d3.json 或本段副本，复制为 p2_d3_from_primary.json。"""
    if not resume_from:
        raise RuntimeError("--part sensitivity 需要 --resume-from 指向主集输出目录（含 p2_d3.json；K11：D3 先于观测量写盘）")
    for nm in ("p2_d3.json", "p2_d3_from_primary.json"):
        src = os.path.join(resume_from, nm)
        if os.path.exists(src):
            with open(src, encoding="utf-8") as f:
                d3 = json.load(f)
            dst = os.path.join(out_dir, "p2_d3_from_primary.json")
            if os.path.abspath(src) != os.path.abspath(dst):
                shutil.copyfile(src, dst)
            return {"source": src, "sha256": sha256_file(src), "pass": d3.get("pass"), "verdict": d3.get("verdict"),
                    "written_utc": d3.get("written_utc")}
    raise RuntimeError(f"--part sensitivity：{resume_from} 下没有 p2_d3.json／p2_d3_from_primary.json（主集 D3 未写盘，K11）")


def run_full(args, out_dir, log, t0):
    part = getattr(args, "part", "all") or "all"
    do_prim, do_sens = part in ("all", "primary"), part in ("all", "sensitivity")
    d3_ref = primary_d3_ref(args.resume_from, out_dir) if part == "sensitivity" else None   # 先查，早失败
    gv, gdetail = read_p1b_gates(args.p1b_dir)
    stations, ev, val, prim, extra, plan = build_all(args, out_dir, log)
    jdump(dict(plan, mode="full", part=part),
          os.path.join(out_dir, "p2_plan_sensitivity.json" if part == "sensitivity" else "p2_plan.json"))
    sbn = {s["name"]: s for s in stations}
    if args.smoke:
        keep = ev[PRIMARY_THRESHOLD][:3]
        ev = {PRIMARY_THRESHOLD: keep, 5.0: [], 20.0: []}
        prim, extra = sorted(needed_hours(keep)), []
    if not do_sens:                                # primary：只主集事件、只主集小时（K22 主集小时）
        ev = {PRIMARY_THRESHOLD: ev[PRIMARY_THRESHOLD]}
        extra = []
    # GLODAP（小，先做；敏感性段不需要 β）
    glo, betas = None, {}
    if do_prim:
        glo = glodap_values(stations, out_dir, log, args.resume_from)
        for s in stations:
            g = glo["stations"].get(s["name"], {})
            if all(isnum(g.get(k)) for k in GLODAP_VARS):
                betas[s["name"]] = beta_model(g["TAlk"], g["TCO2"], g["temperature"], g["salinity"])
            else:
                betas[s["name"]] = None
    # CMORPH
    pix = PixelCache(os.path.join(out_dir, "cmorph_pixels.csv"))
    pix.load(os.path.join(args.resume_from, "cmorph_pixels.csv") if args.resume_from else None, pix.path)
    dl = run_downloads(prim + extra, pix, stations, out_dir, log, args.max_conc)
    needed = set(prim) | set(extra)
    dl["missing_hours"] = sum(1 for h in needed if h in pix.missing)
    dl["missing_reasons"] = {}
    for h in needed:
        if h in pix.missing:
            dl["missing_reasons"][pix.missing[h].split(":")[0]] = dl["missing_reasons"].get(pix.missing[h].split(":")[0], 0) + 1
    # 模型（先于任何观测比值）
    wins = [PRIMARY_WIN] + SENS_WINDOWS if do_sens else [PRIMARY_WIN]
    recs = {}
    for thr, evs in ev.items():
        variants = ("cmorph", "gauge") if thr == PRIMARY_THRESHOLD and do_sens else ("cmorph",)
        windows = wins if thr == PRIMARY_THRESHOLD else (PRIMARY_WIN,)
        recs[thr] = [event_record(sbn[e["_stn"]], e, pix, variants, windows, k27=thr == PRIMARY_THRESHOLD and do_prim)
                     for e in evs]
    R = recs[PRIMARY_THRESHOLD]
    avail = {}
    for r in R:
        k = r["reason"]["cmorph"] or "ok"
        avail[k] = avail.get(k, 0) + 1
    k29 = model_availability(R)
    if do_prim and not k29["pass"]:                # K29：D3 之前中止，不写 d3/summary，不判出口
        jdump({"status": "aborted_K29", "script": "p2_rim_test.py", "version": VERSION, "part": part,
               "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               "n_total": k29["n_total"], "n_available": k29["n_available"], "fraction": k29["fraction"],
               "threshold": K29_MIN_AVAIL, "missing_reason_counts": k29["missing_reason_counts"],
               "forcing_availability_primary": avail, "cmorph": dl,
               "note": "主集模型可用率 <50%（K29）：未写 p2_d3.json 与 p2_summary.json，未做出口判定"},
              os.path.join(out_dir, "p2_abort.json"))
        log.log(f"K29 中止：可用 {k29['n_available']}/{k29['n_total']}，原因 {k29['missing_reason_counts']}", echo=True)
        return 3
    d3 = None
    if do_prim:
        d3 = d3_gate(R, sbn)
        d3.update({"written_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                   "note": "D3 只用模型预测与对照窗观测；本文件写盘先于 P5/β/任何比值计算（K11）"})
        jdump(d3, os.path.join(out_dir, "p2_d3.json"))
        print(f"D3：{d3['verdict']}", flush=True)
    # 观测侧（事件 ΔS 与 ΔpCO₂ 在 D3 写盘之后才计算；sensitivity 段 D3 已由主集写盘）
    for thr, rs in recs.items():
        for r in rs:
            fill_obs(r, sbn[r["station"]])
    sens = sensitivity_block(R, recs) if do_sens else None
    if part == "sensitivity":
        out = {"script": "p2_rim_test.py", "version": VERSION, "part": part, "smoke": bool(args.smoke),
               "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               "runtime_s": round(time.monotonic() - t0, 1),
               "status": "四项敏感性（K21）；不做出口判定（K26），不写 p2_summary.json",
               "inputs": {"p1_events_csv": args.p1_events, "p1_events_sha256": sha256_file(args.p1_events),
                          "p1b_dir": args.p1b_dir, "resume_from": args.resume_from},
               "primary_d3_reference": d3_ref, "rebuild_validation": val,
               "cmorph": dict(dl, unique_hours_primary=len(prim), unique_hours_extra=len(extra)),
               "model_availability_primary_info": k29, "forcing_availability_primary": avail,
               "wind": {"factor_4m_to_10m": rnd(WIND_FACTOR, 4), "floor_mps": WIND_FLOOR,
                        "hours_floored": sum(r["diag"].get("cmorph", {}).get("wind_floor", 0) for r in R)},
               "sensitivity": sens}
        jdump(out, os.path.join(out_dir, "p2_sensitivity.json"))
        print("敏感性已写 p2_sensitivity.json（不判出口）", flush=True)
        return 0
    prim_m = block_metrics(R, betas)
    ex = decide_exit(gv, d3, prim_m["R_RIM_1m"], prim_m["R_S20_05m"])
    strata = {}
    for grp in (GROUP_A, GROUP_B):
        sub = [r for r in R if r["group"] == grp]
        strata[grp] = block_metrics(sub, betas)
        strata[grp]["n_events"] = len(sub)
        strata[grp]["n_s1_substitute"] = sum(1 for r in sub if r["s1_substitute"])
    if not do_sens:
        sens = {"status": "未跑（--part primary）；四项敏感性见 --part sensitivity 的 p2_sensitivity.json（不进出口，K26）"}
    k27 = k27_block(R)                             # K27：出口已判定之后才算，只报告
    k28 = k28_layer_points(R)                      # K28：分层点估计，只报告
    ctrl1 = [c for r in R for c in r["ctrl_ds1"] if isnum(c)]
    summary = {
        "script": "p2_rim_test.py", "version": VERSION, "part": part, "smoke": bool(args.smoke),
        "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "runtime_s": round(time.monotonic() - t0, 1),
        "implementation_choices": "K1–K30（脚本头）",
        "inputs": {"p1_events_csv": args.p1_events, "p1_events_sha256": sha256_file(args.p1_events),
                   "p1b_dir": args.p1b_dir,
                   "p1b_summary_sha256": sha256_file(os.path.join(args.p1b_dir, "p1b_summary.json")),
                   "resume_from": args.resume_from},
        "gates_carried_from_p1b": {"verdict": gv, "detail": gdetail, "rule": "K16：不重判"},
        "rebuild_validation": val,
        "cmorph": dict(dl, unique_hours_primary=len(prim), unique_hours_extra=len(extra)),
        "model_availability_primary": avail,
        "model_availability_k29": k29,
        "model_nan_with_valid_forcing_primary": {
            "events_cmorph_halfstep_gt_200mmh": sum(1 for r in R if r["diag"].get("cmorph", {}).get("p_gt_d0_range", 0) > 0),
            "rim_1m_nan_events": sum(1 for r in R if r["reason"].get("cmorph") is None
                                     and not isnum(g_rim(1.0)(r))),
            "note": "强迫齐全但模型 NaN：CMORPH 半步 >200 mm/h 落在 d0 表外（K2，照 Witte 不外推）或该深度雨前观测缺"},
        "wind": {"factor_4m_to_10m": rnd(WIND_FACTOR, 4), "floor_mps": WIND_FLOOR,
                 "hours_floored": sum(r["diag"].get("cmorph", {}).get("wind_floor", 0) for r in R)},
        "glodap": glo, "beta_model_by_station": {k: (None if v is None else {"fCO2_slope": rnd(v[0], 4),
                                                                            "pCO2_slope": rnd(v[1], 4)})
                                                 for k, v in betas.items()},
        "D3": d3,
        "primary": prim_m,
        "exit": ex,
        "strata_by_group": strata,
        "per_station_points": per_station_points(R),
        "k28_R_RIM_1m_by_layer": k28,
        "sensitivity": sens,
        "exploratory_k27": k27,
        "descriptive_info": {"control_dS1_0_6h_mean": rnd(nanmean(ctrl1)), "n_control_windows": len(ctrl1),
                             "note": "K25：比值不减对照；此处只描述"},
        "disclosure": ["扩样在 P1 D4 不过之后决定（P1 原始 D4=14）", "D4=34 只比阈值多 4，新增强淡化事件全在 O 组（1 m 为替代层）",
                       "O 组为中纬度/副热带站，与 P1 暖池站气候背景不同 → 见 strata_by_group",
                       "S20 是海面式，在 0.5 m 只是上界对照（深度不对等）", "0.5 m SSS 在部分站可能即 1 m 传感器（P1 I17）",
                       "β_obs 主集按 ΔS0.5≤−0.2 筛选（分母筛选）会使 β_obs 偏向 0（K18）",
                       "K28：O 组「1 m」为 ≤1.5 m 最浅替代层，模型按实际深度算 → 见 k28_R_RIM_1m_by_layer",
                       "K30：观测 0.5 m 约 3 小时一次而模型用满 6 小时（增噪）；模型小时相对观测超前约 15–45 分钟"],
    }
    if args.smoke:
        summary["exit"] = {"exit": None, "label": "SMOKE：只验证流程，结论无效"}
    jdump(summary, os.path.join(out_dir, "p2_summary.json"))
    with open(os.path.join(out_dir, "p2_events.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=EVENT_FIELDS)
        w.writeheader()
        w.writerows(event_row(r) for r in R)
    try:
        plot(R, os.path.join(out_dir, "p2_fig_ratio.png"))
    except Exception as e:
        log.log(f"出图失败（不影响判定）：{e}")
    rr, rs = prim_m["R_RIM_1m"], prim_m["R_S20_05m"]
    print(f"R_RIM(1 m)={rr.get('point')} CI95={rr.get('ci95')} TOST={rr.get('tost', {}).get('equivalent')}；"
          f"R_S20(0.5 m)={rs.get('point')} CI95={rs.get('ci95')}")
    print(f"出口：{summary['exit']['label']}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="P2 系泊检验")
    ap.add_argument("--out")
    ap.add_argument("--p1-events", default=P1_EVENTS_DEFAULT)
    ap.add_argument("--p1b-dir", default=P1B_DIR_DEFAULT)
    ap.add_argument("--resume-from")
    ap.add_argument("--max-conc", type=int, default=MAX_CONC)
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--part", choices=PARTS, default="all")
    args = ap.parse_args(argv)
    if args.selftest:
        return selftest()
    if args.resume_from and not os.path.isdir(args.resume_from):
        print(f"--resume-from 目录不存在：{args.resume_from}", file=sys.stderr)
        return 3
    args.max_conc = max(1, min(MAX_CONC, args.max_conc))
    out_dir = args.out or os.environ.get("REPRO_OUTPUT_DIR")
    if not out_dir:
        print("需要 --out 或环境变量 REPRO_OUTPUT_DIR", file=sys.stderr)
        return 3
    os.makedirs(out_dir, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(out_dir, "tmp_mpl"))
    log = p1.Log(os.path.join(out_dir, "p2_log.txt"))
    t0 = time.monotonic()
    log.log(f"=== start {VERSION} part={args.part} plan={args.plan} smoke={args.smoke} resume={args.resume_from} "
            f"out={out_dir}", echo=True)
    try:
        rc = run_plan(args, out_dir, log) if args.plan else run_full(args, out_dir, log, t0)
    except (DataSourceError, p1.FetchError) as e:
        log.log(f"FATAL 数据源故障（可 --resume-from 续跑）：{e}", echo=True)
        return 2
    except Exception:
        log.log("FATAL 未预期异常：\n" + traceback.format_exc(), echo=True)
        return 3
    shutil.rmtree(os.path.join(out_dir, "tmp_mpl"), ignore_errors=True)
    log.log("=== done")
    log.close()
    return rc


# ======================================================================== 自测（合成数据，无网络）
def _rim_literal(P, U, S):
    """Witte RIMv3 的逐行移植（numpy 版、手写双线性插值），只用于自测对照。"""
    np = _np()

    def bil(u, r):
        ug, rg, T = np.array(U_GRID, float), np.array(R_GRID, float), np.array(DL, float)
        out = []
        for a, b in zip(np.atleast_1d(u), np.atleast_1d(r)):
            if not (0 <= a <= 200 and 0 <= b <= 200):
                out.append(np.nan)
                continue
            i = min(np.searchsorted(ug, a, side="right") - 1, len(ug) - 2)
            j = min(np.searchsorted(rg, b, side="right") - 1, len(rg) - 2)
            tu = (a - ug[i]) / (ug[i + 1] - ug[i])
            tr = (b - rg[j]) / (rg[j + 1] - rg[j])
            out.append((1 - tu) * (1 - tr) * T[i, j] + tu * (1 - tr) * T[i + 1, j] + (1 - tu) * tr * T[i, j + 1]
                       + tu * tr * T[i + 1, j + 1])
        return np.array(out)

    tc = np.arange(600, 1800, 25)
    ti = -np.arange(0, 24 * 3600, 1800) + 24 * 3600
    res = []
    for t in np.arange(0, len(P) - 48) + 48:
        IRRi = P[t - 48:t] / 1000 / 3600
        Kzi = 2.5 * (10 ** (-5)) * (U[t - 48:t] ** 2)
        d0i = bil(U[t - 48:t], P[t - 48:t])
        prior = d0i / (d0i + 5.7 * IRRi * tc / (np.sqrt(Kzi * ti)))
        IRR = P[t] / 1000 / 3600
        Kz = 2.5 * (10 ** (-5)) * U[t] ** 2
        d0c = bil(U[t], P[t])[0]
        cur = d0c / (d0c + (1 / 11) * IRR * 1800 / (np.sqrt(Kz)))
        res.append(S[t] * (np.prod(prior) * cur))
    return np.array(res)


def _syn_station(name, group, z_obs, n, h0, onsets_main, onsets_extra, seed, true_ratio=0.5):
    """合成站：主事件 4 h×8 mm/h，副事件 2 h×3 mm/h；观测「1 m」层按 z_obs 深度生成（K28 自测）。"""
    np = _np()
    rng = np.random.default_rng(seed)
    ser = {k: np.full(n, np.nan) for k in ("rain", "wind", "s1", "s5", "sss05", "pco2", "sst")}
    ser["wind"][:] = 5.0
    ser["rain"][:] = 0.0
    Pfull = np.zeros(2 * n)
    for i0, nh, rate in [(i, 4, 8.0) for i in onsets_main] + [(i, 2, 3.0) for i in onsets_extra]:
        for k in range(nh):
            Pfull[2 * (i0 + k):2 * (i0 + k) + 2] = rate
            ser["rain"][i0 + k] = rate
    base = 35.0 + rng.normal(0, 0.003, n)
    Uh = np.repeat(ser["wind"] * WIND_FACTOR, 2)
    for zt, key in ((0.5, "sss05"), (z_obs, "s1"), (5.0, "s5")):
        F = np.ones(2 * n)
        F[48:] = rim_factor(Pfull, Uh, zt)
        ser[key] = base * (1 - true_ratio * (1 - F.reshape(-1, 2).mean(axis=1)))
    ser["sss05"] = np.where(np.arange(n) % 3 == 0, ser["sss05"], np.nan)
    ser["sst"][:] = 28.0
    ser["pco2"] = np.where(np.arange(n) % 3 == 0, 400.0 + 10.0 * (ser["sss05"] - 35.0), np.nan)
    return {"name": name, "group": group, "h0": h0, "n": n, "ser": ser, "lat": 0.0, "lon": 180.0,
            "s1_depth": None}, Pfull


def _e2e_fixture(td, z_b=1.5, seed=11):
    """端到端合成夹具：两站（A 真 1 m、B 替代层 z_b），10/5/20 mm 事件集，假下载/GLODAP/P1b 门。"""
    n, h0 = 24 * 800, 24 * 16000
    main_on = list(range(24 * 40 + 5, n - 24 * 40, 24 * 11))
    extra_on = [i + 24 * 5 for i in main_on]
    seasons = ["DJF", "MAM", "JJA", "SON"]
    stations, cm = [], {}
    ev = {PRIMARY_THRESHOLD: [], 5.0: [], 20.0: []}
    for j, (nm, grp, z1) in enumerate((("SYNA", GROUP_A, 1.0), ("SYNB", GROUP_B, z_b))):
        stn, Pfull = _syn_station(nm, grp, z1, n, h0, main_on, extra_on, seed + j)
        stations.append(stn)
        for i0 in range(n):
            cm.setdefault(h0 + i0, {})[nm] = (Pfull[2 * i0], Pfull[2 * i0 + 1], 0.0, 0.0, 9, 9)

        def mk(i0, k, acc):
            ctrl = [i0 + 72, i0 - 72]
            return {"_stn": nm, "i": i0, "hour": h0 + i0, "onset_utc": p1.hour_to_iso(h0 + i0),
                    "season": seasons[k % 4], "acc24": acc, "n_ctrl": 2, "ctrl_hours": [h0 + c for c in ctrl],
                    "ctrl_ds1": [obs_delta(stn["ser"]["s1"], c, PRIMARY_WIN) for c in ctrl],
                    "s1_depth_m": z1, "s1_substitute": abs(z1 - 1.0) > 0.05}
        for k, i0 in enumerate(main_on):
            ev[PRIMARY_THRESHOLD].append(mk(i0, k, 32.0))
            ev[5.0].append(mk(i0, k, 32.0))
            if k % 2 == 0:
                ev[20.0].append(mk(i0, k, 32.0))
        for k, i0 in enumerate(extra_on):
            ev[5.0].append(mk(i0, k, 6.0))
    hrs = {thr: needed_hours(v) for thr, v in ev.items()}
    prim = sorted(hrs[PRIMARY_THRESHOLD])
    extra = sorted((hrs[5.0] | hrs[20.0]) - hrs[PRIMARY_THRESHOLD])
    p1csv = os.path.join(td, "p1_events.csv")
    p1bdir = os.path.join(td, "p1b")
    os.makedirs(p1bdir, exist_ok=True)
    for pth in (p1csv, os.path.join(p1bdir, "p1b_summary.json")):
        with open(pth, "w") as f:
            f.write("{}\n")
    fx = {"stations": stations, "ev": ev, "prim": prim, "extra": extra, "cm": cm, "p1csv": p1csv, "p1bdir": p1bdir,
          "downloaded": [], "missing_hours": set(), "val": {"n_rebuilt": len(ev[PRIMARY_THRESHOLD]), "synthetic": True}}

    def fake_build_all(args, out_dir, log):
        import copy
        plan = {"events": {str(k): len(v) for k, v in fx["ev"].items()},
                "unique_hours": {"10mm_primary": len(prim), "extra_for_5mm_20mm": len(extra)},
                "parts": plan_parts(len(prim), len(extra))}
        return fx["stations"], copy.deepcopy(fx["ev"]), dict(fx["val"]), list(prim), list(extra), plan

    def fake_gates(p1b_dir):
        return {"D1": True, "D2": True, "D4": True, "D5": True}, {"synthetic": True}

    def fake_glodap(stations_, out_dir, log, resume_dir=None):
        return {"stations": {s["name"]: {"TAlk": NAN, "TCO2": NAN, "salinity": NAN, "temperature": NAN}
                             for s in stations_}}

    def fake_downloads(hours, cache, stations_, out_dir, log, max_conc):
        todo = [h for h in hours if h not in cache.data and h not in cache.missing]
        fx["downloaded"].append(list(todo))
        for h in todo:
            if h in fx["missing_hours"]:
                cache.add_missing(h, "HTTP404")
            else:
                cache.add(h, fx["cm"][h])
        return {"hours_needed": len(hours), "downloaded_this_run": len(todo), "bytes_this_run": 0}

    fx["patches"] = {"build_all": fake_build_all, "read_p1b_gates": fake_gates, "glodap_values": fake_glodap,
                     "run_downloads": fake_downloads}
    return fx


def _e2e_run(g, fx, out, part, resume=None):
    """在模块全局 g 上打补丁跑 main(argv)，返回退出码。"""
    saved = {k: g[k] for k in fx["patches"]}
    g.update(fx["patches"])
    argv = ["--out", out, "--p1-events", fx["p1csv"], "--p1b-dir", fx["p1bdir"]]
    if part is not None:
        argv += ["--part", part]
    if resume:
        argv += ["--resume-from", resume]
    try:
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):
            return g["main"](argv)
    finally:
        g.update(saved)


def selftest():
    import tempfile
    np = _np()
    ok, skipped = [], []

    def check(cond, msg):
        if not cond:
            raise AssertionError(msg)
        ok.append(msg)

    rng = np.random.default_rng(7)
    # 1) d0 插值与 Witte 表
    check(len(DL) == 7 and all(len(r) == 7 for r in DL) and len(TC) == 48 and TC[-1] == 1775 and TI[0] == 86400
          and TI[-1] == 1800, "d0 表 7×7、tc 48 个（600…1775）、ti 86400…1800")
    check(abs(d0_interp(np.array([6.0]), np.array([10.0]))[0] - 2.9) < 1e-12
          and abs(d0_interp(np.array([5.0]), np.array([7.5]))[0] - (2.0 + 2.5 + 1.9 + 2.9) / 4) < 1e-12
          and np.isnan(d0_interp(np.array([250.0]), np.array([1.0]))[0]), "d0 双线性插值（格点、格心、界外 NaN）")
    # 2) RIM z=0 与逐行移植一致
    L = 48 + 60
    P = np.where(rng.random(L) < 0.3, rng.gamma(1.2, 6.0, L), 0.0)
    U = rng.uniform(1.0, 12.0, L)
    S = np.full(L, 35.0)
    lit = _rim_literal(P, U, S)
    mine = 35.0 * rim_factor(P, U, 0.0)
    check(np.allclose(lit, mine, rtol=1e-12, atol=1e-12), f"RIM z=0 与 Witte RIMv3 逐行移植一致（max|Δ|={np.max(np.abs(lit - mine)):.2e}）")
    F05, F1, F5, F50 = (rim_factor(P, U, z) for z in (0.5, 1.0, 5.0, 300.0))
    check(np.all(F05 <= F1 + 1e-15) and np.all(F1 <= F5 + 1e-15) and np.allclose(F50, 1.0),
          "深度因子：F 随 z 单调增、z→∞ 时 F→1")
    check(np.all((rim_factor(np.zeros(L), U, 1.0) == 1.0)), "无雨时 F≡1")
    check(abs(s20_delta(np.array([10.0]), np.array([5.0]))[0] - (-0.35 * 10 * 5 ** -0.77)) < 1e-12, "S20 公式")
    # 2b) K27：只改当前项深度因子；默认参数＝主路径；手算单步
    check(all(np.array_equal(rim_factor(P, U, z), rim_factor(P, U, z, t_cur_depth=T_CURRENT_S), equal_nan=True)
              for z in DEPTHS) and np.array_equal(rim_factor(P, U, 0.0), rim_factor(P, U, 0.0, K27_T_CUR_DEPTH_S)),
          "K27：默认参数与主路径逐值相同；z=0 时 K27 与主路径相同")
    Pk, Uk = np.zeros(49), np.full(49, 5.0)
    Pk[-1] = 10.0
    kz = KZ_COEF * 25.0
    d0h = (2.0 + 2.9) / 2
    cmag = C2 * (10.0 / 3.6e6) * 1800.0 / math.sqrt(kz * 1.0)
    f_main = d0h / (d0h + cmag * math.exp(-1.0 / (4 * kz * 1.0)))
    f_k27 = d0h / (d0h + cmag * math.exp(-1.0 / (4 * kz * 1800.0)))
    check(abs(rim_factor(Pk, Uk, 1.0)[0] - f_main) < 1e-12 and abs(rim_factor(Pk, Uk, 1.0, K27_T_CUR_DEPTH_S)[0] - f_k27) < 1e-12
          and f_k27 < f_main - 1e-3 and f_main > 1 - 1e-12, f"K27 手算单步（z=1 m，U=5，当前 10 mm/h）：主 {f_main:.6f}，K27 {f_k27:.6f}")
    check(all(np.all(rim_factor(P, U, z, K27_T_CUR_DEPTH_S) <= rim_factor(P, U, z) + 1e-15) for z in (0.5, 1.0, 5.0)),
          "K27 的 F 不大于主路径（当前项深度衰减更弱 → 淡化更多）")
    # 3) 风插补
    w = fill_wind([5.0, np.nan, np.nan, 8.0, np.nan] + [np.nan] * 7 + [1.0], 6)
    check(abs(w[1] - 6.0) < 1e-12 and abs(w[2] - 7.0) < 1e-12 and np.isnan(w[5]), "风缺口插补（≤6 h 插、>6 h 不插）")
    check(abs(WIND_FACTOR - 1.0865) < 1e-3, "4 m→10 m 风速换算系数")
    # 4) 像元缓存：写、断点、重载
    td = tempfile.mkdtemp(prefix="p2_selftest_")
    pc = PixelCache(os.path.join(td, "px.csv"))
    pc.load()
    pc.add(100, {"A": (1.0, 2.0, 1.5, 2.5, 9, 9), "B/C": (NAN, 0.0, 0.0, 0.0, 8, 9)})
    pc.add_missing(101, "HTTP404")
    with open(pc.path, "a") as f:
        f.write("102,A,1.0,1.0,1.0,1.0,9,9\n")               # 无 DONE 的半截块
    pc2 = PixelCache(os.path.join(td, "px2.csv"))
    pc2.load(pc.path)
    check(set(pc2.data) == {100} and pc2.missing == {101: "HTTP404"} and np.isnan(pc2.get(100, "B/C")[0])
          and pc2.get(100, "A")[1] == 2.0, "像元缓存：DONE 块才算完成、半截块忽略、缺测记录、续跑合并")
    # 5) 合成站：事件观测/模型与比值恢复
    n = 24 * 800
    h0 = 24 * 16000
    ser = {k: np.full(n, np.nan) for k in ("rain", "wind", "s1", "s5", "sss05", "pco2", "sst")}
    ser["wind"][:] = 5.0
    ser["rain"][:] = 0.0
    base = 35.0 + rng.normal(0, 0.003, n)
    true_ratio = 0.5
    onsets = list(range(24 * 40 + 5, n - 24 * 40, 24 * 11))
    cm = {}
    Pfull = np.zeros(2 * n)
    for i0 in onsets:
        for k in range(4):
            Pfull[2 * (i0 + k):2 * (i0 + k) + 2] = 8.0
            ser["rain"][i0 + k] = 8.0
    ser_s = {}
    for z, key in ((0.5, "sss05"), (1.0, "s1"), (5.0, "s5")):
        Uh = np.repeat(ser["wind"] * WIND_FACTOR, 2)
        F = np.ones(2 * n)
        F[48:] = rim_factor(Pfull, Uh, z)
        Fh = F.reshape(-1, 2).mean(axis=1)
        ser_s[key] = base * (1 - true_ratio * (1 - Fh))
    ser["s1"], ser["s5"] = ser_s["s1"], ser_s["s5"]
    ser["sss05"] = np.where(np.arange(n) % 3 == 0, ser_s["sss05"], np.nan)
    ser["sst"][:] = 28.0
    ser["pco2"] = np.where(np.arange(n) % 3 == 0, 400.0 + 10.0 * (ser["sss05"] - 35.0), np.nan)
    stn = {"name": "SYN", "group": GROUP_A, "h0": h0, "n": n, "ser": ser}
    stn2 = dict(stn, name="SYN2", group=GROUP_B)
    for i0 in range(n):
        cm[(h0 + i0)] = {"SYN": (Pfull[2 * i0], Pfull[2 * i0 + 1], 0, 0, 9, 9),
                         "SYN2": (Pfull[2 * i0], Pfull[2 * i0 + 1], 0, 0, 9, 9)}
    pixc = PixelCache(os.path.join(td, "px3.csv"))
    pixc.data = cm
    recs = []
    seasons = ["DJF", "MAM", "JJA", "SON"]
    for s_ in (stn, stn2):
        for j, i0 in enumerate(onsets):
            ctrl = [i0 + 24 * 5, i0 - 24 * 5]
            e = {"i": i0, "hour": h0 + i0, "onset_utc": p1.hour_to_iso(h0 + i0), "season": seasons[j % 4],
                 "acc24": 32.0, "n_ctrl": 2, "ctrl_hours": [h0 + c for c in ctrl],
                 "ctrl_ds1": [obs_delta(ser["s1"], c, PRIMARY_WIN) for c in ctrl], "s1_depth_m": 1.0}
            recs.append(event_record(s_, e, pixc, ("cmorph", "gauge"), [PRIMARY_WIN] + SENS_WINDOWS, k27=True))
    d3 = d3_gate(recs, {"SYN": stn, "SYN2": stn2})
    check(all(not r["obs"] for r in recs) and "z=0.5m" in d3 and "z=1.0m" in d3 and isinstance(d3["pass"], bool),
          f"D3 在事件观测填入前可算（{d3['verdict']}）")
    r0 = event_record(stn, dict(e, i=onsets[0], hour=h0 + onsets[0]), pixc, ("cmorph", "gauge"),
                      [PRIMARY_WIN] + SENS_WINDOWS, k27=False)
    r1 = event_record(stn, dict(e, i=onsets[0], hour=h0 + onsets[0]), pixc, ("cmorph", "gauge"),
                      [PRIMARY_WIN] + SENS_WINDOWS, k27=True)
    check({k: v for k, v in r1["model"].items() if k[0] != K27_TAG} == r0["model"]
          and all(isnum(r1["model"][(K27_TAG, "rim", z, PRIMARY_WIN)]) for z in K27_DEPTHS),
          "K27 开关不改主路径任何模型值，K27 另存独立键")
    saved = {id(r): {k: v for k, v in r["model"].items() if k[0] == K27_TAG} for r in recs}
    for r in recs:
        for k in saved[id(r)]:
            r["model"][k] = -999.0
    d3b = d3_gate(recs, {"SYN": stn, "SYN2": stn2})
    for r in recs:
        r["model"].update(saved[id(r)])
    check(d3b == d3, "K27 值被篡改时 D3 结果逐字不变（D3 不读 K27）")
    for r in recs:
        fill_obs(r, stn if r["station"] == "SYN" else stn2)
    m = block_metrics(recs, {"SYN": (10.0, 10.0), "SYN2": (10.0, 10.0)})
    rr = m["R_RIM_1m"]
    check(rr["evaluable"] and abs(rr["point"] - true_ratio) < 0.03 and rr["ci95"][0] < true_ratio < rr["ci95"][1] + 0.03,
          f"合成数据恢复 R_RIM(1 m)≈{true_ratio}（得 {rr['point']}，CI {rr['ci95']}）")
    rg = ratio_stats(recs, g_obs(1.0), g_rim(1.0, "gauge"), "gauge")
    check(abs(rg["point"] - true_ratio) < 0.03, f"雨量计驱动同样恢复（{rg['point']}）")
    p5 = m["P5"]
    check(abs(p5["point"] - p5["P5_RIM3"]["point"]) < 0.02, f"P5 观测≈模型（{p5['point']} vs {p5['P5_RIM3']['point']}）")
    b = m["beta_all_info"]
    check(b["point"] is not None and abs(b["point"] - 10.0) < 0.5, f"β_obs 恢复 10 µatm/psu（{b['point']}）")
    kb = k27_block(recs)
    k1, k05 = kb["R_RIM_1.0m"], kb["R_RIM_0.5m"]
    check("探索性" in kb["status"] and "不参与出口" in kb["status"] and k1["evaluable"] and k05["evaluable"]
          and "ci95" in k1 and "tost" not in k1 and "ci95_outside_sesoi" not in k1
          and k1["point"] <= rr["point"] + 1e-9 and k1["n_clusters"] == rr["n_clusters"]
          and set(kb["strata_by_group_points"]) == {GROUP_A, GROUP_B}
          and all(kb["strata_by_group_points"][g]["R_RIM_1.0m"]["point"] is not None for g in (GROUP_A, GROUP_B)),
          f"K27 汇总块：探索性标注、点估计＋CI、无 TOST/SESOI 标志、分层点估计（R_RIM,K27(1 m)={k1['point']} ≤ 主 {rr['point']}）")
    # 5b) 计数补丁：雨量计负值置 0 计数；CMORPH >200 mm/h → 模型 NaN 且计数
    ser["rain"][onsets[0] - 10] = -0.3
    Pg, _Ug, dg = event_forcing(stn, h0 + onsets[0], "gauge", pixc)
    check(dg["gauge_neg"] == 1 and np.all(Pg >= 0), "雨量计驱动：负值置 0 并计数（K21）")
    hot = h0 + onsets[0] - 3
    pixc.data[hot] = dict(pixc.data[hot], SYN=(250.0, 0.0, 0, 0, 9, 9))
    rh = event_record(stn, dict(e, i=onsets[0], hour=h0 + onsets[0]), pixc, ("cmorph",), [PRIMARY_WIN], k27=True)
    check(rh["diag"]["cmorph"]["p_gt_d0_range"] == 1 and rh["reason"]["cmorph"] is None
          and not isnum(rh["model"][("cmorph", "rim", 1.0, PRIMARY_WIN)]),
          "CMORPH 半步 >200 mm/h：d0 表外 → 模型 NaN（照 Witte），p_gt_d0_range 计数")
    # 6) 出口判定逻辑
    G = {"D1": True, "D2": True, "D4": True, "D5": True}
    dp, df = {"pass": True}, {"pass": False}

    def R(lo, hi, eq=False):
        return {"evaluable": True, "ci95": [lo, hi], "ci95_outside_sesoi": hi < SESOI_LO or lo > SESOI_HI,
                "ci95_inside_sesoi": lo >= SESOI_LO and hi <= SESOI_HI, "tost": {"equivalent": eq}}

    check(decide_exit(G, df, R(0.9, 1.1, True), R(0.9, 1.1))["exit"] == 3, "D3 不确定 → 出口 3")
    check(decide_exit(dict(G, D5=False), dp, R(0.2, 0.4), R(0.9, 1.1))["exit"] == 3, "D5 不过 → 出口 3")
    check(decide_exit(G, dp, R(0.2, 0.5), R(0.9, 1.1))["exit"] == 1, "R_RIM CI 全在界外 → 出口 1")
    check(decide_exit(G, dp, R(0.8, 1.2, True), R(0.1, 0.3))["exit"] == 1
          and decide_exit(G, dp, R(0.8, 1.2, True), R(0.1, 0.3))["also_satisfied"] == [1, 2], "出口 1、2 同时成立 → 主标 1、两者都报")
    check(decide_exit(G, dp, R(0.8, 1.2, True), R(0.7, 1.5))["exit"] == 2, "CI 在界内且 TOST → 出口 2")
    check(decide_exit(G, dp, R(0.5, 1.1), R(0.7, 1.5))["exit"] is None, "CI 跨边界 → 未落入出口")
    # 7) 可选依赖：netCDF 像元抽取、PyCO2SYS 斜率
    try:
        import netCDF4
        path = os.path.join(td, "cm.nc")
        with netCDF4.Dataset(path, "w") as ds:
            ds.createDimension("time", 2)
            ds.createDimension("lat", 11)
            ds.createDimension("lon", 20)
            ds.createVariable("time", "i4", ("time",))[:] = [3600 * 500000, 3600 * 500000 + 1800]
            ds.createVariable("lat", "f8", ("lat",))[:] = np.linspace(-5, 5, 11)
            ds.createVariable("lon", "f8", ("lon",))[:] = np.arange(20) * 18.0 + 9.0
            v = ds.createVariable("cmorph", "i2", ("time", "lat", "lon"), fill_value=np.int16(-999))
            v.scale_factor = np.float32(0.01)
            v.missing_value = np.int16(-999)
            arr = np.zeros((2, 11, 20))
            arr[0, 5, 9] = 3.0
            arr[1, 5, 9] = 4.5
            v[:] = arr
            v[1, 4, 8] = np.ma.masked
        exr = Extractor([{"name": "X", "lat": 0.2, "lon": -189.0}])   # 171°E → 格点 171（j=9）
        rows, why = exr.extract(path, 500000)
        c = rows["X"]
        check(why is None and abs(c[0] - 3.0) < 1e-6 and abs(c[1] - 4.5) < 1e-6 and c[4] == 9 and c[5] == 8,
              "netCDF 像元抽取（最近格点、经度环绕、掩码→NaN、3×3 计数）")
        check(exr.extract(path, 500001)[0] is None, "时间轴不符 → 记缺测")
    except ImportError:
        skipped.append("netCDF4 不可用：跳过像元抽取自测")
    try:
        fs, ps = beta_model(2300.0, 1990.0, 28.0, 35.0)
        check(5.0 < fs < 20.0 and abs(ps / fs - 1) < 0.01, f"β_model（PyCO2SYS，Witte 口径）={fs:.3f} µatm/psu（fCO2），pCO2 {ps:.3f}")
    except ImportError:
        skipped.append("PyCO2SYS 不可用：跳过 β_model 自测")
    # 8) 事件定义阈值切换不泄漏
    old = p1.RAIN_EVENT_MM
    fake = {"series_raw": {k: [NAN] * 3000 for k in ("rain", "wind", "s1", "s5", "sss05", "pco2")}, "n": 3000,
            "h0": 0, "name": "Z", "regime": "z", "lon": 0.0, "group": GROUP_A, "s1_depth": None}
    fr = fake["series_raw"]
    for i in range(3000):
        fr["rain"][i] = 0.0
        fr["wind"][i] = 5.0
        fr["s1"][i] = 35.0
        if i % 3 == 0:
            fr["sss05"][i] = 35.0
            fr["pco2"][i] = 400.0
    for i in (500, 1500):
        fr["rain"][i] = 7.0                                   # 7 mm：5 mm 阈值是事件、10 mm 不是
    fr["rain"][2500] = 25.0
    n5, n10, n20 = (len(events_at(fake, t)) for t in (5.0, 10.0, 20.0))
    check((n5, n10, n20) == (3, 1, 1) and p1.RAIN_EVENT_MM == old, f"阈值 5/10/20 mm 事件数 {(n5, n10, n20)}，常量已复原")
    # 9) K28：「1 m」层模型深度＝事件 s1_depth_m；缺失/越界报错
    eb = dict(e, i=onsets[0], hour=h0 + onsets[0], s1_depth_m=1.5, s1_substitute=True)
    rb = event_record(stn2, eb, pixc, ("cmorph",), [PRIMARY_WIN], k27=True)
    ra = event_record(stn2, dict(eb, s1_depth_m=1.0, s1_substitute=False), pixc, ("cmorph",), [PRIMARY_WIN], k27=True)
    Pb, Ub, _db = event_forcing(stn2, eb["hour"], "cmorph", pixc)

    def by_hand(z, t=T_CURRENT_S):
        return rb["s0"][1.0] * win_delta(rim_factor(Pb, Ub, z, t).reshape(-1, 2).mean(axis=1), PRIMARY_WIN)
    W = PRIMARY_WIN
    check(rb["zmodel"] == {0.5: 0.5, 1.0: 1.5, 5.0: 5.0} and rb["model"][("cmorph", "rim", 1.0, W)] == by_hand(1.5)
          and rb["model"][("cmorph", "rim", 1.0, W)] != ra["model"][("cmorph", "rim", 1.0, W)]
          and ra["model"][("cmorph", "rim", 1.0, W)] == by_hand(1.0)
          and rb["model"][("cmorph", "rim", 0.5, W)] == ra["model"][("cmorph", "rim", 0.5, W)]
          and rb["model"][("cmorph", "rim", 5.0, W)] == ra["model"][("cmorph", "rim", 5.0, W)]
          and rb["model"][(K27_TAG, "rim", 1.0, W)] == by_hand(1.5, K27_T_CUR_DEPTH_S)
          and event_row(dict(rb, obs={(z, W): NAN for z in DEPTHS}))["model_depth_s1_layer_m"] == 1.5,
          f"K28：O 组事件「1 m」层按 s1_depth_m=1.5 算（{rb['model'][('cmorph', 'rim', 1.0, W)]:.5f} vs 1.0 m "
          f"{ra['model'][('cmorph', 'rim', 1.0, W)]:.5f}），0.5/5 m 不变，K27 同步，CSV 记深度")
    bad = []
    for v in ("__del__", NAN, None, 0.0, 2.0):
        e_bad = {k: x for k, x in eb.items() if k != "s1_depth_m"} if v == "__del__" else dict(eb, s1_depth_m=v)
        try:
            event_record(stn2, e_bad, pixc, ("cmorph",), [PRIMARY_WIN])
            bad.append(v)
        except ValueError:
            pass
    check(not bad, f"K28：s1_depth_m 缺失/NaN/None/0/2.0 m 均报错，不回退 1.0（未报错：{bad}）")
    # 10) K29 可用率函数
    mk = [dict(ra, model={("cmorph", "rim", 1.0, W): 1.0, ("cmorph", "s20", 0.0, W): 1.0}, reason={"cmorph": None})
          for _ in range(2)]
    mk += [dict(ra, model={("cmorph", "rim", 1.0, W): NAN, ("cmorph", "s20", 0.0, W): NAN},
                reason={"cmorph": "cmorph_missing"}) for _ in range(2)]
    mk += [dict(ra, model={("cmorph", "rim", 1.0, W): NAN, ("cmorph", "s20", 0.0, W): 1.0}, reason={"cmorph": None},
                diag={"cmorph": {"p_gt_d0_range": 1}})]
    a4, a5 = model_availability(mk[:4]), model_availability(mk)
    check(a4["pass"] and a4["fraction"] == 0.5 and not a5["pass"] and a5["n_available"] == 2
          and a5["missing_reason_counts"] == {"cmorph_missing": 2, "d0_out_of_range_gt_200mmh": 1},
          "K29 可用率：恰 50% 通过、<50% 不通过，缺测原因分类计数")
    # 11) 端到端（合成数据，假下载）：all／primary／sensitivity／K29 中止／参数保护
    g = globals()
    fx = _e2e_fixture(td, z_b=1.5)
    outs = {k: os.path.join(td, "e2e_" + k) for k in ("all", "prim", "sens", "sens2", "sens_nores", "k29")}

    def jl(d, nm):
        with open(os.path.join(d, nm), encoding="utf-8") as f:
            return json.load(f)

    def has(d, nm):
        return os.path.exists(os.path.join(d, nm))
    rc_all = _e2e_run(g, fx, outs["all"], None)
    S_all = jl(outs["all"], "p2_summary.json")
    rr_all = S_all["primary"]["R_RIM_1m"]
    with open(os.path.join(outs["all"], "p2_events.csv"), encoding="utf-8") as f:
        rows_all = list(csv.DictReader(f))
    st = S_all["strata_by_group"]
    check(rc_all == 0 and has(outs["all"], "p2_d3.json") and not has(outs["all"], "p2_abort.json")
          and S_all["part"] == "all" and "window_0-3h" in S_all["sensitivity"] and "threshold_5mm" in S_all["sensitivity"]
          and abs(rr_all["point"] - 0.5) < 0.03
          and all(abs(st[gp]["R_RIM_1m"]["point"] - 0.5) < 0.03 for gp in (GROUP_A, GROUP_B))
          and "wind_floor_hours" in rows_all[0] and "wind_floor_halfsteps" not in rows_all[0]
          and {r["model_depth_s1_layer_m"] for r in rows_all} == {"1.0", "1.5"}
          and rows_all[0]["dS_rim_1_gauge"] != ""
          and S_all["k28_R_RIM_1m_by_layer"]["substitute_layer"]["n_ratio"] > 0
          and abs(S_all["k28_R_RIM_1m_by_layer"]["substitute_layer"]["point"] - 0.5) < 0.03
          and "hours_floored" in S_all["wind"] and S_all["model_availability_k29"]["pass"],
          f"端到端 all：rc=0、summary/d3/敏感性齐全；B 组观测在 1.5 m 生成时两群 R_RIM(1 m) 都恢复 0.5"
          f"（全 {rr_all['point']}，A {st[GROUP_A]['R_RIM_1m']['point']}，B {st[GROUP_B]['R_RIM_1m']['point']}）")
    fx["downloaded"].clear()
    rc_p = _e2e_run(g, fx, outs["prim"], "primary")
    S_p = jl(outs["prim"], "p2_summary.json")
    with open(os.path.join(outs["prim"], "p2_events.csv"), encoding="utf-8") as f:
        rows_p = list(csv.DictReader(f))
    dl_p = sorted(h for x in fx["downloaded"] for h in x)
    check(rc_p == 0 and has(outs["prim"], "p2_d3.json") and not has(outs["prim"], "p2_sensitivity.json")
          and dl_p == fx["prim"] and "status" in S_p["sensitivity"] and "window_0-3h" not in S_p["sensitivity"]
          and all(r["dS_rim_1_gauge"] == "" for r in rows_p)
          and S_p["primary"] == S_all["primary"] and S_p["exit"] == S_all["exit"]
          and S_p["strata_by_group"] == S_all["strata_by_group"] and S_p["exploratory_k27"] == S_all["exploratory_k27"]
          and jl(outs["prim"], "p2_d3.json")["z=1.0m"] == jl(outs["all"], "p2_d3.json")["z=1.0m"]
          and [{k: v for k, v in r.items() if "gauge" not in k} for r in rows_p]
          == [{k: v for k, v in r.items() if "gauge" not in k} for r in rows_all],
          f"--part primary：只下主集 {len(dl_p)} 小时、不跑敏感性、主结果/出口/K27/D3/事件 CSV（雨量计列除外）与 all 逐值相同")
    rc_nr = _e2e_run(g, fx, outs["sens_nores"], "sensitivity")
    check(rc_nr == 3 and not has(outs["sens_nores"], "p2_sensitivity.json")
          and not has(outs["sens_nores"], "p2_summary.json"),
          "--part sensitivity 无 --resume-from（主集 D3 未写盘）→ 退出 3，不产出")
    fx["downloaded"].clear()
    rc_s = _e2e_run(g, fx, outs["sens"], "sensitivity", resume=outs["prim"])
    S_s = jl(outs["sens"], "p2_sensitivity.json")
    dl_s = sorted(h for x in fx["downloaded"] for h in x)
    check(rc_s == 0 and not has(outs["sens"], "p2_summary.json") and not has(outs["sens"], "p2_d3.json")
          and has(outs["sens"], "p2_d3_from_primary.json") and has(outs["sens"], "p2_plan_sensitivity.json")
          and dl_s == fx["extra"] and S_s["sensitivity"] == S_all["sensitivity"] and "exit" not in S_s
          and S_s["primary_d3_reference"]["pass"] == jl(outs["prim"], "p2_d3.json")["pass"]
          and not has(outs["prim"], "p2_sensitivity.json") and not has(outs["prim"], "p2_d3_from_primary.json"),
          f"--part sensitivity（续主集）：只补下 {len(dl_s)} 个敏感性独有小时，不写 summary/d3，四项敏感性与 all 逐值相同，"
          "不改主集目录")
    fx["downloaded"].clear()
    rc_s2 = _e2e_run(g, fx, outs["sens2"], "sensitivity", resume=outs["sens"])
    check(rc_s2 == 0 and sum(len(x) for x in fx["downloaded"]) == 0
          and jl(outs["sens2"], "p2_sensitivity.json")["sensitivity"] == S_s["sensitivity"],
          "--part sensitivity 自身续跑（续上次敏感性目录）：零下载、结果相同")
    rc_bad = _e2e_run(g, fx, os.path.join(td, "e2e_bad"), "primary", resume=os.path.join(td, "不存在"))
    check(rc_bad == 3, "--resume-from 目录不存在 → 退出 3")
    ev10 = fx["ev"][PRIMARY_THRESHOLD]
    fx["missing_hours"] = needed_hours(sorted(ev10, key=lambda x: x["hour"])[:int(0.6 * len(ev10))])
    for pt in ("primary", None):
        od = outs["k29"] + (pt or "all")
        rc_k = _e2e_run(g, fx, od, pt)
        ab = jl(od, "p2_abort.json") if has(od, "p2_abort.json") else {}
        check(rc_k == 3 and ab.get("status") == "aborted_K29" and ab["n_available"] < 0.5 * ab["n_total"]
              and ab["missing_reason_counts"].get("cmorph_missing", 0) > 0
              and not has(od, "p2_summary.json") and not has(od, "p2_d3.json") and not has(od, "p2_events.csv"),
              f"K29（--part {pt or 'all'}）：可用 {ab.get('n_available')}/{ab.get('n_total')} <50% → 退出 3、写 p2_abort.json、"
              "无 d3/summary/事件 CSV")
    fx["missing_hours"] = set()
    # 12) build_all 在下载前校验全部事件 s1_depth_m（K28 早失败），并在 plan 里记深度计数
    import types
    saved = {k: g[k] for k in ("rebuild_stations", "events_at", "validate_rebuild")}
    evmap = {PRIMARY_THRESHOLD: [dict(x) for x in fx["ev"][PRIMARY_THRESHOLD]], 5.0: [dict(x) for x in fx["ev"][5.0]],
             20.0: [dict(x) for x in fx["ev"][20.0]]}
    g["rebuild_stations"] = lambda http, cache, log: fx["stations"]
    g["events_at"] = lambda stn_, thr, log=None: [dict(x) for x in evmap[thr] if x["_stn"] == stn_["name"]]
    g["validate_rebuild"] = lambda *a: {"d4_count": 0}
    ba = types.SimpleNamespace(p1b_dir=fx["p1bdir"], p1_events=fx["p1csv"])
    try:
        bo = os.path.join(td, "e2e_build")
        os.makedirs(bo, exist_ok=True)
        with open(os.devnull, "w") as dn:
            import contextlib
            with contextlib.redirect_stdout(dn):
                plan_ok = build_all(ba, bo, _QuietLog())[-1]
                evmap[5.0][-1].pop("s1_depth_m")
                try:
                    build_all(ba, bo, _QuietLog())
                    raised = False
                except ValueError:
                    raised = True
    finally:
        g.update(saved)
    kc = plan_ok["k28_model_depth_counts"]
    check(raised and kc["10mm"] == {"SYNA@1.000": 66, "SYNB@1.500": 66} and "parts" in plan_ok
          and plan_ok["parts"]["sensitivity_beyond_primary"]["unique_hours"] == len(fx["extra"]),
          f"build_all：下载前校验 s1_depth_m（5 mm 集缺一条即报错），plan 记深度计数 {kc['10mm']} 与分段小时数")
    shutil.rmtree(td, ignore_errors=True)
    for msg in ok:
        print("PASS", msg)
    for msg in skipped:
        print("SKIP", msg)
    print(f"selftest：{len(ok)} 项通过，{len(skipped)} 项跳过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
