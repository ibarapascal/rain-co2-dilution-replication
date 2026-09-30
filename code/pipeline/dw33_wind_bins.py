#!/usr/bin/env python3
"""C3（阶段 dw33）：全球 RIM-3 界面／稀释通量按像元 U10 细档累计，
供「系泊三档比值逐档缩放」的全球 U（后处理 dw33_post.py）。

只含所需累计量：整体缩放 s＝0、s＝1 两个情形＋一个逐像元检查情形（s 按像元风速所在系泊档取三档点估计，L25 展开），
全球按细档（主：当前半步 U；敏感性：最近 12 个半步 U 均值）累计；不算 s 网格、历史项缩放、S20、z5、G、纬带／洋盆、格点场。
引擎＝import p3_global（不改一字），子类化 Engine 只重写逐步累计（逐式照 Engine.step／_chunk）；
静态场只读沿用 P3 的 p3-work/static（先核 version／code_sha256／N）。

用法（产物写 --out，未给则 $REPRO_OUTPUT_DIR）：
  dw33_wind_bins.py --mode bench --cache-root C --static-dir S --work-dir W
        2000-01-01 一天（含起转）：本引擎 → 原 Engine，同日全局合计逐项比对（相对差 ≤1e-9，不过退出 4），
        报每步耗时、峰值 RSS、全年外推 → dw33_bench.json（同时写 1 天的 dw33_bins.json 供后处理测通）
  dw33_wind_bins.py --mode full --cache-root C --static-dir S --work-dir W --workers 2
        按月并行（每月前补 24 h 起转）→ W/month_YYYYMM.{npz,json}（可续跑）→ dw33_bins.json、dw33_run.json
退出码：0 完成；2 数据缺失；3 其他异常（含静态场核对不符）；4 bench 等价比对不过。
依赖：同 p3_global（numpy、scipy、netCDF4、gsw、PyCO2SYS）。
实现选择：
  W1 细档右闭：bin＝searchsorted(UE, U, 'left')，UE＝细档有限边界去 0；U≤1 入 0 档；U>30 入末档；NaN 另档。
  W2 系泊档由细档上界判定（上界 ≤6.031179→T1，≤9.325534→T2，否则 T3），与逐像元检查情形用 U 直接判定一致（边界在细档边上）。
  W3 6 h 均值：float32 环形 12 槽（含起转），逐步求均值；任一槽 NaN→NaN 档。只用于分档，不进通量。
  W4 逐像元检查情形只在主分档累计；s 由 U 直接查三档（NaN U→T3，其通量本为 NaN、不入配对和）。
  W5 月结果 npz 只存 bacc（量×分档变体×细档）；续跑判据＝参数哈希＋天数（同 P3 L40）。
Change Log:
  2026-09-29：初版。
"""
import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import argparse
import concurrent.futures as cf
import json
import multiprocessing as mp
import sys
import time
import traceback
from datetime import date, datetime, timedelta

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p3_global as g  # noqa: E402

VERSION = "dw33-2026-09-29a"
import repro_paths as _rp  # [repro] 路径集中配置（code/config.*.toml）
# [repro] 核 import 的 p3_global.py 等于 pipeline/ 里当前文件的 sha256（p3-global 阶段写的静态场 meta.code_sha256 也是它），
# [repro] 保证静态场与本引擎来自同一份 p3_global.py
P3_CODE_SHA = _rp.module_sha("p3_global.py")
P3_VERSION = "p3-2026-09-27b"
STATIC_N = 5884455
CUT1, CUT2 = 6.031179, 9.325534
R_T = (0.37105, 0.26521, 0.15673)                          # P4S D_U.cells[k].point（逐像元检查情形）
UE = np.array([1, 1.641173, 2, 2.581536, 3, 4, 5, 6, 6.031179, 7, 8, 9, 9.325534, 10, 11, 12, 13, 14, 15, 15.071225,
               16, 18, 20, 22, 24, 26, 26.633166, 30], float)
NB = len(UE) + 2                                            # 细档 0..len(UE)（末档 >30）＋ NaN 档
NAN_BIN = NB - 1
TERT = np.array([0 if (k < len(UE) and UE[k] <= CUT1) else 1 if (k < len(UE) and UE[k] <= CUT2) else 2
                 for k in range(NB - 1)] + [-1], int)
QS = ["n_W", "W_all", "skin", "n_D", "dep", "turb", "int0", "dil0", "int1", "dil1", "nm_int1", "intP", "dilP"]
QI = {q: i for i, q in enumerate(QS)}
VARIANTS = ("cur", "mean6h")


def bins_of(U):
    U = np.asarray(U, float)
    b = np.searchsorted(UE, U, side="left")
    b[~np.isfinite(U)] = NAN_BIN
    return b.astype(np.int32)


def s_of(U):
    return np.where(U <= CUT1, R_T[0], np.where(U <= CUT2, R_T[1], R_T[2]))


class BinEngine(g.Engine):
    """逐步：稠密风驱与 L30 常数项、D 像元上沉降／湍流／s＝0、1 与逐像元 s 的稀释／界面配对和，按风速细档累计。"""

    def __init__(self, dom, chunk=g.CHUNK):
        super().__init__(dom, chunk=chunk, keep=False)
        self.bacc = np.zeros((len(QS), len(VARIANTS), NB))
        self.uring = np.full((12, dom.N), np.nan, np.float32)
        self.un = 0
        self.b = [None, None]

    def _push_u(self, U):
        self.uring[self.un % 12] = np.asarray(U, np.float32)
        self.un += 1

    def _addb(self, q, x, variants=(0, 1), idx=None):
        x = np.asarray(x, float)
        for v in variants:
            b = self.b[v] if idx is None else self.b[v][idx]
            self.bacc[QI[q], v] += np.bincount(b, weights=x, minlength=NB)

    def push_only(self, P, U):
        super().push_only(P, U)
        self._push_u(U)

    def step(self, P, U):
        dom, N = self.dom, self.N
        P = np.asarray(P, float)
        U = np.asarray(U, float)
        t0 = time.perf_counter()
        self._push_u(U)
        u6 = np.zeros(N)
        with np.errstate(all="ignore"):
            if self.un >= 12:
                for i in range(12):
                    u6 += self.uring[i]
                u6 /= 12.0
            else:
                u6[:] = np.nan
        self.b = [bins_of(U), bins_of(u6)]
        del u6
        t0 = self._t("bins", t0)
        ent = self.ring.entries(P, U)
        logH, logH5, Dmask = self.ring.history()
        t0 = self._t("history", t0)
        Dmask[ent["pos"]] = True
        c = self.co
        A = dom.area
        with np.errstate(all="ignore"):
            kw = c["g_kw"] * (U ** 2)                                 # 同 Engine.step（TB L183）
            Trw = 12 * (kw * g.TCONV_K)
            Fwind = Trw * self.dC0 * self.icef
            okW = np.isfinite(Fwind)
            Fdil0 = Trw * self.dC0d * self.icef
            nonD = ~Dmask
            okS = nonD & okW & np.isfinite(Fdil0)
            self._addb("n_W", okW)
            self._addb("W_all", np.where(okW, A * Fwind, 0.0))
            self._addb("skin", np.where(okS, A * (Fdil0 - Fwind), 0.0))
        t0 = self._t("dense", t0)
        D = np.flatnonzero(Dmask)
        cz0 = np.ones(len(D))
        k = np.searchsorted(D, ent["pos"])
        cz0[k] = ent["cur0"]
        for s in range(0, len(D), self.chunk):
            self._chunk_b(D[s:s + self.chunk], cz0[s:s + self.chunk], P, U, logH, Trw, Fwind)
        self.ring.push(ent)
        self.diag["D_frac"].append(len(D) / max(N, 1))
        self.diag["steps"] += 1

    def _chunk_b(self, idx, cz0, P, U, logH, Trw_f, Fwind_f):
        dom, c = self.dom, self.co
        t0 = time.perf_counter()
        A = dom.area[idx]
        Fw = Fwind_f[idx]
        okW = np.isfinite(Fw)
        S0, T, icef = self.S0[idx], self.T[idx], self.icef[idx]
        sf, fa, beta = self.sf[idx], self.fa[idx], dom.beta[idx]
        aw0, aa0d, dC0 = c["aw0"][idx], c["aa0d"][idx], self.dC0[idx]
        p, u, Trw = P[idx], U[idx], Trw_f[idx]
        lon, lat = dom.lon1d[dom.col[idx]], dom.lat1d[dom.row[idx]]
        with np.errstate(all="ignore"):
            needH = ~((p == 0) & (u > 0))
            TrH = Trw.copy()
            if needH.any():
                TrH[needH] = 12 * (g.k_harrison(c["sck"][idx][needH], u[needH], p[needH]) * g.TCONV_K)
            Fturb = TrH * dC0 * icef
            Fdep = 12 * (p * g.TCONV_R) * c["af"][idx] * fa
            okT = okW & np.isfinite(Fturb)
            okP = okW & np.isfinite(Fdep)
            self._addb("n_D", okW, idx=idx)
            self._addb("turb", np.where(okT, A * (Fturb - Fw), 0.0), idx=idx)
            self._addb("dep", np.where(okP, A * Fdep, 0.0), idx=idx)
            S_R = S0 * (np.exp(logH[idx]) * cz0)
            okR = np.isfinite(S_R)
            missR = ~okR
        t0 = self._t("d_base", t0)
        aw_R = np.full(len(idx), np.nan)
        aa_R = np.full(len(idx), np.nan)
        if okR.any():
            aw_R[okR], aa_R[okR] = g.sol_pair(T[okR], S_R[okR], lon[okR], lat[okR])
        t0 = self._t("d_exact", t0)
        co = {kk: c[kk][idx] for kk in ("aw0", "aa0d", "Kw1", "Kw2", "Ka1", "Ka2")}
        for tag in ("0", "1", "P"):
            with np.errstate(all="ignore"):
                if tag == "1":
                    S, aw, aa = S_R, aw_R, aa_R
                elif tag == "0":
                    S, aw, aa = S0, aw0, aa0d
                else:
                    S = S0 + s_of(u) * (S_R - S0)
                    S = np.where(S < 0, 0.0, S)                         # L28（s≤1 时不会触发）
                    aw, aa, _fb = g.alpha_expand(S, S0, co, T, lon, lat)
                dS = S - S0
                Xc = ((sf + beta * dS) * aw) - (fa * aa)
                Fdil = Trw * Xc * icef
                Fint = TrH * Xc * icef
                okd = okW & np.isfinite(Fdil) & ~missR
                oki = okW & np.isfinite(Fint) & ~missR
                vs = (0,) if tag == "P" else (0, 1)
                self._addb(f"dil{tag}", np.where(okd, A * (Fdil - Fw), 0.0), vs, idx)
                self._addb(f"int{tag}", np.where(oki, A * (Fint - Fw), 0.0), vs, idx)
                if tag == "1":
                    self._addb("nm_int1", okW & ~oki, vs, idx)
        self._t("d_cases", t0)


def check_static(sd):
    with open(os.path.join(sd, "meta.json"), encoding="utf-8") as f:
        meta = json.load(f)
    bad = [k for k, want in (("version", P3_VERSION), ("code_sha256", P3_CODE_SHA), ("N", STATIC_N)) if meta.get(k) != want]
    if bad:
        raise RuntimeError(f"静态场核对不符：{bad}")
    return meta


def run_month(cfg, engine_cls, max_days=None):
    """一个月（或前 max_days 天）＋前一日起转；返回引擎与计时。"""
    ym = cfg["ym"]
    y, m = divmod(ym, 100)
    lay = g.Layout(cfg["cache_root"])
    g.wait_ready(lay, ym, cfg["wait_s"], need_next=not max_days)
    st, meta = g.load_static(cfg["static_dir"])
    dom = g.domain_from_static(st, meta)
    eng = engine_cls(dom, chunk=cfg["chunk"])
    fld = g.Fields(lay, dom)
    first = date(y, m, 1)
    prev = first - timedelta(days=1)
    t_step = 0.0
    nstep = 0
    for hh in range(24):
        dt = datetime(prev.year, prev.month, prev.day, hh, tzinfo=g.UTC)
        P0, P1 = fld.cmorph(dt)
        U0, U30 = fld.wind_pair(g.epoch_hour(dt))
        eng.push_only(P0, U0)
        eng.push_only(P1, U30)
    eng.set_month(*g.watson_month(st, m), pco2=g.watson_is_pco2(meta))
    days = g.days_of_month(y, m)
    if max_days:
        days = days[:max_days]
    for d in days:
        sst, ice, s0 = fld.day(d)
        eng.set_day(sst, ice, s0)
        del sst, ice, s0
        for hh in range(24):
            dt = datetime(d.year, d.month, d.day, hh, tzinfo=g.UTC)
            P0, P1 = fld.cmorph(dt)
            U0, U30 = fld.wind_pair(g.epoch_hour(dt))
            for P, U in ((P0, U0), (P1, U30)):
                ts = time.perf_counter()
                eng.step(P, U)
                t_step += time.perf_counter() - ts
                nstep += 1
        g.log(f"[{ym}] {d} 完成（{engine_cls.__name__}，累计 {nstep} 步，峰值 RSS {g.peak_rss_bytes() / 2**30:.2f} GiB）")
    fld.close()
    return eng, {"days": len(days), "steps": nstep, "step_seconds": t_step, "fields": fld.diag, "meta": meta}


def process_month(cfg):
    ym = cfg["ym"]
    g._LOG_PATH = os.path.join(cfg["work_dir"], f"log_{ym}.txt")
    base = os.path.join(cfg["work_dir"], f"month_{ym}")
    try:
        if os.path.exists(base + ".json") and os.path.exists(base + ".npz"):
            with open(base + ".json", encoding="utf-8") as f:
                old = json.load(f)
            y, m = divmod(ym, 100)
            if old.get("config_hash") == cfg["config_hash"] and old.get("days") == len(g.days_of_month(y, m)):
                g.log(f"[{ym}] 已有月结果，续用")
                return {"ym": ym, "ok": True, "resumed": True}
        t0 = time.perf_counter()
        eng, info = run_month(cfg, BinEngine)
        out = {"ym": ym, "version": VERSION, "config_hash": cfg["config_hash"], "days": info["days"], "steps": info["steps"],
               "seconds": round(time.perf_counter() - t0, 1), "step_seconds": round(info["step_seconds"], 1),
               "peak_rss_bytes": g.peak_rss_bytes(), "timing_s": {k: round(v, 2) for k, v in eng.tim.items()},
               "D_frac_mean": float(np.mean(eng.diag["D_frac"])) if eng.diag["D_frac"] else None,
               "cmorph_missing": info["fields"]["cmorph_missing"], "wind_persist_last": info["fields"]["wind_persist_last"]}
        if os.path.exists(base + ".json"):
            os.remove(base + ".json")
        np.savez(base + ".npz", bacc=eng.bacc)
        g.jdump(out, base + ".json")
        return {"ym": ym, "ok": True, "resumed": False}
    except g.DataError as e:
        g.log(f"[{ym}] 数据错误：{e}")
        return {"ym": ym, "ok": False, "kind": "data", "error": str(e)}
    except Exception as e:
        g.log(f"[{ym}] 异常：{e}\n{traceback.format_exc()}")
        return {"ym": ym, "ok": False, "kind": "error", "error": f"{type(e).__name__}: {e}"}


def bins_doc(bacc, extra):
    return {"version": VERSION, "p3_code_sha256": P3_CODE_SHA, "upper_edges": UE.tolist(), "n_bins": NB, "nan_bin": NAN_BIN,
            "tercile_of_bin": TERT.tolist(), "cuts": [CUT1, CUT2], "s_check_by_tercile": list(R_T), "variants": list(VARIANTS),
            "quantities": QS, "units": "gC（面积×通量逐半步求和；n_* 为像元×半步计数）；百分比口径由后处理按 P3 定义换算",
            "bacc": {v: {q: bacc[QI[q], vi].tolist() for q in QS} for vi, v in enumerate(VARIANTS)}, **extra}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("bench", "full"), required=True)
    ap.add_argument("--cache-root", required=True)
    ap.add_argument("--static-dir", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--chunk", type=int, default=g.CHUNK)
    ap.add_argument("--wait-hours", type=float, default=1.0)
    a = ap.parse_args()
    out_dir = a.out or os.environ.get("REPRO_OUTPUT_DIR")
    if not out_dir:
        print("需要 $REPRO_OUTPUT_DIR 或 --out", file=sys.stderr)
        return 3
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(a.work_dir, exist_ok=True)
    g._LOG_PATH = os.path.join(out_dir, "dw33_log.txt")
    code_sha = g.sha256_file(os.path.abspath(__file__))
    p3_sha = g.sha256_file(os.path.abspath(g.__file__))
    g.log(f"{VERSION}（sha256 {code_sha[:12]}…；p3_global {p3_sha[:12]}…）mode={a.mode}")
    try:
        if p3_sha != P3_CODE_SHA:
            raise RuntimeError(f"p3_global.py sha256 {p3_sha} ≠ 预期 {P3_CODE_SHA}")
        check_static(a.static_dir)
        chash = g.config_hash({"version": VERSION, "code_sha256": code_sha, "chunk": a.chunk,
                               "cache_root": os.path.abspath(a.cache_root), "static_N": STATIC_N})
        base_cfg = {"cache_root": os.path.abspath(a.cache_root), "static_dir": a.static_dir, "work_dir": a.work_dir,
                    "chunk": a.chunk, "wait_s": a.wait_hours * 3600, "config_hash": chash}
        if a.mode == "bench":
            cfg = dict(base_cfg, ym=200001)
            t0 = time.perf_counter()
            eng, info = run_month(cfg, BinEngine, max_days=1)
            rss_new = g.peak_rss_bytes()
            wall_new = time.perf_counter() - t0
            tim_new = dict(eng.tim)
            bacc = eng.bacc.copy()
            del eng
            t0 = time.perf_counter()
            ref, info_ref = run_month(cfg, g.Engine, max_days=1)
            wall_ref = time.perf_counter() - t0
            qi = {q: i for i, q in enumerate(ref.qn)}
            T = lambda q: float(ref.acc[qi[q]].sum())
            B0 = lambda q: float(bacc[QI[q], 0].sum())
            B1 = lambda q: float(bacc[QI[q], 1].sum())
            pairs = [("W_all", "W_all"), ("skin", "skin"), ("dep", "dep"), ("turb", "turb"), ("n_W", "n_W"), ("n_D", "n_D"),
                     ("int1", "int|whole_1"), ("dil1", "dil|whole_1"), ("int0", "int|whole_0"), ("dil0", "dil|whole_0")]
            cmp_ = []
            for mine, theirs in pairs:
                v0, v1, r = B0(mine), B1(mine), T(theirs)
                rel = abs(v0 - r) / max(abs(r), 1e-30)
                rel6 = abs(v1 - r) / max(abs(r), 1e-30)
                cmp_.append({"q": mine, "p3": theirs, "ours_cur": v0, "ours_mean6h": v1, "p3_value": r, "rel_cur": rel,
                             "rel_mean6h": rel6, "pass": bool(rel <= 1e-9 and rel6 <= 1e-9)})
            ok = all(c["pass"] for c in cmp_)
            spd = info["step_seconds"] / max(info["steps"], 1)
            day_over = (wall_new - info["step_seconds"]) / 1.0            # 起转＋读场＋日系数（1 天）
            month_s = [len(g.days_of_month(2000, m)) * (48 * spd) + day_over * len(g.days_of_month(2000, m)) for m in range(1, 13)]
            rep = {"version": VERSION, "code_sha256": code_sha, "p3_code_sha256": p3_sha, "equivalence_pass": ok, "compare": cmp_,
                   "ours": {"wall_s": round(wall_new, 1), "steps": info["steps"], "sec_per_step": spd, "peak_rss_bytes": rss_new,
                            "timing_s": {k: round(v, 2) for k, v in tim_new.items()}},
                   "p3_engine": {"wall_s": round(wall_ref, 1), "sec_per_step": info_ref["step_seconds"] / max(info_ref["steps"], 1),
                                 "peak_rss_after_both": g.peak_rss_bytes()},
                   "extrapolation": {"serial_hours": sum(month_s) / 3600,
                                     **{f"workers_{w}_hours": g.lpt_makespan(month_s, w) / 3600 for w in (1, 2, 3, 4)}},
                   "note": "外推＝每天 48 步×实测每步＋（bench 墙钟−步时）按天计；未计静态场（沿用，0）与月间并行的 I/O 竞争"}
            g.jdump(rep, os.path.join(out_dir, "dw33_bench.json"))
            g.jdump(bins_doc(bacc, {"mode": "bench", "days": "2000-01-01 only", "code_sha256": code_sha}),
                    os.path.join(out_dir, "dw33_bins.json"))
            g.log(f"bench：等价 {'过' if ok else '不过'}；每步 {spd:.2f} s；峰值 RSS {rss_new / 2**30:.2f} GiB；"
                  f"2 worker 外推 {rep['extrapolation']['workers_2_hours']:.2f} h")
            return 0 if ok else 4
        months = [200000 + m for m in range(1, 13)]
        cfgs = [dict(base_cfg, ym=ym) for ym in months]
        t0 = time.perf_counter()
        results = []
        ctx = mp.get_context("spawn")
        with cf.ProcessPoolExecutor(max_workers=max(1, a.workers), mp_context=ctx, max_tasks_per_child=1) as ex:
            futs = {ex.submit(process_month, c): c["ym"] for c in cfgs}
            for fu in cf.as_completed(futs):
                try:
                    r = fu.result()
                except Exception as e:
                    r = {"ym": futs[fu], "ok": False, "kind": "error", "error": f"worker 异常退出：{type(e).__name__}: {e}"}
                results.append(r)
                g.log(f"月 {r['ym']}：{'完成' if r['ok'] else '失败 ' + r.get('error', '')}")
        run = {"version": VERSION, "code_sha256": code_sha, "p3_code_sha256": p3_sha, "config_hash": chash,
               "results": sorted(results, key=lambda r: r["ym"]), "wall_seconds": round(time.perf_counter() - t0, 1)}
        bad = [r for r in results if not r["ok"]]
        if bad:
            g.jdump(run, os.path.join(out_dir, "dw33_run.json"))
            g.log(f"{len(bad)} 个月失败，不汇总（原样重提可续跑）")
            return 2 if all(r.get("kind") == "data" for r in bad) else 3
        bacc = np.zeros((len(QS), len(VARIANTS), NB))
        infos = {}
        for ym in months:
            bacc += np.load(os.path.join(a.work_dir, f"month_{ym}.npz"))["bacc"]
            with open(os.path.join(a.work_dir, f"month_{ym}.json"), encoding="utf-8") as f:
                infos[ym] = json.load(f)
        run["month_info"] = infos
        g.jdump(bins_doc(bacc, {"mode": "full", "months": months, "code_sha256": code_sha,
                                "cmorph_missing_halfhours": sum(len(i["cmorph_missing"]) * 2 for i in infos.values()),
                                "wind_persist_last": sum(i["wind_persist_last"] for i in infos.values())}),
                os.path.join(out_dir, "dw33_bins.json"))
        g.jdump(run, os.path.join(out_dir, "dw33_run.json"))
        g.log("完成：dw33_bins.json 已写")
        return 0
    except g.DataError as e:
        g.log(f"数据错误：{e}")
        return 2
    except Exception as e:
        g.log(f"异常：{type(e).__name__}: {e}\n{traceback.format_exc()}")
        return 3


if __name__ == "__main__":
    sys.exit(main())
