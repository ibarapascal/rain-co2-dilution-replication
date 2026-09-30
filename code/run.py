#!/usr/bin/env python3
"""run.py — single entry point of the pipeline: expands the commands of stages.py with the central configuration and runs
the stages in pipeline/ with uniform logging and run records.

Usage
  python3 run.py --config CFG --list                     list stages, dependencies, network/token needs, completion status
  python3 run.py --config CFG --stage p2 [--mode plan]   run one stage (default mode in stages.py)
  python3 run.py --config CFG --all                      run all in_all stages in topological order; completed stages are skipped
  python3 run.py --config CFG --stage p2 --dry-run       print the expanded commands and environment only
  --resume         allow continuing in a non-empty output directory (stages with built-in resume: p2, p3-fetch, p5-sss, ...)
  --force-deps     do not check upstream completion markers (e.g. when [upstream] points to existing product directories)
Each stage: output directory <out_root>/<stage id>; child environment REPRO_CONFIG = absolute config path,
  REPRO_OUTPUT_DIR = output directory, PYTHONUNBUFFERED=1;
  cwd = pipeline/. stdout/stderr are also written to <out>/_repro_stdout.txt; <out>/_repro_run.json records the command,
  config SHA-256, script SHA-256, SHA-256 of staged inputs, options, start/end and exit code. The host name is recorded only if
  the environment variable REPRO_RECORD_HOST=1 is set. Nothing is deleted; a non-empty output directory is refused without --resume.
Anemometer-height switch: [options] wind_height = uniform_4m (default, as in the released data) / per_station / per_source
  (per-source heights, implemented only by stage dw13a); see docs/reproduction.md. Use a separate out_root for non-default settings.
Multi-step modes ({"steps": [...]}) run in order and stop at the first non-zero exit.
Exit code: the child exit code (first non-zero with --all); 2 = configuration or argument error.
Change Log
  2026-09-27 first version; 2026-09-29 options and wind switch in run records; 2026-09-29c multi-step modes (dw13a).
  2026-09-30 host name recorded only on request.
  2026-10-01 environment variables REPRO_*; the run record no longer holds an assembly manifest SHA-256.
"""

import argparse
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.join(HERE, "pipeline")
sys.path.insert(0, HERE)
sys.path.insert(0, PIPE)

import stages as S  # noqa: E402
import repro_paths as rp  # noqa: E402

OPT_DEFAULTS = {"raw_rate_mbps": 2.0}


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def out_dir(sid):
    return os.path.join(rp.path("out_root"), sid)


def expand(tok, sid, py):
    """Expand all {...} placeholders in one argv fragment."""
    out, i = [], 0
    while i < len(tok):
        j = tok.find("{", i)
        if j < 0:
            out.append(tok[i:])
            break
        k = tok.find("}", j)
        if k < 0:
            raise ValueError(f"unclosed placeholder: {tok}")
        out.append(tok[i:j])
        key = tok[j + 1:k]
        if key == "out":
            out.append(out_dir(sid))
        elif key == "py":
            out.append(py)
        elif key == "tag":
            out.append(rp.run_tag(sid))
        elif key.startswith("path:"):
            out.append(rp.path(key[5:]))
        elif key.startswith("up:"):
            out.append(rp.upstream(key[3:]))
        elif key.startswith("raw:"):
            out.append(rp.raw_sub(key[4:]))
        elif key.startswith("opt:"):
            out.append(str(rp.get(key[4:], OPT_DEFAULTS.get(key[4:]))))
        else:
            raise ValueError(f"unknown placeholder {{{key}}} (stage {sid})")
        i = k + 1
    return "".join(out)


def build(sid, mode=None):
    st = S.BY_ID[sid]
    mode = mode or st["default"]
    if mode not in st["modes"]:
        raise SystemExit(f"stage {sid} has no mode {mode} (available: {sorted(st['modes'])})")
    py = os.path.expanduser(os.path.expandvars(rp.get("python") or sys.executable))
    spec = st["modes"][mode]
    steps = spec["steps"] if isinstance(spec, dict) else [spec]
    argvs = [[py, os.path.join(PIPE, st["script"])] + [expand(t, sid, py) for t in s] for s in steps]
    inputs = {rel: expand(src, sid, py) for rel, src in st.get("inputs", {}).items()}
    ensure = [expand(d, sid, py) for d in st.get("ensure_dirs", [])]
    return st, mode, argvs, inputs, ensure


def wind_conflict(st, wmode):
    """Mutual exclusions between the anemometer-height switch and stages (returns the reason for refusal or None)."""
    if wmode == "per_source" and st.get("wind") == "direct":
        return ("wind_height=per_source is implemented only by stage dw13a (heights per source file/variable); "
                "this stage has no per-source tracking and would produce 4 m results labelled as corrected -- use uniform_4m or per_station")
    if wmode == "per_station" and st.get("wind") == "source":
        return "this stage applies per-source heights itself; per_station would apply a second correction -- use uniform_4m or per_source"
    return None


def is_done(sid):
    st = S.BY_ID[sid]
    d = out_dir(sid)
    return bool(st.get("done")) and all(os.path.exists(os.path.join(d, f)) for f in st["done"])


def upstream_ok(sid):
    missing = []
    for dep in S.BY_ID[sid]["deps"]:
        if not is_done(dep):
            missing.append(dep)
    return missing


def run_stage(sid, mode, cfg_path, dry=False, resume=False, force_deps=False):
    st, mode, argvs, inputs, ensure = build(sid, mode)
    od = out_dir(sid)
    env = dict(os.environ, REPRO_CONFIG=cfg_path, REPRO_OUTPUT_DIR=od, PYTHONUNBUFFERED="1")
    print(f"[run] {sid}（{mode}）：{st['title']}", flush=True)
    wmode = rp.wind_height_mode()
    if wmode != "uniform_4m" and st.get("wind"):
        print(f"[run] note: wind_height={wmode} (not the default 4 m setting); wind effect class of this stage: {st['wind']}"
              f"{' -- constants hard-coded from upstream results do not follow (docs/reproduction.md, known issues)' if st['wind'] == 'const' else ''}", flush=True)
    for argv in argvs:
        print("      " + " ".join(shlex.quote(a) for a in argv), flush=True)
    why = wind_conflict(st, wmode)
    if why:
        print(f"[run] refused: {why}", file=sys.stderr)
        return 2
    if dry:
        for rel, src in inputs.items():
            print(f"      input {src} -> {st.get('stage_dir', 'out')}/{rel}")
        return 0
    if not force_deps:
        miss = upstream_ok(sid)
        if miss:
            print(f"[run] upstream not complete: {miss} (or use --force-deps and point [upstream] to existing products)", file=sys.stderr)
            return 2
    if st["token"]:
        tf = rp.path("token_file")
        if not os.path.isfile(tf):
            print(f"[run] stage {sid} needs an Earthdata token file; the configured path does not exist (existence checked only, not read)", file=sys.stderr)
            return 2
    if os.path.isdir(od) and os.listdir(od) and not resume:
        print(f"[run] output directory not empty: {od} (add --resume to continue, or use a new out_root; nothing is deleted)", file=sys.stderr)
        return 2
    os.makedirs(od, exist_ok=True)
    for d in ensure:
        os.makedirs(d, exist_ok=True)
    staged = {}
    if inputs:
        base = od if st.get("stage_dir") == "out" else os.path.join(od, st.get("stage_dir", "_inputs"))
        for rel, src in inputs.items():
            dst = os.path.join(base, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)
            staged[rel] = {"src": src, "sha256": sha256(dst)}
    rec = {"stage": sid, "mode": mode, "argv": argvs[0] if len(argvs) == 1 else argvs, "cwd": PIPE, "config": cfg_path, "config_sha256": sha256(cfg_path),
           "script_sha256": sha256(os.path.join(PIPE, st["script"])), "inputs_staged": staged,
           "scheduler_attempt": os.environ.get("REPRO_ATTEMPT_ID"),
           "host": platform.node() if os.environ.get("REPRO_RECORD_HOST") == "1" else None, "start": now(),
           "options": rp.load().get("options", {}), "wind_height": wmode, "wind_affected": st.get("wind"),
           "wind_heights_used": rp.wind_heights() if wmode == "per_station" else None}
    t0 = time.monotonic()
    rc = 0
    with open(os.path.join(od, "_repro_stdout.txt"), "ab") as logf:
        for argv in argvs:
            p = subprocess.Popen(argv, cwd=PIPE, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            for line in p.stdout:
                logf.write(line)
                logf.flush()
                sys.stdout.buffer.write(line)
                sys.stdout.flush()
            rc = p.wait()
            if rc != 0:
                break
    rec.update({"end": now(), "seconds": round(time.monotonic() - t0, 1), "rc": rc,
                "done_files": {f: os.path.exists(os.path.join(od, f)) for f in st.get("done", [])}})
    runs = os.path.join(od, "_repro_run.json")
    hist = []
    if os.path.exists(runs):
        with open(runs, encoding="utf-8") as f:
            hist = json.load(f)
    hist.append(rec)
    with open(runs, "w", encoding="utf-8") as f:
        json.dump(hist, f, ensure_ascii=False, indent=1)
    print(f"[run] {sid} finished rc={rc}, {rec['seconds']} s", flush=True)
    return rc


def main():
    ap = argparse.ArgumentParser(description="Pipeline entry point (rain-freshening test of RIM-3 and global recalculation)")
    ap.add_argument("--config", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--stage")
    g.add_argument("--all", action="store_true")
    g.add_argument("--list", action="store_true")
    ap.add_argument("--mode")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--force-deps", action="store_true")
    a = ap.parse_args()
    cfg = os.path.abspath(a.config)
    if not os.path.isfile(cfg):
        print(f"config not found: {cfg}", file=sys.stderr)
        return 2
    os.environ[rp.ENV] = cfg
    rp.load(cfg, force=True)
    for k in ("raw_root", "out_root", "fast_root"):
        if rp.path(k).startswith(rp.PLACEHOLDER):
            print(f"config lacks [paths] {k}", file=sys.stderr)
            return 2
    if a.list:
        print(f"# wind_height={rp.wind_height_mode()}; flags: N = network, T = Earthdata token, A = part of --all; "
              "W/w/c/s = affected by the anemometer-height switch directly / via upstream / via hard-coded constants / own per-source heights (dw13a)")
        for sid in S.topo():
            st = S.BY_ID[sid]
            flags = ("N" if st["net"] else "-") + ("T" if st["token"] else "-") + ("A" if st["in_all"] else "-") \
                + {"direct": "W", "inherit": "w", "const": "c", "source": "s"}.get(st.get("wind"), "-")
            print(f"{sid:15s} {flags} {'done' if is_done(sid) else '    '}  deps {','.join(st['deps']) or '-':34s} {st['title']}")
        return 0
    if a.stage:
        if a.stage not in S.BY_ID:
            print(f"unknown stage {a.stage}", file=sys.stderr)
            return 2
        return run_stage(a.stage, a.mode, cfg, a.dry_run, a.resume, a.force_deps)
    for sid in S.topo([s["id"] for s in S.STAGES if s["in_all"]]):
        if not a.dry_run and is_done(sid):
            print(f"[run] {sid} already complete, skipped")
            continue
        rc = run_stage(sid, None, cfg, a.dry_run, a.resume, a.force_deps)
        if rc != 0:
            print(f"[run] --all stopped at {sid} (rc={rc})", file=sys.stderr)
            return rc
    return 0


if __name__ == "__main__":
    sys.exit(main())
