"""Python 环境自测：打印解释器与 P2/P3 依赖版本（P3 起含 gsw、xarray），可选对数据 binding 做一次写读删，
写一份 JSON 到 REPRO_OUTPUT_DIR。用法：venv_selftest.py [--rw-dir 绝对路径]（P3：fast_root 下 p3-cache）。"""
import json, os, sys, time
import numpy, pandas, scipy, netCDF4, h5py, PyCO2SYS, matplotlib, gsw, xarray
mods = (numpy, pandas, scipy, netCDF4, h5py, PyCO2SYS, matplotlib, gsw, xarray)
info = {"version": sys.version.split()[0], "executable": sys.executable, "prefix": sys.prefix,
        "base_prefix": sys.base_prefix,
        "pkgs": {m.__name__: getattr(m, "__version__", "?") for m in mods},
        # gsw 冒烟：SP=35、p=0、经纬 0 → SA≈35.165
        "gsw_SA_from_SP_35": float(gsw.SA_from_SP(35.0, 0.0, 0.0, 0.0))}
if "--rw-dir" in sys.argv:
    d = sys.argv[sys.argv.index("--rw-dir") + 1]
    p = os.path.join(d, f".selftest_{os.environ.get('REPRO_ATTEMPT_ID', 'x')}.tmp")
    payload = f"fast_root rw {time.time()}".encode()
    with open(p, "wb") as f:
        f.write(payload)
    with open(p, "rb") as f:
        ok = f.read() == payload
    os.remove(p)
    info["rw"] = {"dir": d, "write_read_ok": ok, "deleted": not os.path.exists(p)}
print(json.dumps(info, ensure_ascii=False))
with open(os.path.join(os.environ["REPRO_OUTPUT_DIR"], "selftest.json"), "w") as f:
    json.dump(info, f, indent=1)
if info.get("rw") and not (info["rw"]["write_read_ok"] and info["rw"]["deleted"]):
    sys.exit(1)
