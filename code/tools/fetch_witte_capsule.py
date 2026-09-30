#!/usr/bin/env python3
"""fetch_witte_capsule.py — download the two files of the Witte et al. (2026a) Code Ocean capsule used by this pipeline
(equivalence self-test of p3-global and the conservation check p4c) and verify their SHA-256.

The capsule (DOI 10.24433/CO.9378898.v1; MIT License according to its DataCite record) is not redistributed in this
repository. The files are saved to code/pipeline/witte_capsule/ (git-ignored). p3-global downloads them itself when absent;
p4c needs them in place. Standard library only.

Usage: python3 code/tools/fetch_witte_capsule.py [--dest DIR]
Exit code: 0 both files present with the expected SHA-256; 1 download failed or checksum mismatch.
Change Log: 2026-09-30 first version (endpoint and checksums as used by p3_global.py).
"""

import argparse
import hashlib
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DEST = os.path.normpath(os.path.join(HERE, "..", "pipeline", "witte_capsule"))
FILES = {"CO2_Rain_Flux_Toolbox.py": "cbce94db03fc4a3bc988146ceea71509bd3c85980e7c6b1827126eb2d863df7b",
         "main.py": "53017ba313bc094c423ece48e22e3a0ef94ac0db912c2947e1f8c69d804002b1"}
API = ("https://codeocean.com/api/capsules/6b8891ef-0ca8-43d0-ba62-edfdad8e65d5/blob"
       "?owner_id=verified&path=code/{name}&commit=HEAD&version=1")
UA = "rain-co2-dilution-repro/1.0 (research code; python)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default=DEST)
    a = ap.parse_args()
    os.makedirs(a.dest, exist_ok=True)
    ok = True
    for name, want in FILES.items():
        p = os.path.join(a.dest, name)
        if os.path.exists(p) and hashlib.sha256(open(p, "rb").read()).hexdigest() == want:
            print(f"present  {name}")
            continue
        try:
            req = urllib.request.Request(API.format(name=name), headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
        except Exception as e:  # noqa: BLE001
            print(f"FAILED   {name}: {e}  (download it manually from https://doi.org/10.24433/CO.9378898.v1, folder code/)")
            ok = False
            continue
        got = hashlib.sha256(data).hexdigest()
        if got != want:
            print(f"MISMATCH {name}: {got[:16]}... != {want[:16]}... (not saved)")
            ok = False
            continue
        with open(p, "wb") as f:
            f.write(data)
        print(f"saved    {name}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
