#!/usr/bin/env python3
"""build_event_table.py — merge the per-stage event tables of the mooring test into one 646-row event-level table.

Inputs (released copies in data/products/, all keyed by station + onset_utc):
  p1/p1_events.csv, p1b/p1b_events_new.csv   event definition (round 1 = P1, round 2 = P1b expansion)      -> def_*
  p2/p2_events.csv                            observed and modeled freshening, pCO2 (main test)             -> p2_*
  p4-mech/p4_mech_events.csv                  first-round mechanism diagnostics                              -> p4_*
  p4b/p4b_events.csv                          second-round mechanism diagnostics                             -> p4b_*
  p4e/p4e_events.csv                          freshwater budget, group A only, one row per window            -> p4e_<window>_*
  p5-beta/p5_beta_events.csv                  chemical-slope candidates                                      -> p5b_*
  p7c/p7c_events.csv                          wind-dependence candidates                                     -> p7c_*
  dw13a/measured/dw13a_events_full.csv        anemometer-height sensitivity (measured heights)               -> dw13a_*
Output: data/events/events_merged.csv (one row per event; empty cell = not available for that event or stage).
The 204 P1 events (group A) and the 442 P1b events (group B) are disjoint and together form the 646 events of p2_events.csv.
Every input must match this event set exactly; the script stops if a key is duplicated, unmatched or extra. Standard library only; runs in about a second.

Usage: python3 code/tools/build_event_table.py [--data DATA_DIR] [--out OUT_CSV]
Change Log: 2026-09-30 first version.
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "..", "data"))
KEY = ("station", "onset_utc")
DROP = {"station", "onset_utc", "group", "season"}


def read(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def index(rows, name):
    out = {}
    for r in rows:
        k = (r["station"], r["onset_utc"])
        if k in out:
            sys.exit(f"[build_event_table] duplicate key {k} in {name}")
        out[k] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=DATA)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    prod = os.path.join(a.data, "products")
    out_path = a.out or os.path.join(a.data, "events", "events_merged.csv")

    base = read(os.path.join(prod, "p2", "p2_events.csv"))
    keys = [(r["station"], r["onset_utc"]) for r in base]
    if len(set(keys)) != len(keys):
        sys.exit("[build_event_table] duplicate keys in p2_events.csv")

    # event definitions: P1 (round 1, group A) and P1b (round 2, group B); the two sets are disjoint
    p1 = index(read(os.path.join(prod, "p1", "p1_events.csv")), "p1_events")
    p1b = index(read(os.path.join(prod, "p1b", "p1b_events_new.csv")), "p1b_events_new")
    defs = {}
    for k, r in p1.items():
        defs[k] = dict(r, round="1")
    for k, r in p1b.items():
        defs[k] = dict(r, round="2")

    tables = [
        ("p4", index(read(os.path.join(prod, "p4-mech", "p4_mech_events.csv")), "p4_mech_events")),
        ("p4b", index(read(os.path.join(prod, "p4b", "p4b_events.csv")), "p4b_events")),
        ("p5b", index(read(os.path.join(prod, "p5-beta", "p5_beta_events.csv")), "p5_beta_events")),
        ("p7c", index(read(os.path.join(prod, "p7c", "p7c_events.csv")), "p7c_events")),
        ("dw13a", index(read(os.path.join(prod, "dw13a", "measured", "dw13a_events_full.csv")), "dw13a_events_full")),
    ]
    p4e = {}
    for r in read(os.path.join(prod, "p4e", "p4e_events.csv")):
        k = (r["station"], r["onset_utc"])
        w = r["window"].replace("-", "_")
        if (k, w) in p4e:
            sys.exit(f"[build_event_table] duplicate p4e key {k} {w}")
        p4e[(k, w)] = r
    windows = sorted({w for (_k, w) in p4e})

    def_cols = [c for c in next(iter(p1b.values())).keys() if c not in DROP] + ["round"]
    def_cols = ["round"] + [c for c in def_cols if c != "round"]
    base_cols = [c for c in base[0].keys() if c not in DROP]
    cols = ["event_id", "station", "group", "onset_utc", "season"] + [f"def_{c}" for c in def_cols] + [f"p2_{c}" for c in base_cols]
    for pre, t in tables:
        cols += [f"{pre}_{c}" for c in next(iter(t.values())).keys() if c not in DROP]
    p4e_cols = [c for c in next(iter(p4e.values())).keys() if c not in DROP | {"window"}]
    for w in windows:
        cols += [f"p4e_{w}_{c}" for c in p4e_cols]

    unmatched = {pre: 0 for pre, _ in tables}
    missing_def = 0
    rows = []
    for i, r in enumerate(base):
        k = (r["station"], r["onset_utc"])
        o = {"event_id": f"E{i + 1:03d}", "station": r["station"], "group": r["group"], "onset_utc": r["onset_utc"],
             "season": r["season"]}
        d = defs.get(k)
        if d is None:
            missing_def += 1
        for c in def_cols:
            o[f"def_{c}"] = d.get(c, "") if d else ""
        for c in base_cols:
            o[f"p2_{c}"] = r[c]
        for pre, t in tables:
            m = t.get(k)
            if m is None:
                unmatched[pre] += 1
                continue
            if m.get("group") and m["group"] != r["group"]:
                sys.exit(f"[build_event_table] group mismatch {k} in {pre}")
            for c, v in m.items():
                if c not in DROP:
                    o[f"{pre}_{c}"] = v
        for w in windows:
            m = p4e.get((k, w))
            if m:
                for c in p4e_cols:
                    o[f"p4e_{w}_{c}"] = m[c]
        rows.append(o)
    extra = {pre: len(set(t) - set(keys)) for pre, t in tables}
    if missing_def or any(unmatched.values()) or any(extra.values()):
        sys.exit(f"[build_event_table] join incomplete: missing definitions {missing_def}, unmatched {unmatched}, extra {extra}")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"[build_event_table] {len(rows)} events x {len(cols)} columns -> {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
