#!/usr/bin/env python3
"""元数据核查（阶段 cal-dw13-meta；不计算，只读公开元数据）：
③ MAPCO2 盐度传感器实际深度与来源——NCEI OCADS 各站各部署 PI_OME／元数据 XML 的 Sea_Surface_Salinity 段（Location、Model、
   Other_Comments）与 Depth_of_Sea_Water_Intake；按 P4 第一轮 646 事件的 onset 映射到部署，统计事件所在部署的传感器来源类别。
④ B 群（KEO、Papa、MOSEAN/WHOTS、SOFS、Stratus）风速计高度——P1b 站单里提供风的 OceanSITES 文件，经 NDBC THREDDS OPeNDAP
   只取高度坐标变量（.dds 定名、.ascii 取值）与 sensor_height 属性。
网络请求只带通用 User-Agent，不带任何身份信息；下载的 XML／DAS 原样存 <out>/raw_meta/（元数据原件，原样保留）。

用法：python3 cal_dw13_meta.py --p1b-stations p1b_stations.json --p4e p4_mech_events.csv --out DIR
输出：dw13_meta.json（逐文件解析表、事件映射计数、风高表）、raw_meta/（原件＋SHA256SUMS）。
Change Log:
  2026-09-29：初版；同日补标量高度变量（首跑漏 WHOTS／Stratus 的 HEIGHT_WND，只影响风高表）。
"""
import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import os
import re
import urllib.request

UA = {"User-Agent": "rain-co2-dilution-replication-metadata-check/1.0"}
NCEI = "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/"
ACCESSIONS = {"Stratus": "0100075", "TAO140W": "0100077", "TAO165E": "0113238", "TAO8S165E": "0117073", "WHOTS_a": "0100073",
              "MOSEAN/WHOTS": "0100080", "TAO125W": "0100076", "TAO170W": "0100078", "TAO155W": "0100084", "TAO110W": "0112885",
              "SOFS": "0118546", "KEO": "0100071", "BOBOA": "0162473", "Papa": "0100074"}
GDAC = "https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites/"
HNAMES = ("HEIGHTWIND", "HEIGHT_WIND", "HEIGHT_WND", "WSPD_H", "WND_measurementHeight", "HEIGHT_WIND2")
MON = {m: i + 1 for i, m in enumerate("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split())}


def get(url, path=None):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
        b = r.read()
    if path:
        with open(path, "wb") as f:
            f.write(b)
    return b.decode("utf-8", "ignore")


def flat(t):
    return re.sub(r"\s+", " ", re.sub(r"<([^/>][^>]*)>", r" [\1] ", re.sub(r"</[^>]+>", "", t)))


def parse_ome(t):
    s = flat(t)
    m = s.find("[Sea_Surface_Salinity]")
    sec = None
    if m >= 0:
        e = min([k for k in (s.find(x, m + 20) for x in ("[Atmospheric_Pr", "[Other_Sensors]", "[Method_Description]")) if k > 0] or [len(s)])
        sec = s[m:e]
    g1 = lambda pat, x: (re.search(pat, x or "") or [None, None])[1]
    dep = re.search(r"\[Depth_of_Sea_Water_Intake\] ([^\[]*)", s)
    return {"loc": (g1(r"\[Location\] ([^\[]*)", sec) or "").strip() or None, "model": (g1(r"\[Model\] ([^\[]*)", sec) or "").strip() or None,
            "comment": (g1(r"\[Other_Comments\] ([^\[]*)", sec) or "").strip()[:300] or None,
            "intake": dep.group(1).strip() if dep else None}


def category(r):
    c = r["comment"] or ""
    if (r["model"] or "").startswith("ATLAS") or c.startswith("Sea Surface Salinity collected by NOAA/NDBC/TAO"):
        return "TAO 1 m 传感器（ATLAS 或注明 NDBC/TAO 提供）"
    if any(k in c for k in ("CSIRO", "Integrated Marine", "Open Access", "another group")):
        return "外方提供（CSIRO／IMOS／另一组）"
    if r["loc"] is None:
        return "元数据未解析"
    return "MAPCO2 自带 SBE（Location 记 1 m）"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p1b-stations", required=True)
    ap.add_argument("--p4e", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    raw = os.path.join(a.out, "raw_meta")
    os.makedirs(os.path.join(raw, "ncei_ome"), exist_ok=True)
    os.makedirs(os.path.join(raw, "oceansites_dds_das"), exist_ok=True)
    res = {"ome": [], "wind": []}
    # ③
    deps = []
    for st, acc in ACCESSIONS.items():
        idx = get(NCEI + acc + "/")
        for x in sorted(set(re.findall(r'href="([^"?/]+\.xml)"', idx))):
            p = os.path.join(raw, "ncei_ome", f"{acc}__{x}")
            r = parse_ome(get(NCEI + acc + "/" + x, p))
            m = re.match(r"([A-Za-z]+)_?(\d+[EW])?_(\d+[NS])_([A-Z][a-z]{2})(\d{4})_([A-Z][a-z]{2})(\d{4})", x)
            name = m.group(1) if m else x
            if name == "TAO" and m:
                name = "TAO" + m.group(2)
            if name.startswith("TAO") and m and m.group(3) != "0N":
                name = "TAO8S165E"
            name = {"WHOTS": "MOSEAN/WHOTS"}.get(name, name)
            r.update({"file": f"{acc}/{x}", "station": name, "category": category(r)})
            if m:
                r["from"] = dt.date(int(m.group(5)), MON[m.group(4)], 1).isoformat()
                r["to"] = dt.date(int(m.group(7)), MON[m.group(6)], 28).isoformat()
                deps.append(r)
            res["ome"].append(r)
    cnt, per = collections.Counter(), collections.defaultdict(collections.Counter)
    for e in csv.DictReader(open(a.p4e, encoding="utf-8")):
        t = e["onset_utc"][:10]
        hit = [d for d in deps if d["station"] == e["station"] and d["from"] <= t <= d["to"]]
        k = hit[0]["category"] if hit else "无对应部署元数据（NCEI 未列该期 XML）"
        cnt[k] += 1
        per[e["station"]][k] += 1
    res["events_by_category"] = dict(cnt)
    res["events_by_station"] = {k: dict(v) for k, v in per.items()}
    res["ome_location_counts"] = dict(collections.Counter(r["loc"] for r in res["ome"]))
    res["ome_intake_counts"] = dict(collections.Counter(r["intake"] for r in res["ome"]))
    # ④
    j = json.load(open(a.p1b_stations, encoding="utf-8"))
    for x in j["o_group"]["candidates"]:
        if x["name"] not in ("KEO", "Papa", "MOSEAN/WHOTS", "SOFS", "Stratus"):
            continue
        for f in x["files"]:
            if not any(p[0] == "wind" for p in f["provides"]):
                continue
            base = os.path.join(raw, "oceansites_dds_das", f["path"].replace("/", "__"))
            dds = get(GDAC + f["path"] + ".dds", base + ".dds")
            das = get(GDAC + f["path"] + ".das", base + ".das")
            vals = {}
            for n in HNAMES:
                mm = re.search(r"\b" + n + r"\[(\w+) = (\d+)\]", dds)
                if mm:
                    q = n if int(mm.group(2)) <= 5 else n + "[0:1:0]"
                elif re.search(r"\b(?:Float64|Float32|Int32|Int16) " + n + r";", dds):
                    q = n                                        # 标量高度变量（WHOTS／Stratus 的 D_M 文件）
                else:
                    continue
                txt = get(GDAC + f["path"] + ".ascii?" + q)
                nums = re.findall(r"-?\d+\.?\d*(?:e-?\d+)?", txt.split("-" * 45)[-1])
                vals[n] = float(nums[-1]) if nums else None
            res["wind"].append({"station": x["name"], "file": f["path"], "start": f["start"][:10], "end": f["end"][:10],
                                "heights_m": vals, "sensor_height_attrs": sorted(set(re.findall(r"sensor_height ([\d.]+)", das)))})
    summ = {}
    for st in ("KEO", "Papa", "MOSEAN/WHOTS", "SOFS", "Stratus"):
        prim = [v for r in res["wind"] if r["station"] == st for k, v in r["heights_m"].items() if k != "HEIGHT_WIND2" and v is not None]
        summ[st] = {"primary_min": min(prim) if prim else None, "primary_max": max(prim) if prim else None, "n_files": len(prim)}
    res["wind_summary"] = summ
    with open(os.path.join(raw, "SHA256SUMS"), "w") as fo:
        for root, _, fs in os.walk(raw):
            for fn in sorted(fs):
                if fn == "SHA256SUMS":
                    continue
                p = os.path.join(root, fn)
                fo.write(hashlib.sha256(open(p, "rb").read()).hexdigest() + "  " + os.path.relpath(p, raw) + "\n")
    json.dump(res, open(os.path.join(a.out, "dw13_meta.json"), "w"), ensure_ascii=False, indent=1)
    print(json.dumps({"events_by_category": res["events_by_category"], "ome_location_counts": res["ome_location_counts"],
                      "ome_intake_counts": res["ome_intake_counts"], "wind_summary": summ}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
