#!/usr/bin/env python3
"""p7e_fetch.py — P7e 原始数据下载（只下载，不分析）：SPURS-1 WHOI 中心浮标＋Wave Glider（PO.DAAC）、
PISTON 2018 R/V Thompson／2019 R/V Sally Ride 船载近表层与气象相关文件（NASA ASDC），落 workstation bulk 原始档。

Description：
  原始数据一律保留；先下载好原始数据，再分析。
  下载、核对、心跳接替全部复用 `raw_fetch_list.py`（一行未改）：本文件只做三件事——
    1) `--asdc-list`（laptop，公开 CMR 元数据，无 token）：按产品名正则从 ASDC 两个 PISTON 船载 collection 生成清单
       `raw_archive_lists/piston.json`（格式同 `raw_fetch_list.cmr_list`；ASDC 的 CMR 给 Size＋SizeUnit 而非 SizeInBytes，
       换算为字节并要求是整数，否则拒绝）；SEA-POL 雷达网格（体量大）只列在清单的 `listed_not_downloaded`，不下。
       SPURS-1 清单直接用 `raw_fetch_list.py --cmr-list ... --list-name spurs1` 生成（PO.DAAC，已有路径）。
    2) 认证主机扩展：ASDC 受保护文件在 `data.asdc.earthdata.nasa.gov`（Earthdata Cloud TEA），不在 p5_sss_sat.AUTH_HOSTS 里；
       本进程内把 `raw_archive.auth_bits()` 返回的 auth_host 包一层，只额外放行这一个 Earthdata 主机（S3 预签名、登录页仍不带 token）。
       ASDC 若要求额外授权（被重定向到 urs 登录页或 401/403）→ raw_store.AuthFail → 本清单退出码 5，停下报告，不注册、不申请。
    3) `--run`：按顺序（spurs1 → piston）对每份清单调用 `raw_fetch_list.run`：等 `.p5_fetch_heartbeat` 空出 → 写心跳接替 2 MB/s 名额
       （与其他下载合计 ≤4 MB/s）→ 下载、逐文件核字节与 CMR MD5 → 写 done。两份清单各自一个输出子目录。
  认证：token 只由 raw_fetch_list.run 从 workstation 既有 token 文件读进内存，只进 Earthdata 主机请求头；本文件不读 token。

Usage：
  laptop：python p7e_fetch.py --asdc-list --list-out raw_archive_lists/piston.json
       python raw_fetch_list.py --cmr-list SPURS1_MOORING_WHOI=C2491772311-POCLOUD,SPURS1_WAVEGLIDER=C2491772321-POCLOUD \
              --list-name spurs1 --list-out raw_archive_lists/spurs1.json
  workstation：python p7e_fetch.py --run --lists raw_archive_lists/spurs1.json,raw_archive_lists/piston.json \
              --raw-root <raw_root> --rate 2 --tag raw-p7e:$REPRO_ATTEMPT_ID --out <dir>
  离线自测：python p7e_fetch.py --selftest（认证主机扩展、两份清单合法、raw_fetch_list 自测）
Dependencies：标准库；同目录 raw_fetch_list.py、raw_store.py、raw_archive.py、p5_sss_sat.py（后者仅 --run 时经 raw_archive 取认证函数）。
退出码：各清单 raw_fetch_list 退出码的最大值（0 全部完成且核对一致；2 有失败；3 异常；4 MD5 不符；5 认证失败；6 等待超时）。

Change Log：
  2026-09-27 初版（P7e）。
"""

import argparse
import json
import os
import re
import sys
import types
import urllib.parse
import urllib.request

import raw_fetch_list as rfl
import raw_store as rs

VERSION = "p7e-fetch-2026-09-27a"
CMR_G = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"
EXTRA_AUTH_HOSTS = ("data.asdc.earthdata.nasa.gov",)
ASDC_PREFIX = "https://data.asdc.earthdata.nasa.gov/asdc-prod-protected/"
UNIT = {"B": 1, "KB": 2 ** 10, "MB": 2 ** 20, "GB": 2 ** 30}

PISTON_COLLS = {  # 短名 → CMR concept id（LARC_CLOUD；collection DOI 见清单）
    "PISTON-ONR-NOAA_RVThompson_2018": "C3880791367-LARC_CLOUD",
    "PISTON-ONR-NOAA_RVSallyRide_2019": "C3880791409-LARC_CLOUD",
}
# 下载：SurfOtter（近表层拖曳剖面）、nav-met-sea（船位／光学雨量计／风／TSG，1 min、10 min、60 min 通量）、
# 雨滴谱仪、ROSR 皮温、uCTD、CTD、Chameleon（微结构剖面）
KEEP = re.compile(r"(SurfOtter|nav-met-sea|disdrometer|rosr|uctd|piston-CTD|Chameleon)", re.I)
LIST_ONLY = re.compile(r"SEAPOL", re.I)   # SEA-POL C 波段雷达：只列不下


def asdc_size_bytes(a):
    s, u = a.get("Size"), (a.get("SizeUnit") or "").upper()
    if s is None or u not in UNIT:
        return None
    b = s * UNIT[u]
    return int(round(b)) if abs(b - round(b)) < 1e-3 else None


def asdc_list(name="piston"):
    files, colls, listed = [], {}, []
    for short, cid in PISTON_COLLS.items():
        q = urllib.parse.urlencode({"collection_concept_id": cid, "page_size": 2000})
        req = urllib.request.Request(f"{CMR_G}?{q}", headers={"User-Agent": rs.UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            hits = int(r.headers.get("CMR-Hits", "-1"))
            d = json.load(r)
        items = d.get("items", [])
        if hits != len(items):
            raise RuntimeError(f"{cid}：CMR-Hits {hits} ≠ 返回 {len(items)}（需翻页）")
        colls[short] = {"concept_id": cid, "granules_in_collection": len(items), "granules": 0, "files": 0,
                        "bytes": 0}
        for it in items:
            u = it["umm"]
            g = u["GranuleUR"]
            keep, lonly = bool(KEEP.search(g)), bool(LIST_ONLY.search(g))
            if not (keep or lonly):
                continue
            url = next((x["URL"] for x in u.get("RelatedUrls", []) if x.get("Type") == "GET DATA"
                        and x["URL"].startswith(ASDC_PREFIX)), None)
            if url is None:
                raise RuntimeError(f"{g}：没有 ASDC https 下载链接")
            fn = url.rsplit("/", 1)[1]
            a = next((x for x in u["DataGranule"]["ArchiveAndDistributionInformation"]
                      if x.get("Name") == fn and "Size" in x), None)
            if a is None:
                raise RuntimeError(f"{g}：ArchiveAndDistributionInformation 无大小")
            b = asdc_size_bytes(a)
            ck = a.get("Checksum") or {}
            md5 = ck.get("Value", "").lower() if ck.get("Algorithm", "").upper() == "MD5" else None
            t = (u.get("TemporalExtent") or {}).get("RangeDateTime") or {}
            rec = {"short_name": short, "concept_id": cid, "granule": g, "kind": "data", "url": url,
                   "rel": f"{name}/{short}/{fn}", "bytes": b, "md5": md5,
                   "production": u["DataGranule"].get("ProductionDateTime"),
                   "time": [t.get("BeginningDateTime"), t.get("EndingDateTime")]}
            if lonly:
                rec["bytes_approx"] = b if b is not None else a.get("Size") * UNIT.get((a.get("SizeUnit") or "").upper(), 1)
                listed.append(rec)
                continue
            if b is None or not md5:
                raise RuntimeError(f"{g}：大小不能换算为整数字节或缺 MD5（{a}）")
            files.append(rec)
            colls[short]["granules"] += 1
            colls[short]["files"] += 1
            colls[short]["bytes"] += b
    files.sort(key=lambda f: (f["short_name"], f["granule"]))
    return {"name": name, "built": rs.now_iso(),
            "source": "CMR granules.umm_json（公开元数据，无 token）；ASDC Size×SizeUnit(2^n) 换算整数字节",
            "rel_prefix": name, "collections": colls, "n_files": len(files),
            "total_bytes": sum(f["bytes"] for f in files), "files": files,
            "keep_regex": KEEP.pattern, "listed_not_downloaded": {
                "regex": LIST_ONLY.pattern, "reason": "SEA-POL C 波段雷达网格与快视图体量大（约 10 GB），先只列清单不下",
                "n_files": len(listed), "bytes_approx": int(sum(x["bytes_approx"] for x in listed)),
                "files": listed}}


def patch_auth():
    """把 raw_archive.auth_bits 返回的 auth_host 扩展到 ASDC 主机（只在本进程内）。"""
    import raw_archive as ra
    orig = ra.auth_bits

    def auth_bits():
        read_token, auth_host, redact, AuthError = orig()

        def ah(host):
            h = (host or "").lower()
            return auth_host(h) or h in EXTRA_AUTH_HOSTS
        return read_token, ah, redact, AuthError
    ra.auth_bits = auth_bits
    return ra


def run(args):
    patch_auth()
    rcs = {}
    for lp in [x.strip() for x in args.lists.split(",") if x.strip()]:
        L = rfl.load_list(lp)
        sub = os.path.join(args.out, L["name"])
        os.makedirs(sub, exist_ok=True)
        log = rfl._Log(os.path.join(sub, "raw_fetch_log.txt"))
        log.log(f"=== start {VERSION} (raw_fetch_list {rfl.VERSION}) list={lp} raw={args.raw_root} rate={args.rate} MB/s")
        ns = types.SimpleNamespace(list=lp, raw_root=args.raw_root, rate=args.rate, conc=args.conc,
                                   wait_max_s=args.wait_max_s, settle_s=args.settle_s, tag=args.tag)
        try:
            rc = rfl.run(ns, sub, log)
        except Exception:
            import traceback
            log.log("FATAL 未预期异常：\n" + traceback.format_exc())
            rc = 3
        log.log(f"=== done rc={rc}")
        rcs[L["name"]] = rc
        print(f"[p7e_fetch] {L['name']} rc={rc}", flush=True)
    with open(os.path.join(args.out, "p7e_fetch_rc.json"), "w") as f:
        json.dump({"version": VERSION, "rc": rcs, "t": rs.now_iso()}, f, ensure_ascii=False)
    return max(rcs.values()) if rcs else 3


def selftest():
    res = []

    def chk(n, ok, d=""):
        res.append(bool(ok))
        print(f"  {'ok ' if ok else 'BAD'}  {n}  {d}", flush=True)

    import raw_archive as ra
    fake = types.SimpleNamespace(read_token=None, auth_host=lambda h: h in ("archive.podaac.earthdata.nasa.gov",)
                                 or h.endswith(".earthdatacloud.nasa.gov"), redact=str, AuthError=Exception)
    ra_orig = ra.auth_bits
    ra.auth_bits = lambda: (fake.read_token, fake.auth_host, fake.redact, fake.AuthError)
    try:
        patch_auth()
        _, ah, _, _ = ra.auth_bits()
        chk("认证主机：ASDC TEA 放行", ah("data.asdc.earthdata.nasa.gov"))
        chk("认证主机：PO.DAAC 仍放行", ah("archive.podaac.earthdata.nasa.gov"))
        chk("认证主机：S3 预签名不带 token", not ah("asdc-prod-protected.s3.us-west-2.amazonaws.com"))
        chk("认证主机：登录页不带 token", not ah("urs.earthdata.nasa.gov"))
        chk("认证主机：仿冒后缀不放行", not ah("data.asdc.earthdata.nasa.gov.evil.com"))
    finally:
        ra.auth_bits = ra_orig
    here = os.path.dirname(os.path.abspath(__file__))
    for nm in ("spurs1", "piston"):
        lp = os.path.join(here, "raw_archive_lists", f"{nm}.json")
        if not os.path.exists(lp):
            chk(f"清单 {nm}.json 存在", False)
            continue
        L = rfl.load_list(lp)
        chk(f"清单 {nm}.json 合法", True, f"{L['n_files']} 文件 {L['total_bytes']} B")
        chk(f"清单 {nm}：每个文件有 CMR MD5 与整数字节", all(x["md5"] and re.fullmatch(r"[0-9a-f]{32}", x["md5"])
                                                    and isinstance(x["bytes"], int) for x in L["files"]))
        chk(f"清单 {nm}：rel 都在 {nm}/ 下", all(x["rel"].startswith(nm + "/") for x in L["files"]))
        if nm == "piston":
            chk("清单 piston：URL 都是 ASDC https", all(x["url"].startswith(ASDC_PREFIX) for x in L["files"]))
            chk("清单 piston：SEA-POL 不在下载清单", not any(LIST_ONLY.search(x["granule"]) for x in L["files"]))
    n = len(res)
    print(f"p7e_fetch selftest {sum(res)}/{n}", flush=True)
    rc = rfl.selftest(rfl._Log())
    return 0 if all(res) and rc == 0 else 4


def main(argv=None):
    ap = argparse.ArgumentParser(description="P7e 原始数据下载（SPURS-1＋PISTON）")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--asdc-list", action="store_true")
    g.add_argument("--run", action="store_true")
    g.add_argument("--selftest", action="store_true")
    ap.add_argument("--list-out")
    ap.add_argument("--lists")
    ap.add_argument("--raw-root", default=rs.RAW_ROOT_DEFAULT)
    ap.add_argument("--out")
    ap.add_argument("--rate", type=float, default=2.0)
    ap.add_argument("--conc", type=int, default=2)
    ap.add_argument("--wait-max-s", type=float, default=6 * 3600)
    ap.add_argument("--settle-s", type=float, default=20.0)
    ap.add_argument("--tag", default="raw-p7e")
    a = ap.parse_args(argv)
    if a.asdc_list:
        L = asdc_list()
        with open(a.list_out, "w", encoding="utf-8") as f:
            json.dump(L, f, ensure_ascii=False, indent=1)
            f.write("\n")
        print(json.dumps({"collections": L["collections"], "n_files": L["n_files"], "total_bytes": L["total_bytes"],
                          "listed_not_downloaded": {k: v for k, v in L["listed_not_downloaded"].items()
                                                    if k != "files"}}, ensure_ascii=False, indent=1))
        return 0
    if a.selftest:
        return selftest()
    if a.rate > 2.0 or a.rate <= 0:
        print("限速越界：本任务只接替 2 MB/s 名额", file=sys.stderr)
        return 3
    a.out = a.out or os.environ.get("REPRO_OUTPUT_DIR")
    if not a.out or not a.lists:
        print("需要 --lists 与 --out", file=sys.stderr)
        return 3
    os.makedirs(a.out, exist_ok=True)
    return run(a)


if __name__ == "__main__":
    sys.exit(main())
