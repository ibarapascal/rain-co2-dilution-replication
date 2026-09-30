"""repro_io.py — 阶段脚本共用的原件取数小工具（只做 I/O；一切下载经 raw_store 落 bulk 原始档，不删原件）。

  cmorph_fetcher(exc_cls, log, owner)  → raw_store.CmorphRaw：接口同 p2_rim_test.Fetcher.fetch(h)（返回 ("ok", 原件路径) 或
                                         ("missing", 原因)），原件已在 raw 档就直接用，否则按 [options] raw_rate_mbps 限速下载进
                                         raw 档并记 raw/manifest.jsonl；调用方不得删除返回的路径。P2 由此取代「下载到临时目录→抽像元→删」。
  raw_file(url, rel, dataset, exc_cls) → 单个原件（如 GLODAP tar）落 raw 档后返回本地路径；P2 的 GLODAP 由「流式解包不落原件」改为
                                         先落原件再从本地文件流式解包（解出的成员仍是派生物，照原逻辑用完即删）。
  refuse_overwrite(dst)               → raw 档只增不改：目标已存在就抛 FileExistsError（raw_archive.Copier.copy 用；
                                         原为「同名不同大小即 os.replace 覆盖」）。
  dup_if_exists(root, rel)            → 已有同名原件时返回另存的 (dst, rel)：rel＋.dup<unix 秒>[_n]（与 p7a keep() 同一约定）；
                                         p5_sss_sat 的 OPeNDAP 子集落档 keep() 用（原为无条件覆盖）。

规则：原始数据一律保留。
Change Log：
  2026-09-27 初版。
  2026-10-01 refuse_overwrite、dup_if_exists（raw 档不覆盖）。
"""

import os
import time

import raw_store as rs
import repro_paths as _rp


def _rate():
    return float(_rp.get("raw_rate_mbps", 2.0))


def cmorph_fetcher(exc_cls, log=None, owner="p2"):
    return rs.CmorphRaw(_rp.path("raw_root"), _rate(), owner=owner, source_job=_rp.run_tag(owner),
                        exc_cls=exc_cls, log=log, heartbeat=False)


def raw_file(url, rel, dataset, exc_cls, owner="p2"):
    root = rs.ensure_root(_rp.path("raw_root"))
    dl = rs.Downloader(root, rs.TokenBucket(_rate() * rs.MB), rs.Manifest(root), owner, _rp.run_tag(owner))
    try:
        kind, val, _rec = dl.get(url, rel, dataset)
    except (rs.DownloadError, rs.AuthFail) as e:
        raise exc_cls(f"{url} → {e}")
    if kind == "missing":
        raise exc_cls(f"{url} → {val}")
    return val


def refuse_overwrite(dst):
    if os.path.exists(dst):
        raise FileExistsError(f"raw archive is append-only: {dst} exists with a different size; not replaced "
                              "(check which version is right, then move the old file aside or use a new raw_root)")


def dup_if_exists(root, rel):
    dst = os.path.join(root, rel)
    if not os.path.exists(dst):
        return dst, rel
    base, n = rel + f".dup{int(time.time())}", 0
    rel2 = base
    while os.path.exists(os.path.join(root, rel2)):
        n += 1
        rel2 = f"{base}_{n}"
    return os.path.join(root, rel2), rel2
