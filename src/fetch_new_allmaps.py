"""Catalogue metadata for maps added to Allmaps after the round-1 snapshot (PLAN2 part F).

Reads data/raw/maps_2026-10-23.ndjson, keeps maps whose annotation was created after
2026-10-13 with at least 6 GCPs and a 1-degree extent, and fetches their metadata with
the round-1 fetchers into separate files, so round 1 stays reproducible.
"""
import json
import time

from common import RAW, domain, gcp_span, iter_maps, rumsey_mid
from fetch_iiif_meta import SKIP, fetch_cdm, fetch_manifest, source_of
from fetch_rumsey_meta import fetch_one

NEW = RAW / "maps_2026-10-23.ndjson"
OUT = RAW / "r2_new_meta.jsonl"
SNAPSHOT = "2026-10-13T05:10"


def main():
    old = {m["id"] for m in iter_maps()}
    todo = []
    for m in iter_maps(NEW):
        if m["id"] in old or m.get("created", "") <= SNAPSHOT:
            continue
        if len(m["gcps"]) >= 6 and gcp_span(m) >= 1:
            todo.append(m)
    print(f"{len(todo)} new maps", flush=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        for m in todo:
            mid = rumsey_mid(m)
            rec = {"map_id": m["id"]}
            if mid:
                rec.update({"kind": "rumsey", **fetch_one(mid)})
            elif domain(m) not in SKIP and source_of(m):
                kind, key = source_of(m)
                try:
                    rec.update({"kind": kind, **(fetch_cdm(key) if kind == "cdm" else fetch_manifest(key))})
                except Exception as e:
                    rec.update({"kind": kind, "error": str(e)[:200]})
            else:
                rec["kind"] = "none"
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
            fh.flush()
            time.sleep(0.3)


if __name__ == "__main__":
    main()
