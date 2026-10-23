"""Second source of georeferenced maps: MapWarper.net (open API, volunteer GCPs).

Phase 1 lists every map (about 83,000; 50 per page) into data/raw/mapwarper_maps.jsonl.
Phase 2 fetches the GCPs of warped maps that have a usable date (date_depicted, else a
year in the title) between 1450 and 1950 and a bounding box at least 2 degrees across,
into data/raw/mapwarper_gcps.jsonl. One request at a time; both phases resume.
"""
import json
import socket
import sys
import time
import urllib.request

from common import RAW, USER_AGENT
from dates import parse_year

LIST = RAW / "mapwarper_maps.jsonl"
GCPS = RAW / "mapwarper_gcps.jsonl"
API = "https://mapwarper.net/api/v1"
PAGE = 50
socket.setdefaulttimeout(60)


def get(url):
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode("utf-8", errors="replace"))
        except Exception as e:
            err = e
            time.sleep(5 * (attempt + 1))
    raise err


def dated(a):
    for src, text in (("date_depicted", a.get("date_depicted")), ("title", a.get("title"))):
        y = parse_year(text)
        if y and 1450 <= y[0] <= 1950:
            return y, src
    return None, None


def extent(bbox):
    try:
        x0, y0, x1, y1 = map(float, bbox.split(","))
    except (AttributeError, ValueError):
        return 0.0
    return max(x1 - x0, y1 - y0)


def phase1():
    done = set()
    if LIST.exists():
        for line in open(LIST, encoding="utf-8"):
            done.add(json.loads(line)["page"])
    first = get(f"{API}/maps?page=1&per_page={PAGE}")  # JSON:API page[number] is ignored by this server
    pages = first["meta"]["total_pages"]
    print(f"phase 1: {pages} pages, {len(done)} done", flush=True)
    with open(LIST, "a", encoding="utf-8") as fh:
        for p in range(1, pages + 1):
            if p in done:
                continue
            j = first if p == 1 else get(f"{API}/maps?page={p}&per_page={PAGE}")
            for d in j["data"]:
                a = d["attributes"]
                fh.write(json.dumps({"page": p, "id": d["id"], **{k: a.get(k) for k in (
                    "title", "date_depicted", "bbox", "status", "map_type", "width", "height",
                    "source_uri", "created_at")}, "description": (a.get("description") or "")[:300]},
                    ensure_ascii=False) + "\n")
            if not j["data"]:
                fh.write(json.dumps({"page": p, "id": None}) + "\n")
            fh.flush()
            if p % 100 == 0:
                print("page", p, flush=True)
            time.sleep(0.3)


def phase2():
    maps = {}
    for line in open(LIST, encoding="utf-8"):
        r = json.loads(line)
        if r.get("id") and r.get("status") == "warped" and r.get("map_type") == "is_map":
            y, _ = dated(r)
            if y and extent(r.get("bbox")) >= 2:
                maps[r["id"]] = r
    done = set()
    if GCPS.exists():
        for line in open(GCPS, encoding="utf-8"):
            done.add(json.loads(line)["id"])
    todo = [i for i in maps if i not in done]
    print(f"phase 2: {len(maps)} dated regional maps, {len(todo)} to fetch", flush=True)
    with open(GCPS, "a", encoding="utf-8") as fh:
        for k, i in enumerate(todo):
            try:
                j = get(f"{API}/maps/{i}/gcps")
                pts = [[float(g["attributes"]["x"]), float(g["attributes"]["y"]),
                        float(g["attributes"]["lon"]), float(g["attributes"]["lat"])] for g in j["data"]]
                rec = {"id": i, "gcps": pts}
            except Exception as e:
                rec = {"id": i, "error": str(e)[:200]}
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if k % 200 == 0:
                print("gcps", k, flush=True)
            time.sleep(0.3)


if __name__ == "__main__":
    phase1()
    phase2()
    sys.exit(0)
