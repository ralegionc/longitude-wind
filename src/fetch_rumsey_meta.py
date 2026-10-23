"""Fetch catalogue metadata (date, place, author, title) for Rumsey maps.

Uses the public Luna JSON endpoint `as/fetchMediaSearch`, one request per map,
three at a time. Results are appended to data/raw/rumsey_meta.jsonl so the
script can be stopped and resumed.
"""
import json
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

from common import RAW, USER_AGENT, gcp_span, iter_maps, rumsey_mid

OUT = RAW / "rumsey_meta.jsonl"
API = "https://www.davidrumsey.com/luna/servlet/as/fetchMediaSearch?mid={}&fullData=true"
MIN_GCPS, MIN_SPAN = 6, 1.0
WORKERS = 3


def parse(payload):
    """Flatten Luna's doubly encoded field list into {display name: value}."""
    rec = payload[0]
    out = {"displayName": rec.get("displayName")}
    for key in ("fieldValues", "relatedItemFieldValues"):
        raw = rec.get(key)
        if not raw:
            continue
        for fv in json.loads(raw):
            field = json.loads(fv["field"])
            vals = json.loads(fv.get("valueList") or "[]") or [fv.get("value")]
            val = vals if len(vals) > 1 else vals[0]
            if isinstance(val, str) and len(val) > 400:  # long notes: keep the start only
                val = val[:400]
            out[field["displayName"]] = val
    return out


def fetch_one(mid):
    for attempt in range(4):
        try:
            req = urllib.request.Request(API.format(mid), headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=90) as r:
                payload = json.loads(r.read().decode("utf-8"))
            return {"mid": mid, **(parse(payload) if payload else {"missing": True})}
        except Exception as e:  # network hiccup or server error: back off
            err = str(e)[:200]
            time.sleep(5 * (attempt + 1))
    return {"mid": mid, "error": err}


def main():
    done = set()
    if OUT.exists():
        for line in open(OUT, encoding="utf-8"):
            rec = json.loads(line)
            if "error" not in rec:
                done.add(rec["mid"])
    todo = []
    for m in iter_maps():
        mid = rumsey_mid(m)
        if mid and mid not in done and len(m["gcps"]) >= MIN_GCPS and gcp_span(m) >= MIN_SPAN:
            todo.append(mid)
            done.add(mid)
    print(f"{len(todo)} to fetch with {WORKERS} workers", flush=True)
    with open(OUT, "a", encoding="utf-8") as fh, ThreadPoolExecutor(WORKERS) as pool:
        for i, rec in enumerate(pool.map(fetch_one, todo)):
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
            fh.flush()
            if i % 200 == 0:
                print(i, rec["mid"], rec.get("Date"), flush=True)


if __name__ == "__main__":
    sys.exit(main())
