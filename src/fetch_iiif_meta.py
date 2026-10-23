"""Fetch catalogue metadata for non-Rumsey maps from their IIIF manifests.

Each Allmaps record points at an image service. Where the record names its
manifest (resource.partOf) we fetch that. CONTENTdm images get their parent
record through the CONTENTdm web services instead, since page manifests there
carry no dates. Gallica (Cloudflare challenge) and NYPL (no manifest link) are
skipped. Output: data/raw/iiif_meta.jsonl, one line per manifest or record.
"""
import json
import re
import socket
import ssl
import urllib.error
import sys
import threading
import time
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlparse

from common import RAW, USER_AGENT, domain, gcp_span, iter_maps

OUT = RAW / "iiif_meta.jsonl"
LINKS = RAW / "iiif_links.jsonl"
MIN_GCPS, MIN_SPAN = 6, 1.0
socket.setdefaulttimeout(30)  # a half-closed connection once hung a worker for good
_LOCK = threading.Lock()
_FH = None
SKIP = {"gallica.bnf.fr", "openapi.bnf.fr", "iiif.nypl.org", "www.davidrumsey.com"}
CDM = re.compile(r"https?://(cdm\d+\.contentdm\.oclc\.org)/(?:digital/)?iiif/(?:2/)?([^:/]+)[:/](\d+)")


def manifest_of(m):
    for c in m["resource"].get("partOf") or []:
        if c.get("type") == "Manifest":
            return c["id"]
        for mm in c.get("partOf") or []:
            if mm.get("type") == "Manifest":
                return mm["id"]
    return None


def source_of(m):
    """Return (kind, key) naming where this map's metadata lives, or None."""
    rid = m["resource"]["id"]
    hit = CDM.match(rid)
    if hit:
        return ("cdm", f"{hit.group(1)}|{hit.group(2)}|{hit.group(3)}")
    man = manifest_of(m)
    return ("manifest", man) if man else None


try:  # Windows Python lacks some intermediate certificates (Stanford, Princeton)
    import certifi
    _SSL = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    _SSL = None


def get_json(url):
    # no Accept header: Digital Commonwealth answers 404 to "Accept: application/json"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=30, context=_SSL) as r:
        return json.loads(r.read().decode("utf-8", errors="replace"))


def text_of(v):
    """Flatten IIIF v2/v3 language maps and value lists into one string."""
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    if isinstance(v, list):
        return " ; ".join(text_of(x) for x in v)
    if isinstance(v, dict):
        if "@value" in v:
            return text_of(v["@value"])
        return " ; ".join(text_of(x) for x in v.values())
    return str(v)


def fetch_manifest(url):
    j = get_json(url)
    meta = [{"label": text_of(kv.get("label")), "value": text_of(kv.get("value"))[:400]}
            for kv in (j.get("metadata") or [])]
    return {"label": text_of(j.get("label"))[:400], "navDate": j.get("navDate"), "metadata": meta}


_cdm_fields = {}


def fetch_cdm(key):
    host, alias, ptr = key.split("|")
    base = f"https://{host}/digital/bl/dmwebservices/index.php?q="
    if (host, alias) not in _cdm_fields:
        info = get_json(f"{base}dmGetCollectionFieldInfo/{alias}/json")
        _cdm_fields[(host, alias)] = {f["nick"]: f["name"] for f in info}
    names = _cdm_fields[(host, alias)]
    parent = get_json(f"{base}GetParent/{alias}/{ptr}/json").get("parent", "-1")
    rec_ptr = ptr if str(parent) in ("-1", "") else str(parent)
    item = get_json(f"{base}dmGetItemInfo/{alias}/{rec_ptr}/json")
    meta = [{"label": names.get(k, k), "value": v[:400]} for k, v in item.items()
            if isinstance(v, str) and v.strip()]
    return {"label": item.get("title", "") if isinstance(item.get("title"), str) else "",
            "navDate": None, "metadata": meta, "record": rec_ptr}


def run_host(tasks):
    """Fetch one host's records in turn, writing each as soon as it arrives."""
    n = 0
    for kind, key in tasks:
        rec = {"kind": kind, "key": key}
        for attempt in range(3):
            try:
                rec.update(fetch_cdm(key) if kind == "cdm" else fetch_manifest(key))
                rec.pop("error", None)
                break
            except urllib.error.HTTPError as e:
                rec["error"] = str(e)[:200]
                if 400 <= e.code < 500:  # permanent: missing or forbidden, do not retry
                    break
                time.sleep(3 * (attempt + 1))
            except Exception as e:
                rec["error"] = str(e)[:200]
                time.sleep(3 * (attempt + 1))
        with _LOCK:
            _FH.write(json.dumps(rec, ensure_ascii=False) + "\n")
            _FH.flush()
        n += 1
        time.sleep(0.3)
    return n, tasks[0][1] if tasks else ""


def main():
    links, keys = [], {}
    for m in iter_maps():
        d = domain(m)
        if d in SKIP or len(m["gcps"]) < MIN_GCPS or gcp_span(m) < MIN_SPAN:
            continue
        src = source_of(m)
        links.append({"map": m["id"], "source": src})
        if src:
            keys[src] = urlparse(src[1] if src[0] == "manifest" else "https://" + src[1].split("|")[0]).netloc
    with open(LINKS, "w", encoding="utf-8") as fh:
        for l in links:
            fh.write(json.dumps(l) + "\n")
    done = set()
    if OUT.exists():
        for line in open(OUT, encoding="utf-8"):
            r = json.loads(line)
            if "error" not in r:
                done.add((r["kind"], r["key"]))
    by_host = defaultdict(list)
    for src, host in keys.items():
        if src not in done:
            by_host[host].append(src)
    print(f"{len(links)} maps, {len(keys)} sources, {sum(map(len, by_host.values()))} to fetch "
          f"across {len(by_host)} hosts", flush=True)
    global _FH
    with open(OUT, "a", encoding="utf-8") as _FH, ThreadPoolExecutor(8) as pool:
        for n, first in pool.map(run_host, by_host.values()):
            print("host done", n, first[:60], flush=True)


if __name__ == "__main__":
    sys.exit(main())
