"""Shared paths and loaders for the Longitude Wind pipeline."""
import json
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
DOCS = ROOT / "docs"

MAPS_NDJSON = RAW / "maps.ndjson"
USER_AGENT = "longitude-wind research script (historical map accuracy study)"

for d in (DERIVED, RESULTS, FIGURES, DOCS):
    d.mkdir(parents=True, exist_ok=True)


def iter_maps(path=MAPS_NDJSON):
    """Yield every georeferenced map record in the Allmaps dump."""
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def domain(m):
    return urlparse(m["resource"]["id"]).netloc


def gcp_span(m):
    """Largest of the longitude and latitude extents of a map's GCPs, in degrees."""
    g = m["gcps"]
    lons = [p["geo"][0] for p in g]
    lats = [p["geo"][1] for p in g]
    return max(max(lons) - min(lons), max(lats) - min(lats))


def rumsey_mid(m):
    """Luna media id (RUMSEY~8~1~x~y) for a Rumsey image service URL."""
    rid = m["resource"]["id"]
    if "davidrumsey.com" not in rid:
        return None
    return rid.rstrip("/").split("/")[-1]
