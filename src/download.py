"""Download the inputs: the Allmaps open dataset, the GeoNames gazetteer and Natural
Earth land. The Allmaps file is regenerated twice a day upstream, so a new download
will hold more maps than the 2026-10-13 snapshot used for the README numbers."""
import urllib.request
import zipfile

from common import RAW, USER_AGENT

FILES = {
    "maps.ndjson": "https://files.allmaps.org/maps.ndjson",
    "cities15000.zip": "https://download.geonames.org/export/dump/cities15000.zip",
    "ne/ne_50m_land.zip": "https://naciscdn.org/naturalearth/50m/physical/ne_50m_land.zip",
}

for name, url in FILES.items():
    out = RAW / name
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        print("have", name)
        continue
    print("get", url, flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=600) as r, open(out, "wb") as fh:
        fh.write(r.read())
zipfile.ZipFile(RAW / "cities15000.zip").extractall(RAW)
