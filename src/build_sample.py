"""Geometric sample: every Allmaps map with enough GCPs and a regional extent.

Writes data/derived/maps_geom.csv (one row per kept map) and data/derived/gcps.pkl
({map id: (pixels, lon/lat)}). Dates are joined later in metadata.py.
"""
import hashlib
import pickle

import numpy as np
import pandas as pd

from common import DERIVED, domain, iter_maps, rumsey_mid
from fit import centroid, dedupe

MIN_GCPS = 6
MIN_EXTENT = 1.0          # degrees; PLAN primary threshold (4) and sensitivity (2, 8) applied later
MAX_LON_SPAN, MAX_LAT_SPAN = 120.0, 70.0


def geo_signature(geo):
    """Order-free fingerprint of a set of control-point coordinates."""
    o = np.lexsort((geo[:, 1], geo[:, 0]))
    return hashlib.md5(np.round(geo[o], 4).tobytes()).hexdigest()


def spans(lon, lat):
    """Longitude span measured around the centroid, so maps across 180 E/W are not
    mistaken for world maps; latitude span is plain."""
    lon0, lat0 = centroid(lon, lat)
    rel = (lon - lon0 + 180.0) % 360.0 - 180.0
    return float(rel.max() - rel.min()), float(lat.max() - lat.min()), lon0, lat0


def main():
    rows, gcps, best = [], {}, {}
    for m in iter_maps():
        pix = np.array([g["resource"] for g in m["gcps"]], float)
        geo = np.array([g["geo"] for g in m["gcps"]], float)
        if len(pix) < MIN_GCPS:
            continue
        pix, geo = dedupe(pix, geo)
        if len(pix) < MIN_GCPS:
            continue
        lon_span, lat_span, lon0, lat0 = spans(geo[:, 0], geo[:, 1])
        extent = max(lon_span, lat_span)
        if extent < MIN_EXTENT or lon_span > MAX_LON_SPAN or lat_span > MAX_LAT_SPAN:
            continue
        image = m["resource"]["id"]
        # one record per image: keep the georeference with the most GCPs
        if image in best and best[image][1] >= len(pix):
            continue
        if image in best:
            gcps.pop(best[image][0], None)
        best[image] = (m["id"], len(pix))
        gcps[m["id"]] = (pix, geo)
        rows.append({
            "map_id": m["id"], "image": image, "provider": domain(m), "rumsey_mid": rumsey_mid(m),
            "width": m["resource"].get("width"), "height": m["resource"].get("height"), "n_gcps": len(pix), "lon_span": lon_span, "lat_span": lat_span, "extent_deg": extent,
            "lon0": lon0, "lat0": lat0, "transformation": m["transformation"]["type"],
            "created": m.get("created"),
        })
    df = pd.DataFrame(rows)
    df = df[df["map_id"].isin(gcps)].reset_index(drop=True)
    # the same scan can sit under two image URLs (VU Amsterdam serves each map twice), so its
    # georeference appears twice; keep one record per identical set of control points
    df["geo_sig"] = [geo_signature(gcps[m][1]) for m in df.map_id]
    dup = df.duplicated("geo_sig")
    for m in df.map_id[dup]:
        gcps.pop(m, None)
    print(f"dropped {int(dup.sum())} records whose control points repeat another record's exactly")
    df = df[~dup].reset_index(drop=True)
    df.to_csv(DERIVED / "maps_geom.csv", index=False)
    with open(DERIVED / "gcps.pkl", "wb") as fh:
        pickle.dump(gcps, fh)
    print(f"{len(df)} maps kept; extent>=4: {(df.extent_deg >= 4).sum()}; providers:")
    print(df.provider.value_counts().head(15).to_string())


if __name__ == "__main__":
    main()
