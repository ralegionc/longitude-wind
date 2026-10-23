"""PLAN2 part F: replication on maps round 1 never saw.

MapWarper.net maps (independent volunteers, a different clicking tool) and Allmaps maps
created after the round-1 snapshot, under the round-1 filters and fit. Clusters are
single maps (MapWarper has no maker field).
"""
import json

import numpy as np
import pandas as pd

from analyze import region
from build_sample import MAX_LAT_SPAN, MAX_LON_SPAN, spans
from common import RAW, RESULTS, iter_maps
from dates import parse_year
from fit import dedupe, fit_map
from metadata import DATE_LABEL, FACSIMILE, NOT_DATE, VIEW_GLOBE
from paper import contrasts, decompose, paper_residuals
from stats import boot


def mapwarper_records():
    meta = {}
    for line in open(RAW / "mapwarper_maps.jsonl", encoding="utf-8"):
        r = json.loads(line)
        if r.get("id"):
            meta[r["id"]] = r
    for line in open(RAW / "mapwarper_gcps.jsonl", encoding="utf-8"):
        g = json.loads(line)
        if "error" in g or not g.get("gcps"):
            continue
        m = meta[g["id"]]
        y = None
        for text in (m.get("date_depicted"), m.get("title")):
            y = parse_year(text)
            if y:
                break
        text = f"{m.get('title') or ''} {m.get('description') or ''}"
        pts = np.array(g["gcps"], float)
        yield {"source": "mapwarper", "map_id": f"mapwarper:{g['id']}", "title": m.get("title"),
               "year": y[0] if y else np.nan, "prec": y[1] if y else np.nan,
               "excluded": bool(FACSIMILE.search(text) or VIEW_GLOBE.search(m.get("title") or "")),
               "pix": pts[:, :2], "geo": pts[:, 2:]}


def new_allmaps_records():
    meta = {}
    for line in open(RAW / "r2_new_meta.jsonl", encoding="utf-8"):
        r = json.loads(line)
        meta[r["map_id"]] = r
    for m in iter_maps(RAW / "maps_2026-10-23.ndjson"):
        r = meta.get(m["id"])
        if r is None or "error" in r:
            continue
        y = None
        if r["kind"] == "rumsey":
            for k in ("Date", "Pub Date"):
                v = r.get(k)
                y = parse_year(v[0] if isinstance(v, list) else v)
                if y:
                    break
            text = " ".join(str(r.get(k) or "") for k in ("Full Title", "Pub Title", "Pub Note", "Note", "Type"))
        else:
            pairs = [(kv["label"], kv["value"]) for kv in r.get("metadata", [])]
            for lab, val in pairs:
                if DATE_LABEL.search(lab) and not NOT_DATE.search(lab):
                    y = parse_year(val)
                    if y:
                        break
            if not y:
                y = parse_year(r.get("label"))
            text = " ".join([r.get("label") or ""] + [v for _, v in pairs])
        yield {"source": "allmaps-new", "map_id": m["id"], "title": r.get("Full Title") or r.get("label"),
               "year": y[0] if y else np.nan, "prec": y[1] if y else np.nan,
               "excluded": bool(FACSIMILE.search(text)),
               "pix": np.array([g["resource"] for g in m["gcps"]], float),
               "geo": np.array([g["geo"] for g in m["gcps"]], float)}


def mean_a(d):
    return d.A.mean()


def f2(d):
    o = d[d.region != "Europe"]
    return o[o.year < 1850].A.mean() - o[o.year >= 1850].A.mean()


def est(d):
    return tuple(decompose(d.p.values, d.q.values, d.theta.values))


def fmt(p, ci):
    return f"{p:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]"


def main():
    rows, pts = [], []
    counts = {"mapwarper": 0, "allmaps-new": 0}
    for rec in list(mapwarper_records()) + list(new_allmaps_records()):
        counts[rec["source"]] += 1
        if rec["excluded"] or not (1450 <= rec["year"] <= 1950) or rec["prec"] > 25:
            continue
        pix, geo = dedupe(rec["pix"], rec["geo"])
        if len(pix) < 6:
            continue
        lon_span, lat_span, lon0, lat0 = spans(geo[:, 0], geo[:, 1])
        if max(lon_span, lat_span) < 4 or lon_span > MAX_LON_SPAN or lat_span > MAX_LAT_SPAN:
            continue
        r = fit_map(pix, geo)
        if r is None:
            continue
        reg = region(r["lon0"], r["lat0"])
        rows.append({"source": rec["source"], "map_id": rec["map_id"], "title": rec["title"], "year": rec["year"],
                     "region": reg, "A": r["A"], "S": r["S"], "n": r["n"], "extent_deg": max(lon_span, lat_span)})
        if abs(r["S"]) < 0.15:
            res, theta = paper_residuals(r)
            p, q = contrasts(res, theta)
            pts.append(pd.DataFrame({"map_id": rec["map_id"], "year": rec["year"], "p": p, "q": q, "theta": theta}))
    d = pd.DataFrame(rows)
    d.to_csv(RESULTS / "r2_replication_maps.csv", index=False)
    pts = pd.concat(pts, ignore_index=True)
    L = ["ROUND 2, PART F: replication on maps round 1 never saw", "",
         f"records read: MapWarper {counts['mapwarper']} (dated regional maps with GCPs), new Allmaps {counts['allmaps-new']}",
         f"eligible maps (round-1 filters, extent >= 4 degrees): {len(d)} "
         f"({(d.source == 'mapwarper').sum()} MapWarper, {(d.source == 'allmaps-new').sum()} new Allmaps)",
         "by century: " + ", ".join(f"{c}s {((d.year >= c) & (d.year < c + 100)).sum()}" for c in range(1400, 2000, 100)),
         "by region: " + ", ".join(f"{k} {v}" for k, v in d.region.value_counts().items()), ""]
    late = d[d.year >= 1850]
    p, ci = boot(late, mean_a, late.map_id)
    L.append(f"F1 1850-1950 mean A (n={len(late)}): {fmt(p, ci)} -> {'REPLICATED' if ci[0] > 0 else 'NOT REPLICATED'}")
    p, ci = boot(d, f2, d.map_id)
    o = d[d.region != "Europe"]
    L.append(f"F2 outside Europe, mean A before 1850 (n={(o.year < 1850).sum()}) minus after (n={(o.year >= 1850).sum()}): "
             f"{fmt(p, ci)} -> {'REPLICATED' if ci[0] > 0 else 'NOT REPLICATED'}")
    eu = d[(d.region == "Europe") & (d.year < 1700)]
    if len(eu) >= 10:
        p, ci = boot(eu, mean_a, eu.map_id)
        L.append(f"F3 Europe before 1700 (n={len(eu)}): mean A {fmt(p, ci)} -> "
                 f"{'REPLICATED (no excess)' if ci[0] <= 0 <= ci[1] else 'NOT REPLICATED'}")
    else:
        L.append(f"F3 Europe before 1700: only {len(eu)} maps, not tested")
    lp = pts[pts.year >= 1850]
    (ap, ac), cis = boot(lp, est, lp.map_id)
    L.append(f"F4 decomposition 1850-1950 ({lp.map_id.nunique()} maps): a_p {fmt(ap, cis[0])}, a_c {fmt(ac, cis[1])} -> "
             f"{'REPLICATED' if cis[0][0] > 0 else 'NOT REPLICATED'}")
    text = "\n".join(L)
    (RESULTS / "r2_replication.txt").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
