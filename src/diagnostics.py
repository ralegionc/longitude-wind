"""Post-hoc diagnostics, written after the pre-registered run (see PLAN.md addendum).

K2 failed: maps from 1900 to 1950 still show east-west errors about 14% larger than
north-south ones, in every region alike. One candidate is the book itself: a map
printed across two pages has a vertical fold (the gutter) that hides or stretches a
strip of the map, which can only displace places east-west. This script measures the
residuals as a function of position across the scanned page.

D1  step test: mean east-west residual left of the page centre minus right of it,
    in units of the map's own RMS error, with three controls (north-south residual
    across the same line, both residuals across the horizontal centre line).
D2  profile: |east-west| and |north-south| residual against position across the page.
"""
import json
import pickle

import numpy as np
import pandas as pd

from analyze import load
from common import DERIVED, RESULTS
from fit import aeqd_forward, design, fit_map
from stats import boot


def kind(t):
    if not isinstance(t, str) or t in ("null", "nan"):
        return "unknown"
    try:
        v = json.loads(t)
    except Exception:
        v = t
    v = " ".join(map(str, v if isinstance(v, list) else [v]))
    return "atlas" if "Atlas Map" in v else ("sheet" if any(k in v for k in
                                                               ("Separate Map", "Case Map", "Wall Map", "Pocket Map")) else "other")


def gutter_rows(df, gcps, fit_kwargs=None):
    rows, prof = [], []
    for g in df.itertuples():
        pix, geo = gcps[g.map_id]
        r = fit_map(pix, geo, **(fit_kwargs or {}))
        if r is None or not g.width or not g.height:
            continue
        x, y = r["gcp_u"] / g.width, r["gcp_v"] / g.height
        rms = r["rms_km"]
        ew, ns = r["e_ew"] / rms, r["e_ns"] / rms
        row = {"map_id": g.map_id, "year": g.year, "kind": kind(g.type), "maker": g.maker, "A": r["A"]}
        for name, pos, line in [("v", x, 0.5), ("h", y, 0.5), ("v30", x, 0.3)]:
            left, right = pos < line, pos >= line
            if left.sum() >= 3 and right.sum() >= 3:
                row[f"ew_{name}"] = ew[left].mean() - ew[right].mean()
                row[f"ns_{name}"] = ns[left].mean() - ns[right].mean()
        rows.append(row)
        prof.append(pd.DataFrame({"map_id": g.map_id, "kind": row["kind"], "x": x, "y": y,
                                  "aew": np.abs(ew), "ans": np.abs(ns), "ew": ew, "ns": ns}))
    return pd.DataFrame(rows), pd.concat(prof)


def frames(r):
    """Residuals on the paper, split two ways: along the page's own axes (horizontal,
    vertical) and along the drawn compass directions (local east, local north).
    Returns (A_page, A_compass, mean tilt of drawn north from page-up in degrees)."""
    lon, lat, u, v = r["gcp_lon"], r["gcp_lat"], r["gcp_u"], r["gcp_v"]
    X, Y = aeqd_forward(lon, lat, r["lon0"], r["lat0"])
    mu, sc = np.array([X.mean(), Y.mean()]), np.std(np.r_[X, Y]) or 1.0
    order = r["order"]

    def D(lo, la):
        x, y = aeqd_forward(lo, la, r["lon0"], r["lat0"])
        return design((x - mu[0]) / sc, (y - mu[1]) / sc, order)

    coef, *_ = np.linalg.lstsq(D(lon, lat), np.column_stack([u, v]), rcond=None)
    res = np.column_stack([u, v]) - D(lon, lat) @ coef
    d = 0.05
    north = D(lon, lat + d) @ coef - D(lon, lat) @ coef
    east = D(lon + d / np.cos(np.radians(lat)), lat) @ coef - D(lon, lat) @ coef
    north /= np.linalg.norm(north, axis=1, keepdims=True)
    east /= np.linalg.norm(east, axis=1, keepdims=True)
    a_page = 0.5 * np.log(np.sum(res[:, 0] ** 2) / np.sum(res[:, 1] ** 2))
    a_comp = 0.5 * np.log(np.sum(np.sum(res * east, 1) ** 2) / np.sum(np.sum(res * north, 1) ** 2))
    tilt = np.degrees(np.arctan2(north[:, 0], -north[:, 1]))  # 0 = north straight up the page
    return a_page, a_comp, float(np.degrees(np.arctan2(np.mean(np.sin(np.radians(tilt))), np.mean(np.cos(np.radians(tilt))))))


def frame_rows(df, gcps):
    rows = []
    for g in df.itertuples():
        pix, geo = gcps[g.map_id]
        r = fit_map(pix, geo)
        if r is None:
            continue
        try:
            a_page, a_comp, tilt = frames(r)
        except (FloatingPointError, ZeroDivisionError, np.linalg.LinAlgError):
            continue
        rows.append({"map_id": g.map_id, "year": g.year, "maker": g.maker, "A": r["A"],
                     "A_page": a_page, "A_compass": a_comp, "tilt": tilt})
    return pd.DataFrame(rows)


def main():
    gcps = pickle.load(open(DERIVED / "gcps.pkl", "rb"))
    df = load()
    df = df[df.eligible & (df.extent_deg >= 4)]
    fr = frame_rows(df, gcps)
    fr["cluster"] = fr.maker.fillna(fr.map_id).astype(str)
    fr.to_csv(RESULTS / "frames.csv", index=False)
    d, prof = gutter_rows(df, gcps)
    d["cluster"] = d.maker.fillna(d.map_id).astype(str)
    prof.to_csv(DERIVED / "page_profile.csv", index=False)
    d.to_csv(RESULTS / "gutter_steps.csv", index=False)
    L = ["POST-HOC DIAGNOSTICS (not pre-registered)", "",
         "D1 step across a line on the page, in units of each map's RMS error (mean, 95% cluster bootstrap):"]
    for k in ["atlas", "sheet", "unknown"]:
        g = d[d.kind == k]
        L.append(f"  {k} maps (n={len(g)})")
        for col, label in [("ew_v", "east-west residual across the vertical centre line"),
                           ("ns_v", "north-south residual across the vertical centre line"),
                           ("ew_h", "east-west residual across the horizontal centre line"),
                           ("ns_h", "north-south residual across the horizontal centre line"),
                           ("ew_v30", "east-west residual across a line at 30% of the width")]:
            gg = g[g[col].notna()]
            p, ci = boot(gg, lambda z, c=col: z[c].mean(), gg.cluster, n=1000)
            L.append(f"    {label:55s} n={len(gg):5d} {p:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]")
    L.append("")
    L.append("D2 mean |residual| / RMS by position across the page (x = 0 left edge, 1 right edge):")
    prof["xb"] = (prof.x.clip(0, 0.999) * 10).astype(int) / 10 + 0.05
    t = prof[prof.x.between(0, 1)].groupby(["kind", "xb"])[["aew", "ans"]].mean().unstack(0).round(3)
    L.append(t.to_string())
    L.append("")
    L.append("D3 compass or page? Anisotropy measured on the paper along the page axes (A_page) and along")
    L.append("   the drawn east and north (A_compass), by how far drawn north is tilted from page-up:")
    fr["abs_tilt"] = fr.tilt.abs()
    for lo, hi in [(0, 5), (5, 15), (15, 30), (30, 60), (60, 120), (120, 181)]:
        g = fr[(fr.abs_tilt >= lo) & (fr.abs_tilt < hi)]
        if len(g) < 5:
            continue
        p1, c1 = boot(g, lambda z: z.A_page.mean(), g.cluster, n=1000)
        p2, c2 = boot(g, lambda z: z.A_compass.mean(), g.cluster, n=1000)
        L.append(f"   tilt {lo:3d}-{hi:3d} deg  n={len(g):5d}  A_page {p1:+.3f} [{c1[0]:+.3f}, {c1[1]:+.3f}]   "
                 f"A_compass {p2:+.3f} [{c2[0]:+.3f}, {c2[1]:+.3f}]")
    L.append("")
    L.append("Descriptive: mean A (and map count) by region and period, primary sample:")
    st = pd.read_csv(RESULTS / "map_stats.csv")
    st = st[st.extent_deg >= 4].copy()
    st["period"] = pd.cut(st.year, [1449, 1649, 1749, 1849, 1951],
                          labels=["1450-1649", "1650-1749", "1750-1849", "1850-1950"])
    st["group"] = np.where(st.region == "Europe", "Europe",
                           np.where(st.region == "North America", "North America", "Asia, Africa, Lat. Am., Oceania"))
    t = st.groupby(["group", "period"], observed=True).A.agg(["mean", "count"]).round(3).unstack()
    L.append(t.to_string())
    text = "\n".join(L)
    (RESULTS / "diagnostics.txt").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
