"""PLAN2 part B: separate the page term a_p from the compass term a_c on every map.

B0 synthetic known answers first, then B1-B3 on the real primary maps with |S| < 0.15.
Also writes the per-point table (data/derived/r2_points.csv) used by part F.
"""
import pickle
from functools import partial

import numpy as np
import pandas as pd

from analyze import load
from common import DERIVED, RESULTS
from fit import fit_map
from paper import contrasts, decompose, forward_fit, paper_residuals
from stats import boot


def point_table(df, gcps):
    rows = []
    for g in df.itertuples():
        pix, geo = gcps[g.map_id]
        r = fit_map(pix, geo)
        if r is None or abs(r["S"]) >= 0.15:
            continue
        res, theta = paper_residuals(r)
        p, q = contrasts(res, theta)
        rows.append(pd.DataFrame({"map_id": g.map_id, "year": g.year, "cluster": g.cluster,
                                  "region": g.region, "dist_km": g.dist_km, "p": p, "q": q, "theta": theta}))
    return pd.concat(rows, ignore_index=True)


def est(d):
    return tuple(decompose(d.p.values, d.q.values, d.theta.values))


def diff_ac(d, sel_a, sel_b):
    return est(d[sel_a(d)])[1] - est(d[sel_b(d)])[1]


def early(d):
    return (d.year >= 1500) & (d.year < 1750)


def late(d):
    return d.year >= 1850


def far(d):
    return (d.year < 1750) & (d.dist_km > 2000)


def near(d):
    return (d.year < 1750) & (d.dist_km < 500)


def ap_minus_ac(d):
    a = est(d)
    return a[0] - a[1]


def synthetic(df, gcps, a_p, a_c, rng):
    """Rebuild each map from its own fitted forward model, real layout and real local tilts,
    with residual covariance sigma^2 (I + a_p J0 + a_c J(theta))."""
    rows = []
    for g in df.itertuples():
        pix, geo = gcps[g.map_id]
        r = fit_map(pix, geo)
        if r is None or abs(r["S"]) >= 0.15:
            continue
        D, coef = forward_fit(r)
        pred = D(r["gcp_lon"], r["gcp_lat"]) @ coef
        _, theta = paper_residuals(r)
        sig = np.sqrt(np.mean(paper_residuals(r)[0] ** 2))
        c, s = np.cos(2 * theta), np.sin(2 * theta)
        S = np.empty((len(theta), 2, 2))
        S[:, 0, 0] = 1 + a_p + a_c * c
        S[:, 1, 1] = 1 - a_p - a_c * c
        S[:, 0, 1] = S[:, 1, 0] = a_c * s
        noise = np.einsum("nij,nj->ni", np.linalg.cholesky(S), rng.standard_normal((len(theta), 2))) * sig
        r2 = fit_map(pred + noise, np.column_stack([r["gcp_lon"], r["gcp_lat"]]))
        if r2 is None:
            continue
        res, th = paper_residuals(r2)
        p, q = contrasts(res, th)
        rows.append(pd.DataFrame({"p": p, "q": q, "theta": th}))
    d = pd.concat(rows, ignore_index=True)
    return est(d)


def fmt(p, ci):
    return f"{p:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]"


def main():
    gcps = pickle.load(open(DERIVED / "gcps.pkl", "rb"))
    st = pd.read_csv(RESULTS / "map_stats.csv")
    st = st[st.extent_deg >= 4][["map_id", "region", "dist_km", "cluster"]]
    df = load().merge(st, on="map_id")
    df = df[df.eligible]
    L = ["ROUND 2, PART B: page term a_p and compass term a_c", ""]

    rng = np.random.default_rng(23)
    ok = True
    L.append("B0 synthetic known answers (real layouts and tilts; tolerance 0.04, leak below 0.04):")
    for a_p, a_c in [(0.2, 0.0), (0.0, 0.2), (0.2, 0.2)]:
        e_p, e_c = synthetic(df, gcps, a_p, a_c, rng)
        good = abs(e_p - a_p) <= 0.04 and abs(e_c - a_c) <= 0.04
        ok &= good
        L.append(f"   injected a_p {a_p:.1f}, a_c {a_c:.1f}: recovered a_p {e_p:+.3f}, a_c {e_c:+.3f} -> {'PASS' if good else 'FAIL'}")
    L.append(f"   B0 overall: {'PASS' if ok else 'FAIL'}")
    L.append("")

    pts = point_table(df, gcps)
    pts.to_csv(DERIVED / "r2_points.csv", index=False)
    tilt = np.degrees(np.abs(np.arctan2(np.sin(pts.theta), np.cos(pts.theta))))
    L.append(f"points {len(pts)} on {pts.map_id.nunique()} maps with |S| < 0.15; local tilt of drawn east: "
             f"median {np.median(tilt):.1f} deg, 90th percentile {np.percentile(tilt, 90):.1f}, "
             f"share over 10 deg {np.mean(tilt > 10):.2f}")
    L.append("")
    L.append("By period (a_p page, a_c compass; ratio form exp(atanh a)):")
    rows = []
    for lo, hi in [(1450, 1650), (1650, 1750), (1750, 1850), (1850, 1951)]:
        g = pts[(pts.year >= lo) & (pts.year < hi)]
        (ap, ac), cis = boot(g, est, g.cluster, parallel=True, seed=lo)
        rows.append({"lo": lo, "hi": hi, "n_points": len(g), "n_maps": g.map_id.nunique(),
                     "a_p": ap, "a_p_lo": cis[0][0], "a_p_hi": cis[0][1],
                     "a_c": ac, "a_c_lo": cis[1][0], "a_c_hi": cis[1][1]})
        L.append(f"   {lo}-{hi - 1}: maps {g.map_id.nunique():5d}  a_p {fmt(ap, cis[0])}  a_c {fmt(ac, cis[1])}")
    for reg_name, sel in [("Europe", pts.region == "Europe"), ("elsewhere", pts.region != "Europe")]:
        for lo, hi in [(1450, 1750), (1750, 1850), (1850, 1951)]:
            g = pts[sel & (pts.year >= lo) & (pts.year < hi)]
            if g.map_id.nunique() < 20:
                continue
            (ap, ac), cis = boot(g, est, g.cluster, parallel=True, seed=lo + 7)
            rows.append({"lo": lo, "hi": hi, "group": reg_name, "n_points": len(g), "n_maps": g.map_id.nunique(),
                         "a_p": ap, "a_p_lo": cis[0][0], "a_p_hi": cis[0][1],
                         "a_c": ac, "a_c_lo": cis[1][0], "a_c_hi": cis[1][1]})
            L.append(f"   {reg_name:9s} {lo}-{hi - 1}: maps {g.map_id.nunique():5d}  a_p {fmt(ap, cis[0])}  a_c {fmt(ac, cis[1])}")
    for dlo, dhi in [(0, 500), (500, 2000), (2000, 6000), (6000, 25000)]:
        g = pts[(pts.year < 1750) & (pts.dist_km >= dlo) & (pts.dist_km < dhi)]
        if g.map_id.nunique() < 15:
            continue
        (ap, ac), cis = boot(g, est, g.cluster, parallel=True, seed=dlo + 3)
        rows.append({"lo": 1450, "hi": 1750, "group": f"dist {dlo}-{dhi}", "n_points": len(g),
                     "n_maps": g.map_id.nunique(), "a_p": ap, "a_p_lo": cis[0][0], "a_p_hi": cis[0][1],
                     "a_c": ac, "a_c_lo": cis[1][0], "a_c_hi": cis[1][1]})
        L.append(f"   before 1750, {dlo}-{dhi} km: maps {g.map_id.nunique():4d}  a_p {fmt(ap, cis[0])}  a_c {fmt(ac, cis[1])}")
    pd.DataFrame(rows).to_csv(RESULTS / "r2_decomposition.csv", index=False)
    L.append("")

    lt = pts[late(pts)]
    (ap, ac), cis = boot(lt, est, lt.cluster, parallel=True, seed=1)
    dpc, dci = boot(lt, ap_minus_ac, lt.cluster, parallel=True, seed=2)
    L.append(f"B1 1850-1950: a_p {fmt(ap, cis[0])}, a_c {fmt(ac, cis[1])}, a_p - a_c {fmt(dpc, dci)} -> "
             f"{'SUPPORTED' if cis[0][0] > 0 and dci[0] > 0 else 'NOT SUPPORTED'}")
    b2, b2ci = boot(pts, partial(diff_ac, sel_a=early, sel_b=late), pts.cluster, parallel=True, seed=3)
    L.append(f"B2 a_c 1500-1749 minus a_c 1850-1950: {fmt(b2, b2ci)} -> {'SUPPORTED' if b2ci[0] > 0 else 'NOT SUPPORTED'}")
    b3, b3ci = boot(pts, partial(diff_ac, sel_a=far, sel_b=near), pts.cluster, parallel=True, seed=4)
    L.append(f"B3 before 1750, a_c over 2,000 km minus within 500 km: {fmt(b3, b3ci)} -> "
             f"{'SUPPORTED' if b3ci[0] > 0 else 'NOT SUPPORTED'}")
    text = "\n".join(L)
    (RESULTS / "r2_decomposition.txt").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
