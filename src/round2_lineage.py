"""PLAN2 part E: copy lineages from shared error patterns, and H3 with better statistics.

E1 two maps sharing at least 5 places (within 15 km) get the correlation of their
stacked east and north residuals at those places. The null is pairs at least 150 years
apart; a pair is a copy if it beats the null's 99th percentile. Lineages are connected
components of the copy graph.
E2 H3 with clusters = lineages merged with makers.
E3 H3 as a mixed model (random intercepts for maker, variance component for atlas), plus
a dof-weighted regression with two-way clustered errors because statsmodels' mixed
model takes no weights (see PLAN2 addendum).
"""
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from common import DERIVED, RAW, RESULTS
from fit import R
from stats import boot


def unit(lon, lat):
    lam, phi = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)])


def pair_correlations(err, radius_km=15.0, min_shared=5):
    err = err.reset_index(drop=True)
    codes, names = pd.factorize(err.map_id)
    tree = cKDTree(unit(err.lon.values, err.lat.values))
    pairs = tree.query_pairs(radius_km / R, output_type="ndarray")
    a, b = pairs[:, 0], pairs[:, 1]
    ca, cb = codes[a], codes[b]
    keep = ca != cb
    a, b, ca, cb = a[keep], b[keep], ca[keep], cb[keep]
    swap = ca > cb  # order each map pair consistently
    a, b = np.where(swap, b, a), np.where(swap, a, b)
    ca, cb = np.where(swap, cb, ca), np.where(swap, ca, cb)
    pa, pb = unit(err.lon.values[a], err.lat.values[a]), unit(err.lon.values[b], err.lat.values[b])
    df = pd.DataFrame({"ca": ca, "cb": cb, "a": a, "b": b, "d": np.linalg.norm(pa - pb, axis=1)})
    df = df.sort_values("d").drop_duplicates(["ca", "cb", "a"]).drop_duplicates(["ca", "cb", "b"])
    ew_a, ew_b = err.e_ew.values[df.a], err.e_ew.values[df.b]
    ns_a, ns_b = err.e_ns.values[df.a], err.e_ns.values[df.b]
    g = pd.DataFrame({"ca": np.r_[df.ca, df.ca], "cb": np.r_[df.cb, df.cb],
                      "x": np.r_[ew_a, ns_a], "y": np.r_[ew_b, ns_b]})
    g["xx"], g["yy"], g["xy"] = g.x ** 2, g.y ** 2, g.x * g.y
    s = g.groupby(["ca", "cb"]).agg(n=("x", "size"), sx=("x", "sum"), sy=("y", "sum"),
                                    sxx=("xx", "sum"), syy=("yy", "sum"), sxy=("xy", "sum"))
    s = s[s.n >= 2 * min_shared]
    cov = s.sxy - s.sx * s.sy / s.n
    vx, vy = s.sxx - s.sx ** 2 / s.n, s.syy - s.sy ** 2 / s.n
    out = pd.DataFrame({"r": cov / np.sqrt(vx * vy), "n": s.n // 2}).reset_index()
    out["ma"], out["mb"] = names[out.ca.values], names[out.cb.values]
    return out.drop(columns=["ca", "cb"])


def h3_diff(d):
    def slope(x):
        if len(x) < 10:
            return np.nan
        return np.polyfit(np.log10(x.dist_km.clip(lower=10)), x.A, 1)[0]
    return slope(d[d.year < 1750]) - slope(d[d.year > 1850])


def main():
    st = pd.read_csv(RESULTS / "map_stats.csv")
    P = st[st.extent_deg >= 4].copy()
    err = pd.read_csv(DERIVED / "gcp_errors.csv")
    err = err[err.map_id.isin(P.map_id)]
    L = ["ROUND 2, PART E: copy lineages and better statistics", ""]
    pc = pair_correlations(err)
    yr = P.set_index("map_id").year
    pc["gap"] = (pc.ma.map(yr) - pc.mb.map(yr)).abs()
    null = pc[pc.gap >= 150].r
    thr = float(np.nanpercentile(null, 99))
    pc["is_copy"] = pc.r > thr
    pc.to_csv(DERIVED / "r2_pair_correlations.csv", index=False)
    L.append(f"E1 map pairs sharing >= 5 places: {len(pc)}; null pairs (>= 150 years apart): {len(null)}; "
             f"null 99th percentile r = {thr:.3f}")
    L.append(f"   copy pairs: {int(pc.is_copy.sum())} ({pc.is_copy.mean():.1%}); median r of all pairs {pc.r.median():.3f}; "
             f"pairs with r > 0.9: {int((pc.r > 0.9).sum())}")
    ids = {m: i for i, m in enumerate(P.map_id)}
    cp = pc[pc.is_copy & pc.ma.isin(ids) & pc.mb.isin(ids)]
    rows = np.r_[cp.ma.map(ids).values, cp.mb.map(ids).values]
    cols = np.r_[cp.mb.map(ids).values, cp.ma.map(ids).values]
    _, lineage = connected_components(coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(ids), len(ids))))
    P["lineage"] = lineage[P.map_id.map(ids).values]
    sizes = P.groupby("lineage").size()
    L.append(f"   lineages: {len(sizes)}; maps in lineages of 2+: {int(sizes[sizes > 1].sum())}; largest {int(sizes.max())}")
    # merge lineage with maker: maps join the same cluster if they share either
    mk = {m: i for i, m in enumerate(P.cluster.unique())}
    n_l = int(lineage.max()) + 1
    r2 = np.r_[P.lineage.values, P.cluster.map(mk).values + n_l]
    edges = coo_matrix((np.ones(len(P)), (P.lineage.values, P.cluster.map(mk).values + n_l)),
                       shape=(n_l + len(mk), n_l + len(mk)))
    _, comp = connected_components(edges + edges.T)
    P["cluster2"] = comp[P.lineage.values]
    del r2
    L.append(f"   clusters after merging lineages with makers: {P.cluster2.nunique()} (makers alone: {P.cluster.nunique()})")
    big = P.groupby("cluster2").size().sort_values(ascending=False).head(3)
    L.append("   largest merged clusters (maps): " + ", ".join(str(int(v)) for v in big.values))
    L.append("")

    hz = P[P.dist_km.notna()]
    p0, c0_ = boot(hz, h3_diff, hz.cluster)
    p1, c1 = boot(hz, h3_diff, hz.cluster2)
    L.append(f"E2 H3 slope difference (before 1750 minus after 1850), maker clusters: {p0:+.3f} [{c0_[0]:+.3f}, {c0_[1]:+.3f}]")
    L.append(f"   with lineages merged into clusters: {p1:+.3f} [{c1[0]:+.3f}, {c1[1]:+.3f}] -> "
             f"{'SUPPORTED' if c1[0] > 0 else 'NOT SUPPORTED'}")
    # post hoc (PLAN2 addendum): chaining at the pre-registered threshold makes one giant cluster,
    # so repeat with near-identical pairs only (r > 0.9)
    strict = pc[(pc.r > 0.9) & pc.ma.isin(ids) & pc.mb.isin(ids)]
    rr = np.r_[strict.ma.map(ids).values, strict.mb.map(ids).values]
    cc = np.r_[strict.mb.map(ids).values, strict.ma.map(ids).values]
    _, lin9 = connected_components(coo_matrix((np.ones(len(rr)), (rr, cc)), shape=(len(ids), len(ids))))
    lin9 = lin9[P.map_id.map(ids).values]
    n9 = int(lin9.max()) + 1
    e9 = coo_matrix((np.ones(len(P)), (lin9, P.cluster.map(mk).values + n9)), shape=(n9 + len(mk), n9 + len(mk)))
    _, comp9 = connected_components(e9 + e9.T)
    P["cluster3"] = comp9[lin9]
    hz = P[P.dist_km.notna()]
    big9 = P.groupby("cluster3").size().max()
    p3, c3 = boot(hz, h3_diff, hz.cluster3)
    L.append(f"   post hoc, lineages from r > 0.9 pairs merged with makers ({P.cluster3.nunique()} clusters, "
             f"largest {int(big9)} maps): {p3:+.3f} [{c3[0]:+.3f}, {c3[1]:+.3f}]")
    L.append("")

    atlas = {}
    for line in open(RAW / "rumsey_meta.jsonl", encoding="utf-8"):
        r = json.loads(line)
        v = r.get("Pub List No")
        atlas[r["mid"]] = str(v[0] if isinstance(v, list) else v) if v else None
    m = hz[(hz.year < 1750) | (hz.year > 1850)].copy()
    m["early"] = (m.year < 1750).astype(int)
    m["ld"] = np.log10(m.dist_km.clip(lower=10))
    m["atlas"] = [atlas.get(x) if isinstance(x, str) else None for x in m.rumsey_mid]
    m["atlas"] = m.atlas.fillna(m.cluster)
    try:
        md = smf.mixedlm("A ~ ld * early", m, groups=m.cluster, vc_formula={"atlas": "0 + C(atlas)"}).fit(reml=True)
        b, se, pv = md.params["ld:early"], md.bse["ld:early"], md.pvalues["ld:early"]
        L.append(f"E3 mixed model (maker random intercept, atlas variance component): slope difference {b:+.3f} "
                 f"(SE {se:.3f}, p = {pv:.3g}) -> {'SUPPORTED' if b > 0 and pv < 0.05 else 'NOT SUPPORTED'}")
    except Exception as e:  # report rather than hide a failed fit
        L.append(f"E3 mixed model failed to fit: {e}")
    X = sm.add_constant(pd.DataFrame({"ld": m.ld, "early": m.early, "ld_early": m.ld * m.early}))
    groups = np.column_stack([pd.factorize(m.cluster)[0], pd.factorize(m.atlas)[0]])
    w = sm.WLS(m.A, X, weights=m.dof).fit(cov_type="cluster", cov_kwds={"groups": groups})
    L.append(f"   dof-weighted regression, two-way clustered (maker, atlas): slope difference {w.params['ld_early']:+.3f} "
             f"(SE {w.bse['ld_early']:.3f}, p = {w.pvalues['ld_early']:.3g})")
    P[["map_id", "lineage", "cluster2"]].to_csv(RESULTS / "r2_lineages.csv", index=False)
    text = "\n".join(L)
    (RESULTS / "r2_lineage.txt").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
