"""Per-map fits, known-answer checks K2-K4, hypotheses H1-H3 and robustness (PLAN.md).

Writes results/map_stats.csv, results/summary.txt, results/robustness.csv,
results/curve.csv and data/derived/gcp_errors.csv.
"""
import pickle
import re
from functools import partial

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import DERIVED, RESULTS
from fit import R, fit_map
from stats import B, boot, fit_sigmoid

EUROPE = dict(lon=(-25, 45), lat=(34, 72))


def region(lon, lat):
    if EUROPE["lon"][0] <= lon <= EUROPE["lon"][1] and EUROPE["lat"][0] <= lat <= EUROPE["lat"][1]:
        return "Europe"
    if -170 <= lon <= -50 and lat >= 15:
        return "North America"
    if -120 <= lon <= -30 and lat < 15:
        return "Latin America"
    if -20 <= lon <= 55 and -40 <= lat < 34:
        return "Africa & Middle East"
    if 45 < lon <= 180 and lat >= -10:
        return "Asia"
    if (lon >= 110 and lat < -10) or lon <= -120:
        return "Oceania & Pacific"
    return "Other"


def great_circle_km(lat1, lon1, lat2, lon2):
    p1, p2, dl = np.radians(lat1), np.radians(lat2), np.radians(lon2 - lon1)
    c = np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(dl)
    return R * np.arccos(np.clip(c, -1, 1))


def load():
    geom = pd.read_csv(DERIVED / "maps_geom.csv")
    meta = pd.read_csv(DERIVED / "meta.csv")
    df = geom.merge(meta, on="map_id", how="left")
    df["eligible"] = (df.year.between(1450, 1950) & (df.year_prec <= 25)
                      & ~df.facsimile.fillna(False).astype(bool) & ~df.view_globe.fillna(False).astype(bool))
    return df


def run_fits(df, gcps, force_order=None, keep_errors=False):
    rows, errs = [], []
    for g in df.itertuples():
        pix, geo = gcps[g.map_id]
        r = fit_map(pix, geo, force_order=force_order)
        if r is None:
            continue
        rows.append({"map_id": g.map_id, **{k: v for k, v in r.items() if not isinstance(v, np.ndarray)}})
        if keep_errors:
            errs.append(pd.DataFrame({"map_id": g.map_id, "lon": r["gcp_lon"], "lat": r["gcp_lat"],
                                      "e_ew": r["e_ew"], "e_ns": r["e_ns"]}))
    return pd.DataFrame(rows), (pd.concat(errs) if errs else None)


def mean_of(d, col="A"):
    return d[col].mean()


def fmt(point, ci, digits=3):
    return f"{point:+.{digits}f} [{ci[0]:+.{digits}f}, {ci[1]:+.{digits}f}]"


def h1_stat(col):
    def f(d):
        early = d[(d.year >= 1500) & (d.year < 1700)][col]
        late = d[(d.year >= 1850) & (d.year <= 1950)][col]
        return early.mean() - late.mean()
    return f


def h1_stat_weighted(col):
    def f(d):
        e = d[(d.year >= 1500) & (d.year < 1700)]
        l = d[(d.year >= 1850) & (d.year <= 1950)]
        return np.average(e[col], weights=e.dof) - np.average(l[col], weights=l.dof)
    return f


def h1_stat_regions(col):
    """Region fixed effects: within-region early-minus-late differences, weighted by
    the region's number of early maps (regions with fewer than 5 maps in either era skipped)."""
    def f(d):
        num = den = 0.0
        for _, g in d.groupby("region"):
            e = g[(g.year >= 1500) & (g.year < 1700)][col]
            l = g[(g.year >= 1850) & (g.year <= 1950)][col]
            if len(e) >= 5 and len(l) >= 5:
                num += len(e) * (e.mean() - l.mean())
                den += len(e)
        return num / den if den else np.nan
    return f


def tmid(d, col="A", weighted=False):
    return fit_sigmoid(d.year, d[col], d.dof if weighted else None)[2]


def h2_tmids(d, col="A", weighted=False):
    eu, rest = d[d.region == "Europe"], d[d.region != "Europe"]
    a = fit_sigmoid(eu.year, eu[col], eu.dof if weighted else None)[2]
    b = fit_sigmoid(rest.year, rest[col], rest.dof if weighted else None)[2]
    return (a, b, a - b)


def tmid_stat(col, weighted=False):
    return partial(tmid, col=col, weighted=weighted)


def h2_stat(col, weighted=False):
    return partial(h2_tmids, col=col, weighted=weighted)


def slope(d, col="A"):
    x = np.log10(d.dist_km.clip(lower=10))
    if len(d) < 10 or x.std() == 0:
        return np.nan
    return np.polyfit(x, d[col], 1)[0]


def binned_curve(d, col, label, step=25):
    out = []
    for lo in range(1450, 1950, step):
        g = d[(d.year >= lo) & (d.year < lo + step)]
        if len(g) < 8:
            continue
        p, ci = boot(g, lambda x: x[col].mean(), g.cluster, n=1000, seed=lo)
        out.append({"series": label, "stat": col, "bin_lo": lo, "mid": lo + step / 2, "n": len(g),
                    "mean": p, "lo": ci[0], "hi": ci[1], "ratio": np.exp(p)})
    return out


def main(placebo=False):
    """placebo=True shuffles publication years across maps before any test: a null run
    in which every time-dependent hypothesis should come out empty."""
    gcps = pickle.load(open(DERIVED / "gcps.pkl", "rb"))
    df = load()
    tag = "placebo_" if placebo else ""
    if placebo:
        rng = np.random.default_rng(99)
        ok = df.eligible.values
        df.loc[ok, "year"] = rng.permutation(df.loc[ok, "year"].values)
    base = df[df.eligible & (df.extent_deg >= 2)].copy()
    stats_, errs = run_fits(base, gcps, keep_errors=True)
    aff, _ = run_fits(base, gcps, force_order=1)
    m = base.merge(stats_, on="map_id", suffixes=("", "_fit"))
    m["region"] = [region(a, b) for a, b in zip(m.lon0_fit, m.lat0_fit)]
    m["cluster"] = m.maker.fillna(m.map_id).astype(str)
    m["dist_km"] = great_circle_km(m.place_lat, m.place_lon, m.lat0_fit, m.lon0_fit)
    if not placebo:
        m.to_csv(RESULTS / "map_stats.csv", index=False)
        errs.merge(m[["map_id", "year", "region"]], on="map_id").to_csv(DERIVED / "gcp_errors.csv", index=False)
    aff = base[["map_id", "year", "maker", "extent_deg"]].merge(aff, on="map_id")
    aff["cluster"] = aff.maker.fillna(aff.map_id).astype(str)
    aff["region"] = [region(a, b) for a, b in zip(aff.lon0, aff.lat0)]

    P = m[m.extent_deg >= 4]           # primary sample
    L = []
    L.append("LONGITUDE WIND: " + ("PLACEBO RUN, years shuffled across maps" if placebo else "results (PLAN.md tests in order)"))
    L.append(f"Allmaps maps in geometric sample: {len(df)}; with catalogue metadata: {int(df.has_meta.fillna(False).sum())}")
    L.append(f"eligible (dated 1450-1950, precise, not facsimile/view): {int(df.eligible.sum())}; "
             f"fitted with extent >= 2 deg: {len(m)}; primary (extent >= 4 deg): {len(P)}")
    L.append(f"primary sample by century: " + ", ".join(
        f"{c}s {((P.year >= c) & (P.year < c + 100)).sum()}" for c in range(1400, 2000, 100)))
    L.append(f"primary sample by region: " + ", ".join(f"{k} {v}" for k, v in P.region.value_counts().items()))
    L.append(f"makers (clusters) in primary sample: {P.cluster.nunique()}")
    L.append("")

    # --- known answers on real maps ---
    L.append("KNOWN-ANSWER CHECKS")
    k2 = P[(P.year >= 1900) & (P.S.abs() < 0.15)]
    p, ci = boot(k2, lambda d: d.A.mean(), k2.cluster)
    L.append(f"K2 modern floor (1900-1950, |S|<0.15, n={len(k2)}): mean A {fmt(p, ci)}  "
             f"-> {'PASS' if abs(p) <= 0.10 else 'FAIL'} (|mean| <= 0.10)")

    band = errs[errs.lat.between(30, 46)]
    west = set(band[band.lon < 0].map_id)
    east = set(band[band.lon > 28].map_id)
    med_ids = west & east
    med = m[m.map_id.isin(med_ids)]
    early_med, late_med = med[(med.year >= 1477) & (med.year <= 1650)], med[(med.year >= 1750) & (med.year <= 1950)]
    k3 = lambda d: d[(d.year >= 1477) & (d.year <= 1650)].S.mean() - d[(d.year >= 1750) & (d.year <= 1950)].S.mean()
    p, ci = boot(med, k3, med.cluster)
    L.append(f"K3 long Mediterranean (maps spanning it: {len(med)}; 1477-1650 n={len(early_med)}, 1750-1950 n={len(late_med)}): "
             f"mean S early {early_med.S.mean():+.3f}, late {late_med.S.mean():+.3f}; difference {fmt(p, ci)} "
             f"-> {'PASS' if p >= 0.15 else 'FAIL'} (>= 0.15)")
    pto = m[m.title.fillna("").str.contains("ptolem", case=False) & m.year.between(1477, 1600)]
    if len(pto) >= 10:
        p, ci = boot(pto, lambda d: d.A.mean(), pto.cluster)
        L.append(f"K4 Ptolemy maps 1477-1600 (n={len(pto)}, extent >= 2): mean A {fmt(p, ci)} "
                 f"-> {'PASS' if p > 0 else 'FAIL'} (> 0); mean S {pto.S.mean():+.3f}")
    else:
        L.append(f"K4 Ptolemy maps: only {len(pto)} found, not reported (needs 10)")
    L.append("")

    # --- hypotheses ---
    L.append("HYPOTHESES (primary sample, 95% cluster-bootstrap CIs over makers)")
    early = P[(P.year >= 1500) & (P.year < 1700)]
    late = P[(P.year >= 1850) & (P.year <= 1950)]
    pe, cie = boot(early, lambda d: d.A.mean(), early.cluster)
    pl, cil = boot(late, lambda d: d.A.mean(), late.cluster)
    pd_, cid = boot(P, h1_stat("A"), P.cluster)
    h1 = pe > 0 and cie[0] > 0 and cid[0] > 0
    L.append(f"H1 mean A 1500-1699 (n={len(early)}): {fmt(pe, cie)}  (EW/NS RMS ratio {np.exp(pe):.2f})")
    L.append(f"   mean A 1850-1950 (n={len(late)}): {fmt(pl, cil)}  (ratio {np.exp(pl):.2f})")
    L.append(f"   difference: {fmt(pd_, cid)} -> {'SUPPORTED' if h1 else 'NOT SUPPORTED'}")
    p_all = fit_sigmoid(P.year, P.A)
    pt, cit = boot(P, tmid_stat("A"), P.cluster, parallel=True)
    # PLAN: the anisotropy is lost, so the fitted step must be a decline (dA > 0) to count
    h1t = p_all[1] > 0 and 1680 <= pt <= 1830
    L.append(f"H1t sigmoid: A_late {p_all[0]:+.3f}, dA {p_all[1]:+.3f}, w {p_all[3]:.0f} yr; "
             f"t_mid {pt:.0f} [{cit[0]:.0f}, {cit[1]:.0f}] -> "
             + ("SUPPORTED (a decline, midpoint in 1680-1830)" if h1t else
                ("NOT SUPPORTED: the fitted change is a rise, not a decline" if p_all[1] <= 0
                 else "NOT SUPPORTED: midpoint outside 1680-1830")))
    n_eu = ((P.region == "Europe") & (P.year < 1700)).sum()
    n_rest = ((P.region != "Europe") & (P.year < 1700)).sum()
    if n_eu >= 40 and n_rest >= 40:
        (a, b, dlt), cis = boot(P, h2_stat("A"), P.cluster, parallel=True)
        eu_fit = fit_sigmoid(P[P.region == "Europe"].year, P[P.region == "Europe"].A)
        rest_fit = fit_sigmoid(P[P.region != "Europe"].year, P[P.region != "Europe"].A)
        declines = eu_fit[1] > 0 and rest_fit[1] > 0  # PLAN: both groups must lose anisotropy
        verdict = ("SUPPORTED" if cis[2][1] < 0 and declines else
                   "NOT SUPPORTED" + ("" if declines else
                                      f" (Europe's fitted change is a {'rise' if eu_fit[1] <= 0 else 'decline'}, "
                                      f"elsewhere a {'rise' if rest_fit[1] <= 0 else 'decline'})"))
        L.append(f"H2 t_mid Europe {a:.0f} [{cis[0][0]:.0f}, {cis[0][1]:.0f}] (pre-1700 n={n_eu}); "
                 f"rest {b:.0f} [{cis[1][0]:.0f}, {cis[1][1]:.0f}] (pre-1700 n={n_rest}); "
                 f"difference {dlt:+.0f} [{cis[2][0]:+.0f}, {cis[2][1]:+.0f}] -> {verdict}")
        for reg, g in P.groupby("region"):
            pr = fit_sigmoid(g.year, g.A)
            L.append(f"      {reg:22s} n={len(g):5d} pre-1700 {((g.year < 1700)).sum():4d} "
                     f"A_late {pr[0]:+.3f} dA {pr[1]:+.3f} t_mid {pr[2]:.0f} w {pr[3]:.0f}")
    else:
        L.append(f"H2 untestable: pre-1700 maps Europe {n_eu}, rest {n_rest} (needs 40 each)")
    hz = P[P.dist_km.notna()]
    pre, post = hz[hz.year < 1750], hz[hz.year > 1850]
    ps, cis_ = boot(pre, slope, pre.cluster)
    qs, ciq = boot(post, slope, post.cluster)
    ds, cids = boot(hz, lambda d: slope(d[d.year < 1750]) - slope(d[d.year > 1850]), hz.cluster)
    L.append(f"H3 slope of A on log10 distance from place of publication: before 1750 (n={len(pre)}) {fmt(ps, cis_)}; "
             f"after 1850 (n={len(post)}) {fmt(qs, ciq)}; difference {fmt(ds, cids)} -> "
             f"{'SUPPORTED' if cis_[0] > 0 and cids[0] > 0 else 'NOT SUPPORTED'}")
    rho, pv = spearmanr(P.year, P.relerr)
    L.append(f"S0 Spearman rho(year, RelErr) = {rho:+.3f} (p={pv:.1e}) -> {'PASS' if rho < 0 else 'FAIL'}")
    L.append("")

    # --- robustness ---
    L.append("ROBUSTNESS (H1 difference early-late; H2 t_mid Europe minus rest)")
    variants = [
        ("primary", P, "A", False, h1_stat("A")),
        ("R1 |S|<0.15", P[P.S.abs() < 0.15], "A", False, h1_stat("A")),
        ("R2 A_paper", P, "A_paper", False, h1_stat("A_paper")),
        ("R3 region fixed effects", P, "A", False, h1_stat_regions("A")),
        ("R4 extent >= 2", m, "A", False, h1_stat("A")),
        ("R4 extent >= 8", m[m.extent_deg >= 8], "A", False, h1_stat("A")),
        ("R5 weighted by dof", P, "A", True, h1_stat_weighted("A")),
        ("R6 affine fits only", aff[aff.extent_deg >= 4], "A", False, h1_stat("A")),
    ]
    rob = []
    for name, d, col, w, f in variants:
        p1, c1 = boot(d, f, d.cluster, seed=1)
        if name.startswith("R3"):
            p2, c2 = np.nan, (np.nan, np.nan)
        else:
            (a, b, p2), cs = boot(d, h2_stat(col, w), d.cluster, seed=2, parallel=True)
            c2 = cs[2]
        rob.append({"variant": name, "n": len(d), "h1_diff": p1, "h1_lo": c1[0], "h1_hi": c1[1],
                    "h2_diff": p2, "h2_lo": c2[0], "h2_hi": c2[1]})
        L.append(f"{name:26s} n={len(d):5d}  H1 {fmt(p1, c1)}   H2 {p2:+.0f} [{c2[0]:+.0f}, {c2[1]:+.0f}]")
    pd.DataFrame(rob).to_csv(RESULTS / f"{tag}robustness.csv", index=False)
    r1 = next(r for r in rob if r["variant"].startswith("R1"))
    r3 = next(r for r in rob if r["variant"].startswith("R3"))
    robust = h1 and r1["h1_lo"] > 0 and r3["h1_lo"] > 0
    L.append(f"H1 robust (primary + R1 + R3 all exclude 0): {'YES' if robust else 'NO'}")

    curve = binned_curve(P, "A", "all") + binned_curve(P[P.region == "Europe"], "A", "Europe") \
        + binned_curve(P[P.region != "Europe"], "A", "rest") + binned_curve(P, "relerr", "all") \
        + binned_curve(med, "S", "mediterranean")
    if not placebo:
        pd.DataFrame(curve).to_csv(RESULTS / "curve.csv", index=False)
        # the H3 picture in full: anisotropy by era and distance from the place of publication
        rows = []
        for lo, hi in [(1450, 1750), (1750, 1850), (1850, 1951)]:
            for dlo, dhi in [(0, 500), (500, 2000), (2000, 6000), (6000, 25000)]:
                g = hz[(hz.year >= lo) & (hz.year < hi) & (hz.dist_km >= dlo) & (hz.dist_km < dhi)]
                if len(g) < 8:
                    continue
                p, ci = boot(g, partial(mean_of, col="A"), g.cluster, n=1000, seed=dlo + lo)
                rows.append({"era_lo": lo, "era_hi": hi, "dist_lo": dlo, "dist_hi": dhi, "n": len(g),
                             "A": p, "lo": ci[0], "hi": ci[1], "ratio": np.exp(p)})
        pd.DataFrame(rows).to_csv(RESULTS / "distance_table.csv", index=False)
    text = "\n".join(L)
    (RESULTS / f"{tag}summary.txt").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    import sys
    main(placebo="--placebo" in sys.argv)
