"""PLAN2 part C: longitude stretch beyond any honest projection, for maps >= 20 degrees.

For each of seven projections with conventional parameters, longitudes are scaled about
the map centre by k before projecting, and a similarity transform (scale, rotation,
shift) maps the result to the page. The projection with the smallest residual wins and
ln k is the map's longitude stretch.
"""
import pickle
import warnings
from functools import partial

import numpy as np
import pandas as pd
from pyproj import Proj
from scipy.optimize import minimize_scalar

from analyze import load
from common import DERIVED, RESULTS
from fit import R, fit_map
from stats import boot

FAMILIES = ["equirect", "merc", "lcc", "eqdc", "bonne", "sinu", "stere"]
warnings.filterwarnings("ignore")


def projector(family, lon0, lat0, lat_min, lat_max):
    span = lat_max - lat_min
    p1, p2 = lat_min + span / 6, lat_max - span / 6
    if abs(p1 + p2) < 1e-6:  # conic with standard parallels symmetric about the equator is singular
        p2 += 1.0
    base = f"+lon_0=0 +R={R * 1000} +no_defs"
    if family == "equirect":
        c = np.cos(np.radians(lat0))
        return lambda lam, phi: (R * np.radians(lam) * c, R * np.radians(phi))
    s = {"merc": f"+proj=merc {base}",
         "lcc": f"+proj=lcc +lat_1={p1} +lat_2={p2} +lat_0={lat0} {base}",
         "eqdc": f"+proj=eqdc +lat_1={p1} +lat_2={p2} +lat_0={lat0} {base}",
         "bonne": f"+proj=bonne +lat_1={lat0 if abs(lat0) > 1 else 1.0} {base}",
         "sinu": f"+proj=sinu {base}",
         "stere": f"+proj=stere +lat_0={lat0} {base}"}[family]
    P = Proj(s)
    return lambda lam, phi: P(lam, phi)


def similarity_rss(x, y, u, v):
    """Least-squares similarity (with a mirror for the downward pixel axis); returns RSS."""
    vu = -v
    A = np.zeros((2 * len(x), 4))
    A[0::2] = np.column_stack([x, -y, np.ones_like(x), np.zeros_like(x)])
    A[1::2] = np.column_stack([y, x, np.zeros_like(x), np.ones_like(x)])
    b = np.empty(2 * len(x))
    b[0::2], b[1::2] = u, vu
    coef, *_ = np.linalg.lstsq(A, b, rcond=None)
    return float(np.sum((A @ coef - b) ** 2))


def fit_stretch(lon_rel, lat, u, v, lat0):
    """lon_rel: longitudes relative to the map centre. Returns best family, ln k, rss by family."""
    lat_min, lat_max = lat.min(), lat.max()
    out = {}
    for fam in FAMILIES:
        P = projector(fam, 0.0, lat0, lat_min, lat_max)

        def rss(lnk):
            x, y = P(lon_rel * np.exp(lnk), lat)
            x, y = np.asarray(x), np.asarray(y)
            if not (np.isfinite(x).all() and np.isfinite(y).all()):
                return np.inf
            return similarity_rss(x, y, u, v)
        r = minimize_scalar(rss, bounds=(-0.7, 0.7), method="bounded", options={"xatol": 1e-4})
        out[fam] = (r.fun, r.x)
    best = min(out, key=lambda f: out[f][0])
    return best, out[best][1], {f: out[f][0] for f in out}


def map_inputs(r):
    lon0 = r["lon0"]
    lon_rel = (r["gcp_lon"] - lon0 + 180) % 360 - 180
    return lon_rel, r["gcp_lat"], r["gcp_u"], r["gcp_v"], r["lat0"]


def c0(df, gcps, rng, per_cell=150):
    rows = []
    sample = df.sample(min(len(df), per_cell * 3), random_state=1)
    for fam in FAMILIES:
        for k in (0.8, 1.0, 1.25):
            for g in sample.sample(per_cell, random_state=hash((fam, k)) % 2**31).itertuples():
                pix, geo = gcps[g.map_id]
                r = fit_map(pix, geo)
                if r is None:
                    continue
                lon_rel, lat, _, _, lat0 = map_inputs(r)
                P = projector(fam, 0.0, lat0, lat.min(), lat.max())
                x, y = P(lon_rel * k, lat)
                x, y = np.asarray(x), np.asarray(y)
                if not (np.isfinite(x).all() and np.isfinite(y).all()):
                    continue
                th = rng.uniform(-0.2, 0.2)
                X, Y = np.cos(th) * x - np.sin(th) * y, np.sin(th) * x + np.cos(th) * y
                sc = 6000 / max(np.ptp(X), np.ptp(Y))
                u, v = 3000 + sc * X, 3000 - sc * Y
                u, v = u + rng.normal(0, 60, len(u)), v + rng.normal(0, 60, len(v))
                best, lnk, _ = fit_stretch(lon_rel, lat, u, v, lat0)
                rows.append({"family": fam, "k": k, "best": best, "lnk": lnk, "err": lnk - np.log(k)})
    return pd.DataFrame(rows)


def mean_abs(d):
    return d.abs_lnk.mean()


def diff_abs(d):
    e = d[(d.year >= 1500) & (d.year < 1700)].abs_lnk.mean()
    l_ = d[d.year >= 1850].abs_lnk.mean()
    return e - l_


def first_below(d, thr=0.05, step=25):
    bins = []
    for lo in range(1450, 1950, step):
        g = d[(d.year >= lo) & (d.year < lo + step)]
        if len(g) >= 8:
            bins.append((lo, g.abs_lnk.mean()))
    for i, (lo, m) in enumerate(bins):
        if all(mm < thr for _, mm in bins[i:]):
            return lo, bins
    return None, bins


def main():
    gcps = pickle.load(open(DERIVED / "gcps.pkl", "rb"))
    st = pd.read_csv(RESULTS / "map_stats.csv")
    st = st[st.extent_deg >= 20][["map_id", "region", "cluster", "title"]]
    df = load().drop(columns=["title"]).merge(st, on="map_id")
    df = df[df.eligible]
    rng = np.random.default_rng(31)
    L = ["ROUND 2, PART C: longitude stretch beyond the projection (maps >= 20 degrees)", ""]

    syn = c0(df, gcps, rng)
    syn.to_csv(RESULTS / "r2_stretch_synthetic.csv", index=False)
    L.append("C0 synthetic known answers (tolerance: mean ln k error within 0.03 per family; right family >= 70%):")
    ok = True
    for fam, g in syn.groupby("family"):
        acc = (g.best == fam).mean()
        good = abs(g.err.mean()) <= 0.03
        ok &= good
        L.append(f"   {fam:9s} n={len(g):4d} mean error {g.err.mean():+.4f} (median |err| {g.err.abs().median():.4f}) "
                 f"right family {acc:.0%}")
    acc_all = (syn.best == syn.family).mean()
    L.append(f"   overall right family {acc_all:.0%}; C0 -> {'PASS' if ok and acc_all >= 0.7 else 'FAIL'}")
    L.append("")

    rows = []
    for g in df.itertuples():
        pix, geo = gcps[g.map_id]
        r = fit_map(pix, geo)
        if r is None:
            continue
        lon_rel, lat, u, v, lat0 = map_inputs(r)
        best, lnk, rss = fit_stretch(lon_rel, lat, u, v, lat0)
        plane = -np.log(np.cos(np.radians(lat0)))
        rows.append({"map_id": g.map_id, "year": g.year, "region": g.region, "cluster": g.cluster,
                     "title": g.title, "lat0": lat0, "best": best, "lnk": lnk, "abs_lnk": abs(lnk),
                     "plane_chart": bool(plane > 0.1 and abs(lnk - plane) < 0.05)})
    d = pd.DataFrame(rows)
    d.to_csv(RESULTS / "r2_stretch.csv", index=False)
    L.append(f"maps fitted: {len(d)}; best family counts: " + ", ".join(f"{k} {v}" for k, v in d.best.value_counts().items()))
    L.append(f"plane-chart-like maps (ln k within 0.05 of -ln cos centre latitude): {d.plane_chart.sum()}")
    for lo, hi in [(1450, 1650), (1650, 1750), (1750, 1850), (1850, 1951)]:
        g = d[(d.year >= lo) & (d.year < hi)]
        p, ci = boot(g, mean_abs, g.cluster, seed=lo)
        L.append(f"   {lo}-{hi - 1}: maps {len(g):4d} mean |ln k| {p:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]  mean ln k {g.lnk.mean():+.3f}")
    L.append("")
    p, ci = boot(d, diff_abs, d.cluster)
    L.append(f"C1 mean |ln k| 1500-1699 minus 1850-1950: {p:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}] -> "
             f"{'SUPPORTED' if ci[0] > 0 else 'NOT SUPPORTED'}")
    eu, eub = first_below(d[d.region == "Europe"])
    rest, restb = first_below(d[d.region != "Europe"])
    L.append(f"C2 first 25-year bin after which mean |ln k| stays below 0.05: Europe {eu}, elsewhere {rest} -> "
             + ("SUPPORTED" if eu is not None and (rest is None or eu < rest) else "NOT SUPPORTED"))
    L.append("   Europe bins: " + " ".join(f"{lo}:{m:.3f}" for lo, m in eub))
    L.append("   elsewhere bins: " + " ".join(f"{lo}:{m:.3f}" for lo, m in restb))
    err = pd.read_csv(DERIVED / "gcp_errors.csv")
    band = err[err.lat.between(30, 46)]
    med = set(band[band.lon < 0].map_id) & set(band[band.lon > 28].map_id)
    m = d[d.map_id.isin(med) & (d.year < 1650)]
    if len(m) >= 5:
        p, ci = boot(m, lambda z: z.lnk.mean(), m.cluster)
        L.append(f"C3 Mediterranean maps before 1650 (n={len(m)}): mean ln k {p:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}] -> "
                 f"{'PASS' if ci[0] > 0.1 else 'FAIL'} (needs > 0.1)")
    else:
        L.append(f"C3 Mediterranean maps before 1650: only {len(m)}")
    text = "\n".join(L)
    (RESULTS / "r2_stretch.txt").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
