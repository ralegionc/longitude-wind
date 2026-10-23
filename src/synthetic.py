"""Known-answer checks K1a-c (PLAN.md): rebuild every map as a synthetic map whose
error structure is known, run the real pipeline, and see if it reads the truth back.

K1a  conformal projection + isotropic paper noise      -> mean A      = 0
K1b  non-conformal projection + isotropic paper noise  -> mean A_paper = 0
K1c  conformal projection + ground error 2x east-west  -> mean A      = ln 2
"""
import pickle

import numpy as np
import pandas as pd
from pyproj import Proj

from common import DERIVED, RESULTS
from fit import R, fit_map

CONFORMAL = ["merc", "lcc", "stere", "tmerc"]
NONCONFORMAL = ["platecarree", "equirect", "eqdc", "aea", "bonne", "sinu", "cass"]


def proj_for(name, lon0, lat0, lat_span):
    q = max(lat_span / 4, 1.0)
    lat1, lat2 = lat0 - q, lat0 + q
    if lat1 * lat2 <= 0:  # standard parallels must not straddle the equator symmetrically
        lat1, lat2 = (0.5, lat2) if lat0 >= 0 else (lat1, -0.5)
    base = f"+lon_0={lon0} +R={R * 1000}"
    s = {
        "merc": f"+proj=merc {base}",
        "lcc": f"+proj=lcc +lat_1={lat1} +lat_2={lat2} +lat_0={lat0} {base}",
        "stere": f"+proj=stere +lat_0={lat0} {base}",
        "tmerc": f"+proj=tmerc +lat_0={lat0} {base}",
        "platecarree": f"+proj=eqc +lat_ts=0 {base}",
        "equirect": f"+proj=eqc +lat_ts={lat0} {base}",
        "eqdc": f"+proj=eqdc +lat_1={lat1} +lat_2={lat2} {base}",
        "aea": f"+proj=aea +lat_1={lat1} +lat_2={lat2} {base}",
        "bonne": f"+proj=bonne +lat_1={lat0 if abs(lat0) > 1 else 1.0} {base}",
        "sinu": f"+proj=sinu {base}",
        "cass": f"+proj=cass +lat_0={lat0} {base}",
    }[name]
    return Proj(s)


def to_paper(x, y, rng):
    th = rng.uniform(-0.26, 0.26)
    c, s = np.cos(th), np.sin(th)
    px, py = c * x - s * y, s * x + c * y
    scale = 6000.0 / max(np.ptp(px), np.ptp(py))
    return np.column_stack([3000 + scale * px, 3000 - scale * py])


def one(lon, lat, lon0, lat0, lat_span, family, rng, mode):
    if mode == "ground2x":
        sg = 0.01 * R * np.radians(max(np.ptp(lat), np.ptp(lon) * np.cos(np.radians(lat0))))
        e, n = rng.normal(0, 2 * sg, len(lon)), rng.normal(0, sg, len(lon))
        lon_d = lon + np.degrees(e / (R * np.cos(np.radians(lat))))
        lat_d = lat + np.degrees(n / R)
    else:
        lon_d, lat_d = lon, lat
    x, y = proj_for(family, lon0, lat0, lat_span)(lon_d, lat_d)
    x, y = np.asarray(x), np.asarray(y)
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        return None  # this projection cannot draw this map (e.g. Mercator near a pole)
    pix = to_paper(x, y, rng)
    if mode == "paper":
        pix = pix + rng.normal(0, 0.01 * 6000.0, pix.shape)
    return fit_map(pix, np.column_stack([lon, lat]))


def main(seed=11):
    rng = np.random.default_rng(seed)
    from analyze import load
    geom = load()
    geom = geom[geom.eligible & (geom.extent_deg >= 2)]  # PLAN: every eligible map
    gcps = pickle.load(open(DERIVED / "gcps.pkl", "rb"))
    rows = []
    for _, g in geom.iterrows():
        _, geo = gcps[g.map_id]
        lon = g.lon0 + ((geo[:, 0] - g.lon0 + 180) % 360 - 180)
        lat = geo[:, 1]
        for check, fams, mode in [("K1a", CONFORMAL, "paper"), ("K1b", NONCONFORMAL, "paper"),
                                  ("K1c", CONFORMAL, "ground2x")]:
            fam = fams[rng.integers(len(fams))]
            r = one(lon, lat, g.lon0, g.lat0, g.lat_span, fam, rng, mode)
            if r:
                rows.append({"check": check, "family": fam, "map_id": g.map_id, "extent_deg": g.extent_deg,
                             "order": r["order"], "A": r["A"], "S": r["S"], "A_paper": r["A_paper"]})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "synthetic_checks.csv", index=False)
    lines = []
    for check, stat, target, tol in [("K1a", "A", 0.0, 0.03), ("K1b", "A_paper", 0.0, 0.03), ("K1c", "A", np.log(2), 0.10)]:
        for lo, hi in [(4, 999), (2, 4), (4, 8), (8, 20), (20, 999)]:
            d = df[(df.check == check) & (df.extent_deg >= lo) & (df.extent_deg < hi)]
            m = d[stat].mean()
            se = d[stat].std() / np.sqrt(len(d))
            tag = "primary" if (lo, hi) == (4, 999) else f"extent {lo}-{hi}"
            verdict = ("PASS" if abs(m - target) <= tol else "FAIL") if (lo, hi) == (4, 999) else ""
            lines.append(f"{check} {stat:8s} {tag:14s} n={len(d):5d} mean={m:+.3f} (se {se:.3f}) "
                         f"target {target:+.3f}+-{tol} {verdict}")
        for fam, d in df[(df.check == check) & (df.extent_deg >= 4)].groupby("family"):
            lines.append(f"    {fam:12s} n={len(d):5d} mean {stat}={d[stat].mean():+.3f}  mean S={d.S.mean():+.3f}")
    text = "\n".join(lines)
    (RESULTS / "synthetic_checks.txt").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
