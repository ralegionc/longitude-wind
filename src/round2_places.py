"""PLAN2 part D: the same places through time.

D1 event study at five places with a documented first astronomical longitude.
D2 error trajectories for well-mapped ports and capitals (figure in figures_r2.py).
"""
import numpy as np
import pandas as pd

from analyze import great_circle_km
from common import DERIVED, RESULTS
from stats import boot

RADIUS_KM = 15.0
EVENTS = {  # place: (lat, lon, year of first astronomical longitude, by whom)
    "Cayenne": (4.937, -52.326, 1672, "Richer"),
    "Alexandria": (31.200, 29.919, 1694, "Chazelles"),
    "Lima": (-12.046, -77.043, 1709, "Feuillee"),
    "Cape Town": (-33.925, 18.424, 1751, "Lacaille"),
    "Tahiti": (-17.495, -149.495, 1769, "Green and Cook"),
}
PORTS = {
    "London": (51.507, -0.128), "Paris": (48.857, 2.352), "Amsterdam": (52.368, 4.904),
    "Lisbon": (38.722, -9.139), "Rome": (41.903, 12.496), "Istanbul": (41.008, 28.978),
    "Cairo": (30.044, 31.236), "Alexandria": (31.200, 29.919), "Cape Town": (-33.925, 18.424),
    "Mumbai": (18.939, 72.835), "Jakarta": (-6.135, 106.813), "Guangzhou": (23.129, 113.264),
    "Beijing": (39.904, 116.407), "Nagasaki": (32.750, 129.878), "Manila": (14.600, 120.984),
    "Lima": (-12.046, -77.043), "Mexico City": (19.433, -99.133), "Havana": (23.114, -82.367),
    "Boston": (42.360, -71.059), "Quebec": (46.814, -71.208), "Rio de Janeiro": (-22.907, -43.173),
    "Buenos Aires": (-34.604, -58.382), "Cayenne": (4.937, -52.326), "Tahiti": (-17.495, -149.495),
}


def place_hits(err, lat, lon):
    """One control point per map: the nearest within RADIUS_KM of the place."""
    near = err[(err.lat - lat).abs() < 0.5]
    d = great_circle_km(near.lat.values, near.lon.values, lat, lon)
    near = near.assign(d=d)[d <= RADIUS_KM]
    return near.sort_values("d").groupby("map_id").head(1)


def event_stat(d):
    """Mean over places of [log(post/pre median |EW|) - log(post/pre median |NS|)]."""
    vals = []
    for place, g in d.groupby("place"):
        pre, post = g[g.side == "pre"], g[g.side == "post"]
        if len(pre) >= 5 and len(post) >= 5:
            ew = np.log(post.e_ew.abs().median() / pre.e_ew.abs().median())
            ns = np.log(post.e_ns.abs().median() / pre.e_ns.abs().median())
            vals.append(ew - ns)
    return np.mean(vals) if vals else np.nan


def main():
    err = pd.read_csv(DERIVED / "gcp_errors.csv")
    st = pd.read_csv(RESULTS / "map_stats.csv")[["map_id", "cluster"]]
    err = err.merge(st, on="map_id")
    L = ["ROUND 2, PART D: the same places through time", "",
         "D1 event study: 60 years before vs after the first astronomical longitude (one point per map)"]
    rows = []
    for place, (lat, lon, yr, who) in EVENTS.items():
        h = place_hits(err, lat, lon)
        h = h[(h.year >= yr - 60) & (h.year <= yr + 60)].copy()
        h["side"] = np.where(h.year < yr, "pre", "post")
        h["place"] = place
        pre, post = h[h.side == "pre"], h[h.side == "post"]
        L.append(f"   {place:10s} {yr} ({who}): maps before {len(pre):3d}, after {len(post):3d}; median |EW| "
                 f"{pre.e_ew.abs().median():6.1f} -> {post.e_ew.abs().median():6.1f} km, |NS| "
                 f"{pre.e_ns.abs().median():6.1f} -> {post.e_ns.abs().median():6.1f} km")
        rows.append(h)
    d = pd.concat(rows, ignore_index=True)
    used = [p for p, g in d.groupby("place") if (g.side == "pre").sum() >= 5 and (g.side == "post").sum() >= 5]
    p, ci = boot(d, event_stat, d.cluster)
    L.append(f"   places with at least 5 maps each side: {', '.join(used) if used else 'none'}")
    L.append(f"   D1 pooled difference in log ratios (EW minus NS): {p:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}] -> "
             f"{'PASS' if ci[1] < 0 else 'FAIL'} (prediction: negative)")
    L.append("")
    traj = []
    for place, (lat, lon) in PORTS.items():
        h = place_hits(err, lat, lon).assign(place=place)
        traj.append(h)
        L.append(f"   D2 {place:15s} maps {len(h):4d}  years {h.year.min():.0f}-{h.year.max():.0f}" if len(h) else
                 f"   D2 {place:15s} maps    0")
    pd.concat(traj, ignore_index=True).to_csv(RESULTS / "r2_place_trajectories.csv", index=False)
    text = "\n".join(L)
    (RESULTS / "r2_places.txt").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
