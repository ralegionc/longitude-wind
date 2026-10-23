"""Data for the interactive page in docs/: land outline, error streaks per era, one
point per map (with a link to open it in the Allmaps viewer), and the curves."""
import json
from urllib.parse import quote

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import mapping

from common import DERIVED, DOCS, RAW, RESULTS
from figures import ERAS, segments

OUT = DOCS / "data"
WEB_EXAG = 10          # streaks drawn ten times the error, as the page says
WEB_SAMPLE = 30000     # places per period, to keep each file small


def rnd(x, k=3):
    return float(round(x, k))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    land = gpd.read_file(f"zip://{RAW / 'ne' / 'ne_50m_land.zip'}")
    land["geometry"] = land.geometry.simplify(0.08)
    feats = [{"type": "Feature", "properties": {}, "geometry": mapping(g)} for g in land.geometry if not g.is_empty]
    (OUT / "land.json").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                              separators=(",", ":")), encoding="utf-8")

    stats = pd.read_csv(RESULTS / "map_stats.csv")
    stats = stats[stats.extent_deg >= 4]
    err = pd.read_csv(DERIVED / "gcp_errors.csv")
    err = err[err.map_id.isin(stats.map_id)]
    for i, (lo, hi, label) in enumerate(ERAS):
        e = err[(err.year >= lo) & (err.year < hi)]
        e = e.sample(min(len(e), WEB_SAMPLE), random_state=0)
        seg, orient, tot = segments(e, exag=WEB_EXAG)
        feats = [{"type": "Feature", "properties": {"o": rnd(o, 2), "k": int(t)},
                  "geometry": {"type": "LineString", "coordinates": [[rnd(a), rnd(b)], [rnd(c), rnd(d)]]}}
                 for (a, b), (c, d), o, t in zip(seg[:, 0], seg[:, 1], orient, tot)]
        (OUT / f"streaks_{i}.json").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                                         separators=(",", ":")), encoding="utf-8")

    stretch = {}
    if (RESULTS / "r2_stretch.csv").exists():  # round 2: longitude stretch for maps >= 20 degrees
        stretch = pd.read_csv(RESULTS / "r2_stretch.csv").set_index("map_id").lnk.to_dict()
    pts = []
    for s in stats.itertuples():
        title = str(s.title) if isinstance(s.title, str) else ""
        pts.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [rnd(s.lon0_fit, 2), rnd(s.lat0_fit, 2)]},
                    "properties": {"t": title[:140], "y": int(s.year), "r": rnd(np.exp(s.A), 2),
                                   "e": int(round(s.rms_km)), "n": int(s.n),
                                   **({"k": rnd(np.exp(stretch[s.map_id]), 2)} if s.map_id in stretch else {}),
                                   "u": "https://viewer.allmaps.org/?url=" + quote(s.map_id, safe="")}})
    (OUT / "maps.json").write_text(json.dumps({"type": "FeatureCollection", "features": pts},
                                              separators=(",", ":"), ensure_ascii=False), encoding="utf-8")
    curve = pd.read_csv(RESULTS / "curve.csv")
    c = curve[curve.stat == "A"][["series", "mid", "n", "mean", "lo", "hi"]]
    (OUT / "curves.json").write_text(c.to_json(orient="records"), encoding="utf-8")

    eras = []
    for lo, hi, label in ERAS:
        s = stats[(stats.year >= lo) & (stats.year < hi)]
        e = err[(err.year >= lo) & (err.year < hi)]
        eras.append({"label": label, "n_maps": int(len(s)), "n_places": int(len(e)),
                     "ratio": rnd(np.exp(s.A.mean()), 3), "median_km": rnd(float(np.median(np.hypot(e.e_ew, e.e_ns))), 1)})
    table = pd.read_csv(RESULTS / "distance_table.csv")
    early = table[table.era_lo == 1450]
    far = early[early.dist_lo == 6000].ratio.iloc[0]
    near = early[early.dist_lo == 0].ratio.iloc[0]
    fr = pd.read_csv(RESULTS / "frames.csv")
    side = fr[fr.tilt.abs().between(60, 120)]
    st2 = pd.read_csv(RESULTS / "r2_stretch.csv")
    pre, post = st2[st2.year < 1750].abs_lnk.mean(), st2[st2.year >= 1750].abs_lnk.mean()
    note = (f"Before 1750, maps drew longitudes about {100 * (np.exp(pre) - 1):.0f}% too wide or too narrow; after "
            f"1750, about {100 * (np.exp(post) - 1):.0f}%. Places over 6,000 km from where a map was printed were also "
            f"misplaced {far:.2f} times as far east-west as north-south, against {near:.2f} within 500 km. After 1850 a "
            f"small excess remains on every map, in two independent volunteer tools; whether it follows the page or "
            f"the compass is still open.")
    del side
    (OUT / "summary.json").write_text(json.dumps({"eras": eras, "note": note}, indent=1), encoding="utf-8")
    sizes = {p.name: p.stat().st_size for p in OUT.iterdir()}
    print({k: f"{v / 1e6:.2f} MB" for k, v in sizes.items()})


if __name__ == "__main__":
    main()
