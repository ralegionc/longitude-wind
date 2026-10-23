"""Round-2 figures: page vs compass by period, stretch over time, the same places
through time, and error portraits drawn on the scans themselves."""
import io
import json
import urllib.request

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch

from common import DERIVED, FIGURES, RAW, RESULTS, USER_AGENT
from figures import COAST, EUROPE_C, EW, INK, MUTED, NS, REST_C, TEXT, ratio_axis


def style(ax):
    for s_ in ax.spines.values():
        s_.set_color(COAST)


def ratio(a):
    """Log-variance contrast a to an RMS ratio: exp(atanh a)."""
    return np.exp(np.arctanh(np.clip(a, -0.95, 0.95)))


def decomposition():
    d = pd.read_csv(RESULTS / "r2_decomposition.csv")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), dpi=150, gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    for grp, col, mk in [("Europe", EUROPE_C, "o"), ("elsewhere", REST_C, "s")]:
        g = d[d.group == grp].sort_values("lo")
        x = (g.lo + g.hi.clip(upper=1950)) / 2
        off = -4 if grp == "Europe" else 4
        ax.errorbar(x + off, ratio(g.a_c), yerr=[ratio(g.a_c) - ratio(g.a_c_lo), ratio(g.a_c_hi) - ratio(g.a_c)],
                    fmt="-" + mk, color=col, lw=2.2, ms=6, capsize=0, label=f"compass term, maps of {grp}")
    g = d[d.group.isna()].sort_values("lo")
    x = (g.lo + g.hi.clip(upper=1950)) / 2
    ax.errorbar(x, ratio(g.a_p), yerr=[ratio(g.a_p) - ratio(g.a_p_lo), ratio(g.a_p_hi) - ratio(g.a_p)],
                fmt="--D", color=NS, lw=1.6, ms=5, capsize=0, label="page term, all maps")
    ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ratio_axis(ax, (0.9, 1.0, 1.25, 1.5))
    ax.set_xlim(1450, 1950)
    ax.set_xlabel("year the map was published")
    ax.set_ylabel("east-west / north-south (or across / down the page)")
    ax.legend(frameon=False, fontsize=9.5, loc="upper right")
    ax.set_title("The compass term peaks in 1750-1849, then falls; the page term grows late", loc="left",
                 fontsize=12.5, weight="bold")
    style(ax)
    ax = axes[1]
    g = d[d.group.fillna("").str.startswith("dist")].copy()
    g["mid"] = g.group.str.extract(r"dist (\d+)-(\d+)").astype(float).apply(lambda r: np.sqrt(max(r[0], 100) * min(r[1], 15000)), axis=1)
    g = g.sort_values("mid")
    for col_, name, c in [("a_c", "compass term", EW), ("a_p", "page term", NS)]:
        ax.errorbar(g.mid, ratio(g[col_]), yerr=[ratio(g[col_]) - ratio(g[col_ + "_lo"]), ratio(g[col_ + "_hi"]) - ratio(g[col_])],
                    fmt="-o", color=c, lw=2, capsize=0, label=name)
    ax.set_xscale("log")
    ax.set_xticks([200, 1000, 4000, 10000])
    ax.set_xticklabels(["under 500 km", "500-2,000", "2,000-6,000", "over 6,000"])
    ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ratio_axis(ax, (0.8, 1.0, 1.25, 1.5, 2.0))
    ax.set_xlabel("distance from where the map was printed (maps before 1750)")
    ax.legend(frameon=False, fontsize=9.5, loc="upper left")
    ax.set_title("Before 1750, only the farthest places carry an excess", loc="left", fontsize=12.5, weight="bold")
    style(ax)
    fig.text(0.01, 0.01, "Round 2 (PLAN2.md part B). Bars are 95% cluster-bootstrap intervals over map makers.",
             fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(FIGURES / "r2_decomposition.png")
    plt.close(fig)


def stretch():
    d = pd.read_csv(RESULTS / "r2_stretch.csv")
    fig, ax = plt.subplots(figsize=(12, 5.6), dpi=150)
    for grp, col, sel in [("Europe", EUROPE_C, d.region == "Europe"), ("elsewhere", REST_C, d.region != "Europe")]:
        g = d[sel]
        ax.scatter(g.year, np.exp(g.lnk), s=6, color=col, alpha=0.25, lw=0)
        rows = []
        for lo in range(1450, 1950, 25):
            gg = g[(g.year >= lo) & (g.year < lo + 25)]
            if len(gg) >= 8:
                rows.append((lo + 12.5, np.exp(gg.lnk.median())))
        if rows:
            x, y = zip(*rows)
            ax.plot(x, y, "-o", color=col, lw=2.2, ms=4, label=f"maps of {grp}: 25-year median")
    ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ax.text(1455, 1.02, "longitudes drawn true", color=MUTED, fontsize=9)
    ratio_axis(ax, (0.7, 0.8, 1.0, 1.25, 1.5))
    ax.set_ylim(0.65, 1.6)
    ax.set_xlim(1450, 1950)
    ax.set_xlabel("year the map was published")
    ax.set_ylabel("longitude stretch k (beyond the best projection)")
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    ax.set_title("Maps drew longitudes too wide until about 1720, then true", loc="left", fontsize=14, weight="bold")
    style(ax)
    fig.text(0.01, 0.01, "Round 2 (PLAN2.md part C): maps spanning at least 20 degrees, each fitted with seven projections.",
             fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(FIGURES / "r2_stretch.png")
    plt.close(fig)


def places(n=16):
    t = pd.read_csv(RESULTS / "r2_place_trajectories.csv")
    counts = t.groupby("place").size()
    keep = counts[counts >= 25].sort_values(ascending=False).index[:n]
    cols = 4
    rows = int(np.ceil(len(keep) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(16, 3.2 * rows), dpi=150, sharex=True, sharey=True)
    for ax, place in zip(axes.flat, keep):
        g = t[t.place == place].sort_values("year")
        ax.scatter(g.year, g.e_ew.abs().clip(lower=0.3), s=7, color=EW, alpha=0.5, lw=0)
        ax.scatter(g.year, g.e_ns.abs().clip(lower=0.3), s=7, color=NS, alpha=0.5, lw=0)
        for c, col in [(EW, "e_ew"), (NS, "e_ns")]:
            s = g.set_index("year")[col].abs().rolling(15, center=True, min_periods=6).median()
            ax.plot(s.index, s.values, color=c, lw=2)
        ax.set_yscale("log")
        ax.set_ylim(0.3, 1000)
        ax.set_yticks([1, 10, 100, 1000])
        ax.set_yticklabels(["1", "10", "100", "1,000"])
        ax.minorticks_off()
        ax.set_title(f"{place} ({len(g)} maps)", fontsize=11, loc="left", color=TEXT)
        style(ax)
    for ax in axes.flat[len(keep):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("error (km)")
    fig.suptitle("The same places through time: east-west error (cyan) and north-south error (amber)",
                 x=0.01, ha="left", fontsize=15, weight="bold", color=TEXT)
    fig.text(0.01, 0.005, "One control point per map within 15 km of the place; lines are rolling medians of 15 maps.",
             fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(FIGURES / "r2_places.png")
    plt.close(fig)


def iiif_image(image_id, width=1600):
    url = f"{image_id.rstrip('/')}/full/{width},/0/default.jpg"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        return plt.imread(io.BytesIO(r.read()), format="jpg")


def portrait(ax, map_id, gcps, geom, meta, exag=1.0):
    """The scan, with an arrow at each place from where its true position falls on this
    map (by the map's own best fit) to where the map actually drew it."""
    from fit import fit_map
    from paper import forward_fit
    pix, geo = gcps[map_id]
    r = fit_map(pix, geo)
    D, coef = forward_fit(r)
    pred = D(r["gcp_lon"], r["gcp_lat"]) @ coef
    drawn = np.column_stack([r["gcp_u"], r["gcp_v"]])
    g = geom.loc[map_id]
    img = iiif_image(g.image)
    s = img.shape[1] / g.width
    ax.imshow(img, alpha=0.9)
    for (u0, v0), (u1, v1) in zip(pred * s, (pred + exag * (drawn - pred)) * s):
        a = FancyArrowPatch((u0, v0), (u1, v1), arrowstyle="-|>", mutation_scale=8,
                            color=EW if abs(u1 - u0) > abs(v1 - v0) else NS, lw=1.6)
        a.set_path_effects([pe.Stroke(linewidth=3.2, foreground="#0b0e13"), pe.Normal()])
        ax.add_patch(a)
    ax.set_axis_off()
    m = meta.loc[map_id]
    title = str(m.title)[:70]
    ax.set_title(f"{int(m.year)}: {title}", fontsize=10, loc="left", color=TEXT)


def portraits(choices, exag=3.0):
    import pickle
    gcps = pickle.load(open(DERIVED / "gcps.pkl", "rb"))
    geom = pd.read_csv(DERIVED / "maps_geom.csv").set_index("map_id")
    meta = pd.read_csv(RESULTS / "map_stats.csv").set_index("map_id")
    fig, axes = plt.subplots(1, len(choices), figsize=(8 * len(choices), 7.6), dpi=150)
    for ax, mid in zip(np.atleast_1d(axes), choices):
        portrait(ax, mid, gcps, geom, meta, exag=exag)
    fig.suptitle(f"Where the maps put places: each arrow runs from a place's true position on the map's own frame "
                 f"to where the map drew it, lengthened {exag:g} times", x=0.01, ha="left", fontsize=13, color=TEXT)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURES / "r2_portraits.png")
    plt.close(fig)


PORTRAITS = [  # a distant region early, Europe early, the same distant region late
    "https://annotations.allmaps.org/maps/9b1eab3793de8bb0",  # 1708, South America (Rumsey)
    "https://annotations.allmaps.org/maps/4d8edae715ba7e08",  # 1693, the German circles (Rumsey)
    "https://annotations.allmaps.org/maps/2af39943bf95fb58",  # 1886, South America (Rumsey)
]

if __name__ == "__main__":
    decomposition()
    stretch()
    places()
    portraits(PORTRAITS)
    print("round-2 figures written")
