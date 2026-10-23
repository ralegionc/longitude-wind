"""Figures: the night poster, the anisotropy curve, the long Mediterranean, the horizon
plot and the test summary. Reads results/ and data/derived/; writes figures/."""
import re

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

from common import DERIVED, FIGURES, RAW, RESULTS
from fit import R
from stats import boot, fit_sigmoid, sigmoid

INK = "#0b0e13"
LAND = "#151b25"
COAST = "#2a3442"
TEXT = "#ebe6da"
MUTED = "#8f96a1"
EW = "#4fd6e8"      # east-west error: cyan
NS = "#f2a65a"      # north-south error: amber
EUROPE_C = "#9ad0a3"
REST_C = "#d7a1e8"
ORIENT = LinearSegmentedColormap.from_list("orient", [EW, "#e9e4d8", NS])

ERAS = [(1450, 1650, "1450 to 1649"), (1650, 1750, "1650 to 1749"),
        (1750, 1850, "1750 to 1849"), (1850, 1951, "1850 to 1950")]

mpl.rcParams.update({
    "font.family": ["Georgia", "DejaVu Serif"], "text.color": TEXT, "axes.labelcolor": TEXT,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": COAST, "figure.facecolor": INK,
    "axes.facecolor": INK, "savefig.facecolor": INK, "font.size": 11,
})


def land():
    return gpd.read_file(f"zip://{RAW / 'ne' / 'ne_50m_land.zip'}")


def segments(err, exag=1.0, cap_km=1500.0):
    """Line from each place's true position towards where its map puts it (degrees).
    exag lengthens every streak by the same factor; cap_km limits the drawn length."""
    e = err
    tot = np.hypot(e.e_ew, e.e_ns)
    f = np.minimum(exag, cap_km / tot.clip(lower=1e-9))
    dlon = np.degrees(e.e_ew * f / (R * np.cos(np.radians(e.lat))))
    dlat = np.degrees(e.e_ns * f / R)
    seg = np.stack([np.column_stack([e.lon, e.lat]), np.column_stack([e.lon + dlon, e.lat + dlat])], axis=1)
    orient = np.degrees(np.arctan2(np.abs(e.e_ns), np.abs(e.e_ew))) / 90.0   # 0 = east-west, 1 = north-south
    return seg, orient.values, tot.values


def ratio_axis(ax, ticks=(0.8, 1.0, 1.25, 1.5, 2.0)):
    ax.set_yscale("log")
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:g}×" for t in ticks])
    ax.minorticks_off()


POSTER_EXAG = 25
POSTER_SAMPLE = 10000


def poster(stats, err, curve, table, summary):
    import textwrap
    proj = ccrs.Robinson(central_longitude=10)
    fig = plt.figure(figsize=(22, 16), dpi=150)
    fig.text(0.035, 0.962, "The Longitude Wind", fontsize=42, weight="bold", color=TEXT, va="center")
    fig.text(0.035, 0.925, "\n".join(textwrap.wrap(summary["subtitle"], 150)), fontsize=15, color=MUTED,
             style="italic", va="center", linespacing=1.5)
    fig.text(0.035, 0.885, f"Each streak starts at a real place and points to where a dated map put it, drawn "
             f"{POSTER_EXAG} times longer than the error; at most {POSTER_SAMPLE:,} random places per period. "
             f"Cyan streaks are east-west errors, amber north-south.",
             fontsize=12, color=MUTED, va="center")
    lnd = land()
    pos = [(0.015, 0.56, 0.48, 0.31), (0.505, 0.56, 0.48, 0.31), (0.015, 0.245, 0.48, 0.31), (0.505, 0.245, 0.48, 0.31)]
    for (lo, hi, label), p in zip(ERAS, pos):
        ax = fig.add_axes(p, projection=proj)
        ax.set_global()
        ax.set_facecolor(INK)
        ax.spines["geo"].set_edgecolor(COAST)
        ax.add_geometries(lnd.geometry, ccrs.PlateCarree(), facecolor=LAND, edgecolor=COAST, linewidth=0.3)
        e_all = err[(err.year >= lo) & (err.year < hi)]
        e = e_all.sample(min(len(e_all), POSTER_SAMPLE), random_state=0)  # equal density across periods
        seg, orient, tot = segments(e, exag=POSTER_EXAG)
        order = np.argsort(tot)  # long streaks on top
        cols = ORIENT(orient[order])
        for lw, alpha in [(2.4, 0.05), (0.9, 0.2), (0.35, 0.8)]:
            c = cols.copy()
            c[:, 3] = alpha
            ax.add_collection(LineCollection(seg[order], colors=c, linewidths=lw, transform=ccrs.PlateCarree()))
        s = stats[(stats.year >= lo) & (stats.year < hi)]
        _, _, tot_all = segments(e_all)
        ax.text(0.0, 0.25, label, transform=ax.transAxes, fontsize=19, weight="bold", color=TEXT, va="top")
        ax.text(0.0, 0.16, f"{s.map_id.nunique():,} maps, {len(e_all):,} places\nmedian error {np.median(tot_all):,.0f} km\n"
                f"east-west / north-south {np.exp(s.A.mean()):.2f}",
                transform=ax.transAxes, fontsize=11, color=MUTED, va="top", linespacing=1.45)
    # bottom strip, left: the curve by region
    ax = fig.add_axes([0.06, 0.06, 0.4, 0.14])
    for series, col, label in [("Europe", EUROPE_C, "maps of Europe"), ("rest", REST_C, "maps of everywhere else")]:
        c = curve[(curve.series == series) & (curve.stat == "A") & (curve.n >= 15)]
        ax.fill_between(c.mid, np.exp(c.lo), np.exp(c.hi), color=col, alpha=0.15, lw=0)
        ax.plot(c.mid, np.exp(c["mean"]), color=col, lw=2, label=label)
    ax.axhline(1, color=MUTED, lw=0.7, ls=":")
    ax.set_xlim(1525, 1950)
    ratio_axis(ax, (1.0, 1.25, 1.5))
    ax.set_ylim(0.8, 1.75)
    ax.legend(frameon=False, fontsize=10, loc="upper right", ncol=2)
    ax.set_title("east-west / north-south error by year", fontsize=12, color=TEXT, loc="left")
    ax.tick_params(labelsize=9)
    # bottom strip, right: the frontier
    ax = fig.add_axes([0.56, 0.06, 0.4, 0.14])
    for lo, label, col in [(1450, "printed before 1750", NS), (1850, "1850 and later", EW)]:
        t = table[table.era_lo == lo].sort_values("dist_lo")
        x = np.arange(len(t))
        ax.fill_between(x, np.exp(t.lo), np.exp(t.hi), color=col, alpha=0.15, lw=0)
        ax.plot(x, np.exp(t.A), "-o", color=col, lw=2, ms=4, label=label)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["under 500 km", "500 to 2,000", "2,000 to 6,000", "over 6,000 km"])
    ax.axhline(1, color=MUTED, lw=0.7, ls=":")
    ratio_axis(ax, (1.0, 1.25, 1.5))
    ax.set_ylim(0.8, 1.75)
    ax.legend(frameon=False, fontsize=10, loc="upper left")
    ax.set_title("east-west / north-south error by distance from where the map was printed", fontsize=12,
                 color=TEXT, loc="left")
    ax.tick_params(labelsize=9)
    for a in fig.axes[-2:]:
        for s_ in a.spines.values():
            s_.set_color(COAST)
    fig.text(0.965, 0.012, "Data: Allmaps open georeferencing (CC0), David Rumsey Map Collection and library IIIF catalogues; "
             "control points placed by volunteers. Bands: 95% intervals over map makers.", fontsize=9, color=MUTED, ha="right")
    fig.savefig(FIGURES / "poster.png", dpi=150)
    plt.close(fig)


def anisotropy_curve(stats, curve):
    fig, ax = plt.subplots(figsize=(12, 6.2), dpi=150)
    P = stats[stats.extent_deg >= 4]
    t = np.linspace(1450, 1950, 400)
    for series, col, sub, label in [("Europe", EUROPE_C, P[P.region == "Europe"], "maps of Europe"),
                                    ("rest", REST_C, P[P.region != "Europe"], "maps of everywhere else")]:
        c = curve[(curve.series == series) & (curve.stat == "A")]
        ax.fill_between(c.mid, np.exp(c.lo), np.exp(c.hi), color=col, alpha=0.16, lw=0)
        ax.plot(c.mid, np.exp(c["mean"]), "-o", color=col, ms=4, lw=2, label=label)
        p = fit_sigmoid(sub.year, sub.A)
        ax.plot(t, np.exp(sigmoid(t, *p)), color=col, lw=1, ls="--", alpha=0.6)
    ax.plot([], [], color=MUTED, lw=1, ls="--", label="pre-registered sigmoid fit")
    ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ymax = ax.get_ylim()[1]
    for yr, txt in [(1668, "Cassini's tables of\nJupiter's moons"), (1767, "Nautical Almanac;\nchronometers spread")]:
        ax.axvline(yr, color=COAST, lw=1)
        ax.text(yr + 3, ymax * 0.985, txt, fontsize=9, color=MUTED, va="top")
    ratio_axis(ax, (0.6, 0.8, 1.0, 1.25, 1.5, 2.0))
    ax.set_xlim(1450, 1950)
    ax.set_ylabel("east-west / north-south RMS error")
    ax.set_xlabel("year the map was published (25-year bins, 95% intervals)")
    ax.legend(frameon=False, loc="lower left", fontsize=10)
    ax.set_title("Predicted: the east-west excess fades. Found: maps of everywhere converge on one excess",
                 loc="left", fontsize=14, weight="bold")
    for s_ in ax.spines.values():
        s_.set_color(COAST)
    fig.tight_layout()
    fig.savefig(FIGURES / "anisotropy_curve.png")
    plt.close(fig)


def mediterranean(stats, med_ids, curve):
    d = stats[stats.map_id.isin(med_ids)]
    fig, ax = plt.subplots(figsize=(11, 5.6), dpi=150)
    ax.scatter(d.year, np.exp(d.S), s=10, color=EW, alpha=0.35, lw=0)
    c = curve[(curve.series == "mediterranean") & (curve.stat == "S")]
    ax.fill_between(c.mid, np.exp(c.lo), np.exp(c.hi), color=EW, alpha=0.25, lw=0)
    ax.plot(c.mid, np.exp(c["mean"]), color=EW, lw=2)
    ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ax.axhline(62 / 42, color=NS, lw=1, ls="--")
    ax.text(1460, 62 / 42 * 1.02, "Ptolemy: 62 degrees of Mediterranean for a true 42", fontsize=9, color=NS)
    ax.text(1460, 1.02, "drawn true to scale", fontsize=9, color=MUTED)
    ratio_axis(ax, (0.8, 1.0, 1.25, 1.5, 2.0))
    ax.set_xlim(1450, 1950)
    ax.set_ylabel("east-west stretch (drawn / true)")
    ax.set_xlabel("year the map was published")
    ax.set_title("The long Mediterranean shrinks", loc="left", fontsize=16, weight="bold")
    for s_ in ax.spines.values():
        s_.set_color(COAST)
    fig.tight_layout()
    fig.savefig(FIGURES / "mediterranean.png")
    plt.close(fig)


def frontier(table):
    """H3 in full: east-west excess against distance from the city the map was printed in."""
    fig, ax = plt.subplots(figsize=(10.5, 5.8), dpi=150)
    eras = [(1450, "printed before 1750", NS), (1750, "1750 to 1849", "#e9e4d8"), (1850, "1850 and later", EW)]
    for lo, label, col in eras:
        t = table[table.era_lo == lo].sort_values("dist_lo")
        x = np.sqrt(t.dist_lo.clip(lower=100) * t.dist_hi.clip(upper=15000))
        ax.fill_between(x, np.exp(t.lo), np.exp(t.hi), color=col, alpha=0.14, lw=0)
        ax.plot(x, np.exp(t.A), "-o", color=col, lw=2.2, ms=5, label=label)
    ax.set_xscale("log")
    ratio_axis(ax, (0.9, 1.0, 1.25, 1.5))
    ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ax.set_xticks([200, 1000, 4000, 10000])
    ax.set_xticklabels(["under 500 km", "500 to 2,000", "2,000 to 6,000", "over 6,000 km"])
    ax.set_xlabel("distance from the city where the map was published to the place it shows")
    ax.set_ylabel("east-west / north-south error")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Longitude doubt grew with distance from home", loc="left", fontsize=15, weight="bold")
    ax.text(0.0, -0.2, "Bands are 95% cluster-bootstrap intervals over map makers. Map counts per point are in results/distance_table.csv.",
            transform=ax.transAxes, fontsize=8.5, color=MUTED)
    for s_ in ax.spines.values():
        s_.set_color(COAST)
    fig.tight_layout()
    fig.savefig(FIGURES / "frontier.png")
    plt.close(fig)


def compass_or_page(fr):
    """D3: does the excess turn with the drawn compass, or stay on the page?"""
    fr = fr.copy()
    fr["t"] = fr.tilt.abs()
    bins = [(0, 5), (5, 15), (15, 30), (30, 60), (60, 120)]
    mid, page, comp, n, pci, cci = [], [], [], [], [], []
    for lo, hi in bins:
        g = fr[(fr.t >= lo) & (fr.t < hi)]
        cl = g.maker.fillna(g.map_id).astype(str)
        mid.append(np.clip(g.t.median(), 0, 90))
        p_, ci_p = boot(g, lambda z: z.A_page.mean(), cl, n=1000)
        c_, ci_c = boot(g, lambda z: z.A_compass.mean(), cl, n=1000)
        page.append(p_)
        comp.append(c_)
        pci.append(ci_p)
        cci.append(ci_c)
        n.append(len(g))
    fig, ax = plt.subplots(figsize=(10, 5.6), dpi=150)
    th = np.linspace(0, 90, 200)
    a0 = page[0]
    # what the page measure would do if the excess were east-west on the ground: it turns with the map
    ax.plot(th, a0 * np.cos(np.radians(2 * th)), color=NS, lw=1, ls="--", alpha=0.7)
    ax.text(47, -0.135, "page measure expected if the excess\nwere east-west on the ground",
            color=NS, fontsize=8.5, alpha=0.9)
    off = 1.2  # nudge the two series apart so their interval bars do not overlap
    ax.errorbar(np.array(mid) - off, page, yerr=np.abs(np.array(pci).T - page), fmt="-o", color=NS, lw=2.2,
                elinewidth=1, capsize=0, label="measured along the page (horizontal vs vertical)")
    ax.errorbar(np.array(mid) + off, comp, yerr=np.abs(np.array(cci).T - comp), fmt="-o", color=EW, lw=2.2,
                elinewidth=1, capsize=0, label="measured along the drawn compass (east vs north)")
    for x, k in zip(mid, n):
        ax.text(x, -0.19, f"n={k}", color=MUTED, fontsize=8, ha="center")
    ax.axhline(0, color=MUTED, lw=0.8, ls=":")
    ax.set_xlim(-3, 95)
    ax.set_xlabel("how far the map's north is tilted from the top of the page (degrees, median of each group)")
    ax.set_ylabel("anisotropy A (log ratio)")
    ax.legend(frameon=False, loc="upper left", fontsize=9.5)
    ax.set_title("The modern excess follows the page more than the compass", loc="left", fontsize=15, weight="bold")
    ax.text(0.0, -0.2, "Post-hoc diagnostic. Bars are 95% cluster-bootstrap intervals; n is the number of maps per tilt group.",
            transform=ax.transAxes, fontsize=8.5, color=MUTED)
    for s_ in ax.spines.values():
        s_.set_color(COAST)
    fig.tight_layout()
    fig.savefig(FIGURES / "compass_or_page.png")
    plt.close(fig)


def tests_figure():
    """Known-answer checks: each estimate against its pre-registered target band."""
    syn = (RESULTS / "synthetic_checks.txt").read_text(encoding="utf-8")
    summ = (RESULTS / "summary.txt").read_text(encoding="utf-8")
    plac = (RESULTS / "placebo_summary.txt").read_text(encoding="utf-8") if (RESULTS / "placebo_summary.txt").exists() else ""
    rows = []
    for check, stat, target, tol in [("K1a", "A", 0.0, 0.03), ("K1b", "A_paper", 0.0, 0.03), ("K1c", "A", np.log(2), 0.10)]:
        m = re.search(rf"{check} {stat}\s+primary\s+n=\s*(\d+) mean=([+-][\d.]+) \(se ([\d.]+)\)", syn)
        v, se = float(m.group(2)), float(m.group(3))
        rows.append((f"{check} synthetic ({stat})", v, v - 1.96 * se, v + 1.96 * se, target - tol, target + tol))
    m = re.search(r"K2 .*?mean A ([+-][\d.]+) \[([+-][\d.]+), ([+-][\d.]+)\]", summ)
    if m:
        rows.append(("K2 modern maps, mean A", *map(float, m.groups()), -0.10, 0.10))
    m = re.search(r"K3 .*?difference ([+-][\d.]+) \[([+-][\d.]+), ([+-][\d.]+)\]", summ)
    if m:
        rows.append(("K3 Mediterranean stretch, early - late S", *map(float, m.groups()), 0.15, 1.0))
    m = re.search(r"K4 .*?mean A ([+-][\d.]+) \[([+-][\d.]+), ([+-][\d.]+)\]", summ)
    if m:
        rows.append(("K4 Ptolemy maps, mean A", *map(float, m.groups()), 0.0, 1.5))
    m = re.search(r"H1 .*?\n.*?\n\s+difference: ([+-][\d.]+) \[([+-][\d.]+), ([+-][\d.]+)\]", plac)
    if m:
        rows.append(("Placebo (years shuffled), H1 difference", *map(float, m.groups()), -0.05, 0.05))
    fig, ax = plt.subplots(figsize=(10, 0.65 * len(rows) + 1.2), dpi=150)
    for i, (name, v, lo, hi, a, b) in enumerate(rows[::-1]):
        ax.add_patch(plt.Rectangle((a, i - 0.3), b - a, 0.6, color="#1f3a2a", lw=0))
        ok = a <= v <= b
        ax.plot([lo, hi], [i, i], color=TEXT, lw=1.4)
        ax.plot(v, i, "o", color="#8fe39b" if ok else "#ff7a6b", ms=8)
        ax.text(-0.42, i, name, ha="right", va="center", fontsize=10)
    ax.set_yticks([])
    ax.set_xlim(-0.4, 1.2)
    ax.axvline(0, color=COAST, lw=0.8)
    ax.set_xlabel("estimate with 95% interval; green band = pre-registered pass region")
    ax.set_title("Known-answer checks", loc="left", fontsize=14, weight="bold")
    for s_ in ax.spines.values():
        s_.set_color(COAST)
    fig.subplots_adjust(left=0.38, right=0.97, top=0.88, bottom=0.12)
    fig.savefig(FIGURES / "tests.png")
    plt.close(fig)


def main():
    stats = pd.read_csv(RESULTS / "map_stats.csv")
    err = pd.read_csv(DERIVED / "gcp_errors.csv")
    curve = pd.read_csv(RESULTS / "curve.csv")
    table = pd.read_csv(RESULTS / "distance_table.csv")
    fr = pd.read_csv(RESULTS / "frames.csv")
    early = table[table.era_lo == 1450]
    far, near = early[early.dist_lo == 6000].ratio.iloc[0], early[early.dist_lo == 0].ratio.iloc[0]
    late = np.exp(stats[(stats.extent_deg >= 4) & (stats.year >= 1850)].A.mean())
    st2 = pd.read_csv(RESULTS / "r2_stretch.csv")
    pre, post = st2[st2.year < 1750].abs_lnk.mean(), st2[st2.year >= 1750].abs_lnk.mean()
    subtitle = (f"Before 1750, maps drew longitudes about {100 * (np.exp(pre) - 1):.0f}% too wide or too narrow, and put "
                f"places over 6,000 km from where they were printed {far:.2f} times as far out east-west as north-south "
                f"(near home, {near:.2f}). After 1750 the stretch fell to {100 * (np.exp(post) - 1):.0f}%, but every map "
                f"kept a small east-west excess, about {late:.2f}.")
    band = err[err.lat.between(30, 46)]
    med_ids = set(band[band.lon < 0].map_id) & set(band[band.lon > 28].map_id)
    prim = stats[stats.extent_deg >= 4]  # the pre-registered primary sample, as on the web page
    poster(prim, err[err.map_id.isin(prim.map_id)], curve, table, {"subtitle": subtitle})
    anisotropy_curve(stats, curve)
    mediterranean(stats, med_ids, curve)
    frontier(table)
    compass_or_page(fr)
    tests_figure()
    print("figures written")


if __name__ == "__main__":
    main()
