"""PLAN2 part A: what causes the page excess? A1 clicking, A2 printing method."""
import json
import pickle
import re
from collections import defaultdict
from functools import partial

import numpy as np
import pandas as pd

from analyze import load
from common import DERIVED, RAW, RESULTS, iter_maps
from fit import fit_map
from paper import paper_residuals
from stats import boot

KM = 111.2


def click_pairs():
    """Same-feature control points placed independently on the same scan."""
    by_img = defaultdict(list)
    for m in iter_maps():
        if len(m["gcps"]) >= 4:
            by_img[m["resource"]["id"]].append(m)
    rows = []
    for img, recs in by_img.items():
        if len(recs) < 2:
            continue
        W = recs[0]["resource"]["width"]
        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                a, b = recs[i]["gcps"], recs[j]["gcps"]
                ga, gb = np.array([g["geo"] for g in a]), np.array([g["geo"] for g in b])
                pa, pb = np.array([g["resource"] for g in a], float), np.array([g["resource"] for g in b], float)
                coslat = np.cos(np.radians(ga[:, 1]))[:, None]
                dkm = KM * np.hypot((ga[:, None, 0] - gb[None, :, 0]) * coslat, ga[:, None, 1] - gb[None, :, 1])
                k = dkm.argmin(1)
                for n in range(len(a)):
                    if dkm[n, k[n]] > 0.2:
                        continue
                    d = pa[n] - pb[k[n]]
                    if np.hypot(*d) <= 0.5 or np.hypot(*d) > 0.03 * W:  # copied point, or a different feature
                        continue
                    rows.append({"image": img, "W": W, "du": d[0], "dv": d[1]})
    return pd.DataFrame(rows).drop_duplicates()


def a_click(d):
    return 0.5 * np.log(np.sum(d.du ** 2) / np.sum(d.dv ** 2))


def click_excess(d):
    """Horizontal minus vertical variance of a single click, in units of image width squared."""
    return np.mean((d.du / d.W) ** 2 - (d.dv / d.W) ** 2) / 2


def page_rows(df, gcps):
    rows = []
    for g in df.itertuples():
        pix, geo = gcps[g.map_id]
        r = fit_map(pix, geo)
        if r is None or not g.width:
            continue
        res, _ = paper_residuals(r)
        rows.append(pd.DataFrame({"map_id": g.map_id, "cluster": g.cluster, "ru": res[:, 0] / g.width,
                                  "rv": res[:, 1] / g.width}))
    return pd.concat(rows)


def page_excess(d):
    return np.mean(d.ru ** 2 - d.rv ** 2)


METHODS = [("wax", r"\bwax"), ("engraved", r"engrav|copper ?-?plate|intaglio|etched|etching"),
           ("lithograph", r"lithograph|\blitho"), ("manuscript", r"manuscript|hand-drawn|hand drawn|pen and ink")]


def method_of(text):
    t = text.lower()
    hits = [name for name, pat in METHODS if re.search(pat, t)]
    if "wax" in hits:
        return None  # relief printing: neither class
    hits = [h for h in hits if h != "wax"]
    return hits[0] if len(hits) == 1 else None


def stratified_diff(d, a, b, lo=1450, hi=1951, step=50, weight="total"):
    num = den = 0.0
    for s in range(lo, hi, step):
        g = d[(d.year >= s) & (d.year < s + step)]
        x, y = g[g.method.isin(a)].A_page, g[g.method.isin(b)].A_page
        if len(x) >= 3 and len(y) >= 3:
            w = len(x) + len(y) if weight == "total" else len(x)
            num += w * (x.mean() - y.mean())
            den += w
    return num / den if den else np.nan


def main():
    L = ["ROUND 2, PART A: what causes the page excess?", ""]
    # A1 clicking
    pairs = click_pairs()
    pairs.to_csv(RESULTS / "r2_click_pairs.csv", index=False)
    p, ci = boot(pairs, a_click, pairs.image)
    ce, ce_ci = boot(pairs, click_excess, pairs.image)
    df = load()
    df = df[df.eligible & (df.extent_deg >= 4)].copy()
    df["cluster"] = df.maker.fillna(df.map_id).astype(str)
    gcps = pickle.load(open(DERIVED / "gcps.pkl", "rb"))
    pr = page_rows(df, gcps)
    pe, pe_ci = boot(pr, page_excess, pr.cluster)
    share = ce / pe
    L.append(f"A1 click pairs: {len(pairs)} on {pairs.image.nunique()} scans")
    L.append(f"   A_click (horizontal vs vertical disagreement) {p:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]")
    L.append(f"   single-click excess variance {ce:+.2e} [{ce_ci[0]:+.2e}, {ce_ci[1]:+.2e}] (image widths^2)")
    L.append(f"   page excess variance of round-1 residuals {pe:+.2e} [{pe_ci[0]:+.2e}, {pe_ci[1]:+.2e}]")
    L.append(f"   share of the page excess that clicking could explain: {share:.2f}")
    verdict = ("clicking explains it" if ci[0] > 0 and share >= 0.5 else
               "clicking does not explain it" if (ci[0] <= 0 <= ci[1]) or ci[1] < 0.05 else
               "clicking is anisotropic but too small to explain it" if ci[0] > 0 else "inconclusive")
    L.append(f"   -> {verdict}")
    L.append("")

    # A2 printing method
    fr = pd.read_csv(RESULTS / "frames.csv")
    st = pd.read_csv(RESULTS / "map_stats.csv")[["map_id", "rumsey_mid", "extent_deg"]]
    notes = {}
    for line in open(RAW / "rumsey_meta.jsonl", encoding="utf-8"):
        r = json.loads(line)
        notes[r["mid"]] = " ".join(str(r.get(k) or "") for k in ("Pub Note", "Note", "Type", "Full Title"))
    d = fr.merge(st, on="map_id")
    d = d[d.rumsey_mid.notna()].copy()
    d["method"] = [method_of(notes.get(m, "")) for m in d.rumsey_mid]
    d = d[d.method.notna()]
    d.to_csv(RESULTS / "r2_print_methods.csv", index=False)
    L.append("A2 printing method (Rumsey notes), page anisotropy A_page:")
    for m_, g in d.groupby("method"):
        L.append(f"   {m_:11s} n={len(g):5d} median year {g.year.median():.0f} mean A_page {g.A_page.mean():+.3f}")
    c19 = d[(d.year >= 1800) & (d.year < 1900)]
    p1, c1 = boot(c19, partial(stratified_diff, a=["engraved"], b=["lithograph"], lo=1800, hi=1900),
                  c19.cluster)
    p2, c2 = boot(d, partial(stratified_diff, a=["manuscript"], b=["engraved", "lithograph"], weight="first"),
                  d.cluster)
    L.append(f"   engraved minus lithograph, 1800-1899 within 50-year bins: {p1:+.3f} [{c1[0]:+.3f}, {c1[1]:+.3f}]")
    L.append(f"   manuscript minus printed, within 50-year bins: {p2:+.3f} [{c2[0]:+.3f}, {c2[1]:+.3f}]")
    paper = c1[0] > 0 and c2[1] < 0
    lettering = (c1[0] <= 0 <= c1[1]) and (c2[0] <= 0 <= c2[1])
    L.append("   -> " + ("consistent with paper shrinkage" if paper else
                         "consistent with lettering (no difference by method)" if lettering else
                         "neither pattern cleanly"))
    text = "\n".join(L)
    (RESULTS / "r2_cause.txt").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
