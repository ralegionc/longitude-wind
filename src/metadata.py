"""Join catalogue metadata to the geometric sample: year, maker, place, flags.

Writes data/derived/meta.csv with one row per map in maps_geom.csv.
"""
import json
import re
import unicodedata

import numpy as np
import pandas as pd

from common import DERIVED, RAW
from dates import parse_year

FACSIMILE = re.compile(r"facsim|fac-sim|reproduc|reprint|herdruk|nachdruck|r[ée]impression|re-issue of", re.I)
VIEW_GLOBE = re.compile(r"\bview\b|globe|gores", re.I)
DATE_LABEL = re.compile(r"date|datum|jaar|year|created|issued|publication|erschein|entstehung|fecha|ann[ée]e|dating|datering|period", re.I)
NOT_DATE = re.compile(r"availab|record|process|digiti|modif|online|catalog|upload|accession|added|epoch|copyright", re.I)
PLACE_LABEL = re.compile(r"place|plaats|\bort\b|lieu|origin", re.I)
IMPRINT_LABEL = re.compile(r"uitgave|imprint|published|publisher|publication|verlag", re.I)
NOT_PLACE = re.compile(r"digital|holding|physical|record|signatuur|land van|ppn|identifier|repository", re.I)
MAKER_LABEL = re.compile(r"creator|cartograph|author|vervaardiger|maker|auteur|kartograph|contributor|publisher|uitgever", re.I)

# Latin, Dutch and French imprint forms seen on early maps -> modern GeoNames name.
ALIASES = {
    "amstelodami": "Amsterdam", "amstelaedami": "Amsterdam", "amstelodamum": "Amsterdam",
    "lugduni batavorum": "Leiden", "leide": "Leiden", "leyde": "Leiden", "leyden": "Leiden",
    "parisiis": "Paris", "lutetiae": "Paris", "lutetiae parisiorum": "Paris",
    "londini": "London", "norimbergae": "Nuremberg", "nurnberg": "Nuremberg", "nuernberg": "Nuremberg",
    "augustae vindelicorum": "Augsburg", "venetiis": "Venice", "venezia": "Venice", "venise": "Venice",
    "romae": "Rome", "roma": "Rome", "antverpiae": "Antwerp", "anvers": "Antwerp", "antwerpen": "Antwerp",
    "coloniae": "Cologne", "coloniae agrippinae": "Cologne", "koln": "Cologne", "basileae": "Basel",
    "argentorati": "Strasbourg", "ulmae": "Ulm", "lugduni": "Lyon", "lyons": "Lyon",
    "s-gravenhage": "The Hague", "s gravenhage": "The Hague", "den haag": "The Hague", "la haye": "The Hague",
    "hagae comitis": "The Hague", "haag": "The Hague", "francofurti": "Frankfurt am Main",
    "frankfurt": "Frankfurt am Main", "frankfurt a. m": "Frankfurt am Main", "frankfurt a.m": "Frankfurt am Main",
    "frankfurt a.m.": "Frankfurt am Main", "frankfurt am main": "Frankfurt am Main", "lipsiae": "Leipzig", "vindobonae": "Vienna", "wien": "Vienna",
    "st. petersburg": "Saint Petersburg", "st petersburg": "Saint Petersburg", "washington d.c.": "Washington",
    "washington": "Washington", "new york": "New York City", "new-york": "New York City",
    "nurenberg": "Nuremberg", "mainz": "Mainz", "moguntiae": "Mainz", "duisburg": "Duisburg",
    "edinburgh": "Edinburgh", "dublin": "Dublin", "madrid": "Madrid", "lisboa": "Lisbon", "lisbonne": "Lisbon",
    "berlin": "Berlin", "gotha": "Gotha", "weimar": "Weimar", "milano": "Milan", "firenze": "Florence",
    "glogau": "Glogow", "berolini": "Berlin", "n. york": "New York City", "newyork": "New York City",
    "phiadelphia": "Philadelphia", "philadelphiae": "Philadelphia", "londres": "London", "londra": "London", "cornhill": "London", "patis": "Paris", "v praze": "Prague", "praha": "Prague",
}
# publishing towns below the gazetteer's 15,000 population floor
SMALL_TOWNS = {"hildburghausen": (50.426, 10.729), "gotha": (50.949, 10.701)}


def ascii_key(s):
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower().replace("'", "").replace('"', "")
    return re.sub(r"\s+", " ", re.sub(r"[\[\]\(\)\?\.]+$|^[\[\(]+", "", s.strip())).strip(" ,;:")


def load_gazetteer():
    cols = ["id", "name", "ascii", "alt", "lat", "lon", "fc", "fcode", "cc", "cc2", "a1", "a2", "a3", "a4",
            "pop", "elev", "dem", "tz", "mod"]
    g = pd.read_csv(RAW / "cities15000.txt", sep="\t", names=cols, usecols=["name", "ascii", "alt", "lat", "lon", "pop"],
                    quoting=3, keep_default_na=False, dtype={"alt": str})
    by_name, by_alt = {}, {}
    for r in g.itertuples():
        for key in {ascii_key(r.name), ascii_key(r.ascii)}:
            if key and (key not in by_name or by_name[key][2] < r.pop):
                by_name[key] = (r.lat, r.lon, r.pop)
        if r.pop >= 50000 and r.alt:
            for a in r.alt.split(","):
                key = ascii_key(a)
                if key and (key not in by_alt or by_alt[key][2] < r.pop):
                    by_alt[key] = (r.lat, r.lon, r.pop)
    return by_name, by_alt


def geocode(place, gaz):
    if not place:
        return None
    by_name, by_alt = gaz
    first = re.split(r"[;,:/]| and | & |\bund\b|\bet\b", str(place))[0]
    key = ascii_key(first)
    key = re.sub(r"^((prostant|apud|chez|bij|by|published|gedruckt|impensis|sumptibus|in|at|a|te|zu)\s+)+", "", key)
    if not key or key in ("s.l", "s. l", "sine loco", "n.p", "o.o", "null"):
        return None
    if key in SMALL_TOWNS:
        return SMALL_TOWNS[key]
    key = ascii_key(ALIASES.get(key, key))
    hit = by_name.get(key) or by_alt.get(key)
    return (hit[0], hit[1]) if hit else None


def first(v):
    if isinstance(v, list):
        return v[0] if v else None
    return v


def norm_maker(s):
    if not s:
        return None
    s = first(s)
    s = re.split(r";| and | & |\|", str(s))[0]
    s = re.sub(r"\d{3,4}(-\d{2,4})?|\(.*?\)|\[.*?\]", " ", s)
    s = re.sub(r"[^a-z ]", " ", ascii_key(s))
    s = " ".join(s.split()[:3])
    return s or None


def rumsey_rows():
    out = {}
    for line in open(RAW / "rumsey_meta.jsonl", encoding="utf-8"):
        r = json.loads(line)
        if "error" in r or r.get("missing"):
            continue
        date, src = parse_year(first(r.get("Date"))), "Date"
        if not date:
            date, src = parse_year(first(r.get("Pub Date"))), "Pub Date"
        text = " ".join(str(r.get(k) or "") for k in ("Full Title", "Pub Title", "Pub Note", "Note", "Type", "Pub Type"))
        typ = json.dumps(r.get("Type"))
        out[r["mid"]] = {
            "year": date[0] if date else np.nan, "year_prec": date[1] if date else np.nan,
            "date_src": src if date else None, "title": first(r.get("Full Title")) or r.get("displayName"),
            "maker": norm_maker(r.get("Author") or r.get("Publisher") or r.get("Publication Author")),
            "place": first(r.get("Publisher Location")), "type": typ,
            "facsimile": bool(FACSIMILE.search(text)), "view_globe": bool(VIEW_GLOBE.search(typ)),
        }
    return out


def iiif_rows():
    meta = {}
    for line in open(RAW / "iiif_meta.jsonl", encoding="utf-8"):
        r = json.loads(line)
        if "error" not in r:
            meta[(r["kind"], r["key"])] = r
    out = {}
    for line in open(RAW / "iiif_links.jsonl", encoding="utf-8"):
        link = json.loads(line)
        src = tuple(link["source"]) if link["source"] else None
        r = meta.get(src)
        if not r:
            continue
        pairs = [(kv["label"], kv["value"]) for kv in r.get("metadata", [])]
        date, dsrc = None, None
        for lab, val in pairs:
            if DATE_LABEL.search(lab) and not NOT_DATE.search(lab):
                date = parse_year(val)
                if date:
                    dsrc = lab
                    break
        if not date and r.get("navDate"):
            date, dsrc = parse_year(r["navDate"][:4]), "navDate"
        if not date:
            date, dsrc = parse_year(r.get("label")), "title"
        place = next((v for lab, v in pairs if PLACE_LABEL.search(lab) and not NOT_PLACE.search(lab) and v), None) \
            or next((v for lab, v in pairs if IMPRINT_LABEL.search(lab) and not NOT_PLACE.search(lab) and v), None)
        if place and ":" in place:  # "Place : Publisher, year" imprints
            place = place.split(":")[0]
        maker = next((v for lab, v in pairs if MAKER_LABEL.search(lab) and v), None)
        text = " ".join([r.get("label") or ""] + [v for _, v in pairs])
        out[link["map"]] = {
            "year": date[0] if date else np.nan, "year_prec": date[1] if date else np.nan,
            "date_src": dsrc if date else None, "title": r.get("label"), "maker": norm_maker(maker),
            "place": place, "type": None, "facsimile": bool(FACSIMILE.search(text)),
            "view_globe": bool(re.search(r"globe|gores|\bgezicht\b|\bvue\b|bird's-eye", text, re.I)),
        }
    return out


def main():
    geom = pd.read_csv(DERIVED / "maps_geom.csv")
    rum, iiif = rumsey_rows(), iiif_rows()
    gaz = load_gazetteer()
    rows = []
    for g in geom.itertuples():
        rec = rum.get(g.rumsey_mid) if isinstance(g.rumsey_mid, str) else iiif.get(g.map_id)
        rec = dict(rec) if rec else {}
        loc = geocode(rec.get("place"), gaz)
        rec.update({"map_id": g.map_id, "has_meta": bool(rec),
                    "place_lat": loc[0] if loc else np.nan, "place_lon": loc[1] if loc else np.nan})
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(DERIVED / "meta.csv", index=False)
    print(f"{len(df)} maps; metadata {df.has_meta.sum()}; dated {df.year.notna().sum()}; "
          f"place geocoded {df.place_lat.notna().sum()} of {df.place.notna().sum()}; "
          f"facsimile {df.facsimile.sum()}; view/globe {df.view_globe.sum()}")
    print(df.date_src.value_counts().head(12).to_string())
    unmatched = df[df.place.notna() & df.place_lat.isna()].place.value_counts().head(40)
    (DERIVED / "places_unmatched.txt").write_text(unmatched.to_string(), encoding="utf-8")


if __name__ == "__main__":
    main()
