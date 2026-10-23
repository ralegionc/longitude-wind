# Data

The raw inputs live in `data/raw/` and are not tracked in git (about 470 MB). Every file there
is either a public download or a metadata harvest that the scripts rebuild.
`data_manifest.csv` records each file's exact byte length and a fingerprint. That lets
you check a re-download rather than assume it matches. Files of 32 MiB or less have a
full `sha256`. Larger files have `headtail256`: a sha256 over the first mebibyte, the
last mebibyte and the exact byte length.

| File | What it is | Source | Rebuilt by |
|---|---|---|---|
| `maps.ndjson` | every map georeferenced in Allmaps, with ground control points (snapshot of 2026-10-13, 65,552 maps) | https://files.allmaps.org/maps.ndjson (CC0) | `src/download.py` |
| `rumsey_meta.jsonl` | catalogue records for 5,735 David Rumsey maps | Luna JSON endpoint `as/fetchMediaSearch` | `src/fetch_rumsey_meta.py` |
| `iiif_links.jsonl`, `iiif_meta.jsonl` | which manifest or CONTENTdm record holds each other map's metadata, and the records themselves | each library's IIIF or CONTENTdm service | `src/fetch_iiif_meta.py` |
| `cities15000.zip`, `cities15000.txt` | gazetteer for places of publication | GeoNames (CC BY 4.0) | `src/download.py` |
| `ne/ne_50m_land.zip` | land outlines for the figures | Natural Earth (public domain) | `src/download.py` |
| `maps_2026-10-23.ndjson` | a second Allmaps dump, for maps created after the round-1 snapshot (round 2, part F) | https://files.allmaps.org/maps.ndjson (CC0) | downloaded by hand on 2026-10-23, see PLAN2.md |
| `r2_new_meta.jsonl` | catalogue records for those new maps | as above, per library | `src/fetch_new_allmaps.py` |
| `mapwarper_maps.jsonl`, `mapwarper_gcps.jsonl` | every MapWarper.net map record, and control points for the dated regional ones | https://mapwarper.net/api/v1 | `src/fetch_mapwarper.py` |

Allmaps regenerates its dump twice a day, so a fresh download holds more maps than the
snapshot used for the README. To reproduce the exact numbers, use a file matching the
manifest fingerprint. Otherwise expect small differences in counts and intervals.

The metadata harvest appends to its files and resumes where it stopped. Some records
appear twice, once as an error and once as a later success. `src/metadata.py` keeps the
successful one.

`data/derived/` is also untracked. It holds the GCP arrays, the per-place errors and the
joined metadata, all rebuilt by `src/run_all.py`.
