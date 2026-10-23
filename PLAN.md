# Longitude Wind: pre-registered plan

Written 2026-10-14, before any residual, anisotropy or stretch statistic was computed.
The only things seen at this point are the Allmaps record counts per provider and a
handful of catalogue dates printed by the metadata fetcher.

## Question

Before the 18th century, latitude could be measured from the sun and stars but
longitude had to be estimated from travel times. If that asymmetry shaped maps, the
positional errors on old maps should run mostly east-west, and the excess should
disappear once longitude could be measured: first on land in Europe (Jupiter's moons,
Paris Observatory, about 1680 to 1720), later on distant coasts (chronometers and
lunar-distance tables, about 1765 to 1830).

## Data

- Allmaps open dataset, `maps.ndjson` snapshot of 2026-10-13 (CC0): every map
  georeferenced in Allmaps, with ground control points (GCPs) linking image pixels to
  modern longitude and latitude.
- Dates, titles, authors and places of publication: David Rumsey Map Collection via
  its public Luna JSON endpoint; other libraries via their IIIF manifests or CONTENTdm
  records. Gallica (Cloudflare challenge) and NYPL (no manifest link in the record)
  are not used.

## Unit and inclusion

One georeferenced map. If an image was georeferenced more than once, keep the record
with the most GCPs. A map is eligible if all of these hold:

1. at least 6 GCPs after removing duplicates (same pixel or same coordinates);
2. GCP extent (the larger of its longitude and latitude span) of at least 4 degrees;
3. longitude span at most 120 degrees and latitude span at most 70 degrees (world and
   hemisphere maps are excluded);
4. a publication year between 1450 and 1950 with precision of 25 years or better
   (ranges are reduced to their midpoint);
5. not a facsimile or reproduction (title, type or notes containing "facsim",
   "fac-sim", "reproduc" or "reprint") and not a globe gore or a view.

## Per-map fit

1. Project the true GCP positions to an azimuthal equidistant frame centred on the GCP
   centroid (sphere, R = 6371 km).
2. Fit pixel coordinates to that frame with a polynomial: order 2 if the map has at
   least 12 GCPs, otherwise order 1 (affine).
3. Drop GCPs whose total residual exceeds 4 times the median total residual (only if at
   least 6 remain) and refit once.
4. Scale residuals by 1/sqrt(1 - leverage). Both components share the same design, so
   this cannot create anisotropy.
5. Convert each predicted position back to longitude and latitude. East-west error is
   R * dLon * cos(lat); north-south error is R * dLat.

Per-map statistics:

- **A** (primary), anisotropy: 0.5 * ln(sum EW^2 / sum NS^2). exp(A) is the ratio of
  east-west to north-south RMS error. Under isotropic error, E[A] = 0 by symmetry.
- **S**, aspect: ln(paper scale per km east / paper scale per km north), from an affine
  fit of the ground frame to pixels. Zero for a correctly drawn conformal map; positive
  when east-west distances are drawn too long (or for a plate carree "plane chart").
- **A_paper** = A + S, the anisotropy measured on the paper.
- **RelErr**: RMS total error divided by the largest GCP-to-GCP distance.

## Hypotheses

Confidence intervals are 95% cluster bootstraps (2,000 resamples) over map makers
(normalised first author, else publisher, else the map itself).

- **H1, early anisotropy.** Mean A for maps dated 1500 to 1699 is above zero, and above
  the mean for 1850 to 1950 (difference CI excludes 0).
- **H1t, timing.** A sigmoid A(t) = A_late + dA / (1 + exp((t - t_mid) / w)), with w
  between 5 and 100 years, puts t_mid between 1680 and 1830.
- **H2, two waves.** Maps of Europe (GCP centroid in longitude -25 to 45, latitude 34
  to 72) lose their anisotropy earlier than maps of everywhere else:
  t_mid(Europe) < t_mid(rest), difference CI excludes 0. Needs at least 40 maps dated
  before 1700 in each group, otherwise H2 is reported as untestable.
- **H3, knowledge horizon.** For maps dated before 1750 with a geocoded place of
  publication, A rises with the log distance from that place to the map's centroid
  (slope > 0). For maps after 1850 the slope is smaller than before 1750.
- **S0, sanity.** RelErr falls with date (Spearman rho < 0).

## Known-answer checks (run first; failures are reported, not hidden)

- **K1a, synthetic null.** Replace every eligible map's pixels with a synthetic map:
  true GCP positions drawn through a conformal projection (Mercator, Lambert conformal
  conic, stereographic, transverse Mercator), random rotation, plus isotropic paper
  noise of 1% of the map's extent. Mean A must lie within 0.03 of 0.
- **K1b, identity.** Same with non-conformal projections (plate carree, equirectangular,
  equidistant conic, Albers, Bonne, sinusoidal, Cassini). Mean A_paper must lie within
  0.03 of 0.
- **K1c, power.** Conformal synthetic maps with ground error twice as large east-west
  as north-south. Mean A must recover ln 2 = 0.69 within 0.10.
- **K2, modern floor.** Real maps dated 1900 to 1950 with |S| < 0.15: mean A within
  0.10 of 0.
- **K3, the long Mediterranean.** Maps with GCPs west of 0 E and east of 28 E between
  latitudes 30 and 46: mean S for 1477 to 1650 exceeds mean S for 1750 to 1950 by at
  least 0.15. (Ptolemy's Mediterranean spans about 62 degrees against a true 42, so
  ln 1.48 = 0.39; Delisle's maps of about 1700 corrected it.) Caveat: plane charts also
  give S > 0.
- **K4, Ptolemy.** Maps whose title names Ptolemy, dated 1477 to 1600: mean A > 0
  (Ptolemy's longitudes came from itinerary distances, his latitudes partly from
  observation). Reported only if at least 10 such maps exist.

## Robustness (reported for H1 and H2)

R1 only maps with |S| < 0.15; R2 A_paper instead of A; R3 region fixed effects; R4
extent thresholds 2 and 8 degrees; R5 maps weighted by residual degrees of freedom;
R6 affine fits only. H1 is called robust if the primary test passes and R1 and R3
agree in sign with CIs excluding 0.

## Outputs

A night poster of every GCP error vector by century, the anisotropy curve with region
fits, the Mediterranean aspect series, the horizon plot, a test summary figure, and an
interactive map in `docs/`.

## Addendum, 2026-10-14, written after the pre-registered run

Nothing above was changed after the results were seen. The notes below record what
happened next, so that a reader can tell the planned tests from the later ones.

**Fixes to the scoring code, not to the plan.** The first run printed "IN WINDOW" for
H1t because the script only checked that the sigmoid midpoint fell between 1680 and
1830. The plan asks for the anisotropy to be *lost*, so H1t and H2 now also require the
fitted change to be a decline. On the real data it is a rise, so both are scored as not
supported. The automatic one-line headline, which described the midpoint of that rise
as if it were the predicted decline, was removed.

**Metadata coverage.** Gallica sits behind a Cloudflare challenge and NYPL records carry
no manifest link, so neither is used, as planned. Of 2,516 other IIIF and CONTENTdm
records, 98 could not be read (mostly Library of Congress, which answers 403). Digital
Commonwealth answers 404 to an `Accept: application/json` header; dropping the header
fixed it.

**Post-hoc diagnostics (not pre-registered), prompted by the K2 failure.** Maps from
1900 to 1950 show east-west errors about 14% larger than north-south ones in every
region, which the synthetic checks say the pipeline does not create.

- D1, binding fold: a double-page atlas map has a vertical gutter that could only
  displace places east-west. Test: the jump in east-west residual across the page's
  vertical centre line, against three control lines.
- D2, page profile: residual size against position across the page.
- D3, compass or page: for maps whose drawn north is tilted on the page, compare the
  anisotropy along the page's own axes with the anisotropy along drawn east and north.
  An effect of longitude knowledge follows the compass; an effect of paper, lettering,
  scanning or clicking follows the page.

Descriptive tables of anisotropy by region and era, and by era and distance from the
place of publication (the full picture behind H3), were added at the same time.

## Addendum 2, 2026-10-23: a duplicate bug, fixed

The plan keeps one record per scanned image. The code compared image URLs, but VU
Amsterdam serves each scan under two URLs, so 424 maps in the primary sample (848
records) were counted twice, with identical control points. `src/build_sample.py` now
also drops records whose control points repeat another record's exactly; 1,187 records
went across all extents. Every round-1 test was re-run on the corrected sample. The
outputs from before the fix are kept in `results/before_dedupe/`.

What changed: the primary sample went from 5,187 to 4,763 maps. Every verdict is the
same except K3, which passed narrowly before (+0.152 against a 0.15 bar) and now fails
narrowly (+0.134 [+0.093, +0.219]). H1 still fails (-0.083), H3 still holds (+0.139),
and K2 still fails (+0.133).
