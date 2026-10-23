# Longitude Wind, round 2: pre-registered plan

Written 2026-10-23, after round 1 (PLAN.md, results in README) and before any round-2
statistic was computed. Round 1 stays frozen. Round 2 asks the questions round 1 left
open, and every test below is fixed here before it runs.

What round 1 showed, in one line each: the pre-registered decline in east-west excess
did not happen (H1); the excess grew with distance from the place of publication before
1750 (H3); maps from 1900 to 1950 still carry a 14% excess (K2 failed); and on maps with
tilted north that excess followed the page, not the compass (post hoc).

Shared conventions: the round-1 per-map fit (`src/fit.py`), the round-1 primary sample
unless stated, 95% cluster bootstraps over map makers with 2,000 resamples, and a
"pass" that needs the whole interval on the predicted side unless a tolerance is given.

## A. What causes the page excess?

**A1, clicking.** Some scans were georeferenced twice, in separate Allmaps records. Pairs
of control points whose modern coordinates lie within 200 m of each other (the same
feature) and whose scan positions differ (not copied) measure how much two volunteers
disagree when they click the old map. Disagreement anisotropy on the page:
A_click = 0.5 ln(sum du^2 / sum dv^2), u horizontal, v vertical, after dropping pairs
more than 3% of the image width apart (different features). Bootstrap over images; all
map extents are used, because clicking does not depend on extent.
- If clicking causes the page excess: A_click > 0, and the click disagreement can
  account for at least half the page excess in variance terms.
- Otherwise A_click is near 0 (interval includes 0, or upper limit below 0.05).

**A2, printing.** Rumsey notes name the printing method. Intaglio printing presses damp
paper, which then shrinks unevenly; lithographs are printed with less water, and
manuscripts never went through a press. Page anisotropy A_page (round-1 D3 measure) for
engraved, lithographed and manuscript maps, compared within 50-year bins and pooled by
the number of maps per bin.
- If paper shrinkage causes it: engraved > lithograph within 1800 to 1899, and
  manuscript < printed (engraved and lithograph together).
- If lettering causes it, all three are alike.

## B. Separating page from compass on every map

For each control point, take the residual on the paper (u, v), scaled by the map's RMS,
and the local direction of drawn east, at angle theta from page-right. Model the
residual covariance as proportional to I + a_p J0 + a_c J(theta), with J0 = diag(1, -1)
along the page and J(theta) the same contrast along drawn east and north. Then for each
point, with s^2 the map's mean squared residual per axis:

- p = (r_u^2 - r_v^2) / (2 s^2), with E[p] = a_p + a_c cos 2theta
- q = 2 r_u r_v / (2 s^2), with E[q] = a_c sin 2theta

Stacking p and q, a_p and a_c come from least squares. Tilted maps and the turning
meridians of large conic maps separate the two terms. Only maps with |S| < 0.15 are used,
so paper and ground units agree.

**Known answers first (B0).** Rebuild the round-1 primary maps as synthetic maps, with
real layouts and real local tilts, and inject (i) page anisotropy 0.2 only, (ii) compass
anisotropy 0.2 only, (iii) both. The estimator must recover each term within 0.04, and
leak less than 0.04 into the other term.

Tests:
- **B1**: in 1850 to 1950, a_p > 0, and a_p > a_c.
- **B2** (round-1 H1 asked of the compass term alone): a_c for 1500 to 1749 exceeds a_c
  for 1850 to 1950.
- **B3** (round-1 H3 asked of the compass term alone): before 1750, a_c for places over
  2,000 km from the place of publication exceeds a_c for places within 500 km.

## C. Stretch, not just scatter

For maps spanning at least 20 degrees, fit each of seven projections: equirectangular
true at the centre latitude, Mercator, Lambert conformal conic and equidistant conic
(standard parallels at 1/6 and 5/6 of the latitude range), Bonne (centre latitude),
sinusoidal, and stereographic. All are centred on the map. Before projecting, longitudes
are scaled about the centre by a free factor k, and a similarity transform maps the
result to the page. Keep the projection with the smallest residual; ln k is the
longitude stretch beyond any honest projection. A "plane chart" (equal degrees on the
page) reads as ln k near -ln cos(centre latitude). Those maps are flagged and reported
separately; they stay in the main test.

- **C0 known answers**: synthetic maps drawn in each family with k in {0.8, 1, 1.25}, plus
  1% isotropic paper noise. ln k must be recovered within 0.03 on average for every
  family, and the right family chosen at least 70% of the time.
- **C1**: mean |ln k| for 1500 to 1699 exceeds mean |ln k| for 1850 to 1950.
- **C2**: Europe loses its stretch before everywhere else, judged by the year a 25-year
  binned mean |ln k| first falls below 0.05 and stays there.
- **C3**: replication of K3 with this fit: Mediterranean maps before 1650 have mean ln k
  above 0.1.

## D. The same places through time

Places are reference points: the round-1 control points within 15 km of a reference
city's coordinates. Error is that place's east-west and north-south residual in km.

- **D1 known answer, the event study.** Five places have a documented first
  astronomical longitude: Cayenne (Richer, 1672), Alexandria (Chazelles, 1694), Lima
  (Feuillee, 1709), Cape Town (Lacaille, 1751) and Tahiti (1769). For each, take
  median |east-west error| and median |north-south error| in maps dated within 60 years
  before and within 60 years after. The prediction is that the east-west error falls more
  than the north-south error: a negative difference in log ratios, pooled over places
  with at least 5 maps on each side.
- **D2 figure** (descriptive): error trajectories for about 20 well-mapped ports and
  capitals.

## E. Copying and better statistics

**E1, lineages.** Two maps with at least 5 places within 15 km of each other get a
residual correlation, computed over the stacked east and north components. The null is
pairs of maps of the same region at least 150 years apart. A pair counts as a copy if
its correlation exceeds the null's 99th percentile. Lineages are connected components of
the copy graph.

**E2.** H3 re-estimated with lineages as bootstrap clusters. Prediction: the
before-1750 minus after-1850 slope difference stays above 0.

**E3.** H3 re-estimated as a mixed model of per-map A on log distance by era, with random
intercepts for maker and for atlas (Rumsey publication number). Each map is weighted by
its residual degrees of freedom. Prediction: the slope difference is above 0, with
p < 0.05.

## F. Replication on maps round 1 never saw

**Data.** MapWarper.net maps (open API, independent volunteers, a different clicking
interface), dated by `date_depicted` or else a year in the title, warped, at least
4 degrees across, 1450 to 1950. Plus Allmaps maps first created after 2026-10-13. MapWarper
has no place of publication, so H3 cannot be tested there.

- **F1**: 1850 to 1950 MapWarper maps show A > 0 (the page excess, through another tool).
- **F2**: maps outside Europe lose east-west excess after 1850: mean A before 1850 >
  mean A after 1850.
- **F3**: Europe before 1700 shows no excess (interval includes 0, as in round 1).
- **F4**: the B decomposition on MapWarper maps gives a_p > 0 for 1850 to 1950.

Clusters for MapWarper are single maps (no maker field). Results are reported whatever
they are.

## Addendum, 2026-10-23, written after the round-2 runs

- **The duplicate bug** (PLAN.md addendum 2) was found while choosing maps for the
  portraits, after parts A to E had run once. All parts were re-run on the corrected
  sample, and the numbers reported are from that re-run. The earlier outputs are in
  `results/before_dedupe/`.
- **C0** fails on the family criterion (66% right, against 70%): conic, Bonne and
  stereographic maps look alike over 20 degrees. The quantity used in C1 to C3, ln k, is
  recovered to within 0.003 in every family. C1 to C3 are reported as planned, with this
  caveat.
- **E2**: chaining copy pairs at the pre-registered threshold joins 3,974 of 4,763 maps
  into one cluster, which makes the bootstrap meaningless (interval +0.08 to +0.67). A
  post-hoc repeat uses only near-identical pairs (r > 0.9), giving 1,182 clusters with
  the largest at 718 maps.
- **E3**: statsmodels' mixed model takes no weights. The pre-registered weighted mixed
  model is therefore replaced by two fits: an unweighted mixed model (maker random
  intercept, atlas variance component), and a dof-weighted regression with errors
  clustered two ways (by maker and by atlas).
- **D1** is underpowered: only Alexandria has 5 or more dated maps on each side of its
  event. It is reported as run.
- **F data**: the MapWarper API ignores JSON:API `page[number]` and returns page 1 every
  time. The first crawl was discarded and redone with `page=` and `per_page=`.
