"""Run the whole pipeline in order. The two fetch steps are slow (about an hour each,
politely rate-limited) and resume where they stopped."""
import subprocess
import sys
from pathlib import Path

STEPS = [
    "download.py",            # Allmaps dump, GeoNames, Natural Earth
    "fetch_rumsey_meta.py",   # catalogue metadata for Rumsey maps (needs data/raw/maps.ndjson)
    "fetch_iiif_meta.py",     # catalogue metadata for other libraries
    "build_sample.py",        # geometric sample and GCP arrays
    "metadata.py",            # dates, makers, places
    "synthetic.py",           # K1a-c known-answer checks
    "analyze.py",             # fits, K2-K4, H1-H3, robustness
    "analyze.py --placebo",   # shuffled-year null run
    "diagnostics.py",         # post-hoc: binding fold, page profile, compass or page
    "round2_cause.py",        # round 2 (PLAN2.md): A, what causes the page excess
    "round2_decompose.py",    # B, page term vs compass term
    "round2_stretch.py",      # C, longitude stretch beyond the projection
    "round2_places.py",       # D, the same places through time
    "round2_lineage.py",      # E, copy lineages and better statistics
    "fetch_mapwarper.py",     # F data: MapWarper.net maps (slow; resumes)
    "fetch_new_allmaps.py",   # F data: Allmaps maps added after the round-1 snapshot
    "round2_replication.py",  # F, replication on maps round 1 never saw
    "figures.py",
    "figures_r2.py",
    "docs_data.py",
]

here = Path(__file__).parent
for step in STEPS:
    script, *args = step.split()
    print(f"== {step}", flush=True)
    subprocess.run([sys.executable, str(here / script), *args], check=True)
