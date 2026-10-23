import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from stats import boot, fit_sigmoid, sigmoid  # noqa: E402


def test_sigmoid_recovers_known_midpoint():
    rng = np.random.default_rng(1)
    t = rng.uniform(1500, 1950, 3000)
    a = sigmoid(t, 0.0, 0.4, 1760.0, 25.0) + rng.normal(0, 0.3, len(t))
    p = fit_sigmoid(t, a)
    assert p[2] == pytest.approx(1760, abs=15)
    assert p[1] == pytest.approx(0.4, abs=0.08)


def test_cluster_bootstrap_is_wider_than_naive_when_clusters_share_errors():
    rng = np.random.default_rng(2)
    cl = np.repeat(np.arange(40), 25)
    x = rng.normal(0, 1, 40)[cl] + rng.normal(0, 0.2, len(cl))
    df = pd.DataFrame({"x": x, "c": cl})
    _, (lo_c, hi_c) = boot(df, lambda d: d.x.mean(), df.c, n=500)
    _, (lo_n, hi_n) = boot(df, lambda d: d.x.mean(), np.arange(len(df)), n=500)
    assert (hi_c - lo_c) > 2.5 * (hi_n - lo_n)


def test_bootstrap_ci_covers_truth_mostly():
    hits = 0
    for s in range(60):
        rng = np.random.default_rng(100 + s)
        df = pd.DataFrame({"x": rng.normal(0.5, 1, 200)})
        _, (lo, hi) = boot(df, lambda d: d.x.mean(), np.arange(200), n=300, seed=s)
        hits += lo <= 0.5 <= hi
    assert hits >= 51  # nominal 95% of 60 = 57; allow sampling slack
