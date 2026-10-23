import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from paper import contrasts, decompose  # noqa: E402


def simulate(a_p, a_c, theta, rng):
    """Residuals with covariance I + a_p diag(1,-1) + a_c J(theta), one point per theta."""
    c, s = np.cos(2 * theta), np.sin(2 * theta)
    S = np.empty((len(theta), 2, 2))
    S[:, 0, 0] = 1 + a_p + a_c * c
    S[:, 1, 1] = 1 - a_p - a_c * c
    S[:, 0, 1] = S[:, 1, 0] = a_c * s
    L = np.linalg.cholesky(S)  # batched
    return np.einsum("nij,nj->ni", L, rng.standard_normal((len(theta), 2)))


@pytest.mark.parametrize("a_p, a_c", [(0.2, 0.0), (0.0, 0.2), (0.15, 0.1), (0.0, 0.0)])
def test_decomposition_recovers_both_terms(a_p, a_c):
    rng = np.random.default_rng(3)
    theta = rng.uniform(-0.8, 0.8, 40000)  # tilts up to about 45 degrees
    res = simulate(a_p, a_c, theta, rng)
    p, q = contrasts(res, theta)
    est = decompose(p, q, theta)
    assert est[0] == pytest.approx(a_p, abs=0.02)
    assert est[1] == pytest.approx(a_c, abs=0.02)


def test_untilted_maps_only_identify_the_sum():
    rng = np.random.default_rng(4)
    theta = np.zeros(20000)
    res = simulate(0.1, 0.1, theta, rng)
    p, q = contrasts(res, theta)
    assert np.mean(p) == pytest.approx(0.2, abs=0.02)
