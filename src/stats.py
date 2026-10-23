"""Cluster bootstrap and the sigmoid timing model used by the hypothesis tests."""
import os
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import curve_fit

B = 2000
JOBS = max(1, (os.cpu_count() or 2) - 4)
_POOL = None


def _pool():
    global _POOL
    if _POOL is None:
        _POOL = ProcessPoolExecutor(JOBS)
    return _POOL


def cluster_indices(clusters):
    keys, inv = np.unique(np.asarray(clusters, dtype=str), return_inverse=True)
    groups = [[] for _ in keys]
    for i, k in enumerate(inv):
        groups[k].append(i)
    return [np.array(g) for g in groups]


def _draws(args):
    df, stat, groups, n, seed = args
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        pick = rng.integers(len(groups), size=len(groups))
        idx = np.concatenate([groups[k] for k in pick])
        out.append(stat(df.iloc[idx]))
    return out


def boot(df, stat, clusters, n=B, seed=0, parallel=False):
    """Point estimate and 95% percentile CI of stat(df) under cluster resampling.
    stat may return a scalar or a tuple; NaN draws are dropped. With parallel=True the
    draws are split across worker processes (stat must then be picklable)."""
    groups = cluster_indices(clusters)
    point = stat(df)
    if parallel and JOBS > 1:
        chunks = np.array_split(np.arange(n), JOBS)
        tasks = [(df, stat, groups, len(c), seed * 1000 + k) for k, c in enumerate(chunks) if len(c)]
        draws = [d for part in _pool().map(_draws, tasks) for d in part]
    else:
        draws = _draws((df, stat, groups, n, seed))
    draws = np.array(draws, dtype=float)
    if draws.ndim == 1:
        d = draws[np.isfinite(draws)]
        return point, (np.percentile(d, 2.5), np.percentile(d, 97.5)) if len(d) else (np.nan, np.nan)
    out = []
    for j in range(draws.shape[1]):
        d = draws[:, j][np.isfinite(draws[:, j])]
        out.append((np.percentile(d, 2.5), np.percentile(d, 97.5)) if len(d) else (np.nan, np.nan))
    return point, out


def sigmoid(t, a_late, d_a, t_mid, w):
    return a_late + d_a / (1.0 + np.exp((t - t_mid) / w))


BOUNDS = ([-2.0, -3.0, 1450.0, 5.0], [2.0, 3.0, 1950.0, 100.0])


def fit_sigmoid(t, a, weights=None):
    """Least-squares sigmoid A(t); returns (a_late, d_a, t_mid, w) or NaNs."""
    t, a = np.asarray(t, float), np.asarray(a, float)
    if len(t) < 20:
        return np.full(4, np.nan)
    late = a[t >= 1850].mean() if (t >= 1850).any() else a.mean()
    early = a[t < 1700].mean() if (t < 1700).any() else a.mean()
    p0 = [np.clip(late, -1.9, 1.9), np.clip(early - late, -2.9, 2.9), 1750.0, 30.0]
    sigma = None if weights is None else 1.0 / np.sqrt(np.asarray(weights, float))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p, _ = curve_fit(sigmoid, t, a, p0=p0, bounds=BOUNDS, sigma=sigma, maxfev=5000)
        return p
    except Exception:
        return np.full(4, np.nan)
