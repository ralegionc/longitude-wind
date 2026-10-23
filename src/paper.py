"""Residuals on the paper, with the local drawn compass, for the round-2 analyses.

`paper_residuals(r)` takes a round-1 fit (fit.fit_map output) and refits ground ->
pixel with the same polynomial order. It returns, for each kept control point, the
residual in pixels (u to the right, v down the page) and the angle theta of drawn east
from page-right at that point.
"""
import numpy as np

from fit import aeqd_forward, design


def forward_fit(r):
    lon, lat, u, v = r["gcp_lon"], r["gcp_lat"], r["gcp_u"], r["gcp_v"]
    X, Y = aeqd_forward(lon, lat, r["lon0"], r["lat0"])
    mu = np.array([X.mean(), Y.mean()])
    sc = np.std(np.r_[X, Y]) or 1.0

    def D(lo, la):
        x, y = aeqd_forward(lo, la, r["lon0"], r["lat0"])
        return design((x - mu[0]) / sc, (y - mu[1]) / sc, r["order"])

    coef, *_ = np.linalg.lstsq(D(lon, lat), np.column_stack([u, v]), rcond=None)
    return D, coef


def paper_residuals(r, step=0.05):
    lon, lat = r["gcp_lon"], r["gcp_lat"]
    D, coef = forward_fit(r)
    X = D(lon, lat)
    pred = X @ coef
    h = np.clip(np.einsum("ij,jk,ik->i", X, np.linalg.pinv(X.T @ X), X), 0, 0.99)
    res = (np.column_stack([r["gcp_u"], r["gcp_v"]]) - pred) / np.sqrt(1 - h)[:, None]  # same for both axes
    east = D(lon + step / np.cos(np.radians(lat)), lat) @ coef - pred
    theta = np.arctan2(east[:, 1], east[:, 0])  # 0 = drawn east points right; pixel v runs down
    return res, theta


def contrasts(res, theta):
    """Per-point p and q of PLAN2 part B, scaled by the map's mean squared residual per axis."""
    s2 = np.mean(res ** 2)  # mean over points and both axes = per-axis variance
    p = (res[:, 0] ** 2 - res[:, 1] ** 2) / (2 * s2)
    q = (2 * res[:, 0] * res[:, 1]) / (2 * s2)
    return p, q


def decompose(p, q, theta):
    """Least squares for a_p and a_c from stacked E[p] = a_p + a_c cos2t, E[q] = a_c sin2t."""
    c, s = np.cos(2 * theta), np.sin(2 * theta)
    X = np.r_[np.column_stack([np.ones_like(c), c]), np.column_stack([np.zeros_like(s), s])]
    y = np.r_[p, q]
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return coef  # a_p, a_c
