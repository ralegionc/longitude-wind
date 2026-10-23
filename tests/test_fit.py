import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fit import R, aeqd_forward, aeqd_inverse, aspect, centroid, dedupe, errors, fit_map  # noqa: E402

RNG = np.random.default_rng(7)


def grid(lon0=10.0, lat0=45.0, span=12.0, n=20):
    lon = lon0 + RNG.uniform(-span / 2, span / 2, n)
    lat = lat0 + RNG.uniform(-span / 3, span / 3, n)
    return lon, lat


def to_pixels(x, y, theta=0.3, scale=2.5, off=(4000, 3000)):
    c, s = np.cos(theta), np.sin(theta)
    return np.column_stack([off[0] + scale * (c * x - s * y), off[1] - scale * (s * x + c * y)])


def test_aeqd_round_trip():
    lon, lat = grid(span=60)
    x, y = aeqd_forward(lon, lat, 10, 45)
    lon2, lat2 = aeqd_inverse(x, y, 10, 45)
    assert np.allclose(lon, lon2, atol=1e-9) and np.allclose(lat, lat2, atol=1e-9)


def test_aeqd_preserves_distance_from_centre():
    x, y = aeqd_forward(np.array([10.0]), np.array([45.0 + 30.0]), 10.0, 45.0)
    assert np.hypot(x, y)[0] == pytest.approx(R * np.radians(30.0), rel=1e-9)


def test_perfect_map_has_no_error_and_no_stretch():
    lon, lat = grid()
    lon0, lat0 = centroid(lon, lat)
    x, y = aeqd_forward(lon, lat, lon0, lat0)
    pix = to_pixels(x, y)
    e_ew, e_ns, _ = errors(pix, lon, lat, lon0, lat0, 2)
    assert np.abs(e_ew).max() < 1e-6 and np.abs(e_ns).max() < 1e-6
    S, ortho = aspect(pix, lon, lat, lon0, lat0)
    assert abs(S) < 1e-6
    assert ortho == pytest.approx(90.0, abs=1e-6)


def test_misplaced_town_shows_up_east():
    lon, lat = grid(n=30)
    lon0, lat0 = centroid(lon, lat)
    x, y = aeqd_forward(lon, lat, lon0, lat0)
    x = x.copy()
    x[5] += 50.0  # the map draws town 5 fifty km too far east
    r = fit_map(to_pixels(x, y), np.column_stack([lon, lat]), outlier_factor=1e9)
    i = int(np.argmax(np.abs(r["e_ew"])))
    assert i == 5
    assert 35 < r["e_ew"][5] < 60
    assert abs(r["e_ns"][5]) < 10


def test_plate_carree_aspect():
    lon, lat = grid(lon0=0, lat0=50, span=8)
    pix = np.column_stack([lon * 100, -lat * 100])  # equal degrees east and north on paper
    lon0, lat0 = centroid(lon, lat)
    S, _ = aspect(pix, lon, lat, lon0, lat0)
    assert S == pytest.approx(-np.log(np.cos(np.radians(50))), abs=0.03)


def anisotropy_trials(ratio, trials=400, n=15):
    out = []
    for _ in range(trials):
        lon, lat = grid(n=n)
        lon0, lat0 = centroid(lon, lat)
        x, y = aeqd_forward(lon, lat, lon0, lat0)
        x = x + RNG.normal(0, 10.0 * ratio, n)
        y = y + RNG.normal(0, 10.0, n)
        r = fit_map(to_pixels(x, y, theta=RNG.uniform(-0.3, 0.3)), np.column_stack([lon, lat]))
        out.append(r["A"])
    return float(np.mean(out))


def test_isotropic_error_gives_zero_anisotropy():
    assert abs(anisotropy_trials(1.0)) < 0.04


def test_east_west_error_recovered():
    assert anisotropy_trials(2.0) == pytest.approx(np.log(2.0), abs=0.06)


def test_dedupe_drops_repeats():
    pix = np.array([[0, 0], [0, 0], [1, 1], [2, 2]], float)
    geo = np.array([[0, 0], [1, 1], [1, 1], [3, 3]], float)
    p2, g2 = dedupe(pix, geo)
    # row 1 repeats a pixel and goes; row 2 repeats only row 1's coordinates, so it stays
    assert len(p2) == 3 and np.allclose(p2[:, 0], [0, 1, 2])


def test_too_few_points_rejected():
    lon, lat = grid(n=5)
    x, y = aeqd_forward(lon, lat, 10, 45)
    assert fit_map(to_pixels(x, y), np.column_stack([lon, lat])) is None


def test_geo_signature_ignores_order_and_catches_copies():
    from build_sample import geo_signature
    geo = np.array([[10.0, 45.0], [11.0, 46.0], [12.5, 44.0]])
    assert geo_signature(geo) == geo_signature(geo[::-1])
    assert geo_signature(geo) != geo_signature(geo + [0.01, 0])
