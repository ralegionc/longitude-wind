"""Per-map fit: how far is each place from where the map puts it, east-west vs north-south.

See PLAN.md, "Per-map fit". Everything is spherical (R = 6371 km); errors of a
few km from ellipsoidal effects are far below the errors studied here.
"""
import numpy as np

R = 6371.0


def aeqd_forward(lon, lat, lon0, lat0):
    """Spherical azimuthal equidistant projection, km."""
    lam, phi = np.radians(lon) - np.radians(lon0), np.radians(lat)
    phi0 = np.radians(lat0)
    cos_c = np.clip(np.sin(phi0) * np.sin(phi) + np.cos(phi0) * np.cos(phi) * np.cos(lam), -1, 1)
    c = np.arccos(cos_c)
    k = np.where(c < 1e-12, 1.0, c / np.where(c < 1e-12, 1.0, np.sin(c)))
    x = R * k * np.cos(phi) * np.sin(lam)
    y = R * k * (np.cos(phi0) * np.sin(phi) - np.sin(phi0) * np.cos(phi) * np.cos(lam))
    return x, y


def aeqd_inverse(x, y, lon0, lat0):
    phi0 = np.radians(lat0)
    rho = np.hypot(x, y)
    c = rho / R
    safe = np.where(rho < 1e-9, 1.0, rho)
    phi = np.arcsin(np.clip(np.cos(c) * np.sin(phi0) + y * np.sin(c) * np.cos(phi0) / safe, -1, 1))
    lam = np.arctan2(x * np.sin(c), safe * np.cos(phi0) * np.cos(c) - y * np.sin(phi0) * np.sin(c))
    phi = np.where(rho < 1e-9, phi0, phi)
    lam = np.where(rho < 1e-9, 0.0, lam)
    return np.degrees(lam) + lon0, np.degrees(phi)


def centroid(lon, lat):
    """Spherical mean position of a set of points, degrees."""
    lam, phi = np.radians(lon), np.radians(lat)
    v = np.array([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)]).sum(axis=1)
    v /= np.linalg.norm(v)
    return float(np.degrees(np.arctan2(v[1], v[0]))), float(np.degrees(np.arcsin(v[2])))


def design(u, v, order):
    cols = [np.ones_like(u), u, v]
    if order == 2:
        cols += [u * u, u * v, v * v]
    return np.column_stack(cols)


def dedupe(pix, geo):
    """Drop GCPs that repeat a pixel or a coordinate already seen."""
    keep, seen_p, seen_g = [], set(), set()
    for i, (p, g) in enumerate(zip(pix, geo)):
        kp, kg = (round(p[0], 1), round(p[1], 1)), (round(g[0], 6), round(g[1], 6))
        if kp in seen_p or kg in seen_g:
            continue
        seen_p.add(kp)
        seen_g.add(kg)
        keep.append(i)
    return pix[keep], geo[keep]


def max_distance_km(lon, lat):
    lam, phi = np.radians(lon), np.radians(lat)
    d = np.sin(phi[:, None]) * np.sin(phi[None]) + np.cos(phi[:, None]) * np.cos(phi[None]) * np.cos(lam[:, None] - lam[None])
    return float(R * np.arccos(np.clip(d, -1, 1)).max())


def _fit_once(pix, X, Y, order):
    """Least squares pixel -> ground; returns predictions and leverages."""
    mu, sc = pix.mean(axis=0), pix.std() or 1.0
    u, v = ((pix - mu) / sc).T
    D = design(u, v, order)
    coef, *_ = np.linalg.lstsq(D, np.column_stack([X, Y]), rcond=None)
    pred = D @ coef
    H = D @ np.linalg.pinv(D.T @ D) @ D.T
    return pred[:, 0], pred[:, 1], np.clip(np.diag(H), 0, 1)


def order_for(n):
    return 2 if n >= 12 else 1


def errors(pix, lon, lat, lon0, lat0, order):
    """Ground errors (east, north) in km at each GCP and leverages."""
    X, Y = aeqd_forward(lon, lat, lon0, lat0)
    Xh, Yh, h = _fit_once(pix, X, Y, order)
    lon_h, lat_h = aeqd_inverse(Xh, Yh, lon0, lat0)
    dlon = (lon_h - lon + 180.0) % 360.0 - 180.0
    e_ew = R * np.radians(dlon) * np.cos(np.radians(lat))
    e_ns = R * np.radians(lat_h - lat)
    return e_ew, e_ns, h


def aspect(pix, lon, lat, lon0, lat0):
    """S = ln(paper per km east / paper per km north) from an affine ground -> pixel fit,
    plus the angle between the drawn east and north directions (90 = orthogonal)."""
    X, Y = aeqd_forward(lon, lat, lon0, lat0)
    D = np.column_stack([np.ones_like(X), X, Y])
    coef, *_ = np.linalg.lstsq(D, pix, rcond=None)
    east, north = coef[1], coef[2]
    ne, nn = np.linalg.norm(east), np.linalg.norm(north)
    ang = np.degrees(np.arccos(np.clip(abs(east @ north) / (ne * nn), 0, 1)))
    return float(np.log(ne / nn)), float(ang)


def fit_map(pix, geo, outlier_factor=4.0, force_order=None):
    """Full per-map procedure from PLAN.md. pix: (n,2) pixels; geo: (n,2) lon, lat.

    Returns None if the map has too few usable GCPs, else a dict of statistics and
    per-GCP errors (km, leverage-scaled)."""
    pix, geo = np.asarray(pix, float), np.asarray(geo, float)
    if not (np.isfinite(pix).all() and np.isfinite(geo).all()):
        return None
    pix, geo = dedupe(pix, geo)
    if len(pix) < 6:
        return None
    lon, lat = geo[:, 0], geo[:, 1]
    lon0, lat0 = centroid(lon, lat)
    order = force_order or order_for(len(pix))
    e_ew, e_ns, h = errors(pix, lon, lat, lon0, lat0, order)
    tot = np.hypot(e_ew, e_ns)
    med = np.median(tot)
    keep = tot <= outlier_factor * med if med > 0 else np.ones(len(tot), bool)
    n_out = int((~keep).sum())
    if n_out and keep.sum() >= 6:
        pix, lon, lat = pix[keep], lon[keep], lat[keep]
        lon0, lat0 = centroid(lon, lat)
        order = force_order or order_for(len(pix))
        e_ew, e_ns, h = errors(pix, lon, lat, lon0, lat0, order)
    else:
        n_out = 0
    p = 3 if order == 1 else 6
    dof = len(pix) - p
    if dof < 3:
        return None
    ok = h < 0.99
    scale = 1.0 / np.sqrt(1.0 - h[ok])
    e_ew, e_ns = e_ew[ok] * scale, e_ns[ok] * scale
    ss_ew, ss_ns = float(np.sum(e_ew ** 2)), float(np.sum(e_ns ** 2))
    if ss_ew <= 0 or ss_ns <= 0:
        return None
    A = 0.5 * np.log(ss_ew / ss_ns)
    S, ortho = aspect(pix, lon, lat, lon0, lat0)
    ext = max_distance_km(lon, lat)
    rms = float(np.sqrt((ss_ew + ss_ns) / len(e_ew)))
    return {
        "n": int(len(pix)), "n_out": n_out, "order": order, "dof": dof,
        "lon0": lon0, "lat0": lat0, "extent_km": ext,
        "A": float(A), "S": S, "A_paper": float(A + S), "ortho_deg": ortho,
        "rms_km": rms, "relerr": rms / ext if ext > 0 else np.nan,
        "rms_ew_km": float(np.sqrt(ss_ew / len(e_ew))), "rms_ns_km": float(np.sqrt(ss_ns / len(e_ns))),
        "gcp_lon": lon[ok], "gcp_lat": lat[ok], "e_ew": e_ew, "e_ns": e_ns,
        "gcp_u": pix[ok, 0], "gcp_v": pix[ok, 1],
    }
