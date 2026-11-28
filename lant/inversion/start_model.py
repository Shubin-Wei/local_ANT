"""1-D models: Rayleigh phase and group velocities of a node model, the average curves, depth nodes, and the 1-D
inversion that gives the starting model of the 3-D inversion.

Node model: Vs at depth nodes, linear in between, the last node a half-space (DSurfTomo's model; Fang et al. 2015);
Vp and density from Vs by Brocher (2005), as DSurfTomo does.  The phase velocities come from DSurfTomo's own
surfdisp96 (lant.inversion.dsurf_disp96), so the starting model, the reference curve from a model, and the 3-D
inversion use one forward calculation.
"""

import numpy as np

MAX_PERIODS = 80  # NP of surfdisp96
MAX_LAYERS = 200  # NL of refineGrid2LayerMdl: every node interval becomes minthk + 1 sublayers
MINTHK = 3  # sublayers per node interval for the depth kernels (DSurfTomo's manual default; its range 2-5)


def phase_velocity(depths, vs, periods, minthk=MINTHK, igr=0):
    """Fundamental Rayleigh phase (igr = 0) or group (igr = 1) velocity (km/s) of the node model at the periods (s),
    in the order given; 0 where surfdisp96 finds no root."""
    from .dsurf_disp96 import column_phase

    z = np.asarray(depths, np.float32)
    v = np.asarray(vs, np.float32)
    t = np.asarray(periods, np.float64)
    if z.ndim != 1 or z.shape != v.shape or z.size < 2 or z[0] != 0 or np.any(np.diff(z) <= 0):
        raise ValueError("a node model needs depths from 0, increasing, and one Vs per depth (at least two nodes)")
    if not (np.all(np.isfinite(v)) and np.all(v > 0) and np.all(t > 0)):
        raise ValueError("Vs and periods must be positive")
    if t.size > MAX_PERIODS:
        raise ValueError("at most %d periods at a time (surfdisp96)" % MAX_PERIODS)
    if (z.size - 1) * (minthk + 1) + 1 > MAX_LAYERS:
        raise ValueError("at most %d depth nodes (DSurfTomo's layer limit)" % ((MAX_LAYERS - 1) // (minthk + 1) + 1))
    order = np.argsort(t, kind="stable")  # surfdisp96 solves the periods in increasing order
    out = np.empty(t.size)
    out[order] = column_phase(z, v, np.float32(minthk), np.ascontiguousarray(t[order]), igr)
    return out


def velocity(depths, vs, periods, group):
    """Phase velocities where group is False, group velocities where it is True (in the order given)."""
    group = np.broadcast_to(np.asarray(group, bool), np.shape(periods))
    t = np.asarray(periods, float)
    out = np.empty(t.size)
    for igr in (0, 1):
        m = group == bool(igr)
        if m.any():
            out[m] = phase_velocity(depths, vs, t[m], igr=igr)
    return out


def read_node_model(path):
    """Two columns depth_km Vs_km/s -> (depths, vs)."""
    a = np.loadtxt(path, ndmin=2, comments="#")
    if a.shape[1] < 2 or a[0, 0] != 0 or np.any(np.diff(a[:, 0]) <= 0) or np.any(a[:, 1] <= 0):
        raise ValueError("%s: two columns depth_km Vs_km/s, depths from 0 and increasing, Vs > 0" % path)
    return a[:, 0], a[:, 1]


def average_curve(period, velocity, min_pairs):
    """Median velocity of each period over the pairs that have one -> rows (period, velocity, pairs)."""
    rows = []
    ok = np.isfinite(velocity)
    for t in np.unique(period[ok]):
        v = velocity[ok & (period == t)]
        if v.size >= min_pairs:
            rows.append((float(t), float(np.median(v)), int(v.size)))
    return np.array(rows).reshape(-1, 3)


def depth_nodes(avg):
    """Depth nodes from 0 in steps of max(a quarter of the shortest wavelength, a quarter of the depth), to the first
    node past half the longest wavelength (the last node: the fixed half-space).

    The shortest wavelength limits the thinnest layer the data can resolve and the vertical resolution of surface
    waves decreases with depth; half the longest wavelength is the usual estimate of the investigation depth
    (Foti et al. 2018, sect. 4.2.3).  Steps are rounded to two significant digits and to whole metres."""
    lam = avg[:, 0] * avg[:, 1]
    dz0, half = float(lam.min()) / 4, float(lam.max()) / 2
    z = [0.0]
    while z[-1] <= half + 1e-9:
        nxt = z[-1] + max(dz0, 0.25 * z[-1])
        z.append(max(round(float("%.2g" % nxt), 3), round(z[-1] + 0.001, 3)))
    return np.array(z)


def curvature_operator(z):
    """Rows R with |R m|^2 ~ the integral of (d^2 m / d zeta^2)^2 over zeta = z / z_max on the actual node spacing,
    each row weighted by the depth interval it stands for, so that the weight does not depend on the nodes."""
    z = np.asarray(z, float)
    Z = float(z[-1] - z[0])
    R = np.zeros((len(z) - 2, len(z)))
    for i in range(1, len(z) - 1):
        h1, h2 = z[i] - z[i - 1], z[i + 1] - z[i]
        row = np.array([2 / (h1 * (h1 + h2)), -2 / (h1 * h2), 2 / (h2 * (h1 + h2))]) * Z**2
        R[i - 1, i - 1 : i + 2] = row * np.sqrt(0.5 * (h1 + h2) / Z)
    return R


SMOOTH = 0.01  # weight of the curvature of the 1-D model (km/s of misfit per unit curvature)
DAMP = 0.05  # damping of each Gauss-Newton step
ITERATIONS = 14
MAX_STEP = 0.15  # largest change of ln Vs per step
CONVERGED = 1e-5  # stop when no step changes ln Vs by more


def invert1d(T, d, z, log=print, group=False):
    """Gauss-Newton on m = ln Vs of the node model at depths z fitting the velocities d at periods T (phase
    velocities, or group velocities where group is True).

    Each step minimises |G dm - r|^2 + SMOOTH^2 |R (m + dm)|^2 + DAMP^2 |dm|^2: the curvature of the model itself is
    smoothed ("jumping", Shaw & Orcutt 1985; Constable et al. 1987), so the result does not drift with the number of
    iterations.  Start: Vs = 1.087 x the velocity of the shortest period (a Poisson half-space has c = U = 0.92 Vs).
    The Jacobian is a forward difference of 1 % in Vs.  -> (vs, predicted velocities, rms per iteration)."""
    z = np.asarray(z, float)
    group = np.broadcast_to(np.asarray(group, bool), np.shape(T))
    nz = len(z)
    vs = np.full(nz, 1.087 * d[np.argmin(T)])
    R = curvature_operator(z)
    hist = []
    for it in range(ITERATIONS):
        g = velocity(z, vs, T, group)
        if np.any(g <= 0):
            raise RuntimeError("1-D forward calculation failed at the periods %s" % T[g <= 0])
        r = d - g
        G = np.zeros((len(T), nz))
        for j in range(nz):
            v2 = vs.copy()
            v2[j] *= 1.01
            G[:, j] = (velocity(z, v2, T, group) - g) / 0.01
        A = np.vstack([G, SMOOTH * R, DAMP * np.eye(nz)])
        b = np.concatenate([r, -SMOOTH * (R @ np.log(vs)), np.zeros(nz)])
        dm = np.clip(np.linalg.lstsq(A, b, rcond=None)[0], -MAX_STEP, MAX_STEP)
        vs = vs * np.exp(dm)
        hist.append(float(np.sqrt(np.mean(r**2))))
        if np.abs(dm).max() < CONVERGED:
            break
    log("  1-D inversion: %d iterations, rms misfit %.4f km/s, Vs %.3f-%.3f km/s" % (len(hist), hist[-1], vs.min(),
                                                                                    vs.max()))  # fmt: skip
    return vs, velocity(z, vs, T, group), hist
