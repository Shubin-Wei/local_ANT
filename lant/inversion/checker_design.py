"""checks.checkerboards / checks.vertical "auto": three horizontal checkerboards and one depth-varying one, sized from
the data of the built main project (its surfdata.dat: the pairs, periods and phase velocities DSurfTomo inverts).

  cells (half the period of the pattern, km) = 1, 2 and 3 x the median measured wavelength (user's rule 2026-11-13)
    - wavelength: lambda = c T of every measurement; the median is the wavelength the data carry most.  The lateral
      resolution of surface-wave tomography is limited by the predominant wavelengths (Yao et al. 2006, who set the
      correlation length of their phase-velocity maps to about the wavelength), and the ray theory DSurfTomo uses
      (Fang et al. 2015) assumes structure larger than a wavelength.  Rawlinson & Spakman (2016) recommend tests at
      several scales; 1, 2 and 3 wavelengths are the user's choice.
    - a note is written when the median nearest-station distance exceeds a cell: the stations then sample the
      pattern more coarsely than its sign changes (sampling theorem).
  depth-varying test: the middle cell (2 wavelengths), one sign reversal at lambda_med / 3, the depth of the strongest
    sensitivity of the fundamental-mode Rayleigh phase velocity at the median wavelength (rule of thumb, Fang et al.
    2015): the shallow part sensed mainly by the shorter wavelengths has the opposite sign to the deeper part.

Lines: the cell runs along the line, the pattern is constant across it (as before).  Other tests can be listed in
checks.checkerboards / checks.vertical.
"""

import os

import numpy as np

FACTORS = (1, 2, 3)  # cells in median measured wavelengths
VERTICAL_FACTOR = 2


def _two_digits(x):
    return float("%.2g" % x)


def measured(main_dir):
    """(wavelengths km of every travel time, station lon, lat) of the main project's data."""
    from .run import read_data, read_in

    periods = read_in(os.path.join(main_dir, "DSurfTomo.in"))[4]
    d = read_data(os.path.join(main_dir, "surfdata.dat"), periods)
    pts = np.unique(np.r_[d[:, :2], d[:, 2:4]], axis=0)  # (lat, lon) of every station with data
    return d[:, 5] * d[:, 4], pts[:, 1], pts[:, 0]


def design(main_dir, z, array_type):
    """-> (horizontal tests, depth-varying tests, derivation): the "auto" tests and the numbers they come from."""
    from ..auto import geometry

    lam, lon, lat = measured(main_dir)
    g = geometry({i: dict(lon=float(a), lat=float(b)) for i, (a, b) in enumerate(zip(lon, lat))})
    lam_med = float(np.median(lam))
    nn = float(np.nanmedian(g["nn"]))
    cells = [_two_digits(f * lam_med) for f in FACTORS]
    across = 0.0 if array_type == "line" else 1.0
    info = dict(
        median_wavelength_km=round(lam_med, 4),
        wavelength_range_km=[round(float(lam.min()), 4), round(float(lam.max()), 4)],
        median_nearest_station_km=round(nn, 4),
        aperture_km=round(g["aperture"], 3),
        stations_with_data=int(len(lon)),
        cells_km=cells,
        rule="cells = %s x median measured wavelength %.3g km = %s km"
        % (", ".join(map(str, FACTORS)), lam_med, ", ".join("%g" % c for c in cells)),
    )
    notes = []
    if nn > cells[0]:
        notes.append("the median nearest-station distance %.3g km exceeds the smallest cell %g km: the stations "
                     "sample that pattern more coarsely than its sign changes" % (nn, cells[0]))  # fmt: skip
    if g["aperture"] < 2 * cells[-1]:
        notes.append("the aperture %.2f km holds fewer than two cells of %g km: that pattern changes sign at most once "
                     "across the stations" % (g["aperture"], cells[-1]))  # fmt: skip
    if notes:
        info["note"] = "; ".join(notes)
    horizontal = [{"name": "cb_%gkm" % c, "cell_km": [c, c * across]} for c in cells]
    cell = _two_digits(VERTICAL_FACTOR * lam_med)
    zu = np.asarray(z[:-1], float)  # DSurfTomo does not update the last node
    rev = round(lam_med / 3.0, 3)
    if np.any(zu < rev) and np.any(zu > rev):
        vertical = [{"name": "cb_vertical", "cell_km": [cell, cell * across], "reversal_depths_km": [rev]}]
        info["reversal_km"] = rev
        info["reversal_rule"] = "median measured wavelength / 3 = %.3g km (Fang et al. 2015)" % rev
    else:
        vertical = []
        info["vertical_skipped"] = (
            "no depth-varying test: the reversal depth %.3g km (median wavelength / 3) has no updated depth node on "
            "both sides (nodes %s km)" % (rev, ", ".join("%g" % v for v in zu))
        )
    return horizontal, vertical, info
