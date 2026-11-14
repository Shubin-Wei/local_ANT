"""Read a DSurfTomo project exactly as DSurfTomo's main.f90 does.

Mirrors DSurfTomo a651faf (MIT): list-directed reads of DSurfTomo.in and MOD,
the measurement grouping of surfdata.dat (a new period block restarts the source
counter, so a period that reappears later overwrites its earlier sources), the
single-precision conversion to colatitude/longitude radians, and the
single-precision haversine of delsph.f90 for observed times. Every real value
the Fortran stores as REAL is kept as float32 here so the compatible forward
core starts from identical numbers.
"""

import ctypes
import ctypes.util
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

F32 = np.float32
# REAL parameter pi of main.f90/CalSurfG.f90 (single precision).
PI32 = F32(3.1415926535898)


class _Records:
    """Fortran list-directed reads: each READ starts a new record and may span lines."""

    def __init__(self, lines):
        self.lines = list(lines)
        self.at = 0

    def skip(self):
        self.at += 1

    def read(self, count):
        values = []
        while len(values) < count:
            if self.at >= len(self.lines):
                raise ValueError("input ended while %d values were still expected" % (count - len(values)))
            values.extend(self.lines[self.at].replace(",", " ").split())
            self.at += 1
        return values[:count]


@dataclass
class Control:
    datafile: str
    nx: int
    ny: int
    nz: int
    goxd: np.float32
    gozd: np.float32
    dvxd: np.float32
    dvzd: np.float32
    nsrc: int
    minthk: np.float32
    spfra: float
    periods_rc: np.ndarray  # real*8 periods of Rayleigh phase velocity
    kmax_rg: int
    kmax_lc: int
    kmax_lg: int
    # main.f90 inversion controls: REAL weight0, damp, Minvel, Maxvel, noiselevel, threshold0; integers.
    weight: np.float32 = F32(0)
    damp: np.float32 = F32(0)
    minvel: np.float32 = F32(0)
    maxvel: np.float32 = F32(0)
    maxiter: int = 0
    ifsyn: int = 0
    noiselevel: np.float32 = F32(0)
    threshold: np.float32 = F32(0)

    @property
    def kmax(self):
        return len(self.periods_rc) + self.kmax_rg + self.kmax_lc + self.kmax_lg


def read_control(path):
    r = _Records(Path(path).read_text().splitlines())
    for _ in range(3):
        r.skip()
    datafile = r.read(1)[0]
    nx, ny, nz = (int(v) for v in r.read(3))
    goxd, gozd = (F32(v) for v in r.read(2))
    dvxd, dvzd = (F32(v) for v in r.read(2))
    nsrc = int(r.read(1)[0])
    weight, damp = (F32(v) for v in r.read(2))
    minthk = F32(r.read(1)[0])
    minvel, maxvel = (F32(v) for v in r.read(2))
    maxiter = int(r.read(1)[0])
    spfra = float(r.read(1)[0])
    counts, periods = [], []
    for _ in range(4):
        n = int(r.read(1)[0])
        counts.append(n)
        periods.append(np.array(r.read(n), float) if n > 0 else np.zeros(0))
    ifsyn = int(r.read(1)[0])
    noiselevel, threshold = F32(r.read(1)[0]), F32(r.read(1)[0])
    if min(nx, ny) < 4 or nz < 2 or nsrc < 1:
        raise ValueError("DSurfTomo grid needs nx, ny >= 4, nz >= 2 and at least one source")
    if counts[0] < 1:
        raise ValueError("no Rayleigh phase periods in %s" % path)
    return Control(
        datafile, nx, ny, nz, goxd, gozd, dvxd, dvzd, nsrc, minthk, spfra, periods[0], *counts[1:],
        weight, damp, minvel, maxvel, maxiter, ifsyn, noiselevel, threshold,
    )  # fmt: skip


def colatitude_radians(lat):
    """(90.0 - lat) * pi / 180.0 in single precision, as main.f90."""
    return (F32(90.0) - F32(lat)) * PI32 / F32(180.0)


def longitude_radians(lon):
    return F32(lon) * PI32 / F32(180.0)


def delsph(flat1, flon1, flat2, flon2):
    """delsph.f90: haversine in REAL arithmetic, R = 6371 km, colatitude inputs."""
    # Python floats would promote to double; every step is rounded to float32 like the Fortran.
    flat1, flon1, flat2, flon2 = F32(flat1), F32(flon1), F32(flat2), F32(flon2)
    pi = F32(3.1415926535898)
    dlat = flat2 - flat1
    dlon = flon2 - flon1
    lat1 = pi / F32(2) - flat1
    lat2 = pi / F32(2) - flat2
    s1 = _sinf(dlat / F32(2))
    s2 = _sinf(dlon / F32(2))
    a = s1 * s1 + s2 * s2 * _cosf(lat1) * _cosf(lat2)
    c = F32(2) * _atan2f(np.sqrt(a), np.sqrt(F32(1) - a))
    return F32(6371.0) * c


def _single(name, nargs):
    """The C library's single-precision function, which gfortran calls for REAL sin/cos/atan2."""
    path = ctypes.util.find_library("m")
    function = getattr(ctypes.CDLL(path), name, None) if path else None
    if function is None:  # double precision rounded to REAL; may differ in the last bit
        double = getattr(math, name[:-1])
        return lambda *x: F32(double(*map(float, x)))
    function.restype = ctypes.c_float
    function.argtypes = [ctypes.c_float] * nargs
    return lambda *x: F32(function(*map(float, x)))


_sinf, _cosf, _atan2f = _single("sinf", 1), _single("cosf", 1), _single("atan2f", 2)


@dataclass
class Measurements:
    """main.f90 arrays (0-based here) plus the order in which CalSurfG visits rows."""

    scx: np.ndarray  # (kmax, nsrc) float32 source colatitude, rad
    scz: np.ndarray  # (kmax, nsrc) float32 source longitude, rad
    rcx: np.ndarray  # (kmax, nsrc, nrc) float32 receiver colatitude, rad
    rcz: np.ndarray
    period: np.ndarray  # (kmax, nsrc) period index of each source line (1-based)
    nrc1: np.ndarray  # (kmax, nsrc) receivers per source
    nsrc1: np.ndarray  # (kmax,) sources per period block
    file_row: np.ndarray  # (kmax, nsrc, nrc) index of the surfdata.dat row that set the entry
    observed: np.ndarray  # (rows in file,) float32 delsph distance / velocity, file order
    distance: np.ndarray  # (rows in file,) float32 km

    def calsurfg_rows(self):
        """(knum, source, receiver) of every predicted row, in CalSurfG order (1-based knum)."""
        order = []
        for k in range(self.nsrc1.size):
            for s in range(self.nsrc1[k]):
                order.extend((k + 1, s, r) for r in range(self.nrc1[k, s]))
        return np.array(order, dtype=np.int64).reshape(-1, 3)


def read_measurements(path, control):
    """Parse surfdata.dat with the grouping and overwrite rules of main.f90."""
    kmax, ns = control.kmax, control.nsrc
    nrc = ns
    scx, scz = np.zeros((kmax, ns), F32), np.zeros((kmax, ns), F32)
    rcx, rcz = np.zeros((kmax, ns, nrc), F32), np.zeros((kmax, ns, nrc), F32)
    period = np.zeros((kmax, ns), np.int64)
    nrc1 = np.zeros((kmax, ns), np.int64)
    nsrc1 = np.zeros(kmax, np.int64)
    file_row = np.full((kmax, ns, nrc), -1, np.int64)
    observed, distance = [], []
    krc = len(control.periods_rc)
    knum, knumo, istep, istep1 = 0, 12345, 0, 0
    source = None
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if line.startswith("#"):
            fields = line.split()
            lat, lon, p, wave, kind = F32(fields[1]), F32(fields[2]), int(fields[3]), int(fields[4]), int(fields[5])
            if (wave, kind) != (2, 0):
                raise NotImplementedError("line %d: only Rayleigh phase velocity (2 0) is supported" % number)
            if not 1 <= p <= krc:
                raise ValueError("line %d: period index %d outside 1..%d" % (number, p, krc))
            knum = p
            if knum != knumo:
                istep = 0
            istep += 1
            istep1 = 0
            if istep > ns:
                raise ValueError("line %d: more sources in a period block than nsrc=%d" % (number, ns))
            source = (colatitude_radians(lat), longitude_radians(lon))
            scx[knum - 1, istep - 1], scz[knum - 1, istep - 1] = source
            period[knum - 1, istep - 1] = p
            nsrc1[knum - 1] = istep
            knumo = knum
        elif line.strip():
            if source is None:
                raise ValueError("line %d: receiver before any source line" % number)
            fields = line.split()
            lat, lon, velocity = F32(fields[0]), F32(fields[1]), F32(fields[2])
            istep1 += 1
            if istep1 > nrc:
                raise ValueError("line %d: more receivers than nsrc=%d" % (number, nrc))
            x, z = colatitude_radians(lat), longitude_radians(lon)
            rcx[knum - 1, istep - 1, istep1 - 1], rcz[knum - 1, istep - 1, istep1 - 1] = x, z
            d = delsph(source[0], source[1], x, z)
            file_row[knum - 1, istep - 1, istep1 - 1] = len(observed)
            distance.append(d)
            observed.append(d / velocity)
            nrc1[knum - 1, istep - 1] = istep1
        else:
            # main.f90 reads a blank line as a receiver record and stops with a read error.
            raise ValueError("line %d: blank line in measurement file" % number)
    return Measurements(
        scx, scz, rcx, rcz, period, nrc1, nsrc1, file_row, np.array(observed, F32), np.array(distance, F32)
    )


def _velocities(r, nx, ny, nz):
    vs = np.empty((nz, ny, nx), F32)
    for k in range(nz):
        for j in range(ny):
            vs[k, j] = np.array(r.read(nx), F32)
    return vs


def read_model(path, nx, ny, nz):
    """MOD: nz depths, then for each depth and longitude one READ of nx Vs values."""
    r = _Records(Path(path).read_text().splitlines())
    depths = np.array(r.read(nz), F32)
    return depths, _velocities(r, nx, ny, nz)


def read_true_model(path, nx, ny, nz):
    """MOD.true of the synthetic test: like MOD but without the depth record (main.f90)."""
    return _velocities(_Records(Path(path).read_text().splitlines()), nx, ny, nz)


@dataclass
class Project:
    control: Control
    measurements: Measurements
    depths: np.ndarray
    vs: np.ndarray


def read_project(directory, model=None):
    """DSurfTomo.in, its data file and a MOD-format model (default: the project's MOD)."""
    directory = Path(directory)
    control = read_control(directory / "DSurfTomo.in")
    measurements = read_measurements(directory / control.datafile, control)
    depths, vs = read_model(Path(model) if model else directory / "MOD", control.nx, control.ny, control.nz)
    return Project(control, measurements, depths, vs)
