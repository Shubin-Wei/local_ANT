"""The station table: names and coordinates.

Accepted forms
  * a text table (comma, tab or space separated) with a header line naming the columns: network (optional),
    station, longitude, latitude, elevation (optional).  Short names are accepted: net, sta, lon, lat, elev;
    case does not matter.  Coordinates in decimal degrees (WGS84).
  * a StationXML file (.xml), read with ObsPy.
Every station is called NET.STA from then on.  When the table has no network column, the network is taken from the
record headers (lant.records matches the station code alone and refuses codes that two networks share).
"""

import csv
import os

COLUMNS = {
    "network": ("network", "net"),
    "station": ("station", "sta", "name"),
    "longitude": ("longitude", "lon"),
    "latitude": ("latitude", "lat"),
    "elevation": ("elevation", "elev", "ele"),
}


class StationError(ValueError):
    pass


def _rows(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise StationError("%s is empty" % path)
    head = lines[0]
    if "," in head:
        rows = list(csv.reader(lines))
    elif "\t" in head:
        rows = [ln.split("\t") for ln in lines]
    else:
        rows = [ln.split() for ln in lines]
    return [[c.strip() for c in r] for r in rows]


def read_table(path):
    """-> {name: dict(network, station, lon, lat, elev)}; name is NET.STA, or STA when there is no network column."""
    if not path or not os.path.exists(path):
        raise StationError("station file %s not found" % path)
    if path.lower().endswith(".xml"):
        return _read_stationxml(path)
    rows = _rows(path)
    header = [h.lower() for h in rows[0]]
    col = {}
    for key, names in COLUMNS.items():
        hit = [i for i, h in enumerate(header) if h in names]
        if hit:
            col[key] = hit[0]
    missing = [k for k in ("station", "longitude", "latitude") if k not in col]
    if missing:
        raise StationError(
            "%s: the header line must name the columns %s (found: %s)" % (path, ", ".join(missing), " ".join(rows[0]))
        )
    out = {}
    for n, r in enumerate(rows[1:], start=2):
        try:
            lon, lat = float(r[col["longitude"]]), float(r[col["latitude"]])
            elev = float(r[col["elevation"]]) if "elevation" in col and r[col["elevation"]] else float("nan")
        except (ValueError, IndexError) as exc:
            raise StationError("%s line %d: %s" % (path, n, exc)) from exc
        if not (-180 <= lon <= 360 and -90 <= lat <= 90):
            raise StationError("%s line %d: longitude %g / latitude %g out of range" % (path, n, lon, lat))
        lon = lon - 360 if lon > 180 else lon
        net = r[col["network"]] if "network" in col else ""
        sta = r[col["station"]]
        name = "%s.%s" % (net, sta) if net else sta
        if name in out:
            raise StationError("%s: station %s appears twice" % (path, name))
        out[name] = dict(network=net, station=sta, lon=lon, lat=lat, elev=elev)
    if len(out) < 2:
        raise StationError("%s: fewer than two stations" % path)
    return out


def _read_stationxml(path):
    import obspy

    inv = obspy.read_inventory(path)
    out = {}
    for net in inv:
        for sta in net:
            name = "%s.%s" % (net.code, sta.code)
            out[name] = dict(network=net.code, station=sta.code, lon=sta.longitude, lat=sta.latitude,
                             elev=sta.elevation)  # fmt: skip
    if len(out) < 2:
        raise StationError("%s: fewer than two stations" % path)
    return out


def write_csv(path, stations, extra=None):
    """stations {NET.STA: dict(lon, lat, elev)} -> CSV (extra: {name: {column: value}})."""
    extra = extra or {}
    cols = sorted({k for v in extra.values() for k in v})
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "longitude", "latitude", "elevation"] + cols)
        for name in sorted(stations):
            s = stations[name]
            el = "" if s["elev"] != s["elev"] else "%.1f" % s["elev"]
            w.writerow(
                [name, "%.7f" % s["lon"], "%.7f" % s["lat"], el] + [extra.get(name, {}).get(c, "") for c in cols]
            )


def read_csv(path):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            out[r["name"]] = dict(lon=float(r["longitude"]), lat=float(r["latitude"]),
                                  elev=float(r["elevation"]) if r["elevation"] else float("nan"))  # fmt: skip
    return out
