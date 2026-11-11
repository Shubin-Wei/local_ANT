"""Which station of the table each record of the clean step belongs to, and whether the array geometry suits the
methods: the rows of the station table (of data.line), records without codes named from their files, header codes
matched to the table names (lant/stations.py) or records assigned by serial number and deployment window, the
coordinates of the headers compared with the table, and the geometry checks.  Writes work/01_clean/stations.csv and
stations_from_headers.csv.
"""

import os
from collections import defaultdict

import numpy as np

from . import data as D


def header_coordinates(traces):
    c = defaultdict(list)
    for t in traces:
        if t.get("coord"):
            c[(t["net"], t["sta"])].append(t["coord"])
    return {k: tuple(np.median(np.array(v, float), axis=0)) for k, v in c.items()}


def table_rows(cfg, traces, out, rep):
    """The rows of the station table (of data.line), and records without a station code named from their file
    names (before the channels are chosen, which goes by station)."""
    from . import stations as ST

    d = cfg["data"]
    hc = header_coordinates(traces)
    if hc:
        coords = {
            "%s.%s" % ((k[0] or d["network_default"] or "XX"), k[1]): dict(lon=v[0], lat=v[1], elev=v[2])
            for k, v in hc.items()
        }
        D.write_table(os.path.join(out, "stations_from_headers.csv"), coords)
        rep.append("  coordinates found in the headers of %d stations -> stations_from_headers.csv" % len(hc))
    if not d["station_file"]:
        raise D.DataError(
            "data.station_file is not set%s"
            % (
                "; the headers hold coordinates: check work/01_clean/stations_from_headers.csv and use it as the station file"
                if hc
                else ""
            )
        )
    try:
        rows, how = ST.read_table(d["station_file"], D.label_offset_s(d.get("day_label")))
    except ST.StationError as exc:
        raise D.DataError(str(exc))
    rep.append("  station table %s: %d rows; %s" % (d["station_file"], len(rows), how))
    lines = sorted({r["line"] for r in rows if r.get("line")})
    if len(lines) > 1 and not d["line"]:
        raise D.DataError(
            "the station table holds %d lines (%s): one project per line.  Set data.line to one of them, "
            "or start with `lant run all --raw DIR --stations FILE --out DIR`, which makes one project "
            "per line" % (len(lines), ", ".join(lines[:20]))
        )
    if d["line"]:
        if d["line"] not in lines:
            raise D.DataError(
                "data.line %r is not a line of the station table (%s)"
                % (d["line"], ", ".join(lines) or "no line column")
            )
        rep.append(
            "  line %s: %d of %d table rows" % (d["line"], sum(r.get("line") == d["line"] for r in rows), len(rows))
        )
        rows = [r for r in rows if r.get("line") == d["line"]]
    nfile, unsplit = ST.names_from_files(traces, rows)
    if nfile:
        rep.append("  %d records without a station code in the header: named from their file or directory name" % nfile)
    if unsplit:
        rep.append(
            "  %d files hold several records without station codes at the same times (one trace per channel?): "
            "not used, since the file name cannot tell which record is which station; write one file per "
            "station, or records with station codes (e.g. %s)" % (len(unsplit), unsplit[0])
        )
    return rows


def match(cfg, traces, out, rep, rows):
    """Header codes -> canonical NET.STA of the station table rows; writes work/01_clean/stations.csv."""
    from . import stations as ST

    d = cfg["data"]
    nets = {t["net"] for t in traces if t["net"]} | {r["network"] for r in rows if r["network"]}
    try:
        timed, other, tcoords, rows, trep = ST.assign_by_serial(rows, traces, nets.pop() if len(nets) == 1 else "XX")
    except ST.StationError as exc:
        raise D.DataError(str(exc))
    if trep:
        if not cfg["clean"]["write_archive"]:
            raise D.DataError("a station table with serial numbers needs clean.write_archive: true")
        rep.append(
            "  serial numbers and deployment times of the table: %d points, %d records assigned, %d records of "
            "these serials outside every window (not used)%s"
            % (
                trep["points"],
                trep["records"],
                trep["records_outside_windows"],
                "; points without data: %s" % ", ".join(trep["points_without_data"][:20])
                if trep["points_without_data"]
                else "",
            )
        )
        traces[:] = other
        if not rows:  # every row has a serial: nothing left for the name matching
            used = {n: tcoords[n] for n in {t["station"] for t in timed}}
            if not used:
                raise D.DataError("no record matches a serial number and time window of the station table")
            D.write_table(D.canonical_table(cfg["project"]["output_dir"]), used)
            traces[:] = timed
            return used, dict(
                matched=len(used),
                renamed=0,
                ambiguous=[],
                table_without_data=trep["points_without_data"],
                data_without_table=sorted({"%s.%s" % (t["net"], t["sta"]) for t in other}),
                coordinate_mismatch_km={},
                by_serial=trep,
            )
    codes = {(t["net"], t["sta"]) for t in traces}
    forced = {}
    for k, v in (d["rename"] or {}).items():
        n, s_ = k.split(".", 1)
        forced[(n, s_)] = v
    if d["network_default"]:
        for t in traces:
            if not t["net"]:
                forced.setdefault(("", t["sta"]), "%s.%s" % (d["network_default"], t["sta"]))
    mapping, coords, mrep = ST.match_names(rows, codes - set(forced))
    for code, name in forced.items():
        if code in codes:
            mapping[code] = name
    table_names = set(coords)
    lost = [n for n in set(forced.values()) if n not in table_names]
    if lost:
        raise D.DataError(
            "data.rename / network_default give names that are not in the station table: %s" % ", ".join(sorted(lost))
        )
    if not mapping:
        ex_t = next(iter(table_names))
        ex_d = "%s.%s" % next(iter(codes))
        raise D.DataError(
            "no station of the table could be matched to the record headers (table e.g. %s, headers e.g. %s); "
            "give data.rename" % (ex_t, ex_d)
        )
    used = {name: coords[name] for name in set(mapping.values())}
    used.update({n: tcoords[n] for n in {t["station"] for t in timed}})
    D.write_table(D.canonical_table(cfg["project"]["output_dir"]), used)
    for line in mrep["renamed"][:50]:
        rep.append("    name matched: " + line)
    if len(mrep["renamed"]) > 50:
        rep.append("    ... %d more matched names" % (len(mrep["renamed"]) - 50))
    for key, txt in (
        ("ambiguous", "WARNING ambiguous names (not used): "),
        ("table_without_data", "WARNING table stations without data: "),
        ("data_without_table", "WARNING stations in the data but not in the table (not used): "),
    ):
        if mrep[key]:
            rep.append(
                "  %s%d: %s"
                % (txt, len(mrep[key]), ", ".join(mrep[key][:30]) + (" ..." if len(mrep[key]) > 30 else ""))
            )
    mism = []
    from .signal import geodesic_km

    for code, (lo, la, el) in header_coordinates(traces).items():
        if code in mapping:
            c = used[mapping[code]]
            dkm = float(geodesic_km(lo, la, c["lon"], c["lat"])[0])
            if dkm > 0.05:
                mism.append((mapping[code], round(dkm, 3)))
    if mism:
        rep.append(
            "  WARNING %d stations differ by more than 50 m between the headers and the table (the table is used): %s"
            % (len(mism), ", ".join("%s %.2f km" % m for m in sorted(mism, key=lambda m: -m[1])[:20]))
        )
    for t in traces:
        t["station"] = mapping.get((t["net"], t["sta"]))
    traces += timed
    return used, dict(
        by_serial=trep,
        matched=len(mapping),
        renamed=len(mrep["renamed"]),
        ambiguous=mrep["ambiguous"],
        table_without_data=mrep["table_without_data"],
        data_without_table=mrep["data_without_table"],
        coordinate_mismatch_km=dict(mism),
    )


MAX_LATITUDE = 80.0  # degrees: beyond, the local plane and the regular longitude-latitude grid are too distorted
MIN_SPACING_M = 1.0  # stations closer than this are taken as one point entered twice


def geometry_check(used):
    """The methods use a local plane around the array and a regular longitude-latitude grid.  Refused: an array
    across the date line, beyond MAX_LATITUDE, or with two stations at the same place."""
    from scipy.spatial import cKDTree
    from .signal import local_xy_km

    names = sorted(used)
    lon = np.array([used[n]["lon"] for n in names])
    lat = np.array([used[n]["lat"] for n in names])
    if lon.max() - lon.min() > 180:
        raise D.DataError(
            "the stations lie on both sides of the date line (longitudes %.3f to %.3f): not supported"
            % (lon.min(), lon.max())
        )
    if np.abs(lat).max() > MAX_LATITUDE:
        raise D.DataError(
            "stations beyond %g degrees latitude (%.3f): the local-plane geometry and the "
            "longitude-latitude grid of the inversion are not valid there" % (MAX_LATITUDE, lat[np.abs(lat).argmax()])
        )
    if len(names) > 1:
        x, y = local_xy_km(lon, lat)
        d, j = cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)
        same = [(names[i], names[j[i, 1]]) for i in range(len(names)) if d[i, 1] * 1000 < MIN_SPACING_M and i < j[i, 1]]
        if same:
            raise D.DataError(
                "stations at the same place (closer than %g m): %s.  Correct the coordinates, or leave one "
                "of each pair out (data.exclude_stations)"
                % (MIN_SPACING_M, ", ".join("%s = %s" % p for p in same[:10]))
            )
