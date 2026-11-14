"""Which station of the table each record of the clean step belongs to, and whether the array geometry suits the
methods: the rows of the station table (of data.line), records without codes named from their files, header codes
matched to the table names (match_names; the names are normalised by lant.stations) or records assigned by serial number and deployment window, the
coordinates of the headers compared with the table, and the geometry checks.  Writes work/01_clean/stations.csv and
stations_from_headers.csv.
"""

import os
import re
from collections import defaultdict

import numpy as np

from . import data as D
from .stations import StationError, _key, split_name


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
        raise D.DataError(str(exc)) from exc
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
    nfile, unsplit = names_from_files(traces, rows)
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
        timed, other, tcoords, rows, trep = assign_by_serial(rows, traces, nets.pop() if len(nets) == 1 else "XX")
    except ST.StationError as exc:
        raise D.DataError(str(exc)) from exc
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
    mapping, coords, mrep = match_names(rows, codes - set(forced))
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


def could_match(cfg):
    """-> predicate(net, sta): False only for a header code that no rule of match_names, and no data.rename, can give
    a station of the table (the clean step need not read such files in full); None when every file must be read
    (the table cannot be read here, or data.network_default renames).

    A superset of assign_by_serial (the header code is a serial of the table, whatever the time) and of match_names:
    rule 1 is contained in 2 (same station code), 3 (network and station joined), 4 (the digits of the table name),
    5 (leading zeros ignored).  A record without a station code always counts."""
    from . import stations as ST

    d = cfg["data"]
    if not d["station_file"] or d["network_default"]:
        return None
    try:
        rows, _ = ST.read_table(d["station_file"], D.label_offset_s(d.get("day_label")))
    except (ST.StationError, OSError, ValueError):
        return None
    if d["line"]:
        rows = [r for r in rows if str(r.get("line")) == str(d["line"])]
    if not rows:
        return None
    codes, padded, joined = set(), set(), set()
    for r in rows:
        if r.get("serial"):
            padded.add(_key(r["serial"]))
            continue
        net, sta = r["network"], r["station"]
        if net is None:
            n2, s2 = split_name(sta)
            if n2:
                net, sta = n2, s2
        codes.add(sta.upper())
        padded.add(_key(sta))
        joined.add(((net or "") + sta).upper())
        m = re.fullmatch(r"([A-Za-z]+)(\d+)", sta)
        if m:
            codes.add(m.group(2))
    forced = {tuple(k.split(".", 1)) for k in (d["rename"] or {})}

    def relevant(net, sta):
        if not sta or (net, sta) in forced:
            return True
        return (
            sta.upper() in codes
            or _key(sta) in padded
            or any((net + sep + sta).upper() in joined for sep in ("", ".", "_"))
        )

    return relevant


def names_from_files(traces, table):
    """Records without a station code in the header (SEG-Y, some SAC files): the code is taken from the file name
    or the directory name when exactly one of their parts (split at characters other than letters and digits) is a
    station name or serial of the table (leading zeros ignored).  A file with several records without a code that
    overlap in time under one channel code (one trace per channel, as in many SEG-2 and SEG-Y files, where the
    channel code is empty) cannot be split by its name: its records are left unnamed.  Records of one file following each other in time (a long record cut into traces) are named.
    -> (number of records named, files left unnamed because their records overlap)."""
    import os
    from collections import defaultdict

    bad = set()
    per_file = defaultdict(list)  # records of one channel code at the same time: not separable
    for t in traces:
        if not t["sta"]:
            per_file[(t["file"], t.get("channel", ""))].append((t["start"], t["end"]))
    for (f, _), iv in per_file.items():
        iv.sort()
        if any(b[0] < a[1] - 1e-6 for a, b in zip(iv, iv[1:])):
            bad.add(f)
    known = {}
    for r in table:
        for nm in (r["station"], r.get("serial")):
            if nm:
                known.setdefault(_key(split_name(nm)[1]), nm)
    n = 0
    for t in traces:
        if t["sta"] or t["file"] in bad:
            continue
        stem = os.path.splitext(os.path.basename(t["file"]))[0]
        parts = (
            re.split(r"[^A-Za-z0-9]+", stem)
            + [stem]
            + re.split(r"[^A-Za-z0-9]+", os.path.basename(os.path.dirname(t["file"])))
        )
        hit = {known[_key(p)] for p in parts if p and _key(p) in known}
        if len(hit) == 1:
            t["sta_header"] = t["sta"]  # the code in the file, to pick the record when it is read
            t["sta"] = split_name(hit.pop())[1]
            n += 1
    return n, sorted(bad)


def assign_by_serial(table, traces, default_net="XX"):
    """Nodes moved during the survey: table rows with a serial (and optionally start / end).  A record whose header
    station code is the serial is given to the point whose window [start, end] it overlaps and is cut to the window
    (key "clip"); a record over two windows is used for both points.  The rows and records used here take no part
    in match_names.  -> (records [with "station" and "clip"], other records, coords {NET.STA: ...}, other table
    rows, report dict)."""
    timed = [r for r in table if r.get("serial")]
    rest_rows = [r for r in table if not r.get("serial")]
    if not timed:
        return [], traces, {}, table, {}
    by_serial, coords = {}, {}
    for r in timed:
        net, sta = (r["network"], r["station"]) if r["network"] else split_name(r["station"])
        name = "%s.%s" % (net or default_net, sta)
        if name in coords:
            raise StationError("the station table has two rows for point %s" % name)
        coords[name] = dict(
            lon=r["lon"], lat=r["lat"], elev=r["elev"], instrument=r.get("instrument"), line=r.get("line")
        )
        by_serial.setdefault(_key(r["serial"]), []).append(
            (
                r["t_start"] if r["t_start"] is not None else -np.inf,
                r["t_end"] if r["t_end"] is not None else np.inf,
                name,
            )
        )
    for key, wins in by_serial.items():
        wins.sort()
        for a, b in zip(wins[:-1], wins[1:]):
            if b[0] < a[1]:
                raise StationError(
                    "serial %s is at %s and at %s at the same time (windows overlap): correct the "
                    "start / end times of the station table" % (key, a[2], b[2])
                )
    used, rest, outside, hit = [], [], 0, set()
    for t in traces:
        wins = by_serial.get(_key(t["sta"]))
        if wins is None:
            rest.append(t)
            continue
        got = False
        for a, b, name in wins:
            lo, hi = max(a, t["start"]), min(b, t["end"])
            if hi > lo:
                used.append(dict(t, station=name, start=lo, end=hi, clip=[max(a, -1e18), min(b, 1e18)]))
                hit.add(name)
                got = True
        outside += not got
    return (
        used,
        rest,
        coords,
        rest_rows,
        dict(
            points=len(coords),
            records=len(used),
            records_outside_windows=int(outside),
            points_without_data=sorted(set(coords) - hit),
        ),
    )


def match_names(table, data_codes):
    """Match the table rows to the stations found in the headers.

    data_codes: set of (network, station) of the headers (network may be '').
    Rules, in order (a rule is used only when it gives one candidate):
      1. same network and station;  2. same station code (the table, or the header, has no or another network);
      3. the table name is network + station of the header joined ("YY1001", "YY.1001", "YY_1001");
      4. the header station code is the table name without its leading letters ("1001" for "YY1001");
      5. as 2, ignoring leading zeros of the numbers ("A001" = "A1", "0012" = "12").
    A station without a network code anywhere gets the network code of the data set when there is only one, else XX.
    Returns (mapping {(net, sta) of the header: canonical NET.STA}, table {NET.STA: coords}, report lines)."""
    by_sta = {}
    for n, s in data_codes:
        by_sta.setdefault(s.upper(), []).append((n, s))
    joined = {}
    for n, s in data_codes:
        for key in (n + s, n + "." + s, n + "_" + s):
            joined.setdefault(key.upper(), []).append((n, s))
    by_pad = {}
    for n, s in data_codes:
        by_pad.setdefault(_key(s), []).append((n, s))
    mapping, coords, rep, unmatched, ambiguous = {}, {}, [], [], []
    rows_of = {}
    nets = {n for n, _ in data_codes if n} | {r["network"] for r in table if r["network"]}
    default_net = nets.pop() if len(nets) == 1 else "XX"  # one network in the whole data set: use it everywhere
    for row in table:
        net, sta = row["network"], row["station"]
        if net is None:
            n2, s2 = split_name(sta)
            if n2:
                net, sta = n2, s2
        cand = []
        if net is not None:
            cand = [c for c in data_codes if c[0].upper() == net.upper() and c[1].upper() == sta.upper()]
        rule = 1
        if not cand:
            rule, cand = 2, by_sta.get(sta.upper(), [])
        if not cand:
            rule, cand = 3, joined.get(((net or "") + sta).upper(), [])
        if not cand:
            m = re.fullmatch(r"([A-Za-z]+)(\d+)", sta)
            if m:
                rule, cand = 4, by_sta.get(m.group(2), [])
                if cand and not net:
                    net = m.group(1)
        if not cand:
            rule, cand = 5, by_pad.get(_key(sta), [])
        if len(cand) > 1:
            ambiguous.append("%s%s: %s" % (net + "." if net else "", sta, ", ".join("%s.%s" % c for c in cand)))
            continue
        if not cand:
            name = "%s.%s" % (net or default_net, sta)
        elif rule == 1:
            name = "%s.%s" % (net, sta)
        elif rule == 2:
            name = "%s.%s" % (net or cand[0][0] or default_net, cand[0][1])
        else:
            name = "%s.%s" % (cand[0][0] or net or default_net, cand[0][1])
        if not cand:
            unmatched.append(name)
            coords.setdefault(
                name,
                dict(
                    lon=row["lon"],
                    lat=row["lat"],
                    elev=row["elev"],
                    instrument=row.get("instrument"),
                    line=row.get("line"),
                ),
            )
            continue
        rows_of.setdefault(cand[0], []).append((name, row, rule))
        if cand[0] in mapping:
            continue
        mapping[cand[0]] = name
        coords[name] = dict(
            lon=row["lon"], lat=row["lat"], elev=row["elev"], instrument=row.get("instrument"), line=row.get("line")
        )
        if rule > 1:
            rep.append("%s.%s (header) -> %s (table) [rule %d]" % (cand[0][0], cand[0][1], name, rule))
    for code, lst in rows_of.items():  # one header station claimed by several table rows: none is used
        if len(lst) > 1:
            name = mapping.pop(code)
            coords.pop(name, None)
            rep[:] = [x for x in rep if not x.startswith("%s.%s (header)" % code)]
            ambiguous.append(
                "header station %s.%s matches %d table rows (%s): not used until the table is fixed or "
                "data.rename says which"
                % (
                    code[0],
                    code[1],
                    len(lst),
                    "; ".join("%s at %.5f %.5f" % (n, r["lon"], r["lat"]) for n, r, _ in lst),
                )
            )
    return (
        mapping,
        coords,
        dict(
            renamed=rep,
            table_without_data=unmatched,
            ambiguous=ambiguous,
            data_without_table=sorted("%s.%s" % c for c in data_codes if c not in mapping),
        ),
    )
