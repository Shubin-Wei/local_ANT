"""Station table: read whatever form it comes in, and match its names to the names in the record headers.

Accepted tables (recognised automatically):
  * CSV / TSV / whitespace / semicolon separated text, with or without a header line;
  * Excel (.xlsx, first sheet), StationXML (.xml);
  * columns recognised by name (case and spaces ignored; English, Chinese and SAC header names), e.g.
      network: network net knetwk 台网 ...     station: station sta name code kstnm 台站 台站名 站名 点号 ...
      latitude: latitude lat stla 纬度 y ...   longitude: longitude lon long lng stlo 经度 x ...
      elevation: elevation elev alt height stel 高程 海拔 ...
    without a header line: the first column that is not a number is the name; of the numeric columns the one within
    +-90 is the latitude and the other the longitude (when both are within +-90 the table must have a header);
  * names: separate network and station columns, or one column such as "YY1001", "YY.1001", "YY_1001", "1001";
  * coordinates in decimal degrees, or degrees-minutes-seconds (35°07'24.2", 35:07:24.2, 35d07m24.2s), with an
    optional N/S/E/W letter; longitudes 0-360 are brought to -180..180;
  * optional columns: instrument (model), and for nodes that were moved during the survey
      serial (the code in the record headers: serial sn nodeid 仪器号 序列号 ...), start / end (deployment at this
      point: a UTC date and time, or a date alone = that whole day in the time zone of data.day_label; parse_time),
      line (line profile 测线 线号 ...):
    a record is given to the point whose serial is its header code and whose time window holds it
    (assign_by_serial); `line` splits the table into one project per line (data.line).
Names are matched to the record headers by match_names (see there).
"""

import csv
import re

import numpy as np


class StationError(ValueError):
    pass


ALIASES = {
    "network": ["network", "net", "knetwk", "networkcode", "netcode", "台网", "台网代码", "台网名"],
    "station": [
        "station",
        "sta",
        "stat",
        "stname",
        "stationname",
        "stationcode",
        "kstnm",
        "name",
        "code",
        "id",
        "site",
        "point",
        "台站",
        "台站名",
        "台站名称",
        "站名",
        "站点",
        "测点",
        "点号",
        "台号",
        "台站代码",
        "名称",
    ],
    "latitude": ["latitude", "lat", "stla", "纬度", "北纬", "y", "lat(deg)", "latitude(deg)", "lat_deg"],
    "longitude": [
        "longitude",
        "lon",
        "long",
        "lng",
        "stlo",
        "经度",
        "东经",
        "x",
        "lon(deg)",
        "longitude(deg)",
        "lon_deg",
    ],
    "elevation": [
        "elevation",
        "elev",
        "ele",
        "alt",
        "altitude",
        "height",
        "stel",
        "高程",
        "海拔",
        "高度",
        "z",
        "elev(m)",
        "elevation(m)",
        "高程(m)",
        "海拔(m)",
    ],
    "instrument": [
        "instrument",
        "sensor",
        "model",
        "instrumentmodel",
        "sensormodel",
        "seismometer",
        "geophone",
        "仪器",
        "仪器型号",
        "型号",
        "传感器",
        "检波器",
        "地震计",
    ],
    "serial": [
        "serial",
        "serialnumber",
        "serialno",
        "sn",
        "nodeid",
        "instrumentid",
        "仪器号",
        "仪器编号",
        "序列号",
        "设备号",
        "节点号",
    ],
    "start": [
        "start",
        "starttime",
        "begin",
        "begintime",
        "deploy",
        "deploytime",
        "开始时间",
        "起始时间",
        "布设时间",
        "布台时间",
        "布设日期",
        "布台日期",
        "开始日期",
        "起始日期",
        "startdate",
        "deploydate",
    ],
    "end": [
        "end",
        "endtime",
        "stop",
        "stoptime",
        "pickup",
        "pickuptime",
        "结束时间",
        "终止时间",
        "回收时间",
        "撤台时间",
        "收台时间",
        "回收日期",
        "撤台日期",
        "结束日期",
        "enddate",
        "pickupdate",
    ],
    "line": ["line", "profile", "linename", "测线", "线号", "线名", "剖面"],
}


def _norm(s):
    return re.sub(r"[\s_\-\.（）()]+", "", str(s).strip().lower()).replace("（", "").replace("）", "")


_ALIAS_N = {k: {_norm(a) for a in v} for k, v in ALIASES.items()}


_FLOAT = re.compile(r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?")


def parse_coord(s):
    """'35.12', '1.0543e2', '35°07\'24.2"', '35:07:24.2', '35d07m24.2s', '105.43E', 'S33.1' -> float degrees (NaN if
    empty).  A plain number (scientific notation included) is read first; otherwise one hemisphere letter N/S/E/W
    at the start or the end, and degrees [minutes [seconds]] with minutes and seconds below 60."""
    s = str(s).strip()
    if s in ("", "-", "nan", "NaN", "None"):
        return np.nan
    if _FLOAT.fullmatch(s):
        return float(s)
    sign = 1.0
    s = re.sub(r"(\d\s*[mM′']\s*\d+(?:\.\d+)?)\s*[sS](?=\s*[NSEWnsew]?$)", r"\1", s)  # '...24.2s': seconds, not south
    m = re.fullmatch(r"([NSEWnsew])?\s*(.*?)\s*([NSEWnsew])?", s)
    head, body, tail = m.group(1), m.group(2), m.group(3)
    if head and tail:
        raise StationError("cannot read the coordinate %r (two hemisphere letters)" % s)
    hemi = (head or tail or "").upper()
    if hemi in ("S", "W"):
        sign = -1.0
    if _FLOAT.fullmatch(body):
        v = [float(body)]
    else:
        parts = [p for p in re.split(r"\s*(?:[°º:]|[dD](?=\s*\d)|[mM′'](?=\s*\d|\s*$)|[sS″\"]$|\s+)\s*", body) if p]
        parts = [re.sub(r"[sS″\"]$", "", p) for p in parts]
        try:
            v = [float(p) for p in parts]
        except ValueError as exc:
            raise StationError("cannot read the coordinate %r" % s) from exc
    if not 1 <= len(v) <= 3 or any(x < 0 for x in v[1:]) or any(x >= 60 for x in v[1:]):
        raise StationError("cannot read the coordinate %r" % s)
    neg = body.lstrip().startswith("-")
    if neg and hemi:
        raise StationError("cannot read the coordinate %r (minus sign and hemisphere letter)" % s)
    d = abs(v[0]) + (v[1] / 60 if len(v) > 1 else 0) + (v[2] / 3600 if len(v) > 2 else 0)
    return -d if (neg or sign < 0) else d


_PLAIN = re.compile(r"[-+]?\d+(\.\d*)?([eE][-+]?\d+)?")
_DMS = re.compile(
    r"[NSEWnsew]?\s*[-+]?\d+(\.\d+)?\s*[°º:dD]\s*\d+(\.\d+)?\s*[mM′':]?\s*(\d+(\.\d+)?\s*[sS″\"]?)?\s*[NSEWnsew]?"
)
_HEMI = re.compile(r"[NSEWnsew]\s*\d+\.\d+|\d+\.\d+\s*[NSEWnsew]")


def _is_num(s):
    """A coordinate (or elevation) value: a plain number, degrees-minutes-seconds, or a decimal with N/S/E/W.
    Station names such as S01 or E12 are not numbers."""
    s = str(s).strip()
    return bool(_PLAIN.fullmatch(s) or _DMS.fullmatch(s) or _HEMI.fullmatch(s))


def _rows_text(path):
    raw = open(path, "rb").read()
    for enc in ("utf-8-sig", "gbk", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    lines = [l for l in text.splitlines() if l.strip() and not l.lstrip().startswith("#")]
    if not lines:
        raise StationError("%s is empty" % path)
    for d in (",", "\t", ";", "|"):
        if all(d in l for l in lines[: min(5, len(lines))]):
            return [[c.strip() for c in r] for r in csv.reader(lines, delimiter=d)]
    return [l.split() for l in lines]


def _rows_excel(path):
    try:
        import openpyxl
    except ImportError as exc:
        raise StationError(
            "%s: reading Excel files needs openpyxl (pip install openpyxl), or save the sheet as CSV" % path
        ) from exc
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).worksheets[0]
    rows = [["" if c is None else str(c).strip() for c in r] for r in ws.iter_rows(values_only=True)]
    return [r for r in rows if any(r)]


def _assign_header(head):
    cols = {}
    for i, h in enumerate(head):
        n = _norm(h)
        for key, names in _ALIAS_N.items():
            if n in names and key not in cols:
                cols[key] = i
    return cols


def parse_time(value, end=False, date_offset_s=0.0):
    """A start / end time of the station table -> UTC timestamp.  A full date and time is UTC.  A date without a
    time (2024-05-01, 20240501, 240501) is a day of field work in the time zone of the day labels (date_offset_s,
    from data.day_label): the start is the beginning of that day, the end is the end of that day."""
    import obspy

    v = str(value).strip()
    m = re.fullmatch(r"(\d{2}|\d{4})[-/.]?(\d{2})[-/.]?(\d{2})", v)
    if m:
        y = int(m.group(1))
        t = obspy.UTCDateTime(y + 2000 if y < 100 else y, int(m.group(2)), int(m.group(3))).timestamp
        return t - date_offset_s + (86400.0 if end else 0.0)
    return obspy.UTCDateTime(v).timestamp


def read_table(path, date_offset_s=0.0):
    """-> list of dict(network or None, station, lat, lon, elev) and a description of how the file was read."""
    if path.lower().endswith((".xml", ".stationxml")):
        import obspy

        inv = obspy.read_inventory(path)
        out = [
            dict(
                network=n.code,
                station=s.code,
                lat=float(s.latitude),
                lon=float(s.longitude),
                elev=float(s.elevation),
                instrument=None,
                serial=None,
                line=None,
                t_start=None,
                t_end=None,
            )
            for n in inv
            for s in n
        ]
        if not out:
            raise StationError("%s: no station in the StationXML file" % path)
        return out, "StationXML"
    rows = _rows_excel(path) if path.lower().endswith((".xlsx", ".xlsm")) else _rows_text(path)
    head = rows[0]
    cols = _assign_header(head)
    how = []
    if "latitude" in cols and "longitude" in cols and ("station" in cols or "network" in cols):
        body = rows[1:]
        how.append("header: " + ", ".join("%s=%r" % (k, head[i]) for k, i in sorted(cols.items(), key=lambda t: t[1])))
    else:
        # no recognisable header: infer from the values (a header line of unknown words is skipped)
        body = rows[1:] if sum(_is_num(c) for c in head) < 2 else rows
        ncol = min(len(r) for r in body)
        numeric = [j for j in range(ncol) if all(_is_num(r[j]) for r in body)]
        text = [j for j in range(ncol) if j not in numeric]
        if len(numeric) < 2:
            raise StationError(
                "%s: cannot find two coordinate columns; give the table a header line such as "
                "network,station,latitude,longitude,elevation" % path
            )
        if not text:
            # a numeric station code (e.g. "1001 35.12 105.43"): the first column is the name if it is an integer
            if all(re.fullmatch(r"\d+", r[0]) for r in body) and len(numeric) >= 3:
                text, numeric = [0], numeric[1:]
            else:
                raise StationError("%s: no station-name column found" % path)
        vals = {j: np.array([parse_coord(r[j]) for r in body]) for j in numeric}
        lat_c = [j for j in numeric if np.nanmax(np.abs(vals[j])) <= 90]
        first_two = numeric[:2]
        a, b = first_two
        if a in lat_c and b not in lat_c:
            cols.update(latitude=a, longitude=b)
        elif b in lat_c and a not in lat_c:
            cols.update(latitude=b, longitude=a)
        else:
            raise StationError(
                "%s: the two coordinate columns %d and %d are both within +-90 degrees, so latitude and "
                "longitude cannot be told apart; add a header line (latitude, longitude)" % (path, a + 1, b + 1)
            )
        if len(numeric) > 2:
            cols["elevation"] = numeric[2]
        if len(text) >= 2:
            cols.update(network=text[0], station=text[1])
        else:
            cols["station"] = text[0]
        how.append(
            "no header: columns inferred from the values (%s)"
            % ", ".join("%s=column %d" % (k, i + 1) for k, i in sorted(cols.items(), key=lambda t: t[1]))
        )
    out = []
    for k, r in enumerate(body):
        try:
            lat, lon = parse_coord(r[cols["latitude"]]), parse_coord(r[cols["longitude"]])
        except (IndexError, StationError) as exc:
            raise StationError("%s row %d: %s" % (path, k + 2, exc)) from exc
        if not np.isfinite(lat) or not np.isfinite(lon):
            continue
        if abs(lat) > 90 or abs(lon) > 360:
            raise StationError(
                "%s row %d: latitude %g / longitude %g are not geographic degrees (projected coordinates "
                "such as UTM must be converted to latitude/longitude first)" % (path, k + 2, lat, lon)
            )
        el = np.nan
        if "elevation" in cols and cols["elevation"] < len(r):
            try:
                el = parse_coord(r[cols["elevation"]])
            except StationError:
                el = np.nan
        net = r[cols["network"]].strip() if "network" in cols else None
        sta = r[cols["station"]].strip() if "station" in cols else None
        if not sta:
            sta, net = net, None
        opt = lambda k: (r[cols[k]].strip() or None) if k in cols and cols[k] < len(r) else None
        if lon > 180:
            lon -= 360.0
        row = dict(
            network=net or None,
            station=sta,
            lat=lat,
            lon=lon,
            elev=el,
            instrument=opt("instrument"),
            serial=opt("serial"),
            line=opt("line"),
            t_start=None,
            t_end=None,
        )
        for key, col in (("t_start", "start"), ("t_end", "end")):
            if opt(col):
                try:
                    row[key] = parse_time(opt(col), col == "end", date_offset_s)
                except Exception as exc:
                    raise StationError(
                        "%s row %d: %s time %r is not a date (2024-05-01) or a date and time "
                        "(2024-05-01T08:00:00, UTC)" % (path, k + 2, col, opt(col))
                    ) from exc
        if row["serial"] is None and (row["t_start"] is not None or row["t_end"] is not None):
            row["serial"] = sta  # a deployment window without a serial column: the station name is
            # the code of the record headers, used only inside the window
        out.append(row)
    if not out:
        raise StationError("%s: no station rows read" % path)
    return out, "; ".join(how)


def _key(s):
    """Name for comparison: upper case, leading zeros of its number dropped ("A001" = "A1", "0012" = "12")."""
    return re.sub(r"(?<![0-9])0+(?=[0-9])", "", str(s).strip().upper())


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


def split_name(s):
    """'YY.1001' / 'YY_1001' -> ('YY', '1001'); 'YY1001' -> (None, 'YY1001')."""
    m = re.fullmatch(r"([A-Za-z0-9]{1,2})[._]([A-Za-z0-9]+)", s)
    return (m.group(1), m.group(2)) if m else (None, s)


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
