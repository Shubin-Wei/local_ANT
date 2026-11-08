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
    optional N/S/E/W letter.
Names are matched to the record headers by match_names (see there).
"""
import csv
import math
import re

import numpy as np


class StationError(ValueError):
    pass


ALIASES = {
    "network": ["network", "net", "knetwk", "networkcode", "netcode", "台网", "台网代码", "台网名"],
    "station": ["station", "sta", "stat", "stname", "stationname", "stationcode", "kstnm", "name", "code", "id",
                "site", "point", "台站", "台站名", "台站名称", "站名", "站点", "测点", "点号", "台号", "台站代码", "名称"],
    "latitude": ["latitude", "lat", "stla", "纬度", "北纬", "y", "lat(deg)", "latitude(deg)", "lat_deg"],
    "longitude": ["longitude", "lon", "long", "lng", "stlo", "经度", "东经", "x", "lon(deg)", "longitude(deg)", "lon_deg"],
    "elevation": ["elevation", "elev", "ele", "alt", "altitude", "height", "stel", "高程", "海拔", "高度", "z", "elev(m)",
                  "elevation(m)", "高程(m)", "海拔(m)"],
}


def _norm(s):
    return re.sub(r"[\s_\-\.（）()]+", "", str(s).strip().lower()).replace("（", "").replace("）", "")


_ALIAS_N = {k: {_norm(a) for a in v} for k, v in ALIASES.items()}


def parse_coord(s):
    """'35.12', '35°07\\'24.2"', '35:07:24.2', '35d07m24.2s', '105.43E', 'S33.1' -> float degrees (NaN if empty)."""
    s = str(s).strip()
    if s in ("", "-", "nan", "NaN", "None"):
        return np.nan
    sign = 1.0
    m = re.search(r"[NSEWnsew]", s)
    if m:
        if m.group(0).upper() in "SW":
            sign = -1.0
        s = re.sub(r"[NSEWnsew]", "", s).strip()
    try:
        return sign * float(s)
    except ValueError:
        pass
    parts = [p for p in re.split(r"[°º:dDmM′'″\"sS\s]+", s) if p]
    try:
        v = [float(p) for p in parts]
    except ValueError:
        raise StationError("cannot read the coordinate %r" % s)
    if not 1 <= len(v) <= 3:
        raise StationError("cannot read the coordinate %r" % s)
    d = abs(v[0]) + (v[1] / 60 if len(v) > 1 else 0) + (v[2] / 3600 if len(v) > 2 else 0)
    return sign * math.copysign(d, v[0] if v[0] != 0 else 1)


_PLAIN = re.compile(r"[-+]?\d+(\.\d*)?([eE][-+]?\d+)?")
_DMS = re.compile(r"[NSEWnsew]?\s*[-+]?\d+(\.\d+)?\s*[°º:dD]\s*\d+(\.\d+)?\s*[mM′':]?\s*(\d+(\.\d+)?\s*[sS″\"]?)?\s*[NSEWnsew]?")
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
        if all(d in l for l in lines[:min(5, len(lines))]):
            return [[c.strip() for c in r] for r in csv.reader(lines, delimiter=d)]
    return [l.split() for l in lines]


def _rows_excel(path):
    try:
        import openpyxl
    except ImportError:
        raise StationError("%s: reading Excel files needs openpyxl (pip install openpyxl), or save the sheet as CSV" % path)
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


def read_table(path):
    """-> list of dict(network or None, station, lat, lon, elev) and a description of how the file was read."""
    if path.lower().endswith((".xml", ".stationxml")):
        import obspy
        inv = obspy.read_inventory(path)
        out = [dict(network=n.code, station=s.code, lat=float(s.latitude), lon=float(s.longitude), elev=float(s.elevation))
               for n in inv for s in n]
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
            raise StationError("%s: cannot find two coordinate columns; give the table a header line such as "
                               "network,station,latitude,longitude,elevation" % path)
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
            raise StationError("%s: the two coordinate columns %d and %d are both within +-90 degrees, so latitude and "
                               "longitude cannot be told apart; add a header line (latitude, longitude)" % (path, a + 1, b + 1))
        if len(numeric) > 2:
            cols["elevation"] = numeric[2]
        if len(text) >= 2:
            cols.update(network=text[0], station=text[1])
        else:
            cols["station"] = text[0]
        how.append("no header: columns inferred from the values (%s)"
                   % ", ".join("%s=column %d" % (k, i + 1) for k, i in sorted(cols.items(), key=lambda t: t[1])))
    out = []
    for k, r in enumerate(body):
        try:
            lat, lon = parse_coord(r[cols["latitude"]]), parse_coord(r[cols["longitude"]])
        except (IndexError, StationError) as exc:
            raise StationError("%s row %d: %s" % (path, k + 2, exc))
        if not np.isfinite(lat) or not np.isfinite(lon):
            continue
        if abs(lat) > 90 or abs(lon) > 360:
            raise StationError("%s row %d: latitude %g / longitude %g are not geographic degrees (projected coordinates "
                               "such as UTM must be converted to latitude/longitude first)" % (path, k + 2, lat, lon))
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
        out.append(dict(network=net or None, station=sta, lat=lat, lon=lon, elev=el))
    if not out:
        raise StationError("%s: no station rows read" % path)
    return out, "; ".join(how)


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
      4. the header station code is the table name without its leading letters ("1001" for "YY1001").
    A station without a network code anywhere gets the network code of the data set when there is only one, else XX.
    Returns (mapping {(net, sta) of the header: canonical NET.STA}, table {NET.STA: coords}, report lines)."""
    by_sta = {}
    for n, s in data_codes:
        by_sta.setdefault(s.upper(), []).append((n, s))
    joined = {}
    for n, s in data_codes:
        for key in (n + s, n + "." + s, n + "_" + s):
            joined.setdefault(key.upper(), []).append((n, s))
    mapping, coords, rep, unmatched, ambiguous = {}, {}, [], [], []
    nets = {n for n, _ in data_codes if n} | {r["network"] for r in table if r["network"]}
    default_net = nets.pop() if len(nets) == 1 else "XX"        # one network in the whole data set: use it everywhere
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
                if not net:
                    net = m.group(1)
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
            coords.setdefault(name, dict(lon=row["lon"], lat=row["lat"], elev=row["elev"]))
            continue
        if cand[0] in mapping:
            ambiguous.append("%s: the header station %s.%s matches more than one table row" % (name, *cand[0]))
            continue
        mapping[cand[0]] = name
        coords[name] = dict(lon=row["lon"], lat=row["lat"], elev=row["elev"])
        if rule > 1:
            rep.append("%s.%s (header) -> %s (table) [rule %d]" % (cand[0][0], cand[0][1], name, rule))
    return mapping, coords, dict(renamed=rep, table_without_data=unmatched, ambiguous=ambiguous,
                                 data_without_table=sorted("%s.%s" % c for c in data_codes if c not in mapping))
