"""Station table: read whatever form it comes in, and normalise its names (split_name, _key) for matching.

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
    (lant.station_match.assign_by_serial); `line` splits the table into one project per line (data.line).
Names are matched to the record headers by lant.station_match.match_names (see there).
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


def split_name(s):
    """'YY.1001' / 'YY_1001' -> ('YY', '1001'); 'YY1001' -> (None, 'YY1001')."""
    m = re.fullmatch(r"([A-Za-z0-9]{1,2})[._]([A-Za-z0-9]+)", s)
    return (m.group(1), m.group(2)) if m else (None, s)
