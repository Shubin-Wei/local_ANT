"""Station tables, coordinates and name matching (lant.stations)."""
import pytest


def test_station_tables_and_names(tmp_path):
    """Station tables in several layouts (CSV, plain text, Chinese headers, degree-minute-second) and name matching
    with and without network codes."""
    from lant import stations as ST
    tables = {
        "a.csv": "network,station,channel,latitude,longitude,elevation\nYY,1001,HHZ,35.1234,105.4321,306\nYY,1002,HHZ,35.2,105.5,483\n",
        "b.txt": "YY1001 105.4321 35.1234 306\nYY1002 105.5 35.2 483\n",
        "c.txt": "台站名\t经度\t纬度\nYY.1001\t105°25'55.56\"E\t35°7'24.24\"N\nYY.1002\t105.5\t35.2\n",
        "d.txt": "S01 100.04 30.07 335\nS02 100.01 30.08 517\n",
    }
    for name, text in tables.items():
        f = tmp_path / name
        f.write_text(text)
        rows, how = ST.read_table(str(f))
        assert len(rows) == 2 and abs(rows[0]["lat"] - (30.07 if name == "d.txt" else 35.1234)) < 1e-3, (name, rows, how)
    rows, _ = ST.read_table(str(tmp_path / "b.txt"))
    mapping, _, rep = ST.match_names(rows, {("YY", "1001"), ("", "1002")})
    assert mapping == {("YY", "1001"): "YY.1001", ("", "1002"): "YY.1002"}
    mapping, _, rep = ST.match_names([dict(network=None, station="A01", lat=1, lon=2, elev=0)], {("X", "A01"), ("Y", "A01")})
    assert not mapping and rep["ambiguous"]
    assert not ST._is_num("S01") and ST._is_num("105.43E") and abs(ST.parse_coord("40:30:00") - 40.5) < 1e-9


def test_coordinates_channels_ambiguous_names():
    """Coordinate notations; the vertical channel is chosen (a horizontal one is refused, an unknown orientation is
    left to confirm); one table name for two points is ambiguous."""
    from lant import stations as ST, data as D
    assert abs(ST.parse_coord("1.054321e2") - 105.4321) < 1e-9
    assert abs(ST.parse_coord("35d07m24.2s") - 35.123389) < 1e-5 and ST.parse_coord("35d07m24.2sS") < 0
    assert ST.parse_coord("S33.1") == -33.1 and abs(ST.parse_coord("105°25'55.56\"E") - 105.4321) < 1e-6
    rows, notes = D.choose_channel([dict(net="X", sta="A", channel="HHN", cmpinc=90.0)], "auto")
    assert not rows and "horizontal" in notes[0]
    rows, notes = D.choose_channel([dict(net="X", sta="B", channel="DP3", cmpinc=None)], "auto")
    assert rows and "confirm" in notes[0]
    tab = [dict(network=None, station="A01", lat=1, lon=1, elev=0), dict(network=None, station="A01", lat=2, lon=2, elev=0)]
    m, c, r = ST.match_names(tab, {("ZZ", "A01")})
    assert not m and not c and r["ambiguous"]


def test_moved_nodes_names_geometry(tmp_path):
    """Station table with serial numbers and deployment times: a record is given to the point of its window and
    cut to it; lines; names with other zero padding; records named from the file name; geometry checks."""
    from lant import stations as ST, station_match as SM, data as D
    (tmp_path / "t.csv").write_text(
        "point,serial,start,end,line,lat,lon\n"
        "P1,N7,2024-01-01T00:00:00,2024-01-02T00:00:00,L1,30.0,100.0\n"
        "P2,N7,2024-01-02T06:00:00,2024-01-03T00:00:00,L1,30.0,100.01\n"
        "Q1,N8,,,L2,30.1,260.0\n")
    rows, _ = ST.read_table(str(tmp_path / "t.csv"))
    assert rows[2]["lon"] == -100.0 and rows[0]["line"] == "L1" and rows[0]["t_end"] == 1704153600.0
    t0 = 1704067200.0
    tr = lambda sta, a, b, f="x": dict(net="", sta=sta, start=t0 + a * 3600, end=t0 + b * 3600, file=f)
    used, rest, coords, left, rep = ST.assign_by_serial(rows, [tr("N7", 12, 40), tr("N7", 25, 29), tr("N9", 0, 5)])
    assert [(u["station"], (u["start"] - t0) / 3600, (u["end"] - t0) / 3600) for u in used] == [
        ("XX.P1", 12.0, 24.0), ("XX.P2", 30.0, 40.0)]
    assert rep["records_outside_windows"] == 1 and len(rest) == 1 and not left and rep["points_without_data"] == ["XX.Q1"]
    tab = [dict(network=None, station="A001", lat=1, lon=1, elev=0), dict(network=None, station="B12", lat=2, lon=2, elev=0)]
    m, c, r = ST.match_names(tab, {("", "A1"), ("", "B012")})
    assert m == {("", "A1"): "XX.A1", ("", "B012"): "XX.B012"}
    rec = lambda f, a, b: dict(sta="", file=f, start=a, end=b)
    recs = [rec("/d/line1/A001_Z.sgy", 0, 10), rec("/d/line1/A001_Z.sgy", 10, 20), rec("/d/B12/rec.sgy", 0, 10),
            rec("/d/x/A001_B12.sgy", 0, 10), rec("/d/B12/shot.sgy", 0, 10), rec("/d/B12/shot.sgy", 0, 10)]
    n, unsplit = ST.names_from_files(recs, tab)        # consecutive traces of one file: named; overlapping: not
    assert n == 3 and [x["sta"] for x in recs] == ["A001", "A001", "B12", "", "", ""] and unsplit == ["/d/B12/shot.sgy"]
    ok = {"A": dict(lon=100.0, lat=30.0), "B": dict(lon=100.01, lat=30.0)}
    SM.geometry_check(ok)
    for bad, msg in (({"A": dict(lon=179.9, lat=0.0), "B": dict(lon=-179.9, lat=0.0)}, "date line"),
                     ({"A": dict(lon=10.0, lat=85.0), "B": dict(lon=10.1, lat=85.0)}, "latitude"),
                     ({"A": dict(lon=100.0, lat=30.0), "B": dict(lon=100.0, lat=30.000001)}, "same place")):
        with pytest.raises(D.DataError, match=msg):
            SM.geometry_check(bad)
