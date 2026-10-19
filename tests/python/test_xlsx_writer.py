"""bin/xlsx_writer.py: a stdlib-only .xlsx writer. These tests parse every part of the package with
xml.etree and check the invariants Excel is strict about (style ids, dimension, sheet names, escaping).
The independent reader check (openpyxl) is described in docs/wiring/xlsx_output.md."""

import re
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest
import xlsx_writer as xw

REPO = Path(__file__).resolve().parents[2]
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
M = "{" + NS["m"] + "}"
HOSTILE = [
    "=cmd|' /C calc'!A0",
    "+1+1",
    "-2+3",
    "@SUM(1)",
    '<script>alert("x")</script>',
    "a & b < c > d \" e ' f",
    "tab\there",
    "line1\nline2",
    "JN.1* (excludes KP.2*, KP.3*)",
    "&amp; &#x41; ]]>",
    "é中\U0001f9a0",
]


def build(tmp_path, sheets, name="w.xlsx"):
    path = tmp_path / name
    xw.write_workbook(path, sheets)
    return path


def parts(path):
    with zipfile.ZipFile(path) as archive:
        assert archive.testzip() is None
        return {n: archive.read(n) for n in archive.namelist()}


def load(path):
    """Parse every part (fails on any malformed XML) -> {part name: root element}."""
    return {n: ET.fromstring(b) for n, b in parts(path).items()}  # ParseError if not well formed


def cells(sheet_tree):
    """{ref: (type, style, value)} with inline strings joined."""
    found = {}
    for c in sheet_tree.iter(M + "c"):
        text = "".join(t.text or "" for t in c.iter(M + "t"))
        v = c.find(M + "v")
        found[c.get("r")] = (c.get("t"), int(c.get("s", 0)), v.text if v is not None else text)
    return found


def simple(**extra):
    sheet = {
        "name": "Data",
        "header": ["sample", "x", "note"],
        "rows": [["S1", "0.5", "a"], ["S2", "", ""], ["S3", "0.02", "b"]],
        "numeric_columns": [1],
    }
    sheet.update(extra)
    return sheet


def test_package_structure_is_complete_and_consistent(tmp_path):
    path = build(tmp_path, [simple(), simple(name="Second")])
    raw = parts(path)
    trees = load(path)
    assert next(iter(raw)) == "[Content_Types].xml"
    assert set(raw) == {
        "[Content_Types].xml",
        "_rels/.rels",
        "xl/workbook.xml",
        "xl/_rels/workbook.xml.rels",
        "xl/styles.xml",
        "xl/worksheets/sheet1.xml",
        "xl/worksheets/sheet2.xml",
    }
    ct = trees["[Content_Types].xml"]
    overrides = {o.get("PartName") for o in ct}
    assert {"/xl/workbook.xml", "/xl/styles.xml", "/xl/worksheets/sheet2.xml"} <= overrides
    root_rel = trees["_rels/.rels"][0]
    assert root_rel.get("Target") == "xl/workbook.xml"
    wb = trees["xl/workbook.xml"]
    sheets = wb.find(M + "sheets")
    assert [s.get("name") for s in sheets] == ["Data", "Second"]
    rels = {r.get("Id"): r.get("Target") for r in trees["xl/_rels/workbook.xml.rels"]}
    rid = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
    for sheet in sheets:
        assert f"xl/{rels[sheet.get(rid)]}" in raw
    assert any(t == "styles.xml" for t in rels.values())


def test_styles_table_covers_every_style_id_used(tmp_path):
    heat = {"name": "H", "header": ["id", "v"], "numeric_columns": [1], "heat_columns": [1]}
    heat["rows"] = [
        [f"r{i}", v] for i, v in enumerate(["0.01", "0.07", "0.15", "0.3", "0.4", "0.9", "0", ""])
    ]
    trees = load(build(tmp_path, [simple(), heat]))
    styles = trees["xl/styles.xml"]
    xfs = styles.find(M + "cellXfs")
    assert int(xfs.get("count")) == len(xfs)
    fills = styles.find(M + "fills")
    assert int(fills.get("count")) == len(fills)
    fonts = styles.find(M + "fonts")
    borders = styles.find(M + "borders")
    custom = {int(n.get("numFmtId")) for n in styles.find(M + "numFmts")}
    for xf in xfs:
        assert int(xf.get("fillId")) < len(fills)
        assert int(xf.get("fontId")) < len(fonts)
        assert int(xf.get("borderId")) < len(borders)
        fmt = int(xf.get("numFmtId"))
        assert fmt == 0 or fmt in custom
    for name in ("xl/worksheets/sheet1.xml", "xl/worksheets/sheet2.xml"):
        for _, style, _ in cells(trees[name]).values():
            assert 0 <= style < len(xfs)
    assert len(xfs) < 20  # the style table stays small


def test_heat_buckets_use_six_distinct_fills(tmp_path):
    values = ["0.01", "0.05", "0.1", "0.2", "0.35", "0.5", "0.99", "0.049999"]
    sheet = {"name": "H", "header": ["id", "v"], "numeric_columns": [1], "heat_columns": [1]}
    sheet["rows"] = [[f"r{i}", v] for i, v in enumerate(values)] + [["z", "0"], ["b", ""]]
    trees = load(build(tmp_path, [sheet]))
    found = cells(trees["xl/worksheets/sheet1.xml"])
    xfs = list(trees["xl/styles.xml"].find(M + "cellXfs"))
    fill = lambda ref: int(xfs[found[ref][1]].get("fillId"))
    # 0.01, 0.049999 -> bucket 0; 0.05 -> 1; 0.1 -> 2; 0.2 -> 3; 0.35 -> 4; 0.5, 0.99 -> 5
    assert fill("B2") == fill("B9") == 3
    assert [fill(f"B{r}") for r in (3, 4, 5, 6, 7)] == [4, 5, 6, 7, 8]
    assert fill("B8") == 8
    assert fill("B10") == 0  # zero is not shaded
    assert "B11" not in found  # blank stays a missing cell
    assert len({fill(f"B{r}") for r in range(2, 9)}) == 6


def test_numbers_strings_and_blanks(tmp_path):
    sheet = {
        "name": "N",
        "header": ["id", "num"],
        "numeric_columns": ["num"],
        "rows": [
            ["007", "0.123456"],
            ["1274", "1e-5"],
            ["x", "NA"],
            ["y", "nan"],
            ["z", "1,5"],
            ["w", "  "],
            ["v", None],
            ["u", " 12 "],
            ["t", 7],
            ["s", -3.5],
            ["r", "inf"],
        ],
    }
    trees = load(build(tmp_path, [sheet]))
    found = cells(trees["xl/worksheets/sheet1.xml"])
    assert found["A2"][0] == "inlineStr" and found["A2"][2] == "007"  # IDs stay text
    assert found["A3"][0] == "inlineStr" and found["A3"][2] == "1274"
    assert found["B2"] == (None, found["B2"][1], "0.123456")
    assert found["B3"][2] == "1e-05"
    assert float(found["B3"][2]) == 1e-5
    for ref in ("B4", "B5", "B6", "B12"):  # unparsable stays text
        assert found[ref][0] == "inlineStr"
    assert "B7" not in found and "B8" not in found  # blanks: no cell
    assert found["B9"][2] == "12" and found["B9"][0] is None
    assert found["B10"][2] == "7" and found["B11"][2] == "-3.5"


def test_hostile_strings_are_escaped_and_never_formulas(tmp_path):
    sheet = {
        "name": "H",
        "header": ["h<&>\"'"] + ["c"],
        "rows": [[s, s] for s in HOSTILE],
    }
    path = build(tmp_path, [sheet])
    raw = parts(path)
    trees = load(path)
    xml = raw["xl/worksheets/sheet1.xml"].decode()
    assert "<f>" not in xml and "<f " not in xml
    found = cells(trees["xl/worksheets/sheet1.xml"])
    for i, text in enumerate(HOSTILE, 2):
        for col in "AB":
            kind, _, value = found[f"{col}{i}"]
            assert kind == "inlineStr"
            assert value == text.replace("\r", "")
    assert "<script>" not in xml  # escaped, not raw
    assert found["A1"][2] == "h<&>\"'"


def test_control_characters_are_stripped_and_line_endings_normalised(tmp_path):
    sheet = {
        "name": "C",
        "header": ["a"],
        "rows": [["x\x00y\x01\x08z\x0b\x0c\x1f￾!"], ["a\r\nb\rc"]],
    }
    trees = load(build(tmp_path, [sheet]))
    found = cells(trees["xl/worksheets/sheet1.xml"])
    assert found["A2"][2] == "xyz!"
    assert found["A3"][2] == "a\nb\nc"


def test_overlong_cell_is_truncated_to_excels_limit(tmp_path):
    trees = load(build(tmp_path, [{"name": "L", "header": ["a"], "rows": [["x" * 40000]]}]))
    assert len(cells(trees["xl/worksheets/sheet1.xml"])["A2"][2]) == 32767


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("Lineage summary", "Lineage summary"),
        ("a[b]c:d*e?f/g\\h", "a b c d e f g h"),
        ("'quoted'", "quoted"),
        ("x" * 50, "x" * 31),
        ("", "Sheet"),
        ("History", "History_"),
    ],
)
def test_sheet_names_follow_excels_rules(name, expected):
    assert xw.sheet_title(name, set()) == expected


def test_duplicate_sheet_names_are_made_unique_case_insensitively(tmp_path):
    long = "y" * 31
    sheets = [simple(name="Data"), simple(name="DATA"), simple(name=long), simple(name=long)]
    trees = load(build(tmp_path, sheets))
    names = [s.get("name") for s in trees["xl/workbook.xml"].find(M + "sheets")]
    assert names == ["Data", "DATA_2", long, long[:29] + "_2"]
    assert len({n.casefold() for n in names}) == 4
    assert all(len(n) <= 31 for n in names)


def test_dimension_cols_freeze_and_filter(tmp_path):
    trees = load(build(tmp_path, [simple(freeze=(1, 2))]))
    sheet = trees["xl/worksheets/sheet1.xml"]
    assert sheet.find(M + "dimension").get("ref") == "A1:C4"
    pane = sheet.find(f"{M}sheetViews/{M}sheetView/{M}pane")
    assert (pane.get("xSplit"), pane.get("ySplit"), pane.get("topLeftCell")) == ("2", "1", "C2")
    assert pane.get("state") == "frozen" and pane.get("activePane") == "bottomRight"
    assert sheet.find(M + "autoFilter").get("ref") == "A1:C4"
    names = trees["xl/workbook.xml"].find(M + "definedNames")
    assert names[0].get("name") == "_xlnm._FilterDatabase"
    assert names[0].text == "'Data'!$A$1:$C$4"
    cols = sheet.find(M + "cols")
    assert [int(c.get("min")) for c in cols] == [1, 2, 3]
    assert all(float(c.get("width")) > 0 for c in cols)
    # row-only / column-only / no freeze, no autofilter
    for freeze, pane_name in (((1, 0), "bottomLeft"), ((0, 1), "topRight"), ((0, 0), None)):
        t = load(build(tmp_path, [simple(freeze=freeze, autofilter=False)], "f.xlsx"))
        s = t["xl/worksheets/sheet1.xml"]
        p = s.find(f"{M}sheetViews/{M}sheetView/{M}pane")
        assert (p.get("activePane") if p is not None else None) == pane_name
        assert s.find(M + "autoFilter") is None
        assert t["xl/workbook.xml"].find(M + "definedNames") is None


def test_element_order_follows_the_schema(tmp_path):
    trees = load(build(tmp_path, [simple()]))
    order = [el.tag.replace(M, "") for el in trees["xl/worksheets/sheet1.xml"]]
    assert order == [
        "dimension",
        "sheetViews",
        "sheetFormatPr",
        "cols",
        "sheetData",
        "autoFilter",
        "pageMargins",
        "ignoredErrors",
    ]


def test_rows_are_ascending_and_cells_in_column_order(tmp_path):
    trees = load(build(tmp_path, [simple()]))
    rows = trees["xl/worksheets/sheet1.xml"].find(M + "sheetData")
    numbers = [int(r.get("r")) for r in rows]
    assert numbers == sorted(numbers) == [1, 2, 3, 4]
    for row in rows:
        refs = [c.get("r") for c in row]
        letters = [re.match("[A-Z]+", r).group() for r in refs]
        assert letters == sorted(letters)
        assert all(r.endswith(str(row.get("r"))) for r in refs)


def test_wrap_and_widths_from_content_with_a_cap(tmp_path):
    long = " ".join(f"XFG.{i}" for i in range(60))
    sheet = {
        "name": "W",
        "header": ["id", "list", "short"],
        "rows": [["a", long, "x"], ["b", "", "y"]],
    }
    trees = load(build(tmp_path, [sheet]))
    sheet_tree = trees["xl/worksheets/sheet1.xml"]
    widths = [float(c.get("width")) for c in sheet_tree.find(M + "cols")]
    assert widths[1] == xw.DEFAULT_MAX_WIDTH and widths[2] < 12
    xfs = list(trees["xl/styles.xml"].find(M + "cellXfs"))
    found = cells(sheet_tree)
    wrap = lambda ref: xfs[found[ref][1]].find(M + "alignment").get("wrapText")
    assert wrap("B2") == "1" and wrap("C2") is None
    assert all(r.get("ht") is None for r in sheet_tree.iter(M + "row"))  # Excel autosizes heights


def test_alignment_numbers_right_text_left(tmp_path):
    trees = load(build(tmp_path, [simple()]))
    xfs = list(trees["xl/styles.xml"].find(M + "cellXfs"))
    found = cells(trees["xl/worksheets/sheet1.xml"])
    align = lambda ref: xfs[found[ref][1]].find(M + "alignment").get("horizontal")
    assert align("A2") == "left" and align("B2") == "right" and align("C2") == "left"
    assert align("B1") == "right" and align("A1") == "left"


def test_number_format_default_and_override(tmp_path):
    sheet = simple(number_formats={"x": "0.0000"})
    trees = load(build(tmp_path, [sheet, simple(decimals=2), simple(decimals=0)]))
    codes = [n.get("formatCode") for n in trees["xl/styles.xml"].find(M + "numFmts")]
    assert codes == ["0.0000", "0.00", "0"]
    trees = load(build(tmp_path, [simple()], "d.xlsx"))
    assert [n.get("formatCode") for n in trees["xl/styles.xml"].find(M + "numFmts")] == ["0.000"]


def test_unknown_column_name_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="no column named"):
        xw.write_workbook(tmp_path / "e.xlsx", [simple(numeric_columns=["nope"])])
    with pytest.raises(ValueError, match="out of range"):
        xw.write_workbook(tmp_path / "e.xlsx", [simple(heat_columns=[9])])
    with pytest.raises(ValueError, match="at least one sheet"):
        xw.write_workbook(tmp_path / "e.xlsx", [])


def test_empty_and_ragged_sheets_are_valid(tmp_path):
    sheets = [
        {"name": "Empty", "header": [], "rows": []},
        {"name": "HeaderOnly", "header": ["a", "b"], "rows": []},
        {"name": "Ragged", "header": ["a", "b", "c"], "rows": [["1"], ["1", "2", "3", "4"], []]},
    ]
    trees = load(build(tmp_path, sheets))
    assert trees["xl/worksheets/sheet1.xml"].find(M + "dimension").get("ref") == "A1"
    assert trees["xl/worksheets/sheet2.xml"].find(M + "dimension").get("ref") == "A1:B1"
    ragged = cells(trees["xl/worksheets/sheet3.xml"])
    assert "D2" not in ragged and "D3" not in ragged and ragged["B3"][2] == "2"


def test_output_is_deterministic_and_zip_metadata_is_fixed(tmp_path):
    sheets = [simple(), simple(name="B", heat_columns=["x"])]
    a = build(tmp_path, sheets, "a.xlsx")
    time.sleep(1.1)
    b = build(tmp_path, sheets, "b.xlsx")
    assert a.read_bytes() == b.read_bytes()
    with zipfile.ZipFile(a) as archive:
        for info in archive.infolist():
            assert info.date_time == (1980, 1, 1, 0, 0, 0)
            assert info.compress_type == zipfile.ZIP_DEFLATED


def test_large_workbook_is_fast(tmp_path):
    header = [f"col {i}" for i in range(40)]
    rows = [
        [f"S{r}"] + [f"{(r * c % 97) / 100:.3f}" for c in range(1, 39)] + ["w " * 20]
        for r in range(1000)
    ]
    sheet = {"name": "Big", "header": header, "rows": rows, "numeric_columns": list(range(1, 39)),
             "heat_columns": list(range(1, 39))}  # fmt: skip
    start = time.perf_counter()
    path = build(tmp_path, [sheet])
    elapsed = time.perf_counter() - start
    assert elapsed < 2
    trees = load(path)
    assert trees["xl/worksheets/sheet1.xml"].find(M + "dimension").get("ref") == "A1:AN1001"


def test_column_letters():
    assert [xw.column_letter(i) for i in (0, 25, 26, 51, 52, 701, 702, 16383)] == [
        "A", "Z", "AA", "AZ", "BA", "ZZ", "AAA", "XFD"
    ]  # fmt: skip


# ---- CLI --------------------------------------------------------------------------------------
def run_cli(*args):
    return subprocess.run(
        [sys.executable, str(REPO / "bin" / "xlsx_writer.py"), *map(str, args)],
        capture_output=True,
        text=True,
        check=False,
    )


def test_cli_converts_csv_and_tsv_files(tmp_path):
    (tmp_path / "a.csv").write_text(
        "sample,depth,frac,code\nS1,100,0.25,007\nS2,200,0.5,010\nS3,,,x\n"
    )
    (tmp_path / "qc.tsv").write_text("﻿sample\tn\tnote\nS1\t3\t=1+1\n")
    (tmp_path / "semi.txt").write_text("a;b\n1;2,5\n")
    out = tmp_path / "o.xlsx"
    proc = run_cli("--from-csv", tmp_path / "a.csv", "--sheet-name", "First",
                   "--from-csv", tmp_path / "qc.tsv", "--from-csv", tmp_path / "semi.txt",
                   "--heat-column", "frac", "-o", out)  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    trees = load(out)
    assert [s.get("name") for s in trees["xl/workbook.xml"].find(M + "sheets")] == [
        "First", "qc", "semi"
    ]  # fmt: skip
    first = cells(trees["xl/worksheets/sheet1.xml"])
    assert first["B2"][0] is None and first["B2"][2] == "100"  # numeric column
    assert first["C2"][2] == "0.25"
    assert first["D2"][0] == "inlineStr" and first["D2"][2] == "007"  # leading zero: text
    qc = cells(trees["xl/worksheets/sheet2.xml"])
    assert qc["A1"][2] == "sample"  # BOM removed
    assert qc["C2"][0] == "inlineStr" and qc["C2"][2] == "=1+1"
    assert cells(trees["xl/worksheets/sheet3.xml"])["B2"][2] == "2,5"  # delimiter sniffed as ;
    codes = [n.get("formatCode") for n in trees["xl/styles.xml"].find(M + "numFmts")]
    assert "0.00" in codes  # frac has two decimals: displayed with two


def test_cli_errors_are_friendly(tmp_path):
    (tmp_path / "empty.csv").write_text("")
    proc = run_cli("--from-csv", tmp_path / "empty.csv", "-o", tmp_path / "o.xlsx")
    assert proc.returncode == 2 and "is empty" in proc.stderr and "Traceback" not in proc.stderr
    proc = run_cli("--from-csv", tmp_path / "missing.csv", "-o", tmp_path / "o.xlsx")
    assert proc.returncode == 2 and "Traceback" not in proc.stderr
    assert (
        run_cli(
            "--from-csv", "a.csv", "--sheet-name", "a", "--sheet-name", "b", "-o", "o.xlsx"
        ).returncode
        == 2
    )


def test_script_is_host_python_compatible():
    """The converter may run on a user's older host Python: parse as 3.8 and avoid newer calls."""
    import ast

    source = (REPO / "bin" / "xlsx_writer.py").read_text()
    ast.parse(source, feature_version=(3, 8))
    assert "strict=" not in source and "removeprefix" not in source


def test_cli_on_the_hostile_and_qc_like_fixtures(tmp_path):
    fixtures = REPO / "tests" / "fixtures" / "xlsx"
    out = tmp_path / "o.xlsx"
    proc = run_cli(
        "--from-csv", fixtures / "hostile.csv", "--sheet-name", "Hostile: [x]/y",
        "--from-csv", fixtures / "qc_summary_like.tsv", "--heat-column", "pct_covered", "-o", out,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    trees = load(out)
    names = [s.get("name") for s in trees["xl/workbook.xml"].find(M + "sheets")]
    assert names == ["Hostile   x  y", "qc_summary_like"]
    hostile = cells(trees["xl/worksheets/sheet1.xml"])
    assert hostile["A2"] == ("inlineStr", hostile["A2"][1], "S<1>")
    assert hostile["A1"][2] == "sample" and hostile["B1"][2] == "group a & b"
    assert hostile["D2"][0] == "inlineStr" and hostile["D2"][2] == "=cmd|' /C calc'!A0"
    assert hostile["A4"][2] == "007" and hostile["D5"][2] == "@SUM(1)"
    assert hostile["D6"][2] == "line1\nline2" and hostile["A6"][2] == "café 中"
    assert hostile["B4"][0] is None and hostile["B4"][2] == "0.9"  # all-number column
    assert hostile["C4"][0] == "inlineStr" and hostile["C5"][2] == "nan"  # 'nan' keeps it text
    qc = cells(trees["xl/worksheets/sheet2.xml"])
    assert qc["C2"][2] == "812.4" and qc["C2"][0] is None
    assert qc["E3"][0] == "inlineStr" and qc["E4"][2] == "ivar variants: empty & aborted"
    assert "C4" not in qc  # blank depth: no cell
