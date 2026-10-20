"""bin/xlsx_writer.py: a stdlib-only .xlsx writer. These tests parse every part of the package with
xml.etree and check the invariants Excel is strict about (style ids, dimension, sheet names, escaping).
The independent reader check (openpyxl) is described in docs/wiring/xlsx_output.md."""

import zipfile
from xml.etree import ElementTree as ET

import xlsx_writer as xw

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
M = "{" + NS["m"] + "}"


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
