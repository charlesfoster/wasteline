#!/usr/bin/env python3
"""Minimal, deterministic .xlsx writer (Office Open XML) using only the standard library.

Library use
    from xlsx_writer import write_workbook
    write_workbook("out.xlsx", [
        {"name": "Lineage summary",             # sheet name (sanitised to Excel's rules)
         "header": ["sample", "JN.1*", "notes"],
         "rows": [["S1", "0.5", "text"], ["S2", "", ""]],
         "numeric_columns": [1],                # indices or header names; cells become numbers
         "number_formats": {"JN.1*": "0.000"},  # optional; default "0.000" (see "decimals")
         "heat_columns": [1],                   # optional; static 6-bucket fill on values > 0
         "freeze": (1, 1),                      # (rows, columns) frozen; default (1, 0)
         "widths": None,                        # optional {column: width}; default from content
         "autofilter": True}                    # default True
    ])

Guarantees
  * Every string is XML-escaped; characters XML 1.0 forbids are dropped. Text is always written as an
    inline string (never a formula), so `=cmd|...`, `+1`, `-1`, `@x` stay text. No formulas, macros,
    merged cells, charts or conditional formatting are ever written.
  * A value is a number cell only in a numeric column (or when it is a Python int/float) and only if it
    is a plain decimal number; blanks are empty cells (nothing is written for them).
  * Identical input gives byte-identical output (fixed zip timestamps, order, compression).

CLI (convert CSV/TSV files, e.g. qc_summary.tsv):
    xlsx_writer.py --from-csv a.csv [--sheet-name A] --from-csv b.tsv [--sheet-name B] -o out.xlsx
Columns whose non-empty cells are all plain numbers (no leading zeros) become numeric, except the
first column (an identifier). The first row is frozen, with the first column.

Standard library only; written to run on Python 3.8+.
"""

import argparse
import csv
import io
import math
import re
import sys
import zipfile
from pathlib import Path

NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"
NS_CT = "http://schemas.openxmlformats.org/package/2006/content-types"
REL_OFFICE_DOC = NS_REL + "/officeDocument"
REL_WORKSHEET = NS_REL + "/worksheet"
REL_STYLES = NS_REL + "/styles"
XML_DECL = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
CT_PREFIX = "application/vnd.openxmlformats-officedocument.spreadsheetml."

MAX_CELL_CHARS = 32767  # Excel's limit per cell
MAX_ROWS = 1048576
MAX_COLS = 16384
DEFAULT_MAX_WIDTH = 50  # cap for content-derived column widths (characters)
DEFAULT_DECIMALS = 3
# light sequential blues for heat shading: 0<x<0.05, <0.10, <0.20, <0.35, <0.50, >=0.50
HEAT_EDGES = (0.05, 0.10, 0.20, 0.35, 0.50)
HEAT_COLOURS = ("E3EEF9", "CFE2F3", "B4D3EC", "93BFE3", "74AADA", "5B96CE")
HEADER_FILL = "E7E6E6"
_FORBIDDEN_SHEET_CHARS = re.compile(r"[\[\]:*?/\\]")
# XML 1.0 forbids C0 controls except tab/newline/CR, and U+FFFE/U+FFFF; lone surrogates cannot be encoded.
_BAD_XML = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f￾￿\ud800-\udfff]")
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")
_LEADING_ZERO = re.compile(r"[-+]?0\d")
_PLAIN_DECIMAL = re.compile(r"-?(?:0|[1-9]\d*)(?:\.(\d+))?")


# --------------------------------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------------------------------
def column_letter(index):
    """0 -> A, 25 -> Z, 26 -> AA."""
    letters = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        letters = chr(65 + rem) + letters
    return letters


def clean_text(value):
    """Text safe for XML 1.0: bad characters dropped, line endings normalised, length capped."""
    text = _BAD_XML.sub("", str(value))
    if "\r" in text:
        text = text.replace("\r\n", "\n").replace("\r", "\n")
    return text[:MAX_CELL_CHARS]


def escape(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def escape_attr(text):
    return escape(text).replace('"', "&quot;").replace("\n", "&#10;").replace("\t", "&#9;")


def sheet_title(name, taken):
    """Valid, unique (case-insensitive) sheet name: <= 31 chars, none of []:*?/\\, no edge apostrophes."""
    title = _FORBIDDEN_SHEET_CHARS.sub(" ", clean_text(name)).replace("\n", " ").strip()
    title = title.strip("'").strip()[:31].strip("'").strip() or "Sheet"
    if title.casefold() == "history":  # reserved by Excel
        title = "History_"
    base, n = title, 1
    while title.casefold() in taken:
        n += 1
        tail = f"_{n}"
        title = base[: 31 - len(tail)] + tail
    taken.add(title.casefold())
    return title


def parse_number(text):
    """float for a plain decimal number string, else None (nan/inf and mixed text stay text)."""
    text = text.strip()
    if not _NUMBER.fullmatch(text):
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def number_text(number):
    """Shortest exact text of a number for a <v> element."""
    if isinstance(number, int):
        return str(number)
    if number == int(number) and abs(number) < 1e15:
        return str(int(number))
    return repr(number)


def heat_bucket(number):
    """0..5 for number > 0, else None."""
    if number <= 0:
        return None
    for i, edge in enumerate(HEAT_EDGES):
        if number < edge:
            return i
    return len(HEAT_EDGES)


def number_format_for(decimals):
    return "0" if decimals <= 0 else "0." + "0" * decimals


# --------------------------------------------------------------------------------------------------
# styles
# --------------------------------------------------------------------------------------------------
class Styles:
    """Builds cellXfs on demand (in order of first use) so the style table stays tiny and valid."""

    def __init__(self):
        self.formats = {}  # format code -> numFmtId (custom, from 164)
        self.xfs = []  # (numFmtId, fontId, fillId, borderId, horizontal, vertical, wrap)
        self.index = {}
        self.get(0, 0, 0, 0, None, None, False)  # xf 0: the default for untouched cells
        self.fills = [HEADER_FILL, *HEAT_COLOURS]  # fill ids 2.. (0 none, 1 gray125 are fixed)

    def format_id(self, code):
        if code in (None, "General"):
            return 0
        if code not in self.formats:
            self.formats[code] = 164 + len(self.formats)
        return self.formats[code]

    def get(self, num_fmt, font, fill, border, horizontal, vertical, wrap):
        key = (num_fmt, font, fill, border, horizontal, vertical, wrap)
        if key not in self.index:
            self.index[key] = len(self.xfs)
            self.xfs.append(key)
        return self.index[key]

    def header(self, numeric):
        return self.get(0, 1, 2, 1, "right" if numeric else "left", "center", True)

    def text(self, wrap):
        return self.get(0, 0, 0, 0, "left", "top", wrap)

    def number(self, code, bucket):
        fill = 0 if bucket is None else 3 + bucket
        return self.get(self.format_id(code), 0, fill, 0, "right", "top", False)

    def xml(self):
        out = [XML_DECL, f'<styleSheet xmlns="{NS_MAIN}">']
        if self.formats:
            out.append(f'<numFmts count="{len(self.formats)}">')
            for code, ident in self.formats.items():
                out.append(f'<numFmt numFmtId="{ident}" formatCode="{escape_attr(code)}"/>')
            out.append("</numFmts>")
        base = '<sz val="11"/><name val="Calibri"/><family val="2"/>'
        out.append(f'<fonts count="2"><font>{base}</font><font><b/>{base}</font></fonts>')
        out.append(f'<fills count="{len(self.fills) + 2}">')
        out.append('<fill><patternFill patternType="none"/></fill>')
        out.append('<fill><patternFill patternType="gray125"/></fill>')
        for colour in self.fills:
            out.append(
                f'<fill><patternFill patternType="solid"><fgColor rgb="FF{colour}"/>'
                '<bgColor indexed="64"/></patternFill></fill>'
            )
        out.append("</fills>")
        out.append(
            '<borders count="2"><border><left/><right/><top/><bottom/><diagonal/></border>'
            '<border><left/><right/><top/><bottom style="thin"><color rgb="FF808080"/></bottom>'
            "<diagonal/></border></borders>"
        )
        out.append(
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>'
            "</cellStyleXfs>"
        )
        out.append(f'<cellXfs count="{len(self.xfs)}">')
        for num_fmt, font_id, fill, border, horizontal, vertical, wrap in self.xfs:
            attrs = f'numFmtId="{num_fmt}" fontId="{font_id}" fillId="{fill}" borderId="{border}" xfId="0"'
            if num_fmt:
                attrs += ' applyNumberFormat="1"'
            if font_id:
                attrs += ' applyFont="1"'
            if fill:
                attrs += ' applyFill="1"'
            if border:
                attrs += ' applyBorder="1"'
            if horizontal or vertical or wrap:
                align = ""
                if horizontal:
                    align += f' horizontal="{horizontal}"'
                if vertical:
                    align += f' vertical="{vertical}"'
                if wrap:
                    align += ' wrapText="1"'
                out.append(f'<xf {attrs} applyAlignment="1"><alignment{align}/></xf>')
            else:
                out.append(f"<xf {attrs}/>")
        out.append("</cellXfs>")
        out.append(
            '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        )
        out.append(
            '<dxfs count="0"/><tableStyles count="0" defaultTableStyle="TableStyleMedium9"/>'
        )
        out.append("</styleSheet>")
        return "".join(out)


# --------------------------------------------------------------------------------------------------
# sheets
# --------------------------------------------------------------------------------------------------
def _column_index(key, header, what):
    """Column index from an index or a header name."""
    if isinstance(key, int) and not isinstance(key, bool):
        if not 0 <= key < len(header):
            raise ValueError(f"{what}: column index {key} is out of range")
        return key
    if key not in header:
        raise ValueError(f"{what}: no column named {key!r}")
    return list(header).index(key)


def _display_length(text):
    return max((len(line) for line in text.split("\n")), default=0)


def _cell_value(value, is_numeric):
    """-> ('n', number) | ('s', text) | None (blank)."""
    if value is None:
        return None
    if isinstance(value, bool):
        return ("s", "TRUE" if value else "FALSE")
    if isinstance(value, (int, float)):
        if isinstance(value, float) and not math.isfinite(value):
            return None
        return ("n", value)
    text = clean_text(value)
    if text == "" or (is_numeric and text.strip() == ""):
        return None
    if is_numeric:
        number = parse_number(text)
        if number is not None:
            return ("n", number)
    return ("s", text)


def _column_widths(sheet, header, longest, numeric, has_newline, max_width, what):
    given = {_column_index(k, header, what): w for k, w in (sheet.get("widths") or {}).items()}
    widths, wrap = [], []
    for c, title in enumerate(header):
        if c in given:
            widths.append(float(given[c]))
            wrap.append(longest[c] > given[c] or has_newline[c])
            continue
        if c in numeric:
            width = max(longest[c] + 2, min(len(title), 12) + 2, 8)
        else:
            width = max(longest[c] + 2, min(len(title), 20) + 2, 6)
        wrap.append(longest[c] + 2 > max_width or has_newline[c])
        widths.append(float(min(width, max_width)))
    return widths, wrap


def _sheet_view(selected, freeze_rows, freeze_cols):
    selected_attr = ' tabSelected="1"' if selected else ""
    view = f'<sheetViews><sheetView{selected_attr} workbookViewId="0">'
    if freeze_rows > 0 or freeze_cols > 0:
        top_left = f"{column_letter(freeze_cols)}{freeze_rows + 1}"
        if freeze_rows and freeze_cols:
            pane = "bottomRight"
        else:
            pane = "bottomLeft" if freeze_rows else "topRight"
        attrs = (f' xSplit="{freeze_cols}"' if freeze_cols else "") + (
            f' ySplit="{freeze_rows}"' if freeze_rows else ""
        )
        view += f'<pane{attrs} topLeftCell="{top_left}" activePane="{pane}" state="frozen"/>'
        if pane == "bottomRight":
            view += '<selection pane="topRight"/><selection pane="bottomLeft"/>'
        view += f'<selection pane="{pane}" activeCell="{top_left}" sqref="{top_left}"/>'
    return view + "</sheetView></sheetViews>"


def build_sheet(sheet, name, styles, selected=False, max_width=DEFAULT_MAX_WIDTH):
    """-> (worksheet xml, auto-filter range or None)."""
    header = [clean_text(h) for h in sheet.get("header") or []]
    ncols = len(header)
    rows = sheet.get("rows") or []
    if ncols > MAX_COLS:
        raise ValueError(f"{name}: more than {MAX_COLS} columns")
    if len(rows) + 1 > MAX_ROWS:
        raise ValueError(f"{name}: more than {MAX_ROWS - 1} rows")
    numeric = {_column_index(k, header, name) for k in sheet.get("numeric_columns") or ()}
    heat = {_column_index(k, header, name) for k in sheet.get("heat_columns") or ()}
    formats = {
        _column_index(k, header, name): v for k, v in (sheet.get("number_formats") or {}).items()
    }
    default_format = number_format_for(sheet.get("decimals", DEFAULT_DECIMALS))
    freeze_rows, freeze_cols = sheet.get("freeze", (1, 0))
    freeze_rows, freeze_cols = int(freeze_rows), min(int(freeze_cols), ncols)

    grid, longest, has_newline = [], [0] * ncols, [False] * ncols
    for row in rows:
        cells = []
        for c in range(ncols):
            cell = _cell_value(row[c] if c < len(row) else None, c in numeric)
            cells.append(cell)
            if cell is None:
                continue
            if cell[0] == "s":
                longest[c] = max(longest[c], _display_length(cell[1]))
                has_newline[c] = has_newline[c] or "\n" in cell[1]
            else:
                shown = len(formats.get(c) or number_text(cell[1])) + 1
                longest[c] = max(longest[c], min(shown, 12))
        grid.append(cells)
    widths, wrap = _column_widths(sheet, header, longest, numeric, has_newline, max_width, name)

    last_row = len(rows) + 1
    letters = [column_letter(c) for c in range(ncols)]
    dimension = f"A1:{letters[-1]}{last_row}" if ncols else "A1"
    out = [XML_DECL, f'<worksheet xmlns="{NS_MAIN}" xmlns:r="{NS_REL}">']
    out.append(f'<dimension ref="{dimension}"/>')
    out.append(_sheet_view(selected, freeze_rows, freeze_cols))
    out.append('<sheetFormatPr defaultRowHeight="15"/>')
    if ncols:
        out.append("<cols>")
        for c in range(ncols):
            out.append(f'<col min="{c + 1}" max="{c + 1}" width="{widths[c]}" customWidth="1"/>')
        out.append("</cols>")
    out.append("<sheetData>")
    if ncols:
        out.append('<row r="1">')
        for c in range(ncols):
            out.append(
                f'<c r="{letters[c]}1" s="{styles.header(c in numeric)}" t="inlineStr">'
                f'<is><t xml:space="preserve">{escape(header[c])}</t></is></c>'
            )
        out.append("</row>")
    text_style = [styles.text(w) for w in wrap]
    for r, cells in enumerate(grid):
        parts = []
        for c, cell in enumerate(cells):
            if cell is None:
                continue
            ref = f"{letters[c]}{r + 2}"
            if cell[0] == "s":
                parts.append(
                    f'<c r="{ref}" s="{text_style[c]}" t="inlineStr">'
                    f'<is><t xml:space="preserve">{escape(cell[1])}</t></is></c>'
                )
            else:
                code = formats.get(c) or (default_format if c in numeric else None)
                bucket = heat_bucket(cell[1]) if c in heat else None
                parts.append(
                    f'<c r="{ref}" s="{styles.number(code, bucket)}"><v>{number_text(cell[1])}</v></c>'
                )
        if parts:
            out.append(f'<row r="{r + 2}">{"".join(parts)}</row>')
    out.append("</sheetData>")
    use_filter = bool(sheet.get("autofilter", True)) and ncols > 0
    if use_filter:
        out.append(f'<autoFilter ref="A1:{letters[-1]}{last_row}"/>')
    out.append(
        '<pageMargins left="0.7" right="0.7" top="0.75" bottom="0.75" header="0.3" footer="0.3"/>'
    )
    if ncols:  # a sample ID such as 1274 stored as text is intentional
        out.append(
            f'<ignoredErrors><ignoredError sqref="A1:{letters[-1]}{last_row}" numberStoredAsText="1"/>'
            "</ignoredErrors>"
        )
    out.append("</worksheet>")
    return "".join(out), (f"$A$1:${letters[-1]}${last_row}" if use_filter else None)


# --------------------------------------------------------------------------------------------------
# workbook
# --------------------------------------------------------------------------------------------------
def build_parts(sheets, max_width=DEFAULT_MAX_WIDTH):
    """-> ordered list of (part name, xml text)."""
    if not sheets:
        raise ValueError("a workbook needs at least one sheet")
    styles = Styles()
    taken = set()
    names, sheet_xml, filters = [], [], []
    for i, sheet in enumerate(sheets):
        title = sheet_title(sheet.get("name", "Sheet"), taken)
        xml, filter_ref = build_sheet(sheet, title, styles, selected=(i == 0), max_width=max_width)
        names.append(title)
        sheet_xml.append(xml)
        filters.append(filter_ref)

    n = len(sheets)
    overrides = "".join(
        f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" ContentType="{CT_PREFIX}worksheet+xml"/>'
        for i in range(n)
    )
    content_types = (
        f'{XML_DECL}<Types xmlns="{NS_CT}">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        f'<Override PartName="/xl/workbook.xml" ContentType="{CT_PREFIX}sheet.main+xml"/>'
        f"{overrides}"
        f'<Override PartName="/xl/styles.xml" ContentType="{CT_PREFIX}styles+xml"/>'
        "</Types>"
    )
    root_rels = (
        f'{XML_DECL}<Relationships xmlns="{NS_PKG_REL}">'
        f'<Relationship Id="rId1" Type="{REL_OFFICE_DOC}" Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    defined = ""
    for i, ref in enumerate(filters):
        if ref:
            quoted = "'" + names[i].replace("'", "''") + "'"
            defined += (
                f'<definedName name="_xlnm._FilterDatabase" localSheetId="{i}" hidden="1">'
                f"{escape(quoted)}!{ref}</definedName>"
            )
    sheet_tags = "".join(
        f'<sheet name="{escape_attr(names[i])}" sheetId="{i + 1}" r:id="rId{i + 1}"/>'
        for i in range(n)
    )
    workbook = (
        f'{XML_DECL}<workbook xmlns="{NS_MAIN}" xmlns:r="{NS_REL}">'
        '<bookViews><workbookView xWindow="0" yWindow="0" windowWidth="28800" windowHeight="17000"/></bookViews>'
        f"<sheets>{sheet_tags}</sheets>"
        + (f"<definedNames>{defined}</definedNames>" if defined else "")
        + "</workbook>"
    )
    rels = "".join(
        f'<Relationship Id="rId{i + 1}" Type="{REL_WORKSHEET}" Target="worksheets/sheet{i + 1}.xml"/>'
        for i in range(n)
    )
    workbook_rels = (
        f'{XML_DECL}<Relationships xmlns="{NS_PKG_REL}">{rels}'
        f'<Relationship Id="rId{n + 1}" Type="{REL_STYLES}" Target="styles.xml"/>'
        "</Relationships>"
    )
    parts = [
        ("[Content_Types].xml", content_types),
        ("_rels/.rels", root_rels),
        ("xl/workbook.xml", workbook),
        ("xl/_rels/workbook.xml.rels", workbook_rels),
        ("xl/styles.xml", styles.xml()),
    ]
    parts += [(f"xl/worksheets/sheet{i + 1}.xml", sheet_xml[i]) for i in range(n)]
    return parts


def workbook_bytes(sheets, max_width=DEFAULT_MAX_WIDTH):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, xml in build_parts(sheets, max_width):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o644 << 16
            archive.writestr(info, xml.encode("utf-8"), zipfile.ZIP_DEFLATED, 6)
    return buffer.getvalue()


def write_workbook(path, sheets, max_width=DEFAULT_MAX_WIDTH):
    """Write `sheets` (list of dicts, see the module docstring) to `path` as an .xlsx file."""
    data = workbook_bytes(sheets, max_width)
    Path(path).write_bytes(data)


# --------------------------------------------------------------------------------------------------
# CSV/TSV conversion
# --------------------------------------------------------------------------------------------------
def sniff_delimiter(path, first_line):
    if Path(path).suffix.lower() in (".tsv", ".tab"):
        return "\t"
    best, best_n = ",", 0
    for delim in (",", "\t", ";"):
        n = len(next(csv.reader([first_line], delimiter=delim), []))
        if n > best_n:
            best, best_n = delim, n
    return best


def auto_numeric(rows, ncols):
    """Columns (except the first) whose non-empty cells are all plain numbers without leading zeros."""
    found = []
    for c in range(1, ncols):
        values = [row[c].strip() for row in rows if c < len(row) and row[c].strip()]
        if values and all(
            parse_number(v) is not None and not _LEADING_ZERO.match(v) for v in values
        ):
            found.append(c)
    return found


def decimals_of(rows, c):
    most = 0
    for row in rows:
        if c < len(row):
            match = _PLAIN_DECIMAL.fullmatch(row[c].strip())
            if match and match.group(1):
                most = max(most, len(match.group(1)))
    return min(most, 10)


def sheet_from_table(path, name=None, freeze=(1, 1), heat=()):
    """Sheet dict from a CSV/TSV file (delimiter sniffed; UTF-8, BOM tolerated)."""
    text = Path(path).read_bytes().decode("utf-8-sig", errors="replace")
    first = next((ln for ln in text.splitlines() if ln.strip()), "")
    reader = csv.reader(io.StringIO(text, newline=""), delimiter=sniff_delimiter(path, first))
    table = list(reader)
    table = [r for r in table if any(c.strip() for c in r)]
    if not table:
        raise ValueError(f"{path} is empty")
    header, rows = table[0], table[1:]
    numeric = auto_numeric(rows, len(header))
    formats = {}
    for c in numeric:
        places = decimals_of(rows, c)
        formats[c] = number_format_for(places) if places else "General"
    return {
        "name": name or Path(path).stem,
        "header": header,
        "rows": rows,
        "numeric_columns": numeric,
        "number_formats": formats,
        "heat_columns": [header.index(h) for h in heat if h in header],
        "freeze": freeze,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert CSV/TSV files into one formatted .xlsx workbook."
    )
    parser.add_argument(
        "--from-csv",
        action="append",
        required=True,
        metavar="FILE",
        help="CSV/TSV file (repeatable; one sheet each)",
    )
    parser.add_argument(
        "--sheet-name",
        action="append",
        default=[],
        metavar="NAME",
        help="sheet name for the matching --from-csv (default: the file name without extension)",
    )
    parser.add_argument(
        "--heat-column",
        action="append",
        default=[],
        metavar="HEADER",
        help="shade this numeric column by magnitude (repeatable; matched by header text)",
    )
    parser.add_argument("--freeze-rows", type=int, default=1)
    parser.add_argument("--freeze-cols", type=int, default=1)
    parser.add_argument("-o", "--output", required=True, help="output .xlsx")
    args = parser.parse_args(argv)
    if len(args.sheet_name) > len(args.from_csv):
        parser.error("more --sheet-name than --from-csv")
    names = args.sheet_name + [None] * (len(args.from_csv) - len(args.sheet_name))
    freeze = (args.freeze_rows, args.freeze_cols)
    try:
        sheets = [
            sheet_from_table(f, n, freeze, args.heat_column) for f, n in zip(args.from_csv, names)
        ]
        write_workbook(args.output, sheets)
    except (OSError, ValueError) as err:
        print(f"xlsx_writer.py: error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
