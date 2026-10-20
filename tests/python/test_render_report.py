import json
import re
from html.parser import HTMLParser
from pathlib import Path

import build_report_data as brd
import render_report as rr

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures" / "report"
TEMPLATE = ROOT / "assets" / "templates" / "report.html.j2"

VOID = {"meta", "link", "br", "hr", "img", "input", "col", "area", "base", "source", "wbr"}


class Balance(HTMLParser):
    """Tracks tag balance; self-closed SVG shapes (`<rect/>`) do not open anything."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack, self.errors, self.tags = [], [], []

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))
        if tag not in VOID:
            self.stack.append(tag)

    def handle_startendtag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if not self.stack or self.stack[-1] != tag:
            self.errors.append(f"unexpected </{tag}> with stack tail {self.stack[-3:]}")
            if tag in self.stack:
                while self.stack and self.stack.pop() != tag:
                    pass
        else:
            self.stack.pop()


def parse(html):
    parser = Balance()
    parser.feed(html)
    parser.close()
    return parser


def assert_valid(html):
    parser = parse(html)
    assert parser.errors == []
    assert parser.stack == []
    assert html.lstrip().lower().startswith("<!doctype html>")
    return parser


def build_fixture(name, tmp_path, depths=False):
    base = FIX / name
    mapping = {
        "--qc": "qc_summary.tsv",
        "--amplicons": "amplicon_coverage.tsv",
        "--aggregated": "freyja_aggregated.tsv",
        "--demix-status": "demix_status.tsv",
        "--moi": "mutations_of_interest.tsv",
        "--consensus": "consensus_stats.tsv",
        "--barcode-compare": "barcode_compare.tsv",
    }
    args = ["--run-info", str(base / "run_info.json")]
    for option, fname in mapping.items():
        if (base / fname).exists():
            args += [option, str(base / fname)]
    if depths:
        ddir = tmp_path / "depths"
        ddir.mkdir()
        for sample, scale in (
            ("WW_001", 800),
            ("WW_002", 640),
            ("WW_003.rep1", 400),
            ("NTC_01", 0),
        ):
            lines = [
                f"NC_045512.2\t{p}\tA\t{0 if (sample == 'WW_003.rep1' and 4800 < p < 5200) else scale}"
                for p in range(1, 29904)
            ]
            (ddir / f"{sample}.depths").write_text("\n".join(lines) + "\n")
        args += ["--depths", str(ddir)]
    out = tmp_path / "report_data.json"
    assert brd.main(args + ["--out", str(out)]) == 0
    return json.loads(out.read_text())


def render(data):
    return rr.render(data, TEMPLATE)


def test_full_run_content(tmp_path):
    html = render(build_fixture("a_full", tmp_path, depths=True))
    assert_valid(html)
    # barcode version prominent in the header, before any other section
    header = html[: html.index('id="sec-flags"')]
    assert 'id="barcode-version"' in header and "09_17_2026-00-46" in header
    assert "0.1.0dev" in header and "3f9c2ab71d0e" in header and "2.0.5" in header
    # badges
    for badge in ("b-pass", "b-warn", "b-fail"):
        assert f'class="badge {badge}"' in html
    # every section present
    for sec in (
        "flags",
        "qc",
        "lineage",
        "coverage",
        "amplicons",
        "demix",
        "moi",
        "timeseries",
        "barcode",
        "consensus",
        "provenance",
    ):
        assert f'id="sec-{sec}"' in html, sec
    assert "not available" not in html
    # run flags
    assert "Demixing failed for 1 sample(s): WW_005" in html
    # sample names and lineage names appear (stacked bars and heatmap)
    assert "WW_003.rep1" in html and "BA.2.86" in html and "Omicron" in html
    assert 'id="lin-full" class="lin" hidden' in html
    # inline svg, no external references
    assert "<svg" in html and "<polyline" in html and "<path" in html
    assert "http://" not in html and "https://" not in html
    for tag, attrs in parse(html).tags:
        assert "src" not in attrs and "href" not in attrs or tag == "a", (tag, attrs)
        assert tag not in {"iframe", "object", "embed", "img", "link"}
    # dark mode, print, palette
    assert "prefers-color-scheme: dark" in html and "@media print" in html and "--c9" in html
    # QC table is sortable and filterable
    assert 'class="sortable"' in html and 'data-f="text"' in html

    # dropped amplicons are shaded in the coverage card (WW_003.rep1 has 4)
    card = html[html.index('data-sample="WW_003.rep1" data-qc="WARN" data-ntc="0">\n      <h4>') :]
    card = card[: card.index("</svg>")]
    assert card.count('class="drop"') == 4
    # regression: `r.values` in Jinja resolves to dict.values, not the "values" key
    assert "built-in method" not in html
    assert re.search(r'class="am" data-v="[0-9.,]*"', html)


def test_failed_and_empty_runs_render_status_and_reason(tmp_path):
    html = render(build_fixture("d_all_failed", tmp_path))
    assert_valid(html)
    assert "all_samples_failed" in html and "ntc_reads" in html
    assert html.count('class="badge b-fail"') >= 3
    assert "no reads" in html  # demix reason is shown
    assert "skipped_low_coverage" in html
    # zero samples: still a valid report, no traceback
    empty = tmp_path / "empty"
    empty.mkdir()
    html = render(build_fixture("e_zero_samples", empty))
    assert_valid(html)
    assert "No samples" in html and "Traceback" not in html


def test_xss_in_sample_names_and_reasons_is_escaped(tmp_path):
    evil = "<script>alert(1)</script>"
    evil2 = '"><img src=x onerror=alert(2)>'
    qc = tmp_path / "qc.tsv"
    qc.write_text(
        "sample\tsample_type\tcollection_date\tsite\tqc_status\tqc_reasons\n"
        f"{evil}\tsample\t2026-01-01\t{evil2}\tFAIL\t<b>x</b>\n"
        f"{evil2}\tsample\t2026-02-01\t{evil2}\tPASS\t\n"
    )
    status = tmp_path / "ds.tsv"
    status.write_text(
        "sample\tstatus\treason\tbarcode_version\n"
        f"{evil}\tfailed\t{evil} reason\t{evil2}\n{evil2}\tok\t\t{evil2}\n"
    )
    agg = tmp_path / "agg.tsv"
    agg.write_text(
        "\tsummarized\tlineages\tabundances\tresid\n"
        f"{evil2}\t[('Delta<i>', 1.0)]\t\"['AY.4<u>']\"\t\"[1.0]\"\t1\n"
    )
    info = json.loads((FIX / "a_full" / "run_info.json").read_text())
    info["git_commit"] = "<svg onload=alert(3)>"
    info["freyja_data"]["barcode_version"] = evil2
    ri = tmp_path / "ri.json"
    ri.write_text(json.dumps(info))
    out = tmp_path / "r.json"
    assert (
        brd.main(
            [
                "--run-info",
                str(ri),
                "--qc",
                str(qc),
                "--demix-status",
                str(status),
                "--aggregated",
                str(agg),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    html = render(json.loads(out.read_text()))
    for needle in ("<script>alert", "<img src=x", "<svg onload", "<b>x</b>", "<i>", "<u>"):
        assert needle not in html, needle
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    # no injected elements: the only <script> is our own
    parser = parse(html)
    assert [t for t, _ in parser.tags].count("script") == 1
    assert all(not any(k.startswith("on") for k in attrs) for _, attrs in parser.tags)
    # amplicon names with quotes cannot break out of the data attribute
    data = json.loads(out.read_text())
    data["amplicons"] = {
        data["samples"][0]["sample"]: [
            {
                "amplicon": "a'\"><x onclick=1>",
                "start": 1,
                "end": 10,
                "mean_depth": 5.0,
                "dropped": True,
            }
        ]
    }
    html = render(data)
    assert all(not any(k.startswith("on") for k in attrs) for _, attrs in parse(html).tags)
    assert_valid(html)


A = FIX / "i_anonymize"
A_OPTIONS = {
    "--qc": "qc_summary.tsv",
    "--demix-status": "demix_status.tsv",
    "--aggregated": "freyja_aggregated.tsv",
    "--consensus": "consensus_stats.tsv",
    "--samplesheet": "samplesheet.csv",
    "--site-aggregates": "site_aggregates_freyja.tsv",
    "--site-aggregates-meta": "site_aggregates_meta.tsv",
}
A_IDS = ["NTC_Central_01", "WW_Central_001", "WW_Central_002", "WW_North_001", "WW_North_002.rep1"]
A_NAMES = ["Central_WWTP", "North_Plant", "Metro_Region", "Central", "North_", "WW_", "NTC_"]
A_PATH_LEAKS = [
    "jdoe",
    "j.doe",
    "example.org",
    "/Users",
    "/home",
    "C:\\",
    "Run 1 Plate A",
    "My Results",
    "~/",
]


def test_full_report_html_contains_no_identifying_string(tmp_path):
    args = ["--run-info", str(A / "run_info.json")]
    for option, name in A_OPTIONS.items():
        args += [option, str(A / name)]
    out = tmp_path / "full.json"
    assert brd.main([*args, "--anonymize", "full", "--out", str(out)]) == 0
    html = render(json.loads(out.read_text()))
    assert_valid(html)
    for leak in A_IDS + A_NAMES + A_PATH_LEAKS:
        assert leak not in html, leak
    for pseudo in ("S001", "S005", "SITE1", "SITE2", "AGG1"):
        assert pseudo in html, pseudo
