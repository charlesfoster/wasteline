import json
import random
import re
import time
from html.parser import HTMLParser
from pathlib import Path

import build_report_data as brd
import pytest
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


# --------------------------------------------------------------------------- (a) full
@pytest.fixture(scope="module")
def full_html(tmp_path_factory):
    data = build_fixture("a_full", tmp_path_factory.mktemp("a"), depths=True)
    return render(data)


def test_full_run_content(full_html):
    html = full_html
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


def test_coverage_marks_dropped_amplicons(full_html):
    # WW_003.rep1 has 4 dropped amplicons -> 4 shaded rects in its card
    card = full_html[
        full_html.index('data-sample="WW_003.rep1" data-qc="WARN" data-ntc="0">\n      <h4>') :
    ]
    card = card[: card.index("</svg>")]
    assert card.count('class="drop"') == 4
    assert 'class="thr"' in card


def test_amplicon_strip_data_attribute_is_the_depth_list(full_html):
    # `r.values` in Jinja resolves to the dict.values method, not the "values" key
    assert "built-in method" not in full_html
    assert re.search(r'class="am" data-v="[0-9.,]*"', full_html)


def test_amplicon_dropout_view(full_html):
    assert "data-amps=" in full_html and "midnight_1" in full_html
    assert full_html.count('class="am"') == 6
    assert "dropout frequency" in full_html


def test_mutations_and_time_series_and_comparison(full_html):
    assert "S:N501Y" in full_html or "N501Y" in full_html
    assert 'class="h9"' in full_html or 'class="h8"' in full_html
    assert "SiteA" in full_html and "SiteB" in full_html
    assert "08_01_2026-00-10" in full_html  # comparison label B
    assert "Largest lineage changes" in full_html


def test_template_autoescape_enabled():
    assert rr.make_environment(TEMPLATE).autoescape is True


# --------------------------------------------------------------------------- (b) demix only
def test_demix_only_omits_sections(tmp_path):
    html = render(build_fixture("b_demix_only", tmp_path))
    assert_valid(html)
    assert "not available in mode &#39;demix&#39;" in html
    for sec in ("moi", "timeseries", "barcode"):
        assert f'id="sec-{sec}"' not in html, sec
    assert (
        'id="qc-table"' not in html and 'id="cov-grid"' not in html and 'id="amp-rows"' not in html
    )
    assert "09_17_2026-00-46" in html and "Delta" in html
    assert "previous_run_barcode_unknown" in html
    assert 'class="badge b-fail"' in html  # S2 failed
    assert "https://" not in html and "http://" not in html


# --------------------------------------------------------------------------- (c) legacy barcode
def test_legacy_unknown_barcode_is_prominent(tmp_path):
    html = render(build_fixture("c_legacy_barcode", tmp_path))
    assert_valid(html)
    header = html[: html.index('id="sec-flags"')]
    assert 'class="barcode unknown"' in header and "<b>unknown</b>" in header
    assert "barcode_unknown" in html
    assert "legacy" in header.lower()


# --------------------------------------------------------------------------- (d) all failed
def test_all_failed_samples(tmp_path):
    html = render(build_fixture("d_all_failed", tmp_path))
    assert_valid(html)
    assert "all_samples_failed" in html and "ntc_reads" in html
    assert "No sample has a successful lineage result" in html
    assert html.count('class="badge b-fail"') >= 3
    assert "no reads" in html  # demix reason is shown
    assert "skipped_low_coverage" in html


# --------------------------------------------------------------------------- (e) zero samples
def test_zero_samples(tmp_path):
    html = render(build_fixture("e_zero_samples", tmp_path))
    assert_valid(html)
    assert "No samples" in html and "no_samples" in html
    assert "Traceback" not in html


# --------------------------------------------------------------------------- robustness
@pytest.mark.parametrize(
    "data",
    [
        {},
        {"samples": None, "demix": [], "provenance": "x"},
        {"samples": [{"sample": "x"}]},
        {
            "samples": [{"sample": "x", "qc_status": "PASS", "mean_depth": "NaN"}],
            "demix": {"x": {"status": "ok", "groups": "bad"}},
        },
        {"time_series": {"sites": [{"dates": ["2026-01-01"], "series": []}]}},
    ],
)
def test_minimal_and_malformed_data_still_render(data):
    assert_valid(render(data))


# --------------------------------------------------------------------------- escaping
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


# --------------------------------------------------------------------------- scale
def synthetic_data(n=150, seed=7):
    rng = random.Random(seed)
    groups = ["Omicron", "Delta", "Alpha", "Gamma", "Beta", "Other"]
    lineages = [f"BA.{i}" for i in range(1, 30)]
    qc, status, agg, amps, tracks = {}, {}, {}, {}, {}
    n_bins = 600
    for i in range(n):
        sid = f"NTC_{i:03d}" if i % 50 == 0 else f"WW_{i:03d}_site{i % 5}"
        neg = sid.startswith("NTC")
        rec = brd.blank_sample(sid)
        rec.update(
            sample_type="negative" if neg else "sample",
            collection_date=f"2026-0{1 + i % 9}-1{i % 9}",
            site=f"Site{i % 5}",
            raw_reads=rng.randint(1000, 500000),
            mean_depth=rng.uniform(0, 900),
            pct_genome_ge_10x=rng.uniform(40, 100),
            qc_status=rng.choice(["PASS", "PASS", "WARN", "FAIL"]),
            amplicons_total=98,
            amplicons_dropped=rng.randint(0, 6),
            consensus_n_pct=rng.uniform(0, 5),
        )
        qc[sid] = rec
        status[sid] = {"status": "ok", "reason": "", "barcode_version": "09_17_2026-00-46"}
        fr = [rng.random() for _ in groups]
        tot = sum(fr)
        agg[sid] = {
            "groups": sorted(
                [[g, round(f / tot, 5)] for g, f in zip(groups, fr, strict=True)],
                key=lambda p: -p[1],
            ),
            "lineages": [[lin, round(1 / len(lineages), 5)] for lin in lineages[:18]],
            "resid": 3.0,
            "coverage": 98.0,
        }
        amps[sid] = [
            {
                "amplicon": f"scheme_{k}",
                "start": (k - 1) * 300 + 1,
                "end": (k - 1) * 300 + 400,
                "mean_depth": 0.0 if rng.random() < 0.05 else round(rng.uniform(5, 900), 2),
                "dropped": False,
            }
            for k in range(1, 99)
        ]
        for a in amps[sid]:
            a["dropped"] = a["mean_depth"] < 10
        tracks[sid] = [int(rng.uniform(0, 1500)) for _ in range(n_bins)]
    info = json.loads((FIX / "a_full" / "run_info.json").read_text())
    data = brd.build(
        info,
        qc=qc,
        amplicons=amps,
        aggregated=agg,
        demix_status=status,
        consensus=[
            {
                "sample": s,
                "length": 29903,
                "n_count": 5,
                "n_pct": 0.02,
                "ambiguous_count": 0,
                "status": "ok",
            }
            for s in qc
        ],
        depths={
            "contig": "NC_045512.2",
            "genome_length": 29903,
            "bin_size": 50,
            "n_bins": n_bins,
            "samples": tracks,
        },
        moi=[
            {
                "sample": s,
                "name": f"m{k}",
                "mutation": f"S:X{k}Y",
                "ref_depth": 5,
                "alt_depth": 5,
                "depth": 10,
                "alt_freq": 0.5,
            }
            for s in list(qc)[:150]
            for k in range(8)
        ],
    )
    return data


def test_150_samples_size_and_speed():
    data = synthetic_data(150)
    assert len(data["samples"]) == 150
    start = time.perf_counter()
    html = render(data)
    elapsed = time.perf_counter() - start
    size = len(html.encode())
    assert_valid(html)
    assert size < 3_000_000, f"report is {size / 1e6:.2f} MB"
    assert elapsed < 10, f"render took {elapsed:.1f}s"
    assert html.count("<tr data-sample=") >= 150
    assert "https://" not in html and "http://" not in html


# --------------------------------------------------------------------------- views / CLI
def test_lineage_view_merges_small_segments_and_uses_palette():
    samples = [{"sample": "a", "sample_type": "sample", "qc_status": "PASS"}]
    pairs = [[f"L{i}", 0.01] for i in range(40)] + [["Big", 0.6]]
    demix = {"a": {"status": "ok", "lineages": sorted(pairs, key=lambda p: -p[1])}}
    view = rr.lineage_view(samples, demix, "lineages")
    classes = {s["cls"] for s in view["rows"][0]["segs"]}
    assert "c0" in classes and "co" in classes
    assert all(s["w"] >= 0 for s in view["rows"][0]["segs"])
    assert view["legend"][0]["name"] == "Big"


def test_coverage_view_uses_one_based_inclusive_coordinates():
    data = rr.normalise(
        {
            "samples": [{"sample": "s"}],
            "coverage_tracks": {
                "genome_length": 1000,
                "bin_size": 10,
                "n_bins": 100,
                "samples": {"s": [100] * 100},
            },
            "amplicons": {
                "s": [{"amplicon": "a", "start": 11, "end": 30, "mean_depth": 0, "dropped": True}]
            },
        }
    )
    drop = rr.coverage_view(data)["cards"][0]["drops"][0]
    assert drop["x"] == 1.0 and drop["w"] == 2.0  # positions 11..30 = bins 1 and 2


def test_cli_roundtrip(tmp_path):
    data = build_fixture("a_full", tmp_path)
    data_path = tmp_path / "d.json"
    data_path.write_text(json.dumps(data))
    out = tmp_path / "report.html"
    assert rr.main(["--data", str(data_path), "--template", str(TEMPLATE), "--out", str(out)]) == 0
    assert_valid(out.read_text())
    # coverage tracks absent (no depths): note, not an error
    assert "Per-position depth was not available" in out.read_text()
    assert rr.main(["--data", str(tmp_path / "nope.json"), "--out", str(out)]) == 2


def test_cli_bad_template(tmp_path, capsys):
    (tmp_path / "t.j2").write_text("{{ broken ")
    (tmp_path / "d.json").write_text("{}")
    assert (
        rr.main(
            [
                "--data",
                str(tmp_path / "d.json"),
                "--template",
                str(tmp_path / "t.j2"),
                "--out",
                str(tmp_path / "o.html"),
            ]
        )
        == 2
    )
    assert "template error" in capsys.readouterr().err


def test_negative_controls_block_and_boot_interval_tooltip():
    data = {
        "schema_version": 1,
        "mode": "full",
        "samples": [
            {"sample": "S1", "sample_type": "sample", "qc_status": "PASS"},
            {
                "sample": "NEG_C",
                "sample_type": "negative",
                "qc_status": "FAIL",
                "qc_reasons": ["ntc_reads"],
            },
        ],
        "demix": {
            "S1": {
                "status": "ok",
                "groups": [["Omicron", 1.0]],
                "lineages": [["BA.2", 1.0]],
                "intervals": {"BA.2": [0.95, 1.0]},
            },
            "NEG_C": {"status": "ok", "groups": [], "lineages": []},
        },
        "available": {"demix": True},
        "negative_controls": [
            {
                "sample": "NEG_C",
                "status": "ok",
                "qc_status": "FAIL",
                "qc_reasons": ["ntc_reads"],
                "mapped_reads": 5998,
                "pct_genome_ge_1x": 99.5,
                "lineages": [["BA.2", 1.0]],
                "coverage": 90.0,
                "excluded_from_aggregate": True,
            }
        ],
    }
    html = rr.render(data)
    assert 'id="sec-controls"' in html and "excluded from the aggregate" in html
    assert 'href="#sec-controls"' in html  # in the table of contents
    assert "BA.2 100.0%" in html and "ntc_reads" in html
    assert "95% interval 95.0" in html
    without = rr.render({**data, "negative_controls": []})
    assert 'id="sec-controls"' not in without and 'href="#sec-controls"' not in without
