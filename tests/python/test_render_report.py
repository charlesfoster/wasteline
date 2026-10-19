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


def test_capture_protocol_header_and_amplicon_note():
    data = {
        "schema_version": 1,
        "mode": "full",
        "header": {"protocol": "capture", "duplicates": "mark"},
        "samples": [{"sample": "S1", "sample_type": "sample", "qc_status": "PASS"}],
        "demix": {
            "S1": {"status": "ok", "groups": [["Omicron", 1.0]], "lineages": [["BA.2", 1.0]]}
        },
        "available": {"demix": True, "amplicons": False},
    }
    html = rr.render(data)
    assert "Library protocol" in html and "capture (no primer trimming, duplicates marked)" in html
    assert "Not applicable for protocol <b>capture</b>" in html
    amplicon = rr.render({**data, "header": {"protocol": "amplicon"}})
    assert (
        "Not applicable for protocol" not in amplicon
        and "Amplicon coverage is not available" in amplicon
    )


# =============================================================================================
# Primer scheme columns, dropout impact, unexplained mutations, time series + growth rates
# =============================================================================================
G = FIX / "g_new_analyses"
G_OPTIONS = {
    "--qc": "qc_summary.tsv",
    "--aggregated": "freyja_aggregated.tsv",
    "--samplesheet": "samplesheet.csv",
    "--dropout-impact": "dropout_impact.tsv",
    "--dropout-summary": "dropout_impact_summary.tsv",
    "--unexplained": "unexplained_mutations.tsv",
    "--timeseries": "freyja_timeseries.tsv",
    "--timeseries-lineages": "freyja_timeseries_lineages.tsv",
    "--growth-rates": "freyja_growth_rates.csv",
    "--growth-inputs": "freyja_growth_inputs.json",
}


UE_HEADER = (
    "sample\tposition\tref\talt\tobserved_freq\texpected_freq\texcess\tdepth\talt_depth"
    "\tgene\taa_change\teffect\tn_samples\trecurrent\tp_value"
)


def build_g(tmp_path, drop=(), extra=None, run_info=None):
    args = ["--run-info", str(run_info or G / "run_info.json")]
    for option, name in G_OPTIONS.items():
        if option not in drop:
            args += [option, str(G / name)]
    args += list(extra or []) + ["--out", str(tmp_path / "g.json")]
    assert brd.main(args) == 0
    return json.loads((tmp_path / "g.json").read_text())


@pytest.fixture(scope="module")
def g_html(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("g")
    return render(build_g(tmp, extra=["--dashboard-name", "freyja_dashboard.html"]))


def test_all_new_sections_render(g_html):
    html = g_html
    assert_valid(html)
    for sec in ("dropout-impact", "unexplained", "timeseries"):
        assert f'id="sec-{sec}"' in html, sec
        assert f'href="#sec-{sec}"' in html, sec  # nav entry
    assert html.count('id="sec-timeseries"') == 1
    assert "http://" not in html and "https://" not in html
    assert html.count("<script>") == 2  # the original script plus the one for the new controls


def test_default_report_has_no_new_markup(tmp_path):
    html = render(build_fixture("a_full", tmp_path))
    for needle in ("sec-dropout-impact", "sec-unexplained", 'id="ts-root"', "b-info", "Primer scheme",
                   "showall", "growth-rates", "rtip"):  # fmt: skip
        assert needle not in html, needle
    assert html.count("<script>") == 1 and html.count("<style>") == 1


def test_qc_table_primer_columns_badges_and_tooltips(g_html):
    assert "Primer trimmed %" in g_html and "Primer scheme" in g_html
    assert "informational" in g_html  # tooltip of the primer-trimmed column
    # mismatch: amber badge, best match in the tooltip, reason badge with an explanation
    assert 'class="badge b-warn" title="best match: articv4.1 (declared: midnight)' in g_html
    assert (
        'class="rtip" title="Read starts line up much better with another primer scheme' in g_html
    )
    # ok: green; uninformative: grey with the tagmentation explanation
    assert (
        'class="badge b-pass" title="Read starts line up with the declared primer scheme' in g_html
    )
    assert (
        "no scheme check possible" in g_html
        and 'class="badge b-na" title="Read starts do not line up' in g_html
    )
    assert "primer_scheme_mismatch" in g_html  # run flag code and reason token
    # NA rows have no badge at all: LOWCOV_1 / NEG_1
    row = g_html[g_html.index('data-sample="NEG_1" data-qc="PASS"') :]
    row = row[: row.index("</tr>")]
    assert '<td data-v=""></td>' in row


def test_low_primer_trimming_reason_gets_a_tooltip(tmp_path):
    qc = tmp_path / "qc.tsv"
    qc.write_text(
        "sample\tqc_status\tqc_reasons\tprimer_trimmed_pct\nA\tWARN\tlow_coverage;low_primer_trimming\t0.5\n"
    )
    out = tmp_path / "r.json"
    assert (
        brd.main(["--run-info", str(G / "run_info.json"), "--qc", str(qc), "--out", str(out)]) == 0
    )
    html = render(json.loads(out.read_text()))
    assert_valid(html)
    assert (
        'low_coverage; <span class="rtip" title="Very few reads had primer sequence trimmed' in html
    )
    assert "wrong primer scheme" in html
    assert "Primer scheme</th>" not in html  # no check column: only the trimmed % one
    assert "Primer trimmed %" in html


def test_old_qc_table_shows_no_primer_columns(tmp_path):
    html = render(build_fixture("a_full", tmp_path))
    assert "Primer trimmed" not in html


def test_dropout_impact_section(g_html):
    sec = g_html[g_html.index('<section id="sec-dropout-impact">') :]
    sec = sec[: sec.index("</section>")]
    assert "Amplicon dropout impact" in sec and "lineage-informative" in sec
    assert sec.count('class="hl') >= 4  # two summary rows + two table rows highlighted
    assert "affects detected lineages" in sec and "XBB.1.5 (9), BA.2 (3), BA.5 (3)" in sec
    assert "21,533&ndash;22,612" in sec and "nCoV-2019_22" in sec
    assert "no lineage estimate" in sec  # LOWCOV_1
    assert "850 / 8,006" in sec and "10.62" in sec
    assert "28 sample(s) without dropped amplicons are not listed" in sec
    assert 'id="dx-table"' in sec and 'data-scope="#dx-table"' in sec  # filterable
    # affecting rows come first in the amplicon table
    table = sec[sec.index('id="dx-table"') :]
    assert table.index("SiteA_W03_1") < table.index("LOWCOV_1")


def test_dropout_flag_is_rendered(g_html):
    assert "dropout_affects_lineages" in g_html
    assert (
        "Amplicon dropout hides defining mutations (those that separate a lineage from its parent)"
        in g_html
    )


def test_unexplained_section(g_html):
    sec = g_html[g_html.index('<section id="sec-unexplained">') :]
    sec = sec[: sec.index("</section>")]
    assert "<b>It is not a variant caller</b>" in sec and "p-value \u2264 1.67e-06" in sec
    assert sec.count("recurrent (3)") == 3 and "single sample" in sec
    assert "2.39e-42" in sec and "E180D" in sec and "R177*" in sec
    assert 'class="hl' in sec
    # recurrent rows first (file order is kept)
    assert sec.index("SiteC_W05_1") < sec.index("SiteB_W06_1")
    assert "unexplained_mutations" in g_html  # info flag


def test_unexplained_empty_and_absent(tmp_path):
    header = UE_HEADER
    empty = tmp_path / "ue.tsv"
    empty.write_text(header + "\n")
    html = render(build_g(tmp_path, drop=("--unexplained",), extra=["--unexplained", str(empty)]))
    assert_valid(html)
    assert 'id="sec-unexplained"' in html and "None found with the current thresholds" in html
    assert 'id="ue-table"' not in html
    absent = render(build_g(tmp_path, drop=("--unexplained",)))
    assert 'id="sec-unexplained"' not in absent and 'href="#sec-unexplained"' not in absent


def test_timeseries_views_toggles_and_notes(g_html):
    sec = g_html[g_html.index('<section id="sec-timeseries"') :]
    sec = sec[: sec.index("</section>")]
    assert sec.count('class="tsview"') == 4  # summarized/full x unweighted/weighted
    assert sec.count("hidden>") + sec.count(" hidden") >= 3
    assert 'data-v="w" aria-pressed="false">Viral-load weighted' in sec
    assert "disabled" not in sec  # loads exist, lineages exist
    assert 'id="ts-site"' in sec and "All sites (pooled)" in sec and "SiteB" in sec
    assert "(n=" in sec and "mean viral load" in sec  # tooltips carry n and the load
    assert sec.count('data-site="__all__"') == 4


def test_timeseries_weighting_disabled_without_viral_load(tmp_path):
    tsv = tmp_path / "ts.tsv"
    rows = [
        "site\tcollection_date\tlineage_group\tabundance\tabundance_load_weighted\tn_samples\tmean_viral_load\ttotal_viral_load"
    ]
    for d, v in (("2026-01-05", 0.2), ("2026-01-12", 0.5), ("2026-01-19", 0.7)):
        rows += [f"A\t{d}\tDelta\t{v}\tNA\t1\tNA\tNA", f"A\t{d}\tAlpha\t{1 - v}\tNA\t1\tNA\tNA"]
    tsv.write_text("\n".join(rows) + "\n")
    html = render(
        build_g(
            tmp_path,
            drop=("--timeseries", "--timeseries-lineages"),
            extra=["--timeseries", str(tsv)],
        )
    )
    assert_valid(html)
    sec = html[html.index('<section id="sec-timeseries"') :]
    assert 'disabled title="No viral load recorded' in sec
    assert 'disabled title="No full-lineage time series' in sec
    assert sec.count('class="tsview"') == 1
    assert "All sites (pooled)" not in sec  # a single site


def test_timeseries_falls_back_to_the_derived_panel(tmp_path):
    html = render(
        build_g(
            tmp_path,
            drop=("--timeseries", "--timeseries-lineages", "--growth-rates", "--growth-inputs"),
        )
    )
    assert_valid(html)
    sec = html[html.index('<section id="sec-timeseries"') :]
    assert 'id="ts-root"' not in html and 'class="tsview"' not in html
    assert "One panel per site, shared axes." in sec  # today's panel


def test_timeseries_with_a_single_date_draws_no_chart(tmp_path):
    tsv = tmp_path / "ts.tsv"
    tsv.write_text(
        "site\tcollection_date\tlineage_group\tabundance\tabundance_load_weighted\tn_samples\tmean_viral_load\ttotal_viral_load\n"
        "A\t2026-01-05\tDelta\t1.0\tNA\t1\tNA\tNA\n"
    )
    data = build_g(
        tmp_path, drop=("--timeseries", "--timeseries-lineages"), extra=["--timeseries", str(tsv)]
    )
    assert rr.timeseries_view(rr.normalise(data), {"summarized": {}, "lineages": {}}) is None
    assert_valid(render(data))


def test_growth_rate_table_and_notes(g_html):
    sec = g_html[g_html.index('id="growth-rates"') :]
    sec = sec[: sec.index("</section>")]
    assert "XBB.1.5" in sec and "+60.3" in sec and "-35.1" in sec and "[44.74% , 49.89%]" in sec
    assert (
        sec.index("XBB.1.5") < sec.index("JN.1") < sec.index("BA.1") < sec.index("BA.2")
    )  # by advantage
    assert "How to read this" in sec and "not as a confidence interval" in sec
    assert (
        "5.5 days" in sec
        and "Fitted on 30 sample(s) over 8 collection date(s) spanning 49 days" in sec
    )
    # dashboard: a file name as text, never a link, and the internet caveat
    assert "<code>freyja_dashboard.html</code>" in sec and "needs internet access" in sec
    assert 'href="freyja_dashboard' not in g_html


def test_growth_failed_notice_is_shown(tmp_path):
    failed = tmp_path / "f.txt"
    failed.write_text(
        "freyja relgrowthrate produced no growth-rate table: insufficient data: 3 distinct collection date(s) spanning 14 day(s)\n"
    )
    html = render(
        build_g(
            tmp_path,
            drop=("--growth-rates", "--growth-inputs"),
            extra=["--growth-failed", str(failed)],
        )
    )
    assert_valid(html)
    assert (
        "Growth rates could not be computed." in html
        and "insufficient data: 3 distinct collection date(s)" in html
    )
    assert 'id="growth-table"' not in html


def test_growth_without_a_time_series_table_still_renders_a_section(tmp_path):
    data = build_g(tmp_path, drop=("--timeseries", "--timeseries-lineages", "--samplesheet"))
    data["time_series"] = None
    html = render(data)
    assert_valid(html)
    assert html.count('id="sec-timeseries"') == 1 and 'id="growth-table"' in html
    assert 'id="ts-root"' not in html and "One panel per site" not in html


def test_long_tables_are_capped_with_show_all(tmp_path):
    header = UE_HEADER
    rows = [
        f"SiteA_W01_1\t{1000 + i}\tA\tT\t0.2\t0\t0.2\t300\t60\tS\tX{i}Y\tmissense\t1\tfalse\t1e-9"
        for i in range(1700)
    ]
    big = tmp_path / "ue.tsv"
    big.write_text(header + "\n" + "\n".join(rows) + "\n")
    html = render(build_g(tmp_path, drop=("--unexplained",), extra=["--unexplained", str(big)]))
    assert_valid(html)
    sec = html[html.index('<section id="sec-unexplained">') :]
    sec = sec[: sec.index("</section>")]
    assert sec.count("<tr data-sample=") == rr.MAX_TABLE_ROWS
    assert sec.count(' more"') == rr.MAX_TABLE_ROWS - rr.VISIBLE_ROWS
    assert (
        f"Show all {rr.MAX_TABLE_ROWS:,} rows" in sec
        and "of 1,700" in sec
        and "unexplained_mutations.tsv" in sec
    )
    assert len(html) < 3_000_000


def test_xss_in_the_new_tables_is_escaped(tmp_path):
    evil, evil2 = "<script>alert(1)</script>", '"><img src=x onerror=alert(2)>'
    qc = tmp_path / "qc.tsv"
    qc.write_text(
        "sample\tsample_type\tcollection_date\tsite\tqc_status\tqc_reasons\tprimer_trimmed_pct\tprimer_scheme_check\tprimer_scheme_best\n"
        f"{evil}\tsample\t2026-01-01\t{evil2}\tWARN\tprimer_scheme_mismatch\t5\tmismatch\t{evil2}\n"
        f"{evil2}\tsample\t2026-02-01\t{evil}\tPASS\t\t5\tok\tNA\n"
    )
    agg = tmp_path / "agg.tsv"
    agg.write_text(
        "\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
        f"{evil}\t[('Delta<i>', 1.0)]\tAY.4\t1.0\t1\t99\n{evil2}\t[('Delta<i>', 1.0)]\tAY.4\t1.0\t1\t99\n"
    )
    imp = tmp_path / "imp.tsv"
    imp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tn_informative_sites\tn_sites_in_detected_lineages\taffects_detected_lineage\taffected_lineages\n"
        f"{evil}\t{evil2}\t1\t100\t0\t5\t2\ttrue\t{evil2}:2\n"
    )
    ue = tmp_path / "ue.tsv"
    ue.write_text(
        "sample\tposition\tref\talt\tobserved_freq\texpected_freq\texcess\tdepth\talt_depth\tgene\taa_change\teffect\tn_samples\trecurrent\tp_value\n"
        f"{evil}\t5\tA\tT\t0.3\t0\t0.3\t100\t30\t{evil2}\t{evil}\t{evil2}\t2\ttrue\t1e-9\n"
    )
    tsv = tmp_path / "ts.tsv"
    tsv.write_text(
        "site\tcollection_date\tlineage_group\tabundance\tabundance_load_weighted\tn_samples\tmean_viral_load\ttotal_viral_load\n"
        f"x{evil2}\t2026-01-01\t{evil}\t0.4\t0.4\t1\t10\t10\nx{evil2}\t2026-02-01\t{evil}\t0.6\t0.6\t1\t10\t10\n"
    )
    growth = tmp_path / "g.csv"
    growth.write_text(f'Lineage,Estimated Advantage,Bootstrap 95% interval\n{evil},5%,"{evil2}"\n')
    failed = tmp_path / "f.txt"
    failed.write_text(evil)
    args = ["--run-info", str(G / "run_info.json"), "--qc", str(qc), "--aggregated", str(agg),
            "--dropout-impact", str(imp), "--unexplained", str(ue), "--timeseries", str(tsv),
            "--growth-rates", str(growth), "--dashboard-name", evil2, "--out", str(tmp_path / "r.json")]  # fmt: skip
    assert brd.main(args) == 0
    data = json.loads((tmp_path / "r.json").read_text())
    html = render(data)
    failed_html = render(dict(data, growth_rates={"status": "failed", "reason": evil, "rows": []}))
    for page in (html, failed_html):
        assert_valid(page)
        for needle in ("<script>alert", "<img src=x", "Delta<i>"):
            assert needle not in page, needle
        parser = parse(page)
        assert [t for t, _ in parser.tags].count("script") == 2
        assert all(not any(k.startswith("on") for k in attrs) for _, attrs in parser.tags)
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    assert (
        'id="sec-unexplained"' in html
        and 'id="sec-dropout-impact"' in html
        and 'id="ts-root"' in html
    )


def test_150_samples_with_all_new_tables_stay_small_and_fast():
    rng = random.Random(3)
    data = synthetic_data(150)
    ids = [s["sample"] for s in data["samples"]]
    real = [s for s in ids if not s.startswith("NTC")]
    impact = {}
    for sid in ids:
        amps = [
            {"amplicon": f"scheme_{k}", "start": k * 300, "end": k * 300 + 400, "mean_depth": 1.5,
             "n_informative_sites": 40, "n_sites_in_detected_lineages": 5, "affects_detected_lineage": k % 2 == 0,
             "affected_lineages": [{"lineage": "BA.2", "n_sites": 3}, {"lineage": "XBB.1.5", "n_sites": 2}]}
            for k in range(1, 7)
        ]  # fmt: skip
        impact[sid] = {
            "summary": {"n_dropped_amplicons": 6, "informative_sites_total": 8006, "informative_sites_in_dropped": 240,
                        "pct_informative_sites_lost": 3.0, "n_dropped_affecting_detected_lineages": 3},
            "amplicons": amps,
        }  # fmt: skip
    data["dropout_impact"] = impact
    data["unexplained_mutations"] = [
        {"sample": sid, "position": 1000 + k, "ref": "A", "alt": "T", "observed_freq": 0.2, "expected_freq": 0.0,
         "excess": 0.2, "depth": 300, "alt_depth": 60, "gene": "S", "aa_change": f"X{k}Y", "effect": "missense",
         "n_samples": 1, "recurrent": False, "p_value": 1e-9}
        for sid in real for k in range(8)
    ]  # fmt: skip
    dates = [f"2026-{m:02d}-{d:02d}" for m in (1, 2, 3) for d in (3, 10, 17, 24)]
    groups = ["Omicron", "Delta", "Alpha", "Gamma", "Beta", "XBB", "BA.2"]
    rows = [
        {"site": f"Site{s}", "collection_date": dt, "lineage_group": g, "abundance": rng.random() / 3,
         "abundance_load_weighted": rng.random() / 3, "n_samples": 3, "mean_viral_load": 1e5, "total_viral_load": 3e5}
        for s in range(5) for dt in dates for g in groups
    ]  # fmt: skip
    data["timeseries"] = {"enabled": True, "has_viral_load": True, "min_coverage": 60, "sites": [f"Site{s}" for s in range(5)],
                          "dates": dates, "summarized": rows, "lineages": rows, "lineages_total_groups": None}  # fmt: skip
    data["growth_rates"] = {"status": "ok", "reason": None, "inputs": {}, "params": {},
                            "rows": [{"lineage": g, "advantage_pct": 5.0, "interval_text": "[1% , 2%]"} for g in groups]}  # fmt: skip
    start = time.perf_counter()
    html = render(data)
    elapsed = time.perf_counter() - start
    assert_valid(html)
    assert len(html.encode()) < 3_000_000, len(html)
    assert elapsed < 10
    assert html.count("tsview") >= 4


# --------------------------------------------------------------------------- consensus reason column
def test_consensus_reason_column_and_tooltip(tmp_path):
    html = render(build_g(tmp_path, extra=["--consensus", str(G / "consensus_stats.tsv")]))
    assert_valid(html)
    sec = html[html.index('<section id="sec-consensus">') :]
    sec = sec[: sec.index("</section>")]
    assert "<th data-sort>Reason</th>" in sec
    assert (
        "bcftools consensus failed: The fasta sequence does not match the REF allele at NC_045512.2:11287:"
        in sec
    )
    assert 'class="badge b-fail" title="bcftools consensus failed' in sec
    assert (
        'class="badge b-pass" title="skipped variant at NC_045512.2:22942 (overlaps another variant)"'
        in sec
    )
    assert "consensus_failed" in html and "consensus_warnings" in html
    assert "consensus genome failed for 1 sample(s): LOWCOV_1" in html


def test_consensus_without_reasons_has_no_reason_column(tmp_path):
    html = render(build_g(tmp_path, extra=["--consensus", str(G / "consensus_stats_6col.tsv")]))
    sec = html[html.index('<section id="sec-consensus">') :]
    assert "<th data-sort>Reason</th>" not in sec[: sec.index("</section>")]
    assert "consensus_failed" in html  # the flag does not depend on a reason
    plain = render(build_fixture("a_full", tmp_path))
    assert (
        "<th data-sort>Reason</th>"
        not in plain[plain.index('<section id="sec-consensus">') :].split("</section>")[0]
    )


def test_consensus_reason_is_escaped(tmp_path):
    cons = tmp_path / "c.tsv"
    cons.write_text(
        "sample\tlength\tn_count\tn_pct\tambiguous_count\tstatus\treason\n"
        'LOWCOV_1\t29903\t29903\t100\t0\tfailed\t"><img src=x onerror=alert(1)><script>alert(2)</script>\n'
    )
    html = render(build_g(tmp_path, extra=["--consensus", str(cons)]))
    assert_valid(html)
    assert "<img src=x" not in html and "<script>alert(2)" not in html
    assert all(not any(k.startswith("on") for k in attrs) for _, attrs in parse(html).tags)


# =============================================================================================
# Lineage groups and site aggregates
# =============================================================================================
H = FIX / "h_lineage_groups"


def build_h(tmp_path, drop=(), extra=None, only=None):
    options = {
        "--aggregated": "freyja_aggregated.tsv",
        "--demix-status": "demix_status.tsv",
        "--samplesheet": "samplesheet.csv",
        "--lineage-summary": "lineage_summary.csv",
        "--lineage-groups-resolved": "lineage_groups_resolved.tsv",
        "--site-aggregates": "site_aggregates_freyja.tsv",
        "--site-aggregates-meta": "site_aggregates_meta.tsv",
    }
    args = ["--run-info", str(H / "run_info.json")]
    for option, name in options.items():
        if option not in drop and (only is None or option in only):
            args += [option, str(H / name)]
    args += list(extra or []) + ["--out", str(tmp_path / "h.json")]
    assert brd.main(args) == 0
    return json.loads((tmp_path / "h.json").read_text())


def section(html, sec_id):
    part = html[html.index(f'<section id="{sec_id}"') :]
    return part[: part.index("</section>")]


@pytest.fixture(scope="module")
def h_html(tmp_path_factory):
    return render(build_h(tmp_path_factory.mktemp("h")))


def test_site_aggregates_section(h_html):
    sec = section(h_html, "sec-site-aggregates")
    assert 'href="#sec-site-aggregates"' in h_html
    assert "never from pooled reads" in sec and "weights" in sec.lower()
    assert sec.count('class="sarow"') == 10  # 5 periods x (summarized + full view)
    assert 'id="sa-summ"' in sec and 'id="sa-full" class="sa" hidden' in sec
    assert 'data-show="#sa-full" data-hide="#sa-summ"' in sec  # same toggle as the per-sample bars
    # grouped by aggregate, periods in order
    assert sec.index("Catchment_X") < sec.index("Coastal")
    first_view = sec[: sec.index('id="sa-full"')]
    assert (
        first_view.index("2026-W02") < first_view.index("2026-W03") < first_view.index("2026-W04")
    )
    # badges: contributing sites, weighting basis, weight share, missing sites with reasons
    assert ">2/2 sites<" in sec and ">1/2 sites<" in sec
    assert (
        'class="chip2 weak" title="Weighting basis of the average (weaker than requested)">equal (viral_load missing)'
        in sec
    )
    assert ">50% of weight<" in sec and ">64% of weight<" in sec
    assert 'class="chip2 miss" title="Site_2: not_ok:failed">missing Site_2' in sec
    assert 'title="Site_3: no_sample"' in sec
    assert "<title>NB.1.8: 100.0%</title>" in sec
    assert "site_aggregate_incomplete" in h_html


def test_site_aggregate_bars_reuse_the_per_sample_colours(h_html):
    per_sample = section(h_html, "sec-lineage")

    def legend(part):
        return {
            m.split("</li>")[0].split("</i>")[1]: m.split('class="sw ')[1].split('"')[0]
            for m in part.split("<li><i ")[1:]
        }

    sa = section(h_html, "sec-site-aggregates")
    full_sa = legend(sa[sa.index('id="sa-full"') :])
    full_samples = legend(per_sample[per_sample.index('id="lin-full"') :].split("<h3>")[0])
    shared = set(full_sa) & set(full_samples)
    assert shared and all(full_sa[n] == full_samples[n] for n in shared)


def test_period_without_a_solution_shows_a_status_not_a_bar(tmp_path):
    old = ROOT / "tests" / "fixtures" / "lineage_groups"
    data = build_h(tmp_path, drop=("--site-aggregates", "--site-aggregates-meta"))
    data["site_aggregates"] = brd.load_site_aggregates(
        old / "site_aggregates_freyja.tsv", old / "site_aggregates_meta.tsv", brd.Problems()
    )
    html = render(data)
    assert_valid(html)
    sec = section(html, "sec-site-aggregates")
    assert (
        "no site contributed a result in this period" in sec
        and ">0/2 sites<" in sec
        and 'class="badge b-fail" title="Sites that contributed: none"' in sec
    )


def test_default_report_has_no_lineage_group_markup(tmp_path):
    html = render(build_fixture("a_full", tmp_path))
    for needle in ("sec-lineage-groups", "sec-site-aggregates", "sarow", "Lineage groups:"):
        assert needle not in html, needle
    assert html.count("<script>") == 1 and html.count("<style>") == 1


def test_no_lineage_group_table_section_or_nav_link(h_html):
    assert 'id="sec-lineage-groups"' not in h_html and "#sec-lineage-groups" not in h_html
    assert (
        'table class="sortable capped lg"' not in h_html and "Check your groups file" not in h_html
    )
    assert "BA.2.86* (excludes JN.1*)" not in h_html  # no group column heading anywhere


def test_provenance_line_about_the_lineage_groups(h_html):
    sec = section(h_html, "sec-provenance")
    assert (
        "Lineage groups: <code>groups.csv</code>, sha256 3cbeef2ff042, dataset 01_02_2026-03-04; "
        "table in <code>summary/lineage_summary.csv</code>"
    ) in sec


def test_provenance_line_without_a_configured_file_and_without_audit_facts(tmp_path):
    info = json.loads((H / "run_info.json").read_text())
    del info["params"]["lineage_groups"]
    ri = tmp_path / "ri.json"
    ri.write_text(json.dumps(info))
    out = tmp_path / "r.json"
    assert (
        brd.main(
            [
                "--run-info",
                str(ri),
                "--lineage-summary",
                str(H / "real_lineage_summary.csv"),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    html = render(json.loads(out.read_text()))
    assert_valid(html)
    assert (
        "Lineage groups: <code>bundled default</code>; table in <code>summary/lineage_summary.csv</code>"
        in html
    )


def test_lineage_group_flags_render(h_html):
    assert "lineage_groups_warnings" in h_html and "lineage_groups_ungrouped" in h_html
    assert (
        "add them to your lineage groups file; see summary/lineage_ungrouped_summary.tsv" in h_html
    )


def test_lineage_group_strings_are_escaped(tmp_path):
    evil, evil2 = "<script>alert(1)</script>", 'x"><img src=x onerror=alert(2)>'
    audit = tmp_path / "audit.tsv"
    audit.write_text(
        "group\tinclude\texclude\tin_total\tn_lineages_matched\twarnings\tdataset_version\tgroups_sha256\n"
        f"{evil}\tJN.1*\t\tyes\t2\t{evil2} | second {evil}\t{evil2}\t{evil}\n"
    )
    summary = tmp_path / "summary.csv"
    summary.write_text(
        "sample,status,other,ungrouped_lineages\n"
        f"{evil},ok,0.5,Evil<u>:0.5\n{evil2},ok,0.4,<i>Odd</i>:0.4\n"
    )
    info = json.loads((H / "run_info.json").read_text())
    info["params"]["lineage_groups"] = f"/x/{evil2}.csv"
    ri = tmp_path / "ri.json"
    ri.write_text(json.dumps(info))
    out = tmp_path / "r.json"
    args = ["--run-info", str(ri), "--lineage-summary", str(summary), "--lineage-groups-resolved", str(audit),
            "--out", str(out)]  # fmt: skip
    assert brd.main(args) == 0
    html = render(json.loads(out.read_text()))
    assert_valid(html)
    for needle in ("<script>alert", "<img src=x", "Evil<u>", "<i>Odd"):
        assert needle not in html, needle
    parser = parse(html)
    assert [t for t, _ in parser.tags].count("script") == 1
    assert all(not any(k.startswith("on") for k in attrs) for _, attrs in parser.tags)
    assert "lineage_groups_warnings" in html and "&lt;script&gt;" in html


def test_site_aggregate_strings_are_escaped(tmp_path):
    evil, evil2 = "<script>alert(1)</script>", 'x"><img src=x onerror=alert(2)>'
    meta = tmp_path / "meta.tsv"
    meta.write_text(
        "row_label\taggregate_name\tperiod\tperiod_start\tn_sites_included\tsites_included\tsites_missing\t"
        "missing_reasons\tweight_basis\tweight_fraction_included\tn_samples\n"
        f"Agg_2026-W01\t{evil}\t2026-W01\t2025-12-29\t1\t{evil2}\t{evil}\tnot_ok:{evil2}\t{evil2}\t0.5\t1\n"
    )
    freyja = tmp_path / "f.tsv"
    freyja.write_text(
        "\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
        "Agg_2026-W01\t[('Omicron<i>', 1.0)]\tEvil<u>\t1.0\t1\tNA\n"
    )
    out = tmp_path / "r.json"
    args = ["--run-info", str(H / "run_info.json"), "--site-aggregates", str(freyja),
            "--site-aggregates-meta", str(meta), "--out", str(out)]  # fmt: skip
    assert brd.main(args) == 0
    html = render(json.loads(out.read_text()))
    assert_valid(html)
    for needle in ("<script>alert", "<img src=x", "Evil<u>", "Omicron<i>"):
        assert needle not in html, needle
    parser = parse(html)
    assert [t for t, _ in parser.tags].count("script") == 2
    assert all(not any(k.startswith("on") for k in attrs) for _, attrs in parser.tags)
    assert 'id="sec-site-aggregates"' in html


def test_150_samples_with_site_aggregates_stay_small_and_fast():
    data = synthetic_data(150)
    data["lineage_groups"] = {"n_groups": 25, "groups_sha256": "a" * 64, "dataset_version": "x", "file": "g.csv",
                              "warnings": [], "n_samples_with_ungrouped": 3, "top_ungrouped": [["B.52", 0.01, 3]]}  # fmt: skip
    data["site_aggregates"] = [
        {"row_label": f"Agg{i % 4}_2026-W{i % 20 + 1:02d}", "aggregate_name": f"Agg{i % 4}", "period": f"2026-W{i % 20 + 1:02d}",
         "period_start": f"2026-01-{i % 20 + 1:02d}", "n_sites_included": 3, "sites_included": ["a", "b", "c"],
         "sites_missing": ["d"], "missing_reasons": {"d": "no_sample"}, "weight_basis": "viral_load",
         "weight_fraction_included": 0.8, "n_samples": 5,
         "groups": [["Omicron", 0.6], ["Delta", 0.4]], "lineages": [[f"X{j}", 0.1] for j in range(10)]}
        for i in range(60)
    ]  # fmt: skip
    start = time.perf_counter()
    html = render(data)
    elapsed = time.perf_counter() - start
    assert_valid(html)
    assert len(html.encode()) < 3_000_000, len(html)
    assert elapsed < 10 and html.count('class="sarow"') == 120


# =============================================================================================
# Anonymised reports: nothing identifying may remain anywhere in the HTML
# =============================================================================================
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


def build_anon(tmp_path, mode, extra_args=()):
    args = ["--run-info", str(A / "run_info.json")]
    for option, name in A_OPTIONS.items():
        args += [option, str(A / name)]
    out = tmp_path / f"{mode}.json"
    assert brd.main(args + ["--anonymize", mode, "--out", str(out), *extra_args]) == 0
    return json.loads(out.read_text())


@pytest.fixture(scope="module")
def anon_html(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("anon")
    return {mode: render(build_anon(tmp, mode)) for mode in ("none", "paths", "full")}


def test_full_report_html_contains_no_identifying_string(anon_html):
    html = anon_html["full"]
    assert_valid(html)
    for leak in A_IDS + A_NAMES + A_PATH_LEAKS:
        assert leak not in html, leak
    # pseudonyms are used everywhere instead: tables, flags, charts, tooltips, nav targets
    for pseudo in ("S001", "S002", "S005", "SITE1", "SITE2", "AGG1"):
        assert pseudo in html, pseudo
    assert 'id="anonymized-note"' in html and "Anonymised for sharing (full)" in html
    assert (
        "pseudonyms (S001, SITE1" in html
        and "lineage names, barcode versions and QC values are unchanged" in html
    )
    # structure that must survive
    assert "XBB.1.5* [Omicron (XBB.1.5.X)]" in html and "10_06_2026-05-25" in html
    assert "quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0" in html
    assert "2026-01-05" in html and "http://" not in html and "https://" not in html
    # tooltips and data attributes carry pseudonyms too
    assert 'data-sample="S002"' in html and 'title="S005"' in html
    assert "<title>wasteline report</title>" in html


def test_paths_report_html_has_no_path_user_or_email_but_keeps_ids(anon_html):
    html = anon_html["paths"]
    assert_valid(html)
    for leak in A_PATH_LEAKS:
        assert leak not in html, leak
    for sid in A_IDS:
        assert sid in html
    assert (
        "Anonymised for sharing (paths)" in html
        and "pseudonyms" not in html.split('id="anonymized-note"')[1].split("</p>")[0]
    )
    assert (
        "<code>.../main.nf" not in html and "nextflow run .../main.nf" in html
    )  # sanitised command line (provenance block)
    assert ".../samplesheet.csv" in html and "&lt;user&gt;" in html  # run_name and the reason text


def test_default_report_is_not_anonymised_and_has_no_note(anon_html):
    html = anon_html["none"]
    assert 'id="anonymized-note"' not in html
    assert "/Users/jdoe" in html and "WW_Central_001" in html


def test_hostile_ids_in_full_and_paths_html_are_inert(tmp_path):
    evil = "<script>alert(1)</script>"
    qc = tmp_path / "qc.tsv"
    qc.write_text(
        "sample\tsample_type\tcollection_date\tsite\tqc_status\tqc_reasons\n"
        f"{evil}\tsample\t2026-01-01\tSite<b>\tFAIL\tlow_coverage\nWW_1\tsample\t2026-01-08\tSite<b>\tPASS\t\n"
    )
    for mode in ("full", "paths"):
        out = tmp_path / f"{mode}.json"
        args = [
            "--run-info",
            str(A / "run_info.json"),
            "--qc",
            str(qc),
            "--anonymize",
            mode,
            "--out",
            str(out),
        ]
        assert brd.main(args) == 0
        html = render(json.loads(out.read_text()))
        assert_valid(html)
        parser = parse(html)
        assert [t for t, _ in parser.tags].count("script") == 1
        assert "<script>alert" not in html
        assert ("script&gt;alert(1)" in html) is (
            mode == "paths"
        )  # escaped when kept, gone when pseudonymised
        assert ("Site&lt;b&gt;" in html) is (mode == "paths")


def test_anonymised_150_samples_stay_small(tmp_path):
    data = synthetic_data(150)
    info = json.loads((A / "run_info.json").read_text())
    data["provenance"] = info
    anon = brd.Anonymizer("full", data)
    out = anon.apply(data)
    html = render(out)
    assert_valid(html)
    assert len(html.encode()) < 3_000_000
    assert not any(sid in html for sid in ("NTC_000", "WW_001_site1")) and "S150" in html


def test_full_anonymisation_covers_every_analysis_table(tmp_path):
    data = build_g(
        tmp_path,
        extra=["--consensus", str(G / "consensus_stats.tsv"), "--dashboard-name", "x.html"],
    )
    ids = [s["sample"] for s in data["samples"]]
    sites = {s["site"] for s in data["samples"] if s["site"]}
    assert "SiteA_W03_1" in ids and sites == {"SiteA", "SiteB", "SiteC"}
    anon = brd.Anonymizer("full", data)
    html = render(anon.apply(data))
    assert_valid(html)
    for sid in ids:
        assert sid not in html, sid
    for site in sites:
        assert site not in html, site
    for section_id in ("sec-dropout-impact", "sec-unexplained", "sec-timeseries"):
        assert f'id="{section_id}"' in html
    assert "SITE1" in html and "S0" in html
    # lineage names, p-values and growth rates are untouched
    assert "XBB.1.5* [Omicron (XBB.1.5.X)]" in html and "2.39e-42" in html and "+60.3" in html
    mapped = {row[1]: row[2] for row in anon.mapping_rows()}
    assert mapped["SiteA_W03_1"].startswith("S") and mapped["SiteB"] == "SITE2"


# =============================================================================================
# Defining-sites dropout section and the unexplained `in_dataset` marker
# =============================================================================================
J = FIX / "j_defining_sites"


def build_j(tmp_path):
    args = ["--run-info", str(J / "run_info.json")]
    for option, name in G_OPTIONS.items():
        base = J if option in ("--dropout-impact", "--dropout-summary", "--unexplained") else G
        args += [option, str(base / name)]
    out = tmp_path / "j.json"
    assert brd.main(args + ["--out", str(out)]) == 0
    return json.loads(out.read_text())


@pytest.fixture(scope="module")
def j_html(tmp_path_factory):
    return render(build_j(tmp_path_factory.mktemp("j")))


def test_dropout_section_shows_defining_sites_first(j_html):
    assert_valid(j_html)
    sec = section(j_html, "sec-dropout-impact")
    assert (
        '<abbr title="Defining mutations are the mutations a lineage has that its parent lineage lacks'
        in sec
    )
    assert (
        "defining mutations</abbr>" in sec and "(abundance &ge; 1.0%)" in sec and "&ge; 5.0%" in sec
    )
    # per-sample summary: defining sites lost (primary), pct defining, informative as secondary context
    head = sec[sec.index("Per-sample summary") :]
    head = head[: head.index("</thead>")]
    assert (
        "Defining sites lost" in head
        and "% defining sites lost" in head
        and "Informative sites lost (%)" in head
    )
    assert head.index("% defining sites lost") < head.index("Informative sites lost (%)")
    assert "<b>2.50</b>" in sec and "3 / 120" in sec and "850 / 8,006 (10.6)" in sec
    assert "no lineage estimate" in sec  # LOWCOV_1
    # amplicon table: defining sites hidden is the primary count
    table = sec[sec.index('id="dx-table"') :]
    th = table[: table.index("</thead>")]
    assert (
        th.index("Defining sites hidden")
        < th.index("Lineages losing defining sites")
        < th.index("Differential sites")
        < th.index("Informative sites")
    )
    assert (
        "<b>3</b>" in table and "XBB.1.5 (2), BA.2 (1)" in table and "hides defining sites" in table
    )
    # a lineage under the affect threshold is listed but not highlighted
    quiet = table[table.index("PQ.2.1.17 (1)") - 700 : table.index("PQ.2.1.17 (1)")]
    assert 'class="hl' not in quiet.split("<tr")[-1]
    assert (
        "hides defining sites" in j_html and "Amplicon dropout hides defining mutations" in j_html
    )


def test_old_dropout_tables_keep_the_old_layout(tmp_path):
    html = render(build_g(tmp_path))
    sec = section(html, "sec-dropout-impact")
    assert "Defining sites" not in sec and "lineage-informative</b> sites it covers" in sec
    assert "Affected lineages (informative sites)" in sec and "affects detected lineages" in sec
    assert "XBB.1.5 (9), BA.2 (3), BA.5 (3)" in sec


def test_unexplained_marker_floor_and_note(j_html):
    sec = section(j_html, "sec-unexplained")
    assert (
        sec.count('class="chip" title="No lineage in the Freyja dataset models this mutation') == 4
    )
    assert ">not in dataset<" in sec
    assert "&lt;=1e-300" in sec and "1.41e-17" in sec and "2.28e-14" in sec
    assert (
        "P-values are floored at 1e-300" in sec
        and "dataset blind spot" in sec
        and "(2 mutations here)" in sec
    )
    # recurrent rows that the dataset does model carry no marker
    n501 = sec[sec.index("N501Y") : sec.index("N501Y") + 900]
    assert "not in dataset" not in n501.split("</tr>")[0]
    assert "(1 of them at sites no lineage in the dataset models)" in j_html


def test_unexplained_without_in_dataset_column_has_no_marker_or_floor_note(g_html):
    sec = section(g_html, "sec-unexplained")
    assert "not in dataset" not in sec and "floored" not in sec
    assert "2.39e-42" in sec


def test_pval_filter():
    assert rr.pval(1e-300) == "<=1e-300" and rr.pval(2.5e-310) == "<=1e-300"
    assert rr.pval(2.39e-42) == "2.39e-42" and rr.pval(None) == "–"


def test_new_dropout_and_unexplained_strings_are_escaped(tmp_path):
    evil = "<script>alert(1)</script>"
    imp = tmp_path / "imp.tsv"
    imp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tn_informative_sites\tn_sites_in_detected_lineages\taffects_detected_lineage\t"
        "affected_lineages\tn_defining_sites\tn_differential_sites\n"
        f"{evil}\t{evil}\t1\t100\t0\t5\t2\ttrue\t{evil}:2\t2\t1\n"
    )
    ue = tmp_path / "ue.tsv"
    ue.write_text(
        UE_HEADER
        + "\tin_dataset\n"
        + f"{evil}\t5\tA\tT\t0.3\t0\t0.3\t100\t30\t{evil}\t{evil}\t{evil}\t2\ttrue\t1e-300\tfalse\n"
    )
    data_args = ["--run-info", str(J / "run_info.json"), "--qc", str(G / "qc_summary.tsv"), "--dropout-impact", str(imp),
                 "--unexplained", str(ue), "--out", str(tmp_path / "o.json")]  # fmt: skip
    assert brd.main(data_args) == 0
    html = render(json.loads((tmp_path / "o.json").read_text()))
    assert_valid(html)
    assert "<script>alert" not in html and [t for t, _ in parse(html).tags].count("script") == 2
