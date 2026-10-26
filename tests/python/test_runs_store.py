"""bin/runs_store.py: a cross-run store and its alerts, on four tiny hand-made weekly results directories.

Site A, one sample per week (2026-09-07 .. 09-28, ISO weeks 37-40), plus a low_confidence sample in the last
week. The groups: G.new appears only in the last week (emerging), G.grow rises every week (sustained growth),
G.low is only in the low_confidence sample. C6501A is unexplained in the last three weeks (recurring), and
S:N501Y is on the watchlist.
"""

import csv
import json
import shutil
import sqlite3

import pytest
import runs_store as rs

WEEKS = ["2026-09-07", "2026-09-14", "2026-09-21", "2026-09-28"]
GROW = [0.05, 0.10, 0.20, 0.30]


def make_results(root, name, start, barcode_iso, samples):
    """samples: [(id, date, status, {group: abundance}, has_unexplained, n501y_freq)] -> the results dir."""
    out = root / name
    (out / "summary").mkdir(parents=True)
    (out / "pipeline_info").mkdir()
    (out / "pipeline_info" / "run_info.json").write_text(
        json.dumps(
            {
                "run_name": name,
                "start": start,
                "mode": "full",
                "success": True,
                "freyja_data": {"barcode_version": "x", "barcode_version_iso": barcode_iso},
            }
        )
    )

    def table(file, header, rows, sep="\t"):
        lines = [sep.join(header)] + [sep.join(str(c) for c in r) for r in rows]
        (out / "summary" / file).write_text("\n".join(lines) + "\n")

    table(
        "demix_status.tsv",
        ["sample", "status", "reason", "genome_coverage_pct", "barcode_version"],
        [(s[0], s[2], "", 90, "x") for s in samples],
    )
    table(
        "samplesheet.csv",
        ["sample", "sample_type", "collection_date", "site"],
        [(s[0], "sample", s[1], "A") for s in samples],
        ",",
    )
    table(
        "lineage_summary.csv",
        [
            "sample",
            "status",
            "genome_coverage_pct",
            "G.new",
            "G.grow",
            "G.low",
            "other",
            "proportion_summary",
        ],
        [
            (
                s[0],
                s[2],
                90,
                s[3].get("G.new", ""),
                s[3].get("G.grow", ""),
                s[3].get("G.low", ""),
                "",
                1,
            )
            for s in samples
        ],
        ",",
    )
    table(
        "freyja_aggregated.tsv",
        ["", "lineages", "abundances"],
        [(s[0], "X.1", "1.0") for s in samples],
    )
    table(
        "unexplained_mutations.tsv",
        [
            "sample",
            "position",
            "ref",
            "alt",
            "observed_freq",
            "gene",
            "aa_change",
            "n_samples",
            "in_dataset",
        ],
        [(s[0], 6501, "C", "A", 0.4, "ORF1ab", "P2079Q", 3, "false") for s in samples if s[4]],
    )
    table(
        "mutations_of_interest.tsv",
        ["sample", "name", "mutation", "ref_depth", "alt_depth", "depth", "alt_freq"],
        [(s[0], "S N501Y", "S:N501Y", 100, 10, 110, s[5]) for s in samples],
    )
    return out


def weekly(root, new_in_last=0.20):
    runs = []
    for i, date in enumerate(WEEKS):
        groups = {"G.grow": GROW[i], "G.new": new_in_last if i == 3 else 0.0}
        samples = [(f"S{i + 1}", date, "ok", groups, i >= 1, 0.30 if i == 3 else 0.02)]
        if i == 3:
            samples.append(("S5", date, "low_confidence", {"G.low": 0.40}, False, 0.0))
        runs.append(
            make_results(root, f"run{i + 1}", date + "T09:00:00+10:00", "2026-10-06T05-25", samples)
        )
    return runs


def tsv(path):
    with open(str(path), newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def alerts(store, tmp_path, *extra):
    watch = tmp_path / "watch.csv"
    watch.write_text("name,mutation\nS N501Y,S:N501Y\n")
    out = tmp_path / "alerts.tsv"
    argv = ["alerts", str(store), "--watchlist", str(watch), "--out", str(out)] + list(extra)
    assert rs.main(argv) == 0
    return tsv(out)


def table_counts(store):
    conn = sqlite3.connect(str(store))
    names = [r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type = 'table'")]
    return {n: conn.execute("SELECT COUNT(*) FROM " + n).fetchone()[0] for n in names}


@pytest.fixture
def store(tmp_path):
    dirs = weekly(tmp_path / "results")
    path = tmp_path / "store.sqlite"
    assert rs.main(["add", str(path)] + [str(d) for d in dirs]) == 0
    return path


def test_add_then_alerts_round_trip(store, tmp_path):
    found = alerts(store, tmp_path, "--markdown", str(tmp_path / "alerts.md"))
    assert [
        (
            a["rule"],
            a["level"],
            a["site"],
            a["group_or_mutation"],
            a["first_date"],
            a["last_date"],
            a["samples"],
        )
        for a in found
    ] == [
        ("emerging", "group", "A", "G.new", WEEKS[3], WEEKS[3], "S4"),
        ("sustained_growth", "group", "A", "G.grow", WEEKS[0], WEEKS[3], "S1,S2,S3,S4"),
        (
            "recurring_unexplained",
            "mutation",
            "A",
            "C6501A (ORF1ab:P2079Q)",
            WEEKS[1],
            WEEKS[3],
            "S2,S3,S4",
        ),
        ("watchlist", "mutation", "A", "S N501Y (S:N501Y)", WEEKS[3], WEEKS[3], "S4"),
    ]
    assert "in_dataset=no" in found[2]["detail"]
    assert "not epidemiological conclusions" in (tmp_path / "alerts.md").read_text()


def test_re_adding_is_idempotent_and_queries_use_the_latest_run_per_sample(store, tmp_path):
    before = table_counts(store)
    again = [str(tmp_path / "results" / f"run{i}") for i in (1, 2, 3, 4)]
    assert rs.main(["add", str(store)] + again) == 0
    assert table_counts(store) == before and before["runs"] == 4

    # a re-analysis of the last week with newer barcodes no longer sees G.new: that sample's latest result wins
    redo = make_results(
        tmp_path / "results",
        "run4_reanalysis",
        "2026-10-20T09:00:00+10:00",
        "2026-10-19T00-00",
        [("S4", WEEKS[3], "ok", {"G.grow": 0.30, "G.new": 0.0}, True, 0.30)],
    )
    assert rs.main(["add", str(store), str(redo)]) == 0
    assert table_counts(store)["runs"] == 5
    assert {a["rule"] for a in alerts(store, tmp_path)} == {
        "sustained_growth",
        "recurring_unexplained",
        "watchlist",
    }


def test_low_confidence_samples_are_left_out_of_emerging_unless_asked(store, tmp_path):
    assert "G.low" not in {a["group_or_mutation"] for a in alerts(store, tmp_path)}
    with_low = alerts(store, tmp_path, "--include-low-confidence")
    assert ("emerging", "G.low") in {(a["rule"], a["group_or_mutation"]) for a in with_low}


def test_panel_results_are_stored_per_target_and_alerts_never_mix_targets(tmp_path):
    """Panel results (summary/<t>/, pipeline_info/<t>/run_info.json) are one run per target, idempotent per
    (dir, target). Target TA is the weekly series above, TB has the same sample IDs and dates but nothing
    to alert on: only TA alerts, the TSV gets a target column, and TB does not overwrite TA's samples."""
    stage, panel = tmp_path / "stage", tmp_path / "panel"
    for i, date in enumerate(WEEKS):
        for target, quiet in (("TA", False), ("TB", True)):
            groups = {} if quiet else {"G.grow": GROW[i], "G.new": 0.20 if i == 3 else 0.0}
            src = make_results(
                stage / target,
                f"week{i + 1}",
                date + "T09:00:00+10:00",
                "2026-10-06T05-25",
                [(f"S{i + 1}", date, "ok", groups, i >= 1 and not quiet, 0.02 if quiet else 0.30)],
            )
            (panel / f"week{i + 1}" / "pipeline_info" / target).mkdir(parents=True)
            shutil.copytree(src / "summary", panel / f"week{i + 1}" / "summary" / target)
            shutil.copy(
                src / "pipeline_info" / "run_info.json",
                panel / f"week{i + 1}" / "pipeline_info" / target / "run_info.json",
            )
    store = tmp_path / "store.sqlite"
    dirs = [str(panel / f"week{i}") for i in (1, 2, 3, 4)]
    assert rs.main(["add", str(store), *dirs]) == 0
    assert rs.main(["add", str(store), str(panel)]) == 0  # a folder of panel results, again
    counts = table_counts(store)
    assert counts["runs"] == 8 and counts["samples"] == 8
    found = alerts(store, tmp_path)
    assert list(found[0])[:3] == ["rule", "target", "level"]
    assert {(a["rule"], a["target"]) for a in found} == {
        ("emerging", "TA"),
        ("sustained_growth", "TA"),
        ("recurring_unexplained", "TA"),
        ("watchlist", "TA"),
    }
