import csv

import panel_select
import pytest


def write(path, text):
    path.write_text(text)
    return str(path)


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_panel(tmp_path, with_second_target=True):
    """RSVa (contigs HA + NA, marker HA) and RSVb (one contig that is also called HA)."""
    a = write(tmp_path / "a.fa", ">HA desc\nACGT\nAC\n>NA\nGGGG\n")
    b = write(tmp_path / "b.fa", ">HA\nTTTT\n")
    rows = f"target,fasta,marker_contig,family\nRSVa,{a},HA,RSV\n"
    if with_second_target:
        rows += f"RSVb,{b},,RSV\n"
    return write(tmp_path / "panel_resolved.csv", rows)


def test_prepare_renames_contigs_per_target_but_keeps_a_single_target_as_is(tmp_path):
    fasta, cmap = str(tmp_path / "panel.fasta"), str(tmp_path / "contig_map.tsv")
    assert (
        panel_select.main(
            ["prepare", "--panel", make_panel(tmp_path), "--out-fasta", fasta, "--out-map", cmap]
        )
        == 0
    )
    assert (
        tmp_path / "panel.fasta"
    ).read_text() == ">RSVa|HA\nACGT\nAC\n>RSVa|NA\nGGGG\n>RSVb|HA\nTTTT\n"
    rows = read_tsv(cmap)
    assert [(r["panel_contig"], r["is_marker"], r["family"]) for r in rows] == [
        ("RSVa|HA", "true", "RSV"),
        ("RSVa|NA", "false", "RSV"),
        ("RSVb|HA", "true", "RSV"),  # no marker_contig: every contig is a marker
    ]
    single = [
        "prepare",
        "--panel",
        make_panel(tmp_path, False),
        "--out-fasta",
        fasta,
        "--out-map",
        cmap,
    ]
    assert panel_select.main(single) == 0
    assert [r["panel_contig"] for r in read_tsv(cmap)] == ["HA", "NA"]


# fragments per target (+ aligned bases = 100 x fragments, breadth 40 % when any), panel_fragments
CASES = {
    "pure": ({"RSVa": 500, "RSVb": 0}, {"RSVa": "detected", "RSVb": "not_detected"}, "ok", "false"),
    "mixed": ({"RSVa": 700, "RSVb": 300}, {"RSVa": "detected", "RSVb": "detected"}, "ok", "true"),
    "trace_floor": (
        {"RSVa": 540, "RSVb": 40},
        {"RSVa": "detected", "RSVb": "trace"},
        "ok",
        "false",
    ),
    "trace_share": (
        {"RSVa": 5000, "RSVb": 120},
        {"RSVa": "detected", "RSVb": "trace"},
        "ok",
        "false",
    ),
    "no_reads": (
        {"RSVa": 0, "RSVb": 0},
        {"RSVa": "not_detected", "RSVb": "not_detected"},
        "qc_fail",
        "false",
    ),
}


def test_call_statuses_and_one_row_per_target_even_without_reads(tmp_path):
    for case, spec in CASES.items():
        (tmp_path / case).mkdir()
        check_call(tmp_path / case, case, *spec)


def check_call(tmp_path, case, fragments, expected, sample_status, mixed):
    fasta, cmap = str(tmp_path / "panel.fasta"), str(tmp_path / "contig_map.tsv")
    panel_select.main(
        ["prepare", "--panel", make_panel(tmp_path), "--out-fasta", fasta, "--out-map", cmap]
    )
    header = "\t".join(panel_select.STATS_COLUMNS) + "\n"
    stats = header + "".join(
        f"S\t{t}\t{n}\t0\t{n * 100}\t{60 if n else 0}\t{40 if n else 0}\t{n / 10}\t{sum(fragments.values())}\n"
        for t, n in fragments.items()
    )
    out = str(tmp_path / "S.targets.tsv")
    args = [
        "call",
        "--sample",
        "S",
        "--sample-type",
        "sample",
        "--stats",
        write(tmp_path / "stats.tsv", stats),
    ]
    assert panel_select.main([*args, "--contig-map", cmap, "--output", out]) == 0
    rows = read_tsv(out)
    assert list(rows[0]) == panel_select.TARGET_COLUMNS
    assert {r["target"]: r["status"] for r in rows} == expected
    assert {r["sample_status"] for r in rows} == {sample_status}
    assert {r["mixed"] for r in rows if r["status"] == "detected"} <= {mixed}
    assert {r["analysed"] for r in rows if r["status"] != "detected"} <= {"false"}
    by_target = {r["target"]: r for r in rows}
    if case == "mixed":
        assert by_target["RSVa"]["secondary_targets"] == "RSVb"
        assert float(by_target["RSVa"]["fraction_of_family"]) == pytest.approx(0.7)
    if case == "trace_floor":
        assert by_target["RSVb"]["reason"].startswith("low_fragments(40<100)")
    if case == "trace_share":
        assert by_target["RSVb"]["reason"] == "low_family_share(2.34<5)"
    if case == "no_reads":
        assert {r["reason"] for r in rows} == {"no_panel_reads"}
        # a stats file that is empty (assign died) must still give the rows
        empty = write(tmp_path / "empty.tsv", "")
        assert (
            panel_select.main(
                ["call", "--sample", "S", "--stats", empty, "--contig-map", cmap, "--output", out]
            )
            == 0
        )
        assert len(read_tsv(out)) == 2


def test_prepare_fails_clearly_on_a_duplicated_contig_missing_marker_and_missing_fasta(
    tmp_path, capsys
):
    dup = write(tmp_path / "dup.fa", ">HA\nAC\n>HA\nGG\n")
    ok = write(tmp_path / "ok.fa", ">HA\nAC\n")
    out = ["--out-fasta", str(tmp_path / "p.fa"), "--out-map", str(tmp_path / "m.tsv")]
    for rows, message in [
        (f"{dup},HA", "appears twice"),
        (f"{ok},NA", "marker_contig 'NA'"),
        (f"{tmp_path / 'nope.fa'},HA", "FASTA not found"),
    ]:
        fasta, marker = rows.split(",")
        panel = write(tmp_path / "panel.csv", f"target,fasta,marker_contig\nT,{fasta},{marker}\n")
        assert panel_select.main(["prepare", "--panel", panel, *out]) == 2
        assert message in capsys.readouterr().err
