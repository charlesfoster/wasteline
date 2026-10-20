import json

import depth_stats
import pytest


def write_depths(path, depths, contig="ref"):
    path.write_text("".join(f"{contig}\t{i}\tA\t{d}\n" for i, d in enumerate(depths, start=1)))
    return str(path)


def test_compute_stats(tmp_path):
    depths = write_depths(tmp_path / "d.depths", [0, 5, 10, 20, 100])
    stats = depth_stats.compute_stats(depths, covcut=10)
    assert stats["positions"] == 5
    assert stats["positions_ge_covcut"] == 3
    assert stats["pct_genome_ge_covcut"] == 60.0
    assert stats["mean_depth"] == pytest.approx(27.0)


def test_covcut_is_inclusive(tmp_path):
    depths = write_depths(tmp_path / "d.depths", [9, 10, 11])
    assert depth_stats.compute_stats(depths, covcut=10)["positions_ge_covcut"] == 2


def test_gate_pass_writes_json(tmp_path):
    depths = write_depths(tmp_path / "d.depths", [10] * 8 + [0] * 2)
    out = tmp_path / "stats.json"
    rc = depth_stats.main(
        [
            depths,
            "--covcut",
            "10",
            "--min-genome-coverage",
            "60",
            "--sample",
            "S_1",
            "--output",
            str(out),
        ]
    )
    assert rc == depth_stats.EXIT_OK
    record = json.loads(out.read_text())
    assert record["sample"] == "S_1"
    assert record["passes_gate"] is True
    assert record["pct_genome_ge_covcut"] == 80.0


def test_gate_fail_returns_10_but_still_writes_json(tmp_path):
    depths = write_depths(tmp_path / "d.depths", [10] * 3 + [0] * 7)
    out = tmp_path / "stats.json"
    rc = depth_stats.main([depths, "--min-genome-coverage", "60", "--output", str(out)])
    assert rc == depth_stats.EXIT_BELOW_GATE
    assert json.loads(out.read_text())["passes_gate"] is False


def test_low_confidence_band(tmp_path):
    out = tmp_path / "stats.json"
    in_band = write_depths(
        tmp_path / "a.depths", [10] * 4 + [0] * 6
    )  # 40 %: floor 25 <= 40 < gate 50
    rc = depth_stats.main(
        [
            in_band,
            "--min-genome-coverage",
            "50",
            "--floor-genome-coverage",
            "25",
            "--output",
            str(out),
        ]
    )
    assert rc == depth_stats.EXIT_LOW_CONFIDENCE
    record = json.loads(out.read_text())
    assert record["low_confidence"] is True and record["passes_gate"] is False
    below = write_depths(tmp_path / "b.depths", [10] * 2 + [0] * 8)  # 20 %: below the floor
    assert (
        depth_stats.main([below, "--min-genome-coverage", "50", "--floor-genome-coverage", "25"])
        == depth_stats.EXIT_BELOW_GATE
    )
    for floor in ("0", "50", "60"):  # band off: floor 0 or not below the gate
        assert (
            depth_stats.main(
                [in_band, "--min-genome-coverage", "50", "--floor-genome-coverage", floor]
            )
            == depth_stats.EXIT_BELOW_GATE
        )


def test_gate_boundary_is_inclusive(tmp_path):
    depths = write_depths(tmp_path / "d.depths", [10] * 6 + [0] * 4)
    assert depth_stats.main([depths, "--min-genome-coverage", "60"]) == depth_stats.EXIT_OK


def test_no_gate_always_passes(tmp_path):
    depths = write_depths(tmp_path / "d.depths", [0] * 10)
    assert depth_stats.main([depths]) == depth_stats.EXIT_OK


def test_empty_file_is_an_error(tmp_path):
    path = tmp_path / "empty.depths"
    path.write_text("")
    assert depth_stats.main([str(path), "--min-genome-coverage", "0"]) == depth_stats.EXIT_ERROR


def test_malformed_rows_are_errors(tmp_path):
    short = tmp_path / "short.depths"
    short.write_text("ref\t1\tA\n")
    assert depth_stats.main([str(short)]) == depth_stats.EXIT_ERROR
    bad = tmp_path / "bad.depths"
    bad.write_text("ref\t1\tA\tNA\n")
    assert depth_stats.main([str(bad)]) == depth_stats.EXIT_ERROR


def test_missing_file_is_an_error(tmp_path):
    assert depth_stats.main([str(tmp_path / "nope.depths")]) == depth_stats.EXIT_ERROR


def test_blank_lines_ignored(tmp_path):
    path = tmp_path / "d.depths"
    path.write_text("ref\t1\tA\t10\n\nref\t2\tA\t10\n")
    assert depth_stats.compute_stats(str(path), 10)["positions"] == 2
