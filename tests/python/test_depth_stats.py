import json

import depth_stats
import pytest


def write_depths(path, depths, contig="ref"):
    path.write_text("".join(f"{contig}\t{i}\tA\t{d}\n" for i, d in enumerate(depths, start=1)))
    return str(path)


def test_stats_and_gate_exit_codes(tmp_path):
    stats = depth_stats.compute_stats(write_depths(tmp_path / "s.depths", [0, 5, 10, 20, 100]), 10)
    assert stats["positions_ge_covcut"] == 3  # covcut is inclusive
    assert stats["pct_genome_ge_covcut"] == 60.0
    assert stats["mean_depth"] == pytest.approx(27.0)

    out = tmp_path / "stats.json"
    ok = write_depths(tmp_path / "ok.depths", [10] * 8 + [0] * 2)
    args = [
        "--covcut",
        "10",
        "--min-genome-coverage",
        "60",
        "--sample",
        "S_1",
        "--output",
        str(out),
    ]
    assert depth_stats.main([ok, *args]) == depth_stats.EXIT_OK
    record = json.loads(out.read_text())
    assert (record["sample"], record["passes_gate"]) == ("S_1", True)

    low = write_depths(tmp_path / "low.depths", [10] * 3 + [0] * 7)
    assert depth_stats.main([low, *args]) == depth_stats.EXIT_BELOW_GATE
    assert json.loads(out.read_text())["passes_gate"] is False  # JSON still written


def test_low_confidence_band(tmp_path):
    out = tmp_path / "stats.json"
    in_band = write_depths(
        tmp_path / "a.depths", [10] * 4 + [0] * 6
    )  # 40 %: floor 25 <= 40 < gate 50
    args = ["--min-genome-coverage", "50", "--floor-genome-coverage", "25"]
    assert (
        depth_stats.main([in_band, *args, "--output", str(out)]) == depth_stats.EXIT_LOW_CONFIDENCE
    )
    record = json.loads(out.read_text())
    assert record["low_confidence"] is True and record["passes_gate"] is False
    below = write_depths(tmp_path / "b.depths", [10] * 2 + [0] * 8)  # 20 %: below the floor
    assert depth_stats.main([below, *args]) == depth_stats.EXIT_BELOW_GATE
