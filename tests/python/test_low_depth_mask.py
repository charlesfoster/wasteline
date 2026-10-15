import low_depth_mask as ldm
import pytest


def write_depths(path, depths, contig="ref"):
    path.write_text("".join(f"{contig}\t{i}\tA\t{d}\n" for i, d in enumerate(depths, start=1)))
    return str(path)


def run(tmp_path, depths, *extra, contig="ref"):
    dep = write_depths(tmp_path / "d.depths", depths, contig)
    out = tmp_path / "mask.bed"
    rc = ldm.main(["--depths", dep, "-o", str(out), *extra])
    return rc, out.read_text() if out.exists() else None


def test_merged_intervals_zero_based_half_open(tmp_path):
    # positions 1-based: 1,2 low | 3,4 ok | 5 low | 6 ok | 7,8 low
    rc, text = run(tmp_path, [0, 1, 20, 20, 14, 15, 3, 0], "--min-depth", "15")
    assert rc == 0
    assert text == "ref\t0\t2\nref\t4\t5\nref\t6\t8\n"


def test_depth_equal_to_min_is_not_masked(tmp_path):
    _, text = run(tmp_path, [15, 15, 15], "--min-depth", "15")
    assert text == ""


def test_all_zero_depth_masks_whole_contig(tmp_path):
    _, text = run(tmp_path, [0] * 100, "--min-depth", "15")
    assert text == "ref\t0\t100\n"


def test_start_and_end_edges(tmp_path):
    _, text = run(tmp_path, [0, 0, 30, 30, 0], "--min-depth", "15")
    assert text == "ref\t0\t2\nref\t4\t5\n"


def test_zero_depth_stretch_in_middle(tmp_path):
    _, text = run(tmp_path, [30] * 3 + [0] * 4 + [30] * 3, "--min-depth", "15")
    assert text == "ref\t3\t7\n"


def test_length_masks_unlisted_tail(tmp_path):
    _, text = run(tmp_path, [30] * 10, "--min-depth", "15", "--length", "15")
    assert text == "ref\t10\t15\n"


def test_length_merges_with_trailing_low_depth(tmp_path):
    _, text = run(tmp_path, [30] * 8 + [0, 0], "--min-depth", "15", "--length", "12")
    assert text == "ref\t8\t12\n"


def test_unlisted_gap_is_masked(tmp_path):
    dep = tmp_path / "d.depths"
    dep.write_text("ref\t1\tA\t30\nref\t2\tA\t30\nref\t5\tA\t30\n")
    out = tmp_path / "m.bed"
    assert ldm.main(["--depths", str(dep), "-o", str(out)]) == 0
    assert out.read_text() == "ref\t2\t4\n"


def test_contig_name_comes_from_depths(tmp_path):
    _, text = run(tmp_path, [0, 0], contig="NC_045512.2")
    assert text == "NC_045512.2\t0\t2\n"


def test_keep_bed_subtracts_variant_positions(tmp_path):
    keep = tmp_path / "keep.bed"
    keep.write_text("ref\t3\t4\nref\t0\t1\nref\t8\t12\nother\t0\t100\n")
    _, text = run(tmp_path, [0] * 10, "--min-depth", "15", "--keep-bed", str(keep))
    assert text == "ref\t1\t3\nref\t4\t8\n"


def test_keep_bed_covering_everything_gives_empty_mask(tmp_path):
    keep = tmp_path / "keep.bed"
    keep.write_text("ref\t0\t50\n")
    _, text = run(tmp_path, [0] * 10, "--keep-bed", str(keep))
    assert text == ""


def test_keep_bed_ignores_headers_and_extra_columns(tmp_path):
    keep = tmp_path / "keep.bed"
    keep.write_text("track name=x\n#comment\nref\t2\t3\tvar\t0.9\t+\n")
    _, text = run(tmp_path, [0] * 5, "--keep-bed", str(keep))
    assert text == "ref\t0\t2\nref\t3\t5\n"


@pytest.mark.parametrize(
    "intervals,keep,expected",
    [
        ([(0, 10)], [(2, 4), (6, 8)], [(0, 2), (4, 6), (8, 10)]),
        ([(0, 10)], [], [(0, 10)]),
        ([(5, 10)], [(0, 7)], [(7, 10)]),
        ([(0, 5), (10, 15)], [(4, 11)], [(0, 4), (11, 15)]),
        ([(0, 5)], [(5, 9)], [(0, 5)]),
    ],
)
def test_subtract(intervals, keep, expected):
    assert ldm.subtract(intervals, keep) == expected


def test_merge_adjacent_and_overlapping():
    assert ldm.merge([(5, 8), (0, 2), (2, 3), (7, 9)]) == [(0, 3), (5, 9)]


def test_empty_depths_is_an_error(tmp_path):
    dep = tmp_path / "d.depths"
    dep.write_text("")
    assert ldm.main(["--depths", str(dep), "-o", str(tmp_path / "m.bed")]) == 2


def test_missing_depths_is_an_error(tmp_path):
    assert ldm.main(["--depths", str(tmp_path / "nope"), "-o", str(tmp_path / "m.bed")]) == 2


def test_length_with_multiple_contigs_is_an_error(tmp_path):
    dep = tmp_path / "d.depths"
    dep.write_text("a\t1\tA\t1\nb\t1\tA\t1\n")
    assert ldm.main(["--depths", str(dep), "--length", "5", "-o", str(tmp_path / "m.bed")]) == 2


def test_full_size_genome(tmp_path):
    depths = [30] * 29903
    for i in range(50):
        depths[i] = 0
    for i in range(29903 - 20, 29903):
        depths[i] = 2
    _, text = run(tmp_path, depths, "--min-depth", "15", "--length", "29903", contig="NC_045512.2")
    assert text == "NC_045512.2\t0\t50\nNC_045512.2\t29883\t29903\n"
