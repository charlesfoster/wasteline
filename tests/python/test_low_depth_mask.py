import low_depth_mask as ldm


def write_depths(path, depths, contig="ref"):
    path.write_text("".join(f"{contig}\t{i}\tA\t{d}\n" for i, d in enumerate(depths, start=1)))
    return str(path)


def run(tmp_path, depths, *extra, contig="ref"):
    dep = write_depths(tmp_path / "d.depths", depths, contig)
    out = tmp_path / "mask.bed"
    rc = ldm.main(["--depths", dep, "-o", str(out), *extra])
    return rc, out.read_text() if out.exists() else None


def test_merged_intervals_zero_based_half_open(tmp_path):
    # positions 1-based: 1,2 low | 3,4 ok | 5 low | 6 ok (depth == min is not masked) | 7,8 low
    rc, text = run(tmp_path, [0, 1, 20, 20, 14, 15, 3, 0], "--min-depth", "15")
    assert rc == 0
    assert text == "ref\t0\t2\nref\t4\t5\nref\t6\t8\n"
    # --length masks the unlisted tail and merges with trailing low depth
    _, text = run(tmp_path, [30] * 8 + [0, 0], "--min-depth", "15", "--length", "12")
    assert text == "ref\t8\t12\n"


def test_keep_bed_subtracts_variant_positions(tmp_path):
    keep = tmp_path / "keep.bed"
    keep.write_text("ref\t3\t4\nref\t0\t1\nref\t8\t12\nother\t0\t100\n")
    _, text = run(tmp_path, [0] * 10, "--min-depth", "15", "--keep-bed", str(keep))
    assert text == "ref\t1\t3\nref\t4\t8\n"
