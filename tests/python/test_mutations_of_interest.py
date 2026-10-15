from pathlib import Path

import mutations_of_interest as moi
import pytest

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures" / "qc"
GFF = str(ROOT / "assets" / "references" / "NC_045512.gff3")
REF = str(ROOT / "assets" / "references" / "NC_045512.fasta")
SLICE = str(FIX / "real_variants_slice.tsv")

HEADER = (
    "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\t"
    "GFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\tPOS_AA\n"
)


def write_moi(tmp_path, entries):
    path = tmp_path / "moi.csv"
    path.write_text("name,mutation\n" + "".join(f"{n},{m}\n" for n, m in entries))
    return str(path)


def write_variants(tmp_path, rows):
    """rows: (pos, ref, alt, ref_dp, alt_dp, freq[, feature, ref_aa, alt_aa, pos_aa])"""
    lines = [HEADER]
    for r in rows:
        pos, ref, alt, ref_dp, alt_dp, freq = r[:6]
        ann = r[6:] if len(r) > 6 else ("NA", "NA", "NA", "NA")
        feature, ref_aa, alt_aa, pos_aa = ann
        lines.append(
            f"NC_045512.2\t{pos}\t{ref}\t{alt}\t{ref_dp}\t1\t36\t{alt_dp}\t1\t37\t{freq}\t{ref_dp + alt_dp}\t0\tTRUE\t"
            f"{feature}\tNA\t{ref_aa}\tNA\t{alt_aa}\t{pos_aa}\n"
        )
    path = tmp_path / "v.tsv"
    path.write_text("".join(lines))
    return str(path)


def run(tmp_path, variants, entries, *extra):
    out = tmp_path / "out.tsv"
    rc = moi.main(
        [
            "--variants",
            variants,
            "--moi",
            write_moi(tmp_path, entries),
            "--sample",
            "S_1",
            "-o",
            str(out),
            *extra,
        ]
    )
    if rc != 0:
        return rc, None
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == moi.COLUMNS
    return rc, [dict(zip(moi.COLUMNS, line.split("\t"))) for line in lines[1:]]


def write_depths(tmp_path, overrides, n=29903, default=500):
    depths = [default] * n
    for pos, d in overrides.items():
        depths[pos - 1] = d
    path = tmp_path / "d.depths"
    path.write_text("".join(f"NC_045512.2\t{i}\tA\t{d}\n" for i, d in enumerate(depths, start=1)))
    return str(path)


# ---- parsing --------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "text,expected",
    [
        ("A23063T", {"kind": "snv", "ref": "A", "pos": 23063, "alt": "T"}),
        ("C21995-", {"kind": "del", "bases": "C", "pos": 21995}),
        ("TCT21991-", {"kind": "del", "bases": "TCT", "pos": 21991}),
        ("T22204+GAGCCAGAA", {"kind": "ins", "ref": "T", "pos": 22204, "bases": "GAGCCAGAA"}),
        ("S:N501Y", {"kind": "aa", "gene": "S", "ref": "N", "pos": 501, "alt": "Y"}),
        ("ORF1a:T3255I", {"kind": "aa", "gene": "ORF1a", "ref": "T", "pos": 3255, "alt": "I"}),
        ("N:R203*", {"kind": "aa", "gene": "N", "ref": "R", "pos": 203, "alt": "*"}),
    ],
)
def test_parse_mutation(text, expected):
    assert moi.parse_mutation(text) == expected


@pytest.mark.parametrize("text", ["", "N501Y", "S:501Y", "A23063", "23063T", "foo", "S:N501"])
def test_parse_mutation_rejects(text):
    with pytest.raises(ValueError):
        moi.parse_mutation(text)


def test_codon_table():
    assert moi.translate("ATG") == "M"
    assert moi.translate("TAA") == "*"
    assert moi.translate("TGG") == "W"
    assert moi.translate("GGG") == "G"
    assert moi.translate("NNN") == "X"
    assert len(moi.CODON_TABLE) == 64
    assert sum(1 for aa in moi.CODON_TABLE.values() if aa == "*") == 3


def test_cds_models_from_bundled_gff():
    models = moi.read_cds_models(GFF)
    assert len(models["s"]["positions"]) == 3822
    assert models["orf1a"]["positions"][0] == 266
    assert len(models["orf1a"]["positions"]) == 13218
    # ribosomal slippage: base 13468 is read twice in ORF1ab
    assert models["orf1ab"]["positions"].count(13468) == 2
    assert len(models["orf1ab"]["positions"]) == 21291


# ---- real TSV slice (ivar output from Freyja's bundled test.bam) ------------------------------------------


def test_real_slice_nucleotide_snv_and_absent(tmp_path):
    rc, rows = run(tmp_path, SLICE, [("D614G", "A23403G"), ("N501Y", "A23063T")])
    assert rc == 0
    assert rows[0] == {
        "sample": "S_1",
        "name": "D614G",
        "mutation": "A23403G",
        "ref_depth": "340",
        "alt_depth": "220",
        "depth": "560",
        "alt_freq": "0.392857",
    }
    assert rows[1]["alt_freq"] == "0" and rows[1]["alt_depth"] == "0"
    assert rows[1]["ref_depth"] == "NA" and rows[1]["depth"] == "NA"


def test_real_slice_amino_acid_fallback_from_gff(tmp_path):
    entries = [
        ("D614G", "S:D614G"),
        ("P681R", "S:P681R"),
        ("T478K", "S:T478K"),
        ("T19R", "S:T19R"),
        ("syn", "ORF1a:F924F"),
        ("absent", "S:N501Y"),
    ]
    rc, rows = run(tmp_path, SLICE, entries, "--gff", GFF, "--reference", REF)
    assert rc == 0
    by = {r["name"]: r for r in rows}
    assert (by["D614G"]["alt_depth"], by["D614G"]["depth"], by["D614G"]["alt_freq"]) == (
        "220",
        "560",
        "0.392857",
    )
    assert by["P681R"]["alt_freq"] == "0.490665"
    assert by["T478K"]["alt_freq"] == "0.496124"
    assert by["T19R"]["alt_depth"] == "24"
    assert by["syn"]["alt_freq"] == "0.396923"
    assert by["absent"]["alt_freq"] == "0"


def test_real_slice_indels(tmp_path):
    entries = [
        ("ins", "A14+C"),  # ivar: POS 14 ALT +C
        ("del1", "C75-"),  # ivar: POS 74 ALT -C (first deleted base is 75)
        ("del6", "AGTTCA22029-"),  # ivar: POS 22028 ALT -AGTTCA
        ("del_absent", "T22029-"),
    ]
    rc, rows = run(tmp_path, SLICE, entries)
    assert rc == 0
    assert rows[0]["alt_freq"] == "0.579592" and rows[0]["alt_depth"] == "142"
    assert rows[1]["alt_freq"] == "0.000416667" and rows[1]["alt_depth"] == "1"
    assert rows[2]["alt_freq"] == "0.389892" and rows[2]["depth"] == "277"
    assert rows[3]["alt_freq"] == "0"


def test_deletion_bases_match_the_bundled_reference():
    ref = moi.read_fasta_sequence(REF)
    assert ref[22028:22034] == "AGTTCA"  # first deleted base 22029 (1-based)
    assert ref[74] == "C"


# ---- semantics ----------------------------------------------------------------------------------------------


def test_absent_mutation_depth_from_depths_file(tmp_path):
    dep = write_depths(tmp_path, {23063: 73})
    rc, rows = run(
        tmp_path,
        SLICE,
        [("N501Y", "A23063T"), ("N501Y_aa", "S:N501Y")],
        "--depths",
        dep,
        "--gff",
        GFF,
        "--reference",
        REF,
    )
    assert rc == 0
    assert (rows[0]["ref_depth"], rows[0]["alt_depth"], rows[0]["depth"], rows[0]["alt_freq"]) == (
        "73",
        "0",
        "73",
        "0",
    )
    assert (rows[1]["depth"], rows[1]["alt_freq"]) == (
        "73",
        "0",
    )  # minimum over the codon's three positions


def test_snv_requires_matching_ref_and_alt(tmp_path):
    v = write_variants(tmp_path, [(100, "A", "T", 90, 10, 0.1)])
    _, rows = run(tmp_path, v, [("a", "A100T"), ("b", "A100G"), ("c", "C100T")])
    assert [r["alt_freq"] for r in rows] == ["0.1", "0", "0"]


def test_duplicate_rows_are_not_double_counted(tmp_path):
    v = write_variants(tmp_path, [(100, "A", "T", 90, 10, 0.1), (100, "A", "T", 90, 10, 0.1)])
    _, rows = run(tmp_path, v, [("a", "A100T")])
    assert rows[0]["alt_depth"] == "10"


def test_orf1ab_beyond_the_slippage_site(tmp_path):
    v = write_variants(tmp_path, [(14408, "C", "T", 10, 90, 0.9)])
    rc, rows = run(tmp_path, v, [("P314L", "ORF1ab:P4715L")], "--gff", GFF, "--reference", REF)
    assert rc == 0
    assert rows[0]["alt_freq"] == "0.9"
    # the same protein position does not exist in ORF1a (13218 nt = 4406 aa)
    rc, _ = run(tmp_path, v, [("x", "ORF1a:P4715L")], "--gff", GFF, "--reference", REF)
    assert rc == 2


def test_orf1a_position_in_orf1a(tmp_path):
    v = write_variants(tmp_path, [(10029, "C", "T", 20, 80, 0.8)])
    _, rows = run(
        tmp_path,
        v,
        [("a", "ORF1a:T3255I"), ("b", "ORF1ab:T3255I")],
        "--gff",
        GFF,
        "--reference",
        REF,
    )
    assert rows[0]["alt_freq"] == "0.8"
    assert rows[1]["alt_freq"] == "0.8"  # same CDS coordinates in ORF1ab before the slippage site


def test_stop_gain(tmp_path):
    ref = moi.read_fasta_sequence(REF)
    model = moi.read_cds_models(GFF)["n"]
    # first TGG (W) codon of N: G>A at the third base gives TGA (stop)
    pp = next(
        i
        for i in range(1, 420)
        if "".join(ref[p - 1] for p in moi.codon_positions(model, i)) == "TGG"
    )
    third = moi.codon_positions(model, pp)[2]
    v = write_variants(tmp_path, [(third, "G", "A", 50, 50, 0.5)])
    _, rows = run(tmp_path, v, [("stop", f"N:W{pp}*")], "--gff", GFF, "--reference", REF)
    assert rows[0]["alt_freq"] == "0.5"


def test_two_nucleotide_changes_to_same_aa_are_summed(tmp_path):
    ref = moi.read_fasta_sequence(REF)
    models = moi.read_cds_models(GFF)
    cps = moi.codon_positions(models["s"], 346)
    codon = "".join(ref[p - 1] for p in cps)
    assert moi.translate(codon) == "R"  # AGA
    # AGA -> AGG is synonymous R; -> CGA is R as well: two rows, same aa change R346R
    v = write_variants(tmp_path, [(cps[2], "A", "G", 50, 30, 0.3), (cps[0], "A", "C", 50, 20, 0.2)])
    _, rows = run(tmp_path, v, [("syn", "S:R346R")], "--gff", GFF, "--reference", REF)
    assert rows[0]["alt_depth"] == "50"
    assert rows[0]["alt_freq"] == "0.5"


def test_annotated_tsv_resolves_aa_without_gff(tmp_path):
    v = write_variants(
        tmp_path,
        [
            (23063, "A", "T", 10, 90, 0.9, "S", "N", "Y", 501),
            (23403, "A", "G", 60, 40, 0.4, "S", "D", "G", 614),
            (23063, "A", "C", 90, 5, 0.05),
        ],
    )
    rc, rows = run(tmp_path, v, [("N501Y", "S:N501Y"), ("D614G", "S:D614G"), ("K417N", "S:K417N")])
    assert rc == 0
    assert (rows[0]["alt_freq"], rows[0]["alt_depth"], rows[0]["depth"]) == ("0.9", "90", "100")
    assert rows[1]["alt_freq"] == "0.4"
    assert rows[2]["alt_freq"] == "0" and rows[2]["depth"] == "NA"


def test_annotation_for_another_feature_is_ignored(tmp_path):
    v = write_variants(tmp_path, [(23063, "A", "T", 10, 90, 0.9, "ORF3a", "N", "Y", 501)])
    rc, rows = run(tmp_path, v, [("x", "S:N501Y")], "--gff", GFF, "--reference", REF)
    assert rc == 0
    assert rows[0]["alt_freq"] == "0.9"  # recomputed from nucleotides, which really are S 501


def test_aa_needs_gff_when_tsv_unannotated(tmp_path):
    rc, rows = run(tmp_path, SLICE, [("x", "S:D614G")])
    assert rc == 2 and rows is None


def test_unknown_gene_is_an_error(tmp_path):
    rc, _ = run(tmp_path, SLICE, [("x", "Q:D614G")], "--gff", GFF, "--reference", REF)
    assert rc == 2


def test_protein_position_beyond_cds_is_an_error(tmp_path):
    rc, _ = run(tmp_path, SLICE, [("x", "E:V9999A")], "--gff", GFF, "--reference", REF)
    assert rc == 2


def test_reference_residue_mismatch_warns(tmp_path, capsys):
    rc, _rows = run(tmp_path, SLICE, [("x", "S:G614D")], "--gff", GFF, "--reference", REF)
    assert rc == 0
    assert "reference residue is D" in capsys.readouterr().err


def test_bad_mutation_in_moi_file_is_an_error(tmp_path):
    rc, _ = run(tmp_path, SLICE, [("x", "garbage")])
    assert rc == 2


def test_moi_file_needs_expected_header(tmp_path):
    p = tmp_path / "moi.csv"
    p.write_text("a,b\n1,2\n")
    assert (
        moi.main(["--variants", SLICE, "--moi", str(p), "--sample", "S", "-o", str(tmp_path / "o")])
        == 2
    )


def test_missing_variants_file_gives_na_rows(tmp_path):
    rc, rows = run(tmp_path, str(tmp_path / "nope.tsv"), [("a", "A100T"), ("b", "S:D614G")])
    assert rc == 0
    assert all(r["alt_freq"] == "NA" and r["depth"] == "NA" for r in rows)
    assert [r["name"] for r in rows] == ["a", "b"]


def test_header_only_variants_gives_absent_rows(tmp_path):
    v = write_variants(tmp_path, [])
    _, rows = run(tmp_path, v, [("a", "A100T")])
    assert rows[0]["alt_freq"] == "0" and rows[0]["alt_depth"] == "0"


def test_ivar_tsv_without_pos_aa_column_is_accepted(tmp_path):
    v = tmp_path / "old.tsv"
    v.write_text(
        "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"
        "NC_045512.2\t100\tA\tT\t90\t1\t36\t10\t1\t37\t0.1\t100\t0\tTRUE\n"
    )
    _, rows = run(tmp_path, str(v), [("a", "A100T")])
    assert rows[0]["alt_freq"] == "0.1"
