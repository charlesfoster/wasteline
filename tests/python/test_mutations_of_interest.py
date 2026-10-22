from pathlib import Path

import mutations_of_interest as moi

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures" / "qc"
GFF = str(ROOT / "assets" / "references" / "NC_045512.gff3")
REF = str(ROOT / "assets" / "references" / "NC_045512.fasta")
SLICE = str(FIX / "real_variants_slice.tsv")


def write_moi(tmp_path, entries):
    path = tmp_path / "moi.csv"
    path.write_text("name,mutation\n" + "".join(f"{n},{m}\n" for n, m in entries))
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


def test_real_ivar_slice_snv_aa_and_indels(tmp_path):
    """Nucleotide SNVs, amino-acid changes resolved from the GFF, ivar-style indels and absent
    mutations, against real ivar output."""
    entries = [
        ("D614G", "A23403G"),
        ("N501Y", "A23063T"),
        ("P681R", "S:P681R"),
        ("ins", "A14+C"),  # ivar: POS 14 ALT +C
        ("del6", "AGTTCA22029-"),  # ivar: POS 22028 ALT -AGTTCA
        ("absent_aa", "S:N501Y"),
    ]
    rc, rows = run(tmp_path, SLICE, entries, "--gff", GFF, "--reference", REF)
    assert rc == 0
    by = {r["name"]: r for r in rows}
    assert by["D614G"] == {
        "sample": "S_1",
        "name": "D614G",
        "mutation": "A23403G",
        "ref_depth": "340",
        "alt_depth": "220",
        "depth": "560",
        "alt_freq": "0.392857",
    }
    assert by["N501Y"]["alt_freq"] == "0" and by["N501Y"]["depth"] == "NA"
    assert by["P681R"]["alt_freq"] == "0.490665"
    assert (by["ins"]["alt_freq"], by["ins"]["alt_depth"]) == ("0.579592", "142")
    assert (by["del6"]["alt_freq"], by["del6"]["depth"]) == ("0.389892", "277")
    assert by["absent_aa"]["alt_freq"] == "0"


def test_missing_variants_file_gives_na_rows(tmp_path):
    """A sample with no variants file still gets a row per mutation (NA), not an empty file."""
    rc, rows = run(tmp_path, str(tmp_path / "nope.tsv"), [("a", "A100T"), ("b", "S:D614G")])
    assert rc == 0
    assert all(r["alt_freq"] == "NA" and r["depth"] == "NA" for r in rows)
    assert [r["name"] for r in rows] == ["a", "b"]


def test_contig_selects_fasta_record_gff_seqid_and_depths(tmp_path):
    ref_seq = "".join(Path(REF).read_text().splitlines()[1:])
    fa = tmp_path / "multi.fa"
    fa.write_text(f">segA first\n{ref_seq}\n>segB\n{'A' * len(ref_seq)}\n")
    gff_lines = [
        line.replace("NC_045512.2", "segA", 1) for line in Path(GFF).read_text().splitlines()
    ]
    cds = next(line for line in gff_lines if "\tCDS\t" in line)
    gff = tmp_path / "multi.gff3"
    gff.write_text("\n".join([*gff_lines, cds.replace("segA", "segB", 1)]) + "\n")
    depths = tmp_path / "multi.depths"
    depths.write_text(
        "".join(f"segA\t{p}\tA\t{p}\n" for p in (23063, 23064, 23065)) + "segB\t23063\tA\t9\n"
    )
    entries = [("N501Y", "S:N501Y"), ("D614G", "A23403G")]
    multi = ["--gff", str(gff), "--reference", str(fa), "--depths", str(depths)]
    rc, rows = run(tmp_path, SLICE, entries, *multi, "--contig", "segA")
    assert rc == 0
    # N501Y absent (as in the single-contig test): its depth is the minimum over the codon (23063-23065) in segA, not segB's 9
    assert (rows[0]["alt_freq"], rows[0]["depth"]) == ("0", "23063")
    assert rows[1]["alt_freq"] == "0.392857"
    assert run(tmp_path, SLICE, entries, *multi)[0] == 2  # several contigs, no --contig
