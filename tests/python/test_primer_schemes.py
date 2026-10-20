"""Every bundled primer scheme must use the contig name of the bundled reference.

articv3, articv5.3.2 and qiagenv1 originally used MN908947.3 (same genome and coordinates as NC_045512.2)
and failed the start-up contig check against the bundled reference; the contig is now NC_045512.2.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCHEMES = sorted((ROOT / "assets" / "primer_schemes").glob("*.bed"))
REFERENCE = ROOT / "assets" / "references" / "NC_045512.fasta"


def reference_contig():
    header = REFERENCE.read_text().splitlines()[0]
    return header[1:].split()[0]


def test_every_bed_uses_the_reference_contig():
    assert len(SCHEMES) == 9
    for bed in SCHEMES:
        contigs = {line.split("\t")[0] for line in bed.read_text().splitlines() if line.strip()}
        assert contigs == {reference_contig()}, f"{bed.name}: {contigs}"
