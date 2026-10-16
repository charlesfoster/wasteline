#!/usr/bin/env python
# authors: borrowed from nextflow viralrecon team
"""Convert an iVar / Freyja variants TSV to VCF.

  ivar_variants_to_vcf.py [--reference ref.fa] [-po] [-af X] FILE_IN FILE_OUT

With --reference (the pipeline always passes it) indel alleles are built from the REFERENCE, never
from the bases ivar prints. This matters because ivar run without a reference prints the deleted
bases as N (row `G  -NNNNNNNNN`); building REF from them gives `GNNNNNNNNN`, which does not match
the reference and makes `bcftools consensus` abort. Rules (ivar POS p is the 1-based anchor base,
i.e. the base BEFORE the indel):

  deletion  `-X`  L = len(X)      REF = ref[p-1 : p+L]   ALT = ref[p-1]
                                  (the printed bases are ignored, only their number is used; a
                                  deletion overrunning the contig end is truncated to the contig;
                                  nothing left to delete -> skipped)
  insertion `+X`                  REF = ref[p-1]         ALT = ref[p-1] + X
                                  (X containing anything but A/C/G/T, e.g. N, is skipped: an
                                  unresolved inserted base must not be written into a consensus)
  SNV                             REF/ALT as printed by ivar (a REF differing from the reference is
                                  only counted; CONSENSUS_BCFTOOLS then skips that record)

Skipped records and reference mismatches are counted on stderr. Without --reference the old
behaviour (REF/ALT from the printed bases) is kept and a loud warning is printed, because it
produces invalid records for ivar's `-N...` deletions.
"""

import argparse
import errno
import os
import re
import sys


def parse_args(args=None):
    Description = "Convert iVar variants tsv file to vcf format."
    Epilog = """Example usage: python ivar_variants_to_vcf.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input tsv file.")
    parser.add_argument("FILE_OUT", help="Full path to output vcf file.")
    parser.add_argument(
        "-po",
        "--pass_only",
        dest="PASS_ONLY",
        help="Only output variants that PASS all filters.",
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--reference",
        dest="REFERENCE",
        default=None,
        help="Reference FASTA. Indel REF/ALT alleles are built from it (strongly recommended; "
        "without it ivar's printed bases are used, which breaks for deletions printed as N).",
    )
    parser.add_argument(
        "-af",
        "--allele_freq_thresh",
        type=float,
        dest="ALLELE_FREQ_THRESH",
        default=0,
        help="Only output variants where allele frequency greater than this number (default: 0).",
    )

    return parser.parse_args(args)


def make_dir(path):
    if path:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def read_reference(path):
    """Return {contig name (first word of the header): upper-case sequence} for a FASTA file."""
    seqs = {}
    name = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0] if line[1:].split() else ""
                if name in seqs:
                    raise ValueError(f"{path}: duplicate contig {name}")
                seqs[name] = []
            elif name is None:
                raise ValueError(f"{path}: sequence data before the first '>' header")
            else:
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def build_alleles(contig, pos, ref_base, ivar_alt, reference=None):
    """Return (REF, ALT, var_type, note) for one TSV row, or (None, None, None, reason) to skip it.

    `note` is None, or a short tag counted on stderr ("ref_mismatch", "deletion_truncated").
    `reference` is {contig: sequence} or None (legacy behaviour, see the module docstring).
    """
    kind = ivar_alt[0]
    if kind not in "+-":
        if reference is not None and contig in reference:
            seq = reference[contig]
            if 1 <= pos <= len(seq) and seq[pos - 1] != ref_base.upper():
                return ref_base, ivar_alt, "SNP", "ref_mismatch"
        return ref_base, ivar_alt, "SNP", None
    if reference is None:
        if kind == "+":
            return ref_base, ref_base + ivar_alt[1:], "INS", None
        return ref_base + ivar_alt[1:], ref_base, "DEL", None
    seq = reference.get(contig)
    if seq is None:
        return None, None, None, "unknown_contig"
    if not 1 <= pos <= len(seq):
        return None, None, None, "pos_out_of_range"
    anchor = seq[pos - 1]
    note = None if anchor == ref_base.upper() else "ref_mismatch"
    bases = ivar_alt[1:]
    if kind == "+":
        if not bases or any(c not in "ACGTacgt" for c in bases):
            return None, None, None, "insertion_non_acgt"
        return anchor, anchor + bases.upper(), "INS", note
    length = len(bases)
    if length == 0:
        return None, None, None, "empty_deletion"
    if pos + length > len(seq):
        length = len(seq) - pos
        if length <= 0:
            return None, None, None, "deletion_past_contig_end"
        note = note or "deletion_truncated"
    return seq[pos - 1 : pos + length], anchor, "DEL", note


def ivar_variants_to_vcf(FileIn, FileOut, passOnly=False, minAF=0, reference=None):
    """`reference` is None or {contig: sequence} (see read_reference)."""
    if reference is None:
        print(
            "WARNING: no --reference given: indel alleles are built from ivar's printed bases, "
            "which is invalid for deletions printed as N (bcftools consensus will reject them)",
            file=sys.stderr,
        )
    counts = {}
    filename = os.path.splitext(FileIn)[0]
    header = (
        "##fileformat=VCFv4.2\n"
        "##source=iVar\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">\n'
        '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">\n'
        '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">\n'
        '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">\n'
        '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">\n'
        '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Depth of alternate base on reverse reads">\n'
        '##FORMAT=<ID=ALT_QUAL,Number=1,Type=Integer,Description="Mean quality of alternate base">\n'
        '##FORMAT=<ID=ALT_FREQ,Number=1,Type=Float,Description="Frequency of alternate base">\n'
    )
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + filename + "\n"

    varList = set()
    varCountDict = {"SNP": 0, "INS": 0, "DEL": 0}
    OutDir = os.path.dirname(FileOut)
    make_dir(OutDir)
    with open(FileOut, "w") as fout, open(FileIn) as f:
        fout.write(header)
        for line in f:
            if not re.match("REGION", line) and line.strip():
                line = re.split("\t", line.rstrip("\n"))
                CHROM = line[0]
                POS = line[1]
                ID = "."
                REF, ALT, var_type, note = build_alleles(
                    CHROM, int(POS), line[2], line[3], reference
                )
                if note is not None:
                    counts[note] = counts.get(note, 0) + 1
                if REF is None:
                    continue
                QUAL = "."
                pass_test = line[13]
                if pass_test == "TRUE":
                    FILTER = "PASS"
                else:
                    FILTER = "FAIL"
                INFO = "DP=" + line[11]
                FORMAT = "GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ"
                SAMPLE = (
                    "1:"
                    + line[4]
                    + ":"
                    + line[5]
                    + ":"
                    + line[6]
                    + ":"
                    + line[7]
                    + ":"
                    + line[8]
                    + ":"
                    + line[9]
                    + ":"
                    + line[10]
                )
                oline = (
                    CHROM
                    + "\t"
                    + POS
                    + "\t"
                    + ID
                    + "\t"
                    + REF
                    + "\t"
                    + ALT
                    + "\t"
                    + QUAL
                    + "\t"
                    + FILTER
                    + "\t"
                    + INFO
                    + "\t"
                    + FORMAT
                    + "\t"
                    + SAMPLE
                    + "\n"
                )
                writeLine = True
                if passOnly and FILTER != "PASS":
                    writeLine = False
                if float(line[10]) < minAF:
                    writeLine = False
                if (CHROM, POS, REF, ALT) in varList:
                    writeLine = False
                else:
                    varList.add((CHROM, POS, REF, ALT))
                if writeLine:
                    varCountDict[var_type] += 1
                    fout.write(oline)

    for note, n in sorted(counts.items()):
        print(f"ivar_variants_to_vcf: {n} record(s): {note}", file=sys.stderr)

    print(
        "ivar_variants_to_vcf: written "
        + ", ".join(f"{k}={v}" for k, v in sorted(varCountDict.items())),
        file=sys.stderr,
    )


def main(args=None):
    args = parse_args(args)
    reference = read_reference(args.REFERENCE) if args.REFERENCE else None
    ivar_variants_to_vcf(
        args.FILE_IN, args.FILE_OUT, args.PASS_ONLY, args.ALLELE_FREQ_THRESH, reference
    )


if __name__ == "__main__":
    sys.exit(main())
