#!/usr/bin/env python3
"""Simulate paired-end Illumina reads of a mixture of Freyja barcode lineages (a dev tool for validating a pathogen).

Each lineage is a pseudo-genome: the reference with every mutation of that lineage's row in the Freyja barcodes CSV applied
(single-nucleotide changes only; indels and ambiguous columns are skipped and counted). Read pairs are drawn from the lineages
in the given proportions, shotgun-style (uniform fragments) or amplicon-style (--amplicon-bed: one fragment per amplicon, equal
weight; amplicons are the LEFT/RIGHT primer pairs, named <scheme>_<n>_LEFT_<alt> / _RIGHT_<alt>, as in primalscheme BEDs;
--shear takes a random piece of each amplicon instead of its ends, for long amplicons).
Single-sequence references only. Writes <out>_R1.fastq.gz, <out>_R2.fastq.gz and <out>.truth.tsv (lineage, proportion, mutations).

  simulate_pathogen_mixture.py --reference ref.fasta --barcodes RSVa_barcodes.csv \
      --lineage RSVa-A.D.1 --lineage RSVa-A.3 --proportion 0.3 --proportion 0.7 --out sim/mix --depth 500
"""

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

import numpy as np

COMP = str.maketrans("ACGT", "TGCA")
SNV = re.compile(r"([ACGT])(\d+)([ACGT])")


def read_fasta(path):
    names, seq = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            names.append(line[1:].split()[0])
        else:
            seq.append(line.strip().upper())
    if len(names) != 1:
        sys.exit(f"{path}: need exactly one sequence, found {len(names)}")
    return names[0], "".join(seq)


def lineage_genome(ref, barcodes, lineage):
    """Reference with the lineage's barcode mutations applied; returns (genome, applied, skipped)."""
    with open(barcodes) as fh:
        rows = csv.reader(fh)
        header = next(rows)
        for row in rows:
            if row[0] == lineage:
                break
        else:
            sys.exit(f"{barcodes}: lineage {lineage} not found")
    genome, applied, skipped = list(ref), 0, 0
    for col, value in zip(header[1:], row[1:]):
        if float(value or 0) <= 0:
            continue
        m = SNV.fullmatch(col)
        if not m or ref[int(m[2]) - 1] != m[1]:
            skipped += 1
            continue
        genome[int(m[2]) - 1] = m[3]
        applied += 1
    return "".join(genome), applied, skipped


def amplicons(bed):
    """Amplicon (start, end) 0-based half-open intervals from the primer BED (outer ends of the LEFT/RIGHT primers)."""
    spans = {}
    for line in Path(bed).read_text().splitlines():
        f = line.split("\t")
        m = re.match(r"(.+?_\d+)_(LEFT|RIGHT)", f[3]) if len(f) > 3 else None
        if m:
            lo, hi = spans.get(m[1], (10**9, 0))
            spans[m[1]] = (min(lo, int(f[1])), max(hi, int(f[2])))
    return sorted(spans.values())


def mutate(codes, rng, rate):
    """Substitute each base with probability `rate` (a different base); codes are 0-3 uint8 arrays."""
    hit = rng.random(codes.shape) < rate
    codes[hit] = (codes[hit] + rng.integers(1, 4, hit.sum())) % 4
    return codes


def sample_reads(genome, a, n, rng):
    """n read pairs (R1, R2 strings) from one pseudo-genome: shotgun fragments or one fragment per amplicon."""
    if a.amplicon_bed:
        amps = amplicons(a.amplicon_bed)
        spans = [amps[i] for i in rng.integers(0, len(amps), n)]
        if a.shear:  # long amplicons (2.5 kb): the fragment is a random piece of the amplicon, not the whole amplicon
            sheared = []
            for lo, hi in spans:
                length = int(np.clip(rng.normal(a.insert_mean, 40), a.read_length, hi - lo))
                start = lo + int(rng.random() * (hi - lo - length + 1))
                sheared.append((start, start + length))
            spans = sheared
    else:
        lens = np.clip(rng.normal(a.insert_mean, 40, n).astype(int), a.read_length, len(genome))
        starts = (rng.random(n) * (len(genome) - lens + 1)).astype(int)
        spans = list(zip(starts, starts + lens))
    code = {c: i for i, c in enumerate("ACGT")}
    for s, e in spans:
        frag = genome[s:e]
        rl = min(a.read_length, len(frag))
        pair = []
        for r in (frag[:rl], frag[-rl:].translate(COMP)[::-1]):
            c = mutate(np.array([code.get(b, 0) for b in r], dtype=np.uint8), rng, a.error_rate)
            pair.append("".join("ACGT"[x] for x in c))
        yield s, pair


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--reference", required=True)
    p.add_argument("--barcodes", required=True, help="Freyja barcodes CSV (first column = lineage)")
    p.add_argument("--lineage", action="append", required=True)
    p.add_argument(
        "--proportion", action="append", type=float, required=True, help="one per --lineage"
    )
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--depth", type=float, default=500, help="mean depth (default 500)")
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--insert-mean", type=float, default=350, help="shotgun fragment mean (sd 40)")
    p.add_argument("--error-rate", type=float, default=0.002, help="substitution error per base")
    p.add_argument("--amplicon-bed", help="primer BED: amplicon fragments instead of shotgun ones")
    p.add_argument(
        "--shear",
        action="store_true",
        help="with --amplicon-bed: fragments are random pieces of the amplicons (insert-mean), as in tagmented long-amplicon libraries",
    )
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()
    if len(a.lineage) != len(a.proportion):
        sys.exit("need one --proportion per --lineage")
    rng = np.random.default_rng(a.seed)
    _, ref = read_fasta(a.reference)
    props = np.array(a.proportion) / sum(a.proportion)
    genomes = [lineage_genome(ref, a.barcodes, lin) for lin in a.lineage]
    which = rng.choice(len(genomes), round(a.depth * len(ref) / (2 * a.read_length)), p=props)
    counts = np.bincount(which, minlength=len(genomes))
    truth = ["lineage\tproportion\tmutations_applied\tmutations_skipped\tpairs\n"]
    for lin, pr, (_, na, ns), n in zip(a.lineage, props, genomes, counts):
        truth.append(f"{lin}\t{pr:.4f}\t{na}\t{ns}\t{n}\n")
    Path(a.out + ".truth.tsv").write_text("".join(truth))
    q = chr(33 + round(-10 * np.log10(max(a.error_rate, 1e-6) * 3 / 4)))  # Phred of the error rate
    with (
        gzip.open(f"{a.out}_R1.fastq.gz", "wt") as out1,
        gzip.open(f"{a.out}_R2.fastq.gz", "wt") as out2,
    ):
        for gi, (genome, _, _) in enumerate(genomes):
            for k, (s, (r1, r2)) in enumerate(sample_reads(genome, a, counts[gi], rng)):
                name = f"@sim{gi}_{k}_{s}"
                out1.write(f"{name}/1\n{r1}\n+\n{q * len(r1)}\n")
                out2.write(f"{name}/2\n{r2}\n+\n{q * len(r2)}\n")
    print("".join(truth), end="")


if __name__ == "__main__":
    main()
