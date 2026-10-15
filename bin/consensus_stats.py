#!/usr/bin/env python3
"""Consensus genome statistics (consensus_stats.tsv).

  consensus_stats.py --fasta S.consensus.fa --sample S [--reference-length 29903] -o S.stats.tsv
  consensus_stats.py --concat a.tsv b.tsv ... -o consensus_stats.tsv   (sorted by sample)

Columns: sample, length, n_count, n_pct, ambiguous_count, status.

  length           residues in the (single) FASTA record, whitespace ignored, multi-line FASTA supported
  n_count          N, n or '-' (deletion / masked) characters
  n_pct            n_count / denominator * 100 (2 dp); the denominator is --reference-length when given,
                   otherwise the sequence length. A length different from --reference-length only warns.
  ambiguous_count  characters other than A, C, G, T, N and '-' (IUPAC ambiguity codes and anything else)
  status           ok      sequence present with at least one non-N/- base
                   empty   sequence present but every character is N or '-'
                   failed  file missing/unreadable, not FASTA, no record, more than one record, or zero length
                           (all numeric fields are NA)
"""

import argparse
import contextlib
import sys

COLUMNS = ["sample", "length", "n_count", "n_pct", "ambiguous_count", "status"]
NA = "NA"
_ACGT = set("ACGT")


def read_single_fasta(path):
    """Return the sequence of the single FASTA record in `path`; ValueError otherwise."""
    try:
        with open(path) as handle:
            text = handle.read()
    except OSError as exc:
        raise ValueError(f"cannot read {path}: {exc}") from exc
    seqs = []
    current = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = []
            seqs.append(current)
        elif current is None:
            raise ValueError(f"{path}: sequence data before the first '>' header")
        else:
            current.append(line)
    if not seqs:
        raise ValueError(f"{path}: no FASTA record")
    if len(seqs) > 1:
        raise ValueError(f"{path}: {len(seqs)} FASTA records, expected exactly one")
    return "".join(seqs[0])


def compute_row(sample, fasta, reference_length=None):
    failed = [sample, NA, NA, NA, NA, "failed"]
    if not fasta:
        return failed
    try:
        seq = read_single_fasta(fasta).upper()
    except ValueError as exc:
        print(f"warning: {sample}: {exc}", file=sys.stderr)
        return failed
    length = len(seq)
    if length == 0:
        print(f"warning: {sample}: zero-length sequence", file=sys.stderr)
        return failed
    n_count = seq.count("N") + seq.count("-")
    ambiguous = sum(1 for c in seq if c not in _ACGT and c not in "N-")
    if reference_length and length != reference_length:
        print(
            f"warning: {sample}: consensus length {length} differs from reference length {reference_length}",
            file=sys.stderr,
        )
    denominator = reference_length or length
    status = "empty" if n_count == length else "ok"
    return [
        sample,
        str(length),
        str(n_count),
        f"{100.0 * n_count / denominator:.2f}",
        str(ambiguous),
        status,
    ]


def write_rows(path, rows):
    handle = sys.stdout if path == "-" else open(path, "w")  # noqa: SIM115
    with contextlib.nullcontext(handle) if path == "-" else handle as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


def concat(paths):
    rows = []
    for path in paths:
        with open(path) as handle:
            lines = [line.rstrip("\r\n") for line in handle if line.strip()]
        if not lines:
            raise ValueError(f"{path}: empty file")
        if lines[0].split("\t") != COLUMNS:
            raise ValueError(f"{path}: header does not match the consensus_stats contract")
        for line in lines[1:]:
            fields = line.split("\t")
            if len(fields) != len(COLUMNS):
                raise ValueError(f"{path}: row has {len(fields)} fields, expected {len(COLUMNS)}")
            rows.append(fields)
    samples = [r[0] for r in rows]
    dupes = sorted({s for s in samples if samples.count(s) > 1})
    if dupes:
        raise ValueError(f"duplicate samples: {', '.join(dupes)}")
    rows.sort(key=lambda r: r[0])
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--concat", nargs="*", metavar="TSV")
    p.add_argument("--fasta")
    p.add_argument("--sample")
    p.add_argument("--reference-length", type=int)
    p.add_argument("-o", "--output", default="-")
    args = p.parse_args(argv)

    if args.concat is not None:
        try:
            rows = concat(args.concat)
        except (OSError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        write_rows(args.output, rows)
        return 0
    if not args.sample:
        p.error("--sample is required (or use --concat)")
    write_rows(args.output, [compute_row(args.sample, args.fasta, args.reference_length)])
    return 0


if __name__ == "__main__":
    sys.exit(main())
