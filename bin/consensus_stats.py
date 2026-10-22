#!/usr/bin/env python3
"""Consensus genome statistics (consensus_stats.tsv).

  consensus_stats.py --fasta S.consensus.fa --sample S [--reference-length 29903]
                     [--status-json S.consensus_status.json] [--contig-out S.stats_contigs.tsv] -o S.stats.tsv
  consensus_stats.py --concat a.tsv b.tsv ... -o consensus_stats.tsv   (sorted by sample)

Columns: sample, length, n_count, consensus_coverage_pct, ambiguous_count, status, reason.

  length           residues in the (single) FASTA record, whitespace ignored, multi-line FASTA supported
  n_count          N, n or '-' (deletion / masked) characters
  consensus_coverage_pct
                   % of the reference with a called base (A, C, G, T or an IUPAC code): 100 - n_count /
                   denominator * 100 (2 dp, never below 0); the denominator is --reference-length when given,
                   otherwise the sequence length. A length different from --reference-length only warns.
                   0.00 for an `empty` consensus, NA for a `failed` one.
  ambiguous_count  characters other than A, C, G, T, N and '-' (IUPAC ambiguity codes and anything else)
  status           ok      sequence present with at least one non-N/- base
                   empty   sequence present but every character is N or '-'
                   failed  file missing/unreadable, not FASTA, no record, more than one record, or zero length
                           (all numeric fields are NA)
  reason           NA, or the text explaining a non-ok status or the warnings of an ok one (e.g. variants
                   skipped by CONSENSUS_BCFTOOLS). Tabs and newlines are replaced by spaces.

--contig-out FILE accepts a FASTA with several records (multi-contig references, e.g. the eight influenza
segments; without the flag more than one record stays an error). The sample row is then aggregated over all
records with the same formulas (`length` = total, `n_count` = sum, `ambiguous_count` = sum, consensus_coverage_pct
from the totals, status on the whole) and FILE gets one row per record with columns sample, contig, length,
n_count, consensus_coverage_pct (against that record's own length), ambiguous_count, status. `contig` is the
part of the record's first header word after the last `|` (CONSENSUS_BCFTOOLS writes `>sample|contig`). A
record's status is computed from its sequence (`ok`, `empty`, `failed` for zero length); a `failed` or `empty`
status in --status-json applies to the whole sample and so to every record. When the FASTA is missing or unusable
FILE gets one row with contig NA, NA numbers and status `failed`.

--status-json is the <sample>.consensus_status.json written by CONSENSUS_BCFTOOLS ({"status", "reason",
"warnings"}). When present and readable its status wins over the one recomputed from the sequence for
`failed` and `empty` (a failed sample is stored as an all-N sequence, which alone would read as `empty`);
for `ok` the sequence-based status is used, so an unreadable FASTA is still `failed`. Its reason is
reported; the numeric columns are always computed from the FASTA when it is readable.
"""

import argparse
import contextlib
import json
import sys

COLUMNS = [
    "sample",
    "length",
    "n_count",
    "consensus_coverage_pct",
    "ambiguous_count",
    "status",
    "reason",
]
CONTIG_COLUMNS = [
    "sample",
    "contig",
    "length",
    "n_count",
    "consensus_coverage_pct",
    "ambiguous_count",
    "status",
]
LEGACY_COLUMNS = COLUMNS[:-1]  # pre-`reason` per-sample files are still accepted by --concat
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


def read_fasta_records(path):
    """Return [(header, sequence)] for every record of the FASTA `path`; ValueError when unreadable or empty."""
    try:
        with open(path) as handle:
            text = handle.read()
    except OSError as exc:
        raise ValueError(f"cannot read {path}: {exc}") from exc
    records = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            records.append((line[1:].strip(), []))
        elif not records:
            raise ValueError(f"{path}: sequence data before the first '>' header")
        else:
            records[-1][1].append(line)
    if not records:
        raise ValueError(f"{path}: no FASTA record")
    return [(header, "".join(parts)) for header, parts in records]


def contig_name(header):
    """Part of the first header word after the last `|` (the whole word without one); NA when empty."""
    words = header.split()
    return words[0].rsplit("|", 1)[-1] if words else NA


def clean(text):
    """One-line, tab-free text; NA when empty."""
    text = " ".join(str(text).split()) if text is not None else ""
    return text or NA


def read_status_json(path):
    """Return (status, reason) from a consensus_status.json, or (None, None) when absent/unusable."""
    if not path:
        return None, None
    try:
        with open(path) as handle:
            data = json.load(handle)
        status = str(data["status"]).strip().lower()
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        print(f"warning: cannot use status file {path}: {exc}", file=sys.stderr)
        return None, None
    if status not in ("ok", "empty", "failed"):
        print(f"warning: {path}: unknown status {status!r}, ignored", file=sys.stderr)
        return None, None
    return status, data.get("reason") or ""


def seq_counts(seq):
    """(n_count, ambiguous_count) of an upper-case sequence."""
    n_count = seq.count("N") + seq.count("-")
    ambiguous = sum(1 for c in seq if c not in _ACGT and c not in "N-")
    return n_count, ambiguous


def compute_row(sample, fasta, reference_length=None, status_json=None, contig_rows=None):
    """Return the consensus_stats row of one sample.

    With `contig_rows` (a list) several FASTA records are accepted: the row is aggregated over them and one
    row per record is appended to `contig_rows` (see the module docstring). Without it more than one record
    is a `failed` row.
    """
    json_status, json_reason = read_status_json(status_json)

    def failed_row(why):
        reason = json_reason if json_status == "failed" and json_reason else why
        if contig_rows is not None:
            contig_rows.append([sample, NA, NA, NA, NA, NA, "failed"])
        return [sample, NA, NA, NA, NA, "failed", clean(reason)]

    if not fasta:
        return failed_row("no consensus file")
    try:
        if contig_rows is None:
            records = [(sample, read_single_fasta(fasta))]
        else:
            records = read_fasta_records(fasta)
    except ValueError as exc:
        print(f"warning: {sample}: {exc}", file=sys.stderr)
        return failed_row(str(exc))
    records = [(header, seq.upper()) for header, seq in records]
    length = sum(len(seq) for _, seq in records)
    if length == 0:
        print(f"warning: {sample}: zero-length sequence", file=sys.stderr)
        return failed_row("zero-length consensus sequence")
    n_count = ambiguous = 0
    for _, seq in records:
        n, a = seq_counts(seq)
        n_count += n
        ambiguous += a
    if reference_length and length != reference_length:
        print(
            f"warning: {sample}: consensus length {length} differs from reference length {reference_length}",
            file=sys.stderr,
        )
    denominator = reference_length or length
    status = "empty" if n_count == length else "ok"
    reason = NA
    if status == "empty":
        reason = "consensus contains only N or - characters"
    if json_status in ("failed", "empty"):
        status = json_status
        reason = json_reason or ("consensus failed" if json_status == "failed" else reason)
    elif json_status == "ok" and json_reason and status == "ok":
        reason = json_reason  # e.g. variants skipped
    if contig_rows is not None:
        for header, seq in records:
            c_len = len(seq)
            c_n, c_amb = seq_counts(seq)
            if c_len == 0:
                c_row = [sample, contig_name(header), NA, NA, NA, NA, "failed"]
            else:
                c_status = json_status if json_status in ("failed", "empty") else None
                c_status = c_status or ("empty" if c_n == c_len else "ok")
                c_row = [
                    sample,
                    contig_name(header),
                    str(c_len),
                    str(c_n),
                    f"{max(0.0, 100.0 * (c_len - c_n) / c_len):.2f}",
                    str(c_amb),
                    c_status,
                ]
            contig_rows.append(c_row)
    return [
        sample,
        str(length),
        str(n_count),
        f"{max(0.0, 100.0 * (denominator - n_count) / denominator):.2f}",
        str(ambiguous),
        status,
        clean(reason),
    ]


def write_rows(path, rows, columns=COLUMNS):
    handle = sys.stdout if path == "-" else open(path, "w")  # noqa: SIM115
    with contextlib.nullcontext(handle) if path == "-" else handle as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


def concat(paths):
    rows = []
    for path in paths:
        with open(path) as handle:
            lines = [line.rstrip("\r\n") for line in handle if line.strip()]
        if not lines:
            raise ValueError(f"{path}: empty file")
        header = lines[0].split("\t")
        if header not in (COLUMNS, LEGACY_COLUMNS):
            raise ValueError(f"{path}: header does not match the consensus_stats contract")
        for line in lines[1:]:
            fields = line.split("\t")
            if len(fields) != len(header):
                raise ValueError(f"{path}: row has {len(fields)} fields, expected {len(header)}")
            if header == LEGACY_COLUMNS:
                fields.append(NA)
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
    p.add_argument("--status-json", help="<sample>.consensus_status.json from CONSENSUS_BCFTOOLS")
    p.add_argument(
        "--contig-out",
        help="accept a multi-record FASTA (one record per contig); write the per-contig table here",
    )
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
    contig_rows = [] if args.contig_out else None
    row = compute_row(args.sample, args.fasta, args.reference_length, args.status_json, contig_rows)
    write_rows(args.output, [row])
    if args.contig_out:
        write_rows(args.contig_out, contig_rows, CONTIG_COLUMNS)
    return 0


if __name__ == "__main__":
    sys.exit(main())
