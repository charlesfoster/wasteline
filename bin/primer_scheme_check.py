#!/usr/bin/env python3
"""Does the declared primer scheme match the reads? A statistical check on read 5' ends.

  primer_scheme_check.py --starts S.read_starts.tsv --declared SCHEME_OR_BED --sample S
        [--beds-dir assets/primer_schemes] [--window 2] [--min-enrichment 5] [--min-reads 500]
        [--genome-length 29903] [--max-p-value 0.001] [--min-hit-fraction 0.05]
        [--equivalence-jaccard 0.9]
        -o S.primer_scheme_check.tsv [--summary-out S.primer_scheme_summary.tsv]

Idea. In anchored libraries (ligation/ARTIC-style, QIAseq single-primer extension) a read starts exactly at a
primer's 5' end, on the primer's strand. For every candidate scheme (each `*.bed` in --beds-dir, plus the declared
one if it is a BED path) count the mapped primary reads whose 5' alignment end lies within +-W bp of a primer 5'
end of the SAME strand, and compare with the number expected if reads started at random.

`--starts` is a TSV with a header and the columns `strand` (`+`/`-`) and `five_prime_pos` (1-based reference
coordinate of the read's 5' aligned base: the alignment start for forward reads, the last aligned base for reverse
reads; soft clips are ignored). The READ_STARTS module writes it with samtools. Single-contig references only
(contig names are not compared).

Primer 5' ends (1-based): `+` primer -> BED start + 1; `-` primer -> BED end. The strand is BED column 6, or, when
that column is absent/invalid (xgen), the primer name suffix (`F`/`LEFT` -> `+`, `R`/`RIGHT` -> `-`, an `_alt<N>` tail
is ignored). Other columns (swift's column 7) are ignored.

Statistics, per strand s with n_s reads: the window set W_s = every position within +-W of a primer 5' end of that
strand (its size is counted, so overlapping windows of dense schemes are not double counted);
expected_s = n_s * |W_s| / genome_length; hits_s = reads whose 5' position is in W_s. Over both strands
hits = sum hits_s, expected = sum expected_s, enrichment = hits / expected (NA if expected is 0), hit_fraction =
hits / n_reads. p_value = one-sided P(X >= hits) for X ~ Binomial(n_reads, expected / n_reads): the pooled
binomial, which is slightly conservative compared with the exact sum of the two strands' binomials (equal when both
strands have the same window fraction). It is computed exactly in log space (terms summed outward from `hits` until
they are negligible), reported with 3 significant digits and floored at 1e-300.

Scheme comparison. Two schemes are *equivalent* when the Jaccard index of their sets of (strand, 5' position) sites
is >= --equivalence-jaccard (default 0.9; e.g. midnight and midnightv2 share 58 of 59 sites): reads cannot tell
them apart. A scheme "reaches" the bar when ALL of: enrichment >= --min-enrichment, p_value <= --max-p-value (stops
a handful of chance hits passing on a ratio alone) and hit_fraction >= --min-hit-fraction (an absolute floor: a sparse
scheme has a tiny expected count, so even ~1.5 % of reads at its sites gives a large enrichment ratio; 0 disables).

Summary row (`--summary-out`): sample, declared_scheme, declared_enrichment, best_scheme, best_enrichment,
scheme_check, reason, then the appended declared_hit_fraction, best_hit_fraction (4 dp; consumers read by column name
and ignore these). best_scheme = the highest-enrichment scheme that reaches the bar (NA if none; ties/near-ties of
equivalent schemes are all listed, joined with `/`; best_hit_fraction is that top scheme's); scheme_check:
  ok            the declared scheme reaches the bar (reason `declared_enriched`)
  mismatch      the declared scheme does not reach the bar, the best scheme does (so its hit_fraction is also
                >= --min-hit-fraction), is not equivalent to the declared one and either has >= 3x the declared
                enrichment or >= 3x the declared hit_fraction (reason `best_other_scheme`)
  uninformative no scheme reaches the bar, e.g. fragmentation/tagmentation libraries (reason `no_scheme_enriched`),
                or the declared scheme does not reach it and the best scheme is not decisively better
                (`best_not_decisive`), or the declared scheme passes the enrichment and p-value tests but not the
                hit-fraction floor and nothing else reaches the bar (`declared_sparse_inflated`)
  NA            fewer than --min-reads usable reads (`too_few_reads`), unreadable/missing starts
                (`missing_input`) or an unknown declared scheme (`declared_not_found`)
Never an error: exit 0 always except unusable arguments (2).
"""

import argparse
import glob
import math
import os
import re
import sys
from collections import Counter

NA = "NA"
DETAIL_COLUMNS = [
    "scheme",
    "n_primer_5prime_sites",
    "n_reads",
    "hits",
    "expected_hits",
    "hit_fraction",
    "enrichment",
    "p_value",
    "equivalent_to_declared",
]
SUMMARY_COLUMNS = [
    "sample",
    "declared_scheme",
    "declared_enrichment",
    "best_scheme",
    "best_enrichment",
    "scheme_check",
    "reason",
    "declared_hit_fraction",
    "best_hit_fraction",
]
NAME_STRAND = re.compile(r"(LEFT|RIGHT|F|R)$")
ALT_TAIL = re.compile(r"[_-]?ALT\d*$")


def primer_strand(fields):
    """'+' or '-' for a BED line (already split), or None."""
    if len(fields) > 5 and fields[5] in ("+", "-"):
        return fields[5]
    name = ALT_TAIL.sub("", fields[3].upper()) if len(fields) > 3 else ""
    m = NAME_STRAND.search(name)
    if not m:
        return None
    return "+" if m.group(1) in ("F", "LEFT") else "-"


def read_bed_sites(path):
    """Return (set of (strand, 1-based 5' position), number of skipped lines)."""
    sites = set()
    skipped = 0
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            try:
                start, end = int(fields[1]), int(fields[2])
            except (IndexError, ValueError):
                skipped += 1
                continue
            strand = primer_strand(fields)
            if strand is None:
                skipped += 1
                continue
            sites.add((strand, start + 1 if strand == "+" else end))
    return sites, skipped


def read_starts(path):
    """Return Counter of (strand, pos); bad lines are skipped. ValueError if the file has no usable header."""
    counts = Counter()
    with open(path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if "strand" not in header or "five_prime_pos" not in header:
            raise ValueError(f"{path}: header must contain strand and five_prime_pos")
        i_s, i_p = header.index("strand"), header.index("five_prime_pos")
        top = max(i_s, i_p)
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) <= top or fields[i_s] not in ("+", "-"):
                continue
            try:
                counts[(fields[i_s], int(fields[i_p]))] += 1
            except ValueError:
                continue
    return counts


def log_pmf(k, n, p):
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + k * math.log(p)
        + (n - k) * math.log1p(-p)
    )


def binom_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p), exact in log space. Returns a float in [1e-300, 1]."""
    if k <= 0 or p >= 1:
        return 1.0
    if k > n or p <= 0:
        return 1e-300
    mode = int((n + 1) * p)
    # Terms fall monotonically away from the mode, so sum outward from k and stop once negligible.
    if k > mode:
        indices, complement = range(k, n + 1), False  # upper tail directly
    else:
        indices, complement = range(k - 1, -1, -1), True  # 1 - lower tail (k-1 < mode)
    acc = None
    for i in indices:
        term = log_pmf(i, n, p)
        if acc is None:
            acc = term
        else:
            acc = max(acc, term) + math.log1p(math.exp(-abs(acc - term)))
            if term < acc - 40:
                break
    value = 1.0 - math.exp(acc) if complement else math.exp(acc)
    return min(1.0, max(1e-300, value))


def window_set(positions, window):
    out = set()
    for pos in positions:
        out.update(range(max(1, pos - window), pos + window + 1))
    return out


def score_scheme(name, sites, counts, window, genome_length):
    row = {"scheme": name, "n_primer_5prime_sites": len(sites)}
    n_reads = hits = 0
    expected = 0.0
    for strand in ("+", "-"):
        windows = window_set((pos for s, pos in sites if s == strand), window)
        n_s = sum(c for (s, _), c in counts.items() if s == strand)
        n_reads += n_s
        hits += sum(c for (s, pos), c in counts.items() if s == strand and pos in windows)
        expected += n_s * len(windows) / genome_length
    row.update({"n_reads": n_reads, "hits": hits, "expected": expected})
    row["enrichment"] = hits / expected if expected > 0 else None
    row["hit_fraction"] = hits / n_reads if n_reads else None
    row["p_value"] = (
        binom_sf(hits, n_reads, min(expected / n_reads, 1.0)) if n_reads and expected > 0 else None
    )
    row["sites"] = sites
    return row


def jaccard(a, b):
    union = len(a | b)
    return len(a & b) / union if union else 0.0


def fmt(value, digits):
    return NA if value is None else f"{value:.{digits}f}"


def load_schemes(beds_dir, declared):
    """Return (dict name -> sites, declared_name or None, notes)."""
    schemes, notes = {}, []
    if beds_dir and os.path.isdir(beds_dir):
        for path in sorted(glob.glob(os.path.join(beds_dir, "*.bed"))):
            name = os.path.splitext(os.path.basename(path))[0]
            sites, skipped = read_bed_sites(path)
            if skipped:
                notes.append(f"{name}: {skipped} BED lines without coordinates/strand skipped")
            schemes[name] = sites
    declared_name = None
    if declared:
        if os.path.isfile(declared):
            declared_name = os.path.splitext(os.path.basename(declared))[0]
            schemes[declared_name], skipped = read_bed_sites(declared)
            if skipped:
                notes.append(f"{declared_name}: {skipped} BED lines skipped")
        elif declared in schemes:
            declared_name = declared
    return schemes, declared_name, notes


def enriched(row, args):
    """Enrichment and p-value tests only (the ratio evidence)."""
    return (
        row["enrichment"] is not None
        and row["enrichment"] >= args.min_enrichment
        and row["p_value"] is not None
        and row["p_value"] <= args.max_p_value
    )


def reaches(row, args):
    """The bar: ratio evidence AND an absolute hit-fraction floor."""
    return (
        enriched(row, args)
        and row["hit_fraction"] is not None
        and row["hit_fraction"] >= args.min_hit_fraction
    )


def summarise(args, scored, declared_name, status, reason):
    out = {c: NA for c in SUMMARY_COLUMNS}
    out["sample"] = args.sample
    out["declared_scheme"] = declared_name or args.declared or NA
    out["scheme_check"] = status
    out["reason"] = reason
    if scored and declared_name:
        by_name = {r["scheme"]: r for r in scored}
        out["declared_enrichment"] = fmt(by_name[declared_name]["enrichment"], 1)
        out["declared_hit_fraction"] = fmt(by_name[declared_name]["hit_fraction"], 4)
    return out


def check(args, notes):
    """Return (detail_rows, summary_dict)."""
    schemes, declared_name, load_notes = load_schemes(args.beds_dir, args.declared)
    notes.extend(load_notes)
    if declared_name is None:
        notes.append(
            f"declared scheme {args.declared!r} is neither a BED file nor in {args.beds_dir}"
        )
        return [], summarise(args, None, None, NA, "declared_not_found")
    try:
        counts = read_starts(args.starts)
    except (OSError, ValueError) as exc:
        notes.append(str(exc))
        return [], summarise(args, None, declared_name, NA, "missing_input")

    scored = [
        score_scheme(name, sites, counts, args.window, args.genome_length)
        for name, sites in sorted(schemes.items())
    ]
    declared_sites = schemes[declared_name]
    for row in scored:
        row["equivalent"] = (
            row["scheme"] == declared_name
            or jaccard(row["sites"], declared_sites) >= args.equivalence_jaccard
        )

    n_reads = scored[0]["n_reads"]
    if n_reads < args.min_reads:
        return scored, summarise(args, scored, declared_name, NA, "too_few_reads")

    by_name = {r["scheme"]: r for r in scored}
    declared_row = by_name[declared_name]
    qualifying = sorted(
        (r for r in scored if reaches(r, args)), key=lambda r: r["enrichment"], reverse=True
    )
    best_label, best_enrichment, best_hit_fraction, best_equivalent = NA, None, None, False
    if qualifying:
        top = qualifying[0]
        # the best scheme together with every scheme the reads cannot tell apart from it
        group = [
            r for r in qualifying if jaccard(r["sites"], top["sites"]) >= args.equivalence_jaccard
        ]
        best_label = "/".join(sorted(r["scheme"] for r in group))
        best_enrichment = top["enrichment"]
        best_hit_fraction = top["hit_fraction"]
        best_equivalent = any(r["equivalent"] for r in group)
    if reaches(declared_row, args):
        status, reason = "ok", "declared_enriched"
    elif not qualifying or best_equivalent:
        # the declared scheme's ratio looked fine but too few reads sit at its sites (sparse scheme)
        sparse = enriched(declared_row, args)
        status = "uninformative"
        reason = "declared_sparse_inflated" if sparse else "no_scheme_enriched"
    elif (
        declared_row["enrichment"] is None
        or best_enrichment >= 3 * declared_row["enrichment"]
        or declared_row["hit_fraction"] is None
        or best_hit_fraction >= 3 * declared_row["hit_fraction"]
    ):
        status, reason = "mismatch", "best_other_scheme"
    else:
        status, reason = "uninformative", "best_not_decisive"
    out = summarise(args, scored, declared_name, status, reason)
    out["best_scheme"] = best_label
    out["best_enrichment"] = fmt(best_enrichment, 1)
    out["best_hit_fraction"] = fmt(best_hit_fraction, 4)
    return scored, out


def write_tsv(path, columns, rows):
    text = (
        "\t".join(columns) + "\n" + "".join("\t".join(r[c] for c in columns) + "\n" for r in rows)
    )
    if path == "-":
        sys.stdout.write(text)
    else:
        with open(path, "w") as out:
            out.write(text)


def detail_row(row):
    return {
        "scheme": row["scheme"],
        "n_primer_5prime_sites": str(row["n_primer_5prime_sites"]),
        "n_reads": str(row["n_reads"]),
        "hits": str(row["hits"]),
        "expected_hits": fmt(row["expected"], 2),
        "hit_fraction": fmt(row["hit_fraction"], 4),
        "enrichment": fmt(row["enrichment"], 2),
        "p_value": NA if row["p_value"] is None else f"{row['p_value']:.3g}",
        "equivalent_to_declared": "true" if row["equivalent"] else "false",
    }


def default_beds_dir():
    path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "assets", "primer_schemes"
    )
    return path if os.path.isdir(path) else None


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--starts", required=True, help="read_starts.tsv (strand, five_prime_pos)")
    p.add_argument(
        "--declared", required=True, help="declared scheme name (in --beds-dir) or a BED path"
    )
    p.add_argument("--sample", default=NA)
    p.add_argument("--beds-dir", default=default_beds_dir(), help="directory of <scheme>.bed files")
    p.add_argument("--window", type=int, default=2, help="+- bp around a primer 5' end (default 2)")
    p.add_argument("--min-enrichment", type=float, default=5.0)
    p.add_argument("--min-reads", type=int, default=500)
    p.add_argument("--genome-length", type=int, default=29903, help="default 29903 (SARS-CoV-2)")
    p.add_argument("--max-p-value", type=float, default=0.001)
    p.add_argument(
        "--min-hit-fraction",
        type=float,
        default=0.05,
        help="a scheme reaches the bar only if >= this fraction of reads start at its primers "
        "(default 0.05; 0 disables)",
    )
    p.add_argument("--equivalence-jaccard", type=float, default=0.9)
    p.add_argument("-o", "--output", default="-", help="per-scheme table")
    p.add_argument("--summary-out", help="one-row summary table")
    args = p.parse_args(argv)
    if args.window < 0 or args.genome_length < 1:
        p.error("--window must be >= 0 and --genome-length >= 1")

    notes = []
    scored, summary = check(args, notes)
    for note in notes:
        print(f"note: {args.sample}: {note}", file=sys.stderr)
    write_tsv(args.output, DETAIL_COLUMNS, [detail_row(r) for r in scored])
    if args.summary_out:
        write_tsv(args.summary_out, SUMMARY_COLUMNS, [summary])
    return 0


if __name__ == "__main__":
    sys.exit(main())
