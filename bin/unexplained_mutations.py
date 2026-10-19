#!/usr/bin/env python3
"""Unexplained mutations: high-frequency mutations that the demixed lineage mixture does not account for.

An early-warning signal for a novel / cryptic lineage in a wastewater sample. This is NOT a variant caller:
it only asks, for every SNV that ivar saw at a high frequency, whether the lineages Freyja assigned (and their
abundances) predict that frequency.

  unexplained_mutations.py --variants S.variants.tsv --demix S.demix --dataset DIR|barcodes.csv --sample S
        [--gff NC_045512.gff3 --reference NC_045512.fasta] [thresholds] -o S.unexplained_mutations.tsv
  unexplained_mutations.py --concat a.tsv b.tsv ... [--min-recurrence 2] -o unexplained_mutations.tsv

Method (per sample)
  expected_freq(m) = sum over the demixed lineages l of  a_l * barcode(l, m)
where a_l is the abundance Freyja reports for lineage l (the `lineages` / `abundances` rows of the per-sample
`.demix` file, NOT the `summarized` row, which is grouped) and barcode(l, m) is 1 when the nucleotide SNV m (e.g.
C241T) is in lineage l's barcode, else 0 (a fractional value in a custom CSV is used as given, clipped to [0, 1]).
The expectation is capped at 1. Barcodes come from the dataset: `lineage_mutations.json` of a Freyja dataset
directory (identical content to `usher_barcodes.feather`), or `usher_barcodes.csv`, or a custom barcodes `.csv`
given directly (rows = lineages, columns = mutations). A directory that only has `usher_barcodes.feather` with no
`lineage_mutations.json`, or a bare `.feather` file, is an error (exit 2); pandas is deliberately not required.

Observed frequency is ivar's ALT_FREQ (ALT_DP / TOTAL_DP when ALT_FREQ is not numeric) for each variants-TSV row
whose ALT is a single base. Insertions/deletions (ALT starting with `+`/`-`) in the TSV and non-SNV tokens in the
barcodes are not considered in v1; each is counted in one stderr note. A mutation is reported ("unexplained") when

  TOTAL_DP >= --min-depth (20)   ALT_DP >= --min-alt-depth (10)   ALT_QUAL >= --min-alt-qual (20, if numeric)
  observed >= --min-freq (0.05)  observed - expected >= --min-excess (0.05)

(comparisons are inclusive; a 1e-9 tolerance keeps values that are exactly on a threshold from being lost to
floating-point noise) and, as a depth-aware guard against sampling noise,

  P(X >= ALT_DP) <= --max-pvalue (0.05 / 29903 = 1.67e-6),  X ~ Binomial(TOTAL_DP, max(expected, --background-rate (0.01)))

i.e. the alt count must be significantly more than the mixture predicts even when the expectation is allowed to be
as high as the background floor. Without it, low-depth positions (20-50x) of an ordinary mixture routinely differ
from the expectation by more than 0.05 by chance alone. The default is a Bonferroni-style family-wise 5 % over one
test per genome position (0.05 / 29903; about 30,000 positions are tested per sample, so a per-site 0.001 would still
let ~30 chance calls through). `--max-pvalue 1` switches the test off (the plain rules above). The ivar PASS
column is ignored, as elsewhere in the pipeline. Rows are sorted by excess (descending), at most --max-rows (100)
are written; the header is always written. The p-value of every reported row is the last output column
(`p_value`, 3 significant digits, e.g. 1.2e-07), also when the test is switched off. The tail is summed in
log space, so deep, clear-cut calls get their real (tiny) p-value; values below 1e-300 are written as the floor
`1e-300`, which means "at most 1e-300" (the column is always a plain float, never `0`, `<1e-300` or `NA`).

Lineages named by the demix but absent from the dataset cannot be expanded: they are counted (with the abundance
they carry) in a stderr note and contribute nothing to the expectation, which makes the table more sensitive, never
less. A sample with no usable demix (no --demix, an empty or unparseable file: skipped / failed samples) cannot have an
expectation, so nothing is called: a header-only file is written, with a stderr note. The same holds for an
unreadable variants file.

Annotation: gene, aa_change (`<ref aa><protein position><alt aa>`, e.g. D614G, `*` = stop) and effect
(missense / synonymous / stop_gained / stop_lost / start_lost / intergenic) are derived from the GFF3 CDS models and the reference
FASTA with the translator of mutations_of_interest.py (ivar's own annotation columns are not used; they are NA unless
ivar ran with --annot). A position inside several CDS (ORF1a/ORF1ab, ORF7a/ORF7b) lists the genes joined by `;`
(ORF1a is dropped when ORF1ab covers its change; at the ORF1ab slippage base two codons can be listed joined by
`/`); effect is the most severe. Without --gff/--reference these
three columns are NA.

Per-sample columns: sample, position, ref, alt, observed_freq, expected_freq, excess, depth, alt_depth, gene,
aa_change, effect, p_value, in_dataset.

`in_dataset` (last column) is `true` when at least one lineage in the dataset's barcodes (any lineage, not only
the demixed ones) carries this exact (position, alt) mutation, else `false`. `false` means the barcodes never
model this mutation (a masked / problematic site, e.g. 21987 G>A = S:G142D, which no lineage of Freyja 2.0.6 carries,
or a genuinely novel change): the table cannot tell which, so read a recurrent hit with `in_dataset = false` first as
a dataset blind spot, second as a candidate novel variant. It does not influence the call or the order.

Concat mode inserts `n_samples` (number of distinct samples in the input that have the same position + alt
unexplained: recurrence across samples is the strongest evidence that the signal is real and shared) and
`recurrent` (`true` when n_samples >= --min-recurrence, default 2). Rows are sorted recurrent first, then excess
(descending), then sample and position. Inputs must share the per-sample header (exit 2 otherwise); header-only inputs
are fine and give a header-only output. Tables from before the `in_dataset` column (13 columns) are rejected by
that header check.
"""

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

import mutations_of_interest as moi

COLUMNS = [
    "sample",
    "position",
    "ref",
    "alt",
    "observed_freq",
    "expected_freq",
    "excess",
    "depth",
    "alt_depth",
    "gene",
    "aa_change",
    "effect",
    "p_value",
    "in_dataset",
]
# concat: n_samples and recurrent go before p_value; p_value and in_dataset stay the last two columns
CONCAT_COLUMNS = [*COLUMNS[:-2], "n_samples", "recurrent", *COLUMNS[-2:]]
GENOME_LENGTH = 29903  # one test per position; only used for the default significance level
DEFAULT_MAX_PVALUE = 0.05 / GENOME_LENGTH
P_VALUE_FLOOR = 1e-300  # reported p-values are never smaller; `1e-300` reads "at most 1e-300"
LOG_NEGLIGIBLE = math.log(1e-17)  # stop the tail sum once the next term is this small vs the total
NA = "NA"
EPS = 1e-9
_SNV = re.compile(r"^([ACGT])(\d+)([ACGT])$")
_SEVERITY = {"stop_gained": 5, "start_lost": 4, "stop_lost": 3, "missense": 2, "synonymous": 1}


def warn(message):
    print(f"warning: {message}", file=sys.stderr)


# --------------------------------------------------------------------------- demix


def _tokens(text):
    """Split a demix cell: plain `a b c`, or a Python/numpy repr such as `['a', 'b']` / `[0.7 0.3]`."""
    return [t for t in re.split(r"[\s,\[\]()'\"]+", text) if t]


def read_demix(path):
    """Return [(lineage, abundance)] from a Freyja `.demix` file; ValueError if unusable."""
    rows = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\r\n").split("\t", 1)
            if len(parts) == 2 and parts[0].strip() in ("lineages", "abundances"):
                rows[parts[0].strip()] = parts[1]
    if "lineages" not in rows or "abundances" not in rows:
        raise ValueError("no `lineages` / `abundances` rows (empty or not a demix file)")
    names = _tokens(rows["lineages"])
    try:
        values = [float(t) for t in _tokens(rows["abundances"])]
    except ValueError as exc:
        raise ValueError(f"unparseable abundances: {exc}") from exc
    if len(names) != len(values):
        raise ValueError(f"{len(names)} lineages but {len(values)} abundances")
    if not names:
        raise ValueError("no lineages")
    return list(zip(names, values))


# --------------------------------------------------------------------------- barcodes


def _add(table, token, weight, counters):
    m = _SNV.match(token)
    if m is None:
        counters["non_snv"] += 1
        return
    table[(int(m.group(2)), m.group(1), m.group(3))] = weight


def _note_sites(all_sites, tokens):
    """Add the (position, alt) of every SNV token to `all_sites` (a set, or None to skip)."""
    if all_sites is None:
        return
    for token in tokens:
        m = _SNV.match(token)
        if m is not None:
            all_sites.add((int(m.group(2)), m.group(3)))


def _json_barcodes(path, wanted, counters, all_sites=None):
    with open(path) as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a JSON object lineage -> mutations")  # noqa: TRY004
    out = {}
    if all_sites is not None:
        unique = set()
        for tokens in data.values():
            unique.update(str(t).strip() for t in tokens)
        _note_sites(all_sites, unique)
    for lineage in wanted:
        if lineage not in data:
            continue
        table = {}
        for token in data[lineage]:
            _add(table, str(token).strip(), 1.0, counters)
        out[lineage] = table
    return out


def _csv_barcodes(path, wanted, counters, all_sites=None):
    wanted = set(wanted)
    out = {}
    with open(path, newline="") as handle:
        header = next(csv.reader([handle.readline()]), None)
        if not header:
            raise ValueError(f"{path}: empty barcodes file")
        columns = [c.strip() for c in header[1:]]
        carried = [False] * len(columns)
        for line in handle:
            name, _, rest = line.partition(",")
            name = name.strip().strip('"')
            if all_sites is not None:
                cells = rest.rstrip("\r\n").split(",")
                for i, cell in enumerate(cells[: len(columns)]):
                    if not carried[i] and cell.strip().strip('"') not in ("", "0", "0.0", "0.00"):
                        weight = _float(cell)
                        carried[i] = weight is not None and weight > 0
            if name not in wanted:
                continue
            values = next(csv.reader([rest]))
            if len(values) != len(columns):
                raise ValueError(
                    f"{path}: row {name!r} has {len(values)} values, header has {len(columns)}"
                )
            table = {}
            for col, val in zip(columns, values):
                try:
                    weight = min(max(float(val), 0.0), 1.0)
                except ValueError as exc:
                    raise ValueError(
                        f"{path}: row {name!r}, column {col!r}: {val!r} is not a number"
                    ) from exc
                if weight > 0:
                    _add(table, col, weight, counters)
            out[name] = table
    if all_sites is not None:
        _note_sites(all_sites, (c for c, flag in zip(columns, carried) if flag))
    return out


def load_barcodes(dataset, wanted, all_sites=None):
    """Return ({lineage: {(pos, ref, alt): weight}}, n_non_snv_tokens) for the wanted lineages.

    Lineages absent from the dataset are simply missing from the returned dict. When `all_sites` is a
    set it is filled with the (position, alt) of every SNV carried by ANY lineage of the dataset.
    """
    path = Path(dataset)
    counters = {"non_snv": 0}
    if path.is_dir():
        feather = path / "usher_barcodes.feather"
        csv_file = path / "usher_barcodes.csv"
        json_file = path / "lineage_mutations.json"
        json_ok = json_file.is_file() and json_file.stat().st_size > 0
        if feather.exists() and json_ok:
            table = _json_barcodes(json_file, wanted, counters, all_sites)
        elif feather.exists():
            raise ValueError(
                f"{dataset}: only usher_barcodes.feather found; lineage_mutations.json is required "
                "(.feather barcodes are not read here)"
            )
        elif csv_file.is_file():
            table = _csv_barcodes(csv_file, wanted, counters, all_sites)
        elif json_ok:
            table = _json_barcodes(json_file, wanted, counters, all_sites)
        else:
            raise ValueError(f"{dataset}: no lineage_mutations.json or usher_barcodes.csv")
    elif path.is_file():
        suffix = path.suffix.lower()
        if suffix == ".csv":
            table = _csv_barcodes(path, wanted, counters, all_sites)
        elif suffix == ".json":
            table = _json_barcodes(path, wanted, counters, all_sites)
        elif suffix == ".feather":
            raise ValueError(
                f"{dataset}: .feather barcodes are not supported; give the dataset directory "
                "(with lineage_mutations.json) or a .csv"
            )
        else:
            raise ValueError(f"{dataset}: unknown dataset file type (expected .csv or .json)")
    else:
        raise ValueError(f"{dataset}: dataset not found")
    return table, counters["non_snv"]


def expected_frequencies(lineages, barcodes):
    """{(pos, ref, alt): expected frequency} for a mixture [(lineage, abundance)]."""
    expected = {}
    for lineage, abundance in lineages:
        for key, weight in barcodes.get(lineage, {}).items():
            expected[key] = expected.get(key, 0.0) + abundance * weight
    return {k: min(v, 1.0) for k, v in expected.items()}


# --------------------------------------------------------------------------- variants & annotation


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_variants(path):
    """Return (rows, n_indel). Rows: dict pos, ref, alt, alt_dp, total_dp, alt_qual, freq (SNVs only, deduplicated)."""
    rows = []
    seen = set()
    n_indel = 0
    with open(path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        required = ["POS", "REF", "ALT", "ALT_DP", "ALT_FREQ", "TOTAL_DP"]
        missing = [c for c in required if c not in header]
        if missing:
            raise ValueError(f"{path}: missing columns {', '.join(missing)}")
        idx = {c: i for i, c in enumerate(header)}
        for line in handle:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            f = line.split("\t")
            f += [NA] * (len(header) - len(f))
            ref, alt = f[idx["REF"]], f[idx["ALT"]]
            if len(alt) != 1 or len(ref) != 1 or alt not in "ACGT" or ref not in "ACGT":
                n_indel += 1
                continue
            key = (f[idx["POS"]], ref, alt)
            if key in seen:
                continue
            seen.add(key)
            alt_dp = _float(f[idx["ALT_DP"]])
            total_dp = _float(f[idx["TOTAL_DP"]])
            freq = _float(f[idx["ALT_FREQ"]])
            if freq is None and alt_dp is not None and total_dp:
                freq = alt_dp / total_dp
            if alt_dp is None or total_dp is None or freq is None:
                continue
            rows.append(
                {
                    "pos": int(f[idx["POS"]]),
                    "ref": ref,
                    "alt": alt,
                    "alt_dp": int(alt_dp),
                    "total_dp": int(total_dp),
                    "alt_qual": _float(f[idx["ALT_QUAL"]]) if "ALT_QUAL" in idx else None,
                    "freq": freq,
                }
            )
    return rows, n_indel


def _effect(ref_aa, pos_aa, alt_aa):
    if ref_aa == alt_aa:
        return "synonymous"
    if pos_aa == 1 and ref_aa == "M":
        return "start_lost"
    if alt_aa == "*":
        return "stop_gained"
    if ref_aa == "*":
        return "stop_lost"
    return "missense"


def annotate(pos, ref, alt, models, reference):
    """Return (gene, aa_change, effect) strings for an SNV; NA x3 when no GFF, intergenic when in no CDS."""
    if models is None or reference is None:
        return NA, NA, NA
    row = {
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "feature": NA,
        "ref_aa": NA,
        "alt_aa": NA,
        "pos_aa": NA,
    }
    hits = []  # (gene, [(aa_change, effect), ...])
    for model in models.values():
        changes = []
        for ref_aa, pos_aa, alt_aa in moi.aa_changes_for_row(row, model["name"], model, reference):
            change = (f"{ref_aa}{pos_aa}{alt_aa}", _effect(ref_aa, pos_aa, alt_aa))
            if change not in changes:
                changes.append(change)
        if changes:
            hits.append((model["name"], changes))
    # ORF1a is the same reading frame as ORF1ab up to the slippage site: keep only ORF1ab when it covers its changes
    ab = [c for g, c in hits if g == "ORF1ab"]
    if ab:
        hits = [h for h in hits if not (h[0] == "ORF1a" and set(h[1]) <= set(ab[0]))]
    if not hits:
        return NA, NA, "intergenic"
    effects = [e for _, changes in hits for _, e in changes]
    return (
        ";".join(g for g, _ in hits),
        ";".join("/".join(c for c, _ in changes) for _, changes in hits),
        max(effects, key=lambda e: _SEVERITY[e]),
    )


# --------------------------------------------------------------------------- per-sample call


def binom_logsf(k, n, p):
    """ln P(X >= k) for X ~ Binomial(n, p), exact (stdlib only), computed in log space.

    Exact down to arbitrarily small tails (a plain sum of probabilities underflows to 0.0 below ~1e-308,
    which happens as soon as the alt count is a large share of a depth of a few hundred).
    """
    if k <= 0:
        return 0.0
    if k > n:
        return -math.inf
    p = min(max(p, 1e-12), 1 - 1e-12)
    log_p, log_q = math.log(p), math.log1p(-p)
    log_n = math.lgamma(n + 1)
    total = -math.inf
    for i in range(k, n + 1):
        term = log_n - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * log_p + (n - i) * log_q
        total = (
            term
            if total == -math.inf
            else max(total, term) + math.log1p(math.exp(-abs(total - term)))
        )
        if i > n * p and term < total + LOG_NEGLIGIBLE:
            break
    return min(total, 0.0)


def binom_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p); underflows to 0.0 below ~1e-308 (see binom_logsf)."""
    return math.exp(binom_logsf(k, n, p))


def format_pvalue(log_p):
    """p-value text for ln p: 3 significant digits, floored at P_VALUE_FLOOR (written `1e-300`).

    The floor means "at most 1e-300"; the column stays a plain float for any parser.
    """
    return f"{max(math.exp(log_p), P_VALUE_FLOOR):.3g}"


def call_unexplained(
    rows, expected, thresholds, sample, models=None, reference=None, dataset_sites=frozenset()
):
    """Return the output rows (lists of strings), sorted by excess descending, uncapped.

    `dataset_sites`: (position, alt) of every mutation any lineage of the dataset carries (`in_dataset`).
    """
    called = []
    max_log_p = math.log(thresholds["max_pvalue"]) if thresholds["max_pvalue"] > 0 else -math.inf
    for r in rows:
        if r["total_dp"] < thresholds["min_depth"] or r["alt_dp"] < thresholds["min_alt_depth"]:
            continue
        if r["alt_qual"] is not None and r["alt_qual"] < thresholds["min_alt_qual"]:
            continue
        if r["freq"] < thresholds["min_freq"] - EPS:
            continue
        exp = expected.get((r["pos"], r["ref"], r["alt"]), 0.0)
        excess = r["freq"] - exp
        if excess < thresholds["min_excess"] - EPS:
            continue
        null_rate = max(exp, thresholds["background_rate"])
        log_p = binom_logsf(r["alt_dp"], r["total_dp"], null_rate)
        if log_p > max_log_p:
            continue
        called.append((excess, r, exp, log_p))
    called.sort(
        key=lambda t: (-round(t[0], 9), t[1]["pos"], t[1]["alt"])
    )  # rounding: deterministic ties
    out = []
    for excess, r, exp, log_p in called:
        gene, aa_change, effect = annotate(r["pos"], r["ref"], r["alt"], models, reference)
        out.append(
            [
                sample,
                str(r["pos"]),
                r["ref"],
                r["alt"],
                f"{r['freq']:.6g}",
                f"{exp:.6g}",
                f"{excess:.6g}",
                str(r["total_dp"]),
                str(r["alt_dp"]),
                gene,
                aa_change,
                effect,
                format_pvalue(log_p),
                "true" if (r["pos"], r["alt"]) in dataset_sites else "false",
            ]
        )
    return out


def write_table(path, columns, rows):
    text = "\t".join(columns) + "\n" + "".join("\t".join(r) + "\n" for r in rows)
    if path == "-":
        sys.stdout.write(text)
    else:
        Path(path).write_text(text)


def run_sample(args):
    thresholds = {
        "min_depth": args.min_depth,
        "min_alt_depth": args.min_alt_depth,
        "min_freq": args.min_freq,
        "min_excess": args.min_excess,
        "min_alt_qual": args.min_alt_qual,
        "max_pvalue": args.max_pvalue,
        "background_rate": args.background_rate,
    }
    if bool(args.gff) != bool(args.reference):
        raise ValueError("--gff and --reference must be given together")
    models = moi.read_cds_models(args.gff) if args.gff else None
    reference = moi.read_fasta_sequence(args.reference) if args.reference else None

    def header_only(reason):
        warn(f"{args.sample}: {reason}; no unexplained mutations called")
        write_table(args.output, COLUMNS, [])
        return 0

    if not args.demix:
        return header_only("no demix (sample skipped or failed), no expectation possible")
    try:
        lineages = read_demix(args.demix)
    except (OSError, ValueError) as exc:
        return header_only(f"demix unusable ({exc})")
    try:
        rows, n_indel = read_variants(args.variants)
    except (OSError, ValueError) as exc:
        return header_only(f"variants unusable ({exc})")
    if not args.dataset:
        raise ValueError("--dataset is required")
    dataset_sites = set()
    barcodes, n_non_snv = load_barcodes(args.dataset, [name for name, _ in lineages], dataset_sites)
    missing = [(n, a) for n, a in lineages if n not in barcodes]
    if missing:
        mass = sum(a for _, a in missing)
        warn(
            f"{args.sample}: {len(missing)} of {len(lineages)} demixed lineages not in the dataset "
            f"({', '.join(n for n, _ in missing)}; abundance {mass:.4g}); they contribute no expected mutations"
        )
    if n_indel or n_non_snv:
        warn(
            f"{args.sample}: skipped {n_indel} indel row(s) in the variants TSV and {n_non_snv} non-SNV "
            "barcode mutation(s) (only SNVs are considered)"
        )
    expected = expected_frequencies(lineages, barcodes)
    out = call_unexplained(
        rows, expected, thresholds, args.sample, models, reference, dataset_sites
    )
    if len(out) > args.max_rows:
        warn(
            f"{args.sample}: {len(out)} unexplained mutations; reporting the top {args.max_rows} by excess"
        )
        out = out[: args.max_rows]
    write_table(args.output, COLUMNS, out)
    return 0


# --------------------------------------------------------------------------- concat


def run_concat(args):
    rows = []
    for path in args.concat:
        with open(path) as handle:
            lines = handle.read().splitlines()
        if not lines or lines[0].split("\t") != COLUMNS:
            raise ValueError(
                f"{path}: header is not the per-sample unexplained_mutations header "
                f"(expected {len(COLUMNS)} columns ending in p_value, in_dataset; tables from before "
                "the in_dataset column must be regenerated)"
            )
        for number, line in enumerate(lines[1:], start=2):
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) != len(COLUMNS):
                raise ValueError(
                    f"{path}: line {number} has {len(f)} columns, expected {len(COLUMNS)}"
                )
            rows.append(f)
    samples_by_site = {}
    for f in rows:
        samples_by_site.setdefault((f[1], f[3]), set()).add(f[0])
    out = []
    for f in rows:
        n = len(samples_by_site[(f[1], f[3])])
        out.append([*f[:-2], str(n), "true" if n >= args.min_recurrence else "false", *f[-2:]])

    def sort_key(f):
        excess = _float(f[6])
        return (f[13] != "true", -(excess if excess is not None else 0.0), f[0], int(f[1]), f[3])

    out.sort(key=sort_key)
    write_table(args.output, CONCAT_COLUMNS, out)
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--concat", nargs="+", metavar="TSV", help="merge per-sample tables (concat mode)"
    )
    p.add_argument(
        "--min-recurrence", type=int, default=2, help="concat: samples needed for `recurrent` (2)"
    )
    p.add_argument("--variants", help="ivar/Freyja variants TSV")
    p.add_argument("--demix", help="Freyja .demix file (omit for skipped/failed samples)")
    p.add_argument("--dataset", help="Freyja dataset directory, or a custom barcodes .csv")
    p.add_argument("--sample")
    p.add_argument("--gff", help="GFF3 with CDS features (gene / amino-acid annotation)")
    p.add_argument("--reference", help="reference FASTA (with --gff)")
    p.add_argument("--min-depth", type=int, default=20)
    p.add_argument("--min-alt-depth", type=int, default=10)
    p.add_argument("--min-freq", type=float, default=0.05)
    p.add_argument("--min-excess", type=float, default=0.05)
    p.add_argument("--min-alt-qual", type=float, default=20)
    p.add_argument(
        "--max-pvalue",
        type=float,
        default=DEFAULT_MAX_PVALUE,
        help="binomial test; 1 disables (0.05/29903)",
    )
    p.add_argument(
        "--background-rate", type=float, default=0.01, help="floor for the null rate (0.01)"
    )
    p.add_argument(
        "--max-rows",
        type=int,
        default=100,
        help="write at most this many rows per sample, the highest excess first (ties: position, "
        "then alt); 0 writes the header only; a stderr warning says how many were cut (100)",
    )
    p.add_argument("-o", "--output", default="-")
    args = p.parse_args(argv)
    if args.max_rows < 0:
        p.error("--max-rows must be >= 0")

    try:
        if args.concat:
            return run_concat(args)
        if not args.variants or not args.sample:
            p.error("--variants and --sample are required (or use --concat)")
        return run_sample(args)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
