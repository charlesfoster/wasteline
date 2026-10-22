#!/usr/bin/env python3
"""Frequency and depth of mutations of interest in one sample's ivar/Freyja variants TSV.

  mutations_of_interest.py --variants S.variants.tsv --moi moi.csv --sample S
        [--depths S.depths] [--gff NC_045512.gff3 --reference NC_045512.fasta] -o S.moi.tsv

moi.csv has columns `name,mutation`. `mutation` is one of

  A23063T            nucleotide SNV (ref base, 1-based position, alt base)
  C21995-            deletion: the deleted bases and the 1-based position of the FIRST deleted base
                     (ivar reports it as ALT `-C` at the preceding position; the conversion is done here)
  CTT21991-          multi-base deletion, same convention
  T22204+GAGCCAGAA   insertion after the 1-based position 22204 (ivar ALT `+GAGCCAGAA` at that position)
  S:N501Y            amino-acid substitution: gene (CDS/transcript name, e.g. S, N, ORF1a, ORF1ab, ORF3a),
                     reference residue, protein position, alt residue (`*` = stop)

Amino-acid mutations are resolved from the variants TSV's GFF_FEATURE/REF_AA/POS_AA/ALT_AA columns when the
TSV was annotated (ivar --annot) and the feature equals the gene; otherwise from the nucleotide SNV, the
GFF3 CDS model and the reference FASTA with a minimal codon translator (+ strand CDS, spliced/slippage
CDS such as ORF1ab included). Each SNV row is translated on its own: several SNVs inside one codon
are not combined. If several nucleotide rows give the same amino-acid change their alt depths and
frequencies are summed. All ivar rows are considered regardless of the PASS column.

Output columns: sample, name, mutation, ref_depth, alt_depth, depth, alt_freq (ref_depth/alt_depth/depth =
REF_DP/ALT_DP/TOTAL_DP of the matching row). A mutation absent from the variants file is reported with
alt_depth 0 and alt_freq 0; with --depths its depth is the depth at the position (amino acid: the minimum over
the three codon positions) and ref_depth equals depth, otherwise both are NA. If the variants file is missing or
unreadable every row is NA (with a warning). Unparseable mutations or unresolvable amino-acid mutations are
errors (exit 2), as is an invalid MOI file.

Multi-contig references: `--contig C` selects the contig of the reference FASTA (first header word), the GFF3
(seqid), the depths file (column 1) and, when its REGION column holds several contigs, the variants TSV. An
input with several contigs and no `--contig` (or a `--contig` that is not in the FASTA/GFF3/depths) is an error
(exit 2), never the first record. Without `--contig` single-contig inputs behave as before.
"""

import argparse
import contextlib
import csv
import re
import sys

COLUMNS = ["sample", "name", "mutation", "ref_depth", "alt_depth", "depth", "alt_freq"]
NA = "NA"

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AAS[i]
    for i, (a, b, c) in enumerate((a, b, c) for a in _BASES for b in _BASES for c in _BASES)
}

_NT_SNV = re.compile(r"^([ACGT])(\d+)([ACGT])$")
_NT_DEL = re.compile(r"^([ACGT]+)(\d+)-$")
_NT_INS = re.compile(r"^([ACGT])(\d+)\+([ACGT]+)$")
_AA = re.compile(r"^([A-Za-z0-9_.]+):([A-Z*])(\d+)([A-Z*])$")


def translate(codon):
    return CODON_TABLE.get(codon.upper(), "X")


def parse_mutation(text):
    """Return a dict describing the mutation; ValueError if not understood."""
    text = text.strip()
    m = _AA.match(text)
    if m:
        return {
            "kind": "aa",
            "gene": m.group(1),
            "ref": m.group(2),
            "pos": int(m.group(3)),
            "alt": m.group(4),
        }
    m = _NT_SNV.match(text)
    if m:
        return {"kind": "snv", "ref": m.group(1), "pos": int(m.group(2)), "alt": m.group(3)}
    m = _NT_DEL.match(text)
    if m:
        return {"kind": "del", "bases": m.group(1), "pos": int(m.group(2))}
    m = _NT_INS.match(text)
    if m:
        return {"kind": "ins", "ref": m.group(1), "pos": int(m.group(2)), "bases": m.group(3)}
    raise ValueError(
        f"cannot parse mutation {text!r} (expected e.g. A23063T, C21995-, T22204+GAGCCAGAA or S:N501Y)"
    )


def read_moi(path):
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if (
            not reader.fieldnames
            or "name" not in reader.fieldnames
            or "mutation" not in reader.fieldnames
        ):
            raise ValueError(f"{path}: expected a header with columns name,mutation")
        for line_number, rec in enumerate(reader, start=2):
            name = (rec.get("name") or "").strip()
            mutation = (rec.get("mutation") or "").strip()
            if not name and not mutation:
                continue
            try:
                parsed = parse_mutation(mutation)
            except ValueError as exc:
                raise ValueError(f"{path}: line {line_number}: {exc}") from exc
            rows.append((name or mutation, mutation, parsed))
    return rows


def region_filter(path, header, lines, contig):
    """Return the REGION value to keep, or None (keep all rows).

    Only a variants table whose REGION column holds several contigs is filtered: it needs `contig`
    (ContigError otherwise). A single-region table is used as is, whatever `contig` is.
    """
    if "REGION" not in header:
        return None
    i = header.index("REGION")
    regions = {line.split("\t")[i] for line in lines if len(line.split("\t")) > i}
    if len(regions) <= 1:
        return None
    if contig is None:
        raise ContigError(f"{path}: {len(regions)} contigs in the REGION column, pass --contig")
    return contig


def read_variants(path, contig=None):
    """Return (rows, annotated). Rows are dicts with the columns used here; duplicates are dropped.

    A variants table with several contigs in REGION needs `contig` and keeps that contig's rows only.
    """
    rows = []
    seen = set()
    annotated = False
    with open(path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        required = ["POS", "REF", "ALT", "REF_DP", "ALT_DP", "ALT_FREQ", "TOTAL_DP"]
        missing = [c for c in required if c not in header]
        if missing:
            raise ValueError(f"{path}: missing columns {', '.join(missing)}")
        idx = {c: i for i, c in enumerate(header)}
        lines = [line.rstrip("\r\n") for line in handle]
        lines = [line for line in lines if line.strip()]
        keep_region = region_filter(path, header, lines, contig)
        for line in lines:
            f = line.split("\t")
            if len(f) < len(header):
                f += [NA] * (len(header) - len(f))
            if keep_region is not None and f[idx["REGION"]] != keep_region:
                continue

            def get(col, f=f):
                return f[idx[col]] if col in idx else NA

            key = (get("POS"), get("REF"), get("ALT"))
            if key in seen:
                continue
            seen.add(key)
            row = {
                "pos": int(get("POS")),
                "ref": get("REF"),
                "alt": get("ALT"),
                "ref_dp": get("REF_DP"),
                "alt_dp": get("ALT_DP"),
                "total_dp": get("TOTAL_DP"),
                "alt_freq": get("ALT_FREQ"),
                "feature": get("GFF_FEATURE"),
                "ref_aa": get("REF_AA"),
                "alt_aa": get("ALT_AA"),
                "pos_aa": get("POS_AA"),
            }
            if row["feature"] not in (NA, "", "."):
                annotated = True
            rows.append(row)
    return rows, annotated


class ContigError(ValueError):
    """Several contigs in an input and no (or an unknown) --contig: a hard error, never softened to a warning."""


def read_depths(path, contig=None):
    """Return {position: depth}. A depths file with several contigs needs `contig` (its rows only)."""
    depths = {}
    contigs = set()
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 4:
                raise ValueError(
                    f"{path}: line {line_number} has {len(f)} columns, expected at least 4"
                )
            contigs.add(f[0])
            if contig is not None and f[0] != contig:
                continue
            depths[int(f[1])] = int(float(f[3]))
    if contig is None and len(contigs) > 1:
        raise ContigError(f"{path}: {len(contigs)} contigs in the depths file, pass --contig")
    if contig is not None and contigs and contig not in contigs:
        raise ContigError(
            f"{path}: contig {contig!r} not in the depths file (contigs: {', '.join(sorted(contigs))})"
        )
    return depths


def read_fasta_sequence(path, contig=None):
    """Upper-case sequence of the FASTA record named `contig` (first word of the header).

    Without `contig` the file must hold a single record (ContigError otherwise).
    """
    records = {}
    names = []
    current = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                words = line[1:].split()
                current = []
                names.append(words[0] if words else "")
                records.setdefault(names[-1], current)
            elif line:
                if current is None:  # a header-less file is read as one record
                    current = []
                    names.append("")
                    records[""] = current
                current.append(line)
    if contig is None:
        if len(names) > 1:
            raise ContigError(f"{path}: {len(names)} FASTA records, pass --contig")
        seq = records[names[0]] if names else []
    else:
        if contig not in records:
            raise ContigError(
                f"{path}: contig {contig!r} not in the FASTA (records: {', '.join(names) or 'none'})"
            )
        seq = records[contig]
    if not seq:
        raise ValueError(f"{path}: no sequence")
    return "".join(seq).upper()


def _attrs(text):
    out = {}
    for part in text.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def read_cds_models(path, contig=None):
    """Return dict lower-case name -> {'name', 'positions': [1-based genomic position per CDS base]}.

    The feature name is the Name of the parent transcript (so ORF1a and ORF1ab are distinct although they
    share a gene Name), falling back to the gene Name, then the CDS's own Name/gene attribute.

    Only features on the GFF seqid `contig` are used. Without `contig` the GFF must hold a single seqid
    (ContigError otherwise); with it the seqid must occur in the file.
    """
    features = {}
    cds = {}
    order = []
    seqids = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 9:
                continue
            if f[0] not in seqids:
                seqids.append(f[0])
            if contig is not None and f[0] != contig:
                continue
            attrs = _attrs(f[8])
            if f[2] in ("gene", "mRNA", "transcript") and "ID" in attrs:
                features[attrs["ID"]] = (attrs.get("Name"), attrs.get("Parent"))
            elif f[2] == "CDS":
                parent = attrs.get("Parent") or attrs.get("ID") or "CDS"
                if parent not in cds:
                    cds[parent] = []
                    order.append(parent)
                cds[parent].append((int(f[3]), int(f[4]), f[6], attrs))
    if contig is None and len(seqids) > 1:
        raise ContigError(f"{path}: {len(seqids)} seqids in the GFF3, pass --contig")
    if contig is not None and seqids and contig not in seqids:
        raise ContigError(
            f"{path}: contig {contig!r} not a seqid of the GFF3 (seqids: {', '.join(seqids)})"
        )
    models = {}
    for parent in order:
        parts = sorted(cds[parent])
        if any(strand != "+" for _, _, strand, _ in parts):
            continue  # minus-strand CDS are not supported (none in SARS-CoV-2)
        name = None
        node = parent
        for _ in range(3):
            if node in features:
                name = features[node][0]
                if name:
                    break
                node = features[node][1]
        if not name:
            attrs = parts[0][3]
            name = attrs.get("Name") or attrs.get("gene") or parent
        positions = []
        for start, end, _, _ in parts:
            positions.extend(range(start, end + 1))
        models.setdefault(name.lower(), {"name": name, "positions": positions})
    return models


def codon_positions(model, protein_pos):
    start = 3 * (protein_pos - 1)
    cp = model["positions"][start : start + 3]
    if len(cp) < 3:
        raise ValueError(f"{model['name']}: protein position {protein_pos} is beyond the CDS")
    return cp


def _is_snv(row):
    return (
        len(row["ref"]) == 1
        and len(row["alt"]) == 1
        and row["alt"] in "ACGT"
        and row["ref"] in "ACGT"
    )


def aa_changes_for_row(row, gene, model, reference):
    """List of (ref_aa, protein_pos, alt_aa) for an SNV row in `gene` (annotation first, else computed)."""
    if not _is_snv(row):
        return []
    if (
        row["feature"].lower() == gene.lower()
        and row["alt_aa"] not in (NA, "", ".")
        and row["pos_aa"] not in (NA, "", ".")
    ):
        return [(row["ref_aa"], int(row["pos_aa"]), row["alt_aa"])]
    if model is None or reference is None:
        return []
    out = []
    positions = model["positions"]
    for i, pos in enumerate(positions):
        if pos != row["pos"]:
            continue
        c0 = (i // 3) * 3
        cps = positions[c0 : c0 + 3]
        if len(cps) < 3:
            continue
        ref_codon = "".join(reference[p - 1] for p in cps)
        alt_codon = "".join(row["alt"] if p == row["pos"] else reference[p - 1] for p in cps)
        out.append((translate(ref_codon), c0 // 3 + 1, translate(alt_codon)))
    return out


def _num(value):
    try:
        return float(value)
    except ValueError:
        return None


def _fmt_count(value):
    num = _num(value)
    return NA if num is None else str(int(num))


def _fmt_freq(value):
    num = _num(value)
    return NA if num is None else f"{num:.6g}"


def _row_result(row):
    return [
        _fmt_count(row["ref_dp"]),
        _fmt_count(row["alt_dp"]),
        _fmt_count(row["total_dp"]),
        _fmt_freq(row["alt_freq"]),
    ]


def _absent_result(depth):
    if depth is None:
        return [NA, "0", NA, "0"]
    return [str(depth), "0", str(depth), "0"]


def evaluate(parsed, rows, annotated, models, reference, depths):
    """Return [ref_depth, alt_depth, depth, alt_freq] strings for one parsed mutation."""
    kind = parsed["kind"]
    if kind == "snv":
        for r in rows:
            if (
                r["pos"] == parsed["pos"]
                and r["alt"] == parsed["alt"]
                and r["ref"] == parsed["ref"]
            ):
                return _row_result(r)
        return _absent_result(depths.get(parsed["pos"]) if depths is not None else None)
    if kind == "del":
        anchor = parsed["pos"] - 1
        for r in rows:
            if r["pos"] == anchor and r["alt"] == "-" + parsed["bases"]:
                return _row_result(r)
        return _absent_result(depths.get(parsed["pos"]) if depths is not None else None)
    if kind == "ins":
        for r in rows:
            if (
                r["pos"] == parsed["pos"]
                and r["alt"] == "+" + parsed["bases"]
                and r["ref"] == parsed["ref"]
            ):
                return _row_result(r)
        return _absent_result(depths.get(parsed["pos"]) if depths is not None else None)

    gene = parsed["gene"]
    model = models.get(gene.lower()) if models is not None else None
    if model is None and not annotated:
        if models is None:
            raise ValueError(
                f"{gene}:{parsed['ref']}{parsed['pos']}{parsed['alt']}: needs --gff and --reference "
                "(the variants TSV carries no GFF annotation)"
            )
        raise ValueError(
            f"gene {gene!r} not found in the GFF3 (known: {', '.join(sorted(m['name'] for m in models.values()))})"
        )
    cps = None
    if model is not None:
        cps = codon_positions(model, parsed["pos"])
        if reference is not None:
            ref_aa = translate("".join(reference[p - 1] for p in cps))
            if ref_aa != parsed["ref"]:
                print(
                    f"warning: {gene}:{parsed['ref']}{parsed['pos']}{parsed['alt']}: reference residue is {ref_aa}",
                    file=sys.stderr,
                )
    hits = []
    for r in rows:
        for ref_aa, pos_aa, alt_aa in aa_changes_for_row(r, gene, model, reference):
            if pos_aa == parsed["pos"] and alt_aa == parsed["alt"]:
                hits.append(r)
                break
    if hits:
        first = _row_result(hits[0])
        if len(hits) == 1:
            return first
        alt_dp = sum(_num(h["alt_dp"]) or 0 for h in hits)
        freq = sum(_num(h["alt_freq"]) or 0 for h in hits)
        return [first[0], str(int(alt_dp)), first[2], f"{min(freq, 1.0):.6g}"]
    depth = None
    if depths is not None and cps is not None:
        depth = min(depths.get(p, 0) for p in cps)
    return _absent_result(depth)


def write_rows(path, rows):
    handle = sys.stdout if path == "-" else open(path, "w")  # noqa: SIM115
    with contextlib.nullcontext(handle) if path == "-" else handle as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--variants", required=True, help="ivar/Freyja variants TSV")
    p.add_argument("--moi", required=True, help="CSV with columns name,mutation")
    p.add_argument("--sample", required=True)
    p.add_argument("--depths", help="Freyja depths file (depth of absent mutations)")
    p.add_argument("--gff", help="GFF3 with CDS features (amino-acid fallback)")
    p.add_argument("--reference", help="reference FASTA (amino-acid fallback)")
    p.add_argument(
        "--contig",
        help="multi-contig inputs: use this contig of the reference FASTA, GFF3 seqid, depths and "
        "variants REGION (required when any of them holds several contigs; an error otherwise)",
    )
    p.add_argument("-o", "--output", default="-")
    args = p.parse_args(argv)

    try:
        moi = read_moi(args.moi)
        models = read_cds_models(args.gff, args.contig) if args.gff else None
        reference = read_fasta_sequence(args.reference, args.contig) if args.reference else None
        if models is not None and reference is None:
            raise ValueError("--gff requires --reference")
        depths = None
        if args.depths:
            try:
                depths = read_depths(args.depths, args.contig)
            except ContigError:
                raise
            except (OSError, ValueError) as exc:
                print(
                    f"warning: {args.sample}: depths unusable ({exc}); depth reported as NA",
                    file=sys.stderr,
                )
        try:
            rows, annotated = read_variants(args.variants, args.contig)
        except ContigError:
            raise
        except (OSError, ValueError) as exc:
            print(
                f"warning: {args.sample}: variants unusable ({exc}); all values NA", file=sys.stderr
            )
            write_rows(args.output, [[args.sample, n, m, NA, NA, NA, NA] for n, m, _ in moi])
            return 0
        out = []
        for name, mutation, parsed in moi:
            values = evaluate(parsed, rows, annotated, models, reference, depths)
            out.append([args.sample, name, mutation, *values])
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    write_rows(args.output, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
