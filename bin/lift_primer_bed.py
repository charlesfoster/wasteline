#!/usr/bin/env python3
"""Lift a primer scheme BED from the scheme's reference onto another reference (maintenance tool, PLAN D27).

Not part of the pipeline: run it by hand when a published scheme was designed on a different strain of the pathogen than the reference
Freyja's barcodes (and so wasteline) use. The pipeline never remaps coordinates. Needs `minimap2` on PATH (or --minimap2), e.g.
`CONDA_SUBDIR=osx-64 mamba create -p liftenv -c conda-forge -c bioconda minimap2`.

Method
  1. minimap2 (`-x asm20 -c`, base-level CIGAR) aligns the scheme's reference contig (query) to our reference (target). Every query base
     that sits in an aligned (match/mismatch) column gets a target coordinate. If several alignments cover a base (repeats such as the
     mpox terminal repeats), the one whose diagonal is closest to that of the longest alignment wins, so a base is never lifted to the
     other copy of a repeat.
  2. Each primer's start is the first and its end the last mapped base *inside* the primer (BED half-open, 0-based). A primer with no
     mapped base is "in a gap". Indels inside the primer show up as a length change.
  3. The primer sequence (BED column 7, else the scheme reference at the original coordinates, reverse-complemented for '-') is compared
     with OUR reference at the lifted position (reverse-complemented for RIGHT primers), anchored at the primer's 3' end, IUPAC-aware.
     Identity = matching bases / compared bases.

Safe behaviour (a primer is "bad" if it is in a gap, changes length by more than --max-length-change (2) bp, or has identity below
--min-identity (0.8)): the bad primer is DROPPED from the output and reported. An amplicon (primer name minus `_LEFT_n`/`_RIGHT_n`)
keeps its other alternates, but if one side (LEFT or RIGHT) is left without any primer the WHOLE amplicon is dropped, so the output
never has an unpaired LEFT or RIGHT (a left-only and a right-only name with the same trailing number count as one
amplicon, as in amplicon_coverage.py). Primers whose name cannot be paired are kept if good. Nothing is silently lost: every dropped
primer and amplicon is in the report. A dropped amplicon is a hole in the lifted scheme (a coverage gap to expect in the results).

Limits: a base-level alignment of two strains, not a model of the primers' binding. Lifted coordinates are only as good as the
alignment around indels and the genome ends (minimap2 clips the first and last few dozen bases of divergent ends, so primers there
are reported as gaps). Identity is measured against the single reference sequence; alternate primers for other strains legitimately
differ from it. Check the report, and the amplicon sizes, before using a lifted BED; validate it on reads.

Usage:
  lift_primer_bed.py --bed primer.bed --scheme-ref reference.fasta --scheme-contig NC_038235.1 --target-ref ours.fasta \\
      --out lifted.bed --report lifted.report.tsv   # also writes lifted.report.amplicons.tsv
Output BED: contig = our reference's sequence name, columns contig start end name pool strand [sequence], plus `#` header lines.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}
NAME = re.compile(r"^(.+?)_(LEFT|RIGHT)(?:_.*)?$")


def read_fasta(path):
    seqs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def rc(s):
    return s.upper().translate(COMP)[::-1]


def position_map(query_fa, target_fa, minimap2):
    """query position -> target position (None if not in an aligned column), from minimap2 CIGARs."""
    out = subprocess.run(
        [minimap2, "-x", "asm20", "-c", "-t", "1", str(target_fa), str(query_fa)],
        capture_output=True,
        text=True,
        check=False,
    )
    if out.returncode:
        sys.exit(f"minimap2 failed: {out.stderr[-500:]}")
    recs = []
    for line in out.stdout.splitlines():
        f = line.split("\t")
        cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
        recs.append(
            (int(f[9]), f[4], int(f[2]), int(f[7]), cg)
        )  # matches, strand, query start, target start, cigar
    if not recs:
        sys.exit("minimap2 found no alignment between the two references")
    if any(r[1] == "-" for r in recs):
        sys.stderr.write("WARNING: some alignments are on the reverse strand and are ignored\n")
    recs = [r for r in recs if r[1] == "+"]
    main = max(recs)
    main_diag = main[3] - main[2]
    best, cost = {}, {}
    for _, _, q, t, cg in recs:
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cg):
            n = int(n)
            if op in "M=X":
                c = abs(t - q - main_diag)
                for i in range(n):
                    if q + i not in best or c < cost[q + i]:
                        best[q + i], cost[q + i] = t + i, c
                q, t = q + n, t + n
            elif op == "I":
                q += n
            elif op in "DN":
                t += n
    return best


def identity(primer, ref_slice, strand_minus):
    """Fraction of bases of `primer` (5'->3') that match the reference slice, anchored at the 3' end, IUPAC-aware."""
    target = rc(ref_slice) if strand_minus else ref_slice  # the slice as the primer would read it
    n = min(len(primer), len(target))
    if n == 0:
        return 0.0
    p, t = primer[-n:], target[-n:]
    return sum(set(IUPAC.get(a, "N")) & set(IUPAC.get(b, "N")) != set() for a, b in zip(p, t)) / n


def tsv(fh, cols):
    """One tab-separated line; None becomes NA."""
    fh.write("\t".join("NA" if c is None else str(c) for c in cols) + "\n")


def read_primers(a, qseq, tseq, pmap):
    """The scheme's primers on the contig, lifted: ns/ne = new start/end, dlen = length change, ident, why = reasons to drop."""
    primers = []
    for line in Path(a.bed).read_text().splitlines():
        f = [x.strip() for x in line.split("\t")]
        if not line.strip() or line.startswith("#") or f[0] != a.scheme_contig:
            continue
        s, e, strand = int(f[1]), int(f[2]), f[5]
        have_seq = len(f) > 6 and re.fullmatch(r"[A-Za-z]+", f[6])
        seq = f[6].upper() if have_seq else (rc(qseq[s:e]) if strand == "-" else qseq[s:e].upper())
        m = NAME.match(f[3])
        p = {"name": f[3], "pool": f[4], "strand": strand, "seq": seq, "s": s, "e": e, "why": []}
        p["amp"], p["side"] = (m.group(1), m.group(2)) if m else (None, None)
        p["ns"] = p["ne"] = p["dlen"] = p["ident"] = None
        mapped = [pmap[i] for i in range(s, e) if i in pmap]
        if not mapped:
            p["why"].append("in alignment gap")
        else:
            p["ns"], p["ne"] = min(mapped), max(mapped) + 1
            p["dlen"] = (p["ne"] - p["ns"]) - (e - s)
            p["ident"] = identity(seq, tseq[p["ns"] : p["ne"]], strand == "-")
            if abs(p["dlen"]) > a.max_length_change:
                p["why"].append(f"length change {p['dlen']:+d}")
            if p["ident"] < a.min_identity:
                p["why"].append(f"identity {p['ident']:.0%}")
        primers.append(p)
    return primers


def drop_unpaired(primers):
    """Mark the good primers of amplicons with no good LEFT or no good RIGHT primer left; return ({amplicon: primers}, dropped)."""
    amps = {}
    for p in primers:
        if p["amp"]:
            amps.setdefault(p["amp"], []).append(p)
    # as in amplicon_coverage.py: one left-only and one right-only name with the same trailing number are one amplicon
    # (artic-inrb-mpox amplicon 1 is `2o0fvmwf_1_LEFT_*` + `359ba5b8_1_RIGHT_*`)
    number = lambda k: re.search(r"(\d+)_*$", k)
    only = {
        s: [k for k, ps in amps.items() if {p["side"] for p in ps} == {s}]
        for s in ("LEFT", "RIGHT")
    }
    for k in only["LEFT"]:
        mates = [
            r
            for r in only["RIGHT"]
            if r in amps and number(r) and number(k) and number(r)[1] == number(k)[1]
        ]
        if len(mates) == 1:
            amps[k] += amps.pop(mates[0])
            for p in amps[k]:
                p["amp"] = k
    good = {}
    for k, ps in amps.items():
        good[k] = {p["side"] for p in ps if not p["why"]}
    dropped = {k for k in amps if good.get(k, set()) != {"LEFT", "RIGHT"}}
    for p in primers:
        if p["amp"] in dropped and not p["why"]:
            p["why"].append("amplicon dropped (no good primer on one side)")
    return amps, dropped


def span(ps, new):
    """Amplicon span: leftmost LEFT start to rightmost RIGHT end (None if one side is missing)."""
    ls = [p["ns"] if new else p["s"] for p in ps if p["side"] == "LEFT"]
    rs = [p["ne"] if new else p["e"] for p in ps if p["side"] == "RIGHT"]
    return (min(ls), max(rs)) if ls and rs else None


def size(x):
    return x[1] - x[0] if x else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bed", required=True, help="the scheme's primer.bed")
    ap.add_argument("--scheme-ref", required=True, help="the scheme's reference FASTA")
    ap.add_argument("--scheme-contig", required=True, help="contig of the BED/FASTA to lift")
    ap.add_argument("--target-ref", required=True, help="OUR reference FASTA (first record)")
    ap.add_argument("--out", required=True, help="lifted BED")
    ap.add_argument("--report", required=True, help="per-primer TSV; <report>.amplicons.tsv too")
    ap.add_argument("--comment", action="append", default=[], help="extra '#' BED header line")
    ap.add_argument("--min-identity", type=float, default=0.8)
    ap.add_argument("--max-length-change", type=int, default=2)
    ap.add_argument("--minimap2", default="minimap2")
    a = ap.parse_args()

    sref, tref = read_fasta(a.scheme_ref), read_fasta(a.target_ref)
    if a.scheme_contig not in sref:
        sys.exit(f"{a.scheme_contig} not in {a.scheme_ref} (has {', '.join(sref)})")
    tname = next(iter(tref))
    tseq, qseq = tref[tname], sref[a.scheme_contig]
    q, t = Path(a.out).with_suffix(".query.fa"), Path(a.out).with_suffix(".target.fa")
    q.write_text(f">{a.scheme_contig}\n{qseq}\n")
    t.write_text(f">{tname}\n{tseq}\n")
    try:
        pmap = position_map(q, t, a.minimap2)
    finally:
        q.unlink(missing_ok=True)
        t.unlink(missing_ok=True)

    primers = read_primers(a, qseq, tseq, pmap)
    amps, dropped_amps = drop_unpaired(primers)
    kept = [p for p in primers if not p["why"]]

    with open(a.out, "w") as fh:
        fh.write(
            f"# lifted by bin/lift_primer_bed.py from {a.scheme_contig} ({Path(a.bed).name}) onto {tname}\n"
        )
        fh.write(
            f"# kept {len(kept)} of {len(primers)} primers (dropped: gap, length change > {a.max_length_change} bp, "
        )
        fh.write(f"identity < {a.min_identity:.0%}, unpaired amplicon)\n")
        fh.writelines(f"# {c}\n" for c in a.comment)
        for p in sorted(kept, key=lambda p: (p["ns"], p["ne"], p["name"])):
            tsv(fh, [tname, p["ns"], p["ne"], p["name"], p["pool"], p["strand"], p["seq"]])
    with open(a.report, "w") as fh:
        tsv(
            fh,
            [
                "primer",
                "amplicon",
                "strand",
                "pool",
                "orig_start",
                "orig_end",
                "new_start",
                "new_end",
                "length_change",
                "identity",
                "status",
                "reason",
            ],
        )
        for p in primers:
            ident = None if p["ident"] is None else f"{p['ident']:.3f}"
            row = [
                p["name"],
                p["amp"],
                p["strand"],
                p["pool"],
                p["s"],
                p["e"],
                p["ns"],
                p["ne"],
                p["dlen"],
                ident,
            ]
            tsv(fh, row + ["dropped" if p["why"] else "kept", "; ".join(p["why"]) or "-"])
    with open(a.report.replace(".tsv", "") + ".amplicons.tsv", "w") as fh:
        tsv(
            fh,
            [
                "amplicon",
                "primers",
                "kept_primers",
                "orig_start",
                "orig_end",
                "orig_size",
                "new_start",
                "new_end",
                "new_size",
                "size_change",
                "status",
            ],
        )
        for k, ps in amps.items():
            o = span(ps, False)  # None for an amplicon that is unpaired in the scheme itself
            n = None if k in dropped_amps else span([p for p in ps if not p["why"]], True)
            change = size(n) - size(o) if n and o else None
            row = [
                k,
                len(ps),
                sum(not p["why"] for p in ps),
                o and o[0],
                o and o[1],
                size(o),
                n and n[0],
                n and n[1],
                size(n),
            ]
            tsv(fh, row + [change, "kept" if n else "dropped"])
    idents = [p["ident"] for p in primers if p["ident"] is not None]
    bins = ", ".join(f">={b}%: {sum(i >= b / 100 for i in idents)}" for b in (95, 90, 80))
    sys.stderr.write(
        f"{a.scheme_contig} -> {tname}: {len(primers)} primers, {len(kept)} kept, {len(primers) - len(kept)} dropped; "
        f"{len(amps) - len(dropped_amps)} of {len(amps)} amplicons kept\n"
        f"primer identity at the lifted position (all primers): {bins} of {len(idents)}\n"
    )


if __name__ == "__main__":
    main()
