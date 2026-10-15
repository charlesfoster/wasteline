"""Deterministic synthetic wastewater dataset generator for wasteline (PLAN section 8).

Pure Python standard library only (no numpy / badread / art).

What it does
------------
1. (optional, one-off) ``--extract-barcodes SRC.csv`` : cut a small Freyja ``--barcodes`` .csv out of the
   full ``usher_barcodes.csv`` bundled in the Freyja container (see docs/test_data.md).
2. Builds "lineage genomes" by applying a lineage's barcode SNVs to NC_045512.2.
3. Simulates 2x150 Illumina-like amplicon reads from a primer BED (default ``midnight``), with per-amplicon depth
   jitter, optional amplicon dropout, a small error rate, and primer sequence at read starts.
4. Writes ``fastq/*.fastq.gz``, ``samplesheet.csv``, ``truth.json``, ``barcodes_test.csv``,
   ``primer_scheme.txt`` into ``--out``.

Usage
-----
    python tests/make_test_data.py --out tests/data --seed 42
    python tests/make_test_data.py --protocol capture              # -> tests/data_capture (probe-capture data)

``--protocol capture`` (PLAN D13) writes a SEPARATE dataset (tests/data_capture/) that mimics probe-capture
libraries (e.g. Twist): fragment starts uniform over the genome (reads may overhang the genome ends), no primer
sequences at read starts, mild depth variation, and PCR duplicates (exact copies of ~17% of the read pairs). The
amplicon dataset in tests/data is never touched by it.

Output is byte-identical for the same seed, reference, primer BED and barcode file (gzip with mtime=0; only
``random.Random.random()`` and ``gauss`` are used so the stream does not depend on Python-version-specific helpers).
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import random
import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = REPO / "assets" / "references" / "NC_045512.fasta"
DEFAULT_PRIMER_DIR = REPO / "assets" / "primer_schemes"
DEFAULT_BARCODES = REPO / "tests" / "data" / "barcodes_test.csv"

# Real lineage names present in Freyja's lineages.yml. "B" is the Wuhan-like reference (no barcode mutations).
# Chosen to be well separated (early B.1.x, VOCs, Delta AY.x, Omicron BA.x, XBB, JN.1/KP.x).
BARCODE_LINEAGES = [
    "B",
    "B.1",
    "B.1.1.7",
    "B.1.351",
    "P.1",
    "AY.4",
    "BA.1",
    "BA.2",
    "BA.5",
    "BQ.1.1",
    "XBB.1.5",
    "EG.5.1",
    "JN.1",
    "KP.3",
]

READ_LEN = 150
ERROR_RATE = 0.002
FRAG_MEAN, FRAG_SD, FRAG_MIN, FRAG_MAX = 300, 40, 200, 450
# Fraction of fragments anchored at the amplicon 5' end / 3' end (these start with the primer; the remainder start
# at a uniform random position inside the amplicon, like a tagmentation library).
P_ANCHOR_5, P_ANCHOR_3 = 0.20, 0.20
DEPTH_JITTER_SD = 0.35  # sd of log-normal per-amplicon weight
LOW_COVERAGE_MAX_PAIRS = 1000  # samples below this are expected to fail the coverage gate
QC_MAX_DROPPED_DEFAULT = (
    1  # pipeline default --qc_max_dropped_amplicons: more dropouts than this is a QC reason (WARN)
)

# Sample definitions. ``dropout`` = amplicon names (from the BED) that receive zero reads. ``pct_genome_ge_10x`` is
# the expected [min, max] percent of the genome covered at >= 10x (measured with bwa mem + freyja variants).
SAMPLES = [
    {
        "id": "MIX_A",
        "pct_genome_ge_10x": [95, 100],
        "sample_type": "sample",
        "date": "2026-01-05",
        "site": "SiteA",
        "mix": {"BA.2": 0.70, "AY.4": 0.30},
        "pairs": 6000,
    },
    {
        "id": "MIX_B",
        "pct_genome_ge_10x": [95, 100],
        "sample_type": "sample",
        "date": "2026-01-12",
        "site": "SiteA",
        "mix": {"BA.1": 0.50, "XBB.1.5": 0.30, "KP.3": 0.20},
        "pairs": 7000,
    },
    {
        "id": "PURE_C",
        "pct_genome_ge_10x": [95, 100],
        "sample_type": "sample",
        "date": "2026-01-19",
        "site": "SiteB",
        "mix": {"JN.1": 1.0},
        "pairs": 8000,
    },
    {
        "id": "LOWCOV_D",
        "pct_genome_ge_10x": [0, 5],
        "sample_type": "sample",
        "date": "2026-01-26",
        "site": "SiteB",
        "mix": {"BA.2": 1.0},
        "pairs": 60,
    },
    {
        "id": "NEG_E",
        "pct_genome_ge_10x": [0, 5],
        "sample_type": "negative",
        "date": "2026-01-05",
        "site": "SiteA",
        "mix": {"B": 1.0},
        "pairs": 3,
    },
    {
        # contaminated negative control (PLAN D8): deep enough to PASS the coverage gate, so it is demixed
        # (status ok) but kept out of the aggregate/plot, and the NTC assessment must flag it (ntc_reads)
        "id": "NEG_CONTAM",
        "pct_genome_ge_10x": [90, 100],
        "sample_type": "negative",
        "date": "2026-01-19",
        "site": "SiteB",
        "mix": {"BA.2": 1.0},
        "pairs": 3000,
    },
    {
        # underscores in the ID: guards against ID rewriting (the Snakemake pipeline turned _ into -)
        "id": "Site_1_week_2",
        "pct_genome_ge_10x": [88, 97],
        "sample_type": "sample",
        "date": "2026-02-02",
        "site": "Site_1",
        "mix": {"BA.5": 0.60, "XBB.1.5": 0.40},
        "pairs": 6500,
        "dropout": ["nCoV-2019_7", "nCoV-2019_22"],
    },
]

_SNV_RE = re.compile(r"([ACGT])(\d+)([ACGT])")
_PRIMER_NAME_RE = re.compile(r"^(.+?)_(?:LEFT|RIGHT|F|R)(?:_alt\d*)?$")
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
BASES = "ACGT"


# ----------------------------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------------------------
def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def read_fasta(path: Path) -> str:
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sample_rng(seed: int, name: str) -> random.Random:
    digest = hashlib.sha256(f"{seed}:{name}".encode()).digest()
    return random.Random(int.from_bytes(digest[:8], "big"))


def randint(rng: random.Random, lo: int, hi: int) -> int:
    """Inclusive integer in [lo, hi] using only rng.random() (stable across Python versions)."""
    return lo + int(rng.random() * (hi - lo + 1))


def write_gzip(path: Path, data: bytes) -> None:
    """Deterministic gzip: no filename, mtime 0, fixed level."""
    with (
        open(path, "wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as gz,
    ):
        gz.write(data)


# ----------------------------------------------------------------------------------------------------------------
# barcode extraction (one-off; reads the 100 MB Freyja usher_barcodes.csv)
# ----------------------------------------------------------------------------------------------------------------
def lineages_in_yml(path: Path) -> set[str]:
    names = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("- name:"):
                names.add(line.split(":", 1)[1].strip())
    return names


def extract_barcodes(src: Path, out: Path, lineages: list[str], lineages_yml: Path | None) -> None:
    if lineages_yml is not None:
        known = lineages_in_yml(lineages_yml)
        missing = [lin for lin in lineages if lin not in known]
        if missing:
            raise SystemExit(f"lineages not in lineages.yml: {missing}")
    csv.field_size_limit(sys.maxsize)
    wanted = set(lineages)
    rows: dict[str, list[str]] = {}
    with open(src, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)[1:]
        for row in reader:
            if row[0] in wanted:
                rows[row[0]] = row[1:]
    missing = [lin for lin in lineages if lin not in rows]
    if missing:
        raise SystemExit(f"lineages not in barcodes csv: {missing}")
    keep = [i for i in range(len(header)) if any(float(rows[lin][i]) != 0.0 for lin in lineages)]
    seen: dict[tuple, str] = {}
    for lin in lineages:
        key = tuple(rows[lin][i] for i in keep)
        if key in seen:
            raise SystemExit(
                f"{lin} has an identical barcode to {seen[key]}; choose another lineage"
            )
        seen[key] = lin
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow([""] + [header[i] for i in keep])
        for lin in lineages:
            w.writerow([lin] + [rows[lin][i] for i in keep])
    print(f"wrote {out}: {len(lineages)} lineages x {len(keep)} mutations", file=sys.stderr)


def read_barcodes(path: Path) -> dict[str, list[str]]:
    """Return {lineage: [mutation, ...]} for cells > 0."""
    out: dict[str, list[str]] = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)[1:]
        for row in reader:
            out[row[0]] = [m for m, v in zip(header, row[1:]) if float(v) > 0]
    return out


# ----------------------------------------------------------------------------------------------------------------
# lineage genomes
# ----------------------------------------------------------------------------------------------------------------
def apply_mutations(ref: str, muts: list[str]) -> tuple[str, dict]:
    """Apply Freyja barcode mutations (SNV 'A23063T'; deletions '-21633-9' / insertions '+...' if ever present)
    to ``ref``. Returns (genome, stats). Anything that cannot be applied is skipped and counted."""
    seq = list(ref)
    stats = {
        "applied": 0,
        "skipped_ref_mismatch": 0,
        "skipped_unsupported": 0,
        "skipped_conflict": 0,
    }
    touched: set[int] = set()
    snv_re = _SNV_RE
    for m in muts:
        match = snv_re.fullmatch(m)
        if not match:
            # Deletion / insertion notation. The shipped Freyja barcodes only contain SNVs (verified for
            # barcode 09_17_2026-00-46); the real notation is not guessed here.
            stats["skipped_unsupported"] += 1
            continue
        r, pos, a = match.group(1), int(match.group(2)), match.group(3)
        if pos < 1 or pos > len(ref) or ref[pos - 1] != r:
            stats["skipped_ref_mismatch"] += 1
            continue
        if pos in touched:
            stats["skipped_conflict"] += 1
            continue
        touched.add(pos)
        seq[pos - 1] = a
        stats["applied"] += 1
    return "".join(seq), stats


# ----------------------------------------------------------------------------------------------------------------
# primer scheme / amplicons
# ----------------------------------------------------------------------------------------------------------------
def load_amplicons(bed: Path, ref: str) -> list[dict]:
    """Amplicon = left-most F primer start .. right-most R primer end, grouped by the primer name minus _F/_R."""
    groups: dict[str, dict] = {}
    order: list[str] = []
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or line.startswith(("#", "track", "browser")):
                continue
            start, end, name, strand = int(f[1]), int(f[2]), f[3], f[5]
            m = _PRIMER_NAME_RE.match(name)
            if not m:  # scheme with an unsupported naming style (supported: *_F/_R, *_LEFT/_RIGHT, optional _alt)
                continue
            amp = m.group(1)
            g = groups.setdefault(amp, {"name": amp, "fwd": [], "rev": []})
            if amp not in order:
                order.append(amp)
            (g["fwd"] if strand == "+" else g["rev"]).append((start, end))
    amplicons = []
    for amp in order:
        g = groups[amp]
        if not g["fwd"] or not g["rev"]:
            continue
        fstart, fend = min(g["fwd"])
        rstart, rend = max(g["rev"], key=lambda t: t[1])
        amplicons.append(
            {
                "name": amp,
                "start": fstart,
                "end": rend,
                "fwd_primer": ref[fstart:fend],
                "rev_primer": revcomp(ref[rstart:rend]),  # as sequenced on the reverse strand
                "fwd_len": fend - fstart,
                "rev_len": rend - rstart,
            }
        )
    amplicons.sort(key=lambda a: a["start"])
    return amplicons


def amplicon_template(genome: str, ref: str, amp: dict) -> str:
    """Amplicon sequence of a lineage genome. The primer-binding ends carry the *primer* (= reference) sequence."""
    s, e = amp["start"], amp["end"]
    seq = genome[s:e]
    seq = (
        ref[s : s + amp["fwd_len"]]
        + seq[amp["fwd_len"] : len(seq) - amp["rev_len"]]
        + ref[e - amp["rev_len"] : e]
    )
    return seq


def uncovered_regions(amplicons: list[dict], dropped: set[str], genome_len: int) -> list[list[int]]:
    """0-based half-open genome intervals that no non-dropped amplicon spans (only returned inside dropped amplicons)."""
    covered = [False] * genome_len
    for a in amplicons:
        if a["name"] in dropped:
            continue
        for i in range(a["start"], a["end"]):
            covered[i] = True
    regions: list[list[int]] = []
    for a in amplicons:
        if a["name"] not in dropped:
            continue
        i = a["start"]
        while i < a["end"]:
            if not covered[i]:
                j = i
                while j < a["end"] and not covered[j]:
                    j += 1
                if regions and regions[-1][1] == i:
                    regions[-1][1] = j
                else:
                    regions.append([i, j])
                i = j
            else:
                i += 1
    return regions


# ----------------------------------------------------------------------------------------------------------------
# read simulation
# ----------------------------------------------------------------------------------------------------------------
def _mutate_and_qual(rng: random.Random, seq: str) -> tuple[str, str]:
    out = []
    qual = []
    n = len(seq)
    for i, base in enumerate(seq):
        # 3' tail degradation: slightly lower quality at the end of the read
        q_ok = 37 if i < n - 25 else (32 if rng.random() < 0.7 else 27)
        if rng.random() < ERROR_RATE:
            alts = [b for b in BASES if b != base]
            out.append(alts[int(rng.random() * 3)])
            qual.append(
                12 if rng.random() < 0.6 else 30
            )  # most errors carry a low Q (filtered by ivar -q 20)
        else:
            out.append(base)
            qual.append(q_ok)
    return "".join(out), "".join(chr(33 + q) for q in qual)


def simulate_sample(
    rng: random.Random,
    sample: dict,
    genomes: dict[str, str],
    ref: str,
    amplicons: list[dict],
) -> tuple[bytes, bytes, dict]:
    dropped = set(sample.get("dropout", []))
    names = {a["name"] for a in amplicons}
    unknown = dropped - names
    if unknown:
        raise SystemExit(
            f"{sample['id']}: dropout amplicons not in primer scheme: {sorted(unknown)}"
        )

    # per-amplicon depth weights (log-normal jitter); dropped amplicons get zero weight
    weights = []
    for a in amplicons:
        w = 0.0 if a["name"] in dropped else math.exp(DEPTH_JITTER_SD * rng.gauss(0.0, 1.0))
        weights.append(w)
    cum = []
    tot = 0.0
    for w in weights:
        tot += w
        cum.append(tot)

    lineages = list(sample["mix"])
    lcum = []
    t = 0.0
    for lin in lineages:
        t += sample["mix"][lin]
        lcum.append(t)

    templates = {
        (lin, a["name"]): amplicon_template(genomes[lin], ref, a)
        for lin in lineages
        for a in amplicons
    }

    r1_lines: list[str] = []
    r2_lines: list[str] = []
    counts = {lin: 0 for lin in lineages}
    per_amp = {a["name"]: 0 for a in amplicons}
    for n in range(1, sample["pairs"] + 1):
        ai = min(bisect.bisect_left(cum, rng.random() * tot), len(amplicons) - 1)
        while weights[ai] == 0.0:  # guard against landing exactly on a zero-weight boundary
            ai = (ai + 1) % len(amplicons)
        amp = amplicons[ai]
        li = min(bisect.bisect_left(lcum, rng.random() * t), len(lineages) - 1)
        lin = lineages[li]
        tmpl = templates[(lin, amp["name"])]
        counts[lin] += 1
        per_amp[amp["name"]] += 1

        flen = round(rng.gauss(FRAG_MEAN, FRAG_SD))
        flen = max(FRAG_MIN, min(FRAG_MAX, flen))
        flen = min(flen, len(tmpl))
        u = rng.random()
        if u < P_ANCHOR_5:
            fs = 0
        elif u < P_ANCHOR_5 + P_ANCHOR_3:
            fs = len(tmpl) - flen
        else:
            fs = randint(rng, 0, len(tmpl) - flen)
        frag = tmpl[fs : fs + flen]
        if rng.random() < 0.5:  # which strand is read 1
            frag = revcomp(frag)
        rl = min(READ_LEN, len(frag))
        s1, q1 = _mutate_and_qual(rng, frag[:rl])
        s2, q2 = _mutate_and_qual(rng, revcomp(frag)[:rl])
        name = f"{sample['id']}:{n:06d}"
        r1_lines.append(f"@{name} 1:N:0:ACGTACGT\n{s1}\n+\n{q1}\n")
        r2_lines.append(f"@{name} 2:N:0:ACGTACGT\n{s2}\n+\n{q2}\n")

    info = {"pairs_per_lineage": counts, "pairs_per_amplicon": per_amp}
    return "".join(r1_lines).encode(), "".join(r2_lines).encode(), info


# ----------------------------------------------------------------------------------------------------------------
# main generation
# ----------------------------------------------------------------------------------------------------------------
def generate(
    out: Path,
    seed: int = 42,
    primer_scheme: str = "midnight",
    reference: Path = DEFAULT_REFERENCE,
    barcodes: Path = DEFAULT_BARCODES,
    primer_dir: Path = DEFAULT_PRIMER_DIR,
    path_prefix: str = "tests/data",
    log=lambda msg: print(msg, file=sys.stderr),
) -> dict:
    out = Path(out)
    (out / "fastq").mkdir(parents=True, exist_ok=True)
    ref = read_fasta(reference)
    bed = primer_dir / f"{primer_scheme}.bed"
    amplicons = load_amplicons(bed, ref)
    barcode_muts = read_barcodes(barcodes)

    needed = sorted({lin for s in SAMPLES for lin in s["mix"]})
    genomes: dict[str, str] = {}
    mut_stats: dict[str, dict] = {}
    for lin in needed:
        if lin not in barcode_muts:
            raise SystemExit(f"lineage {lin} not in {barcodes}")
        genomes[lin], mut_stats[lin] = apply_mutations(ref, barcode_muts[lin])
        skipped = sum(v for k, v in mut_stats[lin].items() if k.startswith("skipped"))
        log(
            f"{lin}: applied {mut_stats[lin]['applied']} mutations, skipped {skipped} ({mut_stats[lin]})"
        )

    # barcodes_test.csv is an input; copy it next to the other outputs (unless it already is there)
    dest = out / "barcodes_test.csv"
    if Path(barcodes).resolve() != dest.resolve():
        shutil.copyfile(barcodes, dest)

    truth_samples: dict[str, dict] = {}
    sheet_rows = []
    for s in SAMPLES:
        rng = sample_rng(seed, s["id"])
        b1, b2, info = simulate_sample(rng, s, genomes, ref, amplicons)
        write_gzip(out / "fastq" / f"{s['id']}_R1.fastq.gz", b1)
        write_gzip(out / "fastq" / f"{s['id']}_R2.fastq.gz", b2)
        sheet_rows.append(
            {
                "sample": s["id"],
                "fastq_1": f"{path_prefix}/fastq/{s['id']}_R1.fastq.gz",
                "fastq_2": f"{path_prefix}/fastq/{s['id']}_R2.fastq.gz",
                "sample_type": s["sample_type"],
                "collection_date": s["date"],
                "site": s["site"],
            }
        )
        dropped = set(s.get("dropout", []))
        is_low = s["pairs"] < LOW_COVERAGE_MAX_PAIRS
        is_negative = s["sample_type"] == "negative"
        # Expectations under the pipeline defaults (nextflow.config qc_* / ntc_* params)
        qc_reasons = []
        if len(dropped) > QC_MAX_DROPPED_DEFAULT:
            qc_reasons.append("amplicon_dropout")
        if s["id"] == "LOWCOV_D":
            qc_reasons += [
                "low_coverage",
                "amplicon_dropout",
            ]  # ~0.6x: every amplicon is below the threshold
        if s["id"] == "NEG_CONTAM":
            qc_reasons.append("ntc_reads")
        qc_reasons = sorted(set(qc_reasons))
        demix_ok = not is_low
        truth_samples[s["id"]] = {
            "sample_type": s["sample_type"],
            "true_lineage_fractions": s["mix"],
            "realised_pairs_per_lineage": info["pairs_per_lineage"],
            "n_read_pairs": s["pairs"],
            "expected_demix_status": "ok" if demix_ok else "skipped_low_coverage",
            "expected_qc_reasons": qc_reasons,
            "expected_qc_status": (
                "FAIL"
                if set(qc_reasons) & {"low_coverage", "ntc_reads"}
                else ("WARN" if qc_reasons else "PASS")
            ),
            # negatives are demixed (if they pass the gate) but never aggregated (PLAN D8)
            "expected_in_aggregate": demix_ok and not is_negative,
            # no consensus is built for negatives; an unmappable sample gives an all-N consensus
            "expected_consensus_status": "none" if is_negative else ("empty" if is_low else "ok"),
            "dropout_amplicons": sorted(dropped),
            "expected_zero_coverage_regions_0based_half_open": uncovered_regions(
                amplicons, dropped, len(ref)
            ),
            "expected_pct_genome_ge_10x_range": s["pct_genome_ge_10x"],
            "abundance_tolerance": 0.05 if not is_low else None,
        }
    truth = {
        "seed": seed,
        "reference": reference.name,
        "primer_scheme": primer_scheme,
        "primer_bed_sha256": sha256_file(bed),
        "barcode_file": "barcodes_test.csv",
        "lineage_mutation_application": mut_stats,
        "samples": truth_samples,
    }
    (out / "truth.json").write_text(json.dumps(truth, indent=2, sort_keys=True) + "\n")

    with open(out / "samplesheet.csv", "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["sample", "fastq_1", "fastq_2", "sample_type", "collection_date", "site"],
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(sheet_rows)

    (out / "primer_scheme.txt").write_text(f"{primer_scheme}\n")
    return truth


# ----------------------------------------------------------------------------------------------------------------
# probe-capture dataset (PLAN D13): tests/data_capture
# ----------------------------------------------------------------------------------------------------------------
CAPTURE_SAMPLES = [
    {
        "id": "CAP_MIX",
        "sample_type": "sample",
        "date": "2026-03-02",
        "site": "SiteC",
        "mix": {"BA.2": 0.70, "AY.4": 0.30},
        "pairs": 6000,
        "pct_genome_ge_10x": [93, 100],
    },
    {
        "id": "CAP_MIX3",
        "sample_type": "sample",
        "date": "2026-03-09",
        "site": "SiteC",
        "mix": {"BA.1": 0.50, "XBB.1.5": 0.30, "KP.3": 0.20},
        "pairs": 7000,
        "pct_genome_ge_10x": [93, 100],
    },
    {
        "id": "CAP_PURE",
        "sample_type": "sample",
        "date": "2026-03-16",
        "site": "SiteD",
        "mix": {"JN.1": 1.0},
        "pairs": 8000,
        "pct_genome_ge_10x": [93, 100],
    },
    {
        "id": "CAP_LOWCOV",
        "sample_type": "sample",
        "date": "2026-03-23",
        "site": "SiteD",
        "mix": {"BA.2": 1.0},
        "pairs": 60,
        "pct_genome_ge_10x": [0, 5],
    },
    {
        "id": "CAP_NEG",
        "sample_type": "negative",
        "date": "2026-03-02",
        "site": "SiteC",
        "mix": {"B": 1.0},
        "pairs": 3,
        "pct_genome_ge_10x": [0, 5],
    },
]
CAPTURE_FRAG_MEAN, CAPTURE_FRAG_SD, CAPTURE_FRAG_MIN, CAPTURE_FRAG_MAX = 300, 70, 150, 450
CAPTURE_BIN = 1000  # depth varies per 1 kb bin (log-normal)
CAPTURE_DEPTH_JITTER_SD = 0.15
CAPTURE_DUP_RATE = 0.17  # share of read pairs that are exact PCR copies of an earlier pair
CAPTURE_DUP_FRACTION_RANGE = [
    0.10,
    0.25,
]  # expected share of reads flagged duplicate by samtools markdup


def simulate_capture_sample(
    rng: random.Random, sample: dict, genomes: dict[str, str]
) -> tuple[bytes, bytes, dict]:
    """2x150 reads from fragments starting uniformly over the genome, plus exact PCR duplicates."""
    ref_len = len(next(iter(genomes.values())))
    nbins = -(-ref_len // CAPTURE_BIN)
    weights = [math.exp(CAPTURE_DEPTH_JITTER_SD * rng.gauss(0.0, 1.0)) for _ in range(nbins)]
    cum = []
    tot = 0.0
    for w in weights:
        tot += w
        cum.append(tot)
    lineages = list(sample["mix"])
    lcum = []
    t = 0.0
    for lin in lineages:
        t += sample["mix"][lin]
        lcum.append(t)

    unique: list[
        tuple[str, str, str, str]
    ] = []  # (seq1, qual1, seq2, qual2) of every non-duplicate pair
    r1_lines: list[str] = []
    r2_lines: list[str] = []
    counts = {lin: 0 for lin in lineages}
    n_dups = 0
    for n in range(1, sample["pairs"] + 1):
        if unique and rng.random() < CAPTURE_DUP_RATE:
            s1, q1, s2, q2 = unique[int(rng.random() * len(unique))]
            n_dups += 1
        else:
            li = min(bisect.bisect_left(lcum, rng.random() * t), len(lineages) - 1)
            lin = lineages[li]
            genome = genomes[lin]
            counts[lin] += 1
            flen = round(rng.gauss(CAPTURE_FRAG_MEAN, CAPTURE_FRAG_SD))
            flen = max(CAPTURE_FRAG_MIN, min(CAPTURE_FRAG_MAX, flen))
            b = min(bisect.bisect_left(cum, rng.random() * tot), nbins - 1)
            centre = b * CAPTURE_BIN + int(rng.random() * CAPTURE_BIN)
            start = max(
                0, min(ref_len, centre - flen // 2)
            )  # fragments may overhang the genome ends
            frag = genome[start : start + flen]
            if len(frag) < 60:  # too little genome left: take the terminal stretch instead
                frag = genome[-flen:] if start > ref_len // 2 else genome[:flen]
            if rng.random() < 0.5:
                frag = revcomp(frag)
            rl = min(READ_LEN, len(frag))
            s1, q1 = _mutate_and_qual(rng, frag[:rl])
            s2, q2 = _mutate_and_qual(rng, revcomp(frag)[:rl])
            unique.append((s1, q1, s2, q2))
        name = f"{sample['id']}:{n:06d}"
        r1_lines.append(f"@{name} 1:N:0:ACGTACGT\n{s1}\n+\n{q1}\n")
        r2_lines.append(f"@{name} 2:N:0:ACGTACGT\n{s2}\n+\n{q2}\n")
    info = {"pairs_per_lineage": counts, "duplicate_pairs": n_dups}
    return "".join(r1_lines).encode(), "".join(r2_lines).encode(), info


def generate_capture(
    out: Path,
    seed: int = 42,
    reference: Path = DEFAULT_REFERENCE,
    barcodes: Path = DEFAULT_BARCODES,
    path_prefix: str = "tests/data_capture",
    log=lambda msg: print(msg, file=sys.stderr),
) -> dict:
    out = Path(out)
    (out / "fastq").mkdir(parents=True, exist_ok=True)
    ref = read_fasta(reference)
    barcode_muts = read_barcodes(barcodes)
    genomes: dict[str, str] = {}
    mut_stats: dict[str, dict] = {}
    for lin in sorted({lin for s in CAPTURE_SAMPLES for lin in s["mix"]}):
        if lin not in barcode_muts:
            raise SystemExit(f"lineage {lin} not in {barcodes}")
        genomes[lin], mut_stats[lin] = apply_mutations(ref, barcode_muts[lin])
    truth_samples: dict[str, dict] = {}
    sheet_rows = []
    for s in CAPTURE_SAMPLES:
        rng = sample_rng(seed, "capture:" + s["id"])
        b1, b2, info = simulate_capture_sample(rng, s, genomes)
        write_gzip(out / "fastq" / f"{s['id']}_R1.fastq.gz", b1)
        write_gzip(out / "fastq" / f"{s['id']}_R2.fastq.gz", b2)
        sheet_rows.append(
            {
                "sample": s["id"],
                "fastq_1": f"{path_prefix}/fastq/{s['id']}_R1.fastq.gz",
                "fastq_2": f"{path_prefix}/fastq/{s['id']}_R2.fastq.gz",
                "sample_type": s["sample_type"],
                "collection_date": s["date"],
                "site": s["site"],
            }
        )
        is_low = s["pairs"] < LOW_COVERAGE_MAX_PAIRS
        is_negative = s["sample_type"] == "negative"
        qc_reasons = ["low_coverage"] if s["id"] == "CAP_LOWCOV" else []  # never amplicon_dropout
        truth_samples[s["id"]] = {
            "sample_type": s["sample_type"],
            "true_lineage_fractions": s["mix"],
            "realised_pairs_per_lineage": info["pairs_per_lineage"],
            "n_read_pairs": s["pairs"],
            "n_duplicate_pairs": info["duplicate_pairs"],
            "expected_duplicate_fraction_range": None
            if s["pairs"] < 1000
            else CAPTURE_DUP_FRACTION_RANGE,
            "expected_demix_status": "ok" if not is_low else "skipped_low_coverage",
            "expected_qc_reasons": qc_reasons,
            "expected_qc_status": "FAIL" if qc_reasons else "PASS",
            "expected_in_aggregate": (not is_low) and not is_negative,
            "expected_consensus_status": "none" if is_negative else ("empty" if is_low else "ok"),
            "expected_pct_genome_ge_10x_range": s["pct_genome_ge_10x"],
            "abundance_tolerance": 0.05 if not is_low else None,
            "dropout_amplicons": [],
        }
    truth = {
        "seed": seed,
        "protocol": "capture",
        "reference": reference.name,
        "barcode_file": "../data/barcodes_test.csv",
        "lineage_mutation_application": mut_stats,
        "samples": truth_samples,
    }
    (out / "truth.json").write_text(json.dumps(truth, indent=2, sort_keys=True) + "\n")
    with open(out / "samplesheet.csv", "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["sample", "fastq_1", "fastq_2", "sample_type", "collection_date", "site"],
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(sheet_rows)
    return truth


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--protocol",
        choices=["amplicon", "capture"],
        default="amplicon",
        help="amplicon: tests/data (default); capture: tests/data_capture (probe-capture reads, PLAN D13)",
    )
    ap.add_argument(
        "--out", type=Path, default=None, help="output directory (default per protocol)"
    )
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--primer-scheme", default="midnight", help="name of a BED in assets/primer_schemes"
    )
    ap.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    ap.add_argument(
        "--barcodes", type=Path, default=DEFAULT_BARCODES, help="small Freyja barcodes .csv (input)"
    )
    ap.add_argument("--path-prefix", default=None, help="prefix for FASTQ paths in samplesheet.csv")
    ap.add_argument(
        "--extract-barcodes",
        type=Path,
        metavar="USHER_BARCODES_CSV",
        help="one-off: write --barcodes",
    )
    ap.add_argument(
        "--lineages-yml",
        type=Path,
        help="Freyja lineages.yml, to validate lineage names on extraction",
    )
    args = ap.parse_args(argv)

    if args.extract_barcodes:
        extract_barcodes(args.extract_barcodes, args.barcodes, BARCODE_LINEAGES, args.lineages_yml)
        return 0
    if args.protocol == "capture":
        generate_capture(
            args.out or REPO / "tests" / "data_capture",
            seed=args.seed,
            reference=args.reference,
            barcodes=args.barcodes,
            path_prefix=args.path_prefix or "tests/data_capture",
        )
        return 0
    generate(
        args.out or REPO / "tests" / "data",
        seed=args.seed,
        primer_scheme=args.primer_scheme,
        reference=args.reference,
        barcodes=args.barcodes,
        path_prefix=args.path_prefix or "tests/data",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
