"""Deterministic synthetic time-series fixture for wasteline (no Docker, standard library only).

Writes, into ``--out`` (default ``tests/data_timeseries``):

``demix/<sample>.demix``
    Synthetic ``freyja demix`` result files in the exact format Freyja 2.x writes (a ``<TAB><label>`` first
    line, then ``summarized``, ``lineages``, ``abundances``, ``resid``, ``coverage``, ``pathogen`` rows).
``samplesheet.csv``
    ``sample,demix,sample_type,collection_date,site,viral_load`` (``demix`` is relative to the repository
    root, like the ``fastq_*`` columns of ``tests/data/samplesheet.csv``). The sheet is also the metadata
    file of the time-series scripts.
``truth.json``
    The generating model and the noise-free expected abundances per site and date.

Model: 3 sites x 8 weekly dates. Lineage A (``BA.2``) declines and lineage B (``XBB.1.5``) rises
logistically: the share of B among A+B is ``logistic(R_PER_DAY * (day - t0_site))`` with the same true
log-odds slope ``R_PER_DAY`` at every site and a site-specific midpoint ``t0_site``. A third lineage C
(``JN.1``) is stable at ``C_SHARE`` and a minor noise lineage (``BA.1``) appears at random. About 30 dated
samples (1-2 per site and date), plus deliberate edge cases: a negative control, a sample without a
collection date, a low-coverage sample and several blank viral loads.

The ``freyja_aggregated.tsv`` that the other test inputs use is the real ``freyja aggregate`` of these
files (needs Docker, linux/amd64): ``--aggregate`` runs ``relabel_demix.py`` + ``freyja aggregate`` in the
pinned container + ``clean_aggregate_labels.py --sort``, exactly as FREYJA_AGGREGATE_PLOT does.

Usage::

    python tests/make_timeseries_fixture.py --out tests/data_timeseries --seed 42
    python tests/make_timeseries_fixture.py --out tests/data_timeseries --seed 42 --aggregate

Output is byte-identical for the same seed (only ``random.Random.random()`` and ``gauss`` are used).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FREYJA_IMAGE = "quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0"

# Real lineages of the Freyja hierarchy and the `summarized` group Freyja puts them in.
LINEAGE_A = "BA.2"
LINEAGE_B = "XBB.1.5"
LINEAGE_C = "JN.1"
LINEAGE_NOISE = "BA.1"
GROUPS = {
    "BA.2": "Omicron",
    "BA.1": "Omicron",
    "XBB.1.5": "XBB.1.5* [Omicron (XBB.1.5.X)]",
    "JN.1": "BA.2.86* (BA.2.86.X)",
}

R_PER_DAY = 0.10  # true daily log-odds growth advantage of B over A (the same at every site)
C_SHARE = 0.10  # constant share of the stable lineage C
SITES = {"SiteA": 21, "SiteB": 28, "SiteC": 35}  # site -> day (from the first date) where B = A
SITE_LOAD_MEDIAN = {"SiteA": 2.0e5, "SiteB": 1.0e6, "SiteC": 5.0e5}  # gene copies per litre
START = dt.date(2026, 1, 5)  # a Monday
N_DATES = 8
LOGIT_NOISE_SD = 0.25  # sample-level noise on the log-odds of B vs A
LOAD_SD = 0.9  # lognormal sd of the viral load
PATHOGEN = "SARS-CoV-2"
SHEET_COLUMNS = ["sample", "demix", "sample_type", "collection_date", "site", "viral_load"]


def logistic(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def true_abundances(site: str, day: int) -> dict[str, float]:
    """Noise-free abundance of A, B and C at a site on a day (days since START)."""
    q = logistic(R_PER_DAY * (day - SITES[site]))
    return {
        LINEAGE_A: (1 - C_SHARE) * (1 - q),
        LINEAGE_B: (1 - C_SHARE) * q,
        LINEAGE_C: C_SHARE,
    }


def sample_abundances(rng: random.Random, site: str, day: int) -> dict[str, float]:
    """Noisy abundance: noise on the A/B log-odds, a jittered C, an occasional minor noise lineage."""
    q = logistic(R_PER_DAY * (day - SITES[site]) + rng.gauss(0.0, LOGIT_NOISE_SD))
    c = max(0.0, C_SHARE + rng.gauss(0.0, 0.02))
    noise = rng.random() * 0.03 if rng.random() < 0.5 else 0.0
    rest = 1.0 - c - noise
    out = {LINEAGE_A: rest * (1 - q), LINEAGE_B: rest * q, LINEAGE_C: c, LINEAGE_NOISE: noise}
    out = {name: round(value, 8) for name, value in out.items() if value >= 0.005}
    total = sum(out.values())
    return {name: round(value / total, 8) for name, value in out.items()}


def demix_text(sample: str, abundances: dict[str, float], coverage: float, resid: float) -> str:
    """The text of a Freyja 2.x ``.demix`` file."""
    ordered = sorted(abundances.items(), key=lambda kv: (-kv[1], kv[0]))
    groups: dict[str, float] = {}
    for name, value in ordered:
        groups[GROUPS[name]] = groups.get(GROUPS[name], 0.0) + value
    summarized = sorted(groups.items(), key=lambda kv: (-kv[1], kv[0]))
    return (
        f"\t{sample}.variants.tsv\n"
        f"summarized\t{[(g, v) for g, v in summarized]!r}\n"
        f"lineages\t{' '.join(name for name, _ in ordered)}\n"
        f"abundances\t{' '.join(f'{value:.8f}' for _, value in ordered)}\n"
        f"resid\t{resid!r}\n"
        f"coverage\t{coverage!r}\n"
        f"pathogen\t{PATHOGEN}\n"
    )


def build(seed: int):
    """Return (files {relative name: text}, sheet rows, truth dict)."""
    rng = random.Random(seed)
    files: dict[str, str] = {}
    rows: list[dict[str, str]] = []
    expected: list[dict] = []
    n = 0

    def add(sample, abundances, coverage, sample_type, date, site, load):
        files[f"demix/{sample}.demix"] = demix_text(
            sample, abundances, coverage, round(1.0 + 4.0 * rng.random(), 6)
        )
        rows.append(
            {
                "sample": sample,
                "demix": f"tests/data_timeseries/demix/{sample}.demix",
                "sample_type": sample_type,
                "collection_date": date,
                "site": site,
                "viral_load": "" if load is None else str(load),
            }
        )

    for s_idx, site in enumerate(SITES):
        for d_idx in range(N_DATES):
            day = 7 * d_idx
            date = (START + dt.timedelta(days=day)).isoformat()
            truth = true_abundances(site, day)
            expected.append(
                {
                    "site": site,
                    "collection_date": date,
                    "day": day,
                    "abundance": {k: round(v, 6) for k, v in truth.items()},
                }
            )
            n_samples = 2 if (s_idx * N_DATES + d_idx) % 4 == 1 else 1
            for rep in range(1, n_samples + 1):
                n += 1
                sample = f"{site}_W{d_idx + 1:02d}_{rep}"
                load = round(math.exp(math.log(SITE_LOAD_MEDIAN[site]) + rng.gauss(0.0, LOAD_SD)))
                if n % 7 == 0:
                    load = None  # unknown viral load
                coverage = round(94.0 + 5.0 * rng.random(), 4)
                add(
                    sample,
                    sample_abundances(rng, site, day),
                    coverage,
                    "sample",
                    date,
                    site,
                    load,
                )

    # Edge cases ---------------------------------------------------------------------------------
    # negative control: must never reach a time series
    add("NEG_1", {LINEAGE_A: 1.0}, 99.0, "negative", "2026-01-12", "SiteA", 1)
    # sample whose collection date is unknown
    add("NODATE_1", sample_abundances(rng, "SiteB", 14), 97.0, "sample", "", "SiteB", 300000)
    # low genome coverage (below Freyja's plot --mincov 60): excluded from time series
    add("LOWCOV_1", {LINEAGE_A: 0.5, LINEAGE_B: 0.5}, 42.5, "sample", "2026-01-26", "SiteC", 400000)

    rows.sort(key=lambda r: r["sample"])
    truth = {
        "seed": seed,
        "model": (
            "share of B among A+B = logistic(r * (day - t0_site)); A = (1 - c) * (1 - share), "
            "B = (1 - c) * share, C = c; day counted from first_date"
        ),
        "lineage_declining": LINEAGE_A,
        "lineage_rising": LINEAGE_B,
        "lineage_stable": LINEAGE_C,
        "lineage_noise": LINEAGE_NOISE,
        "true_logodds_growth_per_day": R_PER_DAY,
        "true_relative_growth_serial_interval_5_5": round(math.exp(5.5 * R_PER_DAY) - 1, 6),
        "stable_share": C_SHARE,
        "site_midpoint_day": SITES,
        "first_date": START.isoformat(),
        "n_dates": N_DATES,
        "logit_noise_sd": LOGIT_NOISE_SD,
        "excluded_samples": {
            "NEG_1": "negative control",
            "NODATE_1": "no collection date",
            "LOWCOV_1": "coverage 42.5 < 60",
        },
        "expected_site_date_abundance": expected,
    }
    return files, rows, truth


def write_fixture(out: Path, seed: int) -> None:
    files, rows, truth = build(seed)
    if (out / "demix").exists():
        shutil.rmtree(out / "demix")
    (out / "demix").mkdir(parents=True, exist_ok=True)
    for name, text in files.items():
        (out / name).write_text(text)
    with open(out / "samplesheet.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SHEET_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (out / "truth.json").write_text(json.dumps(truth, indent=2, sort_keys=True) + "\n")


def aggregate_with_docker(out: Path) -> None:
    """Real `freyja aggregate` of the fixture, labelled and cleaned like FREYJA_AGGREGATE_PLOT."""
    bin_dir = REPO / "bin"
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "labelled").mkdir()
        with open(out / "samplesheet.csv", newline="") as handle:
            samples = [row["sample"] for row in csv.DictReader(handle)]
        for sample in samples:
            subprocess.run(
                [
                    sys.executable,
                    str(bin_dir / "relabel_demix.py"),
                    str(out / "demix" / f"{sample}.demix"),
                    "--sample",
                    sample,
                    "--output",
                    str(tmp / "labelled" / f"{sample}.demix"),
                ],
                check=True,
            )
        subprocess.run(
            [
                "docker",
                "run",
                "--rm",
                "--platform",
                "linux/amd64",
                "-v",
                f"{tmp}:/w",
                "-w",
                "/w",
                FREYJA_IMAGE,
                "freyja",
                "aggregate",
                "labelled/",
                "--ext",
                "demix",
                "--output",
                "aggregated.raw.tsv",
            ],
            check=True,
        )
        subprocess.run(
            [
                sys.executable,
                str(bin_dir / "clean_aggregate_labels.py"),
                str(tmp / "aggregated.raw.tsv"),
                str(out / "freyja_aggregated.tsv"),
                "--sort",
            ],
            check=True,
        )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", default=str(REPO / "tests" / "data_timeseries"))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--aggregate",
        action="store_true",
        help="also write freyja_aggregated.tsv with the real Freyja container (needs Docker)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    out = Path(args.out)
    write_fixture(out, args.seed)
    if args.aggregate:
        aggregate_with_docker(out)
    print(f"wrote time-series fixture to {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
