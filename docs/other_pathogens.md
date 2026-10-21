# Other pathogens (`--pathogen`)

wasteline was built for SARS-CoV-2, but Freyja also deconvolves other pathogens (`freyja update --pathogen X`), and
`--pathogen X` now works end to end for the single-sequence ones (PLAN D22). Everything is opt-in: with the default
`--pathogen SARS-CoV-2` nothing changes. Design as built: [wiring/multi_pathogen.md](wiring/multi_pathogen.md).

**Two things decide whether a pathogen works.** The genome must be one sequence (influenza has 8 segments, which is refused),
and you must supply the pathogen's **reference FASTA**: Freyja's barcodes are tied to the `reference.fasta` in the pathogen's
folder of the [Freyja-barcodes](https://github.com/andersen-lab/Freyja-barcodes) repository, and a different reference gives
silently wrong results. A primer BED for your amplicon panel (or `--protocol capture`) is the third ingredient; no primer
schemes are bundled for other pathogens.

## Support matrix

| `--pathogen`                                              | Status                       | What you need                                                                                                        |
| --------------------------------------------------------- | ---------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `SARS-CoV-2` (default)                                    | supported, validated         | nothing: reference, annotation, primer schemes and lineage groups are bundled                                        |
| `MEASLESgenome`                                           | supported, **validated**     | the reference is bundled (`NC_001498.1`); `--primer_bed` (or `--protocol capture`). Validation below                 |
| `MPX`, `RSVa`, `RSVb`, `DENV1`, `DENV2`, `DENV3`, `DENV4` | single sequence, **not run** | `--fasta` = the `reference.fasta` of the pathogen's Freyja-barcodes folder, `--primer_bed` (or `--protocol capture`) |
| `H5Nx`, `H5Nx-cattle`, `H1N1`, `H3N2`, `FLU-B-VIC`        | **unsupported**              | influenza has 8 segments (a multi-contig reference), see below                                                       |
| `MTB`                                                     | **unsupported**              | a 4.4 Mb bacterial genome, see below                                                                                 |

"Not run" means the code path is the same one that was validated on measles (nothing in it is pathogen-specific), but no
data of that pathogen has been run through it: expect to check the first results yourself. MPX (about 197 kb) in particular
is much larger than the genomes the pipeline has been timed on.

**Unsupported pathogens stop at start-up** with the reason, before any work is done:

- **Influenza** (`H5Nx`, `H5Nx-cattle`, `H1N1`, `H3N2`, `FLU-B-VIC`): the genome is 8 segments, so the reference has several
  contigs. wasteline supports single-contig references only: the depth QC, amplicon coverage, low-depth masking, the
  consensus and the primer checks all assume one contig. Supporting it needs per-segment (per-contig) handling in all of
  them, and per-segment Freyja calls (`freyja variants --refname`). The Freyja documentation also says influenza barcodes are
  available on request only. Any other pathogen whose `--fasta` has more than one sequence is refused in the same way.
- **MTB**: _Mycobacterium tuberculosis_ is a 4.4 Mb bacterial genome with feather-only Freyja barcodes. The pipeline is
  built and tested for viral genomes of 10-200 kb (per-position depth tables, QC, consensus and the dropout/mutation tables
  are not sized for 4.4 Mb) and nothing has been validated on bacterial data. Supporting it needs those steps scaled and
  checked on real MTB data.

## What changes for a pathogen other than SARS-CoV-2

| Part of the pipeline                                  | Behaviour                                                                                                                                                                                                                                                                                                           |
| ----------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Reference, annotation                                 | `conf/pathogens/<profile>.config`: `measles.config` bundles `assets/references/NC_001498.1.fasta` (no annotation). Any other pathogen uses `generic.config`: `--fasta` is **required** (the error says it must be the reference of the barcodes), `--gff` is optional.                                              |
| Primer scheme                                         | The bundled schemes are SARS-CoV-2 only. Without `--primer_bed` (or a per-sample `primer_scheme` path, or `--protocol capture` / `--skip_primer_trimming`) the run stops with a message instead of trimming with the wrong BED. The BED's first column must be the reference's sequence name (checked at start-up). |
| Freyja dataset (`--barcodes`)                         | `latest` works (`freyja update --pathogen X`). `container` does not (the Freyja container bundles SARS-CoV-2 only; refused). A dataset directory of `freyja update --pathogen X`, or a barcodes `.csv`, can be given. See "Freyja datasets" below.                                                                  |
| Lineage-group summary (`summary/lineage_summary.csv`) | **Off** unless `--lineage_groups <CSV>` is given (the bundled groups are SARS-CoV-2 Pango lineages); one `INFO` line says so. It also needs the dataset's lineage hierarchy (`<X>_lineages.yml`, only shipped for some pathogens): without it the summary is skipped with one warning. See "Known limits".          |
| Dropout impact, unexplained mutations                 | Work from the dataset's barcodes `.csv` (there is no `lineage_mutations.json`). The unexplained-mutation significance level 0.05 / genome length follows the reference length (a user-set `--unexplained_max_pvalue` is never rescaled).                                                                            |
| Primer scheme check (`--primer_scheme_check`)         | Skipped: it scores the reads against the bundled SARS-CoV-2 schemes, which do not fit another reference. Trimming and per-amplicon coverage use your BED as usual.                                                                                                                                                  |
| `--annotate_variants`, amino-acid columns             | Need `--gff` (an error says so). Amino-acid effects are only as good as the annotation you give; the ORF1ab slippage handling is name-based and does nothing without such a gene.                                                                                                                                   |
| Report                                                | The subtitle names the pathogen (`Wastewater MEASLESgenome lineage deconvolution`).                                                                                                                                                                                                                                 |
| `bin/reanalyse_runs.py`                               | Carries `--pathogen` over; only exercised for SARS-CoV-2.                                                                                                                                                                                                                                                           |

### Freyja datasets of other pathogens

What `freyja update --pathogen MEASLESgenome` (Freyja 2.0.6) actually writes, measured for this work:

```
MEASLESgenome_barcodes.csv          the barcodes (a .feather for MTB)
MEASLESgenome_lineages.yml          lineage hierarchy; only for the pathogens whose Freyja config names one (MPX, H5Nx-cattle, MEASLESgenome)
last_barcode_update_other.txt       one line per pathogen: "MEASLESgenome:2025-05-01"  (a date, no time)
pathogen_config.yml
```

There is no `usher_barcodes.*`, `last_barcode_update.txt`, `lineage_mutations.json`, `curated_lineages.json` or `lineages.yml`, and
no reference FASTA. wasteline's dataset code (`bin/freyja_dataset_info.py`) recognises these names: the barcode version is the
date (`2025-05-01`, ISO form `2025-05-01T00-00`, so two downloads on one day are the same dataset), and the barcode file,
hierarchy file and their checksums go into `pipeline_info/freyja_data.json` as for SARS-CoV-2.

Datasets are stored per pathogen: `<db_dir>/freyja` for SARS-CoV-2 (unchanged) and `<db_dir>/freyja_<pathogen>/<date>T00-00` for
the others, so `--barcodes 2025-05-01` can never pick up another pathogen's data. `--mode update_db --pathogen X` downloads one.

One upstream quirk: Freyja's `pathogen_config.yml` points `MPX` at a GitHub web page instead of the raw `lineages.yml`, so
`freyja update --pathogen MPX` saves an HTML page as `MPX_lineages.yml` (checked today). wasteline ignores a hierarchy file that
is not a YAML list, i.e. MPX has no usable hierarchy until Freyja fixes that.

## Running measles

Measles (`MEASLESgenome`) is the validated example. The Freyja tutorial uses the `NC_001498.1` reference (bundled) and the
[ARTIC measles v1.0.0](https://labs.primalscheme.com/detail/artic-measles/400/v1.0.0/) primer panel (47 amplicons, contig
`NC_001498.1`). The primer BED is **not bundled** (licence unclear): download it (Freyja ships a copy as
`docs/data/artic-measles-v1.0.0.bed`) and pass it.

```bash
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results \
    --pathogen MEASLESgenome \
    --primer_bed artic-measles-v1.0.0.bed --ivar_offset 3
```

- `--barcodes` defaults to `latest`: Freyja's measles barcodes are downloaded to `<db_dir>/freyja_MEASLESgenome/`. Pin them
  later with `--barcodes 2025-05-01` (the barcode date).
- `--ivar_offset 3` is what the Freyja tutorial uses for this panel (`ivar trim -x 3`); the default for a custom BED is 0.
- Reads must be paired (the samplesheet needs `fastq_2`).
- For a panel with another name format, the BED's 4th column must still tell left from right primers (`..._LEFT_0` / `..._RIGHT_0`,
  `..._F` / `..._R`, or the strand column). The primalscheme `<scheme>_<n>_LEFT_<alt>` names of the ARTIC panels were not understood
  before PLAN D22 (all primers collapsed into one amplicon) and now are.
- The results are the usual ones (`summary/freyja_aggregated.tsv`, `qc_summary.tsv`, `amplicon_coverage.tsv`, consensus genomes,
  `report.html`, ...). Lineage names are those of the barcodes, e.g. `MEASLESgenome-H1`.

## Validation (measles, simulated reads)

**Data**: the simulated reads of the Freyja measles tutorial (`docs/data` of the Freyja repository; made with the Bygul amplicon
simulator from the ARTIC measles panel): `GCA_031128185.1-simulated.fastq` (a pure H1 sample) and `measles-mixed-simulated.fastq`
(H1 20 % + D9 80 %). They are **single-end** 150 bp reads, but the pipeline needs pairs, so each read became a pair with its reverse
complement as read 2 (a fully overlapping pair; script not kept in the repo). The reads stay outside the repository.

**Run**: `--mode full --pathogen MEASLESgenome --primer_bed artic-measles-v1.0.0.bed --ivar_offset 3`, `--barcodes latest` (the
dataset of 2025-05-01, 15 lineages), docker, Freyja 2.0.6, 49 tasks, about 1 minute on a laptop. A second run with
`--barcodes 2025-05-01` (the pinned dataset), `--freyja_boot` (20 bootstraps) and `--lineage_groups` gave identical abundances.

| Sample            | Truth            | Recovered (`freyja_aggregated.tsv`)                          | Genome >= 10x | Mean depth | Amplicons (dropped) | QC   |
| ----------------- | ---------------- | ------------------------------------------------------------ | ------------- | ---------- | ------------------- | ---- |
| `H1_pure`         | H1 100 %         | `MEASLESgenome-H1` **100.0 %**                               | 92.7 %        | 196x       | 47 (0)              | PASS |
| `mix_H1_20_D9_80` | H1 20 %, D9 80 % | `MEASLESgenome-D9` **77.9 %**, `MEASLESgenome-H1` **22.1 %** | 93.4 %        | 618x       | 47 (0)              | PASS |

The pure sample is recovered exactly; the mixture is within 2.1 percentage points of the simulated proportions (bootstrap
intervals with 20 draws: D9 0.769-0.775, H1 0.225-0.231; the point estimates lie just outside them, which is what so few draws
give). Consensus genomes are 15,894 bp with 6.6-7.3 % N. The `unexplained_mutations` table lists positions where the sample differs
from the reference at ~100 % but no barcode lists the change (presumably the simulated genomes' private mutations; not investigated).

Not run: the custom barcode file `measles-wg-barcode.csv` of the tutorial; any pathogen other than measles and SARS-CoV-2; real
measles data. The SARS-CoV-2 wiring was re-checked with the 52 stub nf-tests (no snapshot changed), the single docker nf-test (synthetic
data, passes) and the pytest suite; `tests/e2e.sh` was not re-run for this change.

## Known limits

- **Lineage groups for measles do not work with today's Freyja data**: the lineage hierarchy (`MEASLESgenome_lineages.yml`) names
  lineages `MEASLES-H1` while the barcodes call them `MEASLESgenome-H1`. A groups file naming either is accepted, but no
  observed lineage matches, so every abundance lands in `other` (the lineage names remain visible in the `ungrouped_lineages`
  column and in `freyja_aggregated.tsv`). Use the aggregate until Freyja's naming is made consistent.
- Freyja's `summarized` column (what `freyja plot` groups by) is `Other` for every measles sample, so its plots have one group.
- The bundled amplicon, QC and report defaults were tuned on SARS-CoV-2 wastewater (D17 QC thresholds: PASS at 60 % of the genome
  > = 10x); review them for your pathogen and panel.
- Only measles has been run. For the other single-sequence pathogens read the first results critically (reference, primer BED
  contig names, barcodes date).
