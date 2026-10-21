# Other pathogens (`--pathogen`)

wasteline was built for SARS-CoV-2, but Freyja also deconvolves other pathogens (`freyja update --pathogen X`), and
`--pathogen X` now works end to end for the single-sequence ones (PLAN D22). RSV-A and RSV-B (`RSVa`, `RSVb`), measles and the
five influenza pathogens (`H1N1`, `H3N2`, `FLU-B-VIC`, `H5Nx`, `H5Nx-cattle`; the HA segment only, see "Running influenza")
are validated and bundle their reference (PLAN D23, D24). Everything is opt-in: with the default
`--pathogen SARS-CoV-2` nothing changes. Design as built: [wiring/multi_pathogen.md](wiring/multi_pathogen.md).

**Two things decide whether a pathogen works.** The reference must be one sequence (a multi-sequence FASTA is refused; Freyja's
influenza barcodes cover only the HA segment, which is one sequence of 1.7-1.9 kb, so influenza qualifies), and it must be the
**reference FASTA** that Freyja's barcodes are tied to: the `reference.fasta` in the pathogen's folder of the
[Freyja-barcodes](https://github.com/andersen-lab/Freyja-barcodes) repository (a different reference gives silently wrong
results). For RSV-A, RSV-B, measles and influenza it is bundled; for the others you pass `--fasta`. A primer BED for your amplicon
panel (or `--protocol capture`) is the third ingredient; no primer schemes are bundled for other pathogens.

## Support matrix

| `--pathogen`                              | Status                                         | What you need                                                                                                              |
| ----------------------------------------- | ---------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `SARS-CoV-2` (default)                    | supported, validated                           | nothing: reference, annotation, primer schemes and lineage groups are bundled                                              |
| `RSVa`                                    | supported, **validated**                       | reference and annotation are bundled (`PP109421.1`); `--primer_bed` (or `--protocol capture`). Validation below            |
| `RSVb`                                    | supported, **validated**                       | reference and annotation are bundled (`OP975389.1`); `--primer_bed` (or `--protocol capture`). Validation below            |
| `MEASLESgenome`                           | supported, **validated**                       | the reference is bundled (`NC_001498.1`); `--primer_bed` (or `--protocol capture`). Validation below                       |
| `MPX`, `DENV1`, `DENV2`, `DENV3`, `DENV4` | single sequence, **not run**                   | `--fasta` = the `reference.fasta` of the pathogen's Freyja-barcodes folder, `--primer_bed` (or `--protocol capture`)       |
| `H1N1`                                    | supported (HA only), **validated**             | reference and annotation are bundled (`CY121680.1`, HA); `--primer_bed` (or `--protocol capture`). See "Running influenza" |
| `H3N2`                                    | supported (HA only), **validated**             | bundled `CY163680.1` (HA); `--primer_bed` (or `--protocol capture`)                                                        |
| `FLU-B-VIC`                               | supported (HA only), **validated (quick run)** | bundled `KX058884.1` (HA); `--primer_bed` (or `--protocol capture`)                                                        |
| `H5Nx`, `H5Nx-cattle`                     | supported (HA only), **validated (quick run)** | bundled `AF144305.1` (HA, one file for both); `--primer_bed` (or `--protocol capture`)                                     |
| `MTB`                                     | **unsupported, not wanted**                    | a 4.4 Mb bacterial genome, see below                                                                                       |

"Not run" means the code path is the same one that was validated on measles (nothing in it is pathogen-specific), but no
data of that pathogen has been run through it: expect to check the first results yourself. MPX (about 197 kb) in particular
is much larger than the genomes the pipeline has been timed on. "Validated" means simulated reads of known composition (no
real reads of these pathogens were available), see the validation sections.

**Unsupported pathogens stop at start-up** with the reason, before any work is done:

- **Whole-genome influenza** is not a refused pathogen any more (PLAN D24): the five influenza pathogens run on their HA segment.
  What stays refused is a `--fasta` with more than one sequence (any pathogen), because the depth QC, amplicon coverage,
  low-depth masking, consensus and primer checks assume one contig; see "What is supported and what Route B would add" below.
- **MTB**: _Mycobacterium tuberculosis_ is a 4.4 Mb bacterial genome with feather-only Freyja barcodes. The pipeline is
  built and tested for viral genomes of 10-200 kb (per-position depth tables, QC, consensus and the dropout/mutation tables
  are not sized for 4.4 Mb) and nothing has been validated on bacterial data. Supporting it needs those steps scaled and
  checked on real MTB data.

## What changes for a pathogen other than SARS-CoV-2

| Part of the pipeline                                  | Behaviour                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Reference, annotation                                 | `conf/pathogens/<profile>.config`: `measles.config` bundles `assets/references/NC_001498.1.fasta` (no annotation), `rsva.config` / `rsvb.config` bundle `PP109421.1` / `OP975389.1` with a `.gff3` each, and the influenza profiles bundle the HA record of each pathogen (`CY121680.1`, `CY163680.1`, `KX058884.1`, `AF144305.1`) with a `.gff3` each. Any other pathogen uses `generic.config`: `--fasta` is **required** (the error says it must be the reference of the barcodes), `--gff` is optional. |
| Primer scheme                                         | The bundled schemes are SARS-CoV-2 only. Without `--primer_bed` (or a per-sample `primer_scheme` path, or `--protocol capture` / `--skip_primer_trimming`) the run stops with a message instead of trimming with the wrong BED. The BED's first column must be the reference's sequence name (checked at start-up). Lines of contigs that are not in the reference (a whole-genome influenza panel against the HA reference) are dropped with one warning, see "Running influenza".                         |
| Freyja dataset (`--barcodes`)                         | `latest` works (`freyja update --pathogen X`). `container` does not (the Freyja container bundles SARS-CoV-2 only; refused). A dataset directory of `freyja update --pathogen X`, or a barcodes `.csv`, can be given. See "Freyja datasets" below.                                                                                                                                                                                                                                                          |
| Lineage-group summary (`summary/lineage_summary.csv`) | **Off** unless `--lineage_groups <CSV>` is given (the bundled groups are SARS-CoV-2 Pango lineages); one `INFO` line says so. It also needs the dataset's lineage hierarchy (`<X>_lineages.yml`, only shipped for some pathogens): without it the summary is skipped with one warning. See "Known limits".                                                                                                                                                                                                  |
| Dropout impact, unexplained mutations                 | Work from the dataset's barcodes `.csv` (there is no `lineage_mutations.json`). The unexplained-mutation significance level 0.05 / genome length follows the reference length (a user-set `--unexplained_max_pvalue` is never rescaled).                                                                                                                                                                                                                                                                    |
| Primer scheme check (`--primer_scheme_check`)         | Skipped: it scores the reads against the bundled SARS-CoV-2 schemes, which do not fit another reference. Trimming and per-amplicon coverage use your BED as usual.                                                                                                                                                                                                                                                                                                                                          |
| `--annotate_variants`, amino-acid columns             | Need `--gff` (an error says so). Amino-acid effects are only as good as the annotation you give; the ORF1ab slippage handling is name-based and does nothing without such a gene.                                                                                                                                                                                                                                                                                                                           |
| Report                                                | The subtitle names the pathogen (`Wastewater MEASLESgenome lineage deconvolution`).                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `bin/reanalyse_runs.py`                               | Carries `--pathogen` over; only exercised for SARS-CoV-2.                                                                                                                                                                                                                                                                                                                                                                                                                                                   |

### Freyja datasets of other pathogens

What `freyja update --pathogen MEASLESgenome` (Freyja 2.0.6) actually writes, measured for this work:

```
MEASLESgenome_barcodes.csv          the barcodes (a .feather for MTB)
MEASLESgenome_lineages.yml          lineage hierarchy; only for the pathogens whose Freyja config names one (MPX, H5Nx-cattle, MEASLESgenome)
last_barcode_update_other.txt       one line per pathogen: "MEASLESgenome:2025-05-01"  (a date, no time)
pathogen_config.yml
```

The same files (`RSVa_barcodes.csv` / `RSVb_barcodes.csv`, no `lineages.yml`, `last_barcode_update_other.txt` reading
`RSVa:2024-03-15`) are what `freyja update --pathogen RSVa` and `RSVb` write (Freyja 2.0.6, checked for PLAN D23). In the
Freyja-barcodes repository each pathogen has a folder with dated subfolders (RSVa: `2024-03-08`, `2024-03-15`) and `latest/`; the
folders hold `barcode.csv`, `reference.fasta`, `auspice_tree.json` and `barcode.html`, so a barcode set and its reference are
published together. `freyja update` downloads only `barcode.csv` (plus a lineage hierarchy where the config names one), never the
reference.

There is no `usher_barcodes.*`, `last_barcode_update.txt`, `lineage_mutations.json`, `curated_lineages.json` or `lineages.yml`, and
no reference FASTA. wasteline's dataset code (`bin/freyja_dataset_info.py`) recognises these names: the barcode version is the
date (`2025-05-01`, ISO form `2025-05-01T00-00`, so two downloads on one day are the same dataset), and the barcode file,
hierarchy file and their checksums go into `pipeline_info/freyja_data.json` as for SARS-CoV-2.

Datasets are stored per pathogen: `<db_dir>/freyja` for SARS-CoV-2 (unchanged) and `<db_dir>/freyja_<pathogen>/<date>T00-00` for
the others, so `--barcodes 2025-05-01` can never pick up another pathogen's data. `--mode update_db --pathogen X` downloads one.

The influenza pathogens write the same files (`H1N1_barcodes.csv`, ...; checked for PLAN D24 with Freyja 2.0.6; only `H5Nx-cattle` also
writes a `_lineages.yml`, a valid YAML list), see "Running influenza".

One upstream quirk: Freyja's `pathogen_config.yml` points `MPX` at a GitHub web page instead of the raw `lineages.yml`, so
`freyja update --pathogen MPX` saves an HTML page as `MPX_lineages.yml` (checked today). wasteline ignores a hierarchy file that
is not a YAML list, i.e. MPX has no usable hierarchy until Freyja fixes that.

## Running RSV-A and RSV-B

`--pathogen RSVa` and `--pathogen RSVb` bundle the reference and its annotation, so the only extra ingredient is the library
description. The primer BED is **not bundled** (licences of RSV panels are unclear or unknown): pass your panel's BED, whose first
column must be the reference sequence name (`PP109421.1` for RSV-A, `OP975389.1` for RSV-B), or say that the library is not
amplicon-based.

```bash
# amplicon library: your panel's primer BED (set --ivar_offset if your panel needs it)
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results \
    --pathogen RSVa --primer_bed my_rsv_a_panel.bed

# probe-capture or shotgun library: no primers to trim, no amplicon dropout analysis
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results \
    --pathogen RSVb --protocol capture
```

Without `--primer_bed` and without `--protocol capture` / `--skip_primer_trimming` the run stops at start-up with a message saying
so. `--barcodes` defaults to `latest` (downloaded with `freyja update --pathogen RSVa`, stored in `<db_dir>/freyja_RSVa/<date>T00-00`;
the dataset of 2024-03-15 today); the lineages are named `RSVa-A.D.1`, `RSVb-B.D.4.1`, ... `--barcodes container` is refused.
`--annotate_variants` works (gene and amino-acid effects from the bundled GFF3: NS1, NS2, N, P, M, SH, G, F, M2-1, M2-2, L).

- **Reference**: the bundled sequences are GenBank `PP109421.1` (hRSV/A/England/397/2017, 15,225 bp) and `OP975389.1`
  (hRSV/B/Australia/VIC-RCH056/2019, 15,222 bp). They are base-for-base identical to the `reference.fasta` of the `RSVa/latest`
  and `RSVb/latest` folders of Freyja-barcodes, which names them by their GISAID ids (EPI_ISL_412866, EPI_ISL_1653999); the
  public GenBank records were found by strain name and compared. Checked for PLAN D23: every SNV column of the barcodes
  (RSVa: 1499 columns, 23 lineages; RSVb: 915 columns, 16 lineages) has the stated reference base at the stated 1-based position
  of this reference (100 % of columns, no exceptions).
- **Annotation**: the NCBI GFF3 of the same records. NCBI's file has only `gene` and `CDS` lines, which `bcftools csq` cannot use
  (it reported "No usable transcripts"), so the bundled `.gff3` adds one `mRNA` and one `exon` line per gene (same coordinates;
  the file header says so). Gene names are NCBI's.
- **Lineage summary**: off (no bundled lineage groups, and Freyja ships no lineage hierarchy for RSV): one `INFO` line says so. With
  `--lineage_groups` given it is skipped with one `WARN` line (the dataset has no `lineages.yml`). Dropout impact and unexplained
  mutations work from the barcodes CSV. The primer scheme check is skipped (it scores against the SARS-CoV-2 schemes).
- The lineage tree of Freyja's RSV barcodes is nested (`RSVa-A.D` contains `RSVa-A.D.1`, ...): expect abundances of about 0.2 % on a
  few neighbouring lineages in every sample.

## Running measles

Measles (`MEASLESgenome`) is the validated example (RSV is just below). The Freyja tutorial uses the `NC_001498.1` reference (bundled) and the
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

## Validation (RSV-A and RSV-B, simulated reads)

**Why simulated**: no real RSV reads were available. The reads come from `tests/simulate_pathogen_mixture.py` (a small development
tool, not part of the test suite): each lineage is a pseudo-genome made by applying every mutation of that lineage's row of the
dataset's barcodes CSV (RSVa: dataset 2024-03-15, 23 lineages; RSVb: 16 lineages) to the bundled reference, and read pairs
(2 x 150 bp, mean depth about 500x, 0.2 % substitution errors, fragment length 350 +/- 40) are drawn from the lineages in fixed
proportions. Three well separated lineages were used per pathogen (RSVa A.1 / A.3 / A.D.1, RSVb B.1 / B.3 / B.D.4.1; 315-844
mutations apart) and one closely related pair (A.D.1 / A.D.1.2, 41 apart; B.D.4.1 / B.D.4.1.1, 25 apart). Because the
pseudo-genomes are exactly the barcodes plus the reference, this tests the plumbing and the arithmetic of the pipeline (mapping,
variant calling, Freyja, consensus, QC), **not** how Freyja performs on real RSV diversity, where genomes carry private mutations
and recombinants exist.

**Runs** (docker, Freyja 2.0.6, one run per pathogen with six samples, `--mode full --barcodes latest --annotate_variants`): four
shotgun-style samples with `--protocol capture`, and two amplicon-style samples (per-sample `protocol` / `primer_scheme` cells in
the samplesheet) simulated from a toy tiling BED derived from the genome (280 bp amplicons every 240 bp, 24 nt primers, names
`RSV_<n>_LEFT_0` / `RSV_<n>_RIGHT_0`, 63 amplicons, one fragment per amplicon). The toy primers overlap at every junction and the
genome ends are not amplified, so 4.0 % of the genome (8 bp at each of 62 junctions plus the ends) has no coverage after primer
trimming, which is why the amplicon samples reach 96.0 % of the genome at 10x. RSV-B took about 1.5 minutes (117 tasks) on a
laptop that was also running other work.

| Sample (RSV-A)            | Truth                          | Recovered (`freyja_aggregated.tsv`)                      | Mean depth | Genome >= 10x | QC   |
| ------------------------- | ------------------------------ | -------------------------------------------------------- | ---------- | ------------- | ---- |
| `pure_AD1` (shotgun)      | A.D.1 100 %                    | A.D.1 **99.8 %** (+ A.D.2.2.1 0.2 %)                     | 493x       | 99.95 %       | PASS |
| `mix2` (shotgun)          | A.3 30 %, A.D.1 70 %           | A.3 **29.6 %**, A.D.1 **69.4 %** (+ 0.9 % in 3 others)   | 494x       | 99.91 %       | PASS |
| `mix3` (shotgun)          | A.1 20 %, A.3 30 %, A.D.1 50 % | A.1 **19.3 %**, A.3 **30.4 %**, A.D.1 **49.6 %**         | 491x       | 99.96 %       | PASS |
| `close` (shotgun)         | A.D.1 70 %, A.D.1.2 30 %       | A.D.1 **70.7 %**, A.D.1.2 **29.1 %** (+ A.D.2.2.1 0.2 %) | 494x       | 99.93 %       | PASS |
| `amp_pure_AD1` (amplicon) | A.D.1 100 %                    | A.D.1 **100.0 %**                                        | 387x       | 96.00 %       | PASS |
| `amp_mix3` (amplicon)     | A.1 20 %, A.3 30 %, A.D.1 50 % | A.1 **20.2 %**, A.3 **28.7 %**, A.D.1 **50.4 %**         | 387x       | 96.00 %       | PASS |

| Sample (RSV-B)             | Truth                            | Recovered (`freyja_aggregated.tsv`)                    | Mean depth | Genome >= 10x | QC   |
| -------------------------- | -------------------------------- | ------------------------------------------------------ | ---------- | ------------- | ---- |
| `pure_BD41` (shotgun)      | B.D.4.1 100 %                    | B.D.4.1 **99.8 %** (+ B.D.4.1.1 0.2 %)                 | 494x       | 99.92 %       | PASS |
| `mix2` (shotgun)           | B.3 30 %, B.D.4.1 70 %           | B.3 **30.3 %**, B.D.4.1 **69.5 %** (+ B.D.4.1.1 0.2 %) | 493x       | 99.94 %       | PASS |
| `mix3` (shotgun)           | B.1 20 %, B.3 30 %, B.D.4.1 50 % | B.1 **19.9 %**, B.3 **29.4 %**, B.D.4.1 **50.4 %**     | 493x       | 99.89 %       | PASS |
| `close` (shotgun)          | B.D.4.1 70 %, B.D.4.1.1 30 %     | B.D.4.1 **69.7 %**, B.D.4.1.1 **30.3 %**               | 493x       | 99.94 %       | PASS |
| `amp_pure_BD41` (amplicon) | B.D.4.1 100 %                    | B.D.4.1 **100.0 %**                                    | 387x       | 96.02 %       | PASS |
| `amp_mix3` (amplicon)      | B.1 20 %, B.3 30 %, B.D.4.1 50 % | B.1 **20.4 %**, B.3 **29.8 %**, B.D.4.1 **49.6 %**     | 387x       | 96.02 %       | PASS |

- The largest error of any true lineage is 1.3 percentage points (RSV-A `amp_mix3`, A.3 28.7 % for 30 %); most are under 1. The
  closely related pairs are separated correctly. The amplicon samples have 63 amplicons, none dropped, so their dropout tables
  (amplicon list, `dropout_impact*.tsv`) show zero loss.
- **Consensus genomes** (majority base, IUPAC code where the lineages differ): compared with the pseudo-genome of the dominant
  lineage, no position in any of the 12 samples has a base incompatible with it, and the pure samples are identical to it
  (outside the Ns). Mixed samples carry 24-382 IUPAC codes (the positions where the mixed lineages differ). N content is
  0.04-0.11 % for shotgun samples and 4.0 % for amplicon samples (the uncovered positions above).
- `unexplained_mutations`: 7 positions flagged in the RSV-A run and 2 in RSV-B, all in mixed samples, with an observed frequency
  0.11-0.17 above the fitted expectation (for example RSV-A position 5076 in G: 0.62-0.65 observed for 0.50 expected in both
  `mix3` samples). Not investigated. Plausible causes are the sampling noise of the simulation (mates of one pair overlap) and
  reference-mapping bias against the most divergent lineage (RSV-A A.1 is 861 mutations from the reference and three of the nine
  flagged rows are in G, the most variable gene); no pipeline fault was found.
- `annotated_variants.tsv` lists 1462 (RSV-A) and 994 (RSV-B) annotated variants with gene and amino-acid effect (for example
  `NS1 S7S`, `G L142S`, `M2-2 H41Y`). The report subtitle reads `Wastewater RSVa lineage deconvolution` (`RSVb`). One cosmetic
  leftover: the console banner still says "wastewater SARS-CoV-2 lineage abundance" for every pathogen.
- The RSV-A run was repeated with `--lineage_groups` (a one-line groups file): identical aggregate, one `WARN` line, no
  `lineage_summary.csv`.

Not run: real RSV reads (amplicon or capture); a real RSV primer panel (the BED above is a toy tiling, not a published scheme);
`--freyja_boot`; a comparison of `latest/` with the newest dated folder of Freyja-barcodes; genomes with private mutations or
recombinants; `tests/e2e.sh` for these pathogens.

## Running influenza (HA segment only)

Freyja's influenza barcodes (`H1N1`, `H3N2`, `FLU-B-VIC`, `H5Nx`, `H5Nx-cattle`) describe **the HA segment only**. Each pathogen's
folder in Freyja-barcodes holds a single-sequence HA `reference.fasta` (1.7-1.9 kb), which is bundled, so influenza runs on the existing
single-contig code (PLAN D24, "Route A"). Freyja deconvolves the subtype's lineages from HA mutations alone, whatever the other seven
segments do; reassortants and the other segments are not analysed. If you need whole-genome consensus sequences or per-segment QC, see
"What is supported and what Route B would add" below.

| `--pathogen`  | Bundled reference (`assets/references/`)                | Length  | Lineages | Barcode columns checked | Barcode dataset (`latest`, today) |
| ------------- | ------------------------------------------------------- | ------- | -------- | ----------------------- | --------------------------------- |
| `H1N1`        | `CY121680.1` (A/California/07/2009 H1N1 HA) + `.gff3`   | 1752 bp | 15       | 140 / 140 (100 %)       | 2025-07-03                        |
| `H3N2`        | `CY163680.1` (A/Wisconsin/67/2005 H3N2 HA) + `.gff3`    | 1737 bp | 32       | 227 / 227 (100 %)       | 2025-12-09                        |
| `FLU-B-VIC`   | `KX058884.1` (B/Brisbane/60/2008 HA) + `.gff3`          | 1885 bp | 8        | 67 / 67 (100 %)         | 2025-07-03                        |
| `H5Nx`        | `AF144305.1` (A/Goose/Guangdong/1/96 H5N1 HA) + `.gff3` | 1760 bp | 11       | 523 / 523 (100 %)       | 2024-09-12                        |
| `H5Nx-cattle` | `AF144305.1` (the same file as `H5Nx`)                  | 1760 bp | 70       | 653 / 666 (98.0 %)      | 2024-01-01                        |

- **The references are Freyja's own sequences.** Each is the `latest/reference.fasta` of the pathogen's folder in Freyja-barcodes,
  sequence for sequence identical to the GenBank record (compared with NCBI efetch): the headers of `H1N1`, `H3N2` and `FLU-B-VIC`
  name `CY121680.1`, `CY163680.1` and `KX058884.1` (the `FLU-B-VIC` header is the bare accession); the `H5Nx` header has no accession
  (`Influenza A virus (A/Goose/Guangdong/1/96(H5N1)) hemagglutinin (HA) gene, complete cds.`) and its sequence is identical to GenBank
  `AF144305.1`, found by strain name. The bundled FASTA files are the NCBI records, whose header starts with the accession (the
  first word is the contig name; Freyja's `H5Nx` header would have given the contig name `Influenza`). `H5Nx` and `H5Nx-cattle`
  have the identical reference file.
- **Reference bases vs the barcodes.** "Barcode columns checked" counts the single-nucleotide columns of the barcodes whose stated
  reference base equals the base at that 1-based position of the bundled reference. All columns of `H1N1` (one more column, `A163R`,
  is an ambiguity code and not checked), `H3N2`, `FLU-B-VIC` and `H5Nx` match. **`H5Nx-cattle` has 13 of 666 columns (2.0 %) whose
  stated base differs from the reference**: `A486G`, `G1299A`, `G1071A`, `A448G`, `C1338A`, `G483A`, `G1028A`, `T927C`, `T1113C`,
  `A684G`, `A1338C`, `A915G`, `C1143T` (the reference has `G`, `A`, `A`, `G`, `T`, `T`, `C`, `C`, `A`, `G`, `T`, `G`, `T` there).
  `C1338A` is in 33 of the 70 lineages, the others in 1-6. Freyja looks up each barcode column by its name (`<ref><position><alt>`)
  in the sample's variants (not verified in Freyja's code), so on this reference `C1338A` can never be observed (`T1338A` would be) and
  these columns are lost for the lineages that carry them: an inconsistency of Freyja's `H5Nx-cattle` barcodes with the file in their own
  folder, not something wasteline can fix. The simulator skips them too.
- **Annotation**: the NCBI GFF3 of the same records. It has only `gene` and `CDS` lines (the HA CDS, gene name `HA`), which
  `bcftools csq` cannot use, so the bundled `.gff3` adds an `mRNA` and an `exon` line with the same coordinates (as for RSV; the file
  header says so). `--annotate_variants` works: the amino-acid effects are numbered from the HA start codon (sequential, not the H1/H3
  numbering used in the literature).
- **One run per subtype.** H1 and H3 HA are too different to map to one reference, and the samplesheet has no per-sample pathogen. A
  sample with several subtypes needs one run per `--pathogen` (the same samplesheet each time); the reads of the other subtypes, and
  of the seven other segments of every subtype, do not map to the HA reference and are ignored.
- **No primer BED is bundled.** Give your panel's BED with `--primer_bed` (or per sample in the samplesheet), or `--protocol capture`.
  A whole-genome influenza panel lists the primers of all eight segments: the lines of the contigs that are not in the HA reference are
  **dropped** once, with one warning (below), and the filtered BED is what trimming, amplicon coverage and the provenance record use.
  A BED none of whose lines is on the reference's contig is still refused with the usual error.
- **Freyja datasets**: `latest` is downloaded on first use to `~/.wasteline/db/freyja_<pathogen>/<date>T00-00` (see "Freyja datasets of
  other pathogens"). `freyja update --pathogen H1N1|H3N2|FLU-B-VIC|H5Nx` writes `<P>_barcodes.csv`, `last_barcode_update_other.txt`
  and `pathogen_config.yml` (no lineage hierarchy); `H5Nx-cattle` also writes `H5Nx-cattle_lineages.yml`, a valid YAML list of 70
  entries whose names equal the barcode rows, so `--lineage_groups` works for it. In that hierarchy `H5Nx-2.3.4.4` has 65 children and
  `H5Nx-2.3.4.4b.B3.13` 35, but the entry `H5Nx-2.3.4.4b` lists only the alias `H5Nx-2.3.4.4.b` as its child (an inconsistent alias, as
  for measles), so a group `H5Nx-2.3.4.4b*` matches that one lineage only: write `H5Nx-2.3.4.4*`. No lineage groups are bundled for
  influenza: the summary is off unless you give `--lineage_groups`.

```bash
# amplicon panel; the BED may list all eight segments (only the HA lines are used, the others are dropped with a warning)
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results_h3n2 \
    --pathogen H3N2 --primer_bed my_flu_wholegenome_panel.bed --annotate_variants

# probe-capture or shotgun library: no primers to trim, no amplicon dropout analysis
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results_h1n1 \
    --pathogen H1N1 --protocol capture

# samples with H1N1 and H3N2 reads: one run per subtype (each ignores the other subtype's reads)
nextflow run charlesfoster/wasteline -profile docker --input samplesheet.csv --outdir results_h1n1 --pathogen H1N1 --protocol capture
nextflow run charlesfoster/wasteline -profile docker --input samplesheet.csv --outdir results_h3n2 --pathogen H3N2 --protocol capture
```

The BED's first column must be the HA contig for its lines to be kept: the accession of the bundled record (`CY121680.1`,
`CY163680.1`, `KX058884.1`, `AF144305.1`). A panel BED that names segment 4 as `4` or `seg4` has to be renamed to that accession
for the HA lines (the other seven contigs are dropped anyway). The warning looks like this:

```
WARN: Primer BED: 10 line(s) on 2 contig(s) that are not in the reference FASTA 'CY121680.1.fasta' (seg4, seg6) were dropped; the reference's own lines are used for trimming and amplicon coverage. Expected for a whole-genome panel run against one segment (for example influenza HA).
```

The filtered copy (same file name; comment, blank and `track` lines are kept) is written once, at start-up, to
`<workDir>/wasteline_primer_beds/<12 hex>/`, and is the file `ivar trim`, `amplicon_coverage.py` and the provenance checksum in
`pipeline_info/run_info.json` use. A BED whose contigs all are in the reference is passed on untouched (the same file, so SARS-CoV-2
runs are unchanged).

**What an HA-only run reports.**

- `mapped_pct` is informational. It is low whenever the library holds reads of other segments or other subtypes (about 10 % in the
  validation sample with 8-fold more reads of another virus), and **no QC rule uses it** (the sample status comes from the share of the
  HA reference at >= 10x, dropped amplicons and the primer checks). Depth, `Genome >= 10x`, the coverage tracks and the consensus all
  describe the HA segment only; the consensus genome is the HA sequence (1.7-1.9 kb).
- The number of HA reads is small. At 500x a sample is about 2,900 read pairs; at 50x it is 290. See the validation for what that
  does to the accuracy.
- The QC settings that are percentages of the genome (PASS at 60 % of the reference at >= 10x, FAIL below 30 %; the low-confidence
  band 25-50 %), the significance level of the unexplained-mutation test (0.05 / reference length, about 2.9e-5 here) and the report's
  coverage track (binned to the reference length) follow the reference length and were used unchanged. **The negative-control
  rule does not scale**: a negative control FAILs above `--ntc_max_genome_coverage_pct` (10 %) of the genome at >= 1x, which is
  about 175 bases, i.e. one or two read pairs, on a 1.75 kb reference (and `--ntc_max_mapped_reads` 1000 is about 85x of mean depth
  there). In the validation a negative control with 12 mapped reads (1x depth, 59.7 % of HA at >= 1x) failed with `ntc_reads`. Raise
  `--ntc_max_genome_coverage_pct` / `--ntc_max_mapped_reads` for influenza HA if a few stray reads should not fail it. Nothing was
  changed in the code for this.
- With a whole-genome panel the amplicon table, `--qc_max_dropped_amplicons` (more than 1 dropped amplicon = WARN) and the dropout impact
  cover the HA amplicons only; a panel with only two or three HA amplicons turns one more drop into a WARN.

## Validation (influenza HA, simulated reads)

**Why simulated**: no real influenza reads were available. The reads come from `tests/simulate_pathogen_mixture.py`, as for RSV: each
lineage is a pseudo-genome (the reference with every single-nucleotide mutation of that lineage's row of the pathogen's barcodes CSV;
ambiguity columns are skipped, and for `H5Nx-cattle` the 13 columns above), and read pairs (2 x 150 bp, 0.2 % substitution errors,
fragment length 350 +/- 40) are drawn from the lineages in fixed proportions. Because the pseudo-genomes are exactly the barcodes plus
the reference, this tests the plumbing and the arithmetic, **not** how Freyja performs on real influenza diversity (private mutations,
reassortment, real primer behaviour). Lineage names are shown without the `<pathogen>-` prefix.

**Runs** (docker, Freyja 2.0.6, `--mode full --barcodes latest --annotate_variants`, `--protocol capture` as the run default and a
per-sample `protocol` / `primer_scheme` cell where a sample is amplicon-style): H1N1 and H3N2 with 11 samples each (203 tasks, about
1.5 minutes on a laptop), FLU-B-VIC, H5Nx and H5Nx-cattle with 3 samples each (63-65 tasks). Samples: pure, a 30/70 mixture, a
20/30/50 mixture, a close lineage pair (H1N1 5a / 5a.1, 3 mutations apart; H3N2 G.1.1 / G.1.1.1), the 20/30/50 mixture at 50x and 20x,
an amplicon-style pure and 20/30/50 sample simulated from a toy tiling BED (280 bp amplicons every 240 bp, 24 nt primers, 7 amplicons,
`FLU_<n>_LEFT_0` / `_RIGHT_0`), the amplicon 20/30/50 reads listed a second time with a **whole-genome-style BED** that also has dummy
lines for `seg4` and `seg6` (10 extra lines; the run printed the warning above), and the 20/30/50 sample with 8-fold more reads of
another virus mixed in (RSV-A reads standing in for the other segments).

H1N1 (reference `CY121680.1`, dataset 2025-07-03):

| Sample                | Truth                                                      | Recovered (`freyja_aggregated.tsv`)                               | Mean depth | Genome >= 10x | Consensus N | QC   |
| --------------------- | ---------------------------------------------------------- | ----------------------------------------------------------------- | ---------- | ------------- | ----------- | ---- |
| `pure`                | 5a.2a.1 100 %                                              | 5a.2a.1 **99.8 %**                                                | 492x       | 99.4 %        | 0.63 %      | PASS |
| `mix2`                | 6B.2 30 %, 5a.2a.1 70 %                                    | 6B.2 **31.1 %**, 5a.2a.1 **67.4 %** (+ 5a.2a 0.9 %, 5a.2 0.3 %)   | 491x       | 99.5 %        | 0.46 %      | PASS |
| `mix3`                | 6B.2 20 %, 1 30 %, 5a.2a.1 50 %                            | 6B.2 **22.4 %**, 1 **28.1 %**, 5a.2a.1 **47.1 %** (+ 5 1.7 %)     | 493x       | 99.4 %        | 0.57 %      | PASS |
| `close`               | 5a 70 %, 5a.1 30 %                                         | 5a **65.6 %**, 5a.1 **34.2 %**                                    | 493x       | 99.4 %        | 0.63 %      | PASS |
| `pure_50x`            | 5a.2a.1 100 %                                              | 5a.2a.1 **100.0 %**                                               | 49x        | 94.9 %        | 5.08 %      | PASS |
| `mix3_50x`            | as `mix3`                                                  | 6B.2 **20.0 %**, 1 **29.2 %**, 5a.2a.1 **50.8 %**                 | 49x        | 93.1 %        | 6.91 %      | PASS |
| `mix3_20x`            | as `mix3`                                                  | 6B.2 **19.4 %**, 1 **27.1 %**, 5a.2a.1 **53.6 %**                 | 20x        | 85.6 %        | 14.38 %     | PASS |
| `amp_pure`            | 5a.2a.1 100 %                                              | 5a.2a.1 **99.8 %**                                                | 387x       | 92.7 %        | 7.31 %      | PASS |
| `amp_mix3`            | as `mix3`                                                  | 6B.2 **20.0 %**, 1 **29.6 %**, 5a.2a.1 **46.7 %** (+ 5a.2a 3.4 %) | 387x       | 92.7 %        | 7.31 %      | PASS |
| `amp_mix3_wgbed`      | as `amp_mix3` (same reads, BED with `seg4` / `seg6` lines) | identical to `amp_mix3`                                           | 387x       | 92.7 %        | 7.31 %      | PASS |
| `mix3_plus_othersegs` | as `mix3`, plus 8x more reads of another virus             | identical to `mix3`; **Mapped % 10.3**                            | 493x       | 99.4 %        | 0.57 %      | PASS |

H3N2 (reference `CY163680.1`, dataset 2025-12-09):

| Sample                | Truth                                                      | Recovered (`freyja_aggregated.tsv`)                                                          | Mean depth | Genome >= 10x | Consensus N | QC   |
| --------------------- | ---------------------------------------------------------- | -------------------------------------------------------------------------------------------- | ---------- | ------------- | ----------- | ---- |
| `pure`                | G.1 100 %                                                  | G.1 **99.6 %**                                                                               | 482x       | 99.4 %        | 0.58 %      | PASS |
| `mix2`                | D 30 %, G.1 70 %                                           | D **31.1 %**, G.1 **68.4 %** (+ F 0.3 %)                                                     | 485x       | 99.7 %        | 0.35 %      | PASS |
| `mix3`                | D 20 %, G.1 30 %, K 50 %                                   | D **23.4 %**, G.1 **29.1 %**, K **45.2 %** (+ J.2 1.6 %)                                     | 482x       | 99.4 %        | 0.63 %      | PASS |
| `close`               | G.1.1 70 %, G.1.1.1 30 %                                   | G.1.1 **71.6 %**, G.1.1.1 **27.9 %**                                                         | 479x       | 99.5 %        | 0.46 %      | PASS |
| `pure_50x`            | G.1 100 %                                                  | G.1 **100.0 %**                                                                              | 49x        | 95.0 %        | 5.01 %      | PASS |
| `mix3_50x`            | as `mix3`                                                  | D **20.6 %**, G.1 **28.6 %**, K **47.4 %** (+ J.2.4 3.4 %)                                   | 48x        | 94.2 %        | 5.81 %      | PASS |
| `mix3_20x`            | as `mix3`                                                  | D **9.5 %**, G.1 **25.8 %**, K **46.7 %** (+ B.1.2.1 6.6 %, J.2 4.8 %, F 3.9 %, G.1.3 2.7 %) | 19x        | 84.2 %        | 15.77 %     | PASS |
| `amp_pure`            | G.1 100 %                                                  | G.1 **99.8 %** (+ C.1 0.2 %)                                                                 | 389x       | 93.5 %        | 6.51 %      | PASS |
| `amp_mix3`            | as `mix3`                                                  | D **20.5 %**, G.1 **30.8 %**, K **48.3 %**                                                   | 385x       | 93.5 %        | 6.51 %      | PASS |
| `amp_mix3_wgbed`      | as `amp_mix3` (same reads, BED with `seg4` / `seg6` lines) | identical to `amp_mix3`                                                                      | 385x       | 93.5 %        | 6.51 %      | PASS |
| `mix3_plus_othersegs` | as `mix3`, plus 8x more reads of another virus             | identical to `mix3`; **Mapped % 10.2**                                                       | 482x       | 99.4 %        | 0.63 %      | PASS |

Quick runs (3 samples each: pure, 30/70, 20/30/50, shotgun at 500x):

| Pathogen (dataset)         | Sample | Truth                                                   | Recovered                                                                                        | Genome >= 10x | QC   |
| -------------------------- | ------ | ------------------------------------------------------- | ------------------------------------------------------------------------------------------------ | ------------- | ---- |
| `FLU-B-VIC` (2025-07-03)   | `pure` | V1A.3 100 %                                             | V1A.3 **99.8 %**                                                                                 | 99.5 %        | PASS |
|                            | `mix2` | V1B 30 %, V1A.3a.2 70 %                                 | V1B **30.9 %**, V1A.3a.2 **68.1 %** (+ V1A.3 0.7 %)                                              | 99.5 %        | PASS |
|                            | `mix3` | V1B 20 %, V1A.2 30 %, V1A.3a.2 50 %                     | V1B **20.8 %**, V1A.2 **29.9 %**, V1A.3a.2 **49.1 %**                                            | 99.6 %        | PASS |
| `H5Nx` (2024-09-12)        | `pure` | 2.3.4.4b 100 %                                          | 2.3.4.4b **99.8 %**                                                                              | 99.0 %        | PASS |
|                            | `mix2` | 2.2.1.2 30 %, 2.3.4.4b 70 %                             | 2.2.1.2 **30.4 %**, 2.3.4.4b **68.7 %**                                                          | 99.0 %        | PASS |
|                            | `mix3` | 2.2.1.2 20 %, 2.3.4.4b 30 %, Am-nonGsGD 50 %            | 2.2.1.2 **25.4 %**, 2.3.4.4b **31.7 %**, Am-nonGsGD **40.8 %** (+ 2.2.1.1a 2.1 %)                | 99.5 %        | PASS |
| `H5Nx-cattle` (2024-01-01) | `pure` | 2.3.4.4b.A3 100 %                                       | 2.3.4.4b.A3 **99.8 %**                                                                           | 99.0 %        | PASS |
|                            | `mix2` | 2.2.1.2 30 %, 2.3.4.4b.A3 70 %                          | 2.2.1.2 **32.0 %**, 2.3.4.4b.A3 **67.2 %**                                                       | 99.0 %        | PASS |
|                            | `mix3` | 2.2.1.2 20 %, 2.3.4.4b.A3 30 %, 2.3.4.4b.B3.13.A.1 50 % | 2.2.1.2 **22.6 %**, 2.3.4.4b.A3 **27.3 %**, B3.13 **22.5 %** + B3.6 **22.5 %** (B3.13.A.1 0.0 %) | 99.5 %        | PASS |

For `H5Nx-cattle` the first run (before the lineages of `mix3` were changed to the table above) used `Am-nonGsGD` as the third
lineage: 2.2.1.2 **24.4 %**, 2.3.4.4b.A3 **29.8 %**, Am-nonGsGD **39.2 %** for 20 / 30 / 50 %.

Findings:

- **Accuracy at 500x**: pure samples are recovered at 99.6-100 %. The largest errors on a true lineage are 4.4 points (H1N1
  `close`: 5a 65.6 % for 70 %, three mutations apart) and 4.8 points (H3N2 `mix3`, K 45.2 % for 50 %); the other well separated
  mixtures are within 3 points. That is somewhat worse than RSV (1.3 points) and expected: Freyja has only the HA mutations (H1N1
  lineages carry 49-86 barcode mutations, H3N2 73-121, FLU-B-VIC 23-39, against hundreds for the RSV lineages), and a lineage 3
  mutations from its neighbour is a weak signal. The amplicon-style samples perform like the shotgun ones (1.7 points at most for
  H3N2; 3.3 for H1N1, where 3.4 % went to the parent lineage 5a.2a).
- **Whole-genome BED**: `amp_mix3_wgbed` (the BED with the `seg4` / `seg6` lines) gave the same abundances, depths, 7 amplicons
  (0 dropped) and consensus as `amp_mix3` with the plain BED, in both runs, and each run printed the one warning. The provenance
  records the filtered BED (14 lines, 519 bytes).
- **Reads of other sequences**: the sample with 8-fold more reads of another virus has `Mapped %` 10.3 / 10.2, the same abundances
  and consensus as without them, and QC PASS: nothing in QC or the report treats a low mapped percentage as a problem. (They are not
  "reads of the other seven influenza segments", but for the HA reference they behave the same way: they do not map.)
- **Depth**: HA is short, so reads are few. At 50x (about 290 read pairs) the mixtures are still recovered (H1N1 within 0.8 points,
  H3N2 within 2.6); at 20x they are not reliable (H1N1 3.6 points, H3N2 10.5 points with four spurious lineages of 2.7-6.6 %). Lower
  still (3 further H1N1 samples, same lineages, `--protocol capture`): 15x gave a 12.5-point error and QC PASS (83.9 % of HA at >=
  10x); 10x gave 8.3 points and QC WARN (`low_coverage`, 59.8 % at >= 10x, demixed as `ok`); 5x and 3x were `skipped_low_coverage` (below
  the 25 % floor of the low-confidence band) with QC FAIL. So the existing coverage rules behave as designed at 1.7 kb (they are
  percentages of the reference), but **QC PASS does not mean the abundances are accurate**: below about 50x of HA the lineage
  fractions of mixtures get noisy long before the QC status changes. Consensus genomes at 20x have 14-16 % N (positions below 10x are
  masked).
- **Consensus**: compared with the pseudo-genome of the dominant lineage, no position in any H1N1 sample has a base incompatible
  with it (pure samples are identical outside the Ns; mixtures carry 3-83 IUPAC codes). H3N2 has two positions (631, 639) in the two
  amplicon 20/30/50 samples where the consensus keeps the reference base although K (48 %) has the alternative: both lie in a stretch
  with five K-specific mutations within 13 bases near the 3' end of the reverse reads; the BAM has only reference-allele reads there
  (183 reads, all `126M24S`) although the FASTQ has 212 reads with the K allele, which suggests that `bwa mem` soft-clips the
  divergent 3' ends (not investigated further). The other 9 H3N2 samples and FLU-B-VIC agree fully; the 30/70 `H5Nx` and `H5Nx-cattle` samples have one incompatible position
  (base 12 of the sequence; not investigated).
- **Divergent lineages map badly**: the `Am-nonGsGD` lineage of `H5Nx` / `H5Nx-cattle` differs from the reference at 236 positions
  (13 %). A sample that is half this lineage has `Mapped %` of about 88, the lineage is under-estimated (40.8 % and 39.2 % for 50 %),
  the consensus has 37-39 positions incompatible with it and `unexplained_mutations.tsv` lists 26 (`H5Nx`) and 59 (`H5Nx-cattle`)
  rows. The Goose/Guangdong reference is far from the American non-Gs/GD lineage; 2.3.4.4b (147-178 mutations) behaved. If your
  samples are dominated by a lineage that diverged from the reference by a tenth, expect the same.
- **Indistinguishable lineages**: many `H5Nx-cattle` lineages have identical or almost identical barcode rows (at least six have no
  mutations at all; `B3.13.A.1` is 3 mutations from `B3.13` and `B3.6`), so the 50 % of `B3.13.A.1` was reported as `B3.13` 22.5 %
  plus `B3.6` 22.5 %. The same applies to `H1N1-6` and `H1N1-6B.1`, whose barcodes are identical. This is Freyja's resolution on HA,
  not a pipeline fault.
- `unexplained_mutations.tsv`: 0 rows for H1N1 and FLU-B-VIC, 1 for H3N2 (`mix2`: a synonymous change at HA 533, observed 0.41 for
  0.31 expected), 9 for the `H5Nx-cattle` run above. `annotated_variants.tsv` lists gene `HA` with amino-acid effects (for example
  `HA L8M`). The report builds, its coverage track covers the 1.7 kb reference (3-base bins), and the subtitle reads `Wastewater H1N1
lineage deconvolution` (`H3N2`, ...).
- **Lineage groups with `H5Nx-cattle`**: with a groups file (`H5Nx-2.3.4.4*`, `H5Nx-2.3.4.4b.B3.13*`, ...) the summary
  (`lineage_summary.csv`, `.xlsx`, `lineage_groups_resolved.tsv`) was made from the downloaded `lineages.yml`; the groups matched 65
  and 35 lineages.
- **Negative controls**: three negative-control samples of 12, 36 and 116 mapped reads (1x, 3x and 10x of HA) all FAILed with
  `ntc_reads` because of the 10 % genome-coverage rule (59.7 %, 93.4 % and 98.9 % of HA at >= 1x). See the remark above.

Not run: real influenza reads; a real influenza primer panel (the BED is a toy tiling); the other seven segments; `--freyja_boot`;
`FLU-B-VIC`, `H5Nx` and `H5Nx-cattle` at lower depth or with amplicon-style samples; `tests/e2e.sh` for these pathogens; a comparison
of `latest/` with the newest dated folder of Freyja-barcodes.

## What is supported and what Route B would add

**Supported (PLAN D24, "Route A")**: the HA segment of the five influenza pathogens as a single-sequence reference, with everything the
single-contig pipeline does (QC, amplicon coverage, HA consensus, variant annotation, dropout impact, unexplained mutations, report),
and a whole-genome primer panel tolerated by dropping the lines of the contigs that are not in the reference. Not supported: the
other seven segments (no consensus, QC or variants for them), per-segment reassortant deconvolution (Freyja does not offer it), and
several subtypes in one run (one run per `--pathogen`).

**Route B, all eight segments** (needed only for whole-genome consensus, per-segment QC and variants; not implemented). Every item
below assumes one contig today:

| Part                                                         | Where                                                                                                                                                                         | What would change                                                                                                                                                                                                                                               |
| ------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Start-up check                                               | `validateInputParameters` (single-sequence FASTA)                                                                                                                             | allow a multi-sequence FASTA when a pathogen declares which contig carries the barcodes (HA)                                                                                                                                                                    |
| Freyja                                                       | `modules/local/freyja_variants`                                                                                                                                               | `freyja variants --refname <HA contig>` (its help: "Ref name (for bams with multiple sequences)"); demix, boot and aggregate then stay as they are                                                                                                              |
| Depth QC (`qc_summary.py`, `depth_stats.py`)                 | genome-level `pct_genome_ge_10x`, `mean_depth`                                                                                                                                | per-segment statistics and a rule for the sample status (for example the HA segment must pass; how many of the 8 should)                                                                                                                                        |
| Amplicon coverage and dropout                                | `amplicon_coverage.py` (reads depths per contig already), `dropout_impact.py`                                                                                                 | key amplicons by (segment, amplicon); name handling per scheme; dropout impact only for the barcode segment, because the barcodes only know HA positions                                                                                                        |
| Low-depth mask and consensus                                 | `low_depth_mask` module (first header of the FASTA only), `consensus_stats.py` (`read_single_fasta` raises on more than one record), `consensus_bcftools`, `consensus_concat` | one mask and one consensus record per segment (`<sample>                                                                                                                                                                                                        | <segment>`), consensus statistics per segment |
| Report                                                       | `build_report_data.read_depths` keeps the first contig only                                                                                                                   | a coverage track per segment                                                                                                                                                                                                                                    |
| Unexplained mutations, dropout impact, mutations of interest | positions are bare integers with no contig                                                                                                                                    | restrict to the barcode contig (genome length of that contig for the significance level)                                                                                                                                                                        |
| Variant annotation                                           | `annotate_variants` (`bcftools csq` copes with several contigs)                                                                                                               | a GFF3 with gene -> mRNA -> exon/CDS lines per segment (see the RSV files for the conversion NCBI files need)                                                                                                                                                   |
| Primers                                                      | primer BED (contig is column 1), `amplicon_coverage.py` name parsing                                                                                                          | segment-aware panels already fit the BED format (one contig per segment). Whole-segment RT-PCR (one amplicon per segment) has nothing to trim: treat it with `--protocol capture` or add a per-segment "amplicon" notion, to be decided with the intended panel |
| Primer scheme check                                          | already skipped for any non-SARS-CoV-2 reference                                                                                                                              | nothing                                                                                                                                                                                                                                                         |

Route A was chosen because Freyja's influenza barcodes are HA-only, so the extra segments would have nothing to be compared with; Route B
becomes worthwhile only if whole-genome consensus sequences or per-segment QC are wanted.

## Known limits

- **Lineage groups for measles do not work with today's Freyja data**: the lineage hierarchy (`MEASLESgenome_lineages.yml`) names
  lineages `MEASLES-H1` while the barcodes call them `MEASLESgenome-H1`. A groups file naming either is accepted, but no
  observed lineage matches, so every abundance lands in `other` (the lineage names remain visible in the `ungrouped_lineages`
  column and in `freyja_aggregated.tsv`). Use the aggregate until Freyja's naming is made consistent.
- Freyja's `summarized` column (what `freyja plot` groups by) is `Other` for every measles sample, so its plots have one group.
- The bundled amplicon, QC and report defaults were tuned on SARS-CoV-2 wastewater (D17 QC thresholds: PASS at 60 % of the genome
  > = 10x); review them for your pathogen and panel.
- Only RSV-A, RSV-B, measles and the five influenza HA pathogens have been run, all on simulated reads. For the other
  single-sequence pathogens (MPX, dengue) read the first results critically (reference, primer BED contig names, barcodes date).
- Influenza: HA only, one run per subtype; the negative-control thresholds (`--ntc_max_genome_coverage_pct`, `--ntc_max_mapped_reads`)
  do not scale to a 1.75 kb reference (a few stray reads fail a negative control); abundances of mixtures are noisy below about 50x
  of HA even when QC says PASS; a lineage far from the reference (the `Am-nonGsGD` lineage of `H5Nx`, 13 %) maps poorly and is
  under-estimated; the `H5Nx-cattle` barcodes have 13 columns that do not agree with their own reference and many indistinguishable
  lineages. See "Running influenza" and its validation.
- No RSV primer scheme is bundled; real RSV panels, real reads and the QC defaults (tuned on SARS-CoV-2 wastewater) have not been
  tried together. The RSV barcodes are small (23 and 16 lineages, dataset 2024-03-15) and nested; do not read 0.2 % abundances as
  detections.
- Freyja's own `freyja plot` summaries group by its `summarized` column, which is `Other` for RSV as for measles.
- The dataset (`~/.wasteline/db/freyja_RSVa/`) is downloaded on first use, as for measles; `--barcodes <date>` pins it.
- The influenza datasets (`~/.wasteline/db/freyja_H1N1/` etc.) are downloaded on first use the same way.
