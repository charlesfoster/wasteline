# Other pathogens (`--pathogen`)

wasteline was built for SARS-CoV-2, but Freyja also deconvolves other pathogens (`freyja update --pathogen X`), and
`--pathogen X` now works end to end for the single-sequence ones (PLAN D22). RSV-A and RSV-B (`RSVa`, `RSVb`) and measles are
validated and bundle their reference (PLAN D23). Everything is opt-in: with the default
`--pathogen SARS-CoV-2` nothing changes. Design as built: [wiring/multi_pathogen.md](wiring/multi_pathogen.md).

**Two things decide whether a pathogen works.** The genome must be one sequence (influenza has 8 segments, which is refused),
and the **reference FASTA** must be the one Freyja's barcodes are tied to: the `reference.fasta` in the pathogen's folder of the
[Freyja-barcodes](https://github.com/andersen-lab/Freyja-barcodes) repository (a different reference gives silently wrong
results). For RSV-A, RSV-B and measles it is bundled; for the others you pass `--fasta`. A primer BED for your amplicon panel (or
`--protocol capture`) is the third ingredient; no primer schemes are bundled for other pathogens.

## Support matrix

| `--pathogen`                                       | Status                          | What you need                                                                                                        |
| -------------------------------------------------- | ------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `SARS-CoV-2` (default)                             | supported, validated            | nothing: reference, annotation, primer schemes and lineage groups are bundled                                        |
| `RSVa`                                             | supported, **validated**        | reference and annotation are bundled (`PP109421.1`); `--primer_bed` (or `--protocol capture`). Validation below      |
| `RSVb`                                             | supported, **validated**        | reference and annotation are bundled (`OP975389.1`); `--primer_bed` (or `--protocol capture`). Validation below      |
| `MEASLESgenome`                                    | supported, **validated**        | the reference is bundled (`NC_001498.1`); `--primer_bed` (or `--protocol capture`). Validation below                 |
| `MPX`, `DENV1`, `DENV2`, `DENV3`, `DENV4`          | single sequence, **not run**    | `--fasta` = the `reference.fasta` of the pathogen's Freyja-barcodes folder, `--primer_bed` (or `--protocol capture`) |
| `H5Nx`, `H5Nx-cattle`, `H1N1`, `H3N2`, `FLU-B-VIC` | **unsupported (next priority)** | refused at start-up; see "Influenza: what it needs" below                                                            |
| `MTB`                                              | **unsupported, not wanted**     | a 4.4 Mb bacterial genome, see below                                                                                 |

"Not run" means the code path is the same one that was validated on measles (nothing in it is pathogen-specific), but no
data of that pathogen has been run through it: expect to check the first results yourself. MPX (about 197 kb) in particular
is much larger than the genomes the pipeline has been timed on. "Validated" means simulated reads of known composition (no
real reads of these pathogens were available), see the validation sections.

**Unsupported pathogens stop at start-up** with the reason, before any work is done:

- **Influenza** (`H5Nx`, `H5Nx-cattle`, `H1N1`, `H3N2`, `FLU-B-VIC`): refused because the genome has 8 segments (a multi-contig
  reference) and the depth QC, amplicon coverage, low-depth masking, the consensus and the primer checks assume one contig.
  The start-up message still says that; what influenza needs, including a finding that makes the reason less clear-cut than
  stated there (Freyja's influenza barcodes cover the HA segment only), is in "Influenza: what it needs" below. Any other
  pathogen whose `--fasta` has more than one sequence is refused in the same way.
- **MTB**: _Mycobacterium tuberculosis_ is a 4.4 Mb bacterial genome with feather-only Freyja barcodes. The pipeline is
  built and tested for viral genomes of 10-200 kb (per-position depth tables, QC, consensus and the dropout/mutation tables
  are not sized for 4.4 Mb) and nothing has been validated on bacterial data. Supporting it needs those steps scaled and
  checked on real MTB data.

## What changes for a pathogen other than SARS-CoV-2

| Part of the pipeline                                  | Behaviour                                                                                                                                                                                                                                                                                                                                                    |
| ----------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Reference, annotation                                 | `conf/pathogens/<profile>.config`: `measles.config` bundles `assets/references/NC_001498.1.fasta` (no annotation), `rsva.config` / `rsvb.config` bundle `PP109421.1` / `OP975389.1` with a `.gff3` each. Any other pathogen uses `generic.config`: `--fasta` is **required** (the error says it must be the reference of the barcodes), `--gff` is optional. |
| Primer scheme                                         | The bundled schemes are SARS-CoV-2 only. Without `--primer_bed` (or a per-sample `primer_scheme` path, or `--protocol capture` / `--skip_primer_trimming`) the run stops with a message instead of trimming with the wrong BED. The BED's first column must be the reference's sequence name (checked at start-up).                                          |
| Freyja dataset (`--barcodes`)                         | `latest` works (`freyja update --pathogen X`). `container` does not (the Freyja container bundles SARS-CoV-2 only; refused). A dataset directory of `freyja update --pathogen X`, or a barcodes `.csv`, can be given. See "Freyja datasets" below.                                                                                                           |
| Lineage-group summary (`summary/lineage_summary.csv`) | **Off** unless `--lineage_groups <CSV>` is given (the bundled groups are SARS-CoV-2 Pango lineages); one `INFO` line says so. It also needs the dataset's lineage hierarchy (`<X>_lineages.yml`, only shipped for some pathogens): without it the summary is skipped with one warning. See "Known limits".                                                   |
| Dropout impact, unexplained mutations                 | Work from the dataset's barcodes `.csv` (there is no `lineage_mutations.json`). The unexplained-mutation significance level 0.05 / genome length follows the reference length (a user-set `--unexplained_max_pvalue` is never rescaled).                                                                                                                     |
| Primer scheme check (`--primer_scheme_check`)         | Skipped: it scores the reads against the bundled SARS-CoV-2 schemes, which do not fit another reference. Trimming and per-amplicon coverage use your BED as usual.                                                                                                                                                                                           |
| `--annotate_variants`, amino-acid columns             | Need `--gff` (an error says so). Amino-acid effects are only as good as the annotation you give; the ORF1ab slippage handling is name-based and does nothing without such a gene.                                                                                                                                                                            |
| Report                                                | The subtitle names the pathogen (`Wastewater MEASLESgenome lineage deconvolution`).                                                                                                                                                                                                                                                                          |
| `bin/reanalyse_runs.py`                               | Carries `--pathogen` over; only exercised for SARS-CoV-2.                                                                                                                                                                                                                                                                                                    |

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

## Influenza: what it needs (not implemented)

Influenza is the next priority after RSV. It is refused at start-up today (`unsupportedPathogenReason()`), and nothing below has
been built or run. What was checked while preparing PLAN D23:

**1. Freyja's influenza barcodes cover the HA segment only.** The Freyja-barcodes folders of the five influenza pathogens each hold
a single-sequence `reference.fasta` that is the **HA gene**, and every barcode position lies inside it:

| `--pathogen`  | `reference.fasta` (1 sequence)                                | Length  | Lineages | Largest barcode position |
| ------------- | ------------------------------------------------------------- | ------- | -------- | ------------------------ |
| `H1N1`        | `CY121680.1`, A/California/07/2009 (H1N1) HA                  | 1752 bp | 15       | 1738                     |
| `H3N2`        | `CY163680.1`, A/Wisconsin/67/2005 (H3N2) HA                   | 1737 bp | 32       | 1731                     |
| `FLU-B-VIC`   | `KX058884.1` (HA)                                             | 1885 bp | 8        | 1835                     |
| `H5Nx`        | A/Goose/Guangdong/1/96 (H5N1) HA (no accession in the header) | 1760 bp | 11       | 1749                     |
| `H5Nx-cattle` | the same HA; also has a `lineages.yml` hierarchy              | 1760 bp | 70       | 1749                     |

So Freyja deconvolves influenza from HA mutations alone, whatever the other seven segments do, and per-segment (reassortant)
deconvolution is not available from Freyja. The statement in earlier documentation and in the start-up message that influenza is
refused because "the reference has 8 segments" is therefore only half the story: the **reference Freyja needs is one 1.7-1.9 kb
sequence**. The 8-segment question arises only if the whole genome should be mapped, called and reported.

**2. Route A, HA only (smallest change; the pipeline is already single-contig).** Use `--fasta` = Freyja's HA `reference.fasta`
and remove the by-name refusal for the five pathogens. Reads of the other segments simply do not map (`mapped_pct` is
informational, not a QC rule). Things to build or check, none of them exercised yet:

- one subtype per run: H1 and H3 HA are too different to map to one reference, so a sample of mixed subtypes needs one run per
  `--pathogen` (the samplesheet has no per-sample pathogen);
- the primer BED check refuses a BED with contigs absent from the FASTA; a whole-genome panel BED lists all eight segments, so the
  check must ignore (or the user must filter) the lines of the other segments;
- the defaults tuned on 30 kb genomes (QC thresholds, the amplicon and `--demix_floor_coverage` settings, the significance level,
  which already follows the reference length) need a look on a 1.7 kb reference; the plots and the report assume a genome-length axis
  but nothing breaks at that length that is known;
- real data to validate with (none was available; the simulator `tests/simulate_pathogen_mixture.py` already works on any single
  reference plus barcodes CSV, and was written for this).

**3. Route B, all eight segments (needed only for whole-genome consensus, per-segment QC and variants).** Every item below assumes
one contig today:

| Part                                                         | Where                                                                                                                                                                         | What would change                                                                                                                                                                                                                                               |
| ------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Start-up check                                               | `validateInputParameters` (single-sequence FASTA), `unsupportedPathogenReason()`                                                                                              | allow a multi-sequence FASTA when a pathogen declares which contig carries the barcodes (HA)                                                                                                                                                                    |
| Freyja                                                       | `modules/local/freyja_variants`                                                                                                                                               | `freyja variants --refname <HA contig>` (its help: "Ref name (for bams with multiple sequences)"); demix, boot and aggregate then stay as they are                                                                                                              |
| Depth QC (`qc_summary.py`, `depth_stats.py`)                 | genome-level `pct_genome_ge_10x`, `mean_depth`                                                                                                                                | per-segment statistics and a rule for the sample status (for example the HA segment must pass; how many of the 8 should)                                                                                                                                        |
| Amplicon coverage and dropout                                | `amplicon_coverage.py` (reads depths per contig already), `dropout_impact.py`                                                                                                 | key amplicons by (segment, amplicon); name handling per scheme; dropout impact only for the barcode segment, because the barcodes only know HA positions                                                                                                        |
| Low-depth mask and consensus                                 | `low_depth_mask` module (first header of the FASTA only), `consensus_stats.py` (`read_single_fasta` raises on more than one record), `consensus_bcftools`, `consensus_concat` | one mask and one consensus record per segment (`<sample>                                                                                                                                                                                                        | <segment>`), consensus statistics per segment |
| Report                                                       | `build_report_data.read_depths` keeps the first contig only                                                                                                                   | a coverage track per segment                                                                                                                                                                                                                                    |
| Unexplained mutations, dropout impact, mutations of interest | positions are bare integers with no contig                                                                                                                                    | restrict to the barcode contig (genome length of that contig for the significance level)                                                                                                                                                                        |
| Variant annotation                                           | `annotate_variants` (`bcftools csq` copes with several contigs)                                                                                                               | a GFF3 with gene -> mRNA -> exon/CDS lines per segment (see the RSV files for the conversion NCBI files need)                                                                                                                                                   |
| Primers                                                      | primer BED (contig is column 1), `amplicon_coverage.py` name parsing                                                                                                          | segment-aware panels already fit the BED format (one contig per segment). Whole-segment RT-PCR (one amplicon per segment) has nothing to trim: treat it with `--protocol capture` or add a per-segment "amplicon" notion, to be decided with the intended panel |
| Primer scheme check                                          | already skipped for any non-SARS-CoV-2 reference                                                                                                                              | nothing                                                                                                                                                                                                                                                         |

The decision for the next step is therefore Route A (small, can reuse what exists) versus Route B (a per-contig refactor of the QC,
consensus and report code), and whether the user's influenza use needs more than HA.

## Known limits

- **Lineage groups for measles do not work with today's Freyja data**: the lineage hierarchy (`MEASLESgenome_lineages.yml`) names
  lineages `MEASLES-H1` while the barcodes call them `MEASLESgenome-H1`. A groups file naming either is accepted, but no
  observed lineage matches, so every abundance lands in `other` (the lineage names remain visible in the `ungrouped_lineages`
  column and in `freyja_aggregated.tsv`). Use the aggregate until Freyja's naming is made consistent.
- Freyja's `summarized` column (what `freyja plot` groups by) is `Other` for every measles sample, so its plots have one group.
- The bundled amplicon, QC and report defaults were tuned on SARS-CoV-2 wastewater (D17 QC thresholds: PASS at 60 % of the genome
  > = 10x); review them for your pathogen and panel.
- Only RSV-A, RSV-B and measles have been run, all on simulated reads. For the other single-sequence pathogens (MPX, dengue) read
  the first results critically (reference, primer BED contig names, barcodes date).
- No RSV primer scheme is bundled; real RSV panels, real reads and the QC defaults (tuned on SARS-CoV-2 wastewater) have not been
  tried together. The RSV barcodes are small (23 and 16 lineages, dataset 2024-03-15) and nested; do not read 0.2 % abundances as
  detections.
- Freyja's own `freyja plot` summaries group by its `summarized` column, which is `Other` for RSV as for measles.
- The dataset (`~/.wasteline/db/freyja_RSVa/`) is downloaded on first use, as for measles; `--barcodes <date>` pins it.
