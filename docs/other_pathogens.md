# Other pathogens (`--pathogen`)

wasteline was built for SARS-CoV-2, but Freyja also deconvolves other pathogens (`freyja update --pathogen X`), and
`--pathogen X` now works end to end for the single-sequence ones (PLAN D22). RSV-A and RSV-B (`RSVa`, `RSVb`), mpox (`MPX`), measles and the
five influenza pathogens (`H1N1`, `H3N2`, `FLU-B-VIC`, `H5Nx`, `H5Nx-cattle`; the HA segment only, see "Running influenza")
are validated and bundle their reference (PLAN D23, D24, D26). Primer schemes are bundled for measles, mpox and RSV-A/RSV-B only (PLAN D26, D27, see
"Primer schemes"). Everything is opt-in: with the default
`--pathogen SARS-CoV-2` nothing changes. Design as built: [wiring/multi_pathogen.md](wiring/multi_pathogen.md).

**Two things decide whether a pathogen works.** The reference must be one sequence (a multi-sequence FASTA is refused; Freyja's
influenza barcodes cover only the HA segment, which is one sequence of 1.7-1.9 kb, so influenza qualifies), and it must be the
**reference FASTA** that Freyja's barcodes are tied to: the `reference.fasta` in the pathogen's folder of the
[Freyja-barcodes](https://github.com/andersen-lab/Freyja-barcodes) repository (a different reference gives silently wrong
results). For RSV-A, RSV-B, mpox, measles and influenza it is bundled; for the others you pass `--fasta`. A primer scheme for your amplicon
panel (or `--protocol capture`) is the third ingredient: the published ARTIC/BCCDC/Yale schemes for measles, mpox and RSV are bundled (the RSV and one mpox scheme with coordinates lifted onto our
reference), for every other pathogen you give a BED (`--primer_bed`, `--primer_scheme_dir`), see "Primer schemes".

## Support matrix

| `--pathogen`                       | Status                                         | What you need                                                                                                                                                                                                                                 |
| ---------------------------------- | ---------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SARS-CoV-2` (default)             | supported, validated                           | nothing: reference, annotation, primer schemes and lineage groups are bundled                                                                                                                                                                 |
| `RSVa`                             | supported, **validated**                       | reference and annotation are bundled (`PP109421.1`); `--primer_scheme rsva-rsv-pipeline` or `rsva-1000-v1.0.0-lifted` (bundled), `--primer_bed` or `--protocol capture`. Validation below                                                     |
| `RSVb`                             | supported, **validated**                       | reference and annotation are bundled (`OP975389.1`); `--primer_scheme rsvb-rsv-pipeline` or `rsvb-1000-v1.0.0-lifted` (bundled), `--primer_bed` or `--protocol capture`. Validation below                                                     |
| `MEASLESgenome`                    | supported, **validated**                       | the reference is bundled (`NC_001498.1`); `--primer_scheme artic-measles-v1.0.0` (bundled), `--primer_bed` or `--protocol capture`. Validation below                                                                                          |
| `MPX`                              | supported, **validated**                       | reference and annotation are bundled (`NC_063383.1`, 197 kb); three bundled schemes (`bccdc-mpox-v2.3.0`, `yale-mpox-2000-v1.0.0-lifted`, `artic-inrb-mpox-2500-v1.0.0-cladeiib`), `--primer_bed` or `--protocol capture`. See "Running mpox" |
| `DENV1`, `DENV2`, `DENV3`, `DENV4` | single sequence, **not run**                   | `--fasta` = the `reference.fasta` of the pathogen's Freyja-barcodes folder, `--primer_bed` (or `--protocol capture`)                                                                                                                          |
| `H1N1`                             | supported (HA only), **validated**             | reference and annotation are bundled (`CY121680.1`, HA); `--primer_bed` (or `--protocol capture`). See "Running influenza"                                                                                                                    |
| `H3N2`                             | supported (HA only), **validated**             | bundled `CY163680.1` (HA); `--primer_bed` (or `--protocol capture`)                                                                                                                                                                           |
| `FLU-B-VIC`                        | supported (HA only), **validated (quick run)** | bundled `KX058884.1` (HA); `--primer_bed` (or `--protocol capture`)                                                                                                                                                                           |
| `H5Nx`, `H5Nx-cattle`              | supported (HA only), **validated (quick run)** | bundled `AF144305.1` (HA, one file for both); `--primer_bed` (or `--protocol capture`)                                                                                                                                                        |
| `MTB`                              | **unsupported, not wanted**                    | a 4.4 Mb bacterial genome, see below                                                                                                                                                                                                          |

"Not run" means the code path is the same one that was validated on measles (nothing in it is pathogen-specific), but no
data of that pathogen has been run through it: expect to check the first results yourself. "Validated" means simulated reads of
known composition (no real reads of these pathogens were available), see the validation sections.

**Unsupported pathogens stop at start-up** with the reason, before any work is done:

- **Whole-genome influenza** is not a refused pathogen any more (PLAN D24): the five influenza pathogens run on their HA segment.
  What stays refused is a `--fasta` with more than one sequence (any pathogen), because the depth QC, amplicon coverage,
  low-depth masking, consensus and primer checks assume one contig; see "What is supported and what Route B would add" below.
- **MTB**: _Mycobacterium tuberculosis_ is a 4.4 Mb bacterial genome with feather-only Freyja barcodes. The pipeline is
  built and tested for viral genomes of 10-200 kb (per-position depth tables, QC, consensus and the dropout/mutation tables
  are not sized for 4.4 Mb) and nothing has been validated on bacterial data. Supporting it needs those steps scaled and
  checked on real MTB data.

## What changes for a pathogen other than SARS-CoV-2

| Part of the pipeline                                  | Behaviour                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Reference, annotation                                 | `conf/pathogens/<profile>.config`: `measles.config` bundles `assets/references/NC_001498.1.fasta` (no annotation), `mpx.config` `NC_063383.1` with a `.gff3`, `rsva.config` / `rsvb.config` bundle `PP109421.1` / `OP975389.1` with a `.gff3` each, and the influenza profiles bundle the HA record of each pathogen (`CY121680.1`, `CY163680.1`, `KX058884.1`, `AF144305.1`) with a `.gff3` each. Any other pathogen uses `generic.config`: `--fasta` is **required** (the error says it must be the reference of the barcodes), `--gff` is optional.                                                                                                                                                                                                                                                                                                |
| Primer scheme                                         | Bundled schemes are per pathogen (`assets/primer_schemes/<profile>/`): SARS-CoV-2 keeps its nine and the default `qiagenv1`; measles has one, mpox three, RSV-A and RSV-B two each (one native per RSV type, the published RSV scheme and one mpox scheme lifted, PLAN D27); the others none. No other pathogen has a default scheme: without `--primer_scheme` / `--primer_bed` (or a per-sample `primer_scheme` cell, or `--protocol capture` / `--skip_primer_trimming`) the run stops with a message that lists the schemes bundled for the pathogen, instead of trimming with the wrong BED. The BED's first column must be the reference's sequence name (checked at start-up). Lines of contigs that are not in the reference (a whole-genome influenza panel against the HA reference) are dropped with one warning, see "Running influenza". |
| Freyja dataset (`--barcodes`)                         | `latest` works (`freyja update --pathogen X`). `container` does not (the Freyja container bundles SARS-CoV-2 only; refused). A dataset directory of `freyja update --pathogen X`, or a barcodes `.csv`, can be given. See "Freyja datasets" below.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Lineage-group summary (`summary/lineage_summary.csv`) | **Off** unless `--lineage_groups <CSV>` is given (the bundled groups are SARS-CoV-2 Pango lineages); one `INFO` line says so. It also needs the dataset's lineage hierarchy (`<X>_lineages.yml`, only shipped for some pathogens): without it the summary is skipped with one warning. See "Known limits".                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| Dropout impact, unexplained mutations                 | Work from the dataset's barcodes `.csv` (there is no `lineage_mutations.json`). The unexplained-mutation significance level 0.05 / genome length follows the reference length (a user-set `--unexplained_max_pvalue` is never rescaled).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Primer scheme check (`--primer_scheme_check`)         | Skipped: it scores the reads against the bundled SARS-CoV-2 schemes, which do not fit another reference. Trimming and per-amplicon coverage use your BED as usual.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `--annotate_variants`, amino-acid columns             | Need `--gff` (an error says so). Amino-acid effects are only as good as the annotation you give; the ORF1ab slippage handling is name-based and does nothing without such a gene.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| Report                                                | The subtitle names the pathogen (`Wastewater MEASLESgenome lineage deconvolution`).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `bin/reanalyse_runs.py`                               | Carries `--pathogen` over; only exercised for SARS-CoV-2.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |

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
`freyja update --pathogen MPX` saves an HTML page as `MPX_lineages.yml` (checked today, and again in the mpox validation run below).
wasteline ignores a hierarchy file that is not a YAML list, i.e. MPX has no usable hierarchy until Freyja fixes that.

## Primer schemes (bundled and your own)

PLAN D26. Bundled schemes are **per pathogen**: `--pathogen X` offers the schemes of its own directory,
`assets/primer_schemes/` for SARS-CoV-2 (nine schemes, default `qiagenv1`, unchanged) and `assets/primer_schemes/<profile>/` for the others.
Only the pathogen's own schemes can be named; a name that belongs to another pathogen is refused with a hint. **No pathogen except SARS-CoV-2 has a
default scheme**: without `--primer_scheme`, `--primer_bed`, a per-sample `primer_scheme` cell, `--protocol capture` or
`--skip_primer_trimming` the run stops at start-up and lists what is bundled for the pathogen. Sources, versions, licences and md5 of every
bundled BED: [assets/primer_schemes/README.md](../assets/primer_schemes/README.md).

| `--pathogen`             | Bundled `--primer_scheme`              | Published scheme (<https://github.com/quick-lab/primerschemes>)            | Licence      | Reference of the scheme vs. ours                                                                                        |
| ------------------------ | -------------------------------------- | -------------------------------------------------------------------------- | ------------ | ----------------------------------------------------------------------------------------------------------------------- |
| `SARS-CoV-2`             | nine (see usage.md)                    | as before D26                                                              | as before    | `NC_045512.2`                                                                                                           |
| `MEASLESgenome`          | `artic-measles-v1.0.0`                 | `artic-measles/400/v1.0.0` (ARTIC; validated upstream)                     | CC BY-SA 4.0 | `NC_001498.1`, sequence-identical                                                                                       |
| `MPX`                    | `bccdc-mpox-v2.3.0`                    | `bccdc-mpox/2500/v2.3.0` (BCCDC PHL; validated upstream)                   | CC BY-SA 4.0 | `NC_063383.1`, identical coordinates; the scheme's file masks the 3' ITR with N                                         |
| `MPX`                    | `artic-inrb-mpox-2500-v1.0.0-cladeiib` | `artic-inrb-mpox/2500/v1.0.0-cladeiib` (ARTIC, INRB; validated upstream)   | CC BY-SA 4.0 | the scheme's own clade IIb coordinates on `NC_063383.1_masked` (= our reference with the 3' ITR masked); contig renamed |
| `MPX`                    | `yale-mpox-2000-v1.0.0-lifted`         | `yale-mpox/2000/v1.0.0` (Yale; validated upstream)                         | CC BY-SA 4.0 | `MT903345`: **lifted** (326/326 primers kept, 100 % identity)                                                           |
| `RSVa`                   | `rsva-rsv-pipeline`                    | `charlesfoster/rsv_pipeline` `data/RSVA.primer.bed` (400 bp, 50 amplicons) | author's own | **native**: written on `PP109421.1` (GISAID id `EPI_ISL_412866`); only the contig name rewritten                        |
| `RSVa`                   | `rsva-1000-v1.0.0-lifted`              | `rsva-rsvb/1000/v1.0.0` (Quick-lab, BCCDC; tested upstream)                | CC BY-SA 4.0 | `NC_038235.1`: **lifted** onto `PP109421.1` (211 of 215 primers kept)                                                   |
| `RSVb`                   | `rsvb-rsv-pipeline`                    | `charlesfoster/rsv_pipeline` `data/RSVB.primer.bed` (400 bp, 50 amplicons) | author's own | **native**: written on `OP975389.1` (`EPI_ISL_1653999`); contig rewritten, one primer renamed                           |
| `RSVb`                   | `rsvb-1000-v1.0.0-lifted`              | `rsva-rsvb/1000/v1.0.0` (same scheme file, other contig)                   | CC BY-SA 4.0 | `NC_001781.1`: **lifted** onto `OP975389.1` (99 of 100 primers kept)                                                    |
| `H1N1`, `H3N2`           | none                                   | `ethz-eawag-h1n1-h3n2/400/v1.0.0`, `artic-flu-a/800/v1.0.0` exist          | CC BY-SA 4.0 | consensus alignments / `NC_026433.1`, not our HA records: not remapped                                                  |
| `FLU-B-VIC`, `H5Nx`, ... | none                                   | none found                                                                 | -            | -                                                                                                                       |
| `MPX` (other schemes)    | none                                   | `artic-inrb-mpox/400/v1.0.0` (draft), `/2500/v1.0.1`                       | CC BY-SA 4.0 | clade I `KJ642613.1_masked` only: not bundled, see the README                                                           |

The licence of all of them is Creative Commons Attribution-ShareAlike 4.0 (the repository's `LICENSE` and each scheme's `info.json`). The
BED of measles and the BCCDC mpox panel are unmodified copies; the RSV and Yale ones are adapted (lifted coordinates, some alternates dropped) and the
INRB one has a renamed contig; all are credited in the README and in a `#` header of the file, and stay under CC BY-SA 4.0 (data, not
pipeline code). Where the scheme's reference is not, base for base, our reference, the coordinates were **lifted by alignment** with
`bin/lift_primer_bed.py` (PLAN D27; "Adapting a published scheme to our reference" below), with the figures in the README. The scheme check
(`--primer_scheme_check`, "does the declared scheme match the reads?") scores only against the nine SARS-CoV-2 schemes and stays off for every
other pathogen.

### Using your own primer scheme

Three routes, all checked at start-up (details and examples in [usage.md](usage.md#using-your-own-primer-scheme)):

1. `--primer_bed my_scheme.bed`: one BED for the whole run.
2. A BED path in the samplesheet `primer_scheme` cell (contains a `/` or ends in `.bed`): a different panel per sample.
3. `--primer_scheme_dir schemes/`: a directory of `<name>.bed` files, then `--primer_scheme <name>` and/or `<name>` in the samplesheet cell.
   Your file wins over a bundled scheme of the same name (one warning). An unknown name stops the run and lists the available names; the
   scheme is recorded as `custom: <name>.bed` plus its sha256.

The BED: tab-separated, at least six columns (`contig start end name pool strand`; a seventh with the primer sequence, as in primalscheme
`primer.bed` files, is fine; `#` comment lines are fine). Column 1 must be the **reference's sequence name** (the first word of its FASTA header: `NC_063383.1`
for `MPX`; the error shows the expected name); lines on other contigs are dropped with a warning only when some lines remain on the reference. A BED
for another accession of the same sequence needs only column 1 renamed (`sed`); if the sequences differ, the coordinates do too and the BED
has to be remapped first: the pipeline never does that for you. Primer **names** should be primalscheme-style, `<scheme>_<n>_LEFT_<alt>` /
`<scheme>_<n>_RIGHT_<alt>`, so that `amplicon_coverage.py` pairs them into amplicons (`<scheme>_<n>`; alternates merged); `X_LEFT`/`X_RIGHT` and `X_F`/`X_R`
work too. With any other naming no amplicons can be derived (a warning per primer, an empty amplicon table, `NA` amplicon columns for the
sample, no dropout calls) or, when there is no strand column either, the amplicon step fails for that sample (`NA` again); trimming
still works. `--ivar_offset N` (`ivar trim -x N`, default 0, 5 for the bundled `swift`) widens what counts as a primer-derived read start;
it is run-wide, and only needed when your panel's reads show untrimmed primer bases (the Freyja measles tutorial uses 3).

### Adapting a published scheme to our reference

Published schemes are often designed on another strain than the sequence Freyja's barcodes are tied to (that was the case for the RSV scheme and
for two of the mpox ones). The BED's coordinates are then wrong for our reference, and the pipeline never remaps them. For schemes of
collinear genomes with few indels (RSV, mpox clades) `bin/lift_primer_bed.py` (a maintenance tool, not part of any run; PLAN D27) lifts a BED by
alignment and says how good the result is. It needs `minimap2`:

```bash
CONDA_SUBDIR=osx-64 mamba create -y -p liftenv -c conda-forge -c bioconda minimap2     # CONDA_SUBDIR only on Apple Silicon
PATH=$PWD/liftenv/bin:$PATH bin/lift_primer_bed.py \
    --bed primer.bed --scheme-ref reference.fasta --scheme-contig NC_038235.1 \
    --target-ref assets/references/PP109421.1.fasta \
    --out my_scheme.bed --report my_scheme.report.tsv        # also writes my_scheme.report.amplicons.tsv
```

It aligns the scheme's reference contig to ours, maps every primer's start and end through the alignment and writes a BED on our contig name
(names, pools, strands and the sequence column kept). The per-primer report has the lifted coordinates, the length change and the identity of the
primer sequence to our reference at the lifted position (reverse-complemented for RIGHT primers, anchored at the 3' end, IUPAC-aware); the amplicon
report has the sizes before and after. A primer in an alignment gap, with a length change of more than 2 bp (`--max-length-change`) or an identity below 80 %
(`--min-identity`) is **dropped** and reported; an amplicon left without a LEFT or a RIGHT primer is dropped whole, so the output never contains an
unpaired primer (a dropped amplicon is a coverage hole to expect in the results). Check before you use the BED:

- the identity table (the share of primers at >= 95 / 90 / 80 %): alternates designed for other strains are lower than the primary primers by design;
- the amplicons that changed size by more than 2 bp: real indels between the strains inside an amplicon (the RSV G-gene duplication, +72 / +60 bp) are
  fine, a large shrinking (hundreds of bp) means a deletion in our reference relative to the scheme's and deserves a look;
- that nothing is lifted onto the wrong copy of a repeat (the tool prefers the alignment on the main diagonal);
- then validate on reads: simulated reads of a known mixture (`tests/simulate_pathogen_mixture.py --amplicon-bed`) or real ones.

It is a base-level alignment of two sequences, not a model of primer binding; ungapped identity under-rates a primer that differs by a 1-base indel; the
first and last few dozen bases of the genome are often not aligned (primers there are reported as gaps). If the repository holds a version of the scheme
built on a closer reference (as `artic-inrb-mpox/2500/v1.0.0-cladeiib` does for clade IIb), prefer that over lifting. Licence: a lifted scheme is adapted
material; keep the original's licence and credits (CC BY-SA 4.0 for everything in primerschemes) and say what you changed, as the header of the bundled files does.

## Running RSV-A and RSV-B

`--pathogen RSVa` and `--pathogen RSVb` bundle the reference and its annotation, so the only extra ingredient is the library
description. Two RSV panels are bundled (PLAN D27), neither is a default:

- `--primer_scheme rsva-rsv-pipeline` / `rsvb-rsv-pipeline`: the 400 bp amplicon panel (50 amplicons per type) of the author's
  `charlesfoster/rsv_pipeline` repository. It is **native**: written on exactly Freyja's references (the GISAID ids in their FASTA headers;
  only the contig name was rewritten). It has no primer sequences.
- `--primer_scheme rsva-1000-v1.0.0-lifted` / `rsvb-1000-v1.0.0-lifted`: the published ARTIC/primalscheme panel (`rsva-rsvb/1000/v1.0.0`, Quick-lab
  and BCCDC, 1 kb amplicons, 20 / 19 amplicons, CC BY-SA 4.0). It was designed on other strains (`NC_038235.1`, `NC_001781.1`), so its
  **coordinates were lifted by alignment** onto `PP109421.1` / `OP975389.1` (see "Primer schemes" and "Adapting a published scheme to our reference").

For another panel pass its BED (or put it in a
`--primer_scheme_dir`), whose first column must be the reference sequence name (`PP109421.1` for RSV-A, `OP975389.1` for RSV-B), or say that the
library is not amplicon-based.

```bash
# amplicon library made with the 400 bp panel of rsv_pipeline (native on the Freyja reference)
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results \
    --pathogen RSVa --primer_scheme rsva-rsv-pipeline

# amplicon library made with the published ARTIC/primalscheme panel (lifted onto the Freyja reference)
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results \
    --pathogen RSVa --primer_scheme rsva-1000-v1.0.0-lifted

# amplicon library made with another panel: its primer BED (set --ivar_offset if your panel needs it)
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results \
    --pathogen RSVa --primer_bed my_rsv_a_panel.bed

# probe-capture or shotgun library: no primers to trim, no amplicon dropout analysis
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results \
    --pathogen RSVb --protocol capture
```

Without `--primer_scheme` / `--primer_bed` and without `--protocol capture` / `--skip_primer_trimming` the run stops at start-up with a message
that lists the bundled scheme. `--barcodes` defaults to `latest` (downloaded with `freyja update --pathogen RSVa`, stored in `<db_dir>/freyja_RSVa/<date>T00-00`;
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
`NC_001498.1`). The panel is **bundled** (PLAN D26; CC BY-SA 4.0, byte-identical to the scheme's `primer.bed`, see
[Primer schemes](#primer-schemes-bundled-and-your-own)) as `--primer_scheme artic-measles-v1.0.0`; before D26 you had to download it
(Freyja ships a copy as `docs/data/artic-measles-v1.0.0.bed`) and pass it with `--primer_bed`.

```bash
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results \
    --pathogen MEASLESgenome \
    --primer_scheme artic-measles-v1.0.0 --ivar_offset 3
```

- `--barcodes` defaults to `latest`: Freyja's measles barcodes are downloaded to `<db_dir>/freyja_MEASLESgenome/`. Pin them
  later with `--barcodes 2025-05-01` (the barcode date).
- `--ivar_offset 3` is what the Freyja tutorial uses for this panel (`ivar trim -x 3`); the default for a custom BED is 0.
- Reads must be paired (the samplesheet needs `fastq_2`).
- Without `--primer_scheme` / `--primer_bed` the run stops at start-up (measles has no default scheme) and lists `artic-measles-v1.0.0`.
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

**Run**: `--mode full --pathogen MEASLESgenome --primer_bed artic-measles-v1.0.0.bed --ivar_offset 3` (Freyja's copy of the BED; see the D26 re-run
below for the bundled scheme), `--barcodes latest` (the
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

**Re-run with the bundled scheme (PLAN D26).** The same reads, `--primer_scheme artic-measles-v1.0.0 --ivar_offset 3 --barcodes 2025-05-01`
(the scheme for `H1_pure` named in its samplesheet cell, for the mixture by the run default) gave the same results as the table above to the
digit: `H1_pure` 100.0 % H1, 92.71 % of the genome at >= 10x, 196.34x, 47 amplicons (0 dropped); the mixture D9 77.9 % / H1 22.1 %, 93.38 %,
617.52x, 47 amplicons (0 dropped); QC PASS. A third sample (the mixture's reads again, `primer_scheme` cell `mylab-measles`, a copy of the BED
in a `--primer_scheme_dir`) gave identical numbers, recorded as `custom: mylab-measles.bed` with the same BED sha256. `ivar trim` has no problem with
the `# artic-bed-version v3.0` comment line of the bundled BED. 65 tasks.

Not run: the custom barcode file `measles-wg-barcode.csv` of the tutorial; any pathogen other than measles and SARS-CoV-2 in this section; real
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

## Validation (RSV native schemes of `rsv_pipeline`, PLAN D27)

`rsva-rsv-pipeline` and `rsvb-rsv-pipeline` (`assets/primer_schemes/README.md`, "Native RSV schemes") are on Freyja's references already, so nothing was lifted.
By hand, 2026-10-22; **no real RSV reads**.

- **Structure**: RSV-A 100 primers, RSV-B 101; 50 amplicons each, one LEFT + one RIGHT per amplicon, in genome order; sizes 389-420 bp (median 407) and 389-418 bp (median 405);
  neighbours overlap by 82-212 / 83-168 bp, so 34 % / 33 % of the genome is covered twice; 99.19 % / 99.79 % of the genome is inside an amplicon. One quirk: `RSVB_4_RIGHT_2`
  is `RSVB_45_RIGHT_2` (renamed in the bundled file, see the README). `bin/amplicon_coverage.py` derives 50 amplicons for both. There are no primer sequences, so the
  identity to the reference could not be checked.
- **Against the lifted published scheme**: different designs (400 bp vs 1 kb amplicons); 16 of 211 (RSV-A) and 21 of 99 (RSV-B) lifted published primers start or end within 10 bp of a primer end of
  the native panel (66 / 41 within 50 bp). The area covered by the lifted published scheme (14,753 bp / 14,923 bp) lies completely inside the area covered by the native panel (15,102 / 15,190 bp). Lifting the
  source repository's own `*_prototype.mapped.bed` (the same primers on `NC_038235.1` / `NC_001781.1`) onto our references with `bin/lift_primer_bed.py` gives the native coordinates for
  87 of 100 and 95 of 101 primers (the rest are short partial hits in the mapped files).
- **End-to-end** (docker, `--mode full --barcodes latest` = dataset 2024-03-15, `--annotate_variants`, `--primer_scheme rsva-rsv-pipeline` / `rsvb-rsv-pipeline`; simulated reads from the native BED with
  `tests/simulate_pathogen_mixture.py --amplicon-bed` (whole amplicons as fragments), 2 x 250 bp, about 500x of the genome, 0.2 % errors; with 2 x 150 bp reads the middle of every 400 bp amplicon is uncovered,
  which gave 64 % of the genome at >= 10x and is a property of that simulation, not of the scheme):

| Sample          | Truth                            | Recovered                                                                           | Mean depth | Genome >= 10x | Amplicons (dropped) | Consensus N | Primer-trimmed reads | QC   |
| --------------- | -------------------------------- | ----------------------------------------------------------------------------------- | ---------- | ------------- | ------------------- | ----------- | -------------------- | ---- |
| RSV-A `ua_pure` | A.D.1 100 %                      | A.D.1 **99.8 %** (+ A.D.2.2.1 0.2 %)                                                | 356x       | 98.83 %       | 50 (0)              | 178 (1.2 %) | 100.0 %              | PASS |
| RSV-A `ua_mix3` | A.1 20 %, A.3 30 %, A.D.1 50 %   | A.1 **20.6 %**, A.3 **29.7 %**, A.D.1 **49.3 %** (+ A.3.1.1 0.3 %, A.D.2.2.1 0.1 %) | 356x       | 98.83 %       | 50 (0)              | 178         | 99.8 %               | PASS |
| RSV-B `ub_pure` | B.D.4.1 100 %                    | B.D.4.1 **99.8 %** (+ B.D.E.2 0.2 %)                                                | 354x       | 99.45 %       | 50 (0)              | 84 (0.6 %)  | 100.0 %              | PASS |
| RSV-B `ub_mix3` | B.1 20 %, B.3 30 %, B.D.4.1 50 % | B.1 **20.3 %**, B.3 **29.7 %**, B.D.4.1 **49.3 %** (+ B 0.5 %, B.D.E.2 0.2 %)       | 354x       | 99.45 %       | 50 (0)              | 84          | 99.8 %               | PASS |

The largest error of a true lineage is 0.7 percentage points. All 50 amplicons were found for both, none dropped; Genome >= 10x matches the share inside the amplicons (the uncovered genome ends become N).
Primer-trimmed reads are about 100 % because the simulated fragments are whole amplicons. `unexplained_mutations`: 3 rows in RSV-A `ua_mix3` (positions 4837, 4862, 4871, observed 0.34 for 0.21 expected, in the
overlap of amplicons `RSVA_15`/`RSVA_16`; not investigated), none in RSV-B.

## Validation (lifted primer schemes, PLAN D27)

**What was lifted and how**: see "Adapting a published scheme to our reference" and `assets/primer_schemes/README.md` (sources, commit, licence, the dropped
primers). Everything below was run by hand on 2026-10-22 (minimap2 2.31, `bin/lift_primer_bed.py`); there is no test for it (a maintenance tool). **No real
RSV or mpox reads were available**, so the end-to-end runs use simulated reads, as for D23 and D26.

**(a) Primer identity to our reference at the lifted position** (all primers of the scheme, ungapped, 3'-anchored; the alternates are designed for other strains):

| Scheme (contig lifted -> ours)                                                                              | Primers | >= 95 %               | >= 90 % | >= 80 % | Mean (kept primers) | Dropped (gap / length > 2 bp / identity < 80 % / unpaired amplicon) |
| ----------------------------------------------------------------------------------------------------------- | ------- | --------------------- | ------- | ------- | ------------------- | ------------------------------------------------------------------- |
| `rsva-rsvb/1000/v1.0.0` `NC_038235.1` -> `PP109421.1`                                                       | 215     | 133                   | 182     | 211     | 94.6 %              | 4 (0 / 0 / 4 / 0)                                                   |
| `rsva-rsvb/1000/v1.0.0` `NC_001781.1` -> `OP975389.1`                                                       | 100     | 81                    | 91      | 99      | 96.4 %              | 1 (0 / 0 / 1 / 0)                                                   |
| `yale-mpox/2000/v1.0.0` `MT903345` -> `NC_063383.1`                                                         | 326     | 326                   | 326     | 326     | 100 %               | 0                                                                   |
| `artic-inrb-mpox/2500/v1.0.0` `KJ642613.1_masked` -> `NC_063383.1` (lift only checked, native file bundled) | 147     | 146                   | 146     | 146     | 99.9 %              | 1 (0 / 0 / 1 / 0)                                                   |
| `artic-inrb-mpox/400/v1.0.0` (draft; lifted, not bundled)                                                   | 1,179   | 1,164 of 1,175 mapped | 1,170   | 1,170   | 99.8 %              | 22 (4 gap / 7 length change / 3 identity / 8 unpaired)              |

**(b) Amplicon structure** (from the lifter's amplicon report): every amplicon of the bundled lifted files has at least one LEFT and one RIGHT primer
(RSV-A 20 of 20 amplicons kept with 113 LEFT / 102 RIGHT primers in the scheme and 109 / 102 kept; RSV-B 19 of 19, 44 LEFT / 55 RIGHT kept of 44 / 56; Yale 163 of 163,
163 LEFT + 163 RIGHT). Amplicon sizes, original -> lifted: RSV-A 17 of 20 within 2 bp, `RSVA_5` -7, `RSVA_7` and `RSVA_8` +72 (the G-gene duplication of the
reference, inside those amplicons); RSV-B 17 of 19 within 2 bp, `RSVB_7` -6, `RSVB_8` +60 (G-gene duplication); Yale 157 of 163 within 2 bp, `MPXV_1` +3, `MPXV_110` -3,
`MPXV_143`/`MPXV_144` -46, `MPXV_148`/`MPXV_149` +97 (indels between the two clade IIb genomes). No primer changed length. Genome covered by the union of the
amplicons, original scheme on its own reference vs. lifted scheme on ours: RSV-A 96.5 % -> 96.9 % (the 3' 401 bp and the first 71 bp are in no amplicon, on both references), RSV-B 97.7 % -> 98.0 %, Yale 99.4 % -> 99.5 %.
Independent check of the method: lifting `artic-inrb-mpox/2500/v1.0.0` onto `NC_063383.1` gives the scheme's own clade IIb coordinates for 145 of 146 kept primers (one primer end 2 bp off);
for the RSV lifts, no kept primer has a position within +-6 bp of the lifted one whose identity is more than 10 points higher (shift scan: RSV-A 0 of 211, RSV-B 0 of 99; not run for Yale).

**(c) End-to-end runs of the lifted RSV schemes** (docker, Freyja 2.0.6 container of the pipeline, `--mode full --pathogen RSVa|RSVb --barcodes latest` = the dataset of 2024-03-15,
`--annotate_variants`; simulated reads: `tests/simulate_pathogen_mixture.py --amplicon-bed <the lifted BED> --shear`, 2 x 150 bp, about 500x, 0.2 % errors, 25,000 pairs
over the amplicons; lineages as in the RSV validation above). Truth vs. recovered (`freyja_aggregated.tsv`):

| Sample          | Truth                            | Recovered                                                                                    | Mean depth | Genome >= 10x | Amplicons (dropped) | Consensus N | Primer-trimmed reads | QC   |
| --------------- | -------------------------------- | -------------------------------------------------------------------------------------------- | ---------- | ------------- | ------------------- | ----------- | -------------------- | ---- |
| RSV-A `a_pure`  | A.D.1 100 %                      | A.D.1 **99.8 %** (+ A.D.2.2.1 0.2 %)                                                         | 494x       | 96.53 %       | 20 (0)              | 528 (3.5 %) | 5.6 %                | PASS |
| RSV-A `a_mix3`  | A.1 20 %, A.3 30 %, A.D.1 50 %   | A.1 **20.0 %**, A.3 **28.8 %**, A.D.1 **50.6 %** (+ A.3.1 0.3 %, A.D.2.2.1 0.2 %, A.D 0.1 %) | 490x       | 96.53 %       | 20 (0)              | 528         | 7.9 %                | PASS |
| RSV-A `a_close` | A.D.1 70 %, A.D.1.2 30 %         | A.D.1 **71.0 %**, A.D.1.2 **28.8 %** (+ A.D.2.2.1 0.1 %)                                     | 494x       | 96.53 %       | 20 (0)              | 528         | 5.4 %                | PASS |
| RSV-B `b_pure`  | B.D.4.1 100 %                    | B.D.4.1 **99.8 %** (+ B.D.4.1.1 0.2 %)                                                       | 494x       | 97.66 %       | 19 (0)              | 356 (2.3 %) | 4.8 %                | PASS |
| RSV-B `b_mix3`  | B.1 20 %, B.3 30 %, B.D.4.1 50 % | B.1 **19.9 %**, B.3 **29.0 %**, B.D.4.1 **49.6 %** (+ B 1.0 %, B.D.4 0.3 %, B.D.4.1.1 0.2 %) | 494x       | 97.66 %       | 19 (0)              | 356         | 4.8 %                | PASS |

- The largest error of a true lineage is 1.2 percentage points (RSV-A `a_mix3` A.3 28.8 % for 30 %). Genome >= 10x (96.5 % / 97.7 %) is a little below the share covered by the amplicons (96.9 % / 98.0 %): the 5' 71 / 44 bp and the 3' 401 / 255 bp are in no
  amplicon of the published scheme, and a few more bases at the amplicon ends fall below 10x;
  all 20 / 19 amplicons were found (none dropped), and the dropout tables show no loss. A low `primer_trimmed_pct` is expected with sheared 1 kb amplicons.
- `unexplained_mutations`: 2 rows (RSV-A `a_mix3`, positions 5076 T>C and 5609 A>G at 0.62-0.63; 5076 was flagged in the D23 validation too), none for RSV-B.
- **Routes**: the same `a_mix3` reads were run three ways and gave identical numbers (A.D.1 50.6 %, A.3 28.8 %, A.1 20.0 %, 490x, 20 amplicons, QC PASS): `--primer_scheme
rsva-1000-v1.0.0-lifted` for the run (run 1, three samples), the name in the samplesheet `primer_scheme` cell with no `--primer_scheme` (run 2), and
  `--primer_scheme mylab-rsva --primer_scheme_dir <dir>` with a scratch copy of the BED (run 3; recorded as `custom: mylab-rsva.bed` with its sha256). The RSV-B run named its scheme
  for `b_pure` by `--primer_scheme` and for `b_mix3` by cell. `ivar trim` has no trouble with the five `#` header lines of the bundled BEDs. A pathogen cannot use the other's
  scheme: `--pathogen RSVb --primer_scheme rsva-1000-v1.0.0-lifted` stops at start-up and says it is bundled for another pathogen.

**Mpox (the two new bundled schemes)**, one run (docker, `--mode full --pathogen MPX --barcodes latest` = the dataset of 2025-05-01; the scheme named in each
sample's samplesheet cell; simulated 2 x 150 bp reads with `--amplicon-bed <the bundled BED> --shear`, about 150x; truth IIb.B.1 70 % + IIa 30 %), 49 tasks:

| Sample (scheme)                                     | Recovered                                                                                                                     | Mean depth | Genome >= 10x | Amplicons (dropped) | Consensus N   | QC   |
| --------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- | ---------- | ------------- | ------------------- | ------------- | ---- |
| `yale_mix` (`yale-mpox-2000-v1.0.0-lifted`)         | IIb.B.1 **69.0 %**, IIa **29.9 %** (+ IIb.B.1.7 0.6 %, IIb.B.1.19 0.6 %)                                                      | 148x       | 99.63 %       | 163 (0)             | 731 (0.4 %)   | PASS |
| `inrb_mix` (`artic-inrb-mpox-2500-v1.0.0-cladeiib`) | IIb.B.1 33.8 % + IIb.B.1.3 33.8 % (one mutation apart: 67.6 % together), IIa **30.1 %** (+ IIb.B.1.8 1.6 %, IIb.B.1.10 0.6 %) | 144x       | 95.84 %       | 71 (1)              | 8,199 (4.2 %) | PASS |

- Both files are accepted by the start-up contig check and by `ivar trim`; the Yale lift finds all 163 amplicons. In the INRB run the pipeline pairs amplicon 1
  (`2o0fvmwf_1_LEFT_*` with `359ba5b8_1_RIGHT_*`) and reports it as dropped (mean depth 3x): that is the **simulator's** limit (it pairs primers by the exact name prefix, so it
  made no reads for amplicon 1), not a fault of the scheme. The 4.2 % N are consistent with the scheme's own coverage (its amplicons span 95.2 % of the genome, nothing is amplified
  beyond about 188 kb, where the masked 3' repeat begins) plus amplicon 1; not investigated further.
- Not run for mpox: real reads, the BCCDC scheme again, the clade I schemes.

## Running mpox

`--pathogen MPX` bundles the reference and its annotation, and a primer scheme (PLAN D26):

```bash
# amplicon library made with the BCCDC mpox panel (2.5 kb amplicons)
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv --outdir results \
    --pathogen MPX --primer_scheme bccdc-mpox-v2.3.0

# any other amplicon panel: its BED (contig NC_063383.1), or a directory of your schemes
nextflow run charlesfoster/wasteline -profile docker --input samplesheet.csv --outdir results \
    --pathogen MPX --primer_bed my_mpox_panel.bed

# probe-capture or shotgun library
nextflow run charlesfoster/wasteline -profile docker --input samplesheet.csv --outdir results \
    --pathogen MPX --protocol capture --annotate_variants
```

`--barcodes` defaults to `latest` (`freyja update --pathogen MPX`; the dataset of 2025-05-01 today: 49 lineages, 1319 columns, stored in
`<db_dir>/freyja_MPX/2025-05-01T00-00`). The lineages are named `MPX-IIb.B.1`, `MPX-Ib`, ... `--barcodes container` is refused.

- **Reference**: GenBank/RefSeq `NC_063383.1` (MPXV-M5312_HM12_Rivers, clade IIb, 197,209 bp), sequence-identical to the `reference.fasta` of
  `MPX/latest` in Freyja-barcodes, whose header is `MPXV-M5312_HM12_Rivers` (the NCBI record was found by that name and compared base for
  base; the bundled FASTA is the NCBI record, whose first header word `NC_063383.1` is the contig name). Checked for PLAN D26: **all 1319
  SNV columns** of the barcodes (49 lineages) have the stated reference base at the stated 1-based position: 100 %, no exceptions.
- **Inverted terminal repeats**: the genome has identical repeats of 6,422 bp at both ends (positions 1-6,422 and the reverse complement at
  190,788-197,209; checked: `NC_063383.1[1:6422]` equals the reverse complement of `[190788:]` exactly). The barcodes treat them as two places: 47 SNV columns lie in the
  left repeat and 48 in the right one, 45 of them are mirror pairs (complementary bases at `p` and `197210 - p`, the same lineage pattern),
  and two left-repeat columns (`T600G`, `C786T`) have no mirror. So a mutation in a repeat counts twice for Freyja, and reads from either
  repeat map to both copies (`bwa` places them arbitrarily). Consequences seen in the validation: the reads of the two copies are shared between them (in the shotgun
  simulation both repeats have the depth of the rest of the genome, 487x and 500x against 494x; with the amplicon scheme, whose amplicons cover the 5' repeat and only the first 77 bases of the 3' one,
  both copies have 294x and 308x against 505x), and a mutation that only one copy carries in the
  barcodes (`T600G`) appears at the mirror position (`A196610C`) in the other copy at the share of reads mapped there, where it is reported
  as an unexplained mutation (below). Real genomes carry the same mutation in both copies; the simulator applies the barcode columns as they are.
- **Annotation**: the NCBI RefSeq GFF3 (179 genes named `OPG001`, `OPG002`, ...), converted to gene -> mRNA -> exon + CDS as for RSV (the NCBI file
  has only gene and CDS lines, which `bcftools csq` cannot use). The four genes inside the repeats occur twice; the copy near the 3' end is named
  with the suffix `_R` (`OPG003` and `OPG003_R`). Sequence features (25 `sequence_feature` lines) and the repeat regions of the NCBI file are
  not kept. `--annotate_variants` gives gene and amino-acid effects (for example `OPG002 N313T`).
- **Lineage summary**: off (no bundled lineage groups). Freyja's `pathogen_config.yml` points `MPX` at a GitHub web page instead of the raw
  `lineages.yml`, so `MPX_lineages.yml` in the dataset is an HTML page (confirmed in the validation run: `freyja_data.json` has
  `lineages_yml_sha256: null`); the pipeline ignores it, and `--lineage_groups` is skipped with a warning. Dropout impact and unexplained
  mutations work from the barcodes CSV.
- **The 197 kb genome**: the genome is 6.6 times the size of SARS-CoV-2. What was checked: the QC thresholds are percentages of the
  reference (no change needed); the per-position depth tables and the coverage track (600 bins of 329 bp) follow the reference length; the
  unexplained-mutation significance level is rescaled to 0.05 / 197,209; `bcftools csq` annotates the 179 genes; the report axis reads the
  genome length from the data. **One thing did break**: the nf-core `BWA_INDEX` module asks for 7 bytes of memory per reference base, which is
  1.4 MB here, and Docker refuses `--memory 1m` ("Minimum memory limit allowed is 6MB"); a reference of roughly 150-860 kb fails at
  the first process. `conf/modules.config` now sets a floor of 256 MB for that process (SARS-CoV-2 and the smaller genomes are unaffected in
  practice). The negative-control rule (`--ntc_max_mapped_reads` 1000, `--ntc_max_genome_coverage_pct` 10 %) is more tolerant in reads than for SARS-CoV-2
  (10 % of the genome is about 130 mapped reads of 150 bp, against 20 for SARS-CoV-2), see the validation.
- **Primer schemes** (three bundled, CC BY-SA 4.0, see "Primer schemes"): `bccdc-mpox-v2.3.0`, `yale-mpox-2000-v1.0.0-lifted` (163 amplicons of about
  2 kb, lifted from `MT903345`) and `artic-inrb-mpox-2500-v1.0.0-cladeiib` (71 amplicons of about 2.9 kb, the scheme's own clade IIb
  coordinates). `bccdc-mpox-v2.3.0` has 79 amplicons of about 2.5 kb in two pools (158 primers). Libraries of such long amplicons are
  fragmented (tagmentation) before sequencing, so reads start all over the amplicon, only few begin at a primer and `primer_trimmed_pct` is low
  (1.4 % in the simulation); that is expected and no QC rule fires on it with the default `--qc_min_primer_trimmed_pct 0`.

## Validation (mpox, simulated reads)

**Why simulated**: no real mpox reads were available. `tests/simulate_pathogen_mixture.py` (a development tool, not part of the test suite)
builds one pseudo-genome per lineage (the reference plus every mutation of the lineage's row of the dataset's barcodes CSV; all
SNV, none skipped) and draws read pairs (2 x 150 bp, mean depth about 500x, 0.2 % substitution errors, fragment length 350 +/- 40) in fixed
proportions. The amplicon-style samples draw each fragment from a random amplicon of the bundled BCCDC scheme (new option `--shear`: a random
350 bp piece of the 2.5 kb amplicon, as in a tagmented library; without it a 2.5 kb amplicon would only be covered at both ends). As for RSV
and influenza this tests the plumbing and the arithmetic, **not** how Freyja performs on real mpox diversity (private mutations, indels,
recombinants, deletions are absent from the pseudo-genomes).

**Run** (docker, Freyja 2.0.6, `--mode full --pathogen MPX --protocol capture --barcodes latest --annotate_variants`; the amplicon samples have
`protocol` = `amplicon` and `primer_scheme` = `bccdc-mpox-v2.3.0` in their samplesheet cells): eight samples, 139 tasks, 3 min 12 s on a laptop
(the slowest task, `FREYJA_VARIANTS`, 51 s and 470 MB). Lineages are shown without the `MPX-` prefix; only abundances of 0.3 % and more are listed
for the minor lineages, the rest is summed.

| Sample                   | Truth                                           | Recovered (`freyja_aggregated.tsv`)                                                      | Mean depth | Genome >= 10x   | Amplicons (dropped) | Consensus N | QC                     |
| ------------------------ | ----------------------------------------------- | ---------------------------------------------------------------------------------------- | ---------- | --------------- | ------------------- | ----------- | ---------------------- |
| `pure_B1` (shotgun)      | IIb.B.1 100 %                                   | IIb.B.1 **99.4 %** (+ 3 lineages of 0.2 %)                                               | 494x       | 100.00 %        | -                   | 0.00 %      | PASS                   |
| `mix2` (shotgun)         | IIb.A 30 %, IIb.B.1 70 %                        | IIb.B.1 **69.3 %**, IIb.A **29.5 %** (+ A.1 0.3 %, 5 lineages of 0.2 %)                  | 493x       | 100.00 %        | -                   | 0.00 %      | PASS                   |
| `mix3` (shotgun)         | Ib 20 %, IIa 30 %, IIb.B.1 50 %                 | IIb.B.1 **49.1 %**, IIa **29.4 %**, Ib **19.3 %** (+ I 1.1 %, IIb.A.1 0.6 %, 3 of 0.2 %) | 493x       | 99.99 %         | -                   | 0.01 %      | PASS                   |
| `close` (shotgun)        | IIb.B.1 70 %, IIb.F.1 30 % (12 mutations apart) | IIb.B.1 **69.2 %**, IIb.F.1 **29.5 %** (+ 5 lineages of 0.2 %)                           | 494x       | 99.99 %         | -                   | 0.01 %      | PASS                   |
| `amp_pure_B1` (amplicon) | IIb.B.1 100 %                                   | IIb.B.1 **99.3 %** (+ 4 lineages of 0.2 %)                                               | 496x       | 99.92 %         | 79 (0)              | 0.08 %      | PASS                   |
| `amp_mix3` (amplicon)    | Ib 20 %, IIa 30 %, IIb.B.1 50 %                 | IIb.B.1 **48.1 %**, IIa **29.9 %**, Ib **19.3 %** (+ I 0.5 %, 9 lineages of 0.2-0.4 %)   | 496x       | 99.92 %         | 79 (0)              | 0.08 %      | PASS                   |
| `neg5` (negative)        | 5 read pairs of `pure_B1`                       | not demixed (`skipped_low_coverage`)                                                     | 0.01x      | 0.73 % at >=1x  | -                   | -           | PASS                   |
| `neg100` (negative)      | 100 read pairs of `pure_B1`                     | not demixed (`skipped_low_coverage`)                                                     | 0.15x      | 13.80 % at >=1x | -                   | -           | **FAIL** (`ntc_reads`) |

- **Abundances**: the largest error on a true lineage is 1.9 percentage points (`amp_mix3`, IIb.B.1 48.1 % for 50 %); most are under 1. The
  pairs of lineages 12 mutations apart (IIb.B.1 / IIb.F.1) are separated correctly. Lineages one mutation apart (IIb.B.1 and its children
  IIb.B.1.1 ... , which differ from it by a single mutation each) cannot be told apart by their barcodes: the 0.2 % abundances on neighbouring
  `IIb.B.1.*` lineages in every sample are that, not contamination (the nested tree again, as for RSV).
- **Consensus genomes**: 197,209 bp; compared with the pseudo-genome of the dominant lineage, no position in any of the six samples has a base
  incompatible with it, and the pure samples are identical to it outside the Ns. The Ns are the ends of the genome that the simulated fragments
  hardly reach: 8-15 Ns for the shotgun samples (the first 4-9 and the last 3-6 bases), 151 for the amplicon ones (positions 1-83 and
  197,142-197,209, outside the scheme's first and last amplicon). Mixed samples carry 62 (`mix2`), 417 (`mix3`), 12
  (`close`) and 421 (`amp_mix3`) IUPAC codes (the sites where the mixed lineages differ).
- **Amplicon samples**: 79 amplicons, none dropped; 99.92 % of the genome at >= 10x, including the 3' repeat that the scheme does not cover (its
  reads map to the 5' copy and the other way round). `primer_trimmed_pct` 1.43 %.
- **Unexplained mutations**: 4 rows (`mix3`, `amp_mix3`), two sites each: `A196610C` (observed 0.09-0.19, expected 0, in the 3' repeat: the
  mirror of the barcode column `T600G` that only the 5' copy has, see "Inverted terminal repeats"; the report flags it as a recurrent mutation
  "possible novel lineage", which it is not) and `C179082T` (0.10-0.12, six bases from the clade-I SNVs `T179072A` and `G179076T`; presumably a
  misaligned base next to a SNV cluster of the clade-I reads; not investigated further). Neither affects the abundances.
- **Negative controls**: 5 read pairs (0.01x, 0.7 % of the genome) pass; 100 pairs (200 mapped reads, 0.15x, 13.8 % of the genome at >= 1x) fail
  with `ntc_reads` because of the 10 % genome-coverage rule, as designed.
- `annotated_variants.tsv`: 59-69 annotated variants per clade IIb sample and 273-274 for the clade-mixed ones, with gene and amino-acid effect
  (388 missense, 403 synonymous, 2 intergenic over all samples). The report subtitle reads `Wastewater MPX lineage deconvolution`.
- The Freyja dataset the run downloaded (`MPX_barcodes.csv`, 2025-05-01) is content-identical to `MPX/latest/barcode.csv` that the simulation used.

Not run: real mpox reads; a real mpox amplicon library; lineages with indels (the barcodes list SNVs only); the clade I mpox schemes (see "Primer
schemes" and the PLAN D27 validation below); `--freyja_boot`; a comparison of `latest/` with the newest dated folder of Freyja-barcodes; `tests/e2e.sh` for
this pathogen.

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
- Only RSV-A, RSV-B, mpox, measles and the five influenza HA pathogens have been run, all on simulated reads (the measles ones are the
  Bygul-simulated reads of the Freyja tutorial). For the other single-sequence pathogens (dengue) read the first results critically
  (reference, primer BED contig names, barcodes date).
- Mpox: the 6,422 bp inverted terminal repeats are in the reference twice and in the barcodes as separate columns; reads from either copy map to both,
  so a mutation of a repeat is split between two places and can show up as an "unexplained mutation" at the mirror position (see "Running mpox").
- Influenza: HA only, one run per subtype; the negative-control thresholds (`--ntc_max_genome_coverage_pct`, `--ntc_max_mapped_reads`)
  do not scale to a 1.75 kb reference (a few stray reads fail a negative control); abundances of mixtures are noisy below about 50x
  of HA even when QC says PASS; a lineage far from the reference (the `Am-nonGsGD` lineage of `H5Nx`, 13 %) maps poorly and is
  under-estimated; the `H5Nx-cattle` barcodes have 13 columns that do not agree with their own reference and many indistinguishable
  lineages. See "Running influenza" and its validation.
- The bundled RSV schemes: `rsv{a,b}-rsv-pipeline` is native on Freyja's references but has no primer sequences (identity cannot be checked); the published one is
  **lifted** by alignment (not designed on them; 4 of 215 / 1 of 100 alternate primers dropped). Both were checked on simulated reads only. Real RSV
  panels, real reads and the QC defaults (tuned on SARS-CoV-2 wastewater) have not been tried together. The RSV barcodes are small (23 and 16 lineages, dataset 2024-03-15) and nested; do not read 0.2 % abundances as
  detections.
- Freyja's own `freyja plot` summaries group by its `summarized` column, which is `Other` for RSV as for measles.
- The dataset (`~/.wasteline/db/freyja_RSVa/`) is downloaded on first use, as for measles; `--barcodes <date>` pins it.
- The influenza datasets (`~/.wasteline/db/freyja_H1N1/` etc.) are downloaded on first use the same way.
