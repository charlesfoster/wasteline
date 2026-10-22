# Bundled primer schemes

The primer BEDs that `--primer_scheme <name>` selects. They are per pathogen (PLAN D26): the SARS-CoV-2 schemes sit in this directory
(the default `qiagenv1`, unchanged), the schemes of every other pathogen in a subdirectory named after its profile in
`conf/pathogens/` (`measles/`, `mpx/`, `rsva/`, `rsvb/`). `--pathogen X` offers the schemes of its own directory only. Your own schemes: `--primer_scheme_dir`
(docs/usage.md, "Using your own primer scheme").

A scheme's BED must be on the bundled reference of its pathogen, contig name included: `tests/python/test_boolean_params.py` checks it.
Schemes whose coordinates were **lifted** by sequence alignment onto our reference (PLAN D27) have `-lifted` in their name and a `#` header in the
file; see "Adapted schemes" below.

## SARS-CoV-2 (this directory)

`midnight`, `midnightv2`, `articv3`, `articv4.1`, `articv5.3.2`, `eden`, `swift`, `qiagenv1` (default), `xgen`; contig `NC_045512.2`.
Origin: as before PLAN D26 (unchanged).

## Other pathogens: source, version, licence

Source of every file in this section (the two unmodified ones in the table here, the adapted ones under "Adapted schemes" further down): the primer scheme repository <https://github.com/quick-lab/primerschemes> (the machine-readable repository behind
<https://labs.primalscheme.com>), commit `8972878b4da79d0a959bdbc3a21b80ceba2a0fe1` (checked 2026-10-22). The two files in this table are the scheme's
`primer.bed`, **byte for byte** (md5 as the scheme's own `info.json`/README states it; the leading `# artic-bed-version v3.0` comment line
is kept; `ivar trim`, `amplicon_coverage.py` and the contig check skip comment lines). No coordinate and no contig name was changed.

| Bundled name (`--primer_scheme`) | File                               | `--pathogen`    | Scheme (primerschemes path) | Version | Contig        | Primers / amplicons | Status upstream | Licence      | Authors                        |
| -------------------------------- | ---------------------------------- | --------------- | --------------------------- | ------- | ------------- | ------------------- | --------------- | ------------ | ------------------------------ |
| `artic-measles-v1.0.0`           | `measles/artic-measles-v1.0.0.bed` | `MEASLESgenome` | `artic-measles/400/v1.0.0`  | v1.0.0  | `NC_001498.1` | 374 / 47            | validated       | CC BY-SA 4.0 | Chris Kent, Quick Lab (ARTIC)  |
| `bccdc-mpox-v2.3.0`              | `mpx/bccdc-mpox-v2.3.0.bed`        | `MPX`           | `bccdc-mpox/2500/v2.3.0`    | v2.3.0  | `NC_063383.1` | 158 / 79            | validated       | CC BY-SA 4.0 | BCCDC Public Health Laboratory |

- **Licence** (`license` field of each scheme's `info.json`, and the `LICENSE` file of the repository): Creative Commons
  Attribution-ShareAlike 4.0 International (<https://creativecommons.org/licenses/by-sa/4.0/>). It allows redistribution and use, including
  commercial, provided the creators are credited, the licence is linked and changes are indicated, and anything adapted from them is shared
  under the same licence. **These BED files stay under CC BY-SA 4.0 whatever the licence of the rest of the repository is**: they are
  unmodified copies, not part of the pipeline code. If CC BY-SA is not acceptable for your distribution, delete the `measles/`, `mpx/`, `rsva/` and `rsvb/`
  directories; nothing else refers to them (the profile configs then show no bundled schemes).
- **Citation** for `artic-measles-v1.0.0`: <https://doi.org/10.1101/2024.12.20.629611>. For `bccdc-mpox-v2.3.0`: BCCDC Public Health
  Laboratory, protocols.io dx.doi.org/10.17504/protocols.io.n2bvj34nnlk5/v1 (homepage
  <https://github.com/BCCDC-PHL/artic-mpxv2022/tree/main/primer_schemes/mpxv-2022/V2.3>).
- **Is the scheme's reference our reference?**
  - `artic-measles/400/v1.0.0`: the scheme's `reference.fasta` is `NC_001498.1`, sequence-identical to `assets/references/NC_001498.1.fasta`
    (15,894 bp). Nothing to change. (`pan-measles/1000/v1.0.0` has the same reference, but it is a draft: not bundled.)
  - `bccdc-mpox/2500/v2.3.0`: the scheme's `reference.fasta` is `NC_063383.1` (197,209 bp), the same coordinates and bases as
    `assets/references/NC_063383.1.fasta` **except that the 3' copy of the inverted terminal repeat (bases 190,789-197,209, 6,421 bp) is hard-masked
    with `N`** in the scheme's file. Outside that block the two sequences are identical (0 differences); the masked block is a reverse-complement
    copy of positions 1-6,422, which the unmasked reference keeps. The BED needs no coordinate change. One primer (`MPXV_79_RIGHT_1`,
    190,836-190,865) lies in the masked block. Three primers (`MPXV_43_RIGHT_1`, `MPXV_59_RIGHT_1`, `MPXV_61_RIGHT_1`) have sequences 2-3 nt longer
    than their coordinates, as published.
- The measles scheme has alternative primers for other genotypes (`_LEFT_1`, `_RIGHT_2`, ...): about 7 % of them differ from the `NC_001498.1`
  bases by more than 3 bases; this is by design and amplicons are paired across alternatives (`amplicon_coverage.py`).

## Native RSV schemes of `charlesfoster/rsv_pipeline` (PLAN D27)

`rsva/rsva-rsv-pipeline.bed` (`--pathogen RSVa`) and `rsvb/rsvb-rsv-pipeline.bed` (`RSVb`) are the amplicon primer BEDs `data/RSVA.primer.bed` and
`data/RSVB.primer.bed` of the repository `charlesfoster/rsv_pipeline` (owned by the same author as wasteline; commit
`9371277af2ae40e825a142b03e796d5fb03739d7`), written **on exactly our references**: their contigs are the GISAID ids that Freyja's `reference.fasta`
headers use (`EPI_ISL_412866` = hRSV/A/England/397/2017 = our `PP109421.1`; `EPI_ISL_1653999` = hRSV/B/Australia/VIC-RCH056/2019 = our
`OP975389.1`; the sequences are base-for-base identical). **Only the contig column was rewritten** (`EPI_ISL_412866` -> `PP109421.1`,
`EPI_ISL_1653999` -> `OP975389.1`), plus one primer name (below); coordinates, names, scores and strands are the source's. Each file starts with `#` comment
lines saying so. Their md5 before the rewrite: `data/RSVA.primer.bed` `8e8b5d4178f8581768282a6853795b25`, `data/RSVB.primer.bed` `f5b352557812b71b7adf21b5fe2c0de2`.
The licence is not a concern (the author's own repository, no licence statement needed).

| Bundled name (`--primer_scheme`) | `--pathogen` | Contig       | Primers | Amplicons | Amplicon size (min / median / max) | Overlap of neighbours | Genome covered by amplicons                                  |
| -------------------------------- | ------------ | ------------ | ------- | --------- | ---------------------------------- | --------------------- | ------------------------------------------------------------ |
| `rsva-rsv-pipeline`              | `RSVa`       | `PP109421.1` | 100     | 50        | 389 / 407 / 420 bp                 | 82-212 bp             | 99.19 % (59 bp at the 5' end, 64 bp at the 3' end uncovered) |
| `rsvb-rsv-pipeline`              | `RSVb`       | `OP975389.1` | 101     | 50        | 389 / 405 / 418 bp                 | 83-168 bp             | 99.79 % (18 bp and 14 bp uncovered)                          |

- Structure (checked by hand): every amplicon has exactly one LEFT and one RIGHT primer (`RSVA_<n>_LEFT` / `_RIGHT`, `RSVB_<n>_...`; strands match the names), the
  50 amplicons are numbered in genome order and each overlaps its neighbours (so every position inside the scheme is covered by one or two amplicons). Column 5 of the files is a
  primer score (0-60), not a pool number. **There are no primer sequences**, so the sequence identity to the reference cannot be checked; coordinates are what they are.
  `bin/amplicon_coverage.py` derives the 50 amplicons.
- **Quirk, handled and documented**: `RSVB.primer.bed` has 101 lines. The extra line is `RSVB_4_RIGHT_2` at 13691-13721, the same coordinates as `RSVB_45_RIGHT`; its
  name has lost the `5` (it is a second right primer of amplicon 45). Unchanged, `amplicon_coverage.py` adds it to amplicon `RSVB_4` and reports that
  amplicon as 902-13721 (12.8 kb). The bundled file therefore names it **`RSVB_45_RIGHT_2`** (header line in the file). Worth fixing in the source repository too.
- Cross-check against the published scheme: the two panels are different designs (400 bp amplicons in the source repository, 1 kb in `rsva-rsvb/1000/v1.0.0`),
  so their primers do not coincide: of the 211 / 99 lifted published primers, 16 / 21 start or end within 10 bp of a primer end of this panel (66 / 41 within 50 bp);
  all amplicons of the lifted published scheme lie inside the area this panel covers. The source repository's `data/RSV{A,B}_prototype.mapped.bed` hold the same primers mapped to
  `NC_038235.1` / `NC_001781.1`; lifting those with `bin/lift_primer_bed.py` onto our references reproduces the coordinates of 87 of 100 (RSV-A) and 95 of 101 (RSV-B) primers
  exactly (the others are primers whose entries in the mapped files are short partial hits of 8-12 bp elsewhere in the genome). That is an independent check of the tool and of the files.

## Adapted schemes (PLAN D27)

The scheme's reference is a different strain from the reference Freyja's barcodes (and so wasteline) use, so its coordinates are not valid on
ours. RSV strains are collinear with few indels, so the coordinates were **lifted by alignment** with `bin/lift_primer_bed.py` (maintenance
tool, not part of the pipeline; method, limits and options in its docstring): minimap2 2.31 (`-x asm20 -c`) aligns the scheme's reference to
ours, every primer's start and end are mapped through the alignment, and the primer sequence (BED column 7) is compared with OUR reference at the
lifted position (identity = matching bases / bases, ungapped, anchored at the primer's 3' end, IUPAC-aware). A primer that falls in an alignment
gap, changes length by more than 2 bp or has identity below 80 % is **dropped** and listed below; an amplicon that would be left without a LEFT or
a RIGHT primer is dropped whole (none was). Primer names, pools, strands and sequences are unchanged. Every bundled file keeps a `#` header
with source, commit, licence and the change. **These files are adapted material and stay under CC BY-SA 4.0** (same licence as the originals,
credits as below; changes: lifted coordinates or renamed contig, some alternate primers removed).

| Bundled name (`--primer_scheme`)       | `--pathogen` | Scheme (primerschemes path; upstream status)       | Contig                                                 | Primers kept / in scheme | Amplicons | Identity of all the scheme's primers to our reference (>= 95 % / >= 90 % / >= 80 %) | Dropped primers                                                                                                 | Authors                                                   |
| -------------------------------------- | ------------ | -------------------------------------------------- | ------------------------------------------------------ | ------------------------ | --------- | ----------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------- |
| `rsva-1000-v1.0.0-lifted`              | `RSVa`       | `rsva-rsvb/1000/v1.0.0` (tested)                   | lifted `NC_038235.1` -> `PP109421.1`                   | 211 / 215                | 20        | 133 / 182 / 211 of 215                                                              | 4 alternates: `RSVA_8_LEFT_16` (68 %), `RSVA_8_LEFT_9` (78 %), `RSVA_9_LEFT_1` (79 %), `RSVA_10_LEFT_12` (74 %) | Quick-lab, BCCDC                                          |
| `rsvb-1000-v1.0.0-lifted`              | `RSVb`       | `rsva-rsvb/1000/v1.0.0` (tested)                   | lifted `NC_001781.1` -> `OP975389.1`                   | 99 / 100                 | 19        | 81 / 91 / 99 of 100                                                                 | 1 alternate: `RSVB_1_RIGHT_10` (75 %)                                                                           | Quick-lab, BCCDC                                          |
| `yale-mpox-2000-v1.0.0-lifted`         | `MPX`        | `yale-mpox/2000/v1.0.0` (validated)                | lifted `MT903345` -> `NC_063383.1`                     | 326 / 326                | 163       | 326 / 326 / 326 of 326                                                              | none                                                                                                            | N. F. G. Chen, L. Gagne, G. R. Gallagher, C. B. F. Vogels |
| `artic-inrb-mpox-2500-v1.0.0-cladeiib` | `MPX`        | `artic-inrb-mpox/2500/v1.0.0-cladeiib` (validated) | the scheme's own clade IIb coordinates, contig renamed | 147 / 147                | 71        | 146 / 146 / 146 of 147 (native coordinates, nothing lifted)                         | none                                                                                                            | ARTIC network, INRB, Quick Lab                            |

- Source for all four: <https://github.com/quick-lab/primerschemes>, commit `8972878b4da79d0a959bdbc3a21b80ceba2a0fe1` (checked 2026-10-22);
  licence CC BY-SA 4.0 (the scheme's `info.json`). Source `primer.bed` md5: `rsva-rsvb/1000/v1.0.0` `86ade75b8623839e69b0999ef9dbc0a0`,
  `yale-mpox/2000/v1.0.0` `f9fe29b44517167f18f9eaf66a141c44`, `artic-inrb-mpox/2500/v1.0.0-cladeiib` `209d92e5d0354854064f6ef1c2e5e3d8`.
  Citation of the Yale panel: <https://doi.org/10.1371/journal.pbio.3002151>.
- **RSV** (`rsva-rsvb/1000/v1.0.0`, one scheme with two contigs): `NC_038235.1` (RSV-A, 15,222 bp) and `NC_001781.1` (RSV-B, 15,225 bp) share only
  31 % / 45 % of 25-mers with the Freyja references (`PP109421.1`, 15,225 bp; `OP975389.1`, 15,222 bp). Each pair aligns as a single block
  (RSV-A: 5.4 % divergence, indels of 1-6 bp and one 72 bp insertion in ours; RSV-B: 3.6 %, indels of 1-6 bp and one 60 bp insertion). All 315 primers lift
  with **no length change**. The 72 bp / 60 bp insertions are the duplications in the G gene (RSV-A ON1 and RSV-B BA genotypes of the references); they lie
  _inside_ the amplicons `RSVA_7`, `RSVA_8` (+72 bp) and `RSVB_8` (+60 bp), which are that much longer on our references. That is
  the references, not a lifting error. The other amplicon size changes are `RSVA_5` -7 bp and `RSVB_7` -6 bp (deletions inside the amplicon);
  the rest are within 2 bp. As a control, no primer has a position within +-6 bp of its lifted one at which its identity is more than 10 points higher.
  The scheme's alternates (`_LEFT_1` ... `_LEFT_16`) were designed for strains other than the references, so their identity to a single reference
  is lower by design (mean identity of the kept primers 94.6 % RSV-A, 96.4 % RSV-B); the 80 % limit dropped five alternates, each from an
  amplicon that keeps 12-16 other primers. Ungapped identity underestimates primers that differ from the reference by a 1-base indel.
- **Yale mpox**: `MT903345` (197,205 bp) lifts onto `NC_063383.1` (197,209 bp) with all 326 primers at 100 % identity. Six of the 163 amplicons
  change size (`MPXV_1` +3, `MPXV_110` -3, `MPXV_143` and `MPXV_144` -46, `MPXV_148` and `MPXV_149` +97): indels between the two genomes
  inside amplicons, no primer is affected. The amplicons are about 2 kb long.
- **ARTIC/INRB mpox**: `artic-inrb-mpox/2500/v1.0.0` is built on the clade I reference `KJ642613.1_masked` (88 % shared 25-mers with
  `NC_063383.1`), but the repository also holds the same primers with **their own coordinates on clade IIb** (`v1.0.0-cladeiib`, reference
  `NC_063383.1_masked`: our `NC_063383.1` with the 3' ITR (190,789-197,209) masked by `N`, no other difference). That native file is bundled, so
  nothing had to be lifted; only the contig was renamed `NC_063383.1_masked` -> `NC_063383.1` and whitespace removed from the sequence column.
  **Cross-check of the lifting tool**: lifting the clade I file `v1.0.0` onto `NC_063383.1` reproduced the scheme's own clade IIb coordinates exactly for
  145 of the 146 primers it kept (`359ba5b8_3_RIGHT_1` ends at 8,327 instead of 8,325). It dropped 1: `359ba5b8_68_LEFT_1` (47 % ungapped identity: it differs from its
  alternate by a 1-base indel). The scheme's amplicon 1 is named with two prefixes (`2o0fvmwf_1_LEFT_*`, `359ba5b8_1_RIGHT_*`); `amplicon_coverage.py` and the lifting
  tool pair a left-only and a right-only name with the same trailing number, so the scheme has 71 amplicons. Not bundled: `artic-inrb-mpox/2500/v1.0.1`
  (a rebalanced v1.0.0, clade I coordinates only; its lift is identical to v1.0.0's) and `artic-inrb-mpox/400/v1.0.0` (a draft; lifted, it keeps 1,157 of
  1,179 primers and 546 of 552 amplicons with 1,164 of 1,175 mapped primers >= 95 % identical, but 68 amplicons change size by more than 2 bp and it has no clade IIb version).
- The original schemes' own reference sequences are not bundled.

## Schemes that exist for our pathogens but are NOT bundled

Looked up in the same repository (all CC BY-SA 4.0). The influenza panels cannot be used without remapping coordinates (whole-genome influenza will use the scheme's own references later); the rest is explained in the table:

| `--pathogen`                       | Scheme (primerschemes path)                                       | Why not bundled                                                                                                                                                                                                                   |
| ---------------------------------- | ----------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `MPX`                              | `artic-inrb-mpox/400/v1.0.0`, `artic-inrb-mpox/2500/v1.0.1`       | The 400 bp panel is a draft and exists only on the clade I reference (lifted in the validation, not bundled); 2500/v1.0.1 is a rebalanced v1.0.0 without a clade IIb version. The clade IIb version of v1.0.0 is bundled (above). |
| `H1N1`, `H3N2`                     | `ethz-eawag-h1n1-h3n2/400/v1.0.0` (ETH Zurich, Eawag; wastewater) | Contigs are alignment consensus sequences (`ha.align_consensus`, ...), not the HA records `CY121680.1` / `CY163680.1` (HA consensus: 0 % / 26 % shared 25-mers). Not remapped.                                                    |
| `H1N1`                             | `artic-flu-a/800/v1.0.0` (ARTIC)                                  | Draft. For HA it has three universal-primer lines on `NC_026433.1` (1,701 bp, the RefSeq HA of A/California/07/2009, not `CY121680.1`: 95 % of 25-mers shared, different length). Not remapped.                                   |
| `FLU-B-VIC`, `H5Nx`, `H5Nx-cattle` | none found                                                        |                                                                                                                                                                                                                                   |
| `DENV1-4`                          | `artic-pan-dengue/400`, `artic-dezi-pan-denv/1000`                | Not supported pathogens (no bundled reference); not examined further.                                                                                                                                                             |

To use one of these anyway, remap it to the Freyja reference (`bin/lift_primer_bed.py`; `docs/other_pathogens.md`, "Adapting a published scheme to our
reference") and give the BED with `--primer_bed` / `--primer_scheme_dir`; the pipeline refuses a BED whose contig is not the reference's.
