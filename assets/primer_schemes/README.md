# Bundled primer schemes

The primer BEDs that `--primer_scheme <name>` selects. They are per pathogen (PLAN D26): the SARS-CoV-2 schemes sit in this directory
(the default `qiagenv1`, unchanged), the schemes of every other pathogen in a subdirectory named after its profile in
`conf/pathogens/` (`measles/`, `mpx/`). `--pathogen X` offers the schemes of its own directory only. Your own schemes: `--primer_scheme_dir`
(docs/usage.md, "Using your own primer scheme").

A scheme's BED must be on the bundled reference of its pathogen, contig name included: `tests/python/test_boolean_params.py` checks it.

## SARS-CoV-2 (this directory)

`midnight`, `midnightv2`, `articv3`, `articv4.1`, `articv5.3.2`, `eden`, `swift`, `qiagenv1` (default), `xgen`; contig `NC_045512.2`.
Origin: as before PLAN D26 (unchanged).

## Other pathogens: source, version, licence

Source of every file below: the primer scheme repository <https://github.com/quick-lab/primerschemes> (the machine-readable repository behind
<https://labs.primalscheme.com>), commit `8972878b4da79d0a959bdbc3a21b80ceba2a0fe1` (checked 2026-10-22). The files are the scheme's
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
  unmodified copies, not part of the pipeline code. If CC BY-SA is not acceptable for your distribution, delete the `measles/` and `mpx/`
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

## Schemes that exist for our pathogens but are NOT bundled

Looked up in the same repository (all CC BY-SA 4.0), none of them can be used without remapping coordinates:

| `--pathogen`                       | Scheme (primerschemes path)                                       | Why not bundled                                                                                                                                                                                                                             |
| ---------------------------------- | ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `RSVa`, `RSVb`                     | `rsva-rsvb/1000/v1.0.0` (Quick Lab, BCCDC)                        | Its references are `NC_038235.1` (RSV-A, 15,222 bp) and `NC_001781.1` (RSV-B, 15,225 bp), not the sequences of Freyja's barcodes (`PP109421.1`, `OP975389.1`): only 31 % / 45 % of 25-mers are shared, so coordinates differ. Not remapped. |
| `MPX`                              | `artic-inrb-mpox/{400,2500}/*` (ARTIC, INRB)                      | Reference `KJ642613.1_masked` (196,442 bp, clade I): 88 % of 25-mers shared with `NC_063383.1`. Not remapped.                                                                                                                               |
| `MPX`                              | `yale-mpox/2000/*` (Chen, Gagne, Gallagher, Vogels)               | Reference `MT903345` (197,205 bp) is not our record (99.7 % of 25-mers shared; 4 bp shorter). Not remapped. Listed upstream as a community scheme.                                                                                          |
| `H1N1`, `H3N2`                     | `ethz-eawag-h1n1-h3n2/400/v1.0.0` (ETH Zurich, Eawag; wastewater) | Contigs are alignment consensus sequences (`ha.align_consensus`, ...), not the HA records `CY121680.1` / `CY163680.1` (HA consensus: 0 % / 26 % shared 25-mers). Not remapped.                                                              |
| `H1N1`                             | `artic-flu-a/800/v1.0.0` (ARTIC)                                  | Draft. For HA it has three universal-primer lines on `NC_026433.1` (1,701 bp, the RefSeq HA of A/California/07/2009, not `CY121680.1`: 95 % of 25-mers shared, different length). Not remapped.                                             |
| `FLU-B-VIC`, `H5Nx`, `H5Nx-cattle` | none found                                                        |                                                                                                                                                                                                                                             |
| `DENV1-4`                          | `artic-pan-dengue/400`, `artic-dezi-pan-denv/1000`                | Not supported pathogens (no bundled reference); not examined further.                                                                                                                                                                       |

To use one of these anyway, give the BED with `--primer_bed` / `--primer_scheme_dir` after remapping it to the Freyja reference yourself
(`docs/other_pathogens.md`, "Using your own primer scheme"); the pipeline refuses a BED whose contig is not the reference's.
