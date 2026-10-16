# Aligner evaluation: `bwa mem` (default) vs `minibwa`

Date of the research: 2026-10-16. Decision D2 keeps `bwa mem`; this page records what was found when the alternative was
looked at, and what is in the repository to test it.

## Verdict

- **`bwa mem` stays the default.** Nothing here argues for changing it.
- **minibwa exists, is real, and is installable from an exact-pinned container.** `MINIBWA_MEM`
  (`modules/local/minibwa_mem/`) and `tests/compare_aligners.sh` are ready for an opt-in comparison (`--aligner minibwa`; design notes in [wiring/small_options.md](wiring/small_options.md)).
- **It is a young tool.** The first release is four months old and it is at 0.x with a preprint only. That is reason enough
  to run it as an experiment, not to make it the default of a pipeline whose output is used for surveillance.
- **The speed advantage does not matter for this pipeline.** Mapping a 30 kb genome takes well under a second; the pipeline
  spends its time elsewhere (container start-up, `freyja variants`/`demix`, bootstrap). The only possible reason to switch
  would be better results on real data, and that has not been shown for amplicon wastewater reads.
- **Recommended path:** run `tests/compare_aligners.sh` on the synthetic set (CI, cheap), then on a few real
  wastewater runs with both aligners and compare abundances and consensus genomes. Reconsider the default only when minibwa
  is at 1.x (or has a peer-reviewed publication) and the comparison on real data is clean.

## What minibwa is

| Item         | Finding                                                                                                                                                                                                                                                                                                         |
| ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| What         | A short-read (and accurate long-read) aligner described by its author as the "successor of bwa-mem": bwa-mem's FM-index and variable-length seeding, minimap2's chaining and ksw2 base alignment, batched/prefetched seeding like bwa-mem2.                                                                     |
| Maintainers  | Heng Li (`lh3`, author of bwa, minimap2, samtools/htslib; 408 commits) and Nils Homer (`nh13`, 17 commits). Repository <https://github.com/lh3/minibwa> (created 2025-10-13, last push 2026-10-01, 352 stars, 11 open issues).                                                                                  |
| Publication  | Preprint: Li H, Homer N. "Fast genomic read alignment with minibwa", arXiv:2606.15357, 13 June 2026. No journal or peer review recorded on the arXiv page. Claims: about 4x faster than BWA-MEM and over 2x faster than BWA-MEM2 at "comparable accuracy" (benchmark data are in the paper, not examined here). |
| Releases     | v0.1 (2026-06-15), v0.2 (06-18), v0.3 (06-28), v0.4 (07-12), v0.5 (07-26), v0.6 (07-30), v0.7 (2026-08-06, "r421"). Seven releases in under two months; all 0.x. The index format and options are still changing between versions.                                                                              |
| Packaging    | Bioconda `minibwa` 0.1 to 0.7 (linux-64, linux-aarch64, osx-64, osx-arm64 for >= 0.3/0.4, so `-profile conda` works on Apple Silicon). BioContainers `quay.io/biocontainers/minibwa:0.7--h118bc1c_0` (amd64 only).                                                                                              |
| Requirements | SSE4.2 (x86) or NEON (ARM); zlib. Index uses about 18 bytes of RAM per genome base (540 MB is not a concern for 30 kb).                                                                                                                                                                                         |
| Known limits | From the README: no noisy long reads or spliced RNA-seq, no undirectional bisulfite, no alternate haplotypes. "Output is not identical to bwa-mem" and many bwa-mem options of the legacy `mem` interface are ignored. Breaking changes to bwa-mem behaviour are part of the design.                            |

### CLI used by the module

```
minibwa index ref.fa                       # writes ref.fa.l2b and ref.fa.mbw (not bwa's .amb/.ann/.bwt/...)
minibwa map -t N ref.fa R1.fq.gz R2.fq.gz  # SAM to stdout; paired-end by giving two files; default preset adaptive
minibwa mem ...                            # legacy bwa-mem style CLI, "not recommended" by the author (not used)
minibwa version                            # 0.7-r421
```

Useful `map` options (put them in `ext.args` for `MINIBWA_MEM`): `-x sr` short-read preset instead of the adaptive one
(adaptive treats reads shorter than 325 bp as short anyway), `-u` drop unmapped reads, `-R '@RG\tID:..'` read group, `-k`
minimum seed length (19), `-A/-B/-O/-E` scoring (bwa's 1/4/6/1 differ: minibwa defaults to 2/8/12,23/2,1), `-Y`
soft-clip supplementary alignments. `minibwa map` writes unmapped reads by default, as bwa mem does, so the pipeline's
`samtools view -F 4` step works unchanged.

## Licence

- The repository README states **MIT**, plus incorporated code: libsais (Apache 2), mimalloc (MIT).
- The default `make` on `master` also compiles `bwtgen.c` (GPL 2, low-memory BWT construction): "if you compile this file,
  which is the default, the resulting binary will be GPL'd". `make gpl=0` or the `Apache2` branch avoid it. The Bioconda
  recipe declares **GPL-2.0-or-later**, so the binary in the container should be assumed GPL. GitHub's licence detector
  reports "other / NOASSERTION" because of this mix.
- This is no problem for wasteline: it runs minibwa as a separate program in a container and links nothing, exactly as it
  does with bwa (GPL-3.0). Only someone redistributing a modified minibwa binary would need to care.

## Container

| Check                                          | Result                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Single-tool image                              | `quay.io/biocontainers/minibwa:0.7--h118bc1c_0`, digest `sha256:14b7bd38...`, verified with `https://quay.io/api/v1/repository/biocontainers/minibwa/tag/` (tags 0.1 to 0.7). Contains minibwa only, **no samtools**.                                                                                                                                                                                                                                                                                                                                                          |
| Image with minibwa and samtools (the one used) | `quay.io/biocontainers/mulled-v2-842faba7afc332e008f7d8b210ce965379274249:371090825f1a2952ebfb0c8b2ea8dd7a86651fbe-0`, digest `sha256:d7dbaa62...`, built 2026-09-01, an existing BioContainers mulled image. Verified: the repository name is the mulled hash of `minibwa` + `samtools`, the tag the hash of `0.7` + `1.24`, and `docker run` prints `minibwa 0.7-r421` and `samtools 1.24` (same samtools as the other modules). Nothing was built or published for this. Singularity: `https://depot.galaxyproject.org/singularity/mulled-v2-842f...:3710...-0` (HTTP 200). |
| Architecture                                   | amd64 only (single-platform manifest). On Apple Silicon Docker runs it under emulation, as it already does for some other images; the module tests pass that way. Production is Linux amd64.                                                                                                                                                                                                                                                                                                                                                                                   |
| Conda                                          | `modules/local/minibwa_mem/environment.yml` pins `minibwa=0.7`, `samtools=1.24`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |

`bwa-mem2` (2.3) and `bwa-mem3` (0.15.0) containers also exist in BioContainers; they were not evaluated.

## What was tested

1. **Module tests** (`nf-test test modules/local/minibwa_mem --profile docker`): stub, and real runs of the synthetic MIX_A
   pair against the bundled reference: 12000 of 12000 reads mapped, all in proper pairs, coordinate-sorted BAM, versions
   reported.
2. **Direct comparison, synthetic reads** (untrimmed `tests/data/fastq` pairs against `NC_045512.fasta`, each aligner in its
   container, 4 threads; not the full pipeline):

   | Sample | Reads | Mapped bwa / minibwa | Proper pairs bwa / minibwa | MAPQ 0 bwa / minibwa | Primary alignments identical (name, mate, position, CIGAR) |
   | ------ | ----: | -------------------- | -------------------------- | -------------------- | ---------------------------------------------------------- |
   | MIX_A  | 12000 | 12000 / 12000        | 12000 / 12000              | 0 / 0                | 11998 of 12000 (99.98%)                                    |
   | MIX_B  | 14000 | 14000 / 14000        | 14000 / 14000              | 0 / 0                | 13994 of 14000 (99.96%)                                    |
   | PURE_C | 16000 | 16000 / 16000        | 16000 / 16000              | 0 / 0                | 15991 of 16000 (99.94%)                                    |

   minibwa aligned each sample (index excluded) in about 0.15 s wall time under amd64 emulation. bwa was not timed on the
   same footing, so no speed ratio is claimed. These reads are simulated from the SARS-CoV-2 reference and are close to
   error free: they say the two aligners agree on easy data, **not** how they differ on real amplicon reads with
   sequencing errors, primer-derived clipping, deletions in variants of concern, or low-complexity ends.

3. **End to end through the pipeline** (`tests/compare_aligners.sh`, run 2026-10-16 with the wired `--aligner` option, Docker,
   `tests/data`, `--primer_scheme midnight`): the whole pipeline with `--aligner bwa` and with `--aligner minibwa`.

   | sample        | max abs. abundance difference | mapped reads bwa / minibwa | mean depth bwa / minibwa |
   | ------------- | ----------------------------- | -------------------------- | ------------------------ |
   | LOWCOV_D      | not in the aggregate (both)   | 119 / 119                  | 0.55 / 0.55              |
   | MIX_A         | 0.0000                        | 11998 / 11998              | 54.86 / 54.86            |
   | MIX_B         | 0.0000                        | 13998 / 13998              | 64.17 / 64.17            |
   | NEG_CONTAM    | not in the aggregate (both)   | 5998 / 5998                | 27.39 / 27.39            |
   | NEG_E         | not in the aggregate (both)   | 5 / 5                      | 0.02 / 0.02              |
   | PURE_C        | 0.0000                        | 15998 / 15998              | 73.23 / 73.23            |
   | Site_1_week_2 | 0.0000                        | 12998 / 12998              | 59.55 / 59.55            |

   Identical on this data (the gate is 0.02): as expected for close-to-error-free simulated reads. It says nothing about real reads.

## How to compare

```bash
tests/compare_aligners.sh /scratch/aligners              # docker profile, tolerance 0.02
TOLERANCE=0.01 tests/compare_aligners.sh /scratch/aligners -profile docker,locked
```

It runs the pipeline on `tests/data` with `--aligner bwa` and `--aligner minibwa`, then prints one row per sample: the
largest absolute difference in any lineage abundance between the two runs (gate; exit status 1 if above the tolerance), mapped
reads and mean depth for each aligner with the percentage change, and a status. For a real comparison point `--input` at
real data instead (copy the script's two `nextflow run` calls), and also compare `summary/consensus_genomes.fasta` and
`summary/qc_summary.tsv`: a different aligner can move indel placement and soft clipping, which changes `ivar`
variant calls and therefore consensus genomes even when the Freyja abundances agree.

## What would make minibwa the default

1. A 1.0 release or a peer-reviewed publication, and an index format that no longer changes between releases.
2. A clean comparison on real wastewater runs (abundances within the tolerance, same consensus genomes, same amplicon
   dropout calls), ideally including primer-trimmed amplicons with deletions typical of Omicron lineages.
3. A multi-arch (or arm64) container, so development on Apple Silicon does not depend on emulation.
4. A reason: with 30 kb genomes it would have to be better results, not speed.

Until then `--aligner bwa` remains the default, and a run that used minibwa records it (`tool_versions` lists `minibwa`, the
container is in `containers.MINIBWA_MEM`).
