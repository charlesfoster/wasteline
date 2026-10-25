# Panel typing: file formats and channel contract (U1 -> U2)

Standalone building blocks, not wired into any workflow yet. Script: `bin/panel_select.py` (standard library only, Python 3.9
compatible; it runs in the `python:3.12.12` container, so it is NOT in `HOST_SCRIPTS` of `test_host_python_compat.py`).
Modules: `modules/local/panel_{prepare,assign,call,split}`.

## Files

**panel_resolved.csv** (input of PANEL_PREPARE; U0/U2 build it from the registry). Required columns `target`, `fasta`; used when present:
`marker_contig` (original contig name; several separated by `;`; empty = every contig is a marker), `family` (empty = the target is its own
family). Other columns (gff, freyja_*, ...) are ignored. `fasta` may be an absolute host path; inside the task a path that does not exist is
looked up by file name among the staged FASTAs, so stage them (below). Target names: no whitespace, no `|`, unique.

**panel.fasta**: every target's contigs, the FASTA header reduced to its first word. With more than one target each contig is named
`<target>|<contig>`; with exactly one target the original names are kept (legacy naming). One target = one FASTA (a FASTA shared by two
targets is read once per target).

**contig_map.tsv** (header + one row per panel contig, tab separated): `panel_contig, target, contig, is_marker (true|false), family`.

**<id>.panel_stats.tsv** (PANEL_ASSIGN): `sample, target, fragments, ambiguous_fragments, aligned_bases, mean_mapq, breadth_pct, mean_depth,
panel_fragments`; one row per target of the panel (zeros when nothing mapped), numbers on the target's MARKER contigs only.
- A *fragment* = all mapped primary alignments of one read name (secondary/supplementary dropped with `-F 0x900`; the BAM is collated here).
- Assigned to target T only if every mapped primary alignment of it is on contigs of T with MAPQ >= `min_mapq`; otherwise *ambiguous*. Both
  counts are per target and only for fragments with at least one primary alignment on that target's marker contig; `aligned_bases` (M/=/X) and
  `mean_mapq` come from the marker-contig alignments of assigned fragments.
- `breadth_pct` / `mean_depth`: `samtools coverage` over the assigned marker-contig records only (so ambiguous fragments add no depth),
  length-weighted over several marker contigs.
- `panel_fragments` (same in every row): fragments with any mapped primary alignment anywhere in the panel (0 = no panel reads).

**<id>.targets.tsv** (PANEL_CALL), one row per panel target in panel order, never empty, columns exactly:
`sample, sample_type, target, family, status, reason, sample_status, mixed, secondary_targets, fragments, ambiguous_fragments,
fraction_of_family, breadth_pct, mean_depth, mean_mapq, second_target, margin, rescue, analysed`.
- Within a family targets are ranked by fragments, aligned bases, breadth, then name. `fraction_of_family` = fragments / family fragments (0-1).
- `status`: `detected` = fragments >= min_fragments AND breadth_pct >= min_breadth_pct AND (rank 1 in its family OR share >= family_share_pct);
  `trace` = >= 1 fragment but not detected, `reason` lists the failed floors joined by `;`: `low_fragments(37<100)`, `low_breadth(2.10<5)`,
  `low_family_share(2.34<5)`; `not_detected` = 0 fragments, reason `no_fragments`, or `no_panel_reads` when the whole sample has none.
  `reason` is empty for `detected`.
- `sample_status` (same in every row): `ok` (>= 1 detected), `no_target_detected` (panel reads, nothing detected), `qc_fail` (no panel reads).
- `mixed` = `true` on the detected rows of a family with >= 2 detected targets, with `secondary_targets` = the other detected targets of that
  family joined by `;`; `false` and empty elsewhere.
- `second_target` = best other target of the family with fragments (only when this target also has some), `margin` = own minus second
  `fraction_of_family`. `mean_mapq` is `NA` with no fragments. `rescue` is always empty (placeholder). `analysed` = `true` for detected only.
- Defaults (uncalibrated): min_fragments 100, min_breadth_pct 5, family_share_pct 5 (a percent), MAPQ 10.

## Processes (channels)

| process | input | output |
| --- | --- | --- |
| PANEL_PREPARE | `path panel_resolved`, `path fastas` (all FASTAs, flat, unique file names) | `fasta` (panel.fasta), `contig_map` (contig_map.tsv). Run-level, no meta |
| PANEL_ASSIGN | `tuple val(meta), path(bam)` (panel-mapped, any order), `path contig_map`, `val min_mapq` | `stats` = `[meta, <id>.panel_stats.tsv]` |
| PANEL_CALL | `tuple val(meta), path(stats)`, `path contig_map`, `val min_fragments`, `val min_breadth_pct`, `val family_share_pct` | `targets` = `[meta, <id>.targets.tsv]`; `meta.sample_type` (default `sample`) goes into the table |
| PANEL_SPLIT | `tuple val(meta), path(bam)` (any order), `path contig_map`, `val target` | `counted` = `[meta, <id>.bam, <id>.bam.bai, <id>.mapped_reads.txt]`, same shape as `BAM_MAPPED_READS.out.counted` |

All have `topic: versions` tuples and a `stub:` block; the file prefix is `task.ext.prefix ?: meta.id`. The count file is always
`${meta.id}.mapped_reads.txt`. Containers: PREPARE/CALL reuse `python:3.12.12`, ASSIGN/SPLIT the `htslib_samtools:1.24` image (pins and
digests are in `conf/containers.lock.config`; no new image).

PANEL_SPLIT output: only the target's contigs with their ORIGINAL names (`<target>|` removed through contig_map.tsv, so a one-target run and a
split give the same reference naming), header `@SQ` trimmed to them, coordinate-sorted and indexed. A record whose mate lies on another target
(a discordant fragment) is rewritten as a read with an unmapped mate (flags 0x8 set, 0x2 and 0x20 cleared, mate = itself), so no record points
outside the trimmed header. Secondary and supplementary records are dropped. Unmapped mates placed at a mapped mate stay.

`bin/` scripts are found by Nextflow only through the pipeline's own `bin/` directory (a symlinked `bin` and a wrapper outside the repo do not work).

## Errors

`panel_select.py` exits 2 with `ERROR: ...` on: a missing FASTA, a duplicated contig name inside one FASTA, a `marker_contig` that is not a
contig, a duplicated or invalid target name, an empty table, a stats row for an unknown target. An empty/missing stats file is not an error: it
gives the all-zero rows (`qc_fail`).

## Manual check (2026-10-26)

A 2-target mini panel (H1 and H3, HA marker + NA non-marker, NA 3 % diverged), 100 + 30 + 40 + 1 simulated pairs mapped with `bwa mem`,
real `samtools 1.24` container: stats H1 100 fragments / 1 ambiguous / 20000 bases / breadth 99.00 / depth 13.33, H3 30 / 1 / 6000 / 98.60 / 4.00,
`panel_fragments` 171 (the 40 NA pairs count there only); calls H1 `detected`, H3 `trace` `low_fragments(30<100)`; split of H1: 281 mapped reads,
`@SQ SN:HA`/`SN:NA` (original names), `samtools quickcheck` ok, the discordant pair flag 73.

## U2: panel_front (subworkflow, not yet called)

`subworkflows/local/panel_front/main.nf`, `workflow PANEL_FRONT`. Maps every sample to the panel once, types it, applies the negative-control
barrier and cuts one BAM per analysed (sample, target) pair. U3 calls it from `workflows/wasteline.nf`.

**take**
- `ch_reads`: `[ meta, reads ]` after FASTP. `meta.id` and `meta.sample_type` are used; `collection_date`, `site`, `viral_load` are copied into the per-target
  meta. The sample's library cells come from `meta.protocol_cell` (or, failing that, `meta.protocol`) and `meta.primer_scheme_cell`: the RAW samplesheet
  cells (null = the run defaults). The resolved `meta.primer_scheme` / `primer_bed` of a pre-split meta are NOT reused, because they depend on the target;
  the per-target meta is rebuilt with `sampleMeta(sample, mode, target)` (which calls `libraryMeta` with that target's schemes).
- `panel_targets`: the list of target maps (`targetsFromParams()`), in panel order. (Not called `targets`: that name is also an emit and shadows the
  take inside closures.)
- `mode`: `'full'`.

**What it does.** `panel_resolved.csv` (absolute FASTA paths; written under `workflow.workDir/wasteline_panel/<digest>/`, not rewritten when it is
there, so `-resume` holds) + the distinct FASTAs (a shared FASTA once; two different files with one name are refused) -> PANEL_PREPARE -> BWA_INDEX +
BWA_MEM (`params.aligner == 'minibwa'` is refused) -> PANEL_ASSIGN (`params.panel_mapq`) -> PANEL_CALL (`panel_min_fragments`, `panel_min_breadth_pct`,
`panel_family_share_pct`) -> barrier (`toList` of all samples' targets files) -> PANEL_SPLIT per analysed pair.

**Analysed pairs.** A non-negative sample: its `detected` targets. A negative control (`meta.sample_type == 'negative'`): every target detected in at
least one NON-negative sample, whatever its own call (so a control with no panel reads still gets an empty BAM, count 0, per such target; it flows on as
the single-target path's empty-BAM branch). `trace` / `not_detected` rows of non-negative samples are not split. The `analysed` column of negative
controls is rewritten accordingly (a negative control that is detected for a target no other sample carries gets `false`). One log line:
`panel: N samples x M targets: X analysed, Y trace, Z not detected` (Y and Z count the rows that are not analysed; a trailing clause names
negative-control-only detections).

**emit**
- `mapped`: `[ meta_t, bam, bai, mapped_reads.txt ]` (PANEL_SPLIT.out.counted), one per analysed pair. `meta_t` = `sampleMeta(...)`, so `meta.target`,
  `meta.key` (`<id>__<target>`), the library settings (protocol, primer_trimming, duplicates, primer_scheme, primer_bed) are set. BAM contigs carry the target's
  original names.
- `targets`: the run's `summary/targets.tsv` (one item; also stored in `${params.outdir}/summary` through `collectFile storeDir`): header once, samples
  sorted by id, panel order within a sample, columns exactly as above.
- `typing`: `[ meta, <id>.targets.tsv ]` per sample (the pre-split meta). `contig_map`, `panel_fasta`: PANEL_PREPARE's outputs (value-like, one item).
- `analysed_pairs`: one item, a list of `[ sample id, target name ]`.

**PANEL_SPLIT changed**: input is now `tuple val(meta), path(bam)` + `path contig_map`; the target is `meta.target.name` (the separate `val target` is
gone). The BAM is staged under `panel_mapped/` because it is named like the output and `samtools sort -o` writes through a symlink (found in the first
real run: it overwrote the upstream BWA_MEM BAM and the two splits of one sample raced).

Checked (2026-10-26, scratch wrapper outside the repo, real docker run, RSV panel, `panel_min_fragments=50`, ~300 pairs/sample): pure RSVa, pure
RSVb, 70/30 mix (`mixed=true`, 213 / 91 fragments), RSVa + 41 RSVb fragments (`trace`, not split), NTC with 5 RSVa pairs (analysed for RSVa and RSVb).
