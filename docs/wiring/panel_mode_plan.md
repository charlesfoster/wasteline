# Panel mode (PLAN D25 phase 2/3): plan, 2026-10-22 (nothing built)

Several pathogens, subtypes or genotypes per run and per sample. Produced by a read-only planning pass that read the code and the
user's own pipelines (h2seq, quasispecies_pipeline, clustering_workflow). Phase 1 (multi-contig targets) is done.

## Recommendation: architecture B (typing + partition front end, then the existing pipeline once per target)

A new `--mode panel` per sample: FASTP -> competitive `bwa mem` against the concatenated panel -> per-fragment target assignment and
per-target statistics -> detected / mixed / trace call -> partition of reads per target (each fragment once, plus an ambiguous pool) ->
`panel/targets/<T>/samplesheet.csv` + FASTQs. A thin host driver `bin/run_panel.py` (dry run by default, like `reanalyse_runs.py`) then
runs the UNCHANGED pipeline once per target (`--pathogen <T> --input ... --outdir out/targets/<T>`) and builds a merged `index.html`
(matrix sample x target with status, fragments, breadth, top lineage, links to each target's report).

Why not A (in-pipeline `meta.target`): about 150-200 edit sites (26 `params.pathogen`, 21 `params.fasta`, 14 `params.freyja_contig`, 76
publishDir sites, 113 `meta.id` uses, six config-time hidden facts that cannot be per target), mostly in the files that produce SARS
output, and it breaks the rule that code never names a pathogen. B changes no existing workflow/module/bin script (SARS risk about
zero), is about 1300 lines of new files, and gets demix, aggregate, report, time series, site aggregation, lineage summary, the runs store
and re-analysis per target for free. Costs: reads mapped and fastp-trimmed twice, extra disk for partitioned FASTQs, no single combined
report, per-target read totals exclude other targets' reads.

## Panel definition (CSV)

Columns `target,reference,fasta,gff,freyja_pathogen,freyja_contig,marker_contig,family,primer_scheme,barcodes`. `target` = the unit that is
called and the `--pathogen` name (RSVa, H3N2); `reference` = optional candidate id (phase 3: several rows per target = alternative references,
the winner by aligned bases becomes that target's reference.fasta); `fasta` may have several records (segments); `marker_contig` = contigs
used for typing (default `freyja_contig`, else all; flu internal segments cross-map between subtypes at 80-95 %, so type on HA);
`family` = explicit grouping for mixed calls (RSV, fluA, fluB, HCV-1), a sample is `mixed` only when two or more targets of one family
are detected; `barcodes` = optional custom CSV. Built-in `assets/panels/rsv.csv` (RSVa PP109421.1 + RSVb OP975389.1) and `flu.csv`
(H1N1 CY121680.1_wg, H3N2 CY163680.1_wg, H5Nx AF144305.1_wg, B/Vic KX058884.1 HA only); custom `--panel path.csv`. Contig names must be unique
across the panel (`prepare` errors on a clash and writes the contig-to-target map). A pytest consistency check ties built-in rows to the
profile configs.

## Mapper, statistics, decision (adapted from h2seq select_reference_from_bam.py and filter_mixed_reference_alignments.awk)

- `bwa mem` (reuse BWA_INDEX/BWA_MEM/FASTP, same aligner and MAPQ semantics as the per-target run), `sort_bam=false` (query-name grouped);
  `PANEL_`-aliased includes with their own 256 MB BWA_INDEX memory floor. minimap2 only if a large candidate panel needs it (phase 3).
- Fragment assignment (awk in a samtools container): both mates' primary alignments on the same target and every mapped primary MAPQ >=
  `panel_min_mapq` (10), else `ambiguous`. The quasispecies AS-XS delta is NOT used (bwa MAPQ already derives from it; a raw XS cannot tell
  a same-target alternative from a cross-target one); phase 3 may parse the XA tag for between-target ambiguity (untested).
- Per (sample, target): assigned fragments, ambiguous fragments shared with a family member, aligned bases and identity, and from
  `samtools coverage -q <mapq>` on the marker contigs: breadth >= 1x, mean depth, mean MAPQ. Rank within a family by aligned_bases, breadth,
  identity, MAPQ, name; report the second target, the margin and `close_target` (clustering_workflow close-hit idea).
- Status per row (always present, never an empty file): `detected` = marker fragments >= `panel_min_fragments` (100) AND breadth >=
  `panel_min_breadth_pct` (5) AND the target's share of its family >= `panel_min_fraction` (0.05; the top target is exempt); `trace` = at
  least one fragment but below a floor (reason names the floor); `not_detected`; `sample_status = qc_fail` (`no_reads`/`no_panel_reads`)
  when nothing maps; `mixed=true` + `secondary_targets` for two or more detected targets of one family. All thresholds are UNCALIBRATED.
- Partition: a detected target gets its assigned fragments; the only detected target of a family also gets every fragment whose primary
  alignment is on it (recovers internal flu segments); in a mixed family unresolved internal-segment fragments go to the ambiguous pool with
  a WARN flag `mixed_family_internal_segments_unresolved`.

## Change list

New: `bin/panel_select.py` (stdlib; `prepare`, `call`), `bin/run_panel.py` (stdlib, py>=3.8; `--execute`, `--keep-going`, `--merge`,
`--extra-args`, `--target-args`), modules `panel_prepare`, `panel_assign`, `panel_call`, `panel_partition` (each with a stub nf-test),
`workflows/panel.nf` (~150 lines), `assets/panels/{rsv,flu}.csv`. Edited: `main.nf` (dispatch `panel`, typed `panel_*` params),
`nextflow.config`, `nextflow_schema.json` (`mode` enum + `panel*` params), the utils subworkflow (modes list; `sampleMeta` for panel keeps raw
cells; `--panel` required, `--pathogen`/`--fasta` not; `panel*` params filtered out of run_info unless mode == panel), appended `PANEL_*`
selectors in `conf/modules.config`, one line in `bin/reanalyse_runs.py`. LITERALLY UNCHANGED: `workflows/wasteline.nf`, `conf/pathogens/*`,
results_input, all existing modules, report code, `make_samplesheet.py`, `check_samplesheet.py`.

Layout: `out/panel/targets.tsv` (sample,target,family,status,reason,sample_status,mixed,secondary_targets,fragments,ambiguous,
fraction_of_family,breadth_pct,mean_depth,mean_mapq,identity,second_target,margin,close_target), `out/panel/reads/<T>/...`,
`out/panel/targets/<T>/{samplesheet.csv,reference.fasta}`, `out/targets/<T>/` (a normal results directory), `out/index.html`. Negative
controls go into every target samplesheet. Unknown: whether fastp/the pipeline accept a 0-read FASTQ (else a placeholder; the e2e shows it).

## Phase 3 (targets without barcodes: enterovirus, HCV)

`--mode panel` already gives the best-reference call, coverage and the winner reference per target without Freyja. Consensus/QC for such a
target needs one new single-target switch `--skip_demix` (about 6-8 guarded sites in `wasteline.nf`, a `not_applicable` demix row), a separate
later task gated on SARS byte-identity again. Unknown: whether `freyja variants`/the generic profile accept an arbitrary `--pathogen`.

## Verification (minimal)

ONE new e2e step (10), SARS-free and offline: simulate with `tests/simulate_pathogen_mixture.py` (one-row barcodes CSV made inline): S1 pure
RSVa, S2 pure RSVb, S3 70/30 RSVa/RSVb, S4 RSVa + ~40 RSVb pairs (trace), NTC near-empty; `--mode panel --panel rsv --protocol capture`;
assert in a `--panel` section of `check_outputs.py`: S1 detected RSVa and not detected/trace RSVb, S2 the reverse, S3 `mixed=true` with
fractions within 0.05 of truth, S4 RSVb `trace` with no reads in its samplesheet, NTC row present, each fragment in exactly one FASTQ,
every generated samplesheet passes `check_samplesheet.py`, the `run_panel.py` dry run prints one command per target. The per-target runs are
the unchanged pipeline (RSV validated manually; the barcodes need the network). SARS byte-identity: empty `git diff` of the unchanged files,
e2e steps 1-9 `diff -r` against the baseline, zero nf-test snapshot diffs, identical `run_info.json` params. New tests: one
`tests/python/test_panel.py` (one test per script), four stub nf-tests.

## Coder tasks (ordered, strict file ownership)

T1 scripts + panels (`bin/panel_select.py`, `assets/panels/*.csv`, pytest, the profile-consistency check); T2 modules (`panel_*`, needs T1's CLI);
T3 wiring (`workflows/panel.nf`, `main.nf`, config, schema, utils, modules.config selectors); T4 driver, e2e step 10, docs (`bin/run_panel.py`,
reanalyse_runs one-liner, e2e, check_outputs, docs, PLAN D28); T5 later: `--skip_demix` (owns `workflows/wasteline.nf`).

## Top risks

RSV-A reads cross-mapping onto RSV-B inflating a trace into `detected` (floors; measured with S4); MAPQ ties among same-target references in
candidate-rich panels (phase 3, XA-based ambiguity); thinly covered flu internal segments in mixed samples (WARN flag); per-target read
totals exclude other targets' reads; primer scheme names are pathogen-specific (BED path works for every target); uncalibrated thresholds; a
0-read FASTQ may break fastp; double fastp changes read counts slightly.

## Decisions for the user

1. Go with B (no single combined report, a second mapping)? 2. Keep the starting thresholds (100 fragments, 5 % breadth, 5 % fraction) and
calibrate on real samples? 3. Should `trace` samples stay out of the per-target runs (default) or be included? 4. Flu panel = H1N1 + H3N2 +
H5Nx + B/Vic (HA only)? Whole-genome capture or amplicon? 5. Phase 3 (`--skip_demix`) wanted later?
