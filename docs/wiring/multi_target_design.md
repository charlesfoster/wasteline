# Multi-target design (PLAN D25, proposal): several pathogens, genotypes or segments in one run

Status: **proposal, nothing built**. Written 2026-10-22 after the user asked for (a) whole-genome influenza without throwing
information away when Freyja only demixes HA, and (b) runs whose samples contain several pathogens, subtypes or genotypes at once
(RSV-A and RSV-B together, several influenza subtypes, enterovirus or HCV genotypes; Freyja has no barcodes for the last two yet).

## What the pipeline assumes today

One run = one pathogen = one single-sequence reference (`params.pathogen`, `params.fasta`). Reads are mapped once to that reference,
and everything downstream (primer trimming, Freyja variants/demix, depth QC, amplicon coverage, consensus, aggregates, plots, summaries,
report) assumes one contig and one Freyja dataset. D20 already moved the *library* settings (protocol, primer scheme) into `meta` per
sample; D22-D24 made the reference/pathogen a per-run profile.

## Concepts

- **Target**: one thing to measure in a sample. A target has a name, one or more contigs (a segmented genome has several), an optional
  Freyja analysis (`freyja_pathogen` plus the single contig that carries the barcode positions, e.g. HA for influenza), an optional GFF,
  an optional default primer scheme. Examples: `SARS-CoV-2` (1 contig, Freyja SARS-CoV-2), `RSVa` (1 contig, Freyja RSVa), `H3N2` (8 contigs,
  Freyja H3N2 on the HA contig), `HCV-1a` (1 contig, no Freyja analysis: QC and consensus only).
- **Panel**: the list of targets a run looks for. Today's runs are panels of one target (`--pathogen X`, unchanged). New:
  `--panel rsv` (built-in: RSVa + RSVb) or `--panel panel.csv` (user file: `target,fasta,gff,freyja_pathogen,freyja_contig,primer_scheme`)
  or `--pathogen RSVa,RSVb` as shorthand for a panel of the bundled targets.

## Data flow

1. **One mapping per sample against the whole panel**: the targets' contigs are concatenated into one FASTA (contigs named
   `<target>|<contig>` when more than one target, so names never clash; a single-target panel keeps the original contig names so
   SARS-CoV-2 output is byte-identical). Reads of similar targets (RSV-A vs RSV-B, HCV genotypes) compete for their best reference, which
   is the standard way to avoid cross-mapping; mapping each target separately would let RSV-B reads pile up on the RSV-A reference.
2. **Split by target after mapping**: `samtools view -b` of the target's contigs gives one BAM per (sample, target). From here the existing
   chain runs per (sample, target): primer trimming (the target's scheme), depth/coverage QC, amplicon coverage and dropout, consensus
   (per contig), Freyja variants and demix on the Freyja contig only (barcodes of that pathogen), boot, mutations of interest.
3. **Detection**: a (sample, target) pair whose genome breadth is below the existing floor is `not_detected` (the existing
   `skipped_low_coverage` status with a reason, so a failed or absent target is a row, never an empty file). The coverage columns say
   how much of the target was seen, so a sample with a trace of RSV-B shows exactly that.
4. **Aggregation per target**: aggregates, plots, time series, lineage groups, site aggregation, unexplained mutations and the cross-run
   store run once per target over the samples where that target was detected (or analysed).
5. **Outputs**: with one target, the layout is unchanged. With several: `samples/<id>/<target>/...` and `summary/<target>/...`, plus a
   run-level `summary/targets.tsv` (sample x target: detected, genome breadth, mean depth, top lineage, status) and one report: a header
   block with the sample x target matrix and one collapsible section per target (the existing report content, per target).

## Prerequisite: multi-contig support inside one target (whole-genome influenza)

A target with several contigs needs per-contig handling that the code does not have. The table in `docs/other_pathogens.md` ("Route B") lists
the places (depth QC per contig and genome-wide, amplicon coverage, low-depth mask and consensus per contig, report depth tracks per contig,
variants/mutation tables with a contig column, unexplained mutations restricted to the Freyja contig). Influenza then runs whole-genome (all
eight segments mapped, QC'd, consensus per segment) while Freyja demixes HA only, so no information is thrown away.

## Reference selection and what is borrowed (2026-10-22, from the user's own pipelines; non-SARS-CoV-2 only)

SARS-CoV-2 runs never select a reference: a single-target panel is mapped as today. Selection only exists for panels with several
candidate references per target group (HCV genotypes, enterovirus types, RSV-A/B, influenza subtypes). Borrowed, in order of value:

1. `~/Programs/h2seq` (MIT, the user's): `bin/select_reference_from_bam.py` (per-reference stats from a competitive-mapping BAM: mapped
   reads, aligned bases, breadth, mean MAPQ, identity; ranking; dominant group by summed reads; any other group with >= 5 % of the reads
   is a secondary, i.e. a mixed sample; an empty result is a `qc_fail` row with a reason) and `bin/filter_mixed_reference_alignments.awk`
   (a read pair counts for one reference only if its primary alignments agree and MAPQ >= 10; otherwise it is ambiguous and excluded).
   Adapt: group references by an explicit panel table (`target`, `contig`), not by a name prefix, so segmented genomes work.
2. `~/Programs/quasispecies_pipeline` `bin/classify_genotype.py`: a read is ambiguous when the best and second-best alignment scores
   differ by less than a delta (default 20, from the AS and XS tags, which bwa also writes); counts exclude ambiguous reads; a secondary
   subtype needs >= 5 % of the usable reads. Same competitive-panel idea with minimap2 for long reads: here bwa mem is kept (short reads).
3. `~/Programs/clustering_workflow`: only the reporting idea (score margin, list of close competing hits, pass/fail with reason) and the
   238-record HCV panel as a starting point for an HCV panel; its consensus-level clustering is out of scope.
4. `charlesfoster/rsv_pipeline`: nothing copied (no licence file; it takes the top sylph abundance only, no breadth or depth). Sylph stays a
   possible fast screening step for very large panels.

Not borrowed: any HCV-only logic, GUI/plots, clustering, per-genotype haplotype callers.

## Phases (each ends with an e2e pass and a SARS-CoV-2 byte-identity check)

1. **Multi-contig targets** (Route B): per-contig depth/coverage/consensus/report, `freyja_contig` setting, a bundled whole-genome influenza
   reference per subtype (H1N1pdm09, H3N2, B/Victoria, H5Nx: eight segments each, HA = the Freyja reference; GenBank accessions to be
   chosen and verified like the HA references). One subtype per run.
2. **Panels**: `meta.target`, one mapping to the concatenated panel, split by target, per-target fan-out and aggregation, per-target
   output layout, `targets.tsv`, report sections. Built-in panels `rsv` (RSVa+RSVb) and `flu` (all influenza subtypes), user panels by CSV.
3. **Targets without Freyja barcodes** (enterovirus types, HCV genotypes): QC, coverage, consensus and a best-reference call per sample
   (which genotype the reads support, by breadth and depth), no demixing until barcodes exist; a user can supply a barcodes CSV
   (`--barcodes`) per target as Freyja already allows.

## Risks and how they are contained

- **SARS-CoV-2 regression**: single-target panels must produce identical files. The e2e already checks demix/aggregate byte-identity;
  add nothing else: the existing run is the guard.
- **Size**: phase 2 touches most of `workflows/wasteline.nf` and `conf/modules.config`. It is mechanical (closures read `meta.target`;
  run-level processes are called on `groupTuple(by: target)` channels) but large. Do phase 1 first and keep phase 2 behind the same
  `--pathogen X` interface.
- **Cross-mapping between close targets**: competitive mapping to one panel reference avoids most of it; measure with simulated RSV-A
  reads against RSV-B and the reverse, report the residual (reads on the wrong target) and set a MAPQ rule if needed.
- **Tests**: no new test files. One extra e2e block per phase at most; unit tests only for new scripts (target-table builder, BAM split
  helper) with one happy path each.

## Decisions needed from the user

1. Is the three-phase order right (whole-genome influenza first, then panels, then barcode-less targets)?
2. For influenza panels: which subtypes besides H1N1, H3N2, B/Victoria and H5Nx matter (B/Yamagata is extinct; H7, H9?), and are samples
   whole-genome capture or whole-genome amplicons (ARTIC-style flu schemes)? That decides the primer schemes to bundle.
3. For enterovirus/HCV: which genotypes/types and which references (a type-strain set per genus?), and is a best-reference call enough
   at first?
