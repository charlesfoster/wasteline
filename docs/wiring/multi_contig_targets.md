# Multi-contig targets (PLAN D25 phase 1): the build contract

Whole-genome influenza first: a single-target run whose reference has several contigs (eight segments). Freyja still demixes ONE contig
(the HA segment, profile fact `freyja_contig`); every contig gets depth/coverage QC, variant calling, a consensus genome and a place in
the report. SARS-CoV-2 and every existing single-contig run stay BYTE-IDENTICAL (outputs and existing table structure); new columns and
files appear only for multi-contig references. Written 2026-10-22 from a read-only planning pass that ran `freyja variants` and
`bcftools consensus` in the pinned containers.

## Decisions (planner)

- Data model **(b)**: one pass with all-contig files plus a freyja-contig view (`FREYJA_SELECT_CONTIG`), not a per-contig fan-out
  (about four times the tasks, high SARS byte-identity risk). Gate: `multiContig() = fastaContigs(params.fasta).size() > 1`.
- `--whole_genome` (visible boolean, default false) selects the bundled whole-genome influenza reference; HA-only stays the default (D24).
  A custom multi-record reference also works: `--fasta multi.fa --freyja_contig <name>`.
- Whole-genome references to bundle: H1N1 (A/California/07/2009: CY121680.1 HA, CY121681.1 M, CY121682.1 NA, CY121683.1 NP, CY121684.1 NS,
  CY121685.1 PA, CY121686.1 PB1, CY121687.1 PB2), H3N2 (A/Wisconsin/67/2005: CY163680.1 HA ... CY163687.1, same order) and H5N1
  A/Goose/Guangdong/1/96 (AF144300.1 PB2, AF144301.1 PB1, AF144302.1 PA, AF144303.1 NP, AF144304.1 NA, AF144305.1 HA, AF144306.1 M,
  AF144307.1 NS; shared by H5Nx and H5Nx-cattle). **B/Victoria whole genome is NOT built now**: there is no consecutive set for the
  strain (a composite of three submissions); B/Vic stays HA-only until a proper set is chosen. H5Nx-cattle: bundle the Goose/Guangdong
  set with a documented limit (internal genes 10-15 % divergent from 2.3.4.4b cattle viruses: reads map poorly, the consensus is
  reference-biased); a closer composite can replace it later.
- QC rule for multi-contig: `qc_status` stays on the WHOLE genome with the existing thresholds; add a WARN-only reason
  `contig_low_coverage` when any contig has `pct_ge_min_depth` below the PASS threshold (never FAILs by itself, like
  `low_primer_trimming`). The negative-control rule is unchanged. The demix gate uses the freyja contig.

## Verified facts (container tests)

- `freyja variants BAM --ref multi.fa` without `--refname` runs `samtools mpileup -aa ... | tee >(cut -f1-4 > depths) | ivar variants`:
  depths (col 1 = contig, every position of every contig) and the variants TSV (REGION column) cover ALL contigs. `awk '$1==X'` of the
  depths and `REGION==X` of the variants are byte-identical to `--refname X`.
- `freyja demix` cannot take multi-contig depths (positions collide, coverage is diluted): demix, boot and the gate must get the
  freyja-contig-only pair.
- `bcftools consensus -f multi.fa -H I -m mask.bed` emits all records; masked contigs come out all-N; headers keep the full FASTA header.
- Already multi-contig-safe: `ivar_variants_to_vcf.py`, `low_depth_mask.py` (except `--length`), `depth_stats.py`, `qc_summary.py`
  (position-level, so length-weighted whole-genome numbers for free), `annotate_variants`, `IVAR_TRIM`, `SAMTOOLS_*`.
- Still single-contig and SILENTLY wrong (not failing): `mutations_of_interest.read_fasta_sequence` (first record), `read_cds_models`
  (ignores GFF seqid), `read_depths`; `build_report_data.read_depths` keeps the first contig; `consensus_stats.read_single_fasta` raises on
  more than one record; `amplicon_coverage.py` has no contig column; `render_report.amplicon_view` keys amplicons by name.

## Change list

Scripts (`bin/`, task T1):
- `qc_summary.py`: `--contig-table FILE` (columns sample, contig, length, mean_depth, median_depth, pct_ge_1x, pct_ge_10x,
  pct_ge_min_depth, amplicons_total, amplicons_dropped; amplicon counts from the `contig` column of `--amplicons`) and
  `--warn-contig-below-pct X` (default 0 = off). `COLUMNS` unchanged.
- `consensus_stats.py`: with `--contig-out FILE` accept several records. The sample row is aggregated (`length` = total, `n_count` = sum,
  `consensus_coverage_pct` = (total - N)/total, `ambiguous_count` = sum, status on the whole, same formulas); the per-contig file has
  columns sample, contig, length, n_count, consensus_coverage_pct, ambiguous_count, status (contig = part after the last `|`). Without
  the flag, more than one record stays an error.
- `amplicon_coverage.py`: when the depths have more than one contig, append a trailing `contig` column; single contig byte-identical.
- `dropout_impact.py`: `--contig C` keeps that contig's rows of `amplicon_coverage.tsv` (required when the table has a `contig` column).
- `unexplained_mutations.py`, `mutations_of_interest.py`: `--contig C`; `read_fasta_sequence(path, contig=None)` and `read_cds_models`
  (by GFF seqid) become contig-aware; without `--contig` today's behaviour; multi-record input WITHOUT `--contig` must be an error.
Report (task T2): `build_report_data.py`, `render_report.py`, `assets/templates/report.html.j2`:
- `build_coverage_tracks` concatenates contigs into one axis with `segments: [{contig,start_bin,n_bins}]` (global bin_size =
  ceil(total/points)); one contig keeps the old code path. Template: segment separators/labels on the coverage cards, the note "across N bp
  in K segments", a "Segments" table (per-contig mean depth, % >= 10x, consensus coverage %, amplicons dropped), header line "Freyja demixes
  <contig> only; other segments are QC/consensus only".
- `load_amplicons`, multi runs only: `amplicon = "<contig>/<name>"` and contig offsets added to start/end (offsets from the depths /
  `run_info.contigs`), so `render_report.amplicon_view` needs no change; keep the `contig` field.
- Run flag `segments_not_covered` (warn): samples with a contig whose `pct_ge_10x` is below the PASS threshold or whose consensus is
  `empty`, naming the segments. `dropout_impact` shows the freyja contig only (a note says so). The report must not call whole-genome
  `mapped_pct` "HA".
- Inputs: `contig_coverage.tsv`, `consensus_stats_contigs.tsv`, the `contig` column of the amplicon table, `run_info.contigs`.
Pipeline (task T3, later): `workflows/wasteline.nf`, `conf/modules.config`, profiles, schema, utils subworkflow, modules:
`FREYJA_SELECT_CONTIG` (new, stub nf-test; filters variants+depths to the contig; output basenames stay `<id>.variants.tsv` /
`<id>.depths`; fails loudly if empty); demix/boot/MOI/unexplained get the freyja view, QC/amplicon/consensus/annotate/report get the
all-contig view; `low_depth_mask` omits `--length` for several records; `consensus_bcftools` writes `>id|<first word of header>` per contig
(one record keeps `>id`), `consensus_stats` sums lengths and passes `--contig-out`, `consensus_concat` accepts N records per file;
`annotate_variants` keys on contig+position; `ext.contig` for unexplained/MOI/dropout; `freyja_demix` gate text names the contig;
start-up validation (more than one record allowed only if `freyja_contig` is set and present; duplicate first-word names are an error;
GFF seqids vs FASTA names); `freyja_contig` hidden profile fact (null in every existing profile, excluded from run_info like the others);
`--whole_genome`; summaries `contig_coverage.tsv` and `consensus_stats_contigs.tsv` (multi only); run_info gains `contigs` and
`freyja_contig` only when multi.
References (task T0): `assets/references/<HA accession>_wg.fasta` + `.gff3`, HA record FIRST and byte-identical to the existing HA FASTA;
per segment NCBI gene/CDS plus synthetic mRNA/exon as for HA (M2 and NEP/NS2 are two-exon CDSs; PA-X may be skipped and noted).

## Verification (minimal; the user's standing rule)

- pytest: extend EXISTING files with one happy-path assertion per new option (no new test files); one new stub nf-test for
  `freyja_select_contig`; existing snapshots must not change (any `.snap` diff is a SARS regression signal).
- ONE e2e addition (step 9): `multi.fa` = SARS NC_045512.2 + measles NC_001498.1 (a real second contig, offline); extend
  `tests/simulate_pathogen_mixture.py` to multi-record references (`--freyja-contig` applies barcode mutations only there; reads drawn from
  all contigs by length); S1 simulated (B.1.1.7 + measles reads), S2 reuses a step-1 SARS pair (measles uncovered); run with
  `--pathogen SARS-CoV-2 --fasta multi.fa --freyja_contig NC_045512.2 --protocol capture`; a short `--multicontig` check in
  `tests/check_outputs.py` (contig table 2x2, consensus records `S1|NC_045512.2` and `S1|NC_001498.1`, consensus length sum, S2
  `contig_low_coverage` and `segments_not_covered`, S1 demix = truth, S2 `.demix` cmp-identical to the single-contig run's).
- SARS byte-identity: baseline e2e outputs of `dev` at 0ac7897 are in the scratchpad (e2e9); after the change `diff -r` excluding
  `pipeline_info/`, `report.html`, `report_data.json`; compare `report_data.json` with provenance/header keys removed; zero nf-test
  snapshot diffs.

## Top risks

1. Silent wrong contig in scripts that accept multi-record input: always pass `--contig`; error without it.
2. Wiring a consumer to the all-contig view by mistake (demix/boot/MOI/unexplained): add a unit assertion that the demix input is
   freyja-contig only.
3. Amplicon names colliding across segments (whole-genome BED naming conventions unknown).
4. Whole-genome QC passing with a missing segment (hence `contig_low_coverage`).
5. NCBI complete-cds references lack UTRs; H5 Goose/Guangdong internal genes fit 2.3.4.4b cattle viruses poorly.
6. Demix mode on multi-contig results only works through `FREYJA_SELECT_CONTIG`; the e2e runs full mode only.
