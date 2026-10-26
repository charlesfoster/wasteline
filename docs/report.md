# The HTML report

`summary/report.html` is one self-contained file: inline CSS, SVG and a little vanilla JavaScript, no external
fonts, scripts or images (the shipped report contains no URL), so it can be e-mailed, archived with the results,
opened offline and printed (print styles are included; it also follows the browser's light/dark setting). It is
rendered from `summary/report_data.json` by `bin/render_report.py` using the Jinja2 template
`assets/templates/report.html.j2`. Every piece of data in the report comes from the other output files described in
[output.md](output.md); the report adds no results of its own, except the merging of QC, demix and consensus status per
sample.

The screenshots below are of a real `full` run on the synthetic test dataset (7 samples: three mixtures, one pure
sample, one low-coverage sample, one clean and one contaminated negative control) with `--barcodes container`,
`--compare_barcodes`, `--annotate_variants` and `--site_aggregation` (bootstrap intervals, `--freyja_boot`, are not shown in
the screenshots; see the tooltip notes below).

![Report header, run-level flags and QC table](images/report_overview.png)

![Lineage abundance and coverage tracks](images/report_lineage_coverage.png)

The contaminated control `NEG_CONTAM` is a QC `FAIL` (`ntc_reads`) and raises the red run-level flag, but its demix is
`ok`: in the lineage panel it carries an `ok` badge (no bar, because controls are left out of the aggregate) and its
content is in the "Negative controls" section:

![Negative controls and barcode comparison sections](images/report_controls_comparison.png)

(The "Demix status" and "Time series" sections that lie between these two sections in the report are not shown in
the third image.)

Contents: [header](#header-and-barcode-version) | [run-level flags](#run-level-flags) | [QC](#quality-control) |
[lineages](#lineage-abundance) | [coverage tracks](#coverage-tracks) | [amplicon dropout](#amplicon-dropout) |
[dropout impact](#amplicon-dropout-impact) | [negative controls](#negative-controls) |
[site aggregates](#site-aggregates) |
[demix status](#demix-status) | [mutations of interest](#mutations-of-interest) |
[unexplained mutations](#unexplained-mutations) | [time series](#time-series) |
[barcode comparison](#barcode-comparison) | [consensus](#consensus-genomes) | [provenance](#provenance) |
[modes](#what-is-available-in-each-mode) | [rebuild](#rebuilding-the-report) | [known issues](#known-issues)

## Controls used throughout

- **Filter bar** above most panels: a text filter on the sample name, a QC filter, "only WARN/FAIL", and "hide controls"
  (negative controls). Filters work per panel and are meant for runs with 100+ samples; the QC table and the demix and
  consensus tables are sortable by clicking a column header.
- **Negative controls** carry an `NTC` chip next to their name and are printed in italics.
- **Hover** any bar segment, heatmap cell or coverage track for a tooltip (lineage and abundance, amplicon and mean depth,
  allele counts, ...). With `--freyja_boot` the lineage tooltips also show the 95% bootstrap interval, for example
  `BA.2: 71.2% (95% interval 67.0-72.7%)`.
- Colours are a colour-blind-safe categorical palette (10 colours; further categories are grey "Other").

## Header and barcode version

The top of the page shows the title, the mode and pipeline version and a prominent box with the **Freyja barcode
version** the lineage calls were made with (`09_17_2026-00-46`, or `custom` for your own barcode file). Always read
lineage results together with this version: the same reads give different lineage names with different barcodes.
The box also shows the source (`latest`, `container`, `pinned`, `custom`), the ISO date, and for custom barcodes the
start of the file's sha256. Under it: number of samples (and controls), pipeline version, git commit, Freyja and
Nextflow versions, run start time, and the time the report data was generated. "Key parameters" unfolds the main
settings used. A table of contents links to the sections; sections that have no data for the run are omitted from it
(negative controls, dropout impact, mutations of interest, unexplained mutations, time series, barcode comparison) or
show a "not available" note.

In the report of one target of a `--panel` run (`summary/<target>/report.html`) the subtitle also names the target
(`... pipeline 0.1.0 · target RSVa`) and a small "&larr; overview" link above the title opens `../overview.html`, the
[panel overview](usage.md#panels-several-targets-in-one-run). Reports of single-target runs have neither.

## Run-level flags

A short list of things you should know before trusting the run, most severe first: `ERROR` (red), `WARN` (amber),
`INFO` (blue). With no flags, a green "no flags" note is shown. Each flag has a machine code (shown at the right) and
the affected samples where relevant. All codes:

| Code                           | Level | Raised when                                                                                                                                                                                                                                               |
| ------------------------------ | ----- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ntc_reads`                    | error | A negative control has mapped reads above `--ntc_max_mapped_reads` (or QC reason `ntc_reads`): possible contamination.                                                                                                                                    |
| `ntc_fail`                     | error | A negative control failed QC for another reason.                                                                                                                                                                                                          |
| `demix_failed`                 | error | `freyja demix` failed for one or more (non-control) samples.                                                                                                                                                                                              |
| `all_samples_failed`           | error | Every non-control sample failed QC or has no usable lineage result.                                                                                                                                                                                       |
| `freyja_data_missing`          | error | Lineage results exist but `run_info.json` records no `freyja_data.barcode_version`: the barcode dataset is unknown.                                                                                                                                       |
| `mixed_barcode_versions`       | error | Samples were demixed with different barcode versions.                                                                                                                                                                                                     |
| `qc_fail`                      | warn  | One or more (non-control) samples have QC `FAIL`.                                                                                                                                                                                                         |
| `demix_low_confidence`         | info  | One or more samples were demixed with less genome coverage than `--min_genome_coverage_for_demix` (status `low_confidence`): treat their estimates with care. Not counted as skipped.                                                                     |
| `demix_skipped`                | warn  | Demixing was skipped for low genome coverage (`skipped_low_coverage`) for one or more samples.                                                                                                                                                            |
| `demix_missing`                | warn  | No demix result was found for one or more samples (status `missing`).                                                                                                                                                                                     |
| `barcode_compare_failed`       | warn  | With `--compare_barcodes`: the second demix (the other barcode dataset) failed for the listed samples; they have no rows in the comparison.                                                                                                               |
| `barcode_unknown`              | warn  | The barcode version is `unknown` (legacy Snakemake results).                                                                                                                                                                                              |
| `previous_run_barcode_unknown` | warn  | A previous run the results build on recorded no barcode version (legacy results).                                                                                                                                                                         |
| `no_samples`                   | warn  | No samples at all; the report is empty.                                                                                                                                                                                                                   |
| `input_problems`               | warn  | Some input rows or values could not be parsed (with the first messages).                                                                                                                                                                                  |
| `primer_scheme_mismatch`       | warn  | The read starts of the listed samples line up much better with another primer scheme than with the declared `--primer_scheme`; the message names the best match(es) (needs the primer-scheme QC columns).                                                 |
| `dropout_affects_lineages`     | warn  | Amplicon dropout hides lineage-informative sites of lineages detected in the listed samples (see [dropout impact](#amplicon-dropout-impact)); abundance estimates may be biased.                                                                          |
| `unexplained_mutations`        | info  | At least one mutation that the demixed lineages do not explain recurs in several samples (see [unexplained mutations](#unexplained-mutations)); possible novel lineage.                                                                                   |
| `lineage_groups_warnings`      | warn  | The groups audit (`lineage_groups_resolved.tsv`) has warnings: a pattern that is not a lineage of the dataset, a group that matches nothing, or in-total groups that overlap; the first three are in the message (see [lineage groups](#lineage-groups)). |
| `lineage_groups_ungrouped`     | info  | An `ok` sample has 2% or more (`--ungrouped-flag-threshold`) of its lineages outside every in-total group; the message names the top ungrouped lineages and says to add them to your lineage groups file (see `summary/lineage_ungrouped_summary.tsv`).   |
| `site_aggregate_incomplete`    | info  | Site aggregate periods are missing member sites or used a weaker weighting than requested (see [site aggregates](#site-aggregates)).                                                                                                                      |
| `segments_not_covered`         | warn  | Multi-contig reference only: a sample has a segment with fewer than `--qc_pass_min_pct` % of its positions at 10x or an empty consensus; the message names the segments. The whole-genome QC can still pass.                                              |
| `consensus_failed`             | warn  | The consensus genome failed for the listed samples (status `failed`, all-N sequence); the consensus table shows the reason, for example `bcftools consensus failed: The fasta sequence does not match the REF allele at NC_045512.2:11287:`.              |
| `consensus_warnings`           | info  | The consensus was built (`ok`) with warnings for the listed samples, for example `skipped variant at NC_045512.2:22942 (overlaps another variant)`.                                                                                                       |
| `qc_warn`                      | info  | One or more samples have QC `WARN`.                                                                                                                                                                                                                       |
| `barcode_custom`               | info  | A custom barcode file was used; its version is not tied to a Freyja release (use the checksum).                                                                                                                                                           |
| `barcode_container`            | info  | The snapshot bundled in the Freyja container was used; it may be older than the latest Freyja data.                                                                                                                                                       |
| `previous_run_barcode_differs` | info  | A previous run used a different barcode version than this run.                                                                                                                                                                                            |
| `no_ntc`                       | info  | The run has no negative-control sample.                                                                                                                                                                                                                   |
| `qc_not_available`             | info  | No QC summary was provided (for example `aggregate` mode).                                                                                                                                                                                                |
| `input_note`                   | info  | Free-text notes from the report builder.                                                                                                                                                                                                                  |

In the example run: `ntc_reads` (the contaminated control `NEG_CONTAM`, 5998 mapped reads), `qc_fail` and `demix_skipped`
(`LOWCOV_D`), `barcode_custom`, and `qc_warn` (`Site_1_week_2`). A negative control that demixed `ok` raises no demix flag
and is never counted as a failed sample.

## Quality control

Counters first: how many samples PASS, WARN, FAIL, are "not assessed" (`NA`) and how many negative controls there are
(the counters exclude negative controls). Then the sortable table, one row per sample. The table is wide: it scrolls
horizontally inside its own box (the page itself never gets wider), the **Sample column stays in place** while you
scroll, and column headings wrap so that more columns fit:

| Column                   | Meaning                                                                                                                                                                                                                                                                                                                                                                     |
| ------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Sample, Type, Date, Site | From the samplesheet.                                                                                                                                                                                                                                                                                                                                                       |
| QC                       | Badge: green `PASS`, amber `WARN`, red `FAIL`, grey `NA`. See [QC rules](output.md#summaryqc_summarytsv) for how the status is decided.                                                                                                                                                                                                                                     |
| Reasons                  | The `qc_reasons` tokens: `low_coverage`, `amplicon_dropout`, `no_reads`, `no_depths`, `ntc_reads`, `primer_scheme_mismatch`, `low_primer_trimming` (the last two are underlined and explain themselves on hover).                                                                                                                                                           |
| Raw reads, Mapped %      | From fastp and `samtools flagstat` (`-` in `demix` mode without them).                                                                                                                                                                                                                                                                                                      |
| Mean depth               | Mean depth over the whole genome (zeros included).                                                                                                                                                                                                                                                                                                                          |
| >=10x %                  | % of the genome at depth >= 10.                                                                                                                                                                                                                                                                                                                                             |
| >=min depth %            | % of the genome at depth >= `--consensus_min_depth` (what will not be masked in the consensus).                                                                                                                                                                                                                                                                             |
| Amplicons dropped        | `dropped / total` amplicons.                                                                                                                                                                                                                                                                                                                                                |
| Demix                    | Badge for the demix status: `ok` (green), `low_confidence`, `skipped_low_coverage` and `missing` (amber; `low_confidence` is shown as `low confidence`, with the reason as its tooltip), `failed` (red). A `low_confidence` sample keeps its lineage bar and heatmap row, and is left out of the derived time series unless the run used `--include_low_confidence true`.   |
| Top group                | The most abundant summarized lineage group and its %.                                                                                                                                                                                                                                                                                                                       |
| Consensus coverage %     | % of the reference positions with a called base (A, C, G, T or an IUPAC code) in the consensus; the rest is `N` or deleted.                                                                                                                                                                                                                                                 |
| Primer trimmed %         | Only shown when at least one sample has a value. The share of mapped reads from which `ivar trim` removed primer sequence. **Informational**: calibration on the synthetic data showed it does not separate right from wrong primer schemes, so it is not a pass/fail signal; `low_primer_trimming` only appears when you set `--qc_min_primer_trimmed_pct`.                |
| Primer scheme            | Only shown when at least one sample has a value. Badge from the primer-scheme check: green `ok` (read starts line up with the declared scheme), amber `mismatch` (hover: the best-matching scheme; also a run-level flag), grey `uninformative` (read starts line up with no scheme, for example tagmentation libraries, so no check is possible), blank when not assessed. |

A negative control is judged differently from samples: it FAILS if it has real reads/coverage, and it never fails for
low coverage. `NEG_E` (5 mapped reads) is a PASS; `NEG_CONTAM` (thousands of reads of BA.2) is a FAIL. In the Demix
column a clean control such as `NEG_E` shows `skipped_low_coverage`, which is the expected outcome for a control.

## Lineage abundance

Per sample, what Freyja estimated. A two-button toggle selects the **resolution**:

- **Summarized groups** (default): Freyja's `summarized` grouping (`Omicron`, `Delta`, `XBB.1.5*`, `BA.2.86*`, ...).
  Easy to read for a surveillance overview. Note the grouping: JN.1 appears as `BA.2.86* (BA.2.86.X)`.
- **All lineages**: individual (sub)lineages from Freyja's `lineages`/`abundances` (`BA.2`, `AY.4`, `JN.1`, ...).

For each view the panel has a shared colour legend, a stacked bar per sample (hover for the exact value; segments below
0.2% of the bar are merged into grey "Other"), and a sample x lineage heatmap of % abundance (square-root colour scale;
at most the 20 most abundant groups or 25 most abundant lineages are shown as columns). Samples with no lineage result
show their demix status badge and reason instead of a bar (for example `skipped_low_coverage: less than 25% of the genome
covered at >= 10x`), and are counted under the heatmap. A negative control that demixed `ok` (`NEG_CONTAM`) shows an
`ok` badge without a bar and is also counted there: controls are deliberately not in the aggregate this panel is drawn
from (decision D8); their lineages are in [Negative controls](#negative-controls). The filter bar here works on the
sample name. Colours are shared only within a view.

With `--freyja_boot`, the tooltip of every bar segment and heatmap cell of the **All lineages** view ends with the 95%
bootstrap interval of that lineage (`freyja_boot_summary.tsv`, the `intervals` field of
[report_data.md](report_data.md)). The intervals are per individual lineage, so the **Summarized groups** view shows
none (only a group named like a lineage would pick one up). Intervals only appear in the tooltips; the bars themselves
are not drawn with error bars.

Read abundances as a relative composition of the lineages whose barcodes are in the dataset, not as absolute viral load;
low-abundance calls (a few percent) in mixtures are often noise.

## Coverage tracks

One small card per sample (QC badge in the corner): mean depth per bin of about 50 bp (at most 600 bins) over the
genome on a **log scale** (the note gives the top value). Overlays:

- dashed line: the amplicon dropout threshold (`--amplicon_min_depth`, 10x);
- dotted line: the consensus minimum depth (`--consensus_min_depth`, 10x by default);
- shaded bands: amplicons that were dropped (mean depth below the threshold).

Under each track: mean depth, % of the genome at 10x and the number of dropped amplicons. Use it to see at a glance where
a sample lost coverage (`Site_1_week_2` has two dropped amplicons; `LOWCOV_D` and `NEG_E` have all 29 dropped). The data
come from the depths files (full and `demix` mode).

### Multi-contig references (whole-genome influenza)

When the reference has several contigs (for example the eight influenza segments; see
[multi_contig_targets](wiring/multi_contig_targets.md)) the report adds, and a single-contig report stays byte-identical:

- under the title: "Freyja demixes `<contig>` only; other segments are QC/consensus only";
- the coverage cards run over all contigs on one axis, with alternating shading, a separator line at each boundary and a
  label strip under the track (the contig name, or its number in reference order when the names do not fit; hover for the
  name). The note says "across N bp in K segments". A bin on a boundary mixes two contigs;
- a **Segments** table under the cards: per sample and contig, length, mean depth, % of positions at 10x, consensus
  coverage %, consensus status and dropped amplicons. Amber rows have fewer than `--qc_pass_min_pct` % of positions at 10x
  or an empty consensus (the `segments_not_covered` run flag). Bad rows sort first;
- the amplicon dropout heatmap lists the amplicons of every contig as `<contig>/<name>` in genome order;
- the dropout impact section scores the Freyja contig only (a note says so);
- the QC table's "Mapped %" is the whole reference (all segments), as its tooltip says.

It needs `contig_coverage.tsv` and `consensus_stats_contigs.tsv` (`--contig-coverage`, `--consensus-stats-contigs`),
the `contig` column of the amplicon table and `contigs` / `freyja_contig` in `run_info.json`; without them the report still
builds from the depths files alone (segments in file order; no Segments table, and the header line does not name the
contig without `freyja_contig`).

## Amplicon dropout

A heatmap of samples (rows) x amplicons (columns, in genome order) coloured by mean depth: dark red = dropped with no
coverage, orange = dropped (< threshold), then light to dark blue for increasing depth (<100x, <1000x, <10000x,

> =10000x), grey = no data. The top bar shows, per amplicon, the fraction of non-control samples in which it dropped out:
> a tall bar means a systematic problem (primer failure, a mutation under a primer) rather than one bad sample. The
> amplicon number ruler is below it. Hover a cell for the amplicon and its mean depth. Needs a primer scheme (always in
> `full`; in `demix` only with `--primer_scheme`/`--primer_bed`). For `--protocol capture` (and `--skip_primer_trimming`)
> the section is replaced by the note "Not applicable for protocol capture": there are no amplicons. The header's
> "Library protocol" fact shows the protocol (and, for capture, that duplicates were marked); genome-wide depth is under
> Coverage tracks. In a run that mixes amplicon and capture samples (per-sample `protocol` column) the fact reads
> `mixed (N capture, M amplicon; ...)` and the amplicon section lists the amplicon samples only, with a note.

## Amplicon dropout impact

Only when `summary/dropout_impact.tsv` / `dropout_impact_summary.tsv` exist (full and demix mode, amplicon protocol,
`--dropout_impact`). A dropped amplicon hides every mutation inside it; what matters for the lineage estimate are the
**defining mutations** of the detected lineages: the mutations a lineage has that its parent lineage lacks (hover the word
in the report). Mutations inherited from ancestors do not tell a lineage from its parent, so hiding them does not move the
estimate.

- **Dropped amplicons table**: per dropped amplicon its span and mean depth, the number of **defining sites hidden**
  (primary), which lineages lose them (`XBB.1.5 (2), BA.2 (1)`, defining-site counts), then the _differential_ sites
  (carried by some but not all detected lineages; reported for context, not part of the warning) and the
  lineage-informative sites (carried by some but not nearly all lineages of the Freyja dataset; secondary context).
  Rows in amber with the badge "hides defining sites" mean a lineage detected at or above the affect threshold
  (`--dropout_min_affect_abundance`, default 5%) loses defining mutations, so its abundance estimate may be biased; such
  samples also raise the `dropout_affects_lineages` run flag. A lineage between 1% and the threshold is listed but does not
  highlight the row.
- **Per-sample summary** above it: defining sites lost / observable and the **% defining sites lost** (primary), the number
  of dropouts hiding defining sites, and the share of lineage-informative sites lost as secondary context.
- Samples without a lineage estimate (skipped, failed, controls) say "no lineage estimate". The tables are filterable
  and sortable like the others; long tables show the first 100 rows (500 at most) with a "Show all" button.
- Tables written before this refinement (without the `n_defining_sites`, `n_differential_sites` and
  `defining_sites_*` columns) still render, with the older informative-sites layout.

## Negative controls

Appears when at least one negative control (`sample_type=negative`) demixed `ok`. Negative controls are demixed when
they pass the coverage gate, so that a contaminated control is visible, but they are **excluded from the aggregate, the
lineage view and the time series** (decision D8; unless the run used `--include_negatives_in_aggregate`, in which case
the last column reads "included"). This section is where they are shown: a sortable table with one row per control:

| Column                      | Meaning                                                                                                                                  |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| Control                     | Sample ID with the `NTC` chip.                                                                                                           |
| QC, Reasons                 | The control's QC badge and reason tokens (`FAIL` / `ntc_reads` for a contaminated control).                                              |
| Mapped reads, Genome >=1x % | From the QC row: how much of the control mapped to the genome.                                                                           |
| Demix coverage %            | Freyja's genome coverage of the control's demix.                                                                                         |
| Lineages detected           | The five most abundant lineages with their % (`BA.2 100.0%`), with "..." if there are more. Taken from `freyja_aggregated_controls.tsv`. |
| Aggregate                   | `excluded` (default) or `included`.                                                                                                      |

In the screenshot above `NEG_CONTAM` (5,998 mapped reads, 99.4% of the genome at >= 1x) contains BA.2, the lineage of
the mixtures run beside it: the signature of carry-over contamination. A clean control that was skipped for low
coverage (`NEG_E`) has no row here (look at QC and Demix status). The run-level `ntc_reads` flag in the flags list is
what raises the alarm; this section tells you what the contamination is.

## Demix status

Counters (ok / skipped for low coverage / failed, plus "no result" when relevant) and a table with the status, the
reason, the barcode version used, Freyja's genome coverage % and the fit residual of each sample. This is the place
to find out why a sample has no lineage bar. Negative controls are listed too: an `ok` control (`NEG_CONTAM`) is `ok`
here, shown without coverage and residual because those come from the aggregate it is not part of (they are in the
Negative controls table).

## Mutations of interest

Only when `--mutations_of_interest` was given. A heatmap of the alternate-allele frequency (%) per sample and mutation.
Cells with a dashed outline have depth below the consensus minimum depth (treat with care); a dash (`-`) means no
coverage or no data. Hover for the allele counts ("alt 41 of 52 reads"). The values are those of
`summary/mutations_of_interest.tsv`.

## Unexplained mutations

Only when `summary/unexplained_mutations.tsv` exists (full and demix mode, `--unexplained_mutations`). High-frequency
single-nucleotide mutations of a sample that the demixed lineage mixture does not account for: for each site the
frequency **expected** from the detected lineages and their abundances is compared with the **observed** frequency, and
sites with enough depth, reads, frequency, excess and statistical support are listed. A mutation that is unexplained in
several samples (**recurrent**) is an early hint of a novel or cryptic lineage, or of a lineage missing from the barcodes.
**It is not a variant caller**: sequencing errors, mixed-in lineages that are not in the barcodes, low depth and amplicon
dropout all produce false hits, so check interesting rows on the reads. The table lists recurrent rows first with a
"recurrent (n)" badge and an amber highlight, gene and amino-acid change, observed and expected frequency (hover the
expected frequency for the excess), depth, alt reads and the p-value in scientific notation (a one-sided binomial tail
probability of the alt read count under the expected frequency; p-values are floored at 1e-300 and shown as
`<=1e-300`, meaning "at most"). A small **not in dataset** marker on a row (newer tables only, column `in_dataset`) means that
no lineage in the Freyja dataset models this mutation: a dataset blind spot (for example masked sites) or a novel change;
the info flag then says how many of the recurrent mutations are of this kind. The thresholds used are printed above the
table. When
the file exists but has no rows the section says "None found with the current thresholds"; when the analysis did not run
there is no section. A recurrent mutation also raises the info flag `unexplained_mutations`.

## Lineage groups

The lineage-group table (one column per group of **your** groups file, per sample) is deliberately **not** part of the
HTML report: it lives in `summary/lineage_summary.csv` (and, for site aggregates,
`summary/lineage_summary_site_aggregates.csv`), where it opens in a spreadsheet. The report only mentions it:

- **Provenance** (last section) has one line: `Lineage groups: <groups file name>, sha256 <first 12>, dataset <Freyja
dataset version>; table in summary/lineage_summary.csv`. The file name comes from `--lineage_groups` (otherwise
  "bundled default"); checksum and dataset version come from `summary/lineage_groups_resolved.tsv`, so the groups used
  can be identified later.
- **Run-level flags**: `lineage_groups_warnings` (warn) lists the first warnings of the groups audit (a pattern that is
  not a lineage of the dataset, a group that matches nothing, or counted groups that overlap), and
  `lineage_groups_ungrouped` (info) is raised when an `ok` sample has 2% or more of its lineages outside every counted
  group; it names the top ungrouped lineages and points to `summary/lineage_ungrouped_summary.tsv`: add them to your
  lineage groups file.

## Site aggregates

Only when `summary/site_aggregates_freyja.tsv` / `site_aggregates_meta.tsv` exist (`--site_aggregation`). Groups of
sites (a catchment made of several treatment plants) **combined from the per-site results, never from pooled reads**:
each site's samples are demixed on their own and the finished abundances of the member sites are averaged per period
(ISO week by default) with weights (viral load, a weight such as population served from your sites file, or equal).
The section has one stacked bar per aggregate-period, grouped by aggregate with the periods in order, drawn and
coloured like the per-sample bars (the same Summarized groups / All lineages toggle; lineages keep the colour they have
in the per-sample view). Chips beside each bar: **n/total sites** (green when every member site contributed, amber when
some are missing, red when none did), the **weighting basis** (amber when it is weaker than requested, e.g.
`equal (viral_load missing)`), the **share of the total weight** that contributed, one **missing &lt;site&gt;** chip per
missing site (hover for the reason: `not_ok:failed`, `no_sample`, `no_date`) and the number of samples. A period where no
site contributed shows its status instead of a bar. Incomplete periods raise the info flag `site_aggregate_incomplete`.

## Time series

Only when at least two distinct collection dates exist among non-control samples with a lineage result. Shows the mean
abundance of the top summarized groups (up to 6) per collection date, with one panel per `site` (samples with no site
form a "(no site)" panel) on shared axes. Hover a point for its value. Negative controls and samples without a lineage
result are excluded. (This is the report's own panel; `freyja_plot_times.pdf` is Freyja's time-aware plot.)

The panel needs `collection_date` (and optionally `site`) for the samples. In `full` and `demix` mode they come from
the QC table; in **`aggregate` mode** (no QC table) the report reads them from the samplesheet in effect: an `--input`
sheet with those columns, or the `summary/samplesheet.csv` carried in a `--previous_results` directory. With
`--demix_dir` there is no metadata, so no time series. (Verified: an `aggregate` run on a `full` run's results drew the
4-date panel; the `--demix_dir` run did not.)

**From the time-series table.** When `summary/freyja_timeseries.tsv` exists (with `--timeseries`, the default whenever
at least two collection dates exist) the panel is built from that table instead, with controls:

- **Resolution**: summarized groups, or all lineages when `freyja_timeseries_lineages.tsv` exists (the 30 most abundant
  lineages are kept; the top 6 are drawn).
- **Weighting**: unweighted (each sample counts equally) or **viral-load weighted** (each sample counts in proportion to
  its `viral_load`, so low-load samples move the curve less; dates without any recorded load are blank). The weighted
  button is disabled, with a tooltip, when the samplesheet has no `viral_load` column.
- **Site**: each site (default), one site, or all sites pooled (sites combined by sample count, or by total viral load when
  weighted).

Hover a point for the value, the number of samples averaged on that date and the mean viral load.

**Growth advantage.** With `--freyja_growth` a table from `freyja relgrowthrate` follows the chart: per lineage the
estimated relative growth advantage (%) per serial interval (default 5.5 days) against all the other lineages, and
Freyja's bootstrap interval. Read it with care: the estimate is `exp(serial interval x r) - 1` for the logistic growth
rate `r` of the lineage's share, the interval is on the growth-rate scale (so the estimate can lie outside it) and is a
rough indication of uncertainty rather than a confidence interval, all sites are pooled, and the fit needs several
dates spread over weeks (the pipeline writes a note instead of numbers when there are too few; the note's text is shown
as "Growth rates could not be computed", for example "insufficient data: 3 distinct collection dates ..."). If a Freyja
dashboard was written (`--freyja_dash`) its file name is shown as text; the dashboard is a separate file that loads its
styling from the internet, whereas this report is fully offline.

## Barcode comparison

Appears with `--compare_barcodes` (modes `full` and `demix`): every sample is demixed with a second Freyja barcode
dataset and the report compares the two estimates. A is the dataset of `--barcodes`, B the one of `--compare_barcodes`
(the header line names both barcode versions, here `custom` and `09_17_2026-00-46`). It shows:

- **Samples with the largest changes**: per sample, the _turnover_ (half the sum of the absolute abundance changes over
  all lineages: 0 = identical estimates, 100% = completely different) and the largest single change in percentage
  points; the 50 samples with the highest turnover.
- **Largest lineage changes**: the 100 biggest sample x lineage changes, with A and B as small bars, the percentage of
  each, and the change in points (orange = increased in B, blue = decreased); sortable.

In the screenshot, `Site_1_week_2` has the highest turnover (43.2%): its BA.5 falls from 60.7% to 29.1% and a BA.5.5
appears (29.1%) in the container's barcode set; `PURE_C` (JN.1) and the contaminated control are unchanged. Different
barcode sets can assign the same signal to different lineages, so a big turnover tells you the lineage names are
dataset-dependent for that sample, not that the sample changed. Negative controls that demixed are included in the
comparison. Samples skipped by the coverage gate have no rows (they are `skipped_low_coverage` in
`barcode_compare_status.tsv`); a second demix that fails raises the `barcode_compare_failed` flag. The data are in
`summary/barcode_compare.tsv` ([output.md](output.md#summarybarcode_comparetsv-and-summarybarcode_compare_statustsv)).

## Consensus genomes

Table per consensus: length, `N` count, `Coverage %` (positions with a called base), other (IUPAC) ambiguous bases and status badge (`ok` green, `empty` amber,
`failed` red). See [consensus conventions](output.md#consensus-fasta-conventions). Negative controls have none. The
FASTA itself is `summary/consensus_genomes.fasta`.

A **Reason** column (and a tooltip on the status badge) appears when `consensus_stats.tsv` has a `reason` for at least
one sample. For `failed` consensuses (the sequence is all-N) it says why, for example `bcftools consensus failed: The
fasta sequence does not match the REF allele at NC_045512.2:11287:`; for `ok` consensuses it carries warnings, for
example `skipped variant at NC_045512.2:22942 (overlaps another variant)`. Failed consensuses raise the warning flag
`consensus_failed` and consensuses with warnings the info flag `consensus_warnings` (see
[run-level flags](#run-level-flags)). Older 6-column tables without `reason` still load.

## Provenance

- **Freyja dataset**: the full `freyja_data` block (source, barcode version, sha256 of the barcode file and other dataset
  files, Freyja version, dataset path, download time); see [output.md](output.md#the-freyja_data-block) for each field and
  the difference between `freyja_reported_barcode_version` and `freyja_bundled_barcode_version`.
- **Previous runs**: for `demix`/`aggregate` reports, the earlier runs the results build on, with their barcode version
  (an amber `unknown` for legacy Snakemake results, which recorded none).
- Two folds: "Full run provenance" (the whole `run_info.json` without parameters) and "All pipeline parameters".

## What is available in each mode

| Panel                                            |             `full`             |                                     `demix`                                     |                                               `aggregate`                                                |      `report`       |
| ------------------------------------------------ | :----------------------------: | :-----------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------: | :-----------------: |
| Header, barcode box, run-level flags, provenance |              yes               |                                       yes                                       |                                                   yes                                                    | as the original run |
| Quality control                                  |              yes               |              yes (no mapped/raw reads without a fastp JSON / BAM)               |                                                    no                                                    | as the original run |
| Lineage abundance, demix status                  |              yes               |                                       yes                                       |                                    yes (samples found in the inputs)                                     | as the original run |
| Coverage tracks                                  |              yes               |                                       yes                                       |                                                    no                                                    | as the original run |
| Amplicon dropout                                 |              yes               |                   only with `--primer_scheme`/`--primer_bed`                    |                                                    no                                                    | as the original run |
| Mutations of interest                            | with `--mutations_of_interest` |                         with `--mutations_of_interest`                          |                                                    no                                                    | as the original run |
| Amplicon dropout impact                          |    with `--dropout_impact`     |                             with `--dropout_impact`                             |                                                    no                                                    | as the original run |
| Unexplained mutations                            | with `--unexplained_mutations` |                         with `--unexplained_mutations`                          |                                                    no                                                    | as the original run |
| Site aggregates                                  |   with `--site_aggregation`    |                            with `--site_aggregation`                            |                                        with `--site_aggregation`                                         | as the original run |
| Time series                                      |        when >= 2 dates         |                                 when >= 2 dates                                 | when >= 2 dates in the samplesheet in effect (`--input` or `--previous_results`; not with `--demix_dir`) | as the original run |
| Negative controls                                |  when a control demixed `ok`   |                                      same                                       |                              same (a control `.demix` found in the inputs)                               | as the original run |
| Barcode comparison                               |   with `--compare_barcodes`    |                            with `--compare_barcodes`                            |                                                    no                                                    | as the original run |
| Bootstrap intervals (tooltips)                   |      with `--freyja_boot`      | with `--freyja_boot` <!-- UNVERIFIED: demix mode with --freyja_boot not run --> |                                                    no                                                    | as the original run |
| Consensus genomes                                |   unless `--skip_consensus`    |                            unless `--skip_consensus`                            |                                                    no                                                    | as the original run |

Unavailable panels say "not available in mode '<mode>'". A report built by `--mode report` shows the mode label of
the original run (it re-renders the original `report_data.json`); in the verified run the re-rendered `report.html` was
byte-for-byte identical to the original.

## Sharing a report: anonymisation

A report normally contains local paths and the full command line (provenance block, "Key parameters"), the user name
that appears in them, and your sample IDs and site names. Before sending it to someone outside the lab, build the report
data with `--anonymize` (pipeline parameter `report_anonymize`, values `none`, `paths`, `full`):

- `paths`: file system paths become `.../<file name>`, home directories `~`, the user name `<user>`, e-mail addresses
  `<email>`; the command line keeps its shape with base names only. Sample IDs and site names stay.
- `full`: additionally samples become `S001`, `S002`, ... and sites `SITE1`, ... (site-aggregate names `AGG1`, ...), the
  same everywhere: QC table, charts, tooltips, flags, filters, links. The mapping is written to a separate file
  (`--anonymize-map-out sample_id_map.csv`), never into the report; keep it private.

The header shows a short "Anonymised for sharing" note. **Not anonymised**: collection dates, lineage names, barcode
and software versions, QC values and every other number, and the other files in `summary/`. See
[report_data.md](report_data.md#anonymisation) for the exact rules.

## Rebuilding the report

The report can be re-rendered from a finished run at any time, for example after upgrading the pipeline (new template),
without touching the results:

```bash
nextflow run charlesfoster/wasteline -profile docker \
    --mode report --results_dir results --outdir results_report
```

- `--results_dir` is the earlier results directory; it must contain `pipeline_info/run_info.json`. It is read only.
- `--outdir` must be a different directory; the new `summary/report.html` is written there (plus Nextflow's own
  `pipeline_info`).
- If `summary/report_data.json` exists in `--results_dir` it is reused and only the HTML is rendered again. If it is
  missing the data is rebuilt from the other `summary/` files and the depths files, then rendered.
- To rebuild with different QC thresholds you need to re-run the QC, so use `--mode demix --previous_results results`
  (cheap: no mapping) rather than `--mode report`.

Without Nextflow, the two scripts can be run directly (Python >= 3.11, stdlib plus `jinja2`; see
[report_data.md](report_data.md)). This was run on the results of the real `full` run above and reproduced the
published `report_data.json` (apart from the timestamp and the completion keys in `provenance`), the report data and
HTML, and with `--write-qc-tsv qc_summary.tsv` added a byte-identical `summary/qc_summary.tsv`:

```bash
python bin/build_report_data.py --run-info results/pipeline_info/run_info.json \
    --qc results/summary/qc_summary.tsv --amplicons results/summary/amplicon_coverage.tsv \
    --aggregated results/summary/freyja_aggregated.tsv --demix-status results/summary/demix_status.tsv \
    --consensus results/summary/consensus_stats.tsv --depths results/samples/*/freyja/*.depths \
    --barcode-compare results/summary/barcode_compare.tsv \
    --controls-aggregated results/summary/freyja_aggregated_controls.tsv \
    --samplesheet results/summary/samplesheet.csv --boot results/summary/freyja_boot_summary.tsv \
    --compare-status results/summary/barcode_compare_status.tsv \
    --mode full --out report_data.json
python bin/render_report.py --data report_data.json --template assets/templates/report.html.j2 --out report.html
```

The optional analyses add `--dropout-impact`/`--dropout-summary` (`summary/dropout_impact*.tsv`), `--unexplained`
(`summary/unexplained_mutations.tsv`), `--timeseries`/`--timeseries-lineages` (`summary/freyja_timeseries*.tsv`),
`--growth-rates`/`--growth-failed`/`--growth-inputs` (`summary/freyja_growth_*`), `--dashboard-name`,
`--lineage-summary`/`--lineage-groups-resolved` (`summary/lineage_*`; used for the provenance line and two flags only) and
`--site-aggregates`/`--site-aggregates-meta` (`summary/site_aggregates_*.tsv`); see
[report_data.md](report_data.md). Without them the report is identical to one built before they existed.

(`--depths` takes depths files, so let the shell expand the glob **unquoted**: the script does not expand a quoted
`'*.depths'` pattern itself, and then the report has no coverage tracks. A directory also works. Leave out the options
for files your run did not produce.)

## Known issues

Found while documenting 0.1.0dev (real runs on the synthetic dataset); they affect how the report reads, not the
result files. Two earlier issues are fixed and verified in the final run: an `ok` negative control is no longer shown
as `failed` (it keeps `ok` and has its own section), and `missing` rows in `demix_status.tsv` no longer raise an
`input_problems` warning (an `aggregate` run on a `full` run's results only raised `demix_missing`).

- The report in `--mode report` shows the original run's mode label.
- In `aggregate` mode on `--previous_results`, a negative control that was skipped for low coverage (`NEG_E`) has no
  usable `.demix` and is not in the discovered samplesheet, so the report lists it as an ordinary sample flagged
  `demix_missing` rather than as a control.
- With `--include_negatives_in_aggregate` the controls are in the main aggregate, so they also appear in the lineage
  view (the time series always leaves controls out); the Negative controls table then shows them as "included"
  <!-- UNVERIFIED: not run; derived from the report builder and template code -->.
