# Dropout impact refinement and unexplained-mutation p-values: wiring notes

Two changes, both script/module/test only. Nothing outside `bin/dropout_impact.py`, `bin/unexplained_mutations.py`,
their modules, tests and fixtures was edited; the param, schema, report and docs changes below are for the wiring
agent. Supersedes the "affects_detected_lineage" definition in `docs/wiring/dropout_impact.md` (everything else there
still holds: signature, include line, failure handling, workflow sketch).

## 1. DROPOUT_IMPACT: defining sites instead of informative sites

### Why

With the old rule ("a detected lineage carries >= 1 informative site inside a dropped amplicon") nearly every dropout was
flagged: a Freyja barcode lists hundreds of mutations inherited from ancestors, so any dropped amplicon hides some
mutation of some detected lineage. On `2026-10-19_dev-test` 24 of 29 / 6 of 10 / 17 of 22 / 26 of 36 / 19 of 26 dropped
amplicons "affected detected lineages" and the run-level warn fired on every sample.

### New definitions (exact)

- **defining site of lineage l** = an *informative* site that l carries and its **parent** does not:
  `barcode(l) - barcode(parent(l))`, intersected with the informative set (carried by < 95% of all lineages).
  - Parent = `parent:` in the dataset's `lineages.yml`, by name or alias (case-insensitive). If the parent has no
    barcode in `lineage_mutations.json`, the **nearest ancestor with a barcode** is used (cycle-safe).
  - No parent / no usable ancestor (a root) = all its informative sites are defining.
  - **Recombinant roots** (entries with `recombinant_parents:`, e.g. `XFG`, `XGQ`, `XER`, `XFW`, `XBB`): the sites carried by any of their
    recombinant parents (`BA.1*` = `BA.1`; only parents with a barcode) are inherited, not defining. This goes one step beyond the
    brief ("recombinant roots: all mutations defining") because a detected `XGQ` (19.7%) and `XER` (13.9%) counted 129 and
    123 "defining" sites that are all inherited from their parents and made `WW-1290E` flag 18 of 36 amplicons (2 with
    this rule). `--strict-roots` restores the literal brief (default off).
  - No `lineages.yml` (a bare barcodes `.csv`, or `lineages.yml` absent): the parent is guessed from the name by
    dropping trailing `.N` components (`B.1.1.7 -> B.1.1 -> B.1 -> B`); aliased names (`BA.5`) then act as roots and
    over-count. A stderr warning says so. `lineages.yml` is looked up inside the dataset directory, or next to a dataset
    file; `--lineages-yml PATH` overrides.
- **detected lineage**: abundance `>= --min-abundance` (0.01), unchanged.
- **affects_detected_lineage** = `true` when a detected lineage with abundance `>= --min-affect-abundance` (new, default
  **0.05**, inclusive) has `>= 1` defining site in the dropped span. Lineages between 0.01 and 0.05 are still listed in
  `affected_lineages` and counted in `n_defining_sites`, but do not raise the flag.
- **n_differential_sites** = sites in the span carried by some but not all detected lineages (full barcodes; only when
  >= 2 are detected, else 0). **Reported, not part of the flag.** Reason (real data): with the 25 lineages
  detected at >= 0.01 a differential site exists in 22 of 29, 6 of 10, 12 of 22, 19 of 36 and 15 of 26 dropped amplicons,
  i.e. adding it to the rule would recreate the noise (most are mutations inherited by one clade and not another, and
  hundreds of other such sites remain visible). Defining sites are the evidence that separates a lineage from its parent
  and siblings, so hiding them is what makes the abundance move.

### Columns

New columns are **appended**; old columns keep name and position.

`dropout_impact.tsv` (one row per dropped amplicon):
`sample amplicon start end mean_depth n_informative_sites n_sites_in_detected_lineages affects_detected_lineage affected_lineages` **`n_defining_sites n_differential_sites`**

- `affected_lineages` = detected lineages with defining sites in the span, `name:defining_count`, max 5, e.g. `PJ.2.1:2,RE.1.1:1`
  (was: informative-site counts, so this column now lists far fewer lineages).
- `n_defining_sites` = distinct defining sites in the span, union over detected lineages; `n_differential_sites` as above.
  Both `NA` for samples without an aggregate row.

`dropout_impact_summary.tsv` (one row per sample):
`sample n_dropped_amplicons informative_sites_total informative_sites_in_dropped pct_informative_sites_lost n_dropped_affecting_detected_lineages` **`defining_sites_total defining_sites_lost pct_defining_sites_lost`**

- `n_dropped_affecting_detected_lineages` uses the new rule.
- `defining_sites_total` = distinct defining sites of the detected lineages (>= `--min-abundance`) inside the union of the
  sample's amplicon spans; `defining_sites_lost` = those inside the union of the dropped spans (overlaps once);
  `pct_defining_sites_lost` = 100 * lost / total, 2 decimals. All three `NA` without an aggregate row; the
  percentage is `NA` when total is 0.

### Wiring to do

`conf/modules.config`:

```groovy
    withName: 'DROPOUT_IMPACT' {
        ext.args = { "--min-abundance ${params.dropout_min_abundance} --min-affect-abundance ${params.dropout_min_affect_abundance}" }
        // rest unchanged
    }
```

| param | default | meaning |
|---|---|---|
| `dropout_min_affect_abundance` | `0.05` | a dropped amplicon is flagged only when it hides defining sites of a detected lineage at an abundance `>=` this |

`nextflow_schema.json`: number, min 0, max 1, next to `dropout_min_abundance` (description: "Minimum abundance of a
lineage for its defining sites in a dropped amplicon to raise the dropout flag"). Add the param to `nextflow.config`.
The `FREYJA_DATA` dataset directory already holds `lineages.yml` next to `lineage_mutations.json`, so no extra input is
needed; with a user `--barcodes` custom file there is no hierarchy and the name-based fallback (plus a warning) applies.
Docs: `docs/parameters.md` (the new param), `docs/output.md` (the 2 + 3 new columns and the new meaning of
`affects_detected_lineage`/`affected_lineages`/`n_dropped_affecting_detected_lineages`).

### Report suggestions (builder + template)

- `report_data.json` `dropout_impact[sample].summary`: add `defining_sites_total`, `defining_sites_lost`,
  `pct_defining_sites_lost` (NA -> null); `.amplicons[]`: add `n_defining_sites`, `n_differential_sites` (NA -> null).
  `affected_lineages` keeps `{lineage, n_sites}`, now defining counts.
- Table "Amplicon dropout impact": show *defining sites* and *differential sites* instead of (or next to) informative sites; keep the
  warn highlight on `affects_detected_lineage = true`. Caption per sample: "N of M defining sites of detected lineages lost
  (P%) in K dropped amplicons".
- Run-level flag `dropout_affects_lineages` (unchanged code, unchanged trigger: any sample with a true
  `affects_detected_lineage`). New message: *"Amplicon dropout hides defining mutations (those that separate a lineage from
  its parent) of lineages detected at >= 5% in N sample(s); abundance estimates for those lineages may be biased."*
  (print `thresholds.dropout_min_affect_abundance` instead of 5 and add it to `thresholds`.) The old text said
  "lineage-informative sites", which is no longer what the flag tests.
- Per-sample `dropout_pct_informative_lost` can stay; prefer `dropout_pct_defining_lost` (= `pct_defining_sites_lost`).

### Real-data validation (`2026-10-19_dev-test`, Freyja 2.0.6 dataset of 2026-10-06)

Dropped amplicons that affect detected lineages, old rule -> new rule (default `--min-affect-abundance 0.05`; in brackets: new
rule with `--strict-roots`):

| sample | dropped | before | after | defining sites lost / total (%) |
|---|---|---|---|---|
| WW-1285E | 29 | 24 | 12 (12) | 13 / 49 (26.53) |
| WW-1287E | 10 | 6 | 2 (2) | 3 / 18 (16.67) |
| WW-1289E | 22 | 17 | 5 (5) | 5 / 16 (31.25) |
| WW-1290E | 36 | 26 | 2 (18) | 4 / 18 (22.22) |
| WW-1292E | 26 | 19 | 0 (0) | 0 / 11 (0.00) |

The other seven samples (failed or negative control, 75-223 dropped amplicons) have no aggregate row and stay `NA`.
The whole analysis takes 0.2 s.

Hand checks (independent script on `lineage_mutations.json` + `lineages.yml`):

- `WW-1285E` `QIAseq_4` (706-917): 58 informative sites, 3 carried by detected lineages, **0 defining**. LF.7.11
  (parent LF.7), PJ.2.1 (PJ.2), RS.1.2 (RS.1), RE.1.1 (RE.1), MC.10.1.6 (MC.10.1) each carry one mutation in the span, all
  inherited from their parents. Old: flagged; new: not flagged.
- `WW-1285E` `QIAseq_5` (794-1029): PJ.2.1 carries 2 mutations, one of them `C934T`, which PJ.2 does not have: 1 defining
  site, PJ.2.1 is at 22.0% -> `affects = true`, `affected_lineages = PJ.2.1:1`.
- `WW-1285E` `QIAseq_17` (2391-2628): `C2470T` is defining for PJ.2.1, and RS.1 (2.1%) has one defining site of its own -> `PJ.2.1:1,RS.1:1`;
  flagged through PJ.2.1 (RS.1 alone would not flag).
- min-affect rule: `WW-1290E` `QIAseq_163` (21873-22159) hides 1 defining site of PQ.2.1.17, abundance 2.7% (< 0.05):
  listed in `affected_lineages`, `n_defining_sites = 1`, `affects = false`.

## 2. UNEXPLAINED_MUTATIONS: p-values in log space

`p_value` showed `0` for deep calls: `P(X >= alt) ~ 0.01^407` underflows a plain double. The tail is now summed in log space
(`binom_logsf`) and formatted by `format_pvalue`:

- 3 significant digits as before (`2.39e-42`);
- **floored at `1e-300`**: any p-value below 1e-300 is written `1e-300`, which means *"at most 1e-300"* (it cannot be told
  apart from a true 1e-300). The column is always a plain positive float (`float("1e-300")` works; no `0`, no `<`, no `NA`), so the
  report builder, sorting and concat need no change (concat passes the string through; its sort is on `excess`).
- the `--max-pvalue` guard compares in log space, so thresholds below 1e-308 now work and `--max-pvalue 0` rejects everything.

Real data: 53 of the 386 rows of `summary/unexplained_mutations.tsv` had `0`; they now read `1e-300` (the one row with
`1.99e-301` is floored too); every other value is unchanged. Document the floor in `docs/output.md` / the report footnote
("p-values are floored at 1e-300").

`--max-rows` (no pipeline param exists; script default 100): per sample, rows are sorted by excess (descending, ties by
position then alt) and cut at `--max-rows`, with a stderr warning ("N unexplained mutations; reporting the top K by excess").
Checked: exactly at the cap = no warning; one over = lowest excess dropped and the tie at the cut resolved by position;
`0` = header only (+ warning); negative is now refused (was silently dropping the last rows via slicing). Two consequences
worth a line in the docs: the cut is made *before* concat, so `n_samples`/`recurrent` count only rows that survived the cap
(real `WW-1285E`: 128 called, 100 reported); and ties on excess (many `excess = 1`) are cut by position, not by
evidence. Optional wiring: `unexplained_max_rows` (default 100) -> `--max-rows` in `ext.args`.

### S:G142D (position 21987 G>A) reported at 100% in four samples

Not a mixture artefact and not a bug: **no lineage of the dataset carries any mutation at position 21987** (0 of 5,977
lineages in `lineage_mutations.json`, including Delta and BA.1, which are usually described with G142D), so the expected frequency
is 0 for every possible mixture. The demixed lineages of `WW-1285E` (LF.7.11, PJ.2.1, RS.1.2, RE.1.1, MC.10.1.6, ...),
`WW-1287E`, `WW-1289E`, `WW-1292E` all lack it in their barcodes. The signal itself is clean: ivar shows ALT_DP =
TOTAL_DP (407/407, 62/63, 49/49, 139/139) while 21986 and 21988 are reference. So either all circulating lineages carry
G142D (the position is a blind spot of the Freyja 2.0.6 mutation set, e.g. dropped with the 21987-21995 deletion
encoding), or it is systematic; the table cannot tell. This is what the new `in_dataset` column records (below): all four rows
read `in_dataset = false`.

### New column `in_dataset` (appended LAST, after `p_value`)

- Per-sample table (14 columns): `sample position ref alt observed_freq expected_freq excess depth alt_depth gene aa_change effect p_value` **`in_dataset`**.
- Concat table (16 columns): `... effect n_samples recurrent p_value` **`in_dataset`** (n_samples/recurrent stay before p_value; sorting is unchanged and does not use it).
- `in_dataset` = `true` when at least one lineage of the dataset's barcodes (ANY lineage, not only the demixed ones) carries this exact
  (position, alt) SNV, else `false`. Works for `lineage_mutations.json` and for barcodes `.csv` (a column with no non-zero cell = not carried). It does not change
  the call, the p-value or the order.
- Meaning of `false`: the barcodes never model this mutation (a masked / problematic site such as 21987 G>A = S:G142D, which 0 of 5,977
  lineages carry, or a truly novel change); the table cannot tell which. Read a **recurrent hit with `in_dataset = false` first as a
  dataset blind spot, second as a candidate novel variant**.
- Real data (5 samples, 386 rows): 255 `true`, 131 `false`; 8 recurrent rows are `false`, at 3 sites (17506 C S5748P, 21987 G142D, 4092 I1276T).
- **Header check / wiring:** `--concat` accepts only the 14-column per-sample header (it exits 2 with a message saying the table predates `in_dataset` and must be
  regenerated). Old 13-column per-sample tables and old 15-column concat tables are rejected. The wiring needs no new input (the dataset is already an input of
  `UNEXPLAINED_MUTATIONS`), but any resumed/cached run mixing old and new per-sample files fails at `UNEXPLAINED_MUTATIONS_CONCAT`; rerun the per-sample step (`-resume` is
  invalidated automatically because the script changed). Anything that reads `unexplained_mutations.tsv` by position must read `p_value` as second-to-last now.
- **Report:** `report_data.json` unexplained rows get `in_dataset` (boolean; absent/`null` for older tables). Table: a small "not in dataset" badge or column. The INFO
  flag `unexplained_mutations` (recurrent unexplained mutations) message gains, when applicable, "(N of them at sites no lineage in the dataset models)", N = number of
  distinct recurrent (position, alt) with `in_dataset = false`; omit the parenthesis when N is 0. Suggested wording: *"Recurrent unexplained mutations in K sites
  (N of them at sites no lineage in the dataset models); recurrent hits at such sites are most likely a gap in the lineage barcodes, then candidate novel variants."*
- `docs/output.md` / `docs/report_data.md`: document the column and the 14/16-column layouts.
