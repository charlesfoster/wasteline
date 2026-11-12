# SOP: combining several wastewater sites into one estimate (site aggregation)

**Scope.** Why wasteline combines sites the way it does, how to do it, how to read the result, and how it differs from running
`freyja aggregate` on the same samples. Reference material lives elsewhere: options in [parameters.md](parameters.md), the sites
file and output columns in [usage.md](usage.md#site-aggregation), the experiment behind the design in
[wiring/site_aggregation.md](wiring/site_aggregation.md).

**Audience.** Anyone who has per-site wastewater samples and wants one lineage picture for a group of sites (a city, a region, a
catchment cluster) per week or day.

## 1. The idea in one paragraph

Every sample is demixed on its own, exactly as it would be without any aggregation. Site aggregation then takes the **finished
abundance results** of the member sites and computes a weighted average of them, per group and per period. Reads, variant tables and
depth files of different sites are never merged. If you remember one sentence: _we average answers, we do not pool evidence._

## 2. Philosophy

1. **Each sample is analysed alone.** A sample's lineage estimate does not depend on which other samples are in the run, or on
   whether it belongs to a group.
2. **An aggregate may not contain anything that no member site reported.** A weighted mean of non-negative numbers, where a lineage a
   site lacks counts as 0, cannot create a lineage. A single fit on pooled data can (section 6).
3. **Weights are explicit, chosen by a person, and written next to every result.** Sequencing depth is deliberately not a weight
   option: it measures how much was sequenced, not how much virus a site contributes.
4. **Missing data is shown, not hidden.** A week in which a site has no usable sample still gets a row, but the row says which site
   is missing and why, and how much of the group's weight is in it.
5. **Nothing is renormalised and nothing is smoothed.** The aggregate is plain arithmetic on the per-site results, so anyone can
   recompute a row by hand.

## 3. What goes in

| Input                            | Needed   | Notes                                                                                                                                 |
| -------------------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| Samples with a `site`            | yes      | The `site` column of the samplesheet. Matching to the sites file ignores case and extra spaces (`Site 1` = `site 1`; `Site1` is not). |
| Samples with a `collection_date` | yes      | ISO `YYYY-MM-DD`. A sample without a valid date cannot be placed in a week and is left out.                                           |
| A **sites file** (CSV)           | yes      | `aggregate_name,site_name[,weight]`, one row per site per group. A site can be in several groups.                                     |
| `viral_load` in the samplesheet  | optional | Positive number, any consistent unit. Used for weights if every contributing site has one that week.                                  |
| `weight` in the sites file       | optional | A positive number, typically population served. Used only when the weights basis says so.                                             |

A sample contributes to an aggregate only if its demix status is `ok`, it is not a negative control, and it has a site and a date.
`low_confidence` samples are excluded unless you pass `--include_low_confidence`.

## 4. How the average is computed

For each group and each period (ISO week, Monday to Sunday, or `day`):

1. **Replicates first.** The samples of one site in that period are combined into one site result, weighted by their `viral_load` if
   every one has a load, else equally. A site with two samples therefore never counts twice.
2. **Then sites.** The site results are averaged with the weights of `--site_aggregation_weights`:
   - `viral_load`: the site's mean load;
   - `site_weight`: the `weight` column of the sites file;
   - `equal`: 1 for every site;
   - `auto` (default): `viral_load` if every contributing site has one, else `site_weight` if every contributing site has one, else
     `equal`. The basis actually used is written per row, so it can differ between weeks.
3. For every lineage: `abundance = sum(w_i * a_i) / sum(w_i)` over the contributing sites. The same formula is applied to Freyja's
   `summarized` groups, `resid` and `coverage`.
4. Periods with fewer contributing sites than `--site_aggregation_min_sites` (default 1) are omitted.

If you have no weights, which is common, leave the `weight` column out and keep `auto`: you get an equal-weight mean, and
`weight_basis` says `equal`. See the caution in section 8.

## 5. Toy example

Two sites, `A` and `B`, in a group `Metro`. Three lineages. Week 1, then week 2.

**Week 1.** Site A has two samples (replicates); site B has one.

| Sample | Site | BA.2 | JN.1 | XBB  |
| ------ | ---- | ---- | ---- | ---- |
| A-1    | A    | 0.80 | 0.20 | 0    |
| A-2    | A    | 0.60 | 0.40 | 0    |
| B-1    | B    | 0.20 | 0.60 | 0.20 |

_Step 1, replicates._ Neither A sample has a viral load, so they are weighted equally: site A is (0.70, 0.30, 0). Site B is
unchanged (0.20, 0.60, 0.20).

_Step 2, sites, equal weights_ (no loads or weights given):

| Lineage | A    | B    | Metro week 1                 |
| ------- | ---- | ---- | ---------------------------- |
| BA.2    | 0.70 | 0.20 | (0.70 + 0.20) / 2 = **0.45** |
| JN.1    | 0.30 | 0.60 | (0.30 + 0.60) / 2 = **0.45** |
| XBB     | 0    | 0.20 | (0 + 0.20) / 2 = **0.10**    |

XBB appears in the aggregate only because site B has it, at half B's value (B's share of the weight). No lineage can appear that
neither site has.

_The same week with viral loads._ Say the mean loads are A = 100 and B = 300 gene copies per litre. Weights are 100/400 = 0.25 and
300/400 = 0.75, so B counts three times as much:

| Lineage | Calculation               | Metro week 1 |
| ------- | ------------------------- | ------------ |
| BA.2    | 0.25 × 0.70 + 0.75 × 0.20 | **0.325**    |
| JN.1    | 0.25 × 0.30 + 0.75 × 0.60 | **0.525**    |
| XBB     | 0.25 × 0 + 0.75 × 0.20    | **0.150**    |

The weights change the answer, which is why the basis is recorded in the output.

**Week 2.** Only site A has a usable sample (site B was not sampled). With `--site_aggregation_min_sites 1` (the default) the row is
produced from site A alone, and the meta file says so:

| row_label | n_sites_included | sites_included | sites_missing | missing_reasons | weight_basis | weight_fraction_included |
| --------- | ---------------- | -------------- | ------------- | --------------- | ------------ | ------------------------ |
| Metro_W1  | 2                | A;B            |               |                 | equal        | 1.0                      |
| Metro_W2  | 1                | A              | B             | no_sample       | equal        | 0.5                      |

Week 2 looks like a metro result but is really site A. Use `--site_aggregation_min_sites 2` to leave such weeks out. A site missing
because its sample was `low_confidence` shows `not_ok:low_confidence`; one without a collection date shows `no_date`.

## 6. Why not pool the reads, and what `freyja aggregate` is (and is not)

Three different things get called "aggregating" samples. They are not interchangeable.

| Approach                                                | What it does                                                                              | Combines sites?                    |
| ------------------------------------------------------- | ----------------------------------------------------------------------------------------- | ---------------------------------- |
| **`freyja aggregate`**                                  | Collects the `.demix` files of many samples into **one table, one row per sample**.       | No. The rows stay separate.        |
| **Pooling** (merge reads or variant tables, demix once) | Sums the evidence of several sites and fits one mixture.                                  | Yes, implicitly weighted by depth. |
| **Site aggregation** (this SOP)                         | Takes each sample's finished result and computes a weighted average per group and period. | Yes, with explicit weights.        |

### `freyja aggregate` is bookkeeping

`freyja aggregate` takes a directory of `.demix` outputs and writes `freyja_aggregated.tsv`: for each sample, its lineage list,
abundances, `resid` and `coverage`. It is how a batch of results is put in one file for plotting and downstream tools. It does not
average anything and does not know that sample `WW-1271E` and sample `WW-1270E` come from the same city.

Wasteline runs `freyja aggregate` in every run. Site aggregation **starts from its output**: `bin/site_aggregate.py` reads
`freyja_aggregated.tsv` and adds what `freyja aggregate` does not do: grouping by site, grouping by week, QC filtering, weighting
and an audit trail.

If someone says "just run `freyja aggregate` on the Sydney samples", what they get is a table of per-sample rows. To turn that
into a metro estimate they must then average by hand (or in a spreadsheet), and each of these decisions is made implicitly:

- which samples count (negative controls? `low_confidence` ones? samples with no date?);
- which week each belongs to;
- how replicates are combined, so that a site with two samples does not count twice;
- how sites are weighted (usually an unweighted mean of whatever rows happen to be there);
- what happens when a site is missing that week, which silently changes what the average means.

Site aggregation makes each of those a stated rule, applies it the same way every time, and writes the outcome per row.

### Why not pool the reads or variants?

Merging the data of several sites and demixing once has two problems:

1. **It weights each site by sequencing depth, not by how much virus it contributes.** Toy illustration: a lineage-defining mutation is
   at frequency 0.70 in site A (sequenced to 2000x) and 0.20 in site B (200x). The pooled frequency is about
   (0.70 × 2000 + 0.20 × 200) / 2200 ≈ **0.65**, while the equal-weight average is 0.45. Site A dominates only because it was
   sequenced more deeply.
2. **A single fit to a mixture of mixtures can report lineages that no site contains.** Averaging results cannot.

In the experiment behind the design (real Freyja 2.0.6, four finished test samples treated as four sites, equal and unequal depths),
the result-level aggregate stayed within 0.025 absolute abundance of the weighted truth. The pooled fit was further off when depths
differed, and in all eight pooled fits it reported a lineage that no site contains. That was synthetic data with exactly known,
linear lineage signatures, so it shows the direction and rough size of the effect, not what real wastewater gives.

## 7. Procedure

### 7.1 Prepare

1. Make sure every sample has a `site` and a `collection_date` in its samplesheet (`bin/make_samplesheet.py --metadata` can add
   them).
2. Write the sites file. Save as **CSV UTF-8**. A plain two-column file is enough:

   ```csv
   aggregate_name,site_name
   Sydney_Metro,Bondi
   Sydney_Metro,Quakers Hill
   ```

   `aggregate_name` may contain letters, digits, `_` and `-`. It becomes part of the row label (`Sydney_Metro_2026-W41`).

3. Decide the weights (section 8) and whether single-site weeks should be kept.

### 7.2 Run

The sites file is checked in the first minutes of the run; a mistake names the line and stops the run.

**From reads, in one run:**

```bash
nextflow run charlesfoster/wasteline -profile docker \
    --mode full --input samplesheet.csv \
    --site_aggregation sites.csv --outdir results
```

**From results that already exist** (no remapping; one or more earlier runs):

```bash
nextflow run charlesfoster/wasteline -profile docker \
    --mode aggregate --previous_results Run1,Run2 \
    --site_aggregation sites.csv --outdir results_metro
```

Useful options: `--site_aggregation_period week|day` (default `week`), `--site_aggregation_weights
auto|viral_load|site_weight|equal` (default `auto`), `--site_aggregation_min_sites N` (default 1), `--include_low_confidence`.

### 7.3 Check before you use a row

Open `summary/site_aggregates_meta.tsv` and read, for every row:

1. `n_sites_included`, `sites_included`, `sites_missing`, `missing_reasons`: did the sites you expected contribute?
2. `weight_basis`: was it the basis you intended? (`auto` can change from week to week.)
3. `weight_fraction_included`: what share of the group's weight is in the row (`NA` when it cannot be computed).
4. `n_samples`: how many samples went in.

A row with a missing site or a weaker basis than requested is not an error, but it must be read with that in mind. Also read the
run log for the count of samples left out for lacking a site or date, and any warning that a site name in the sites file matched no
sample (usually a typo).

### 7.4 Outputs

| File                                          | Contents                                                                                 |
| --------------------------------------------- | ---------------------------------------------------------------------------------------- |
| `summary/site_aggregates_freyja.tsv`          | One row per group and period in Freyja's aggregate format (`<aggregate_name>_<period>`). |
| `summary/site_aggregates_meta.tsv`            | What went into each row (section 7.3).                                                   |
| `summary/lineage_summary_site_aggregates.csv` | The lineage-group columns for each aggregate row, when the lineage summary is on.        |
| `summary/report.html`                         | Includes a "Site aggregates" section.                                                    |

## 8. Choosing weights

| Situation                                   | Recommendation                                                                                                                         |
| ------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| Viral load measured for every site and week | `viral_load` (or `auto`). Weights follow the amount of virus entering each site. Use concentration times flow if you want total virus. |
| Population served is known, no viral load   | `site_weight` with the `weight` column. A proxy: it assumes similar per-person shedding and dilution.                                  |
| Neither                                     | `equal`. State plainly that the aggregate is an average across **sites**, not across the **population**.                               |

`auto` picks the best basis available for each week, so the basis can change across weeks. If the aggregate is plotted as a time
series, pin the basis (`--site_aggregation_weights equal`, for example) so that a trend is not produced by a change of method.

With equal weights a small site counts as much as a large one. That is fine for "the typical lineage picture across these sites" and
misleading for "the lineage mix across the people in the metro area".

## 9. Limits

- **Validated lightly.** The design rests on one synthetic experiment. No real multi-site dataset has been checked against a known
  answer.
- **A point estimate.** Per-sample bootstrap intervals are not carried into the aggregate, so it has no confidence interval.
- **`resid` and `coverage` are averaged for convenience.** They are not statistically meaningful at group level.
- **Only as good as the member results and weights.** A biased site result is averaged in; a poor weight basis shifts the answer.
- **Not a Freyja method.** It is this pipeline's design. Freyja has no built-in way to combine sites, and Freyja's own tools have not
  been checked against the output file.
- **Single-site weeks.** Without `--site_aggregation_min_sites 2` a week with one contributing site is reported under the group's
  name.

## 10. Quick reference

```
sample result ──► (replicates of a site combined) ──► site result ──► (sites weighted-averaged) ──► group result per period
      ▲                                                                                  ▲
 freyja demix, one sample at a time                                  weights: viral_load | site_weight | equal
```

- Pool reads? **No.** Average results.
- Run `freyja aggregate`? **Yes, automatically**: it is the input table. It does not combine sites.
- Check `site_aggregates_meta.tsv` before trusting a row.
