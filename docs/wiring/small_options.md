# Wiring notes: `--bam_format`, `--aligner`, `replay_run.py`

Three small features, all delivered as new files only. Nothing below is wired in yet.

| Feature             | New files                                                                                                                                                                  |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--bam_format cram` | `modules/local/samtools_cram/{main.nf,environment.yml,tests/main.nf.test,tests/main.nf.test.snap,tests/roundtrip.nf,tests/roundtrip.nf.test}`                              |
| `--aligner`         | `modules/local/minibwa_mem/{main.nf,environment.yml,tests/main.nf.test,tests/main.nf.test.snap,tests/mapping.nf,tests/mapping.nf.test}`, `tests/compare_aligners.sh`, `docs/aligner_evaluation.md` |
| replay              | `bin/replay_run.py`, `tests/python/test_replay_run.py`, `tests/fixtures/replay/**`                                                                                         |

## 1. CRAM output (`--bam_format`)

### Parameter

| Param        | Type   | Default | Enum         |
| ------------ | ------ | ------- | ------------ |
| `bam_format` | string | `'bam'` | `bam`, `cram` |

`nextflow.config` (typed params block, next to `save_bam`):

```groovy
    save_bam                   = true
    bam_format                 = 'bam'       // bam | cram: format of the published alignments (Freyja always works on the BAM)
```

`nextflow_schema.json`, group `read_trimming_and_bam_options`, after `save_bam`:

```json
"bam_format": {
    "type": "string",
    "default": "bam",
    "enum": ["bam", "cram"],
    "description": "Format of the published alignments (needs --save_bam).",
    "help_text": "`cram` publishes a reference-based CRAM and its `.crai` instead of the BAM and `.bai` (about a third of the size on the synthetic data: 23.6 kB vs 71.6 kB). Reading a CRAM needs the same reference FASTA (`samtools view --reference ref.fa`). Freyja and all QC steps still use the BAM, which stays in the work directory.",
    "fa_icon": "fas fa-file-archive"
}
```

### Include and signature

```groovy
include { SAMTOOLS_CRAM } from '../modules/local/samtools_cram/main'

// SAMTOOLS_CRAM(tuple val(meta), path(bam), path(bai), path fasta, path fai)
//   -> cram: tuple val(meta), path("<prefix>.cram"), path("<prefix>.cram.crai")
//   -> versions_samtools (topic: versions)
SAMTOOLS_CRAM(ch_final_bam, ch_fasta, [])     // [] = no .fai; it is built in the task (30 kb: instant)
```

`ch_fasta` is the existing `channel.value(file(params.fasta))` (a plain path). `[]` for `fai` is deliberate: the pipeline
has no FASTA index.

### Placement

After the final BAM (after `SAMTOOLS_MARKDUP` when duplicates are marked), as a leaf next to `FREYJA_VARIANTS`
(`workflows/wasteline.nf`, right after the `SAMTOOLS_MARKDUP` block):

```groovy
if (params.save_bam && params.bam_format == 'cram') {
    SAMTOOLS_CRAM(ch_final_bam, ch_fasta, [])
}
```

Nothing consumes its output, so Freyja, QC, consensus and the report keep using the BAM (a CRAM must never be fed to
`freyja variants`). Empty BAMs (clean negative controls) give an empty, indexed CRAM (tested), so no branching is needed.

### `conf/modules.config`

```groovy
    withName: 'SAMTOOLS_CRAM' {
        // the name follows the BAM it replaces: <id>.primertrimmed.cram, <id>.mapped.cram, <id>.*.markdup.cram
        ext.prefix = { (params.protocol == 'capture' || params.skip_primer_trimming) ? (marked ? "${meta.id}.mapped.markdup" : "${meta.id}.mapped") : (marked ? "${meta.id}.primertrimmed.markdup" : "${meta.id}.primertrimmed") }
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/bam" },
            mode: params.publish_dir_mode
        ]
    }
```

`marked` cannot be a local variable in a config closure; spell it out as the existing `SAMTOOLS_SORT` saveAs does:
`def marked = params.duplicates == 'mark' || (params.duplicates == 'auto' && params.protocol == 'capture')` inside the closure.
And stop publishing the BAM when CRAM is requested, in the two existing blocks:

```groovy
// SAMTOOLS_SORT saveAs:   params.save_bam && params.bam_format == 'bam' && !marked ? filename : null
// SAMTOOLS_MARKDUP saveAs: filename.endsWith('.markdup_stats.txt') || (params.save_bam && params.bam_format == 'bam') ? filename : null
```

The markdup stats file stays published in both formats. The `IVAR_TRIM` log is unchanged.

### Also needed

- `docs/output.md` (the `bam/` table rows and the `--save_bam` note): `bam/<s>.<...>.cram` + `.cram.crai` with `--bam_format cram`;
  `docs/parameters.md` is generated from the schema.
- `tests/check_outputs.py` greps no `.bam` names (checked), so the default path is untouched; add an nf-test (stub is enough,
  or `tests/data` with `--bam_format cram` asserting `*.cram`/`*.crai` exist and no `*.bam` is published).
- `bin/lock_containers.py --update`: `SAMTOOLS_CRAM` uses the same image as the other samtools modules; the lock must gain
  an entry for it (the stale check reports "SAMTOOLS_CRAM: not in the lock").
- `PLAN.md` step 4 ("CRAM is a future option") can be updated.

### Tests delivered

`nf-test test modules/local/samtools_cram --profile docker`: stub (snapshot); real on `tests/fixtures/markdup/MIX_A_1_4000.bam`
(2270 reads): CRAM and CRAI exist, CRAM smaller than the BAM (23.6 kB vs 71.6 kB), versions reported; header-only BAM gives an
empty indexed CRAM; `roundtrip.nf.test`: CRAM decodes (`samtools view --reference`) to 2270 reads whose name, flag,
contig, position, MAPQ, CIGAR, sequence and qualities equal the BAM's, and a region query through the `.crai` works.

## 2. Alternative aligner (`--aligner`)

Evidence and verdict: [../aligner_evaluation.md](../aligner_evaluation.md). `bwa` stays the default.

### Parameter

| Param     | Type   | Default | Enum              |
| --------- | ------ | ------- | ----------------- |
| `aligner` | string | `'bwa'` | `bwa`, `minibwa`  |

```groovy
    aligner                    = 'bwa'       // bwa | minibwa (experimental, see docs/aligner_evaluation.md)
```

Schema, group `read_trimming_and_bam_options` (or `reference_and_primer_options`):

```json
"aligner": {
    "type": "string",
    "default": "bwa",
    "enum": ["bwa", "minibwa"],
    "description": "Read aligner. `bwa` (bwa mem) is the default and the validated choice; `minibwa` is experimental.",
    "help_text": "minibwa 0.7 (Li and Homer, preprint 2026) is a young successor of bwa-mem. Its alignments are not identical to bwa mem's. Use `tests/compare_aligners.sh` to compare before relying on it. See docs/aligner_evaluation.md.",
    "fa_icon": "fas fa-align-left"
}
```

### Include and signature

```groovy
include { MINIBWA_MEM } from '../modules/local/minibwa_mem/main'

// MINIBWA_MEM(tuple val(meta), path(reads), tuple val(meta2), path(fasta), val sort_bam)
//   -> bam: tuple val(meta), path("<prefix>.bam")   (coordinate-sorted when sort_bam is true)
//   -> versions_minibwa, versions_samtools (topic: versions)
```

Deliberate difference from `BWA_MEM(reads, index, fasta, sort_bam)`: there is **no index input**. minibwa's index
(`.l2b`, `.mbw`) is built inside the task from the FASTA (instant for 30 kb), so no `MINIBWA_INDEX` module exists and the
call takes the same `ch_fasta_meta` the BWA path already has. The output is only a BAM (no `csi/crai/sam` emits).

### Placement (`workflows/wasteline.nf`, replacing lines 170-171 and the first input of `SAMTOOLS_VIEW`)

```groovy
def ch_aligned = channel.empty()
if (params.aligner == 'minibwa') {
    MINIBWA_MEM(FASTP.out.reads, ch_fasta_meta, true)
    ch_aligned = MINIBWA_MEM.out.bam
}
else {
    BWA_INDEX(ch_fasta_meta)
    BWA_MEM(FASTP.out.reads, BWA_INDEX.out.index, ch_fasta_meta, true)
    ch_aligned = BWA_MEM.out.bam
}

SAMTOOLS_VIEW(
    ch_aligned.map { meta, bam -> [ meta, bam, [] ] },
    ...
```

### `conf/modules.config`

```groovy
    withName: 'MINIBWA_MEM' {
        ext.prefix = { "${meta.id}" }
        // ext.args = '-x sr'   // optional: force the short-read preset instead of the adaptive one
    }
```

Nothing is published by the aligner (same as `BWA_MEM`).

### Also needed

- `bin/lock_containers.py --update` (the mulled image is an exact tag; `MINIBWA_MEM` needs a lock entry).
- Provenance: `tool_versions` gets `minibwa` automatically (topic `versions`); `containers.MINIBWA_MEM` is recorded at completion.
  `params.aligner` is recorded in `run_info.params`, so `replay_run.py` reproduces it.
- Tests in `tests/` (nf-test, tag `docker`): a run with `--aligner minibwa` on `tests/data` asserting the lineages of
  `truth.json` are still found. `tests/compare_aligners.sh <scratch>` is the manual comparison (assumes `--aligner bwa|minibwa`;
  exit 1 when an abundance differs by more than `TOLERANCE`, default 0.02). Its table logic was exercised on two copies of an
  existing full run (identical: pass; one abundance changed by 0.11: fail); the two real pipeline runs need the parameter.
- `docs/usage.md`/`parameters.md`: one line for `--aligner`, pointing to `docs/aligner_evaluation.md`. The output docs saying
  "bwa mem" (`docs/output.md` `bam/` rows) should say "bwa mem (or minibwa with `--aligner minibwa`)".

## 3. `bin/replay_run.py`

```
bin/replay_run.py <results dir | pipeline_info dir | run_info.json> [--execute] [--pipeline-dir DIR] [--outdir NEW]
                  [--db_dir DIR] [--no-locked] [--force] [--params-out FILE] [--no-nxf-ver] [--nextflow BIN]
```

Stdlib only, Python >= 3.8 syntax (host script; checked with `/usr/bin/python3` 3.9.6 and by the
`ast feature_version=(3, 8)` guard in `tests/python/test_replay_run.py`, which reuses `test_host_python_compat.py`;
**add `replay_run.py` to `HOST_SCRIPTS` in `tests/python/test_host_python_compat.py`** so the guard also lives there).

What it prints, in order: the pipeline version (`-r <commit>` from GitHub, or `git -C <dir> checkout <commit>` with `--pipeline-dir`,
HEAD compared; a null commit warns and falls back to `-r <version>`), the params file it wrote (default `./replay_params.json`;
recorded values minus nulls, `outdir`, barcodes, `db_dir`, previous results, volatile ones such as `trace_report_suffix`, `email*`,
bundled `fasta`/`gff`; with `--pipeline-dir` also minus the schema defaults), the barcodes pin, the input verification, the
`-profile`, and the new `--outdir`. `--execute` prints the command first and refuses on any hard error; `--force` overrides
only checksum/commit mismatches, never a missing file.

Barcodes: `latest`/`pinned` -> `--barcodes <iso>` if `<db_dir>/freyja/<iso>/` exists under `--db_dir`, the recorded `db_dir` or
`~/.wasteline/db` (barcode file sha256 checked); else the recorded `dataset_path` if it holds the dataset; else `container` if the
Freyja container bundles that very version (`freyja_bundled_barcode_version`); else an ERROR "MUST BE RESTORED" (this is stricter
than a warning: `--execute` is refused, because the pipeline would fail at `FREYJA_DATA` anyway). `container` -> `container`.
`custom` -> the recorded path with a sha256 check. `freyja_data_compare` is handled the same way. `report` mode has no
`run_info.json`, only `params_*.json`: the commit and profile are then unknown (warnings) and `--results_dir` is reproduced.

### Documentation to add

- `docs/usage.md`, "Reproducing the barcodes of an old run" and `docs/output.md`, "How to reproduce a run": a short step 0
  "`python3 bin/replay_run.py results --pipeline-dir <checkout>` prints the whole recipe below; add `--execute` to run it",
  with the caveats: it needs `--pipeline-dir` for schema-default trimming, it never overwrites the original (`<outdir>_replay`),
  stored Freyja datasets cannot be re-downloaded, uncommitted changes of the original checkout are not recorded.
- `docs/development.md`: one line in the helper scripts list (host Python, stdlib only, 3.8 syntax) and the fixtures
  (`tests/fixtures/replay/<mode>/{run_info,params}.json`, paths replaced by `@ROOT@`/`@REPO@`).
- `CHANGELOG.md`: "Added `bin/replay_run.py`, `--bam_format`, `--aligner` (experimental)".

### Tests delivered

`pytest tests/python/test_replay_run.py tests/python/test_host_python_compat.py` (46 tests, also green under `/usr/bin/python3` 3.9.6): every mode (full, capture, demix, demix legacy, aggregate with `--previous_results`, aggregate with `--demix_dir`,
report without `run_info.json`), every barcode source (latest, pinned, container, custom file and directory, stored dataset in
`--db_dir`/recorded db_dir/recorded `dataset_path`, missing dataset with and without the container bundling it, tampered
dataset), checksum mismatch of samplesheet/primer BED/reference/earlier results, missing files, null commit, commit
mismatch with a real git repo, schema-default trimming, locked profile, outdir protection, and `--execute` with a fake
`nextflow` (command printed first, refused on error, `--force`, exit status passed on).
The generated `-params-file` was not run through a real `nextflow run` (Docker was busy with a gate); the recorded
params files of real runs are known to work as `-params-file` (docs/output.md), and ours is a subset of one.
