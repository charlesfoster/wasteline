process FREYJA_RELGROWTHRATE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0'}"

    input:
    path aggregated                          // freyja_aggregated.tsv
    path metadata, stageAs: 'metadata.csv'   // samplesheet-like CSV: sample, sample_type, collection_date, site and an optional numeric viral_load
    path dataset                             // Freyja dataset (lineages.yml is passed as --lineageyml if present)

    output:
    path "freyja_growth_rates.csv",        emit: growth_rates, optional: true // Freyja's rel_growth_rates.csv: Lineage, Estimated Advantage, Bootstrap 95% interval
    path "freyja_growth_rates.failed.txt", emit: failed,       optional: true // why there is no growth-rate table (insufficient data, Freyja error)
    path "freyja_growth_inputs.json",      emit: inputs                       // what went into the fit: n_samples, n_dates, n_dates_window, span, n_with_viral_load
    tuple val("${task.process}"), val('freyja'), eval('freyja --version | sed "s/.*version //"'), emit: versions_freyja, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: further `freyja relgrowthrate` options (--nboots, --serial_interval, --grthresh, --scale_by_viral_load, --config, --pathogen, ...)
    // ext.days (56), ext.mincov (60), ext.seed (1, bootstrap seed), ext.min_abundance (0.01, per-sample cutoff, see --thresh below),
    // ext.min_dates (4) and ext.min_span_days (14) are handled here.
    def args = task.ext.args ?: ''
    def days = task.ext.days ?: 56
    def mincov = task.ext.mincov ?: 60
    def seed = task.ext.seed ?: 1
    def cutoff = task.ext.min_abundance ?: 0.01
    def min_dates = task.ext.min_dates ?: 4
    def min_span = task.ext.min_span_days ?: 14
    def require_load = args.contains('--scale_by_viral_load') ? '--require-viral-load' : ''
    def own_thresh = args.contains('--thresh') ? 'yes' : 'no'
    """
    export MPLCONFIGDIR="\$PWD/.mplconfig"

    fail() {
        printf '%s\\n' "freyja relgrowthrate produced no growth-rate table: \$1" > freyja_growth_rates.failed.txt
        echo "WARNING: freyja_growth_rates.csv could not be produced, see freyja_growth_rates.failed.txt" >&2
    }
    tail_log() { tr '\\r' '\\n' < "\$1" | grep -v 'Rate calculations' | tail -n 8 | cut -c1-300; }
    json_get() { python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])' freyja_growth_inputs.json "\$1"; }

    # Everything below runs in a function so that a soft failure can `return` early: an `exit` would skip the
    # trailing version capture Nextflow appends to this script (the eval() output).
    run_all() {
    # Freyja wants an aggregate and a metadata file that list exactly the same samples (it fails on a
    # sample without metadata), so both are rebuilt from our samplesheet-like CSV: known, non-negative,
    # dated samples above --mincov. Freyja runs in fg/ because it drops scratch CSVs into its cwd.
    mkdir fg
    if ! timeseries.py --aggregated ${aggregated} --metadata metadata.csv --min-coverage ${mincov} \\
            --window-days ${days} ${require_load} \\
            --freyja-metadata-out fg/metadata.csv --freyja-aggregated-out fg/aggregated.tsv \\
            --freyja-summary-out freyja_growth_inputs.json 2> prep.log; then
        cat prep.log >&2
        echo '{}' > freyja_growth_inputs.json
        fail "could not prepare the inputs: \$(tail -n 1 prep.log)"
        return 0
    fi
    cat prep.log >&2

    # Freyja does not refuse a hopeless fit: with one date it prints its optimiser start value
    # (+73.3% at the default serial interval) for every lineage. So require enough dates first.
    N_SAMPLES=\$(json_get n_samples)
    N_DATES=\$(json_get n_dates_window)
    SPAN=\$(json_get window_span_days)
    if [ "\$N_SAMPLES" -eq 0 ]; then
        fail "no sample with a collection date, a non-negative sample type and coverage above ${mincov}% (see the samplesheet columns collection_date and sample_type)"
        return 0
    fi
    if [ "\$N_DATES" -lt ${min_dates} ] || [ "\$SPAN" -lt ${min_span} ]; then
        fail "insufficient data: \$N_DATES distinct collection date(s) spanning \$SPAN day(s) within the last ${days} days; need at least ${min_dates} dates spanning at least ${min_span} days"
        return 0
    fi

    # Freyja compares each sample's lineage abundance with --thresh x (number of samples) when it merges
    # minor lineages into 'Other', so the default --thresh 0.01 lumps everything below 0.3 together
    # with 30 samples. Unless --thresh is given explicitly, divide the per-sample cutoff by the sample count.
    THRESH_ARG=""
    if [ '${own_thresh}' = no ]; then
        THRESH_ARG="--thresh \$(python3 -c "import sys; print('%.6g' % (float(sys.argv[1]) / int(sys.argv[2])))" ${cutoff} "\$N_SAMPLES")"
    fi

    LINEAGEYML=\$(freyja_dataset_info.py ${dataset} --print-lineageyml-path)
    LINEAGEYML_ARG=""
    if [ -n "\$LINEAGEYML" ]; then
        LINEAGEYML_ARG="--lineageyml \$(python3 -c 'import os, sys; print(os.path.abspath(sys.argv[1]))' "\$LINEAGEYML")"
    fi

    # PLAN D6 class of bug: Freyja 2.0.x ships pandas 3 and the CLI dies with a dtype TypeError unless
    # pandas' infer_string is switched off before freyja._cli.cli() runs. Freyja has no seed option and
    # draws its bootstrap from numpy's global RNG, so the same wrapper also seeds it (reproducible CI).
    cat > fg/run_freyja.py <<'PY'
    import sys
    import numpy as np
    import pandas as pd
    if sys.argv[1] == "infer-string-off":
        pd.options.future.infer_string = False
    np.random.seed(int(sys.argv[2]))
    from freyja._cli import cli
    sys.argv = ["freyja"] + sys.argv[3:]
    cli()
    PY

    # usage: run_attempt <plain|infer-string-off> <log>
    run_attempt() {
        rm -f fg/rel_growth_rates.csv
        ( cd fg && python3 run_freyja.py "\$1" ${seed} relgrowthrate aggregated.tsv metadata.csv \\
            --days ${days} --mincov ${mincov} \$THRESH_ARG \$LINEAGEYML_ARG --output rel_growth_rates.csv ${args} ) > "\$2" 2>&1 || true
    }
    has_rows() { [ -s fg/rel_growth_rates.csv ] && [ "\$(wc -l < fg/rel_growth_rates.csv)" -gt 1 ]; }

    run_attempt plain growth_attempt1.log
    if ! has_rows; then
        run_attempt infer-string-off growth_attempt2.log
    fi
    if has_rows; then
        cp fg/rel_growth_rates.csv freyja_growth_rates.csv
    else
        {
            echo "freyja relgrowthrate produced no growth-rate table (plain attempt and infer_string workaround both failed or gave no rows;"
            echo "an empty table means no lineage reached --grthresh in the last ${days} days)."
            echo "--- attempt 1, last lines ---"
            tail_log growth_attempt1.log
            echo "--- attempt 2 (infer_string=False), last lines ---"
            tail_log growth_attempt2.log 2>/dev/null || true
        } > freyja_growth_rates.failed.txt
        echo "WARNING: freyja_growth_rates.csv could not be produced, see freyja_growth_rates.failed.txt" >&2
    fi
    }
    run_all
    """

    stub:
    """
    printf '%s\\n' 'Lineage,Estimated Advantage,Bootstrap 95% interval' 'stub,0.0%,"[0.00% , 0.00%]"' > freyja_growth_rates.csv
    printf '{"n_samples": 0, "n_dates": 0, "n_dates_window": 0}\\n' > freyja_growth_inputs.json
    """
}
