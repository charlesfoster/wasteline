process FREYJA_DASH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(aggregated)        // meta = the target (tkey); its freyja_aggregated.tsv
    path metadata, stageAs: 'metadata.csv'   // samplesheet-like CSV: sample, sample_type, collection_date, site and an optional numeric viral_load
    path dataset                             // Freyja dataset (lineages.yml is passed as --lineageyml if present)
    val title                                // dashboard title (Freyja inserts it as HTML)
    val intro                                // introduction text (Freyja inserts it as HTML)
    path plot_config, stageAs: 'plot_config.yml' // optional Freyja plot config (--config): lineage groups and colours ([] = none)

    output:
    tuple val(meta), path("freyja_dashboard.html"),                  emit: dashboard,    optional: true // Plotly is embedded; Bootstrap CSS/JS and a logo are loaded from the internet when the page is opened
    tuple val(meta), path("freyja_dashboard_rel_growth_rates.csv"),  emit: growth_rates, optional: true // the dashboard's own growth-rate table
    tuple val(meta), path("freyja_dashboard.failed.txt"),            emit: failed,       optional: true // why there is no dashboard
    tuple val(meta), path("freyja_dashboard_inputs.json"),           emit: inputs                       // what went into the dashboard: n_samples, n_dates, n_with_viral_load
    tuple val("${task.process}"), val('freyja'), eval('freyja --version | sed "s/.*version //"'), emit: versions_freyja, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: further `freyja dash` options (--nboots, --serial_interval, --grthresh, --scale_by_viral_load,
    // --headerColor, --bodyColor, --config, --pathogen, ...)
    // ext.days (56), ext.mincov (60), ext.seed (1, bootstrap seed), ext.min_abundance (0.01, per-sample cutoff,
    // see --thresh below), ext.min_dates (4) and ext.min_span_days (14) are handled here.
    def args = task.ext.args ?: ''
    def days = task.ext.days ?: 56
    def mincov = task.ext.mincov != null ? task.ext.mincov : 60
    def seed = task.ext.seed ?: 1
    def cutoff = task.ext.min_abundance ?: 0.01
    def min_dates = task.ext.min_dates ?: 4
    def min_span = task.ext.min_span_days ?: 14
    def require_load = args.contains('--scale_by_viral_load') ? '--require-viral-load' : ''
    def own_thresh = args.contains('--thresh') ? 'yes' : 'no'
    def config_arg = plot_config ? '--config ../plot_config.yml' : '' // freyja runs in fg/
    def title_q = "'" + title.toString().replace("'", "'\"'\"'") + "'"
    def intro_q = "'" + intro.toString().replace("'", "'\"'\"'") + "'"
    """
    export MPLCONFIGDIR="\$PWD/.mplconfig"

    fail() {
        printf '%s\\n' "freyja dash produced no dashboard: \$1" > freyja_dashboard.failed.txt
        echo "WARNING: freyja_dashboard.html could not be produced, see freyja_dashboard.failed.txt" >&2
    }
    tail_log() { tr '\\r' '\\n' < "\$1" | grep -v 'Rate calculations' | tail -n 8 | cut -c1-300; }
    json_get() { python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])' freyja_dashboard_inputs.json "\$1"; }

    # Everything below runs in a function so that a soft failure can `return` early: an `exit` would skip the
    # trailing version capture Nextflow appends to this script (the eval() output).
    run_all() {
    # Same preparation as FREYJA_RELGROWTHRATE: Freyja needs the aggregate and the metadata to list the
    # same samples. It runs in fg/ because it drops scratch CSVs into its cwd.
    mkdir fg
    if ! timeseries.py --aggregated ${aggregated} --metadata metadata.csv --min-coverage ${mincov} \\
            --window-days ${days} ${require_load} \\
            --freyja-metadata-out fg/metadata.csv --freyja-aggregated-out fg/aggregated.tsv \\
            --freyja-summary-out freyja_dashboard_inputs.json 2> prep.log; then
        cat prep.log >&2
        echo '{}' > freyja_dashboard_inputs.json
        fail "could not prepare the inputs: \$(tail -n 1 prep.log)"
        return 0
    fi
    cat prep.log >&2

    # The dashboard fits growth rates too, so it has the same minimum-data requirement.
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
    # minor lineages into 'Other' (see FREYJA_RELGROWTHRATE); scale the per-sample cutoff unless --thresh is given.
    THRESH_ARG=""
    if [ '${own_thresh}' = no ]; then
        THRESH_ARG="--thresh \$(python3 -c "import sys; print('%.6g' % (float(sys.argv[1]) / int(sys.argv[2])))" ${cutoff} "\$N_SAMPLES")"
    fi

    LINEAGEYML=\$(freyja_dataset_info.py ${dataset} --print-lineageyml-path)
    LINEAGEYML_ARG=""
    if [ -n "\$LINEAGEYML" ]; then
        LINEAGEYML_ARG="--lineageyml \$(python3 -c 'import os, sys; print(os.path.abspath(sys.argv[1]))' "\$LINEAGEYML")"
    fi

    # `freyja dash` reads its title and introduction from files
    printf '%s' ${title_q} > fg/title.txt
    printf '%s' ${intro_q} > fg/intro.txt

    # PLAN D6 class of bug: pandas 3 in Freyja 2.0.x (see FREYJA_RELGROWTHRATE); the wrapper also seeds numpy
    # because the bootstrap uses its global RNG and Freyja has no seed option.
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
        rm -f fg/dashboard.html fg/dashboard_rel_growth_rates.csv
        ( cd fg && python3 run_freyja.py "\$1" ${seed} dash aggregated.tsv metadata.csv title.txt intro.txt \\
            --days ${days} --mincov ${mincov} \$THRESH_ARG \$LINEAGEYML_ARG --output dashboard.html ${config_arg} ${args} ) > "\$2" 2>&1 || true
    }

    run_attempt plain dash_attempt1.log
    if [ ! -s fg/dashboard.html ]; then
        run_attempt infer-string-off dash_attempt2.log
    fi
    if [ -s fg/dashboard.html ]; then
        cp fg/dashboard.html freyja_dashboard.html
        if [ -s fg/dashboard_rel_growth_rates.csv ]; then cp fg/dashboard_rel_growth_rates.csv freyja_dashboard_rel_growth_rates.csv; fi
    else
        {
            echo "freyja dash produced no HTML file (plain attempt and infer_string workaround both failed)."
            echo "--- attempt 1, last lines ---"
            tail_log dash_attempt1.log
            echo "--- attempt 2 (infer_string=False), last lines ---"
            tail_log dash_attempt2.log 2>/dev/null || true
        } > freyja_dashboard.failed.txt
        echo "WARNING: freyja_dashboard.html could not be produced, see freyja_dashboard.failed.txt" >&2
    fi
    }
    run_all
    """

    stub:
    """
    printf '<!DOCTYPE html><html><body>stub dashboard</body></html>\\n' > freyja_dashboard.html
    printf '%s\\n' 'Lineage,Estimated Advantage,Bootstrap 95% interval' > freyja_dashboard_rel_growth_rates.csv
    printf '{"n_samples": 0, "n_dates": 0, "n_dates_window": 0}\\n' > freyja_dashboard_inputs.json
    """
}
