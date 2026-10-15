process FREYJA_AGGREGATE_PLOT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0'}"

    input:
    tuple val(ids), path(demix_files, stageAs: 'demix/?.demix') // sample IDs and their .demix files, same order
    path dataset                                                // Freyja dataset (lineages.yml is used for lineage grouping if present)
    path times                                                  // optional: CSV `Sample,sample_collection_datetime` for the time-aware plot ([] = none)

    output:
    path "freyja_aggregated.tsv",         emit: aggregated
    path "freyja_plot.pdf",               emit: plot, optional: true
    path "freyja_plot.failed.txt",        emit: plot_failed, optional: true
    path "freyja_plot_times.pdf",         emit: plot_times, optional: true
    path "freyja_plot_times.failed.txt",  emit: plot_times_failed, optional: true
    tuple val("${task.process}"), val('freyja'), eval('freyja --version | sed "s/.*version //"'), emit: versions_freyja, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_times = task.ext.args2 ?: ''
    def id_list = ids.collect { id -> "'${id}'" }.join(' ')
    def has_times = times ? 'yes' : 'no'
    """
    export MPLCONFIGDIR="\$PWD/.mplconfig"

    # Label every sample with its real ID. Freyja would use the variants file name recorded inside
    # each .demix file (MIX_A.variants.tsv from wasteline, MIX_A.freyja_variants.tsv or a path from
    # legacy results); relabelling a COPY makes every source layout give identical labels, with
    # underscores and everything else exactly as the sample ID.
    mkdir labelled
    i=1
    for id in ${id_list}; do
        relabel_demix.py demix/\$i.demix --sample "\$id" --output "labelled/\$id.demix"
        i=\$((i + 1))
    done

    # freyja aggregate needs a directory path WITH a trailing slash. Rows are then sorted by label
    # (the glob order inside Freyja is arbitrary) so the table is reproducible.
    freyja aggregate labelled/ --ext demix --output aggregated.raw.tsv
    clean_aggregate_labels.py aggregated.raw.tsv freyja_aggregated.tsv --sort

    LINEAGEYML=\$(freyja_dataset_info.py ${dataset} --print-lineageyml-path)
    LINEAGEYML_ARG=""
    if [ -n "\$LINEAGEYML" ]; then LINEAGEYML_ARG="--lineageyml \$LINEAGEYML"; fi

    # Stopgap (PLAN D6): `freyja plot` is broken in Freyja 2.0.3-2.0.5 containers (pandas 3
    # TypeError) and also exits 0 without a PDF when no sample exceeds --mincov. Success means
    # "a PDF exists". Attempt 1: plain freyja plot. Attempt 2: same CLI with pandas'
    # infer_string disabled. Otherwise record why and carry on.
    # usage: run_plot <aggregated.tsv> <output.pdf> <failed-note.txt> <log prefix> [extra freyja plot args]
    run_plot() {
        local agg="\$1" out="\$2" failed="\$3" log="\$4"
        shift 4
        freyja plot "\$agg" --output "\$out" \$LINEAGEYML_ARG "\$@" > "\${log}_attempt1.log" 2>&1 || true
        if [ ! -s "\$out" ]; then
            rm -f "\$out"
            python3 -c "
    import sys
    import pandas as pd
    pd.options.future.infer_string = False
    from freyja._cli import cli
    sys.argv = ['freyja', 'plot'] + sys.argv[1:]
    cli()
    " "\$agg" --output "\$out" \$LINEAGEYML_ARG "\$@" > "\${log}_attempt2.log" 2>&1 || true
        fi
        if [ ! -s "\$out" ]; then
            rm -f "\$out"
            {
                echo "freyja plot produced no PDF (plain attempt and infer_string workaround both failed)."
                echo "--- attempt 1 (freyja plot), last lines ---"
                tail -n 8 "\${log}_attempt1.log"
                echo "--- attempt 2 (infer_string=False), last lines ---"
                tail -n 8 "\${log}_attempt2.log" 2>/dev/null || true
            } > "\$failed"
            echo "WARNING: \$out could not be produced, see \$failed" >&2
        fi
    }

    run_plot freyja_aggregated.tsv freyja_plot.pdf freyja_plot.failed.txt plot ${args}

    # Time-aware plot (`freyja plot --times`) when collection dates are known. Freyja needs a date
    # for EVERY sample it plots, so only the samples with a date go into this plot.
    if [ '${has_times}' = yes ]; then
        head -n 1 freyja_aggregated.tsv > aggregated.dated.tsv
        awk -F'\t' 'NR == FNR { if (FNR > 1) { split(\$0, f, ","); keep[f[1]] = 1 } next } FNR > 1 && (\$1 in keep)' ${times ?: ''} freyja_aggregated.tsv >> aggregated.dated.tsv
        run_plot aggregated.dated.tsv freyja_plot_times.pdf freyja_plot_times.failed.txt plot_times --times ${times ?: ''} ${args} ${args_times}
    fi
    """

    stub:
    // the stub table lists the sample IDs it was given, so tests can see which samples reached the aggregate
    def stub_rows = ids.collect { id -> "${id}\\tstub" }.join('\\n')
    """
    printf '\\tsummarized\\n${stub_rows}\\n' > freyja_aggregated.tsv
    touch freyja_plot.pdf
    ${times ? 'touch freyja_plot_times.pdf' : ''}
    """
}
