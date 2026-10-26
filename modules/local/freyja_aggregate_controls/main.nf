process FREYJA_AGGREGATE_CONTROLS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0'}"

    input:
    tuple val(meta), val(ids), path(demix_files, stageAs: 'demix/?.demix') // meta = the target (tkey); negative controls with status ok, and their .demix files, same order

    output:
    tuple val(meta), path("freyja_aggregated_controls.tsv"), emit: aggregated
    tuple val("${task.process}"), val('freyja'), eval('freyja --version | sed "s/.*version //"'), emit: versions_freyja, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def id_list = ids.collect { id -> "'${id}'" }.join(' ')
    """
    # Negative controls (PLAN D8) are demixed but kept out of freyja_aggregated.tsv and the plots. This
    # table, in the same format, shows what they contain; the report lists it under "Negative controls".
    # Labels are the sample IDs (see FREYJA_AGGREGATE_PLOT).
    mkdir labelled
    i=1
    for id in ${id_list}; do
        relabel_demix.py demix/\$i.demix --sample "\$id" --output "labelled/\$id.demix"
        i=\$((i + 1))
    done
    freyja aggregate labelled/ --ext demix --output aggregated.raw.tsv
    clean_aggregate_labels.py aggregated.raw.tsv freyja_aggregated_controls.tsv --sort
    """

    stub:
    def stub_rows = ids.collect { id -> "${id}\\tstub" }.join('\\n')
    """
    printf '\\tsummarized\\n${stub_rows}\\n' > freyja_aggregated_controls.tsv
    """
}
