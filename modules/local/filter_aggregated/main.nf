process FILTER_AGGREGATED {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(aggregated)               // meta = the target (tkey); freyja_aggregated.tsv
    path demix_status, stageAs: 'demix_status.tsv'  // summary/demix_status.tsv

    output:
    tuple val(meta), path("freyja_aggregated.filtered.tsv"), emit: aggregated // same format, without the samples of the excluded status (header only if none is left)
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: --exclude-status <status> (repeatable)
    def args = task.ext.args ?: '--exclude-status low_confidence'
    """
    filter_aggregated.py \\
        ${aggregated} \\
        --demix-status demix_status.tsv \\
        ${args} \\
        --output freyja_aggregated.filtered.tsv
    """

    stub:
    """
    head -n 1 ${aggregated} > freyja_aggregated.filtered.tsv
    """
}
