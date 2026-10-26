process CONSENSUS_STATS_CONCAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(stats, stageAs: 'stats/*') // meta = the target (tkey); one <sample>.consensus_stats.tsv per sample of the target

    output:
    tuple val(meta), path("consensus_stats.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    consensus_stats.py --concat stats/* --output consensus_stats.tsv
    """

    stub:
    """
    touch consensus_stats.tsv
    """
}
