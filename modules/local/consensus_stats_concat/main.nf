process CONSENSUS_STATS_CONCAT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path stats, stageAs: 'stats/*' // one <sample>.consensus_stats.tsv per sample (collected)

    output:
    path "consensus_stats.tsv", emit: tsv
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
