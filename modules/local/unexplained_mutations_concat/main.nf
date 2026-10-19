process UNEXPLAINED_MUTATIONS_CONCAT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path tables, stageAs: 'unexplained/*' // one <sample>.unexplained_mutations.tsv per sample (collected)

    output:
    path "unexplained_mutations.tsv", emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # adds n_samples (samples with the same position + alt unexplained) and recurrent; headers are validated
    unexplained_mutations.py \\
        --concat unexplained/* \\
        ${args} \\
        --output unexplained_mutations.tsv
    """

    stub:
    """
    printf 'sample\\tposition\\tref\\talt\\tobserved_freq\\texpected_freq\\texcess\\tdepth\\talt_depth\\tgene\\taa_change\\teffect\\tn_samples\\trecurrent\\tp_value\\tin_dataset\\n' > unexplained_mutations.tsv
    """
}
