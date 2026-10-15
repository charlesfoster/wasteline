process BARCODE_COMPARE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(demix_a), path(demix_b) // the same sample demixed with dataset A (primary) and dataset B (--compare_barcodes)
    path version_a, stageAs: 'version_a.json'     // barcode_version.json of dataset A
    path version_b, stageAs: 'version_b.json'     // barcode_version.json of dataset B

    output:
    tuple val(meta), path("${prefix}.barcode_compare.tsv"), emit: rows // data rows of summary/barcode_compare.tsv (no header)
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    barcode_compare.py \\
        --sample ${meta.id} \\
        --a ${demix_a} \\
        --b ${demix_b} \\
        --version-a ${version_a} \\
        --version-b ${version_b} \\
        --output ${prefix}.barcode_compare.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '${meta.id}\\tstub\\t0.500000\\t0.400000\\ta\\tb\\n' > ${prefix}.barcode_compare.tsv
    """
}
