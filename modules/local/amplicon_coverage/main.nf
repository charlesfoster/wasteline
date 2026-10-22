process AMPLICON_COVERAGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(depths)
    path bed // primer BED (bundled scheme or --primer_bed)

    output:
    tuple val(meta), path("${prefix}.amplicon_coverage.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: e.g. --min-depth 10 (amplicon dropout threshold)
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    amplicon_coverage.py \\
        --depths ${depths} \\
        --bed ${bed} \\
        --sample ${meta.id} \\
        ${args} \\
        --output ${prefix}.amplicon_coverage.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def contig_column = task.ext.multi_contig ? '\\tcontig' : '' // several reference contigs (PLAN D25): a trailing contig column
    """
    printf 'sample\\tamplicon\\tstart\\tend\\tmean_depth\\tmin_depth\\tpct_ge_min_depth\\tdropped${contig_column}\\n' > ${prefix}.amplicon_coverage.tsv
    """
}
