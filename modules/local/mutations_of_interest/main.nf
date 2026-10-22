process MUTATIONS_OF_INTEREST {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    // depths is optional ([]): without it depth/ref_depth of absent mutations are NA
    tuple val(meta), path(variants), path(depths)
    path moi        // CSV with columns name,mutation
    path gff        // GFF3 for amino-acid mutations when the TSV is not annotated; [] to omit
    path reference  // reference FASTA, used together with the GFF3; [] to omit

    output:
    tuple val(meta), path("${prefix}.mutations_of_interest.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def depths_arg = depths ? "--depths ${depths}" : ''
    def gff_arg = gff && reference ? "--gff ${gff} --reference ${reference}" : ''
    def contig_arg = task.ext.contig ? "--contig ${task.ext.contig}" : '' // several reference contigs (PLAN D25): the contig of the variants and depths
    """
    mutations_of_interest.py \\
        --variants ${variants} \\
        --moi ${moi} \\
        --sample ${meta.id} \\
        ${depths_arg} \\
        ${gff_arg} \\
        ${contig_arg} \\
        ${args} \\
        --output ${prefix}.mutations_of_interest.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tname\\tmutation\\tref_depth\\talt_depth\\tdepth\\talt_freq\\n' > ${prefix}.mutations_of_interest.tsv
    """
}
