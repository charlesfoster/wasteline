process UNEXPLAINED_MUTATIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    // demix is optional ([]): a sample without a demix (skipped / failed) gets a header-only table
    tuple val(meta), path(variants), path(demix)
    path dataset    // the Freyja dataset directory (lineage_mutations.json) or a custom barcodes CSV
    path gff        // GFF3 for gene / amino-acid annotation; [] to omit
    path reference  // reference FASTA, used together with the GFF3; [] to omit

    output:
    tuple val(meta), path("${prefix}.unexplained_mutations.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def demix_arg = demix ? "--demix ${demix}" : ''
    def gff_arg = gff && reference ? "--gff ${gff} --reference ${reference}" : ''
    """
    unexplained_mutations.py \\
        --variants ${variants} \\
        --dataset ${dataset} \\
        --sample ${meta.id} \\
        ${demix_arg} \\
        ${gff_arg} \\
        ${args} \\
        --output ${prefix}.unexplained_mutations.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tposition\\tref\\talt\\tobserved_freq\\texpected_freq\\texcess\\tdepth\\talt_depth\\tgene\\taa_change\\teffect\\tp_value\\tin_dataset\\n' > ${prefix}.unexplained_mutations.tsv
    """
}
