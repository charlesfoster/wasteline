process MINIBWA_MEM {
    tag "${meta.id}"
    label 'process_high'

    // Experimental alternative to BWA_MEM (--aligner minibwa). The image is the existing BioContainers
    // mulled image for minibwa 0.7 + samtools 1.24, exact tag (verified with the quay.io API, see docs/aligner_evaluation.md)
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-842faba7afc332e008f7d8b210ce965379274249:371090825f1a2952ebfb0c8b2ea8dd7a86651fbe-0'
        : 'quay.io/biocontainers/mulled-v2-842faba7afc332e008f7d8b210ce965379274249:371090825f1a2952ebfb0c8b2ea8dd7a86651fbe-0'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    val sort_bam

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('minibwa'), eval("minibwa version | sed 's/-r[0-9]*\$//'"), emit: versions_minibwa, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view -b'
    """
    # The index of a small (viral) genome takes seconds, so it is built here instead of in a separate
    # index process: this keeps MINIBWA_MEM a drop-in for BWA_MEM(reads, fasta, sort_bam).
    # `minibwa map` is the native interface; with the default adaptive preset reads shorter than 325 bp are
    # treated as short reads. Unmapped reads are written (as bwa mem does); the pipeline filters them later.
    minibwa index -t ${task.cpus} ${fasta}

    minibwa map \\
        ${args} \\
        -t ${task.cpus} \\
        ${fasta} \\
        ${reads} \\
        | samtools ${samtools_command} ${args2} --threads ${task.cpus} -o ${prefix}.bam -
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
