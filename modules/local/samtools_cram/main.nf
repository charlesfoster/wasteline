process SAMTOOLS_CRAM {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${prefix}.cram"), path("${prefix}.cram.crai"), emit: cram
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus > 1 ? task.cpus - 1 : 1
    """
    # Reference-based CRAM: the reference is needed again to read the file (its MD5 is stored per
    # contig), so keep the FASTA that made the run. `fai` may be an empty list: it is then built here.
    [ -s ${fasta}.fai ] || samtools faidx ${fasta}

    samtools view \\
        --threads ${threads} \\
        --cram \\
        --reference ${fasta} \\
        ${args} \\
        -o ${prefix}.cram \\
        ${bam}

    samtools index -@ ${threads} ${prefix}.cram
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cram ${prefix}.cram.crai
    """
}
