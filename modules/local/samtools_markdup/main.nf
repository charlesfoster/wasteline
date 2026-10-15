process SAMTOOLS_MARKDUP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.markdup_stats.txt"),              emit: stats
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.markdup"
    def threads = task.cpus > 1 ? task.cpus - 1 : 1
    """
    # PCR/optical duplicates are only FLAGGED (0x400), never removed (no -r): samtools mpileup, which
    # `freyja variants` runs, excludes DUP-flagged reads from depth and variant counts by default.
    # markdup needs mate scores (fixmate -m) on a coordinate-sorted BAM: name-sort -> fixmate -m -> sort.
    samtools sort -@ ${threads} -n -u -T ${meta.id}.namesort ${bam} \\
        | samtools fixmate -@ ${threads} -m -u - - \\
        | samtools sort -@ ${threads} -u -T ${meta.id}.possort - \\
        | samtools markdup -@ ${threads} -s -f ${meta.id}.markdup_stats.txt - ${prefix}.bam

    samtools index -@ ${threads} ${prefix}.bam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.markdup"
    """
    touch ${prefix}.bam ${prefix}.bam.bai ${meta.id}.markdup_stats.txt
    """
}
