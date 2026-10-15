process BAM_MAPPED_READS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    // The count lets the workflow route samples with no mapped reads around `ivar trim`, which
    // crashes (SIGFPE) on an empty BAM.
    tuple val(meta), path(bam), path(bai), path("${meta.id}.mapped_reads.txt"), emit: counted
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools view -c -F 4 ${bam} > ${meta.id}.mapped_reads.txt
    """

    stub:
    """
    echo 1 > ${meta.id}.mapped_reads.txt
    """
}
