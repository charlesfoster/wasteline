// Test workflow for mapping.nf.test: MINIBWA_MEM, then inspect the BAM. Not part of the pipeline.
include { MINIBWA_MEM } from '../main'

process INSPECT_BAM {
    tag "${meta.id}"
    container 'quay.io/biocontainers/mulled-v2-842faba7afc332e008f7d8b210ce965379274249:371090825f1a2952ebfb0c8b2ea8dd7a86651fbe-0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.inspect.txt")

    script:
    """
    echo "sort_order \$(samtools view -H ${bam} | sed -n 's/^@HD.*SO:\\([a-z]*\\).*/\\1/p')" > ${meta.id}.inspect.txt
    echo "total \$(samtools view -c ${bam})" >> ${meta.id}.inspect.txt
    echo "mapped \$(samtools view -c -F 4 ${bam})" >> ${meta.id}.inspect.txt
    echo "proper_pairs \$(samtools view -c -f 2 ${bam})" >> ${meta.id}.inspect.txt
    """
}

workflow MAPPING {
    take:
    ch_reads // [ meta, [ fastq_1, fastq_2 ] ]
    ch_fasta // [ meta, fasta ]

    main:
    MINIBWA_MEM(ch_reads, ch_fasta, true)
    INSPECT_BAM(MINIBWA_MEM.out.bam)

    emit:
    inspect = INSPECT_BAM.out
}
