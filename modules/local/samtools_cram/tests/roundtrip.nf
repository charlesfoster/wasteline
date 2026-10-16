// Test workflow for roundtrip.nf.test: BAM -> SAMTOOLS_CRAM -> decode both and write the columns that
// must agree. Not part of the pipeline.
include { SAMTOOLS_CRAM } from '../main'

process DECODE_BOTH {
    tag "${meta.id}"
    container 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'

    input:
    tuple val(meta), path(bam), path(bai), path(cram), path(crai)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.bam.txt"), path("${meta.id}.cram.txt"), path("${meta.id}.counts.txt")

    script:
    """
    # read name, flag, contig, pos, mapq, cigar, sequence, qualities: everything but the optional tags
    samtools view ${bam} | cut -f1-11 > ${meta.id}.bam.txt
    samtools view --reference ${fasta} ${cram} | cut -f1-11 > ${meta.id}.cram.txt
    echo "bam \$(samtools view -c ${bam})" > ${meta.id}.counts.txt
    echo "cram \$(samtools view -c --reference ${fasta} ${cram})" >> ${meta.id}.counts.txt
    # the .crai must be usable: a region query returns reads
    echo "cram_region \$(samtools view -c --reference ${fasta} ${cram} NC_045512.2:1-4000)" >> ${meta.id}.counts.txt
    """
}

workflow ROUNDTRIP {
    take:
    ch_bam // [ meta, bam, bai ]
    fasta  // path

    main:
    SAMTOOLS_CRAM(ch_bam, fasta, [])
    DECODE_BOTH(ch_bam.join(SAMTOOLS_CRAM.out.cram), fasta)

    emit:
    decoded = DECODE_BOTH.out
}
