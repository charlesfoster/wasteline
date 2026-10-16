process READ_STARTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bam)    // the mapped BAM BEFORE primer trimming (trimming moves the 5' ends)

    output:
    tuple val(meta), path("${prefix}.read_starts.tsv"), emit: tsv
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // -F 2308 = unmapped (4) + secondary (256) + supplementary (2048). Soft clips are ignored: the 5' end is the
    // alignment start (forward) or the last aligned reference base (reverse = POS + reference length - 1,
    // reference length = sum of the M/D/N/=/X CIGAR operations). Coordinates are 1-based.
    """
    samtools view -F 2308 ${bam} \\
        | awk 'BEGIN { OFS = "\\t"; print "strand", "five_prime_pos" }
            {
                if (int(\$2 / 16) % 2 == 0) { print "+", \$4; next }
                cigar = \$6; len = 0
                while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
                    n = substr(cigar, 1, RLENGTH - 1) + 0; op = substr(cigar, RLENGTH, 1)
                    if (op == "M" || op == "D" || op == "N" || op == "=" || op == "X") len += n
                    cigar = substr(cigar, RLENGTH + 1)
                }
                print "-", \$4 + len - 1
            }' > ${prefix}.read_starts.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'strand\\tfive_prime_pos\\n' > ${prefix}.read_starts.tsv
    """
}
