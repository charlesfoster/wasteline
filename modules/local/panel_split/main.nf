process PANEL_SPLIT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    // bam: the sample's reads mapped to panel.fasta (any order). meta is the PER-TARGET meta (sampleMeta(..., target)): the target
    // to cut out is meta.target.name. meta is passed through.
    // The BAM is staged in a subdirectory: it is usually named <id>.bam, the name of the output, and `samtools sort -o` writes
    // through a symlink (it would overwrite the upstream task's BAM in its work directory, and the splits of one sample race).
    tuple val(meta), path(bam, stageAs: 'panel_mapped/*')
    path contig_map

    output:
    // same shape as BAM_MAPPED_READS.out.counted: the BAM holds only this target's contigs under their ORIGINAL names (the legacy
    // single-target reference naming), coordinate-sorted and indexed; the count is the number of mapped reads in it.
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), path("${meta.id}.mapped_reads.txt"), emit: counted
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def target = meta.target.name
    // Not interpolated (''' string). Primary records on the target's contigs are kept (secondary/supplementary are dropped upstream);
    // @SQ lines of other targets are dropped; contig names go back to the original; a mate that lies on another target (a discordant
    // fragment) is flagged as an unmapped mate, as samtools writes it, so the record does not point outside the trimmed header.
    def split_awk = '''
        BEGIN { FS = OFS = "\\t" }
        FNR == NR {
            if (FNR > 1 && $2 == wanted) keep[$1] = $3
            next
        }
        $1 == "@SQ" {
            name = ""
            for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) name = substr($i, 4)
            if (!(name in keep)) next
            line = $1
            for (i = 2; i <= NF; i++) line = line OFS (($i ~ /^SN:/) ? "SN:" keep[name] : $i)
            print line
            next
        }
        $1 ~ /^@(HD|RG|PG|CO)$/ { print; next }
        {
            if (!($3 in keep)) next
            flag = $2 + 0
            $3 = keep[$3]
            if ($7 != "=" && $7 != "*") {
                if ($7 in keep) $7 = keep[$7]
                else {
                    $7 = "="
                    $8 = $4
                    $9 = 0
                    if (int(flag / 8) % 2 == 0) flag += 8
                    if (int(flag / 2) % 2 == 1) flag -= 2
                    if (int(flag / 32) % 2 == 1) flag -= 32
                }
            }
            $2 = flag
            print
        }
    '''
    """
    samtools view -h -F 0x900 ${bam} \\
        | awk -v wanted='${target}' '${split_awk}' ${contig_map} - \\
        | samtools sort -o ${prefix}.bam -
    samtools index ${prefix}.bam
    samtools view -c -F 4 ${prefix}.bam > ${meta.id}.mapped_reads.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    echo 1 > ${meta.id}.mapped_reads.txt
    """
}
