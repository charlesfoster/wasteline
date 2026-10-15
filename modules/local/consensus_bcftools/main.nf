process CONSENSUS_BCFTOOLS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.23--h3a4d415_0'
        : 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(mask) // CONSENSUS_FILTER_VCF `vcf`, LOW_DEPTH_MASK `bed`
    path fasta                                         // reference

    output:
    tuple val(meta), path("${prefix}.consensus.fa"),          emit: fasta
    tuple val(meta), path("${prefix}.consensus_status.json"), emit: status
    tuple val(meta), path("${prefix}.consensus.log"),         emit: log
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // PLAN 5.11: bcftools consensus --mark-del '-' -H I -m <low-depth mask>, FASTA header = sample ID.
    // A sample with no usable coverage (every position masked) or a failed bcftools run yields an
    // all-N sequence of reference length plus a status flag. Never an empty file, never a crash.
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ref_len=\$(awk '/^>/ { n++; next } n == 1 { total += length(\$0) } END { print total + 0 }' ${fasta})
    mask_bp=\$(awk '{ total += \$3 - \$2 } END { print total + 0 }' ${mask})
    n_var=\$(bcftools view -H ${vcf} | wc -l)
    # span of the called variants: CONSENSUS_FILTER_VCF removed it from the low-depth mask, so it does not count as coverage
    keep_bp=\$(bcftools query -f '%POS0\\t%END\\n' ${vcf} | awk '{ total += \$2 - \$1 } END { print total + 0 }')

    write_n_sequence() {
        awk -v id='${meta.id}' -v len="\$ref_len" 'BEGIN {
            print ">" id
            line = sprintf("%60s", ""); gsub(/ /, "N", line)
            for (i = 0; i < int(len / 60); i++) print line
            if (len % 60 > 0) print substr(line, 1, len % 60)
        }' > ${prefix}.consensus.fa
    }
    write_status() {
        reason=\$(printf '%s' "\$2" | tr -d '"\\\\' | tr '\\n\\t' '  ')
        printf '{"sample": "%s", "status": "%s", "reason": "%s"}\\n' '${meta.id}' "\$1" "\$reason" > ${prefix}.consensus_status.json
    }

    if [ "\$ref_len" -le 0 ]; then
        echo "reference has no sequence" >&2
        exit 1
    elif [ "\$((ref_len - mask_bp))" -le "\$keep_bp" ]; then
        echo "all \$ref_len positions are below the minimum depth (\$keep_bp bp are called variants)" > ${prefix}.consensus.log
        write_n_sequence
        write_status empty "no position reaches the minimum consensus depth (no usable coverage)"
    elif bcftools consensus \\
            -f ${fasta} \\
            --mark-del '-' \\
            -H I \\
            -m ${mask} \\
            ${args} \\
            ${vcf} > ${prefix}.consensus.raw.fa 2> ${prefix}.consensus.log \\
        && [ "\$(grep -c '^>' ${prefix}.consensus.raw.fa)" -eq 1 ] \\
        && [ "\$(grep -v '^>' ${prefix}.consensus.raw.fa | tr -d '\\n' | wc -c)" -gt 0 ]; then
        awk -v id='${meta.id}' '/^>/ { print ">" id; next } { print }' ${prefix}.consensus.raw.fa > ${prefix}.consensus.fa
        if [ "\$n_var" -eq 0 ]; then
            write_status ok "no variant passed the consensus filters; sequence equals the reference outside masked regions"
        else
            write_status ok ""
        fi
    else
        write_n_sequence
        write_status failed "bcftools consensus failed: \$(grep -v '^\\s*\$' ${prefix}.consensus.log | head -n 1)"
    fi
    rm -f ${prefix}.consensus.raw.fa
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>${meta.id}\\nNNNNNNNNNN\\n' > ${prefix}.consensus.fa
    printf '{"sample": "${meta.id}", "status": "ok", "reason": ""}\\n' > ${prefix}.consensus_status.json
    : > ${prefix}.consensus.log
    """
}
