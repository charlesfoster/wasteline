process CONSENSUS_STATS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(fasta) // <sample>.consensus.fa from CONSENSUS_BCFTOOLS
    path reference               // reference FASTA: its length is the n_pct denominator

    output:
    tuple val(meta), path("${prefix}.consensus_stats.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ref_len=\$(awk '/^>/ { n++; next } n == 1 { total += length(\$0) } END { print total + 0 }' ${reference})

    # A missing/unparseable consensus is reported as status=failed (never a crash)
    consensus_stats.py \\
        --fasta ${fasta} \\
        --sample ${meta.id} \\
        --reference-length \$ref_len \\
        ${args} \\
        --output ${prefix}.consensus_stats.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tlength\\tn_count\\tn_pct\\tambiguous_count\\tstatus\\n${meta.id}\\t0\\t0\\t0.00\\t0\\tok\\n' > ${prefix}.consensus_stats.tsv
    """
}
