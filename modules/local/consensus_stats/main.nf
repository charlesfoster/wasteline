process CONSENSUS_STATS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(fasta), path(status_json) // <sample>.consensus.fa and <sample>.consensus_status.json from CONSENSUS_BCFTOOLS
                                                    // (status_json may be [] : status/reason are then recomputed from the sequence)
    path reference               // reference FASTA: its length is the consensus_coverage_pct denominator

    output:
    tuple val(meta), path("${prefix}.consensus_stats.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.consensus_stats_contigs.tsv"), emit: contig_tsv, optional: true // several reference contigs only (PLAN D25): one row per contig
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The status JSON wins for failed/empty (a failed sample is an all-N sequence, which alone reads as "empty");
    // its reason (and the warnings text of an ok sample) goes to the last column, `reason`.
    def args = task.ext.args ?: ''
    def status_arg = status_json ? "--status-json ${status_json}" : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # total reference length (all records): the consensus_coverage_pct denominator. Several records (PLAN D25): the sample row is the
    # whole genome and the per-contig rows go to <sample>.consensus_stats_contigs.tsv
    ref_len=\$(awk '/^>/ { next } { total += length(\$0) } END { print total + 0 }' ${reference})
    contig_arg=""
    if [ "\$(grep -c '^>' ${reference})" -gt 1 ]; then contig_arg="--contig-out ${prefix}.consensus_stats_contigs.tsv"; fi

    # A missing/unparseable consensus is reported as status=failed (never a crash)
    consensus_stats.py \\
        --fasta ${fasta} \\
        --sample ${meta.id} \\
        --reference-length \$ref_len \\
        ${status_arg} \\
        \$contig_arg \\
        ${args} \\
        --output ${prefix}.consensus_stats.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tlength\\tn_count\\tconsensus_coverage_pct\\tambiguous_count\\tstatus\\treason\\n${meta.id}\\t0\\t0\\t100.00\\t0\\tok\\tNA\\n' > ${prefix}.consensus_stats.tsv
    """
}
