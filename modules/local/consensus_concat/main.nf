process CONSENSUS_CONCAT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path fastas, stageAs: 'fasta/*' // one <sample>.consensus.fa per sample (collected)

    output:
    path "consensus_genomes.fasta", emit: fasta
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Multi-FASTA sorted by sample ID (the FASTA header; byte order, same as the TSV summaries).
    # Each input holds the sample's consensus: one record whose header is the sample ID, or (several reference contigs, PLAN D25)
    # one record per contig headed <sample ID>|<contig>, which stay together in reference order.
    for f in fasta/*; do
        if [ "\$(grep -c '^>' "\$f")" -lt 1 ]; then
            echo "ERROR: \$f does not contain a FASTA record" >&2
            exit 1
        fi
        printf '%s\\t%s\\n' "\$(grep -m1 '^>' "\$f" | sed 's/^>//; s/[[:space:]].*//; s/|.*//')" "\$f"
    done | LC_ALL=C sort -t "\$(printf '\\t')" -k1,1 | cut -f2 | while read -r f; do
        cat "\$f"
    done > consensus_genomes.fasta
    """

    stub:
    """
    touch consensus_genomes.fasta
    """
}
