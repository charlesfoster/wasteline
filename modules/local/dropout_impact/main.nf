process DROPOUT_IMPACT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path amplicons  // summary/amplicon_coverage.tsv (QC_CONCAT)
    path aggregated // summary/freyja_aggregated.tsv (FREYJA_AGGREGATE_PLOT)
    path dataset    // resolved Freyja dataset: a directory with lineage_mutations.json (+ lineages.yml for the hierarchy), or a barcodes .csv

    output:
    path "dropout_impact.tsv",         emit: impact
    path "dropout_impact_summary.tsv", emit: summary
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: e.g. --min-abundance 0.01 --min-affect-abundance 0.05 --informative-max-fraction 0.95 --strict-roots
    def args = task.ext.args ?: ''
    def contig_arg = task.ext.contig ? "--contig ${task.ext.contig}" : '' // several reference contigs (PLAN D25): the amplicon table has a contig column; only the demixed contig is used
    """
    dropout_impact.py \\
        --amplicons ${amplicons} \\
        --aggregated ${aggregated} \\
        --dataset ${dataset} \\
        ${contig_arg} \\
        ${args} \\
        --out dropout_impact.tsv \\
        --summary-out dropout_impact_summary.tsv
    """

    stub:
    """
    printf 'sample\\tamplicon\\tstart\\tend\\tmean_depth\\tn_informative_sites\\tn_sites_in_detected_lineages\\taffects_detected_lineage\\taffected_lineages\\tn_defining_sites\\tn_differential_sites\\n' > dropout_impact.tsv
    printf 'sample\\tn_dropped_amplicons\\tinformative_sites_total\\tinformative_sites_in_dropped\\tpct_informative_sites_lost\\tn_dropped_affecting_detected_lineages\\tdefining_sites_total\\tdefining_sites_lost\\tpct_defining_sites_lost\\n' > dropout_impact_summary.tsv
    """
}
