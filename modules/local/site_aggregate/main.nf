process SITE_AGGREGATE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path aggregated                           // freyja_aggregated.tsv
    path demix_status, stageAs: 'demix_status.tsv'  // summary/demix_status.tsv: only status ok samples contribute
    path samplesheet, stageAs: 'samplesheet.csv'    // summary/samplesheet.csv: sample, sample_type, collection_date, site, optional viral_load
    path sites, stageAs: 'sites.csv'                // user file: aggregate_name, site_name, optional weight

    output:
    path "site_aggregates_freyja.tsv", emit: aggregated // Freyja aggregate format, row label <aggregate_name>_<period>
    path "site_aggregates_meta.tsv",   emit: meta       // row_label, aggregate_name, period, period_start, n_sites_included, sites_included, sites_missing, missing_reasons, weight_basis, weight_fraction_included, n_samples
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: --period week|day, --weights auto|viral_load|equal|site_weight, --min-sites <n>
    def args = task.ext.args ?: ''
    """
    site_aggregate.py \\
        --aggregated ${aggregated} \\
        --demix-status demix_status.tsv \\
        --samplesheet samplesheet.csv \\
        --sites sites.csv \\
        ${args} \\
        --out site_aggregates_freyja.tsv \\
        --meta-out site_aggregates_meta.tsv
    """

    stub:
    """
    printf '\\tsummarized\\tlineages\\tabundances\\tresid\\tcoverage\\n' > site_aggregates_freyja.tsv
    printf 'row_label\\taggregate_name\\tperiod\\tperiod_start\\tn_sites_included\\tsites_included\\tsites_missing\\tmissing_reasons\\tweight_basis\\tweight_fraction_included\\tn_samples\\n' > site_aggregates_meta.tsv
    """
}
