process LINEAGE_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path aggregated, stageAs: 'freyja_aggregated.tsv'               // summary/freyja_aggregated.tsv
    path demix_status, stageAs: 'demix_status.tsv'                  // summary/demix_status.tsv: every sample gets a row; [] when no sample needs one (aggregate mode)
    path samplesheet, stageAs: 'samplesheet.csv'                    // summary/samplesheet.csv (site, collection_date, sample_type); [] allowed
    path groups, stageAs: 'groups.csv'                              // the user's lineage groups CSV (--lineage_groups)
    path dataset, stageAs: 'dataset'                                // Freyja dataset directory holding lineages.yml (+ last_barcode_update.txt)
    path aggregated_sites, stageAs: 'site_aggregates_freyja.tsv'    // optional site aggregates in Freyja aggregate format; [] when absent
    path aggregated_sites_meta, stageAs: 'site_aggregates_meta.tsv' // optional: its meta TSV; [] when absent

    output:
    path "lineage_summary.csv", emit: summary                                              // one row per sample; one column per group
    path "lineage_groups_resolved.tsv", emit: resolved                                     // audit: what each group matched in the dataset
    path "lineage_groups_used.csv", emit: used                                             // byte-identical copy of the groups file used: a starting point for the next edit
    path "lineage_ungrouped_summary.tsv", emit: ungrouped                                  // lineages in no in_total group (`other`), biggest first, with a suggested pattern
    path "lineage_summary_site_aggregates.csv", optional: true, emit: sites                // only with the site-aggregate inputs
    path "lineage_summary.xlsx", optional: true, emit: xlsx                                // formatted workbook; only when ext.args carries `--xlsx-out lineage_summary.xlsx`
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: --major-threshold <fraction> --min-abundance <fraction> --decimals <n> --zeros blank|zero
    //           [--xlsx-out lineage_summary.xlsx [--xlsx-decimals <n>] [--wasteline-version <v>]]
    //           (the xlsx is written from the same tables as the CSVs; xlsx_writer.py ships next to lineage_groups.py in bin/)
    def args = task.ext.args ?: ''
    def status_arg = demix_status ? '--demix-status demix_status.tsv' : ''
    def sheet_arg = samplesheet ? '--samplesheet samplesheet.csv' : ''
    def sites_arg = aggregated_sites ? '--aggregated-sites site_aggregates_freyja.tsv --sites-out lineage_summary_site_aggregates.csv' : ''
    def sites_meta_arg = aggregated_sites_meta ? '--aggregated-sites-meta site_aggregates_meta.tsv' : ''
    """
    lineage_groups.py \\
        --aggregated freyja_aggregated.tsv \\
        ${status_arg} \\
        --groups groups.csv \\
        --dataset dataset \\
        ${sheet_arg} \\
        ${sites_arg} \\
        ${sites_meta_arg} \\
        ${args} \\
        --output lineage_summary.csv \\
        --resolved-out lineage_groups_resolved.tsv \\
        --ungrouped-out lineage_ungrouped_summary.tsv

    cp groups.csv lineage_groups_used.csv
    """

    stub:
    // an empty placeholder is enough for a stub run (no python needed); the real run writes a valid workbook
    def xlsx_stub = (task.ext.args ?: '').contains('--xlsx-out') ? 'touch lineage_summary.xlsx' : ''
    def sites_stub = aggregated_sites ? "printf 'aggregate_name,period,period_start,n_sites_included,sites_included,sites_missing,weight_basis,other,proportion_summary,major_lineages,major_abundances,minor_lineages,minor_abundances,ungrouped_lineages\\n' > lineage_summary_site_aggregates.csv" : ''
    """
    printf 'sample,site,collection_date,status,genome_coverage_pct,other,proportion_summary,major_lineages,major_abundances,minor_lineages,minor_abundances,ungrouped_lineages,status_reason\\n' > lineage_summary.csv
    printf 'group\\tinclude\\texclude\\tin_total\\tn_lineages_matched\\twarnings\\tdataset_version\\tgroups_sha256\\n' > lineage_groups_resolved.tsv
    printf 'lineage\\tn_samples\\tmean_abundance\\tmax_abundance\\tsuggested_pattern\\n' > lineage_ungrouped_summary.tsv
    cp groups.csv lineage_groups_used.csv
    ${sites_stub}
    ${xlsx_stub}
    """
}
