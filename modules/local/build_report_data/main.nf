process BUILD_REPORT_DATA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    // meta = the target the report is about (tkey). Every file except run_info is optional: pass an empty list [] when it does not exist (demix /
    // aggregate modes have no QC, amplicon, consensus or MOI files). `depths` is a list of *.depths.
    tuple val(meta), path(run_info), path(qc_summary), path(amplicon_coverage), path(freyja_aggregated), path(demix_status), path(consensus_stats), path(mutations_of_interest), path(barcode_compare), path(depths)
    // optional extras, [] when absent (samplesheet: one file or a list): aggregate of the ok negative controls, the samplesheet in effect
    // (sample_type, collection_date, site for modes without a QC table), freyja boot intervals and the
    // status of the second-barcode demix
    tuple path(controls_aggregated), path(samplesheet), path(boot_summary), path(compare_status)
    // optional analysis tables, [] when absent: dropout impact (+ per-sample summary), unexplained mutations, time series
    // (summarized + lineages level), Freyja growth rates (+ why there are none, + what went into the fit) and the
    // name of the dashboard file ('' when none; it is linked, not embedded)
    tuple path(dropout_impact), path(dropout_summary), path(unexplained), path(timeseries), path(timeseries_lineages, stageAs: "lineages/*"), path(growth_rates), path(growth_failed), path(growth_inputs), val(dashboard_name)
    // lineage-group summary (+ the audit of what each group matched; for flags and provenance only) and site aggregates (+ meta), [] when absent
    tuple path(lineage_summary), path(lineage_groups_resolved), path(site_aggregates), path(site_aggregates_meta)
    // several reference contigs only (PLAN D25), [] otherwise: per-contig coverage (QC_SAMPLE) and per-contig consensus stats (CONSENSUS_STATS)
    tuple path(contig_coverage), path(consensus_stats_contigs)

    output:
    tuple val(meta), path("report_data.json"), emit: json // meta = the target (tkey: [ id: target name, target ])
    tuple val(meta), path("sample_id_map.csv"), emit: id_map, optional: true // --report_anonymize full only: original -> pseudonym table (PRIVATE: never share it with the report)
    tuple val(meta), path("qc_summary.merged.tsv"), emit: qc_summary, optional: true // the QC table with demix/consensus/barcode columns filled in
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python // not topic: versions: PROVENANCE needs every topic entry before it writes run_info.json, and these run after it

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def qc_arg = qc_summary ? "--qc ${qc_summary}" : ''
    def amplicons_arg = amplicon_coverage ? "--amplicons ${amplicon_coverage}" : ''
    def aggregated_arg = freyja_aggregated ? "--aggregated ${freyja_aggregated}" : ''
    def demix_arg = demix_status ? "--demix-status ${demix_status}" : ''
    def consensus_arg = consensus_stats ? "--consensus ${consensus_stats}" : ''
    def moi_arg = mutations_of_interest ? "--moi ${mutations_of_interest}" : ''
    def compare_arg = barcode_compare ? "--barcode-compare ${barcode_compare}" : ''
    def depths_arg = depths ? "--depths ${depths}" : ''
    def controls_arg = controls_aggregated ? "--controls-aggregated ${controls_aggregated}" : ''
    def sheet_arg = samplesheet ? "--samplesheet ${samplesheet}" : ''
    def boot_arg = boot_summary ? "--boot ${boot_summary}" : ''
    def compare_status_arg = compare_status ? "--compare-status ${compare_status}" : ''
    def dropout_arg = dropout_impact ? "--dropout-impact ${dropout_impact}" : ''
    def dropout_summary_arg = dropout_summary ? "--dropout-summary ${dropout_summary}" : ''
    def unexplained_arg = unexplained ? "--unexplained ${unexplained}" : ''
    def timeseries_arg = timeseries ? "--timeseries ${timeseries}" : ''
    def timeseries_lineages_arg = timeseries_lineages ? "--timeseries-lineages ${timeseries_lineages}" : ''
    def growth_arg = growth_rates ? "--growth-rates ${growth_rates}" : ''
    def growth_failed_arg = growth_failed ? "--growth-failed ${growth_failed}" : ''
    def growth_inputs_arg = growth_inputs ? "--growth-inputs ${growth_inputs}" : ''
    def dashboard_arg = dashboard_name ? "--dashboard-name '${dashboard_name}'" : ''
    def lineage_summary_arg = lineage_summary ? "--lineage-summary ${lineage_summary}" : ''
    def lineage_resolved_arg = lineage_groups_resolved ? "--lineage-groups-resolved ${lineage_groups_resolved}" : ''
    def site_aggregates_arg = site_aggregates ? "--site-aggregates ${site_aggregates}" : ''
    def site_aggregates_meta_arg = site_aggregates_meta ? "--site-aggregates-meta ${site_aggregates_meta}" : ''
    def contig_coverage_arg = contig_coverage ? "--contig-coverage ${contig_coverage}" : ''
    def consensus_contigs_arg = consensus_stats_contigs ? "--consensus-stats-contigs ${consensus_stats_contigs}" : ''
    def merged_arg = qc_summary ? "--write-qc-tsv qc_summary.merged.tsv" : ''
    """
    build_report_data.py \\
        --run-info ${run_info} \\
        ${qc_arg} \\
        ${amplicons_arg} \\
        ${aggregated_arg} \\
        ${demix_arg} \\
        ${consensus_arg} \\
        ${moi_arg} \\
        ${compare_arg} \\
        ${depths_arg} \\
        ${controls_arg} \\
        ${sheet_arg} \\
        ${boot_arg} \\
        ${compare_status_arg} \\
        ${dropout_arg} \\
        ${dropout_summary_arg} \\
        ${unexplained_arg} \\
        ${timeseries_arg} \\
        ${timeseries_lineages_arg} \\
        ${growth_arg} \\
        ${growth_failed_arg} \\
        ${growth_inputs_arg} \\
        ${dashboard_arg} \\
        ${lineage_summary_arg} \\
        ${lineage_resolved_arg} \\
        ${site_aggregates_arg} \\
        ${site_aggregates_meta_arg} \\
        ${contig_coverage_arg} \\
        ${consensus_contigs_arg} \\
        ${merged_arg} \\
        ${args} \\
        --out report_data.json
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo '{"schema_version":1,"mode":"stub","header":{},"provenance":{},"samples":[]}' > report_data.json
    ${qc_summary ? "cp -L ${qc_summary} qc_summary.merged.tsv" : ''}
    ${args.contains('--anonymize-map-out') ? "printf 'type,original,pseudonym\\n' > sample_id_map.csv" : ''}
    """
}
