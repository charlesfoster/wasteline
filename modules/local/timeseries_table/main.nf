process TIMESERIES_TABLE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path aggregated                           // freyja_aggregated.tsv
    path metadata, stageAs: 'metadata.csv'    // samplesheet-like CSV: sample, sample_type, collection_date, site and an optional numeric viral_load

    output:
    path "freyja_timeseries.tsv", emit: table // long table: site, collection_date, lineage_group, abundance, abundance_load_weighted, n_samples, mean_viral_load, total_viral_load (header always written)
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: --level summarized|lineages, --min-abundance <fraction>, --min-coverage <percent>
    def args = task.ext.args ?: ''
    """
    timeseries.py \\
        --aggregated ${aggregated} \\
        --metadata metadata.csv \\
        ${args} \\
        --output freyja_timeseries.tsv
    """

    stub:
    """
    printf 'site\\tcollection_date\\tlineage_group\\tabundance\\tabundance_load_weighted\\tn_samples\\tmean_viral_load\\ttotal_viral_load\\n' > freyja_timeseries.tsv
    """
}
