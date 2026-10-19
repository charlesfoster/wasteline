process SITE_AGGREGATE_CHECK {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path sites, stageAs: 'sites.csv'             // the user's sites CSV (--site_aggregation)
    path samplesheet, stageAs: 'samplesheet.csv' // the samplesheet in effect: unknown site names are reported as warnings

    output:
    path "site_aggregate_check.txt", emit: report // the check's warnings (one per line, may be empty) and its OK line
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    // No errorStrategy override: a format error in the user's file must stop the run within minutes (the default
    // is terminate), long before the mapping starts.
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    site_aggregate.py --check-sites --sites sites.csv --samplesheet samplesheet.csv > site_aggregate_check.txt 2> check.err || {
        cat check.err >&2
        exit 1
    }
    cat check.err >> site_aggregate_check.txt
    """

    stub:
    """
    echo 'site_aggregate.py: sites file OK (stub)' > site_aggregate_check.txt
    """
}
