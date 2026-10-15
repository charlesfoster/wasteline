process DISCOVER_RESULTS {
    tag "${results_dir.name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(index), path(results_dir, stageAs: 'results_dir')
    val sheet_mode // demix | aggregate

    output:
    path "discovered_samplesheet_${index}.csv", emit: samplesheet
    path "skipped_${index}.tsv",                emit: skipped // sample<TAB>reason of samples found but unusable
    path "discovery_${index}.log",              emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Single source of truth for reading wasteline and legacy Snakemake results: bin/make_samplesheet.py.
    # The directory is only read (never written to). The real path is resolved first so the sheet
    # does not point through the symlink Nextflow stages for this task.
    make_samplesheet.py \\
        --from-results "\$(readlink -f results_dir)" \\
        --mode ${sheet_mode} \\
        --skipped-out skipped_${index}.tsv \\
        -o discovered_samplesheet_${index}.csv 2> discovery_${index}.log || {
            cat discovery_${index}.log >&2
            exit 1
        }
    cat discovery_${index}.log >&2
    """

    stub:
    // Discovery is cheap and the rest of the pipeline needs the real rows, so the stub is the real command
    """
    make_samplesheet.py \\
        --from-results "\$(readlink -f results_dir)" \\
        --mode ${sheet_mode} \\
        --skipped-out skipped_${index}.tsv \\
        -o discovered_samplesheet_${index}.csv 2> discovery_${index}.log
    """
}
