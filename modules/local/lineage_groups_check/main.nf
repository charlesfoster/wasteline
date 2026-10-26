process LINEAGE_GROUPS_CHECK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(groups, stageAs: 'groups.csv') // meta = the target (tkey); its lineage groups CSV (--lineage_groups or the target's bundled one)
    path dataset, stageAs: 'dataset'   // Freyja dataset directory holding lineages.yml; [] when the dataset has no hierarchy (format check only)

    output:
    tuple val(meta), path("lineage_groups_check.txt"), emit: report // the check's warnings (one per line, may be empty) and its OK line
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    // No errorStrategy override: a format error in the user's file must stop the run within minutes (the default
    // is terminate), long before the mapping starts. Unknown lineages are warnings only.
    when:
    task.ext.when == null || task.ext.when

    script:
    def dataset_arg = dataset ? '--dataset dataset' : ''
    """
    lineage_groups.py --check-groups --groups groups.csv ${dataset_arg} > lineage_groups_check.txt 2> check.err || {
        cat check.err >&2
        exit 1
    }
    cat check.err >> lineage_groups_check.txt
    """

    stub:
    """
    echo 'lineage_groups.py: groups file OK (stub)' > lineage_groups_check.txt
    """
}
