process BOOT_GROUP_INTERVALS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path aggregated, stageAs: 'freyja_aggregated.tsv' // summary/freyja_aggregated.tsv (the point estimates)
    path boots                                        // every sample's raw draws: <id>.boot_lineages_boot.csv (freyja boot --rawboots)
    path groups, stageAs: 'groups.csv'                // the user's lineage groups CSV (--lineage_groups)
    path dataset, stageAs: 'dataset'                  // Freyja dataset directory holding lineages.yml

    output:
    path "lineage_group_intervals.tsv", emit: intervals // sample, group, abundance, lo, hi, n_boot (one row per sample x group + other)
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: --decimals <n>
    def args = task.ext.args ?: ''
    """
    boot_group_intervals.py \\
        --aggregated freyja_aggregated.tsv \\
        --groups groups.csv \\
        --dataset dataset \\
        ${args} \\
        --output lineage_group_intervals.tsv \\
        --boots ${boots}
    """

    stub:
    """
    printf 'sample\\tgroup\\tabundance\\tlo\\thi\\tn_boot\\n' > lineage_group_intervals.tsv
    """
}
