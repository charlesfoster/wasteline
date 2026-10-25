process PANEL_CALL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    // stats: PANEL_ASSIGN's <id>.panel_stats.tsv. meta needs id; sample_type (default 'sample') is optional.
    // Thresholds (uncalibrated starting values: 100, 5, 5): a target is detected with at least min_fragments marker
    // fragments AND marker breadth >= min_breadth_pct AND (the top target of its family OR a share of the family's fragments
    // >= family_share_pct percent).
    tuple val(meta), path(stats)
    path contig_map
    val min_fragments
    val min_breadth_pct
    val family_share_pct

    output:
    // one row per panel target, always (also for a sample with no panel reads)
    tuple val(meta), path("${prefix}.targets.tsv"), emit: targets
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def sample_type = meta.sample_type ?: 'sample'
    """
    panel_select.py call \\
        --sample '${meta.id}' \\
        --sample-type '${sample_type}' \\
        --stats ${stats} \\
        --contig-map ${contig_map} \\
        --min-fragments ${min_fragments} \\
        --min-breadth-pct ${min_breadth_pct} \\
        --family-share-pct ${family_share_pct} \\
        --output ${prefix}.targets.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tsample_type\\ttarget\\tfamily\\tstatus\\treason\\tsample_status\\tmixed\\tsecondary_targets\\tfragments\\tambiguous_fragments\\tfraction_of_family\\tbreadth_pct\\tmean_depth\\tmean_mapq\\tsecond_target\\tmargin\\trescue\\tanalysed\\n' > ${prefix}.targets.tsv
    printf '${meta.id}\\t${meta.sample_type ?: 'sample'}\\tT1\\tF\\tdetected\\t\\tok\\tfalse\\t\\t500\\t0\\t1.0000\\t100.00\\t50.00\\t60.0\\t\\t\\t\\ttrue\\n' >> ${prefix}.targets.tsv
    """
}
