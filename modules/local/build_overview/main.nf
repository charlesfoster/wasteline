process BUILD_OVERVIEW {
    tag "overview"
    label 'process_single'

    // CONTAINER: the same image as RENDER_REPORT (python 3.12.12 + jinja2 3.1.6, built from environment.yml with the Seqera
    // Wave CLI; see docs/development.md). Change both modules together.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_jinja2:c7f3fc5ddefd422e'
        : 'community.wave.seqera.io/library/python_jinja2:b2615084ab562e7a'}"

    input:
    path targets_tsv                          // summary/targets.tsv of the panel run
    path manifest                             // tab separated: target, role, path (as staged below); roles aggregated|qc|demix_status|run_info|report
    path inputs, stageAs: 'in/?/*'            // the per-target files the manifest lists (any may be missing for a target)
    path template                             // ${projectDir}/assets/templates/overview.html.j2, staged so it is checksummed

    output:
    path "overview.html"                      , emit: html
    path "overview_sample_id_map.csv"         , emit: id_map, optional: true // --report_anonymize full only: original -> pseudonym (PRIVATE: never share it with the overview)
    tuple val("${task.process}"), val('jinja2'), eval('python -c "import jinja2; print(jinja2.__version__)"'), emit: versions_jinja2 // not topic: versions: it runs after PROVENANCE, which needs every topic entry before it writes run_info.json

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    build_overview.py \\
        --targets ${targets_tsv} \\
        --inputs ${manifest} \\
        --template ${template} \\
        ${args} \\
        --out overview.html
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo '<!DOCTYPE html><html><head><meta charset="utf-8"><title>wasteline panel overview (stub)</title></head><body></body></html>' > overview.html
    ${args.contains('--anonymize-map-out') ? "printf 'type,original,pseudonym\\n' > overview_sample_id_map.csv" : ''}
    """
}
