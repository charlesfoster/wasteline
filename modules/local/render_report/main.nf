process RENDER_REPORT {
    tag "report"
    label 'process_single'

    // CONTAINER: no small biocontainer ships python + jinja2, so this uses the MultiQC biocontainer purely
    // as a pinned python 3.14 + jinja2 3.1.6 image (jinja2 is a hard MultiQC dependency, so it cannot
    // vanish on a rebuild, and it does not move with Freyja auto-bumps). Interim: the environment.yml
    // (python + jinja2 only) is the source of truth; replace this with a Seqera Wave image built from it
    // (see docs/wiring/report_modules.md) to avoid pulling ~440 MB for one library.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_2'
        : 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_2'}"

    input:
    path report_data // report_data.json from BUILD_REPORT_DATA (or an existing results dir, --mode report)
    path template    // ${projectDir}/assets/templates/report.html.j2, staged so it is checksummed

    output:
    path "report.html", emit: html
    tuple val("${task.process}"), val('jinja2'), eval('python -c "import jinja2; print(jinja2.__version__)"'), emit: versions_jinja2 // not topic: versions: PROVENANCE needs every topic entry before it writes run_info.json, and these run after it

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    render_report.py \\
        --data ${report_data} \\
        --template ${template} \\
        ${args} \\
        --out report.html
    """

    stub:
    """
    echo '<!DOCTYPE html><html><head><meta charset="utf-8"><title>wasteline report (stub)</title></head><body></body></html>' > report.html
    """
}
