process RENDER_REPORT {
    tag "${meta.id}"
    label 'process_single'

    // CONTAINER: no biocontainer ships python + jinja2 (jinja2 is conda-forge only), so these images were built
    // from environment.yml (python 3.12.12 + jinja2 3.1.6) with the Seqera Wave CLI (v1.8.2) in freeze mode into
    // Seqera's free community registry; the tag is the content hash of the recipe, so it never moves:
    //   wave --conda-channels conda-forge --conda-package python=3.12.12 --conda-package jinja2=3.1.6 \
    //        --platform linux/amd64 --freeze --await                 (add --singularity for the oras:// image)
    // Rebuild and replace both strings whenever environment.yml changes (see docs/development.md).
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_jinja2:c7f3fc5ddefd422e'
        : 'community.wave.seqera.io/library/python_jinja2:b2615084ab562e7a'}"

    input:
    tuple val(meta), path(report_data) // meta = the target (tkey); report_data.json from BUILD_REPORT_DATA (or an existing results dir, --mode report)
    path template    // ${projectDir}/assets/templates/report.html.j2, staged so it is checksummed

    output:
    tuple val(meta), path("report.html"), emit: html
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
