process QC_SAMPLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    // flagstat, fastp_json and amplicons are optional: pass [] to omit one. Missing optional inputs
    // give NA columns in the QC row (never a failure). meta needs id; sample_type (default
    // 'sample'), collection_date and site are optional.
    tuple val(meta), path(depths), path(flagstat), path(fastp_json), path(amplicons)

    output:
    tuple val(meta), path("${prefix}.qc_row.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: thresholds, e.g. --min-depth 10 --consensus-min-depth 10 --pass-min-pct 80
    //   --fail-below-pct 50 --max-dropped 10 --ntc-max-mapped-reads 1000 --ntc-max-genome-coverage-pct 10
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sample_type = meta.sample_type ?: 'sample'
    // single-quote free text for the shell
    def date_q = meta.collection_date ? meta.collection_date.toString().replace("'", "'\\''") : ''
    def site_q = meta.site ? meta.site.toString().replace("'", "'\\''") : ''
    def date_arg = meta.collection_date ? "--collection-date '${date_q}'" : ''
    def site_arg = meta.site ? "--site '${site_q}'" : ''
    def flagstat_arg = flagstat ? "--flagstat ${flagstat}" : ''
    def fastp_arg = fastp_json ? "--fastp-json ${fastp_json}" : ''
    def amplicons_arg = amplicons ? "--amplicons ${amplicons}" : ''
    """
    qc_summary.py \\
        --sample ${meta.id} \\
        --sample-type ${sample_type} \\
        ${date_arg} \\
        ${site_arg} \\
        --depths ${depths} \\
        ${flagstat_arg} \\
        ${fastp_arg} \\
        ${amplicons_arg} \\
        ${args} \\
        --output ${prefix}.qc_row.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tsample_type\\tcollection_date\\tsite\\traw_reads\\ttrimmed_reads\\tmapped_reads\\tmapped_pct\\tmean_depth\\tdepth_sd\\tmedian_depth\\tpct_genome_ge_1x\\tpct_genome_ge_10x\\tpct_genome_ge_min_depth\\tamplicons_total\\tamplicons_dropped\\tqc_status\\tqc_reasons\\tdemix_status\\tdemix_reason\\tconsensus_length\\tconsensus_n_pct\\tbarcode_version\\n' > ${prefix}.qc_row.tsv
    printf '${meta.id}\\t${meta.sample_type ?: 'sample'}\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\t\\tNA\\tNA\\tNA\\tNA\\tNA\\n' >> ${prefix}.qc_row.tsv
    """
}
