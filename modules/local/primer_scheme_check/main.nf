process PRIMER_SCHEME_CHECK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(read_starts)   // READ_STARTS output: strand, five_prime_pos
    val declared                         // name of a bundled scheme (a BED in beds_dir); ignored when custom_bed is given
    path beds_dir                        // directory of <scheme>.bed files to score (the bundled assets/primer_schemes)
    path custom_bed                      // user primer BED ([] when absent); then it is the declared scheme (name = file stem)

    output:
    tuple val(meta), path("${prefix}.primer_scheme_check.tsv"),   emit: tsv     // one row per scheme
    tuple val(meta), path("${prefix}.primer_scheme_summary.tsv"), emit: summary // one row: scheme_check ok|mismatch|uninformative|NA
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: --window 2 --min-enrichment 5 --min-reads 500 --genome-length 29903 --max-p-value 0.001
    //   --min-hit-fraction 0.05 --equivalence-jaccard 0.9
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def declared_arg = custom_bed ? custom_bed : declared
    """
    primer_scheme_check.py \\
        --starts ${read_starts} \\
        --declared ${declared_arg} \\
        --sample ${meta.id} \\
        --beds-dir ${beds_dir} \\
        ${args} \\
        --output ${prefix}.primer_scheme_check.tsv \\
        --summary-out ${prefix}.primer_scheme_summary.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def declared_name = custom_bed ? custom_bed.baseName : declared
    """
    printf 'scheme\\tn_primer_5prime_sites\\tn_reads\\thits\\texpected_hits\\thit_fraction\\tenrichment\\tp_value\\tequivalent_to_declared\\n' > ${prefix}.primer_scheme_check.tsv
    printf 'sample\\tdeclared_scheme\\tdeclared_enrichment\\tbest_scheme\\tbest_enrichment\\tscheme_check\\treason\\tdeclared_hit_fraction\\tbest_hit_fraction\\n' > ${prefix}.primer_scheme_summary.tsv
    printf '${meta.id}\\t${declared_name}\\tNA\\tNA\\tNA\\tNA\\tstub\\tNA\\tNA\\n' >> ${prefix}.primer_scheme_summary.tsv
    """
}
