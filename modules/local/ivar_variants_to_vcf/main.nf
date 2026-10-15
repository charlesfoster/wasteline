process IVAR_VARIANTS_TO_VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(variants) // Freyja/ivar <sample>.variants.tsv (header-only is fine)

    output:
    tuple val(meta), path("${prefix}.ivar.vcf"), emit: vcf
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: bin/ivar_variants_to_vcf.py options. Do NOT pass --pass_only / --allele_freq_thresh
    // for the consensus: AF, quality and depth thresholds are applied by CONSENSUS_FILTER_VCF.
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar_variants_to_vcf.py ${args} ${variants} ${prefix}.ivar.raw.vcf

    # The converter names the VCF sample column after the input file ("S.variants"); use the sample ID.
    awk -F'\\t' -v OFS='\\t' -v id='${meta.id}' '/^#CHROM/ { \$NF = id } { print }' ${prefix}.ivar.raw.vcf > ${prefix}.ivar.vcf
    rm ${prefix}.ivar.raw.vcf
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' > ${prefix}.ivar.vcf
    """
}
