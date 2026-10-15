process CONSENSUS_FILTER_VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.23--h3a4d415_0'
        : 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'}"

    input:
    tuple val(meta), path(vcf) // IVAR_VARIANTS_TO_VCF output (plain VCF; header-only is fine)

    output:
    tuple val(meta), path("${prefix}.consensus_variants.vcf.gz"), path("${prefix}.consensus_variants.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}.called_variants.bed"), emit: bed // positions of the variants that will be applied
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Thresholds (PLAN 5.11, ported from the Snakemake ivar_bcftools_setGT rule). Set in conf/modules.config from
    // params: ext.consensus_freq, ext.indel_freq, ext.min_alt_qual, ext.min_depth (DP of a call, NOT the mask depth),
    // ext.iupac, ext.iupac_min_af, ext.iupac_min_alt_reads. Freyja runs with --varthresh 0.0, so the TSV holds every
    // ALT. For any call DP >= min_depth and ALT_QUAL >= min_alt_qual, then:
    //   SNV  AF >= consensus_freq                                   -> GT 1/1 (alt allele)
    //   SNV  iupac_min_af < AF < consensus_freq, ALT_DP >= reads    -> GT 0/1 (IUPAC code via bcftools consensus -H I)
    //   SNV  otherwise (incl. everything when iupac is off)         -> dropped (reference base)
    //   indel AF >= indel_freq -> GT 1/1; below -> dropped (never IUPAC)
    def consensus_freq = task.ext.consensus_freq != null ? task.ext.consensus_freq : 0.75
    def indel_freq = task.ext.indel_freq != null ? task.ext.indel_freq : 0.75
    def min_alt_qual = task.ext.min_alt_qual != null ? task.ext.min_alt_qual : 20
    def min_depth = task.ext.min_depth != null ? task.ext.min_depth : 10
    def iupac = task.ext.iupac != null ? task.ext.iupac : true
    def iupac_min_af = task.ext.iupac_min_af != null ? task.ext.iupac_min_af : 0.25
    def iupac_min_alt_reads = task.ext.iupac_min_alt_reads != null ? task.ext.iupac_min_alt_reads : 5
    def het_clause = iupac ? " || (TYPE!=\"indel\" && FORMAT/ALT_FREQ > ${iupac_min_af} && FORMAT/ALT_DP >= ${iupac_min_alt_reads})" : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Normalise: split multiallelic records and atomise (no reference needed). TYPE and the per-sample FORMAT fields
    # are tested directly by the -i expression, so no bcftools plugins are needed. The expression selects every
    # record that is applied or IUPAC-coded; awk then sets the genotype (1/1 alt, 0/1 ambiguous).
    bgzip -c ${vcf} > ${prefix}.ivar.vcf.gz
    bcftools index -f ${prefix}.ivar.vcf.gz

    bcftools norm -a -m - -Ou ${prefix}.ivar.vcf.gz | \\
        bcftools view \\
            -i 'FORMAT/ALT_QUAL >= ${min_alt_qual} && INFO/DP >= ${min_depth} && ((TYPE="indel" && FORMAT/ALT_FREQ >= ${indel_freq}) || (TYPE!="indel" && FORMAT/ALT_FREQ >= ${consensus_freq})${het_clause})' \\
            -Ov | \\
        awk -F'\\t' -v OFS='\\t' -v cf=${consensus_freq} '
            /^#/ { print; next }
            {
                n = split(\$10, f, ":")
                af = f[n] + 0                      # ALT_FREQ is the last FORMAT field written by the converter
                if (length(\$4) == length(\$5) && af < cf) f[1] = "0/1"; else f[1] = "1/1"
                \$10 = f[1]; for (i = 2; i <= n; i++) \$10 = \$10 ":" f[i]
                print
            }' | \\
        bgzip -c > ${prefix}.consensus_variants.vcf.gz
    bcftools index -t -f ${prefix}.consensus_variants.vcf.gz

    # Positions of the called variants, never low-depth masked (replaces the old bedtools subtract)
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' ${prefix}.consensus_variants.vcf.gz > ${prefix}.called_variants.bed
    rm ${prefix}.ivar.vcf.gz ${prefix}.ivar.vcf.gz.csi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consensus_variants.vcf.gz ${prefix}.consensus_variants.vcf.gz.tbi ${prefix}.called_variants.bed
    """
}
