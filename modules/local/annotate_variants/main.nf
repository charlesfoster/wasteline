process ANNOTATE_VARIANTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.23--h3a4d415_0'
        : 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'}"

    input:
    tuple val(meta), path(vcf), path(tbi) // CONSENSUS_FILTER_VCF `vcf`: the variants that are applied to the consensus
    path fasta                            // reference
    path gff                              // reference annotation (GFF3)

    output:
    tuple val(meta), path("${prefix}.annotated.tsv"), emit: tsv
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Port of the Snakemake rule `amino_acid_consequences`. Input choice: the consensus-filter VCF
    // (genotypes 1/1 or 0/1 already set, only the calls that pass the consensus thresholds), so the table
    // lists exactly the variants that were applied to the consensus. bcftools csq needs proper GT, which
    // the raw ivar VCF lacks. Columns as in the old pipeline: reference, sample_id, gene, nt_pos, ref,
    // alt, effect, aa_bcsq (bcftools notation, e.g. 501N>501Y), aa_standard (N501Y), vaf, depth.
    // A variant without a coding consequence has effect `intergenic` and NA gene/aa; only the first
    // consequence of a variant is reported (as before).
    // Adjacent SNVs in one codon: bcftools gives the joint consequence on the FIRST variant and only
    // `@<position of the first>` on the others. The reference is resolved, so every variant of the codon
    // carries the joint amino-acid change, and the column appended at the end, `compound`, lists the
    // positions that make up the joint change (e.g. 23063+23064; NA for an ordinary variant).
    """
    printf 'reference\\tsample_id\\tgene\\tnt_pos\\tref\\talt\\teffect\\taa_bcsq\\taa_standard\\tvaf\\tdepth\\tcompound\\n' > ${prefix}.annotated.tsv

    bcftools csq -f ${fasta} -g ${gff} --force -pm ${vcf} 2> ${prefix}.csq.log \\
        | bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%ALT_FREQ\\t%INFO/DP\\t%TBCSQ\\n]' \\
        | awk -F'\\t' -v OFS='\\t' -v sample='${meta.id}' '
            {
                csq = \$7
                if (csq ~ /^@/) { csq = by_pos[substr(csq, 2)] } else { by_pos[\$2] = csq }
                effect = "intergenic"; gene = "NA"; aa = "NA"; std = "NA"; compound = "NA"
                n = split(csq, f, "|")
                if (n >= 6) {
                    effect = f[1]; gene = f[2]; aa = f[6]
                    split(aa, p, ">")
                    refp = p[1]; altp = (p[2] != "") ? p[2] : p[1]
                    pos = refp; sub(/[^0-9].*\$/, "", pos)
                    raa = refp; sub(/^[0-9]+/, "", raa)
                    aaa = altp; sub(/^[0-9]+/, "", aaa)
                    std = raa pos aaa
                    if (n >= 7 && f[7] ~ /\\+/) {
                        m = split(f[7], parts, "+"); compound = ""
                        for (i = 1; i <= m; i++) { q = parts[i]; sub(/[^0-9].*\$/, "", q); compound = compound (i > 1 ? "+" : "") q }
                    }
                }
                print \$1, sample, gene, \$2, \$3, \$4, effect, aa, std, \$5, \$6, compound
            }' >> ${prefix}.annotated.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'reference\\tsample_id\\tgene\\tnt_pos\\tref\\talt\\teffect\\taa_bcsq\\taa_standard\\tvaf\\tdepth\\tcompound\\n' > ${prefix}.annotated.tsv
    printf 'stub\\t${meta.id}\\tS\\t23063\\tA\\tT\\tmissense\\t501N>501Y\\tN501Y\\t1\\t100\\tNA\\n' >> ${prefix}.annotated.tsv
    """
}
