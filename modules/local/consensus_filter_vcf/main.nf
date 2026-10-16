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
    // ALT. Records are first limited to those with DP >= min_depth and ALT_QUAL >= min_alt_qual (bcftools view), then
    // decided per position by the awk program below (docs/wiring/consensus_multiallelic.md):
    //   indel AF >= indel_freq -> GT 1/1; below -> dropped (never IUPAC, never merged)
    //   SNVs at one position (ivar writes one row per ALT; bcftools consensus can apply only one record per position)
    //   are merged into ONE record. Let f_i = ALT_FREQ of each ALT and f_ref = max(0, 1 - sum f_i).
    //   (a) exactly one ALT with f >= consensus_freq                    -> that ALT alone, GT 1/1
    //   (b) iupac on: S = ALTs with f > iupac_min_af and ALT_DP >= iupac_min_alt_reads, plus the reference if
    //       f_ref > iupac_min_af. |S| >= 2 -> one record, ALT = the S alleles (highest frequency first), GT lists every
    //       member (ref = 0): 0/1, 1/2, 0/1/2, 1/2/3, 0/1/2/3. `bcftools consensus -H I` turns that into the IUPAC code.
    //   (c) |S| == 1 and it is an ALT (the reference is not a member, f_ref <= iupac_min_af): a site with one candidate row
    //       -> GT 0/1 (ALT with the reference, the single-ALT rule of the old pipeline); a site with several candidate
    //       rows -> that ALT alone, GT 1/1 (a reference that is not in the mixture is not paired with it)
    //   otherwise (incl. everything when iupac is off, unless (a)) -> no record (reference base).
    // The comparisons with iupac_min_af are strict (">"), as in the single-ALT rule before this change.
    def consensus_freq = task.ext.consensus_freq != null ? task.ext.consensus_freq : 0.75
    def indel_freq = task.ext.indel_freq != null ? task.ext.indel_freq : 0.75
    def min_alt_qual = task.ext.min_alt_qual != null ? task.ext.min_alt_qual : 20
    def min_depth = task.ext.min_depth != null ? task.ext.min_depth : 10
    def iupac = task.ext.iupac != null ? task.ext.iupac : true
    def iupac_min_af = task.ext.iupac_min_af != null ? task.ext.iupac_min_af : 0.25
    def iupac_min_alt_reads = task.ext.iupac_min_alt_reads != null ? task.ext.iupac_min_alt_reads : 5
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Normalise: split multiallelic records and atomise (no reference needed). The per-record thresholds (DP,
    # ALT_QUAL) are tested directly by the -i expression, so no bcftools plugins are needed; the awk program then
    # decides per position (see above) and sets the genotype. Keep it free of gawk-only features (busybox awk).
    bgzip -c ${vcf} > ${prefix}.ivar.vcf.gz
    bcftools index -f ${prefix}.ivar.vcf.gz

    bcftools norm -a -m - -Ou ${prefix}.ivar.vcf.gz | \\
        bcftools view -i 'FORMAT/ALT_QUAL >= ${min_alt_qual} && INFO/DP >= ${min_depth}' -Ov | \\
        awk -F'\\t' -v OFS='\\t' -v cf=${consensus_freq} -v indelf=${indel_freq} -v iupac=${iupac ? 1 : 0} \\
            -v minaf=${iupac_min_af} -v minreads=${iupac_min_alt_reads} '
            # BEGIN merge_snv_alleles
            function field_index(names, want,    n, k, a) {
                n = split(names, a, ":")
                for (k = 1; k <= n; k++) if (a[k] == want) return k
                return 0
            }
            function emit(k, alts, gt,    f, n, s, i, line, sm) {
                # print record k, optionally with another ALT list and genotype
                n = split(gl[k], f, FS)
                if (alts != "") f[5] = alts
                s = split(f[10], sm, ":")
                sm[1] = gt
                f[10] = sm[1]
                for (i = 2; i <= s; i++) f[10] = f[10] ":" sm[i]
                line = f[1]
                for (i = 2; i <= n; i++) line = line OFS f[i]
                print line
            }
            function snv_site(cnt,    i, j, k, t, sumf, fref, nmajor, best, ns, ord, alts, gt, seen, refin, total) {
                sumf = 0; nmajor = 0; best = 0
                for (i = 1; i <= cnt; i++) {
                    sumf += saf[m[i]]
                    if (saf[m[i]] >= cf) { nmajor++; if (best == 0 || saf[m[i]] > saf[best]) best = m[i] }
                }
                fref = 1 - sumf
                if (fref < 0) fref = 0
                fref = int(fref * 1000000 + 0.5) / 1000000
                if (nmajor == 1 || (nmajor > 1 && !iupac)) { emit(best, "", "1/1"); return }
                if (!iupac) return
                ns = 0
                for (i = 1; i <= cnt; i++) {
                    k = m[i]
                    if (saf[k] > minaf && sdp[k] >= minreads && !(salt[k] in seen)) { seen[salt[k]] = 1; ord[++ns] = k }
                }
                for (i = 2; i <= ns; i++) {            # insertion sort, highest frequency first
                    t = ord[i]; j = i - 1
                    while (j >= 1 && saf[ord[j]] < saf[t]) { ord[j + 1] = ord[j]; j-- }
                    ord[j + 1] = t
                }
                refin = (fref > minaf) ? 1 : 0
                total = ns + refin
                if (total >= 2) {
                    alts = ""; gt = refin ? "0" : ""
                    for (i = 1; i <= ns; i++) {
                        alts = (alts == "" ? "" : alts ",") salt[ord[i]]
                        gt = (gt == "" ? "" : gt "/") i
                    }
                    emit(ord[1], alts, gt)
                } else if (ns == 1) {
                    # one member only, and the reference is not a contributor (f_ref <= minaf). A lone candidate row keeps
                    # the single-ALT rule (ALT/reference, 0/1); with competing rows the ALT stands alone (plain base).
                    emit(ord[1], "", (cnt == 1) ? "0/1" : "1/1")
                }
            }
            function flush(    k, i, f, n, cnt, s, sm) {
                for (k = 1; k <= ng; k++) {
                    n = split(gl[k], f, FS)
                    ityp[k] = (length(f[4]) == length(f[5])) ? "S" : "I"
                    iref[k] = f[4]; salt[k] = f[5]
                    s = split(f[10], sm, ":")
                    saf[k] = sm[field_index(f[9], "ALT_FREQ")] + 0
                    sdp[k] = sm[field_index(f[9], "ALT_DP")] + 0
                    used[k] = 0
                }
                for (k = 1; k <= ng; k++) {
                    if (ityp[k] == "I") {
                        if (saf[k] >= indelf) emit(k, "", "1/1")
                    } else if (!used[k]) {
                        cnt = 0
                        for (i = k; i <= ng; i++)
                            if (ityp[i] == "S" && !used[i] && iref[i] == iref[k]) { used[i] = 1; m[++cnt] = i }
                        snv_site(cnt)
                    }
                }
                ng = 0
            }
            /^#/ { print; next }
            {
                if (ng > 0 && (\$1 != gchrom || \$2 != gpos)) flush()
                ng++; gchrom = \$1; gpos = \$2; gl[ng] = \$0
            }
            END { if (ng > 0) flush() }
            # END merge_snv_alleles
            ' | \\
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
