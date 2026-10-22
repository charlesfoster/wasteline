process CONSENSUS_BCFTOOLS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.23--h3a4d415_0'
        : 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(mask) // CONSENSUS_FILTER_VCF `vcf`, LOW_DEPTH_MASK `bed`
    path fasta                                         // reference

    output:
    tuple val(meta), path("${prefix}.consensus.fa"),          emit: fasta
    tuple val(meta), path("${prefix}.consensus_status.json"), emit: status
    tuple val(meta), path("${prefix}.consensus.log"),         emit: log
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // PLAN 5.11: bcftools consensus --mark-del '-' -H I -m <low-depth mask>, FASTA header = sample ID.
    // A sample with no usable coverage (every position masked) or an unrecoverable bcftools failure yields an
    // all-N sequence of reference length plus a status flag. Never an empty file, never a crash.
    // Robustness: when bcftools aborts with "The fasta sequence does not match the REF allele at <contig>:<pos>"
    // exactly the records at that position are dropped from the VCF and bcftools is re-run (max 10 retries; each
    // retry removes the first offending position bcftools reports). The sample then ends with status ok and the
    // dropped variants listed in the status JSON `warnings` array (also in `reason`). Variants bcftools itself
    // skips because they overlap another variant ("overlaps with another variant, skipping") are listed there too:
    // "skipped variant at <contig>:<pos> (inside a deleted region)" when the site lies inside an applied deletion (the
    // deletion wins, expected), otherwise "(overlaps another variant)".
    // Any other failure, or more than 10 offending positions, gives status failed + all-N as before.
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # total length of the reference (all records) and the number of records: one record = the single-contig behaviour (header
    # >sample), several (PLAN D25) = one record per contig, headed >sample|<contig>
    ref_len=\$(awk '/^>/ { next } { total += length(\$0) } END { print total + 0 }' ${fasta})
    n_rec=\$(grep -c '^>' ${fasta})
    mask_bp=\$(awk '{ total += \$3 - \$2 } END { print total + 0 }' ${mask})
    n_var=\$(bcftools view -H ${vcf} | wc -l)
    # span of the called variants: CONSENSUS_FILTER_VCF removed it from the low-depth mask, so it does not count as coverage
    keep_bp=\$(bcftools query -f '%POS0\\t%END\\n' ${vcf} | awk '{ total += \$2 - \$1 } END { print total + 0 }')

    write_n_sequence() {
        if [ "\$n_rec" -gt 1 ]; then
            # all-N, one record per reference contig at its own length
            awk -v id='${meta.id}' '
                function flush(   i, line) {
                    if (name == "") return
                    print ">" id "|" name
                    line = sprintf("%60s", ""); gsub(/ /, "N", line)
                    for (i = 0; i < int(len / 60); i++) print line
                    if (len % 60 > 0) print substr(line, 1, len % 60)
                }
                /^>/ { flush(); name = substr(\$1, 2); len = 0; next }
                { len += length(\$0) }
                END { flush() }' ${fasta} > ${prefix}.consensus.fa
            return
        fi
        awk -v id='${meta.id}' -v len="\$ref_len" 'BEGIN {
            print ">" id
            line = sprintf("%60s", ""); gsub(/ /, "N", line)
            for (i = 0; i < int(len / 60); i++) print line
            if (len % 60 > 0) print substr(line, 1, len % 60)
        }' > ${prefix}.consensus.fa
    }
    json_clean() { printf '%s' "\$1" | tr -d '"\\\\' | tr '\\n\\t' '  '; }
    # write_status <status> <reason>; the warnings (one per line in warnings.txt, max 20) become the JSON array
    write_status() {
        reason=\$(json_clean "\$2")
        warn_json=\$(awk -v max=20 'NF { gsub(/["\\\\]/, ""); if (++n <= max) a[n] = \$0 }
            END { out = ""; for (i = 1; i <= n && i <= max; i++) out = out (i > 1 ? ", " : "") "\\"" a[i] "\\"";
                  if (n > max) out = out ", \\"(" n - max " more)\\""; print "[" out "]" }' warnings.txt)
        printf '{"sample": "%s", "status": "%s", "reason": "%s", "warnings": %s}\\n' '${meta.id}' "\$1" "\$reason" "\$warn_json" > ${prefix}.consensus_status.json
    }

    : > warnings.txt
    : > ${prefix}.consensus.log

    if [ "\$ref_len" -le 0 ]; then
        echo "reference has no sequence" >&2
        exit 1
    elif [ "\$((ref_len - mask_bp))" -le "\$keep_bp" ]; then
        echo "all \$ref_len positions are below the minimum depth (\$keep_bp bp are called variants)" >> ${prefix}.consensus.log
        write_n_sequence
        write_status empty "no position reaches the minimum consensus depth (no usable coverage)"
    else
        cur_vcf=${vcf}
        attempt=0
        consensus_ok=0
        while true; do
            if bcftools consensus \\
                    -f ${fasta} \\
                    --mark-del '-' \\
                    -H I \\
                    -m ${mask} \\
                    ${args} \\
                    "\$cur_vcf" >| ${prefix}.consensus.raw.fa 2>| attempt.log \\
                && [ "\$(grep -c '^>' ${prefix}.consensus.raw.fa)" -eq "\$n_rec" ] \\
                && [ "\$(grep -v '^>' ${prefix}.consensus.raw.fa | tr -d '\\n' | wc -c)" -gt 0 ]; then
                consensus_ok=1
            fi
            cat attempt.log >> ${prefix}.consensus.log
            [ "\$consensus_ok" -eq 1 ] && break

            # "The fasta sequence does not match the REF allele at <contig>:<pos>:" -> drop that position, retry
            bad=\$(awk '/does not match the REF allele at / { s = \$0; sub(/^.*allele at /, "", s); sub(/:\$/, "", s); print s; exit }' attempt.log)
            [ -n "\$bad" ] && [ "\$attempt" -lt 10 ] || break
            chrom=\${bad%:*}
            pos=\${bad##*:}
            case "\$pos" in ''|*[!0-9]*) break ;; esac
            attempt=\$((attempt + 1))
            n_before=\$(bcftools view -H "\$cur_vcf" | wc -l)
            bcftools view -e "\$(printf 'CHROM="%s" && POS=%s' "\$chrom" "\$pos")" -Oz -o retry_\$attempt.vcf.gz "\$cur_vcf" 2>> ${prefix}.consensus.log
            bcftools index -t -f retry_\$attempt.vcf.gz
            [ "\$(bcftools view -H retry_\$attempt.vcf.gz | wc -l)" -lt "\$n_before" ] || break
            cur_vcf=retry_\$attempt.vcf.gz
            echo "retry \$attempt: dropped the variant(s) at \$chrom:\$pos (REF does not match the reference)" >> ${prefix}.consensus.log
            echo "skipped variant at \$chrom:\$pos (REF mismatch)" >> warnings.txt
        done
        # variants bcftools itself skipped (overlapping records). A site inside the span of a deletion that was applied
        # is an expected note (the deletion wins), any other overlap keeps the generic wording.
        {
            bcftools query -f '%CHROM\\t%POS\\t%END\\t%REF\\t%ALT\\n' "\$cur_vcf" | awk -F'\\t' 'length(\$4) > length(\$5)'
            echo "END_OF_DELETIONS"
            grep -o 'The site [^ ]* overlaps with another variant' ${prefix}.consensus.log || true
        } | awk -F'\\t' '
            /^END_OF_DELETIONS\$/ { listing = 1; next }
            !listing { nd++; dc[nd] = \$1; ds[nd] = \$2; de[nd] = \$3; next }
            {
                site = \$0; sub(/^The site /, "", site); sub(/ overlaps with another variant\$/, "", site)
                if (!(site in skipped)) { skipped[site] = 1; order[++ns] = site }
            }
            END {
                for (k = 1; k <= ns; k++) {
                    site = order[k]
                    p = site; sub(/^.*:/, "", p); c = site; sub(/:[^:]*\$/, "", c)
                    inside = 0
                    # inside an applied deletion (one that was not skipped itself)
                    for (i = 1; i <= nd; i++)
                        if (dc[i] == c && p + 0 > ds[i] + 0 && p + 0 <= de[i] + 0 && !((c ":" ds[i]) in skipped)) inside = 1
                    print "skipped variant at " site (inside ? " (inside a deleted region)" : " (overlaps another variant)")
                }
            }' >> warnings.txt

        if [ "\$consensus_ok" -eq 1 ]; then
            # bcftools keeps the full FASTA header; one record -> >sample, several -> >sample|<first word of the header>
            awk -v id='${meta.id}' -v multi="\$n_rec" '/^>/ { print (multi > 1 ? ">" id "|" substr(\$1, 2) : ">" id); next } { print }' ${prefix}.consensus.raw.fa > ${prefix}.consensus.fa
            n_warn=\$(grep -c . warnings.txt || true)
            n_var_final=\$(bcftools view -H "\$cur_vcf" | wc -l)
            if [ "\$n_warn" -gt 0 ]; then
                write_status ok "\$n_warn variant(s) not applied: \$(paste -sd ';' warnings.txt | sed 's/;/; /g')"
            elif [ "\$n_var_final" -eq 0 ]; then
                write_status ok "no variant passed the consensus filters; sequence equals the reference outside masked regions"
            else
                write_status ok ""
            fi
        else
            write_n_sequence
            rm -f warnings.txt; touch warnings.txt   # the genome is all-N: the skipped variants are only in the log
            write_status failed "bcftools consensus failed: \$(grep -v '^\\s*\$' attempt.log | head -n 1)"
        fi
    fi
    rm -f ${prefix}.consensus.raw.fa attempt.log retry_*.vcf.gz retry_*.vcf.gz.tbi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>${meta.id}\\nNNNNNNNNNN\\n' > ${prefix}.consensus.fa
    printf '{"sample": "${meta.id}", "status": "ok", "reason": "", "warnings": []}\\n' > ${prefix}.consensus_status.json
    : > ${prefix}.consensus.log
    """
}
