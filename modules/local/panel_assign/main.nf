process PANEL_ASSIGN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    // bam: the sample's reads mapped to panel.fasta (any order: it is collated here). contig_map: PANEL_PREPARE's table.
    // min_mapq: a fragment is assigned to a target only if every mapped primary alignment of it is on that target's contigs
    // with MAPQ >= min_mapq (the panel_mapq parameter, default 10); otherwise it is ambiguous.
    tuple val(meta), path(bam)
    path contig_map
    val min_mapq

    output:
    // columns sample, target, fragments, ambiguous_fragments, aligned_bases, mean_mapq, breadth_pct, mean_depth, panel_fragments;
    // one row per target of the panel (zeros when nothing mapped); the numbers describe the target's MARKER contigs only.
    tuple val(meta), path("${prefix}.panel_stats.tsv"), emit: stats
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Not interpolated (''' string): $1 etc. are awk's. No single quote inside, the shell quotes it with '.
    def assign_awk = '''
        BEGIN { FS = OFS = "\\t" }
        function aligned(c,    total, num, op) {
            total = 0
            while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                num = substr(c, 1, RLENGTH - 1) + 0
                op = substr(c, RLENGTH, 1)
                if (op == "M" || op == "=" || op == "X") total += num
                c = substr(c, RLENGTH + 1)
            }
            return total
        }
        # one fragment = all the mapped primary records of one read name
        function flush(    i, t, first, ambiguous) {
            if (n == 0) return
            panel_fragments++
            first = ""
            ambiguous = 0
            for (i = 1; i <= n; i++) {
                t = (rname[i] in target) ? target[rname[i]] : "?"
                if (first == "") first = t
                else if (t != first) ambiguous = 1
                if (t == "?" || mapq_of[i] < min_mapq) ambiguous = 1
            }
            for (i = 1; i <= n; i++) {
                if (!(rname[i] in target) || !marker[rname[i]]) continue
                t = target[rname[i]]
                if (ambiguous) {
                    if (!(t in counted)) ambig[t]++
                } else {
                    if (!(t in counted)) frag[t]++
                    bases[t] += aligned(cigar[i])
                    qsum[t] += mapq_of[i]
                    qn[t]++
                    print record[i]
                }
                counted[t] = 1
            }
            delete counted
            n = 0
        }
        FNR == NR {
            if (FNR > 1) {
                target[$1] = $2
                marker[$1] = ($4 == "true")
                if (!($2 in known)) { known[$2] = 1; order[++ntargets] = $2 }
            }
            next
        }
        $1 ~ /^@(HD|SQ|RG|PG|CO)$/ { print; next }
        {
            if ($1 != current) { flush(); current = $1 }
            if (int($2 / 4) % 2) next
            n++
            record[n] = $0
            rname[n] = $3
            mapq_of[n] = $5 + 0
            cigar[n] = $6
        }
        END {
            flush()
            for (i = 1; i <= ntargets; i++) {
                t = order[i]
                printf "%s\\t%d\\t%d\\t%d\\t%.1f\\t%d\\n", t, frag[t], ambig[t], bases[t], (qn[t] ? qsum[t] / qn[t] : 0), panel_fragments > partial
            }
        }
    '''
    // the breadth and mean depth of each target's marker contigs (length-weighted over several), joined onto the counts
    def stats_awk = '''
        BEGIN { FS = OFS = "\\t" }
        FNR == 1 { file++ }
        file == 1 {
            if (FNR > 1) { target[$1] = $2; marker[$1] = ($4 == "true") }
            next
        }
        file == 2 {
            if (substr($1, 1, 1) == "#" || !($1 in target) || !marker[$1]) next
            t = target[$1]
            len = $3 - $2 + 1
            total[t] += len
            covered[t] += $5
            depth[t] += $7 * len
            next
        }
        {
            t = $1
            rows[++m] = sample OFS t OFS $2 OFS $3 OFS $4 OFS $5 OFS sprintf("%.2f", total[t] ? 100 * covered[t] / total[t] : 0) OFS sprintf("%.2f", total[t] ? depth[t] / total[t] : 0) OFS $6
        }
        END {
            print "sample", "target", "fragments", "ambiguous_fragments", "aligned_bases", "mean_mapq", "breadth_pct", "mean_depth", "panel_fragments"
            for (i = 1; i <= m; i++) print rows[i]
        }
    '''
    """
    # Fragments are assigned from the primary alignments only; the assigned marker-contig records are kept for samtools coverage.
    samtools collate -O -u ${bam} collate_tmp \\
        | samtools view -h -F 0x900 - \\
        | awk -v min_mapq=${min_mapq} -v partial=partial_counts.tsv '${assign_awk}' ${contig_map} - \\
        | samtools sort -o assigned_marker.bam -
    samtools coverage assigned_marker.bam > coverage.tsv
    awk -v sample='${meta.id}' '${stats_awk}' ${contig_map} coverage.tsv partial_counts.tsv > ${prefix}.panel_stats.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\ttarget\\tfragments\\tambiguous_fragments\\taligned_bases\\tmean_mapq\\tbreadth_pct\\tmean_depth\\tpanel_fragments\\n' > ${prefix}.panel_stats.tsv
    printf '${meta.id}\\tT1\\t500\\t0\\t50000\\t60.0\\t100.00\\t50.00\\t500\\n${meta.id}\\tT2\\t0\\t0\\t0\\t0.0\\t0.00\\t0.00\\t500\\n' >> ${prefix}.panel_stats.tsv
    """
}
