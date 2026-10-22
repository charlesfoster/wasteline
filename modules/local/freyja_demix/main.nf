process FREYJA_DEMIX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(variants), path(depths)
    path dataset         // Freyja dataset directory or single barcode file (one value shared by all samples)
    path barcode_version // barcode_version.json: recorded in every status file

    output:
    tuple val(meta), path("${prefix}.demix"),             emit: demix, optional: true // present only when status is ok or low_confidence
    tuple val(meta), path("${prefix}.demix_status.json"), emit: status
    tuple val(meta), path("${prefix}.depth_stats.json"),  emit: depth_stats, optional: true
    tuple val("${task.process}"), val('freyja'), eval('freyja --version | sed "s/.*version //"'), emit: versions_freyja, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def covcut = task.ext.covcut
    def min_cov = task.ext.min_genome_coverage
    // Low-confidence band (D18): only when the floor is above 0 and below the high-confidence gate
    def floor = task.ext.floor_genome_coverage
    def band = floor != null && (floor as double) > 0 && (floor as double) < (min_cov as double)
    def floor_arg = band ? "--floor-genome-coverage ${floor}" : ''
    def skip_below = band ? floor : min_cov // coverage under which the sample is skipped (the floor when the band is on)
    // multi-contig reference (PLAN D25): the depths are those of the demixed contig only, and the gate text says so
    def region = task.ext.contig ? "of ${task.ext.contig}" : 'of the genome'
    """
    # Never swallow a failure: every outcome (ok / low_confidence / failed / skipped_low_coverage) is
    # written to ${prefix}.demix_status.json with a reason and the barcode version, and the task
    # exits 0 so one bad sample does not abort the run. A .demix file exists only when status is ok
    # or low_confidence.
    write_status() {
        demix_status.py --sample ${meta.id} --status "\$1" --reason "\$2" \\
            --barcode-info ${barcode_version} --stats ${prefix}.depth_stats.json \\
            "\${@:3}" --output ${prefix}.demix_status.json
    }

    # Inputs may come from earlier runs (demix mode): legacy Snakemake results name the variants
    # <sample>.freyja_variants.tsv(.gz) and may gzip the depths. Freyja labels the sample with the
    # variants file name, so always demix a plain <sample>.variants.tsv: results are then identical
    # whichever layout they came from.
    VARIANTS=${variants}
    DEPTHS=${depths}
    if [ "\$VARIANTS" != "${meta.id}.variants.tsv" ]; then
        case "\$VARIANTS" in
            *.gz) gzip -dc "\$VARIANTS" > ${meta.id}.variants.tsv ;;
            *)    cp -L "\$VARIANTS" ${meta.id}.variants.tsv ;;
        esac
        VARIANTS=${meta.id}.variants.tsv
    fi
    case "\$DEPTHS" in
        *.gz) gzip -dc "\$DEPTHS" > ${meta.id}.plain.depths; DEPTHS=${meta.id}.plain.depths ;;
    esac

    BARCODES=\$(freyja_dataset_info.py ${dataset} --print-barcodes-path)
    LINEAGEYML=\$(freyja_dataset_info.py ${dataset} --print-lineageyml-path)
    LINEAGEYML_ARG=""
    if [ -n "\$LINEAGEYML" ]; then LINEAGEYML_ARG="--lineageyml \$LINEAGEYML"; fi

    gate=0
    depth_stats.py "\$DEPTHS" --covcut ${covcut} --min-genome-coverage ${min_cov} ${floor_arg} \\
        --sample ${meta.id} --output ${prefix}.depth_stats.json || gate=\$?

    if [ "\$gate" -eq 10 ]; then
        write_status skipped_low_coverage "less than ${skip_below}% ${region} covered at >= ${covcut}x"
    elif [ "\$gate" -ne 0 ] && [ "\$gate" -ne 11 ]; then
        rm -f ${prefix}.depth_stats.json
        write_status failed "could not compute coverage from the depths file (depth_stats.py exit \$gate)"
    elif freyja demix "\$VARIANTS" "\$DEPTHS" \\
            --output ${prefix}.demix \\
            --barcodes "\$BARCODES" \$LINEAGEYML_ARG \\
            ${args} > ${prefix}.demix.log 2>&1 \\
        && [ -s ${prefix}.demix ]; then
        if [ "\$gate" -eq 11 ]; then
            COV=\$(sed -n 's/.*"pct_genome_ge_covcut": *\\([0-9.]*\\).*/\\1/p' ${prefix}.depth_stats.json)
            write_status low_confidence "\$(printf '%.1f' "\$COV")% ${region} covered at >= ${covcut}x: below the high-confidence threshold of ${min_cov}%; treat the estimate with care"
        else
            write_status ok ""
        fi
    else
        rm -f ${prefix}.demix
        write_status failed "freyja demix failed" --reason-file ${prefix}.demix.log
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.demix ${prefix}.depth_stats.json
    cat <<'EOF' > ${prefix}.demix_status.json
    {
      "sample": "${meta.id}",
      "status": "ok",
      "reason": "",
      "genome_coverage_pct": null,
      "covcut": null,
      "min_genome_coverage": null,
      "barcode_version": "01_01_2000-00-00",
      "barcode_version_iso": "2000-01-01T00-00"
    }
    EOF
    """
}
