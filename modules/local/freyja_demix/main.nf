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
    tuple val(meta), path("${prefix}.demix"),             emit: demix, optional: true // present only when status == ok
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
    """
    # Never swallow a failure: every outcome (ok / failed / skipped_low_coverage) is written to
    # ${prefix}.demix_status.json with a reason and the barcode version, and the task exits 0 so
    # one bad sample does not abort the run. A .demix file exists only when status == ok.
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
    depth_stats.py "\$DEPTHS" --covcut ${covcut} --min-genome-coverage ${min_cov} \\
        --sample ${meta.id} --output ${prefix}.depth_stats.json || gate=\$?

    if [ "\$gate" -eq 10 ]; then
        write_status skipped_low_coverage "less than ${min_cov}% of the genome covered at >= ${covcut}x"
    elif [ "\$gate" -ne 0 ]; then
        rm -f ${prefix}.depth_stats.json
        write_status failed "could not compute coverage from the depths file (depth_stats.py exit \$gate)"
    elif freyja demix "\$VARIANTS" "\$DEPTHS" \\
            --output ${prefix}.demix \\
            --barcodes "\$BARCODES" \$LINEAGEYML_ARG \\
            ${args} > ${prefix}.demix.log 2>&1 \\
        && [ -s ${prefix}.demix ]; then
        write_status ok ""
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
