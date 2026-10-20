process FREYJA_BOOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(variants), path(depths), path(demix) // only samples whose demix status is ok
    path dataset                                               // the dataset the demix used (directory or barcode file)

    output:
    tuple val(meta), path("${prefix}.boot_summary.tsv"),   emit: summary,    optional: true // data rows: sample, lineage, abundance, lo, hi
    tuple val(meta), path("${prefix}.boot_lineages.csv"),  emit: lineages,   optional: true
    tuple val(meta), path("${prefix}.boot_summarized.csv"), emit: summarized, optional: true
    tuple val(meta), path("${prefix}.boot_lineages_boot.csv"), emit: raw, optional: true // the raw draws; only with `--rawboots` in ext.args
    tuple val(meta), path("${prefix}.boot_status.json"),   emit: status
    tuple val("${task.process}"), val('freyja'), eval('freyja --version | sed "s/.*version //"'), emit: versions_freyja, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Bootstraps are reproducible: --bootseed is fixed and a single thread is used (--nt 1).
    # Like demix, a failure is recorded (status json), never fatal.
    write_status() {
        python3 -c 'import json, sys; print(json.dumps({"sample": sys.argv[1], "status": sys.argv[2], "reason": sys.argv[3]}))' \\
            ${meta.id} "\$1" "\$2" > ${prefix}.boot_status.json
    }

    BARCODES=\$(freyja_dataset_info.py ${dataset} --print-barcodes-path)
    LINEAGEYML=\$(freyja_dataset_info.py ${dataset} --print-lineageyml-path)
    LINEAGEYML_ARG=""
    if [ -n "\$LINEAGEYML" ]; then LINEAGEYML_ARG="--lineageyml \$LINEAGEYML"; fi

    if freyja boot ${variants} ${depths} \\
            --output_base ${prefix}.boot \\
            --barcodes "\$BARCODES" \$LINEAGEYML_ARG \\
            --nt 1 ${args} > ${prefix}.boot.log 2>&1 \\
        && [ -s ${prefix}.boot_lineages.csv ] \\
        && freyja_boot_summary.py --sample ${meta.id} --lineages ${prefix}.boot_lineages.csv \\
            --demix ${demix} --output ${prefix}.boot_summary.tsv 2>> ${prefix}.boot.log; then
        write_status ok ""
    else
        rm -f ${prefix}.boot_summary.tsv ${prefix}.boot_lineages.csv ${prefix}.boot_summarized.csv ${prefix}.boot_lineages_boot.csv ${prefix}.boot_summarized_boot.csv
        write_status failed "freyja boot failed: \$(tail -n 3 ${prefix}.boot.log | tr '\\n' ' ' | cut -c1-300)"
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def raw_stub = (task.ext.args ?: '').contains('--rawboots') ? "printf ',stub\\n0,0.5\\n1,0.5\\n' > ${prefix}.boot_lineages_boot.csv" : ''
    """
    ${raw_stub}
    printf '${meta.id}\\tstub\\t0.500000\\t0.400000\\t0.600000\\n' > ${prefix}.boot_summary.tsv
    printf ',stub\\n0.025,0.4\\n0.975,0.6\\n' > ${prefix}.boot_lineages.csv
    printf ',stub\\n0.025,0.4\\n0.975,0.6\\n' > ${prefix}.boot_summarized.csv
    echo '{"sample": "${meta.id}", "status": "ok", "reason": ""}' > ${prefix}.boot_status.json
    """
}
