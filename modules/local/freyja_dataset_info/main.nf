process FREYJA_DATASET_INFO {
    tag "${source}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.5--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.5--pyhdfd78af_0'}"

    input:
    path dataset       // Freyja dataset directory, or a single .csv/.feather barcode file
    val source         // latest | container | pinned | custom
    val dataset_path   // where the dataset came from (recorded verbatim)
    val pathogen       // Freyja --pathogen name

    output:
    path "freyja_data.json",     emit: freyja_data     // full provenance block (PLAN 5.13)
    path "barcode_version.json", emit: barcode_version // small record copied into every status file
    tuple val("${task.process}"), val('freyja'), eval('freyja --version | sed "s/.*version //"'), emit: versions_freyja, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # `freyja demix --version` reports the data bundled in the Freyja package (it ignores
    # --outdir datasets), so it identifies the dataset only when source == container.
    freyja demix --version --pathogen ${pathogen} > reported_version.txt 2>&1 || true

    DATASET_PATH='${dataset_path}'
    if [ '${source}' = latest ]; then
        DATASET_PATH="\$DATASET_PATH/\$(freyja_dataset_info.py ${dataset} --print-iso)"
    fi

    freyja_dataset_info.py ${dataset} \\
        --source ${source} \\
        --pathogen ${pathogen} \\
        --freyja-version "\$(freyja --version | sed 's/.*version //')" \\
        --dataset-path "\$DATASET_PATH" \\
        --reported-file reported_version.txt \\
        --output freyja_data.json \\
        --version-output barcode_version.json
    """

    stub:
    """
    cat <<'EOF' > freyja_data.json
    {
      "source": "${source}",
      "barcode_version": "01_01_2000-00-00",
      "barcode_version_iso": "2000-01-01T00-00",
      "barcodes_file": "usher_barcodes.feather",
      "barcodes_sha256": "stub",
      "barcodes_format": "feather",
      "lineage_mutations_sha256": "stub",
      "curated_lineages_sha256": "stub",
      "lineages_yml_sha256": "stub",
      "pathogen": "${pathogen}",
      "freyja_version": "stub",
      "freyja_reported_barcode_version": null,
      "dataset_path": "${dataset_path}",
      "downloaded_at": null
    }
    EOF
    cat <<'EOF' > barcode_version.json
    {
      "barcode_version": "01_01_2000-00-00",
      "barcode_version_iso": "2000-01-01T00-00",
      "barcodes_sha256": "stub",
      "source": "${source}",
      "freyja_version": "stub"
    }
    EOF
    """
}
