process PROVENANCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(base_info) // meta = the target (tkey); JSON assembled by Nextflow: versions, mode, params, containers, start time
    path samplesheet
    path primer_bed, stageAs: 'primer_beds/?/*' // every distinct primer BED of the samples; [] when none
    path fasta           // the target's reference; [] in modes that do not map reads
    path freyja_data     // freyja_data.json from FREYJA_DATASET_INFO (mandatory)
    path versions_yml    // collated software versions
    path freyja_data_compare, stageAs: 'freyja_data_compare.json' // freyja_data.json of --compare_barcodes, [] when not used

    output:
    tuple val(meta), path("run_info.json"), emit: run_info

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    make_run_info.py \\
        --base ${base_info} \\
        --samplesheet ${samplesheet} \\
        ${primer_bed ? "--bed ${primer_bed}" : ''} \\
        ${fasta ? "--fasta ${fasta}" : ''} \\
        --dataset-info ${freyja_data} \\
        --versions ${versions_yml} \\
        ${freyja_data_compare ? "--dataset-info-compare ${freyja_data_compare}" : ''} \\
        --output run_info.json
    """

    stub:
    """
    make_run_info.py \\
        --base ${base_info} \\
        --samplesheet ${samplesheet} \\
        ${primer_bed ? "--bed ${primer_bed}" : ''} \\
        ${fasta ? "--fasta ${fasta}" : ''} \\
        --dataset-info ${freyja_data} \\
        --versions ${versions_yml} \\
        ${freyja_data_compare ? "--dataset-info-compare ${freyja_data_compare}" : ''} \\
        --output run_info.json
    """
}
