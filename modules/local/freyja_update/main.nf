process FREYJA_UPDATE {
    tag "${meta.id}:${source}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0'}"

    input:
    // meta: [ id: <target name>, target: <target map> ]; the target's freyja_pathogen / freyja_container_pathogen steer ext.args and publishDir
    // source: 'latest': download with `freyja update` (ext.args: --pathogen; other pathogens write <X>_barcodes.csv, <X>_lineages.yml,
    //         last_barcode_update_other.txt: see bin/freyja_dataset_info.py); 'container': copy the SARS-CoV-2 snapshot bundled in the Freyja
    //         package (the pipeline refuses `container` for other pathogens, validateInputParameters)
    tuple val(meta), val(source)

    output:
    tuple val(meta), path("*", type: 'dir'), emit: dataset
    tuple val("${task.process}"), val('freyja'), eval('freyja --version | sed "s/.*version //"'), emit: versions_freyja, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (source == 'latest') {
        """
        mkdir -p download
        freyja update --outdir download ${args}

        # The dataset directory is named by the ISO form of last_barcode_update.txt (YYYY-MM-DDTHH-MM)
        iso=\$(freyja_dataset_info.py download --print-iso)
        mv download "\$iso"
        printf '{"downloaded_at": "%s", "source": "latest"}\\n' "\$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "\$iso/wasteline_dataset.json"
        """
    }
    else if (source == 'container') {
        """
        data_dir=\$(python3 -c 'import freyja, os; print(os.path.join(os.path.dirname(freyja.__file__), "data"))')
        mkdir container_dataset
        for f in usher_barcodes.feather last_barcode_update.txt lineage_mutations.json curated_lineages.json lineages.yml pathogen_config.yml; do
            cp "\$data_dir/\$f" container_dataset/
        done
        """
    }
    else {
        error("FREYJA_UPDATE: source must be 'latest' or 'container', got '${source}'")
    }

    stub:
    def dir = source == 'latest' ? '2000-01-01T00-00' : 'container_dataset'
    """
    mkdir -p ${dir}
    touch ${dir}/usher_barcodes.feather ${dir}/lineage_mutations.json ${dir}/curated_lineages.json ${dir}/lineages.yml ${dir}/pathogen_config.yml
    echo "01_01_2000-00-00" > ${dir}/last_barcode_update.txt
    echo '{"downloaded_at": "2000-01-01T00:00:00Z", "source": "${source}"}' > ${dir}/wasteline_dataset.json
    """
}
