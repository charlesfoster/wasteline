process GUNZIP_FREYJA_INPUTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(variants), path(depths) // legacy results may hold <s>.freyja_variants.tsv.gz / <s>.depths.gz

    output:
    // the modules that read the files (QC, consensus, mutations of interest, report) need plain text
    tuple val(meta), path("plain/${meta.id}.variants.tsv"), path("plain/${meta.id}.depths"), emit: plain

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir plain
    unpack() {
        case "\$1" in
            *.gz) gzip -dc "\$1" > "\$2" ;;
            *)    cp -L "\$1" "\$2" ;;
        esac
    }
    unpack ${variants} plain/${meta.id}.variants.tsv
    unpack ${depths} plain/${meta.id}.depths
    """

    stub:
    """
    mkdir plain
    touch plain/${meta.id}.variants.tsv plain/${meta.id}.depths
    """
}
