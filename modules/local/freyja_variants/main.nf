process FREYJA_VARIANTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("${prefix}.variants.tsv"), emit: variants
    tuple val(meta), path("${prefix}.depths"),       emit: depths
    tuple val("${task.process}"), val('freyja'), eval('freyja --version | sed "s/.*version //"'), emit: versions_freyja, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # One call runs `samtools mpileup | ivar variants` and writes both files. ivar appends .tsv
    # to the --variants prefix (verified with Freyja 2.0.5); the glob below also copes with a
    # version that does not.
    freyja variants ${bam} \\
        --ref ${fasta} \\
        --variants ${prefix}.variants.raw \\
        --depths ${prefix}.depths \\
        ${args}

    found=\$(ls ${prefix}.variants.raw ${prefix}.variants.raw.tsv 2>/dev/null || true)
    if [ \$(echo "\$found" | grep -c .) -ne 1 ]; then
        echo "ERROR: expected exactly one of ${prefix}.variants.raw[.tsv] from freyja variants, found: '\$found'" >&2
        exit 1
    fi
    mv "\$found" ${prefix}.variants.tsv

    if [ ! -s ${prefix}.variants.tsv ] || [ ! -s ${prefix}.depths ]; then
        echo "ERROR: freyja variants produced an empty variants or depths file for ${meta.id}" >&2
        exit 1
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variants.tsv ${prefix}.depths
    """
}
