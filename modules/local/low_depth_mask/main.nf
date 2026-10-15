process LOW_DEPTH_MASK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    // keep_bed: BED of the called-variant positions (CONSENSUS_FILTER_VCF `bed`); never masked.
    // Pass [] to mask purely on depth.
    tuple val(meta), path(depths), path(keep_bed)
    path fasta // reference: supplies the contig name and length (positions absent from the depths file have depth 0)

    output:
    tuple val(meta), path("${prefix}.mask.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: --min-depth <consensus_min_depth>
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def keep_arg = keep_bed ? "--keep-bed ${keep_bed}" : ''
    """
    # Reference contig name and length (single-record FASTA)
    contig=\$(awk '/^>/ { sub(/^>/, ""); print \$1; exit }' ${fasta})
    length=\$(awk '/^>/ { n++; next } n == 1 { total += length(\$0) } END { print total + 0 }' ${fasta})

    if [ -s ${depths} ]; then
        low_depth_mask.py \\
            --depths ${depths} \\
            --length \$length \\
            ${keep_arg} \\
            ${args} \\
            --output ${prefix}.mask.bed
    else
        # An empty depths file means the sample has no coverage at all (e.g. no mapped reads):
        # mask the whole reference. (low_depth_mask.py rejects an empty depths file by design.)
        printf '%s\\t0\\t%s\\n' "\$contig" "\$length" > ${prefix}.mask.bed
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mask.bed
    """
}
