process FREYJA_SELECT_CONTIG {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    // The variants TSV and depths file of ALL reference contigs (FREYJA_VARIANTS, or given files in mode demix). Multi-contig
    // references only (PLAN D25): Freyja demixes one contig, and its depths and variants must not hold the others.
    tuple val(meta), path(variants, stageAs: 'all_contigs/*'), path(depths, stageAs: 'all_contigs/*')
    val contig // the contig Freyja demixes (first word of its FASTA header)

    output:
    // the output names stay <id>.variants.tsv / <id>.depths whatever the input names were. Not published: the full files are
    tuple val(meta), path("${meta.id}.variants.tsv"), path("${meta.id}.depths"), emit: view
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # The same rows as 'freyja variants --refname ${contig}' writes (checked): the rows of one contig of the variants TSV (REGION,
    # column 1; the header line is kept) and the depths file (no header, contig in column 1). Fail loudly when the contig is
    # not there, rather than let Freyja demix nothing.
    awk -F'\\t' -v c='${contig}' 'NR == 1 || \$1 == c' ${variants} > ${meta.id}.variants.tsv
    awk -F'\\t' -v c='${contig}' '\$1 == c' ${depths} > ${meta.id}.depths

    if [ ! -s ${meta.id}.variants.tsv ] || [ ! -s ${meta.id}.depths ]; then
        echo "ERROR: contig '${contig}' has no rows in the depths file of ${meta.id} (or the variants file is empty): it is not a contig of the reference these files were made with" >&2
        exit 1
    fi
    """

    stub:
    """
    touch ${meta.id}.variants.tsv ${meta.id}.depths
    """
}
