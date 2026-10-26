process QC_CONCAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(qc_rows, stageAs: 'qc/*'), path(amplicon_tables, stageAs: 'amplicons/*') // meta = the target (tkey); one <sample>.qc_row.tsv and one <sample>.amplicon_coverage.tsv ([] = none) per sample of the target

    output:
    tuple val(meta), path("qc_summary.tsv"),         emit: qc_summary
    tuple val(meta), path("amplicon_coverage.tsv"),  emit: amplicon_coverage, optional: true // absent when no amplicon tables were given (demix mode without a primer scheme)
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # qc_summary.tsv: validates headers, rejects duplicate samples, sorts by sample
    qc_summary.py --concat qc/* --output qc_summary.tsv

    # amplicon_coverage.tsv: one header (all inputs must carry the same one), rows sorted by
    # sample with the amplicon order inside each sample preserved
    if [ -d amplicons ] && [ -n "\$(ls -A amplicons)" ]; then
        first=\$(ls amplicons/* | head -n 1)
        head -n 1 "\$first" > amplicon_coverage.tsv
        for f in amplicons/*; do
            if [ "\$(head -n 1 "\$f")" != "\$(head -n 1 "\$first")" ]; then
                echo "ERROR: amplicon table header of \$f differs from \$first" >&2
                exit 1
            fi
        done
        for f in amplicons/*; do tail -n +2 "\$f"; done | LC_ALL=C sort -s -t "\$(printf '\\t')" -k1,1 >> amplicon_coverage.tsv
    fi
    """

    stub:
    """
    touch qc_summary.tsv
    if [ -d amplicons ] && [ -n "\$(ls -A amplicons)" ]; then touch amplicon_coverage.tsv; fi
    """
}
