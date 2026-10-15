process MUTATIONS_OF_INTEREST_CONCAT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    path tables, stageAs: 'moi/*' // one <sample>.mutations_of_interest.tsv per sample (collected)

    output:
    path "mutations_of_interest.tsv", emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # One header (all inputs must carry the same one); rows sorted by sample with the
    # mutation order inside each sample (= the order of the MOI CSV) preserved.
    first=\$(ls moi/* | head -n 1)
    head -n 1 "\$first" > mutations_of_interest.tsv
    for f in moi/*; do
        if [ "\$(head -n 1 "\$f")" != "\$(head -n 1 "\$first")" ]; then
            echo "ERROR: header of \$f differs from \$first" >&2
            exit 1
        fi
    done
    for f in moi/*; do tail -n +2 "\$f"; done | LC_ALL=C sort -s -t "\$(printf '\\t')" -k1,1 >> mutations_of_interest.tsv
    """

    stub:
    """
    touch mutations_of_interest.tsv
    """
}
