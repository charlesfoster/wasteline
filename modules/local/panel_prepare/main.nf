process PANEL_PREPARE {
    tag "panel"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    // panel_resolved: the resolved panel table (columns target,fasta,marker_contig,family at least; further columns are ignored).
    // fastas: every FASTA the table names, staged flat (the table may hold absolute host paths; a path that does not exist inside the
    // task is looked up by file name among the staged files). A FASTA shared by two targets is staged once.
    path panel_resolved
    path fastas

    output:
    path "panel.fasta", emit: fasta // all targets' contigs, named <target>|<contig> (original names with a single target)
    path "contig_map.tsv", emit: contig_map // panel_contig, target, contig, is_marker, family
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    panel_select.py prepare \\
        --panel ${panel_resolved} \\
        --fasta-dir . \\
        --out-fasta panel.fasta \\
        --out-map contig_map.tsv
    """

    stub:
    // the targets of the table (T1 and T2 when it is not a panel table, as in the module test), one contig c1 each
    """
    targets=\$(awk -F, 'NR > 1 && \$1 != "" { print \$1 }' ${panel_resolved} | head -n 20)
    case "\$(head -n 1 ${panel_resolved})" in target,*) ;; *) targets="T1 T2" ;; esac
    [ -n "\$targets" ] || targets="T1 T2"
    : > panel.fasta
    printf 'panel_contig\\ttarget\\tcontig\\tis_marker\\tfamily\\n' > contig_map.tsv
    for t in \$targets; do
        printf '>%s|c1\\nACGTACGT\\n' "\$t" >> panel.fasta
        printf '%s|c1\\t%s\\tc1\\ttrue\\tF\\n' "\$t" "\$t" >> contig_map.tsv
    done
    """
}
