// Upstream canary for PLAN decision D6: runs plain `freyja plot` (no workaround) on the example
// aggregate that ships inside the pinned Freyja container and reports whether it produces a PDF.
// The image is taken from the production module so a Freyja bump moves the canary with it.
def freyjaImage() {
    return (file("${projectDir}/modules/local/freyja_update/main.nf").text =~ /'(quay\.io\/biocontainers\/freyja:[^']+)'/)[0][1]
}

process FREYJA_PLOT_CANARY {
    container freyjaImage()

    output:
    path "plain_plot_works.txt", emit: result

    script:
    """
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    data_dir=\$(python3 -c 'import freyja, os; print(os.path.join(os.path.dirname(freyja.__file__), "data"))')
    freyja plot "\$data_dir/aggregated_result.tsv" --output canary.pdf > canary.log 2>&1 || true
    if [ -s canary.pdf ]; then echo true; else echo false; fi > plain_plot_works.txt
    """
}
