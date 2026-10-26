//
// Subworkflow with functionality specific to the charlesfoster/wasteline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    mode              //  string: full | demix | aggregate | report | update_db
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet (not needed by every mode)
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = wastelineLogo(monochrome_logs)
    // `--help` lists only the common options: say how to see the rest (`--help <name>` is a single-parameter lookup)
    def after_text = (help && !help_full && !(help instanceof String && help != "true"))
        ? "Only the common options are shown above. Run with --help_full for every option, or --help <name> for one option.\n"
        : ""

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters (nf-schema cannot express "required depending on mode")
    //
    validateInputParameters(mode)

    //
    // Rows of the samplesheet given with --input, validated against the schema of the mode
    // (assets/schema_input*.json). Sample IDs are checked (pattern + uniqueness) and never rewritten.
    //   full:      [ meta, [ fastq_1, fastq_2 ] ]
    //   demix:     [ meta, variants, depths, fastp_json ]
    //   aggregate: [ meta, demix ]
    // Empty when no samplesheet is used (demix/aggregate from results directories, report, update_db):
    // the RESULTS_INPUT subworkflow produces the same rows from --previous_results / --demix_dir.
    //
    // The rows are built eagerly, so a bad `protocol` / `primer_scheme` cell stops the run here, naming the sample.
    def ch_samplesheet = channel.empty()
    def rows = []
    if (input && mode == 'full' && params.panel) {
        // A panel run (--panel): which targets a sample carries is only known after the typing, so the samples carry the pseudo
        // target 'panel' and the raw library cells; the per-target metas are made by PANEL_FRONT (sampleMeta with the target).
        // Everything that does not depend on the typing is checked here.
        def panel_samples = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
            .collect { meta, fastq_1, fastq_2 -> [ meta, [ fastq_1, fastq_2 ] ] }
        checkPanelSamples(panel_samples.collect { sample -> sample[0] }, targetsFromParams())
        rows = panel_samples.collect { sample -> [ panelSampleMeta(sample[0]), sample[1] ] }
    }
    else if (input && mode == 'full') {
        rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
            .collect { meta, fastq_1, fastq_2 -> [ sampleMeta(meta, mode), [ fastq_1, fastq_2 ] ] }
    }
    else if (input && mode == 'demix') {
        rows = samplesheetToList(input, "${projectDir}/assets/schema_input_demix.json")
            .collect { meta, variants, depths, fastp_json -> [ sampleMeta(meta, mode), variants, depths, fastp_json ?: [] ] }
    }
    else if (input && mode == 'aggregate') {
        rows = samplesheetToList(input, "${projectDir}/assets/schema_input_aggregate.json")
            .collect { meta, demix -> [ sampleMeta(meta, mode), demix ] }
    }
    if (rows) {
        if (!(mode == 'full' && params.panel)) {
            checkSampleBeds(rows.collect { row -> row[0] })
        }
        ch_samplesheet = channel.fromList(rows)
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                [],
            )
        }

        finaliseRunInfo(outdir)
        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// ASCII logo (assets/logo.txt) with version line; ANSI colour unless monochrome_logs
//
def wastelineLogo(monochrome_logs) {
    def logo = file("${projectDir}/assets/logo.txt").text
    def version = "wasteline v${workflow.manifest.version}"
    if (monochrome_logs || System.getenv('NO_COLOR')) {
        return "\n${logo}\n  ${version}\n"
    }
    return "\n\033[0;33m${logo}\033[0m\n  \033[1m${version}\033[0m\n"
}

//
// Container image per process. workflow.container is empty for images declared inside modules,
// so fall back to the `container` column of the execution trace.
//
def runContainers(outdir) {
    def containers = workflow.container ? workflow.container.collectEntries { key, value -> [ (key.toString()): value.toString() ] } : [:]
    def trace = file("${outdir}/pipeline_info/execution_trace_${params.trace_report_suffix}.txt")
    if (!containers && trace.exists()) {
        def lines = trace.readLines()
        def header = lines[0].split('\t') as List
        def process_col = header.indexOf('process')
        def container_col = header.indexOf('container')
        if (process_col >= 0 && container_col >= 0) {
            lines.drop(1).each { line ->
                def cols = line.split('\t', -1)
                if (cols.size() > container_col && cols[container_col] && cols[container_col] != '-') {
                    containers[cols[process_col].tokenize(':').last()] = cols[container_col]
                }
            }
        }
    }
    return containers.sort()
}

//
// Add what is only known when the run ends (container images per process, completion time and
// outcome) to pipeline_info/run_info.json, which PROVENANCE wrote during the run
//
def finaliseRunInfo(outdir) {
    // a panel run has one run_info.json per target (pipeline_info/<target>/), a single-target run one in pipeline_info/
    def subs = params.panel ? targetsFromParams().collect { t -> t.summary_sub } : ['']
    subs.each { sub ->
        def run_info = file("${outdir}/pipeline_info${sub}/run_info.json")
        if (run_info.exists()) {
            def info = new groovy.json.JsonSlurper().parseText(run_info.text)
            info.containers = runContainers(outdir)
            info.run_completed = workflow.complete?.toString()
            info.duration = workflow.duration?.toString()
            info.success = workflow.success
            run_info.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
        }
    }
}

//
// The target of the run (the reference a sample is analysed against) from the effective parameters, i.e. after the
// conf/pathogens/*.config defaults, --fasta / --gff and --whole_genome. With one target every sample carries this map.
// A target map has the same keys whether it comes from here (single-target run) or from targetByName() (panel run):
//   name, freyja_pathogen, fasta, gff, freyja_contig, marker_contig, family, primer_scheme, primer_scheme_default,
//   primer_schemes_bundled, lineage_groups_bundled, freyja_container_pathogen, barcodes, panel, sample_sub, summary_sub
// panel is false here; sample_sub / summary_sub are the path suffixes of per-sample / run-level outputs ('' for a
// single-target run, which keeps the legacy layout, and '/<name>' for a panel target).
//
def targetFromParams() {
    def pathogen = params.pathogen ?: params.freyja_container_pathogen
    return [
        name: pathogen,
        freyja_pathogen: pathogen,
        fasta: params.fasta,
        gff: params.gff,
        freyja_contig: params.freyja_contig,
        marker_contig: null,
        family: null,
        primer_scheme: null,
        primer_scheme_default: params.primer_scheme_default,
        primer_schemes_bundled: params.primer_schemes_bundled,
        lineage_groups_bundled: params.lineage_groups_bundled,
        freyja_container_pathogen: params.freyja_container_pathogen,
        barcodes: null,
        panel: false,
        sample_sub: '',
        summary_sub: '',
    ]
}

//
// The key that identifies a sample within its target: the sample ID while a run has one target, ID__target in a panel run
//
def targetKey(meta) {
    return meta.target?.panel ? "${meta.id}__${meta.target.name}" : meta.id
}

//
// The grouping key of the run-level steps (every step that summarises all samples of ONE target): [ id: target name, target ]
//
def tkey(meta) {
    return [ id: meta.target.name, target: meta.target ]
}

//
// The pseudo target of a sample before the typing of a panel run (mode full, --panel): the sample is not yet assigned to a
// target, so fastp, the mapping to the panel and the typing carry this one. Its outputs are published directly under
// samples/<id>/ (sample_sub is empty). PANEL_FRONT replaces it by the real target of each (sample, target) pair.
//
def panelTarget() {
    return [
        name: 'panel',
        freyja_pathogen: null,
        fasta: null,
        gff: null,
        freyja_contig: null,
        marker_contig: null,
        family: null,
        primer_scheme: null,
        primer_scheme_default: null,
        primer_schemes_bundled: null,
        lineage_groups_bundled: null,
        freyja_container_pathogen: null,
        barcodes: null,
        panel: true,
        sample_sub: '',
        summary_sub: '',
    ]
}

//
// Does a results directory come from a panel run (several targets)? Panel runs write summary/targets.tsv.
//
def resultsAreMultiTarget(dir) {
    return file("${dir}/summary/targets.tsv").exists()
}

//
// Panels (several targets in one run, --panel). assets/panels/registry.csv holds the facts of every target wasteline can
// analyse (the Freyja pathogens that have barcodes); a panel file (assets/panels/<name>.csv or the user's own) lists the
// targets of a run and may override any registry column of a target (and set primer_scheme).
//
def panelRegistryColumns() {
    return [ 'target', 'fasta', 'gff', 'freyja_pathogen', 'freyja_contig', 'marker_contig', 'family', 'primer_scheme_default',
        'primer_schemes_bundled', 'lineage_groups_bundled', 'freyja_container_pathogen', 'barcodes' ]
}

def splitCsvLine(line) {
    return line.split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/, -1).collect { cell -> cell.trim().replaceAll(/^"|"$/, '') }
}

def readCsvRows(csv) {
    def lines = file(csv, checkIfExists: true).readLines().findAll { line -> line.trim() }
    def header = splitCsvLine(lines[0])
    return lines.drop(1).collect { line ->
        def cells = splitCsvLine(line)
        header.indexed().collectEntries { i, name -> [ (name): i < cells.size() && cells[i] ? cells[i] : null ] }
    }
}

def registryRows() {
    return readCsvRows("${projectDir}/assets/panels/registry.csv")
}

//
// A path cell of the registry (relative to base, the pipeline directory) or of a panel file (base null: relative to the launch
// directory) as an absolute path. A barcodes cell may also be a keyword (latest, container, a date): only paths are changed.
//
def panelPath(key, value, base) {
    if (!value || value.startsWith('/')) {
        return value
    }
    def is_path = key != 'barcodes' || value.contains('/') || value.endsWith('.csv') || value.endsWith('.feather')
    if (!is_path) {
        return value
    }
    return base ? "${base}/${value}" : file(value).toAbsolutePath().toString()
}

//
// A target of the registry as a target map (see targetFromParams()); overrides replace registry columns (empty = keep).
// Paths of the registry are relative to the pipeline directory, paths of overrides to the launch directory.
//
def targetByName(name, overrides = [:]) {
    def row = registryRows().find { r -> r.target == name }
    if (!row) {
        error("Unknown target '${name}': the targets wasteline can analyse are ${registryRows().collect { r -> r.target }.join(', ')} (assets/panels/registry.csv)")
    }
    def path_keys = [ 'fasta', 'gff', 'primer_schemes_bundled', 'lineage_groups_bundled', 'barcodes' ]
    def values = row.collectEntries { key, value -> [ (key): key in path_keys ? panelPath(key, value, projectDir.toString()) : value ] }
    overrides.findAll { key, value -> value != null && value != '' && key != 'target' }.each { key, value ->
        values[key] = key in path_keys ? panelPath(key, value.toString(), null) : value
    }
    def panel = overrides.containsKey('panel') ? overrides.panel : true
    return [
        name: name,
        freyja_pathogen: values.freyja_pathogen,
        fasta: values.fasta,
        gff: values.gff,
        freyja_contig: values.freyja_contig,
        marker_contig: values.marker_contig,
        family: values.family,
        primer_scheme: values.primer_scheme ?: null,
        primer_scheme_default: values.primer_scheme_default,
        primer_schemes_bundled: values.primer_schemes_bundled,
        lineage_groups_bundled: values.lineage_groups_bundled,
        freyja_container_pathogen: values.freyja_container_pathogen,
        barcodes: values.barcodes,
        panel: panel,
        sample_sub: panel ? "/${name}" : '',
        summary_sub: panel ? "/${name}" : '',
    ]
}

//
// The panel file of --panel: a built-in name (assets/panels/<name>.csv, e.g. rsv, flu) or the path to a panel CSV
//
def panelFile(spec) {
    def name = spec.toString()
    def builtin = file("${projectDir}/assets/panels/${name}.csv")
    if (!name.contains('/') && !name.endsWith('.csv') && name != 'registry' && builtin.exists()) {
        return builtin
    }
    def custom = file(name)
    if (!custom.isFile()) {
        def builtins = file("${projectDir}/assets/panels").listFiles().collect { f -> f.name }.findAll { n -> n.endsWith('.csv') && n != 'registry.csv' }.collect { n -> n.substring(0, n.length() - 4) }.sort()
        error("--panel '${name}' is neither a built-in panel (${builtins.join(', ')}) nor an existing panel CSV file (a relative path is relative to the directory you launch Nextflow from)")
    }
    return custom
}

//
// The targets of the run: [ targetFromParams() ] without --panel, otherwise the rows of the panel file as target maps,
// validated (a panel is a closed set: nothing is looked up later that could be missing)
//
def targetsFromParams() {
    if (!params.panel) {
        return [ targetFromParams() ]
    }
    def panel_file = panelFile(params.panel)
    def rows = readCsvRows(panel_file)
    def known = panelRegistryColumns() + [ 'primer_scheme' ]
    def columns = rows ? rows[0].keySet() : []
    if (!('target' in columns)) {
        error("Panel file '${panel_file}' has no 'target' column (columns: ${columns.join(', ') ?: 'none'}); it lists one target per row, by name as in assets/panels/registry.csv")
    }
    def unknown_columns = columns.findAll { c -> !(c in known) }
    if (unknown_columns) {
        error("Panel file '${panel_file}' has unknown column(s) ${unknown_columns.join(', ')}; allowed: ${known.join(', ')}")
    }
    def names = rows.collect { r -> r.target }
    if (names.any { n -> !n }) {
        error("Panel file '${panel_file}' has a row without a target name")
    }
    def repeated = names.findAll { n -> names.count(n) > 1 }.unique()
    if (repeated) {
        error("Panel file '${panel_file}' lists target(s) ${repeated.join(', ')} more than once")
    }
    def targets = rows.collect { r -> targetByName(r.target, r) }
    validatePanelTargets(targets)
    return targets
}

//
// A panel is mapped as one concatenated reference, so contig names must be unique across it. Two targets may share a FASTA
// only when they are the same family with the same marker contig (H5Nx and H5Nx-cattle: one reference, two sets of barcodes).
//
def validatePanelTargets(targets) {
    def contig_owner = [:]   // contig -> target
    def fasta_owner = [:]    // real path of a FASTA -> [ target, contigs ]
    targets.each { t ->
        if (!t.fasta) {
            error("Panel target '${t.name}' has no reference FASTA")
        }
        def real = file(t.fasta, checkIfExists: true).toRealPath().toString()
        if (!fasta_owner.containsKey(real)) {
            def contigs = fastaContigs(file(t.fasta))
            def dup = contigs.findAll { c -> contigs.count(c) > 1 }.unique()
            if (dup) {
                error("Panel target '${t.name}': the reference FASTA '${t.fasta}' has several sequences named ${dup.join(', ')}; contig names must be unique")
            }
            contigs.each { c ->
                def other = contig_owner[c]
                if (other) {
                    def hint = (t.name.startsWith('H5Nx') && other.name.startsWith('H5Nx')) ? ' (H5Nx and H5Nx-cattle must use one shared FASTA)' : ''
                    error("Panel targets '${other.name}' and '${t.name}' both have a contig named '${c}' (in '${other.fasta}' and '${t.fasta}'): contig names must be unique across a panel because the reads are mapped to all references at once${hint}")
                }
                contig_owner[c] = t
            }
            fasta_owner[real] = [ target: t, contigs: contigs ]
        }
        else if (!t.family || t.family != fasta_owner[real].target.family || !t.marker_contig || t.marker_contig != fasta_owner[real].target.marker_contig) {
            error("Panel targets '${fasta_owner[real].target.name}' and '${t.name}' use the same reference FASTA '${t.fasta}': that is only allowed for targets of one family with the same marker contig (influenza H5Nx and H5Nx-cattle share a reference)")
        }
        [ 'freyja_contig', 'marker_contig' ].each { key ->
            if (t[key] && !(t[key] in fasta_owner[real].contigs)) {
                error("Panel target '${t.name}': ${key} '${t[key]}' is not a sequence of its reference FASTA '${t.fasta}' (${fasta_owner[real].contigs.take(10).join(', ')})")
            }
        }
    }
}

//
// Sample meta map shared by every mode. In modes that map or demix (full, demix) it also carries the library
// settings of the sample (PLAN D20), see libraryMeta().
//
def sampleMeta(meta, mode, target = targetFromParams()) {
    def sample = [
        id: meta.id,
        sample_type: meta.sample_type ?: 'sample',
        collection_date: meta.collection_date ?: null,
        site: meta.site ?: null,
        viral_load: meta.viral_load ?: null,
        target: target,
    ]
    sample.key = targetKey(sample)
    if (mode in ['full', 'demix']) {
        return sample + libraryMeta(meta.id, meta.protocol, meta.primer_scheme, mode, target)
    }
    return sample
}

//
// The meta of a sample before the typing of a panel run: no target-specific library settings (they depend on the target and are
// resolved per (sample, target) by PANEL_FRONT from the raw cells kept here). key = id: one item per sample at this stage.
//
def panelSampleMeta(meta) {
    return [
        id: meta.id,
        sample_type: meta.sample_type ?: 'sample',
        collection_date: meta.collection_date ?: null,
        site: meta.site ?: null,
        viral_load: meta.viral_load ?: null,
        protocol_cell: meta.protocol ?: null,
        primer_scheme_cell: meta.primer_scheme ?: null,
        target: panelTarget(),
        key: meta.id,
    ]
}

//
// Start-up checks of the samples of a panel run (mode full): the library cells (protocol, primer scheme) must resolve for EVERY
// target of the panel, because the typing decides later which targets a sample is analysed for, and the primer BEDs a target
// will use must be on its reference (checkSampleBeds, once per target). A samplesheet primer_scheme cell names one scheme,
// which belongs to one target: with several targets it is refused (give each target's scheme in the panel CSV primer_scheme column).
//
def checkPanelSamples(samples, panel_targets) {
    samples.each { sample ->
        if (sample.primer_scheme && panel_targets.size() > 1) {
            error("Sample '${sample.id}': the samplesheet primer_scheme '${sample.primer_scheme}' cannot be used with a panel of ${panel_targets.size()} targets (${panel_targets.collect { t -> t.name }.join(', ')}): a primer scheme belongs to one reference. Leave the cell empty and give each target its scheme in the panel CSV (column primer_scheme), or run the targets one by one with --pathogen.")
        }
        if (sample.target) {
            log.warn("Sample '${sample.id}': the samplesheet target column is not used by --mode full with --panel (the typing decides which targets a sample is analysed for) and is ignored")
        }
    }
    panel_targets.each { t ->
        checkSampleBeds(samples.collect { sample -> sampleMeta(sample, 'full', t) }, t)
    }
}

//
// THE effective library settings of one sample, computed here and nowhere else (conf/modules.config spells
// out the naming conditions it needs, because a config file cannot call a function). `protocol` and
// `primer_scheme` are the optional samplesheet cells; empty = the run defaults (--protocol, --primer_bed /
// --primer_scheme, --skip_primer_trimming, --duplicates).
//   protocol         amplicon | capture
//   primer_trimming  true only for amplicon reads that are not --skip_primer_trimming
//   duplicates       mark | none (--duplicates auto = mark for capture, none for amplicon)
//   primer_scheme    bundled scheme name, or 'custom' (then primer_bed is the BED path); null when the sample has no
//                    amplicon-level analysis: not primer-trimmed, or demix mode without a named scheme (the scheme
//                    of old results is unknown). Full mode with nothing named uses the bundled default (D14).
//
def libraryMeta(sample, protocol_cell, scheme_cell, mode, target = targetFromParams()) {
    def protocol = protocol_cell ?: params.protocol
    if (!(protocol in ['amplicon', 'capture'])) {
        error("Sample '${sample}': protocol '${protocol}' is not one of amplicon, capture")
    }
    def settings = protocolSettings(protocol)
    def scheme = null
    def bed = null
    if (settings.primer_trimming) {
        // order: the sample's own cell, the primer_scheme of the target in the panel CSV, --primer_bed, --primer_scheme, the target's default
        def cell = scheme_cell ?: target.primer_scheme
        def who = scheme_cell ? "Sample '${sample}': primer_scheme" : "Panel target '${target.name}': primer_scheme"
        if (cell) {
            def spec = cell.toString().trim()
            if (spec.contains('/') || spec.endsWith('.bed')) {
                bed = file(spec).toAbsolutePath().toString()
                if (!file(bed).exists()) {
                    error("${who} BED file '${spec}' does not exist (a relative path is relative to the directory you launch Nextflow from)")
                }
                scheme = 'custom'
            }
            else {
                def chosen = chooseScheme(spec, "${who} '${spec}'", target)
                scheme = chosen.scheme
                bed = chosen.bed
            }
        }
        else if (params.primer_bed) {
            scheme = 'custom'
            bed = params.primer_bed.toString()
        }
        else if (params.primer_scheme) {
            def chosen = chooseScheme(params.primer_scheme.toString(), "--primer_scheme '${params.primer_scheme}'", target)
            scheme = chosen.scheme
            bed = chosen.bed
        }
        else if (mode == 'full') {
            if (!defaultPrimerScheme(target)) {
                error("Sample '${sample}' has no primer scheme: ${noDefaultSchemeMessage(target)}")
            }
            scheme = defaultPrimerScheme(target)
        }
    }
    return settings + [ primer_scheme: scheme, primer_bed: bed ]
}

//
// Primer schemes by name (PLAN D26). A name is looked up in --primer_scheme_dir first (the user's own <name>.bed files win),
// then in the directory of schemes bundled for the pathogen (target.primer_schemes_bundled, set by conf/pathogens/*.config:
// assets/primer_schemes for SARS-CoV-2, assets/primer_schemes/<profile> for the others, null for pathogens without any).
//
def schemeNamesIn(dir) {
    if (!dir) {
        return []
    }
    def d = file(dir)
    return d.isDirectory() ? d.listFiles().findAll { f -> f.isFile() && f.name.endsWith('.bed') }.collect { f -> f.name.substring(0, f.name.length() - 4) }.sort() : []
}

def bundledPrimerSchemes(target = targetFromParams()) {
    return schemeNamesIn(target.primer_schemes_bundled)
}

def userPrimerSchemes() {
    return schemeNamesIn(params.primer_scheme_dir)
}

//
// What the bundled schemes of the other pathogens are called, to point a user who named one of them for the wrong pathogen
//
def schemesOfOtherPathogens(name) {
    def top = file("${projectDir}/assets/primer_schemes")
    return top.listFiles().findAll { d -> d.isDirectory() && file("${d}/${name}.bed").exists() }.collect { d -> d.name }.sort()
}

//
// A scheme name to [ scheme, bed ] as libraryMeta() stores it: a bundled scheme keeps its name (its BED is found by
// sourceBedFor), a scheme of --primer_scheme_dir is a custom BED (the path). Stops with the names available otherwise.
//
def chooseScheme(name, who, target = targetFromParams()) {
    if (name in userPrimerSchemes()) {
        return [ scheme: 'custom', bed: file("${params.primer_scheme_dir}/${name}.bed").toAbsolutePath().toString() ]
    }
    if (name in bundledPrimerSchemes(target)) {
        return [ scheme: name, bed: null ]
    }
    def others = schemesOfOtherPathogens(name)
    def hint = others ? " ('${name}' is a scheme bundled for another pathogen, in assets/primer_schemes/${others.join(', ')}/; it belongs with that pathogen's reference)" : ''
    error("${who} is not a primer scheme available for ${target.panel ? 'panel target' : '--pathogen'} ${target.name} and not a path to a BED file${hint}. Bundled: ${bundledPrimerSchemes(target).join(', ') ?: 'none'}" +
        (params.primer_scheme_dir ? "; in --primer_scheme_dir '${params.primer_scheme_dir}': ${userPrimerSchemes().join(', ') ?: 'no .bed files'}" : "; give your own with --primer_bed <file> or --primer_scheme_dir <dir of <name>.bed files>"))
}

//
// The message for a run that needs a primer scheme and has none (a pathogen without a default scheme)
//
def noDefaultSchemeMessage(target = targetFromParams()) {
    def bundled = bundledPrimerSchemes(target)
    if (target.panel) {
        return "panel target ${target.name} has no default primer scheme. Give it one in the panel CSV (column primer_scheme: a bundled scheme" +
            (bundled ? " of ${bundled.join(', ')}" : " (none is bundled for ${target.name})") + " or a BED path), or pass --primer_bed <BED whose first column is the reference sequence name>, or say that the reads are not primer-trimmed with --protocol capture / --skip_primer_trimming. See docs/other_pathogens.md."
    }
    return "--pathogen ${target.name} has no default primer scheme (only the pipeline's default pathogen has one). Name one with --primer_scheme (or the samplesheet primer_scheme column)" +
        (bundled ? " (bundled for ${target.name}: ${bundled.join(', ')}; " : " (none is bundled for ${target.name}; ") +
        "your own: --primer_scheme_dir <dir of <name>.bed files>), pass --primer_bed <BED whose first column is the reference sequence name>, or say that the reads are not primer-trimmed with --protocol capture / --skip_primer_trimming. See docs/other_pathogens.md."
}

//
// protocolSettings: the part of the library settings that depends only on the protocol
//
def protocolSettings(protocol) {
    def capture = protocol == 'capture'
    return [
        protocol: protocol,
        primer_trimming: !capture && !params.skip_primer_trimming,
        duplicates: params.duplicates == 'auto' ? (capture ? 'mark' : 'none') : params.duplicates,
    ]
}

//
// The primer BED of a sample with a primer_scheme, as given (the user's file or the bundled scheme of the pathogen)
//
def sourceBedFor(meta, target = targetFromParams()) {
    if (meta.primer_bed) {
        return file(meta.primer_bed, checkIfExists: true)
    }
    return file("${target.primer_schemes_bundled}/${meta.primer_scheme}.bed", checkIfExists: true)
}

//
// ... as it is used by the processes: restricted to the reference's contig(s) (reconcileBed, PLAN D24)
//
def primerBedFor(meta, target = targetFromParams()) {
    return reconcileBed(sourceBedFor(meta, target), target).bed
}

//
// A primer BED against the reference (PLAN D24). Whole-genome panels (influenza) list the primers of every segment, while
// the reference may be one segment (HA): the lines of contigs that are not in the reference FASTA are dropped, and the
// filtered copy (the same file name, in <workDir>/wasteline_primer_beds/<hash of its content>/) is what trimming, amplicon
// coverage and provenance use. A BED whose contigs are all in the reference is returned as it is (the same file, byte for
// byte). Returns [bed: Path, no_match: boolean, dropped_lines: n, dropped_contigs: [names]]. Comment, track and browser
// lines are kept. A BED none of whose lines is on the reference is returned unchanged with no_match = true (the callers
// refuse it with their own message).
//
def reconcileBed(bed, target = targetFromParams()) {
    if (!target.fasta) {
        return [ bed: bed, no_match: false, dropped_lines: 0, dropped_contigs: [] ]
    }
    def contigs = fastaContigs(file(target.fasta, checkIfExists: true))
    def keep = []
    def kept = 0
    def dropped_lines = 0
    def dropped_contigs = []
    bed.readLines().each { line ->
        def record = line.trim() && !line.startsWith('#') && !line.startsWith('track') && !line.startsWith('browser')
        def contig = record ? line.split('\t')[0] : null
        if (!record) {
            keep << line
        }
        else if (contig in contigs) {
            keep << line
            kept += 1
        }
        else {
            dropped_lines += 1
            if (!(contig in dropped_contigs)) {
                dropped_contigs << contig
            }
        }
    }
    if (dropped_lines == 0 || kept == 0) {
        return [ bed: bed, no_match: kept == 0 && dropped_lines > 0, dropped_lines: dropped_lines, dropped_contigs: dropped_contigs ]
    }
    def text = keep.join('\n') + '\n'
    def digest = java.security.MessageDigest.getInstance('SHA-256').digest(text.getBytes('UTF-8')).collect { b -> String.format('%02x', b) }.join('').substring(0, 12)
    def filtered = file("${workflow.workDir}/wasteline_primer_beds/${digest}/${bed.name}")
    if (!filtered.exists()) {
        filtered.parent.mkdirs()
        filtered.text = text
    }
    return [ bed: filtered, no_match: false, dropped_lines: dropped_lines, dropped_contigs: dropped_contigs ]
}

//
// The one warning about primer BED lines dropped because their contig is not in the reference (PLAN D24)
//
def warnDroppedBedLines(reports, fasta = params.fasta) {
    def lines = reports.sum { report -> report.dropped_lines } ?: 0
    if (lines) {
        def contigs = reports.collectMany { report -> report.dropped_contigs }.unique()
        def reference_contigs = fastaContigs(file(fasta))
        def expected = reference_contigs.size() > 1
            ? "The reference has ${reference_contigs.size()} sequences (${reference_contigs.take(10).join(', ')}${reference_contigs.size() > 10 ? ', ...' : ''}); lines on any other contig cannot be matched to it."
            : "Expected for a whole-genome panel run against one segment (for example influenza HA)."
        log.warn("Primer BED: ${lines} line(s) on ${contigs.size()} contig(s) that are not in the reference FASTA '${file(fasta).name}' (${contigs.take(10).join(', ')}${contigs.size() > 10 ? ', ...' : ''}) were dropped; the reference's own lines are used for trimming and amplicon coverage. ${expected}")
    }
}

//
// What run_info.json records about the library settings of the samples (modes full and demix; [] elsewhere):
//   protocol / primer_trimming / duplicates   the common value, or 'mixed' when the samples differ
//                                             (primer_trimming and duplicates are only reported for mode full)
//   sample_libraries                          per sample: protocol, primer_scheme ('custom: <BED name>' for a BED
//                                             file, null when there is no amplicon-level analysis), primer_trimming, duplicates
//
def commonValue(metas, key) {
    def values = metas.collect { meta -> meta[key] }.unique()
    return values.size() == 1 ? values[0] : 'mixed'
}

def libraryInfo(mode, metas) {
    if (!(mode in ['full', 'demix']) || !metas) {
        return [ protocol: params.protocol, primer_trimming: null, duplicates: null ]
    }
    def libraries = metas.sort { meta -> meta.id }.collectEntries { meta ->
        [ (meta.id): [
            protocol: meta.protocol,
            primer_scheme: meta.primer_bed ? "custom: ${file(meta.primer_bed).name}" : meta.primer_scheme,
            primer_trimming: meta.primer_trimming,
            duplicates: meta.duplicates,
        ] ]
    }
    def trimming = commonValue(metas, 'primer_trimming')
    return [
        protocol: commonValue(metas, 'protocol'),
        primer_trimming: mode == 'full' && trimming != 'mixed' ? trimming : null,
        duplicates: mode == 'full' ? commonValue(metas, 'duplicates') : null,
        sample_libraries: libraries,
    ]
}

//
// Every distinct primer BED of the samples must be on a contig of the reference or ivar trim silently trims nothing
//
def checkSampleBeds(metas, target = targetFromParams()) {
    def beds = metas.findAll { meta -> meta.primer_scheme }.collect { meta -> sourceBedFor(meta, target) }.unique()
    if (beds) {
        def fasta_contigs = fastaContigs(file(target.fasta, checkIfExists: true))
        def reports = beds.collect { bed ->
            def report = reconcileBed(bed, target)
            if (report.no_match) {
                error("Primer BED '${bed}' refers to contig(s) ${bedContigs(bed)} that are not in the reference FASTA '${target.fasta}' (which has ${fasta_contigs}). " +
                    (target.panel ? "Panel target ${target.name}: give it a matching primer scheme in the panel CSV (column primer_scheme), not a run-wide --primer_bed/--primer_scheme meant for another reference." : "Use a matching primer scheme or --fasta."))
            }
            report
        }
        warnDroppedBedLines(reports, target.fasta)
    }
}

//
// Freyja pathogens (`freyja update --pathogen` choices) that wasteline cannot run, with the reason and what would be
// needed (PLAN D22); null for every other name. The FASTA check in validateReferenceContigs accepts a multi-sequence
// reference only together with --freyja_contig. The influenza pathogens are supported on the HA segment (PLAN D24) or, with
// --whole_genome, on all eight segments with HA demixed (PLAN D25).
//
def unsupportedPathogenReason(pathogen) {
    if (pathogen == 'MTB') {
        return "Mycobacterium tuberculosis is a 4.4 Mb bacterial genome with feather-only Freyja barcodes. wasteline is built and tested for viral genomes of 10-200 kb (per-position depth tables, QC, consensus and the dropout/mutation tables are not sized for 4.4 Mb), and nothing has been validated on bacterial data. Supporting it needs those steps scaled and checked on real MTB data."
    }
    return null
}

//
// Comma-separated path list (--previous_results) as a list of strings
//
def splitPaths(spec) {
    return spec ? spec.toString().split(',').collect { entry -> entry.trim() }.findAll { entry -> entry } : []
}

def requireDirectory(option, value) {
    def dir = file(value)
    if (!dir.exists()) {
        error("${option} '${value}' does not exist")
    }
    if (!dir.isDirectory()) {
        error("${option} '${value}' is not a directory")
    }
    return dir
}

//
// --panel (several targets in one run): the panel replaces --pathogen and the reference options, which are therefore refused
// together with it, and its file is validated here (targetsFromParams) so a mistake stops the run before any work is done.
// --targets selects targets of the results that modes demix, aggregate and report read.
//
def validatePanelParameters(mode) {
    if (!params.panel) {
        if (params.targets) {
            log.warn("--targets selects targets of panel results: it is only used together with results of a --panel run and is ignored here")
        }
        return
    }
    def clashes = [
        '--pathogen': params.pathogen,
        '--fasta': params.fasta,
        '--gff': params.gff,
        '--freyja_contig': params.freyja_contig,
        '--whole_genome': params.whole_genome,
    ].findAll { _name, value -> value }.keySet()
    if (params.pathogen) {
        // --fasta, --gff and --freyja_contig are then (also) set by the pathogen's profile: name the cause
        clashes = clashes.findAll { name -> name in ['--pathogen', '--whole_genome'] }
    }
    if (clashes) {
        error("--panel cannot be combined with ${clashes.join(', ')}: a panel names its targets (assets/panels/registry.csv gives each its reference and Freyja pathogen), so --pathogen, --fasta, --gff, --freyja_contig and --whole_genome are not used with it. Give a panel CSV with the targets you want instead (see docs/usage.md, Panels).")
    }
    def targets = targetsFromParams()
    if (mode in ['demix', 'aggregate', 'report']) {
        error("--panel ${params.panel} (${targets.collect { t -> t.name }.join(', ')}) is valid, but panel results cannot be re-demixed/aggregated/reported yet (--mode ${mode}). Panel runs are available in --mode full; to redo a single target use --mode ${mode} on that target's results with --pathogen <target>.")
    }
    if (mode == 'full') {
        // the Freyja dataset of every target is resolved before any work is done (a registry row's own barcodes beat --barcodes)
        targets.each { t ->
            def spec = t.barcodes ?: params.barcodes
            if (spec.toString().trim() == 'container' && t.freyja_pathogen != t.freyja_container_pathogen) {
                error("Panel target ${t.name}: --barcodes container is only available for the pathogen whose barcodes the Freyja container bundles (SARS-CoV-2); the container has none for ${t.freyja_pathogen}. Use --barcodes latest, or give the target a barcodes file or dataset directory in the panel CSV (column barcodes).")
            }
            resolveBarcodes(spec, dbDir(), '--barcodes', t)
        }
        if (params.compare_barcodes) {
            // --compare_barcodes is one spec for every target: it must be valid for each (version stamps resolve per target pathogen)
            def compare = params.compare_barcodes.toString().trim()
            def wrong_pathogen = targets.findAll { t -> t.freyja_pathogen != t.freyja_container_pathogen }.collect { t -> t.name }
            if (compare == 'container' && wrong_pathogen) {
                error("--compare_barcodes container is only available for the pathogen whose barcodes the Freyja container bundles (SARS-CoV-2), but the panel also has ${wrong_pathogen.join(', ')}. Use --compare_barcodes latest or a barcode version (YYYY-MM-DD) to compare every target.")
            }
            if (targets.size() > 1 && resolveBarcodes(compare, dbDir(), '--compare_barcodes', targets[0]).source == 'custom') {
                error("--compare_barcodes ${compare} is a dataset directory or barcodes file, which belongs to one pathogen, but the panel has ${targets.size()} targets (${targets.collect { t -> t.name }.join(', ')}). Use --compare_barcodes latest or a barcode version (YYYY-MM-DD) to compare every target.")
            }
            targets.each { t -> resolveBarcodes(compare, dbDir(), '--compare_barcodes', t) }
        }
        if (params.annotate_variants) {
            def without_gff = targets.findAll { t -> !t.gff }.collect { t -> t.name }
            if (without_gff) {
                error("--annotate_variants needs the reference annotation: panel target(s) ${without_gff.join(', ')} have no GFF3 (give one in the panel CSV, column gff)")
            }
        }
    }
}

//
// Check and validate pipeline parameters for the mode. Only `outdir` is required by the schema for
// every mode; what else is needed depends on the mode and is checked here with a helpful message.
//
def validateInputParameters(mode) {
    def modes = ['full', 'demix', 'aggregate', 'report', 'update_db']
    if (!(mode in modes)) {
        error("--mode '${mode}' is not one of ${modes.join(', ')}")
    }
    def sources = [
        '--input': params.input,
        '--previous_results': params.previous_results,
        '--demix_dir': params.demix_dir,
    ].findAll { _name, value -> value }

    def accepted = [
        full: ['--input'],
        demix: ['--input', '--previous_results'],
        aggregate: ['--input', '--demix_dir', '--previous_results'],
        report: [],
        update_db: [],
    ][mode]
    def hint = [
        full: "--mode full needs a FASTQ samplesheet: --input samplesheet.csv (columns sample,fastq_1,fastq_2[,sample_type][,collection_date][,site][,protocol][,primer_scheme]; make one with bin/make_samplesheet.py <fastq_dir>)",
        demix: "--mode demix needs existing Freyja variants and depths: --input samplesheet.csv (columns sample,variants,depths[,fastp_json][,sample_type][,collection_date][,site][,protocol][,primer_scheme]) or --previous_results <results dir> (new or legacy Snakemake layout; never modified)",
        aggregate: "--mode aggregate needs existing .demix files: --input samplesheet.csv (columns sample,demix[,sample_type][,collection_date][,site]), --demix_dir <dir of *.demix files>, or --previous_results <dir>[,<dir>...]",
        report: "--mode report needs --results_dir <results directory of a previous run> (it takes no samplesheet)",
        update_db: "--mode update_db takes no input; it downloads the latest Freyja data into --db_dir",
    ][mode]

    // report and update_db take no input at all: a stray --input (e.g. the default of a test profile) is ignored
    if (mode in ['report', 'update_db'] && params.input) {
        log.warn("--input is not used by --mode ${mode} and is ignored")
        sources.remove('--input')
    }
    def unexpected = sources.keySet().findAll { name -> !(name in accepted) }
    if (unexpected) {
        error("${unexpected.join(', ')} cannot be used with --mode ${mode}. ${hint}")
    }
    if (accepted && sources.size() == 0) {
        error("No input given for --mode ${mode}. ${hint}")
    }
    if (sources.size() > 1) {
        error("Give only one of ${sources.keySet().join(', ')} with --mode ${mode}. ${hint}")
    }

    validatePanelParameters(mode)

    // The pathogen (PLAN D22): the Freyja pathogens that need more than a reference are refused before any work is done
    if (mode in ['full', 'demix', 'aggregate', 'update_db']) {
        def unsupported = unsupportedPathogenReason(params.pathogen)
        if (unsupported) {
            error("--pathogen ${params.pathogen} is not supported: ${unsupported} See docs/other_pathogens.md.")
        }
        if (mode in ['full', 'demix'] && !params.panel) {
            if (!params.fasta) {
                error("--fasta is required with --pathogen ${params.pathogen}: wasteline bundles a reference only for the pathogens that have a file in conf/pathogens/. Pass the reference.fasta that the Freyja barcodes of ${params.pathogen} were built on (the pathogen's folder in https://github.com/andersen-lab/Freyja-barcodes): barcodes are tied to their reference and a different one gives silently wrong results. --gff (annotation) is optional. See docs/other_pathogens.md.")
            }
            validateReferenceContigs()
        }
    }
    if (params.whole_genome && !(mode in ['full', 'demix'])) {
        log.warn("--whole_genome selects the reference of the reads: it is only used by --mode full and demix and is ignored in --mode ${mode}")
    }
    if (params.barcodes == 'container' || params.compare_barcodes == 'container') {
        if (params.pathogen != params.freyja_container_pathogen) {
            error("--barcodes container is only available for the default pathogen (the only one the Freyja container bundles barcodes for): the pinned Freyja container bundles no barcodes for ${params.pathogen}. Use --barcodes latest, or a dataset directory or barcode file made by `freyja update --pathogen ${params.pathogen}`.")
        }
    }

    // Run defaults: samples without a protocol / primer_scheme cell use them. A samplesheet can override them per sample (PLAN D20).
    def settings = protocolSettings(params.protocol)
    if ((params.protocol == 'capture' || params.skip_primer_trimming) && (params.primer_scheme || params.primer_bed) && mode in ['full', 'demix']) {
        log.warn("--primer_scheme/--primer_bed are ignored for samples that are not primer-trimmed (" + (params.protocol == 'capture' ? "--protocol capture" : "--skip_primer_trimming") + "): there is no amplicon-level analysis for them")
    }
    if (params.primer_scheme_dir) {
        requireDirectory('--primer_scheme_dir', params.primer_scheme_dir)
        if (mode in ['full', 'demix']) {
            if (!userPrimerSchemes()) {
                log.warn("--primer_scheme_dir '${params.primer_scheme_dir}' contains no <name>.bed files: nothing can be selected from it")
            }
            def clashes = userPrimerSchemes().intersect(bundledPrimerSchemes())
            if (clashes) {
                log.warn("--primer_scheme_dir: ${clashes.join(', ')} ${clashes.size() == 1 ? 'is' : 'are'} also bundled for --pathogen ${params.pathogen}; your file${clashes.size() == 1 ? '' : 's'} from '${params.primer_scheme_dir}' will be used")
            }
        }
        else {
            log.warn("--primer_scheme_dir is not used by --mode ${mode} and is ignored")
        }
    }
    if (settings.primer_trimming && !params.panel && (mode == 'full' || (mode == 'demix' && (params.primer_scheme || params.primer_bed))) && primerBed()) {
        // The default BED must be on the same contig(s) as the reference or ivar trim silently trims nothing
        // (BEDs named in the samplesheet are checked when it is read)
        // (PLAN D24: lines of contigs that are not in the reference are dropped; only a BED with no line on it is refused)
        // (PLAN D26: no BED here when the pathogen has no default scheme and none was given: libraryMeta() then stops
        // for each sample that names none, so a samplesheet that gives every sample its own scheme still works)
        def bed_report = reconcileBed(primerBed())
        def fasta_contigs = fastaContigs(file(params.fasta, checkIfExists: true))
        def unknown = bed_report.no_match ? bedContigs(primerBed()) - fasta_contigs : []
        if (unknown && !params.primer_bed && !params.primer_scheme) {
            error("The bundled default primer scheme '${defaultPrimerScheme()}' is for a different reference (its contig(s) ${unknown} are not in the reference FASTA '${params.fasta}', which has ${fasta_contigs}): the bundled primer schemes belong to one pathogen only. Pass --primer_bed <BED whose first column is ${fasta_contigs.join(' / ')}> (a primer scheme for ${params.pathogen}), or --protocol capture / --skip_primer_trimming for reads that are not primer-trimmed. See docs/other_pathogens.md.")
        }
        if (unknown) {
            error("Primer BED '${primerBed()}' refers to contig(s) ${unknown} that are not in the reference FASTA '${params.fasta}' (which has ${fasta_contigs}). Use a matching --primer_bed/--primer_scheme or --fasta.")
        }
        if (!params.input) {
            warnDroppedBedLines([ bed_report ]) // with a samplesheet, checkSampleBeds warns once for all BEDs
        }
    }
    if (params.lineage_groups && !params.lineage_summary && mode in ['full', 'demix', 'aggregate']) {
        log.warn("--lineage_groups is ignored with --lineage_summary false")
    }
    if (mode in ['report', 'update_db'] && (params.lineage_groups || params.site_aggregation)) {
        log.warn("--lineage_groups/--site_aggregation are not used by --mode ${mode} and are ignored")
    }
    if (mode in ['report', 'update_db', 'aggregate'] && (params.primer_scheme || params.primer_bed)) {
        log.warn("--primer_scheme/--primer_bed are not used by --mode ${mode} and are ignored")
    }
    if (params.previous_results) {
        splitPaths(params.previous_results).each { dir -> requireDirectory('--previous_results', dir) }
    }
    if (params.demix_dir) {
        requireDirectory('--demix_dir', params.demix_dir)
    }
    if (mode == 'report') {
        if (!params.results_dir) {
            error("No --results_dir given. ${hint}")
        }
        def dir = requireDirectory('--results_dir', params.results_dir)
        if (params.outdir && file(params.outdir).toAbsolutePath().normalize() == dir.toAbsolutePath().normalize()) {
            error("--outdir must differ from --results_dir ('${params.results_dir}'): report mode never writes into the results it reads")
        }
        if (!file("${dir}/pipeline_info/run_info.json").exists()) {
            error("--results_dir '${params.results_dir}' has no pipeline_info/run_info.json: report mode needs the provenance record of the run (it never invents one)")
        }
    }
    else if (params.results_dir) {
        error("--results_dir is only used by --mode report (this is --mode ${mode}). ${hint}")
    }
    // Resolve the Freyja dataset specification early so a typo fails before any work is done
    // (a panel run resolves it per target in validatePanelParameters)
    if (mode in ['full', 'demix', 'aggregate'] && !params.panel) {
        resolveBarcodes(params.barcodes, dbDir(), '--barcodes')
    }
    if (params.compare_barcodes) {
        if (mode in ['full', 'demix']) {
            if (!params.panel) {
                resolveBarcodes(params.compare_barcodes, dbDir(), '--compare_barcodes') // a panel's targets are checked in validatePanelParameters
            }
        }
        else {
            log.warn("--compare_barcodes is only used by --mode full and demix and is ignored in --mode ${mode}")
        }
    }
    if (params.freyja_boot && !(mode in ['full', 'demix'])) {
        log.warn("--freyja_boot needs variants and depths: it is only used by --mode full and demix and is ignored in --mode ${mode}")
    }
    if (params.annotate_variants) {
        if (!(mode in ['full', 'demix'])) {
            log.warn("--annotate_variants is only used by --mode full and demix and is ignored in --mode ${mode}")
        }
        else if (params.skip_consensus) {
            error("--annotate_variants annotates the variants the consensus is built from: it cannot be combined with --skip_consensus")
        }
        else if (!params.gff && !params.panel) {
            error("--annotate_variants needs the reference annotation: --pathogen ${params.pathogen} has no bundled GFF3, pass --gff")
        }
    }
}

//
// Directory holding pipeline_info/run_info.json for a path inside a previous run (looks up at most
// 6 levels), or null
//
def findRunRoot(path) {
    return findRunRootFrom(new File(path.toString()).absoluteFile, 6)
}

def findRunRootFrom(dir, levels) {
    if (dir == null || levels == 0) {
        return null
    }
    if (new File(dir, 'pipeline_info/run_info.json').isFile()) {
        return dir
    }
    return findRunRootFrom(dir.parentFile, levels - 1)
}

//
// Provenance of the results a demix/aggregate run builds on: [ path, run_info_sha256, barcode_version ].
// run_info_sha256 is null and barcode_version 'unknown' when there is no run_info.json (legacy
// Snakemake results record neither).
//
def previousRunRecord(dir) {
    def run_info = new File(dir.toString(), 'pipeline_info/run_info.json')
    def record = [ path: dir.toString(), run_info_sha256: null, barcode_version: 'unknown' ]
    if (run_info.isFile()) {
        def digest = java.security.MessageDigest.getInstance('SHA-256').digest(run_info.bytes)
        record.run_info_sha256 = digest.collect { b -> String.format('%02x', b) }.join('')
        try {
            def info = new groovy.json.JsonSlurper().parse(run_info)
            // an aggregate-mode run_info describes the dataset used for plotting only, not a demix
            if (info.mode != 'aggregate') {
                record.barcode_version = info.freyja_data?.barcode_version ?: 'unknown'
            }
        }
        catch (Exception _e) {
            record.barcode_version = 'unknown'
        }
    }
    return record
}

//
// previous_runs list for run_info.json (mode demix/aggregate). --previous_results lists its
// directories; --input and --demix_dir list the previous runs that can be discovered from the
// locations of the files they point to.
//
def previousRuns(mode) {
    if (!(mode in ['demix', 'aggregate'])) {
        return []
    }
    def dirs = []
    if (params.previous_results) {
        dirs = splitPaths(params.previous_results).collect { dir -> file(dir).toAbsolutePath().toString() }
    }
    else {
        def locations = []
        if (params.demix_dir) {
            locations = [ file(params.demix_dir).toAbsolutePath().toString() ]
        }
        else if (params.input) {
            try {
                def lines = file(params.input).readLines().findAll { line -> line.trim() }
                def header = lines[0].split(',', -1).collect { column -> column.trim() }
                def column = header.indexOf(mode == 'demix' ? 'variants' : 'demix')
                locations = lines.drop(1)
                    .collect { line -> line.split(',', -1) }
                    .findAll { cells -> column >= 0 && column < cells.size() && cells[column].trim() }
                    .collect { cells -> file(cells[column].trim()).toAbsolutePath().toString() }
            }
            catch (Exception _e) {
                locations = []
            }
        }
        dirs = locations.collect { location -> findRunRoot(location)?.toString() }.findAll { root -> root }.unique().sort()
    }
    return dirs.collect { dir -> previousRunRecord(dir) }
}

//
// The bundled scheme used when neither --primer_scheme nor --primer_bed is given: the pathogen's default
// (conf/pathogens/*.config primer_scheme_default; SARS-CoV-2: qiagenv1, PLAN D14), null for pathogens without one (PLAN D26)
//
def defaultPrimerScheme(target = targetFromParams()) {
    return target.primer_scheme_default ?: null
}

//
// Primer BED of the run default: --primer_bed, else --primer_scheme (a scheme of --primer_scheme_dir or a bundled one),
// else the pathogen's default scheme; null when none of them is set (samples may name their own)
//
def primerBed(target = targetFromParams()) {
    if (params.primer_bed) {
        return file(params.primer_bed, checkIfExists: true)
    }
    def name = params.primer_scheme ?: defaultPrimerScheme(target)
    if (!name) {
        return null
    }
    def chosen = chooseScheme(name.toString(), "--primer_scheme '${name}'", target)
    return chosen.bed ? file(chosen.bed, checkIfExists: true) : file("${target.primer_schemes_bundled}/${name}.bed", checkIfExists: true)
}

//
// Local Freyja database root (writable): --db_dir or ~/.wasteline/db
//
def dbDir() {
    return params.db_dir ?: "${System.getProperty('user.home')}/.wasteline/db"
}

//
// Do the bundled primer schemes belong to the run's reference (is their contig in the FASTA)? True for the pathogen
// they were made for. Where they do not, the primer-scheme check (which scores a sample against every bundled scheme)
// has nothing to compare with and is skipped; a primer BED of the user is still used for trimming and amplicon QC.
//
def bundledSchemesApply(target = targetFromParams()) {
    if (!target.fasta || !defaultPrimerScheme(target)) {
        return false
    }
    return !(bedContigs(file("${target.primer_schemes_bundled}/${defaultPrimerScheme(target)}.bed", checkIfExists: true)) - fastaContigs(file(target.fasta, checkIfExists: true)))
}

def bedContigs(bed) {
    return bed.readLines()
        .findAll { line -> line.trim() && !line.startsWith('#') && !line.startsWith('track') }
        .collect { line -> line.split('\t')[0] }
        .unique()
}

def fastaContigs(fasta) {
    def lines = fasta.name.endsWith('.gz')
        ? new java.util.zip.GZIPInputStream(fasta.newInputStream()).newReader().readLines()
        : fasta.readLines()
    return lines
        .findAll { line -> line.startsWith('>') }
        .collect { line -> line.substring(1).trim().split(/\s+/)[0] }
}

//
// [ [ name: ..., length: ... ] ... ] of the FASTA records in file order (name = first word of the header), gzipped or not
//
def fastaRecords(fasta) {
    def lines = fasta.name.endsWith('.gz')
        ? new java.util.zip.GZIPInputStream(fasta.newInputStream()).newReader().readLines()
        : fasta.readLines()
    def records = []
    lines.each { line ->
        if (line.startsWith('>')) {
            records << [ name: line.substring(1).trim().split(/\s+/)[0], length: 0 ]
        }
        else if (records) {
            records[-1].length += line.trim().length()
        }
    }
    return records
}

//
// Does the reference FASTA have more than one sequence (PLAN D25)? Every code path for several contigs is guarded by
// this, so single-contig runs are untouched. Modes without a reference have no FASTA to look at. A run only gets past
// start-up with several sequences when --freyja_contig names one of them (validateReferenceContigs).
//
def multiContig(target = targetFromParams()) {
    return target.fasta && fastaContigs(file(target.fasta)).size() > 1
}

//
// The checks on the reference FASTA of a run that maps reads or demixes variants (modes full and demix): one sequence is
// the normal case; several are allowed only with --freyja_contig naming one of them, and the contig names (the first word of
// each header) must be unique (PLAN D25).
//
def validateReferenceContigs() {
    def contigs = fastaContigs(file(params.fasta, checkIfExists: true))
    def shown = "${contigs.take(10).join(', ')}${contigs.size() > 10 ? ', ...' : ''}"
    if (params.whole_genome && !params.freyja_contig) {
        error("--whole_genome is not available for --pathogen ${params.pathogen}: " + (params.pathogen == 'FLU-B-VIC'
            ? "there is no whole-genome B/Victoria reference (no consecutive eight-segment set exists for the strain), so FLU-B-VIC stays HA-only. Run without --whole_genome, or give your own multi-sequence reference with --fasta <multi.fa> --freyja_contig <name of the HA sequence>."
            : "whole-genome references are bundled for the influenza pathogens H1N1, H3N2, H5Nx and H5Nx-cattle only. For another multi-sequence reference give --fasta <multi.fa> --freyja_contig <name of the sequence Freyja demixes>.") + " See docs/other_pathogens.md.")
    }
    if (contigs.isEmpty()) {
        error("The reference FASTA '${params.fasta}' has no sequence (no '>' header line)")
    }
    def duplicates = contigs.findAll { name -> contigs.count(name) > 1 }.unique()
    if (duplicates) {
        error("The reference FASTA '${params.fasta}' has several sequences whose header starts with the same word (${duplicates.join(', ')}): wasteline names contigs by the first word of the FASTA header, so these must be unique.")
    }
    if (contigs.size() == 1) {
        if (params.freyja_contig) {
            error("--freyja_contig '${params.freyja_contig}' is only for a reference FASTA with several sequences, but '${params.fasta}' has one (${contigs[0]})." + (params.whole_genome ? " --whole_genome selects the bundled whole-genome reference; a --fasta of your own replaces it." : "") + " Drop --freyja_contig for a single-sequence reference.")
        }
        return
    }
    if (!params.freyja_contig) {
        error("The reference FASTA '${params.fasta}' has ${contigs.size()} sequences (${shown}): a multi-sequence reference needs --freyja_contig <name>, the sequence Freyja demixes (the one the barcodes were built on). The other sequences get depth QC, variant calls and a consensus genome, but no lineage demixing. For influenza, --whole_genome selects the bundled eight-segment reference of --pathogen H1N1, H3N2, H5Nx or H5Nx-cattle. See docs/other_pathogens.md.")
    }
    if (!(params.freyja_contig in contigs)) {
        error("--freyja_contig '${params.freyja_contig}' is not one of the sequences of the reference FASTA '${params.fasta}' (${shown}): give the first word of a FASTA header.")
    }
    if (params.gff) {
        def gff = file(params.gff, checkIfExists: true)
        def lines = gff.name.endsWith('.gz')
            ? new java.util.zip.GZIPInputStream(gff.newInputStream()).newReader().readLines()
            : gff.readLines()
        def seqids = lines.findAll { line -> line.trim() && !line.startsWith('#') }.collect { line -> line.split('\t')[0] }.unique()
        def unmatched = seqids - contigs
        if (unmatched) {
            log.warn("The annotation '${gff.name}' has sequence names that are not in the reference FASTA (${unmatched.take(10).join(', ')}${unmatched.size() > 10 ? ', ...' : ''}): its features on them are ignored (no amino-acid consequences there). The FASTA sequences are ${shown}.")
        }
    }
}

//
// Where the datasets of the run's pathogen live in the local database. The pathogen whose data the pinned Freyja
// container bundles (conf/pathogens/*.config: freyja_container_pathogen) keeps <db_dir>/freyja; every other pathogen
// gets <db_dir>/freyja_<pathogen>, so a date such as `--barcodes 2025-05-01` can never pick another pathogen's dataset.
//
def datasetRoot(db_dir, target = targetFromParams()) {
    def pathogen = target.freyja_pathogen
    return pathogen == target.freyja_container_pathogen ? "${db_dir}/freyja" : "${db_dir}/freyja_${pathogen}"
}

//
// Dataset directory names already stored in the local database
//
def availableDatasets(db_dir, target = targetFromParams()) {
    def root = file(datasetRoot(db_dir, target))
    return root.isDirectory() ? root.listFiles().findAll { entry -> entry.isDirectory() }.collect { entry -> entry.name }.sort() : []
}

//
// Resolve --barcodes (PLAN 5.6) into [ kind: 'fetch'|'path', source: latest|container|pinned|custom,
// path: <existing path, kind == 'path' only> ]. Never falls back silently to another dataset. A pinned version is looked up in the
// dataset root of `target` (default: the run's target).
//
def resolveBarcodes(barcodes, db_dir, option, target = targetFromParams()) {
    def spec = barcodes.toString().trim()
    if (spec == 'latest' || spec == 'container') {
        return [ kind: 'fetch', source: spec ]
    }
    def stamp = (spec =~ /^(\d{2})_(\d{2})_(\d{4})-(\d{2})-(\d{2})$/)
    def iso_form = (spec =~ /^(\d{4}-\d{2}-\d{2}T\d{2})[-:](\d{2})$/)
    if (stamp.matches() || iso_form.matches()) {
        def iso = stamp.matches()
            ? "${stamp[0][3]}-${stamp[0][1]}-${stamp[0][2]}T${stamp[0][4]}-${stamp[0][5]}"
            : "${iso_form[0][1]}-${iso_form[0][2]}"
        def dir = file("${datasetRoot(db_dir, target)}/${iso}")
        if (!dir.isDirectory()) {
            error("Freyja dataset '${spec}' (directory ${iso}) not found in ${datasetRoot(db_dir, target)}. Available: ${availableDatasets(db_dir, target) ?: 'none'}. Download datasets with ${option} latest (stored there by date and time), or point ${option} at a dataset directory or .csv/.feather file.")
        }
        return [ kind: 'path', source: 'pinned', path: dir ]
    }
    if (spec ==~ /^\d{4}-\d{2}-\d{2}$/) {
        def matches = availableDatasets(db_dir, target).findAll { name -> name.startsWith("${spec}T") }
        if (!matches) {
            error("No Freyja dataset from ${spec} found in ${datasetRoot(db_dir, target)}. Available: ${availableDatasets(db_dir, target) ?: 'none'}. Download one with ${option} latest, or point ${option} at a dataset directory or .csv/.feather file.")
        }
        return [ kind: 'path', source: 'pinned', path: file("${datasetRoot(db_dir, target)}/${matches.last()}") ]
    }
    def path = file(spec)
    if (!path.exists()) {
        error("${option} '${spec}' is not 'latest', 'container', a barcode version (YYYY-MM-DD, MM_DD_YYYY-HH-MM or YYYY-MM-DDTHH-MM) or an existing dataset directory / .csv / .feather file.")
    }
    if (path.isFile() && !(path.name.endsWith('.csv') || path.name.endsWith('.feather'))) {
        error("${option} file '${spec}' must end in .csv or .feather")
    }
    return [ kind: 'path', source: 'custom', path: path ]
}

//
// Can the mutation-level analyses (dropout impact, unexplained mutations) read the Freyja dataset in effect?
// They need the lineage -> mutation lists: a dataset directory with a non-empty lineage_mutations.json, or a
// barcodes .csv/.json. A .feather-only dataset (e.g. --barcodes custom.feather) has none we can read.
// Returns null when supported, otherwise the reason.
//
def mutationTablesUnsupportedReason(barcodes = params.barcodes, target = targetFromParams()) {
    def resolved = resolveBarcodes(barcodes, dbDir(), '--barcodes', target)
    if (resolved.kind == 'fetch') {
        return null // 'latest' / 'container': a Freyja dataset directory, which has lineage_mutations.json
    }
    def path = resolved.path
    if (path.isFile()) {
        return (path.name.endsWith('.csv') || path.name.endsWith('.json')) ? null : "${path.name} is a .feather barcodes file, which carries no mutation lists we can read"
    }
    def json = file("${path}/lineage_mutations.json")
    def has_csv = path.listFiles()?.any { f -> f.name.endsWith('.csv') && f.name.contains('barcode') }
    return (json.exists() && json.size() > 0) || has_csv ? null : "${path} has neither lineage_mutations.json nor a barcodes .csv"
}

//
// Can the lineage-group summary read the lineage hierarchy (lineages.yml) of the Freyja dataset in effect for `target`?
// A dataset directory has it (downloaded, pinned or bundled in the container); a bare barcodes .csv/.feather
// file (custom --barcodes, or the barcodes cell of a panel row) does not. Returns null when supported, otherwise the reason.
//
def lineageHierarchyUnsupportedReason(target = targetFromParams()) {
    def resolved = resolveBarcodes(target.barcodes ?: params.barcodes, dbDir(), '--barcodes', target)
    if (resolved.kind == 'fetch') {
        return null
    }
    def path = resolved.path
    if (path.isFile()) {
        return "${path.name} is a barcodes file without a lineage hierarchy (lineages.yml)"
    }
    // Freyja names the hierarchy of the pathogens other than the default <PATHOGEN>_lineages.yml
    return (file("${path}/lineages.yml").exists() || path.listFiles()?.any { f -> f.name.endsWith('_lineages.yml') }) ? null : "${path} has no lineages.yml"
}

//
// Does a dataset DOWNLOADED during the run (--barcodes latest) possibly lack the lineage hierarchy? Freyja ships a
// hierarchy for the default pathogen only, and for a few others (not RSV, dengue, ...): then the workflow looks at what
// was downloaded (lineages_yml_sha256 of freyja_data.json) instead of assuming (PLAN D22).
//
def hierarchyMayBeMissing(target = targetFromParams()) {
    return resolveBarcodes(target.barcodes ?: params.barcodes, dbDir(), '--barcodes', target).kind == 'fetch' && target.freyja_pathogen != target.freyja_container_pathogen
}

//
// The Freyja plot config (--freyja_plot_config; colours and lineage groups for freyja plot, dash and relgrowthrate),
// or [] (an optional process input left empty) when none is given
//
def plotConfigFile() {
    return params.freyja_plot_config ? file(params.freyja_plot_config, checkIfExists: true) : []
}

//
// The lineage groups file in effect for `target`: --lineage_groups (the same file for every target of a panel), else the bundled
// assets/lineage_groups.csv of the target (conf/pathogens/*.config or the registry: lineage_groups_bundled). null = none: the
// lineage-group summary of that target is then off.
//
def lineageGroupsFile(target = targetFromParams()) {
    if (params.lineage_groups) {
        return file(params.lineage_groups, checkIfExists: true)
    }
    return target.lineage_groups_bundled ? file(target.lineage_groups_bundled, checkIfExists: true) : null
}

//
// What run_info.json records about the lineage groups file: its name ('bundled default' for the bundled one),
// where it came from and its sha256
//
def lineageGroupsInfo(target = targetFromParams()) {
    def groups = lineageGroupsFile(target)
    if (!groups) {
        return [ lineage_groups_file: null, lineage_groups_source: null, lineage_groups_sha256: null ]
    }
    def digest = java.security.MessageDigest.getInstance('SHA-256').digest(groups.bytes)
    def sha = digest.collect { b -> String.format('%02x', b & 0xff) }.join('')
    return [
        lineage_groups_file: params.lineage_groups ? groups.name : 'bundled default',
        lineage_groups_source: params.lineage_groups ? 'custom' : 'bundled',
        lineage_groups_sha256: sha,
    ]
}
