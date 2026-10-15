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
    def after_text = ""

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
    def ch_samplesheet = channel.empty()
    if (input && mode == 'full') {
        ch_samplesheet = channel
            .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
            .map { meta, fastq_1, fastq_2 -> [ sampleMeta(meta), [ fastq_1, fastq_2 ] ] }
    }
    else if (input && mode == 'demix') {
        ch_samplesheet = channel
            .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input_demix.json"))
            .map { meta, variants, depths, fastp_json -> [ sampleMeta(meta), variants, depths, fastp_json ?: [] ] }
    }
    else if (input && mode == 'aggregate') {
        ch_samplesheet = channel
            .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input_aggregate.json"))
            .map { meta, demix -> [ sampleMeta(meta), demix ] }
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
    def run_info = file("${outdir}/pipeline_info/run_info.json")
    if (run_info.exists()) {
        def info = new groovy.json.JsonSlurper().parseText(run_info.text)
        info.containers = runContainers(outdir)
        info.run_completed = workflow.complete?.toString()
        info.duration = workflow.duration?.toString()
        info.success = workflow.success
        run_info.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
    }
}

//
// Sample meta map shared by every mode
//
def sampleMeta(meta) {
    return [
        id: meta.id,
        sample_type: meta.sample_type ?: 'sample',
        collection_date: meta.collection_date ?: null,
        site: meta.site ?: null,
    ]
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
        full: "--mode full needs a FASTQ samplesheet: --input samplesheet.csv (columns sample,fastq_1,fastq_2[,sample_type][,collection_date][,site]; make one with bin/make_samplesheet.py <fastq_dir>)",
        demix: "--mode demix needs existing Freyja variants and depths: --input samplesheet.csv (columns sample,variants,depths[,fastp_json][,sample_type][,collection_date][,site]) or --previous_results <results dir> (new or legacy Snakemake layout; never modified)",
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

    def settings = protocolSettings(mode)
    if ((params.protocol == 'capture' || params.skip_primer_trimming) && (params.primer_scheme || params.primer_bed) && mode in ['full', 'demix']) {
        log.warn("--primer_scheme/--primer_bed are ignored: " + (params.protocol == 'capture' ? "--protocol capture" : "--skip_primer_trimming") + " reads are not primer-trimmed and there is no amplicon-level analysis")
    }
    if (settings.amplicons) {
        // The BED must be on the same contig(s) as the reference or ivar trim silently trims nothing
        def bed_contigs = bedContigs(primerBed())
        def fasta_contigs = fastaContigs(file(params.fasta, checkIfExists: true))
        def unknown = bed_contigs - fasta_contigs
        if (unknown) {
            error("Primer BED '${primerBed()}' refers to contig(s) ${unknown} that are not in the reference FASTA '${params.fasta}' (which has ${fasta_contigs}). Use a matching --primer_bed/--primer_scheme or --fasta.")
        }
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
    if (mode in ['full', 'demix', 'aggregate']) {
        resolveBarcodes(params.barcodes, dbDir(), '--barcodes')
    }
    if (params.compare_barcodes) {
        if (mode in ['full', 'demix']) {
            resolveBarcodes(params.compare_barcodes, dbDir(), '--compare_barcodes')
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
// THE effective library settings, computed here and nowhere else (conf/modules.config spells out the two
// naming conditions it needs, because a config file cannot call a function):
//   protocol         amplicon | capture
//   primer_trimming  true only for amplicon reads that are not --skip_primer_trimming
//   duplicates       mark | none (--duplicates auto = mark for capture, none for amplicon)
//   amplicons        amplicon-level QC (primer BED needed): amplicon protocol with trimming; in demix mode
//                    only when a scheme/BED is given explicitly (the scheme of old results is unknown)
//
def protocolSettings(mode) {
    def capture = params.protocol == 'capture'
    def trimming = !capture && !params.skip_primer_trimming
    return [
        protocol: params.protocol,
        primer_trimming: trimming,
        duplicates: params.duplicates == 'auto' ? (capture ? 'mark' : 'none') : params.duplicates,
        amplicons: trimming && (mode == 'full' || (mode == 'demix' && (params.primer_scheme || params.primer_bed))),
    ]
}

//
// Primer BED: --primer_bed, else the bundled --primer_scheme
//
def primerBed() {
    if (params.primer_bed) {
        return file(params.primer_bed, checkIfExists: true)
    }
    return file("${projectDir}/assets/primer_schemes/${params.primer_scheme ?: 'midnight'}.bed", checkIfExists: true)
}

//
// ivar trim -x: explicit --ivar_offset wins; the bundled swift scheme implies 5; otherwise 0
//
def ivarOffset() {
    if (params.ivar_offset != null) {
        return params.ivar_offset as int
    }
    return (!params.primer_bed && params.primer_scheme == 'swift') ? 5 : 0
}

//
// Local Freyja database root (writable): --db_dir or ~/.wasteline/db
//
def dbDir() {
    return params.db_dir ?: "${System.getProperty('user.home')}/.wasteline/db"
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
// Dataset directory names already stored in the local database
//
def availableDatasets(db_dir) {
    def root = file("${db_dir}/freyja")
    return root.isDirectory() ? root.listFiles().findAll { entry -> entry.isDirectory() }.collect { entry -> entry.name }.sort() : []
}

//
// Resolve --barcodes (PLAN 5.6) into [ kind: 'fetch'|'path', source: latest|container|pinned|custom,
// path: <existing path, kind == 'path' only> ]. Never falls back silently to another dataset.
//
def resolveBarcodes(barcodes, db_dir, option) {
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
        def dir = file("${db_dir}/freyja/${iso}")
        if (!dir.isDirectory()) {
            error("Freyja dataset '${spec}' (directory ${iso}) not found in ${db_dir}/freyja. Available: ${availableDatasets(db_dir) ?: 'none'}. Download datasets with ${option} latest (stored there by date and time), or point ${option} at a dataset directory or .csv/.feather file.")
        }
        return [ kind: 'path', source: 'pinned', path: dir ]
    }
    if (spec ==~ /^\d{4}-\d{2}-\d{2}$/) {
        def matches = availableDatasets(db_dir).findAll { name -> name.startsWith("${spec}T") }
        if (!matches) {
            error("No Freyja dataset from ${spec} found in ${db_dir}/freyja. Available: ${availableDatasets(db_dir) ?: 'none'}. Download one with ${option} latest, or point ${option} at a dataset directory or .csv/.feather file.")
        }
        return [ kind: 'path', source: 'pinned', path: file("${db_dir}/freyja/${matches.last()}") ]
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
