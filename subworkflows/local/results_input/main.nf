//
// Samplesheet rows for modes demix/aggregate when no --input samplesheet is given
//
// --previous_results <dir>[,<dir>...]  results directories of wasteline (samples/<s>/freyja/) or of
//                                      the legacy Snakemake pipeline (<s>/freyja/). They are only read.
//                                      bin/make_samplesheet.py --from-results (run in DISCOVER_RESULTS)
//                                      is the single source of truth for what is in them.
// --demix_dir <dir>                    (aggregate only) every *.demix file below <dir>
//
// Emits the same rows as a samplesheet given with --input (see PIPELINE_INITIALISATION) plus the
// samplesheet that was effectively used, which is checksummed in run_info.json.
//

include { DISCOVER_RESULTS } from '../../../modules/local/discover_results/main'
include { sampleMeta       } from '../utils_nfcore_wasteline_pipeline'
include { splitPaths       } from '../utils_nfcore_wasteline_pipeline'

workflow RESULTS_INPUT {

    take:
    mode             // string: demix | aggregate
    previous_results // string: comma-separated results directories, or null
    demix_dir        // string: directory of *.demix files, or null

    main:
    def ch_sheet = channel.empty()
    def ch_log = channel.empty()
    def ch_skipped = channel.empty()
    def ch_rows = channel.empty()

    if (previous_results) {
        def dirs = splitPaths(previous_results)
        DISCOVER_RESULTS(
            channel.fromList(dirs.withIndex().collect { dir, index -> [ index + 1, file(dir) ] }),
            mode,
        )
        // one header, rows of all results directories
        ch_sheet = DISCOVER_RESULTS.out.samplesheet
            .collectFile(name: 'discovered_samplesheet.csv', keepHeader: true, skip: 1, sort: true)
        ch_log = DISCOVER_RESULTS.out.log
        // samples found but unusable (missing files, zero-byte legacy .demix): one 'sample<TAB>reason' line each
        ch_skipped = DISCOVER_RESULTS.out.skipped
            .splitCsv(header: true, sep: '\t')
            .map { row -> [ row.sample, row.reason, row.sample_type ?: 'sample', row.collection_date ?: '', row.site ?: '' ] }
        ch_rows = ch_sheet
            .splitCsv(header: true)
            .map { row -> resultsRow(mode, row) }
    }
    else {
        // --demix_dir: sample ID = file name without .demix. sample_type is inferred like
        // make_samplesheet.py does (NEG_/NC_ in the ID means a negative control).
        def found = demixDirFiles(demix_dir)
        // zero-byte .demix files are what the legacy pipeline left behind when `freyja demix` failed
        def empty = file("${demix_dir}/**.demix").findAll { path -> path.size() == 0 && !path.name.endsWith('.compare.demix') }
            .collect { path ->
                def id = path.name.substring(0, path.name.length() - '.demix'.length())
                [ id, 'empty demix file (failed in the old pipeline)', (id =~ /NEG_|NC_/).find() ? 'negative' : 'sample', '', '' ]
            }
        ch_skipped = channel.fromList(empty)
        ch_sheet = channel
            .fromList(found.collect { entry -> "${entry.sample},${entry.file},${entry.sample_type},," })
            .collectFile(name: 'demix_dir_samplesheet.csv', seed: 'sample,demix,sample_type,collection_date,site', sort: true, newLine: true)
        ch_rows = channel.fromList(found.collect { entry ->
            [ sampleMeta([ id: entry.sample, sample_type: entry.sample_type ], mode), entry.file ]
        })
    }

    // The same sample twice (e.g. two --previous_results directories that both contain it) would
    // silently double count: refuse.
    ch_rows
        .map { row -> row[0].id }
        .toList()
        .subscribe { ids ->
            def duplicates = ids.countBy { id -> id }.findAll { _id, count -> count > 1 }.keySet()
            if (duplicates) {
                error("Sample ID(s) found more than once across the inputs: ${duplicates.sort().join(', ')}. Sample IDs must be unique; use a single results directory or remove the duplicates.")
            }
        }

    emit:
    rows        = ch_rows    // [ meta, ... ] as for a samplesheet of the mode
    samplesheet = ch_sheet   // value channel (single file): the effective samplesheet
    skipped     = ch_skipped // [ sample, reason, sample_type, collection_date, site ] of samples that were found but not usable
    notes       = ch_log   // discovery notes (one file per results directory)
}

//
// Row of a discovered samplesheet (strings) -> channel item for the mode
//
def resultsRow(mode, row) {
    // protocol / primer_scheme: carried over from the samplesheet of the previous run when it had the columns
    def meta = sampleMeta([
        id: row.sample,
        sample_type: row.sample_type,
        collection_date: row.collection_date,
        site: row.site,
        protocol: row.protocol,
        primer_scheme: row.primer_scheme,
    ], mode)
    if (mode == 'demix') {
        return [ meta, file(row.variants), file(row.depths), row.fastp_json ? file(row.fastp_json) : [] ]
    }
    return [ meta, file(row.demix) ]
}

//
// *.demix files below a directory: [ sample, file, sample_type ] sorted by sample
//
def demixDirFiles(dir) {
    // <sample>.compare.demix are the --compare_barcodes results of a wasteline run, not samples
    def files = file("${dir}/**.demix").findAll { path -> !path.name.endsWith('.compare.demix') }.sort { a, b -> a.name <=> b.name }
    if (!files) {
        error("No *.demix files found below --demix_dir '${dir}'")
    }
    if (files.every { path -> path.size() == 0 }) {
        error("Every *.demix file below --demix_dir '${dir}' is empty")
    }
    def entries = files.findAll { path -> path.size() > 0 }.collect { path ->
        def sample = path.name.substring(0, path.name.length() - '.demix'.length())
        if (!(sample ==~ /^[A-Za-z0-9._-]+$/)) {
            error("Cannot use '${path}' from --demix_dir: the sample ID '${sample}' must match ^[A-Za-z0-9._-]+\$ (IDs are never rewritten; rename the file)")
        }
        [ sample: sample, file: path.toString(), sample_type: (sample =~ /NEG_|NC_/).find() ? 'negative' : 'sample' ]
    }
    return entries
}
