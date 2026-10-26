//
// Resolve --barcodes into ONE Freyja dataset per target, shared by every sample of that target (PLAN 5.6; a run has one target
// unless it is a panel run)
//
// latest     -> freyja update (published to <db_dir>/freyja/<YYYY-MM-DDTHH-MM>/; <db_dir>/freyja_<pathogen>/ for other pathogens, PLAN D22)
// container  -> snapshot bundled in the pinned Freyja container (offline)
// stamp/date -> existing dataset in <db_dir>/freyja (error if absent)
// path       -> dataset directory or .csv/.feather file, used verbatim
// A run never falls back to a different dataset; the version and checksums are recorded.
//

include { FREYJA_UPDATE       } from '../../../modules/local/freyja_update/main'
include { FREYJA_DATASET_INFO } from '../../../modules/local/freyja_dataset_info/main'
include { resolveBarcodes     } from '../utils_nfcore_wasteline_pipeline'
include { dbDir               } from '../utils_nfcore_wasteline_pipeline'
include { datasetRoot         } from '../utils_nfcore_wasteline_pipeline'

workflow FREYJA_DATA {

    take:
    barcodes // string: --barcodes (or --compare_barcodes)
    targets  // list (plain Groovy) of target maps: one Freyja dataset is resolved per target, with that target's Freyja --pathogen
    option   // string: the option the dataset came from (--barcodes or --compare_barcodes), for error messages

    main:
    // The barcodes spec of a target: a registry row's own `barcodes` cell beats --barcodes (a panel run stays offline when its targets
    // ship their barcodes); --compare_barcodes is the user's explicit second dataset and is never overridden by a row.
    def specs = targets.collect { t ->
        def spec = resolveBarcodes(option == '--barcodes' ? (t.barcodes ?: barcodes) : barcodes, dbDir(), option, t)
        [
            meta: [ id: t.name, target: t ],
            spec: spec,
            dataset_path: spec.kind == 'fetch' ? (spec.source == 'latest' ? datasetRoot(dbDir(), t) : 'bundled in the Freyja container') : spec.path.toString(),
        ]
    }
    def fetch = specs.findAll { s -> s.spec.kind == 'fetch' }
    def fixed = specs.findAll { s -> s.spec.kind == 'path' }
    def fetch_by_name = fetch.collectEntries { s -> [ (s.meta.id): s ] }

    // [ meta, source ] of the targets whose dataset is downloaded or copied out of the container
    FREYJA_UPDATE(channel.fromList(fetch.collect { s -> [ s.meta, s.spec.source ] }))
    def ch_fetched = FREYJA_UPDATE.out.dataset.map { meta, dataset ->
        [ meta, dataset, fetch_by_name[meta.id].spec.source, fetch_by_name[meta.id].dataset_path ]
    }
    // the others use an existing dataset verbatim
    def ch_fixed = channel.fromList(fixed.collect { s -> [ s.meta, s.spec.path, s.spec.source, s.dataset_path ] })
    def ch_resolved = ch_fetched.mix(ch_fixed)

    FREYJA_DATASET_INFO(ch_resolved)

    // [ target name, item ]: join onto the samples of a target with combine(..., by: 0) on meta.target.name
    def ch_dataset_by_target = ch_resolved.map { meta, dataset, _source, _dataset_path -> [ meta.id, dataset ] }
    def ch_freyja_data_by_target = FREYJA_DATASET_INFO.out.freyja_data.map { meta, freyja_data -> [ meta.id, freyja_data ] }
    def ch_barcode_version_by_target = FREYJA_DATASET_INFO.out.barcode_version.map { meta, barcode_version -> [ meta.id, barcode_version ] }

    // The run-level consumers (aggregate, plot, report, provenance) take ONE dataset: they get value channels for a one-target
    // run and an empty channel otherwise (they are not available for several targets yet).
    def one_target = targets.size() == 1
    def ch_dataset = one_target ? ch_dataset_by_target.map { _name, dataset -> dataset }.first() : channel.empty()
    def ch_freyja_data = one_target ? ch_freyja_data_by_target.map { _name, freyja_data -> freyja_data }.first() : channel.empty()
    def ch_barcode_version = one_target ? ch_barcode_version_by_target.map { _name, barcode_version -> barcode_version }.first() : channel.empty()

    emit:
    dataset                   = ch_dataset                    // value channel (one target only): dataset directory or barcode file
    freyja_data               = ch_freyja_data                // value channel (one target only): full provenance block
    barcode_version           = ch_barcode_version            // value channel (one target only): small version record
    dataset_by_target         = ch_dataset_by_target          // [ target name, dataset ] per target
    freyja_data_by_target     = ch_freyja_data_by_target      // [ target name, freyja_data.json ] per target
    barcode_version_by_target = ch_barcode_version_by_target  // [ target name, barcode_version.json ] per target
}
