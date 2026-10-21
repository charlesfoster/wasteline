//
// Resolve --barcodes into ONE Freyja dataset shared by every sample of the run (PLAN 5.6)
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
    barcodes // string: --barcodes
    pathogen // string: Freyja --pathogen name
    option   // string: the option the dataset came from (--barcodes or --compare_barcodes), for error messages

    main:
    def spec = resolveBarcodes(barcodes, dbDir(), option)

    def ch_dataset = channel.empty()
    def dataset_path = ''
    if (spec.kind == 'fetch') {
        FREYJA_UPDATE(channel.value(spec.source))
        ch_dataset = FREYJA_UPDATE.out.dataset
        dataset_path = spec.source == 'latest' ? datasetRoot(dbDir()) : 'bundled in the Freyja container'
    }
    else {
        ch_dataset = channel.value(spec.path)
        dataset_path = spec.path.toString()
    }

    FREYJA_DATASET_INFO(ch_dataset, spec.source, dataset_path, pathogen)

    emit:
    dataset         = ch_dataset                               // value channel: dataset directory or barcode file
    freyja_data     = FREYJA_DATASET_INFO.out.freyja_data      // value channel: full provenance block
    barcode_version = FREYJA_DATASET_INFO.out.barcode_version  // value channel: small version record
}
