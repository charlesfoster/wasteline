//
// Panel front end (--panel, several targets in one run): every sample's trimmed reads are mapped ONCE to the concatenated panel
// reference, each sample is typed (which targets are present), and the mapped reads of every (sample, target) pair that is to be
// analysed are cut out as a BAM of that target alone, in the target's own reference naming. The result feeds the per-sample
// chain of the single-target path (primer trimming, Freyja, QC, consensus), once per pair. Formats and the channel contract of
// the building blocks: docs/wiring/panel_typing_interface.md.
//
include { PANEL_PREPARE } from '../../../modules/local/panel_prepare/main'
include { PANEL_ASSIGN  } from '../../../modules/local/panel_assign/main'
include { PANEL_CALL    } from '../../../modules/local/panel_call/main'
include { PANEL_SPLIT   } from '../../../modules/local/panel_split/main'
include { BWA_INDEX     } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM       } from '../../../modules/nf-core/bwa/mem/main'
include { sampleMeta    } from '../utils_nfcore_wasteline_pipeline'

//
// panel_resolved.csv (target,fasta,marker_contig,family; absolute FASTA paths) for PANEL_PREPARE, written once under the work
// directory (named by a digest of its text, so identical panels share it and a changed panel gets a new file)
//
def panelResolvedFile(targets) {
    def lines = [ 'target,fasta,marker_contig,family' ] + targets.collect { t ->
        [ t.name, file(t.fasta, checkIfExists: true).toAbsolutePath().toString(), t.marker_contig ?: '', t.family ?: '' ].join(',')
    }
    def text = lines.join('\n') + '\n'
    def digest = java.security.MessageDigest.getInstance('SHA-256').digest(text.getBytes('UTF-8')).collect { b -> String.format('%02x', b) }.join('').substring(0, 12)
    def resolved = file("${workflow.workDir}/wasteline_panel/${digest}/panel_resolved.csv")
    if (!resolved.exists() || resolved.text != text) {
        // not rewritten when it is there: a new modification time would invalidate -resume
        resolved.parent.mkdirs()
        resolved.text = text
    }
    return resolved
}

//
// The distinct reference FASTAs of the panel, to be staged flat in one task directory: a FASTA shared by several targets
// (H5Nx and H5Nx-cattle) once; two different files with the same name cannot be staged together and are refused.
//
def panelFastas(targets) {
    def by_name = [:]   // file name -> real path
    def fastas = []
    targets.each { t ->
        def fasta = file(t.fasta, checkIfExists: true)
        def real = fasta.toRealPath().toString()
        def other = by_name[fasta.name]
        if (other && other != real) {
            error("Panel references '${other}' and '${real}' have the same file name '${fasta.name}': give them different names (they are staged together in one directory)")
        }
        if (!other) {
            by_name[fasta.name] = real
            fastas << fasta
        }
    }
    return fastas
}

//
// The barrier of the panel run. items = [ [ meta, <id>.targets.tsv ], ... ] of ALL samples. A non-negative sample is analysed for
// the targets it has `detected`; a negative control for every target detected in at least one non-negative sample (so that a
// control shows what it carries of the targets of the run), whatever its own call. Returns
//   header  the columns of the table
//   rows    all table rows (lists of cells), by sample then panel order; `analysed` rewritten for negative controls
//   pairs   [ [ meta, target name ], ... ] the (sample, target) pairs to analyse
//
def panelDecision(items) {
    def sorted = items.sort(false) { a, b -> a[0].id <=> b[0].id }
    def header = null
    def parsed = sorted.collect { meta, tsv ->
        def lines = tsv.readLines().findAll { line -> line.trim() }
        def cols = lines[0].split('\t', -1) as List
        if (header == null) {
            header = cols
        }
        else if (cols != header) {
            error("Sample '${meta.id}': ${tsv.name} has other columns than the other samples' targets tables")
        }
        [ meta, lines.drop(1).collect { line -> line.split('\t', -1) as List } ]
    }
    if (header == null) {
        return [ header: [], rows: [], pairs: [] ]
    }
    def i_target = header.indexOf('target')
    def i_status = header.indexOf('status')
    def i_analysed = header.indexOf('analysed')
    def detected_somewhere = parsed.findAll { meta, _rows -> meta.sample_type != 'negative' }
        .collectMany { _meta, rows -> rows.findAll { row -> row[i_status] == 'detected' }.collect { row -> row[i_target] } }
        .toSet()
    def rows_out = []
    def pairs = []
    parsed.each { meta, rows ->
        rows.each { row ->
            def analysed = meta.sample_type == 'negative' ? row[i_target] in detected_somewhere : row[i_status] == 'detected'
            row[i_analysed] = analysed ? 'true' : 'false'
            rows_out << row
            if (analysed) {
                pairs << [ meta, row[i_target] ]
            }
        }
    }
    return [ header: header, rows: rows_out, pairs: pairs ]
}

//
// The one log line of the typing: "panel: 5 samples x 2 targets: 6 analysed, 1 trace, 3 not detected"
//
def panelSummaryLine(decision) {
    def i_target = decision.header.indexOf('target')
    def i_status = decision.header.indexOf('status')
    def i_analysed = decision.header.indexOf('analysed')
    def n_samples = decision.rows.collect { row -> row[0] }.toSet().size()
    def n_targets = decision.rows.collect { row -> row[i_target] }.toSet().size()
    def trace = decision.rows.findAll { row -> row[i_status] == 'trace' && row[i_analysed] != 'true' }.size()
    def not_detected = decision.rows.findAll { row -> row[i_status] == 'not_detected' && row[i_analysed] != 'true' }.size()
    def unanalysed = decision.rows.findAll { row -> row[i_status] == 'detected' && row[i_analysed] != 'true' }.size()
    def line = "panel: ${n_samples} samples x ${n_targets} targets: ${decision.pairs.size()} analysed, ${trace} trace, ${not_detected} not detected"
    return unanalysed ? "${line}, ${unanalysed} detected in a negative control only (not analysed)" : line
}

workflow PANEL_FRONT {

    take:
    ch_reads // channel: [ meta, reads ] after fastp; meta.id and meta.sample_type as given, optionally collection_date, site, viral_load, protocol_cell (or protocol) and primer_scheme_cell (the raw samplesheet cells)
    panel_targets // list of target maps (targetsFromParams()), in panel order
    mode     // string: full

    main:
    if (params.aligner == 'minibwa') {
        error("--aligner minibwa is not supported with --panel: the panel reference is mapped with bwa mem. Use --aligner bwa (the default)")
    }

    //
    // Panel reference: contigs of all targets, named <target>|<contig>
    //
    PANEL_PREPARE(channel.value(panelResolvedFile(panel_targets)), channel.value(panelFastas(panel_targets)))
    def ch_panel_meta = PANEL_PREPARE.out.fasta.map { fasta -> [ [ id: 'panel' ], fasta ] }

    //
    // One competitive mapping of every sample to the whole panel (kept unfiltered: unmapped reads are ignored by the typing)
    //
    BWA_INDEX(ch_panel_meta)
    BWA_MEM(ch_reads, BWA_INDEX.out.index, ch_panel_meta, true)

    //
    // Typing: which targets does each sample carry? <id>.targets.tsv has a row per target, also for a sample without panel reads
    //
    PANEL_ASSIGN(BWA_MEM.out.bam, PANEL_PREPARE.out.contig_map, params.panel_mapq)
    PANEL_CALL(
        PANEL_ASSIGN.out.stats,
        PANEL_PREPARE.out.contig_map,
        params.panel_min_fragments,
        params.panel_min_breadth_pct,
        params.panel_family_share_pct,
    )

    //
    // Barrier: negative controls need what the other samples carry, so the decision waits for every sample's typing
    //
    def ch_decision = PANEL_CALL.out.targets.toList().map { items -> panelDecision(items) }

    ch_decision.subscribe { decision -> log.info(panelSummaryLine(decision)) }

    // the run's typing table: one header, every sample's rows
    def ch_table = ch_decision
        .filter { decision -> decision.header }
        .map { decision -> ([ decision.header.join('\t') ] + decision.rows.collect { row -> row.join('\t') }).join('\n') + '\n' }
        .collectFile(name: 'targets.tsv', storeDir: "${params.outdir}/summary")

    //
    // Split: the reads of each analysed (sample, target) pair, mapped to that target alone. The per-target meta carries the
    // target and the sample's library settings for it (sampleMeta -> libraryMeta); the raw samplesheet cells come from the meta.
    //
    def ch_by_id = BWA_MEM.out.bam.map { meta, bam -> [ meta.id, bam ] }
    def ch_pairs = ch_decision.flatMap { decision -> decision.pairs }.map { meta, target_name ->
        def target = panel_targets.find { t -> t.name == target_name }
        def sample = [
            id: meta.id,
            sample_type: meta.sample_type,
            collection_date: meta.collection_date,
            site: meta.site,
            viral_load: meta.viral_load,
            protocol: meta.protocol_cell ?: meta.protocol,
            primer_scheme: meta.primer_scheme_cell,
        ]
        [ meta.id, sampleMeta(sample, mode, target) ]
    }
    PANEL_SPLIT(
        ch_pairs.combine(ch_by_id, by: 0).map { _id, meta_t, bam -> [ meta_t, bam ] },
        PANEL_PREPARE.out.contig_map,
    )

    emit:
    mapped         = PANEL_SPLIT.out.counted                                              // [ meta_t, bam, bai, mapped_reads.txt ], one per analysed (sample, target)
    targets        = ch_table                                                             // summary/targets.tsv (one item)
    typing         = PANEL_CALL.out.targets                                               // [ meta, <id>.targets.tsv ] per sample
    contig_map     = PANEL_PREPARE.out.contig_map                                         // contig_map.tsv
    panel_fasta    = PANEL_PREPARE.out.fasta                                              // panel.fasta
    analysed_pairs = ch_decision.map { decision -> decision.pairs.collect { meta, name -> [ meta.id, name ] } } // one item: [ [ sample, target ], ... ]
}
