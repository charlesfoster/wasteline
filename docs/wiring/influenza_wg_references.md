# Whole-genome influenza references (D25 phase 1, task T0)

Three multi-FASTA references (all eight segments) with a GFF3 each, for `--whole_genome` and for custom `--fasta` multi-record runs.
Built 2026-10-22 from NCBI nuccore with E-utilities (`efetch rettype=fasta|gb|gff3`); the downloads were treated as untrusted data and
parsed with `python -I`. Contract: `docs/wiring/multi_contig_targets.md`.

| Pathogen profile  | Files                                             | Strain                        | Segments | Total     |
| ----------------- | ------------------------------------------------- | ----------------------------- | -------- | --------- |
| H1N1              | `assets/references/CY121680.1_wg.fasta` + `.gff3` | A/California/07/2009 (H1N1)   | 8        | 13,417 bp |
| H3N2              | `assets/references/CY163680.1_wg.fasta` + `.gff3` | A/Wisconsin/67/2005 (H3N2)    | 8        | 13,367 bp |
| H5Nx, H5Nx-cattle | `assets/references/AF144305.1_wg.fasta` + `.gff3` | A/Goose/Guangdong/1/96 (H5N1) | 8        | 13,590 bp |

**B/Victoria is not built.** There is no consecutive eight-segment set for the strain (a composite of three submissions), so FLU-B-VIC stays
HA-only (`KX058884.1`) until a proper set is chosen.

## Record order and contig names

HA first, then a fixed order: **HA, PB2, PB1, PA, NP, NA, M, NS**. The contig name is the first header word (accession.version), as in the
HA-only files; the header text is NCBI's, unchanged. The HA contig name (`CY121680.1`, `CY163680.1`, `AF144305.1`) is the profile's
`freyja_contig`. Sequence lines are 70 columns, one trailing newline, no blank line.

| Segment | H1N1       | length | H3N2       | length | H5N1       | length |
| ------- | ---------- | ------ | ---------- | ------ | ---------- | ------ |
| HA      | CY121680.1 | 1752   | CY163680.1 | 1737   | AF144305.1 | 1760   |
| PB2     | CY121687.1 | 2316   | CY163687.1 | 2286   | AF144300.1 | 2341   |
| PB1     | CY121686.1 | 2317   | CY163686.1 | 2311   | AF144301.1 | 2341   |
| PA      | CY121685.1 | 2192   | CY163685.1 | 2195   | AF144302.1 | 2233   |
| NP      | CY121683.1 | 1541   | CY163683.1 | 1537   | AF144303.1 | 1565   |
| NA      | CY121682.1 | 1432   | CY163682.1 | 1435   | AF144304.1 | 1458   |
| M       | CY121681.1 | 1002   | CY163681.1 | 1001   | AF144306.1 | 1027   |
| NS      | CY121684.1 | 865    | CY163684.1 | 865    | AF144307.1 | 865    |

## Checks done on every record (all passed)

- the title names the strain and the segment (the AF144306 and AF144307 titles, truncated in esummary, read in full: "matrix proteins M1
  and M2 (M) gene" and "nonstructural proteins NS1 and NS2 (NS) gene", both "alternatively spliced products, complete cds");
- the first header word equals the accession.version asked for; no accession and no segment appears twice; sequences are A-Z only (in fact
  only ACGT, no ambiguity codes); the length equals the `##sequence-region` of NCBI's GFF3;
- the HA record is byte-identical to the bundled single-segment `<acc>.fasta` (`head -n <lines of the old file>` of the new file `cmp`s
  equal, all three), and the HA feature lines of the GFF3 are identical to the lines of the old `<acc>.gff3` (`diff` empty);
- `samtools faidx` reports 8 contigs with the lengths above and `bwa index` works (container
  `community.wave.seqera.io/library/bwa_htslib_samtools`);
- every CDS translates cleanly (length a multiple of 3, starts with M, ends with a stop, no internal stop); M2 is 97 aa and NEP 121 aa
  after joining the two exons.

## GFF3

Same convention as the HA-only files: NCBI's gene and CDS coordinates unchanged; a synthetic mRNA and an exon line per exon so that
`bcftools csq` has gene -> transcript -> exon/CDS (IDs `gene:<name>`, `transcript:<name>`, `CDS:<name>`, `biotype=protein_coding`).
seqids equal the FASTA contig names. One `##sequence-region` line per contig.

- Gene names are NCBI's gene symbols: HA, PB2, PB1, PA, NP, NA, M1, M2, NEP, NS1 (H3N2 also PB1-F2). NCBI gave the H5N1 M and NS records
  one gene (`M`, `NS`) with two products, so there M1/M2 and NS1/NEP are named from the product, the same names as in H1N1 and H3N2.
- **M2 and NEP (NS2) are spliced**: one mRNA, two exon lines (`M2.exon1`, `M2.exon2`) and two CDS lines with the same ID, phases 0 and 1
  for M2 and 0 and 0 for NEP. NCBI's GFF3 gives the joined coordinates, so nothing was inferred.
- **Skipped**: PA-X (ribosomal slippage; NCBI annotates it on the H1N1 and H3N2 PA records only); the H1N1 PB1-F2 (NCBI marks it as a
  non-functional `misc_feature`, no CDS); H5N1 has no PB1-F2 or PA-X annotation to begin with.
- H3N2 PB1-F2 (a CDS in another frame inside PB1) is kept as its own gene.

## `bcftools csq` verification

`quay.io/biocontainers/bcftools:1.23--h3a4d415_0` (the pin of `annotate_variants`), `bcftools csq -f <fasta> -g <gff3> --force -pm`, with
one synthetic missense SNV (GT 1/1) per gene in a region owned by that gene, plus, for the spliced genes, SNVs in exon 1, across the
exon junction and near the end of exon 2. All 49 SNVs gave the expected gene and amino-acid change (16 in H1N1, 17 in H3N2, 16 in H5N1),
for every segment, including M2 and NEP in both exons and across the junction (a codon split by the junction is reported as `splice_region`
for the exon-1 base and `missense` for the joined codon). The per-gene results are in the T0 report.

## Known limits

- **Overlapping genes and the pipeline's "first consequence only" rule.** `annotate_variants` reports only the first consequence of a
  variant. In the shared first exon of M1/M2 and NS1/NEP (and in PB1-F2 inside PB1) csq lists both genes; the table shows M1, NS1 or PB1
  (csq's first), so a variant that changes M2 or NEP in the shared region is attributed to M1/NS1. Variants outside the overlap are exact.
- `bcftools csq` prints "Ignoring GFF feature with unknown phase" for the exon lines (phase `.`), exactly as for the existing HA-only
  files; it is harmless.
- NCBI complete-cds records have no (or few) UTRs: the 5' and 3' ends of each segment are short, so reads from the conserved termini map
  poorly; the consensus there is reference-biased.
- These are single historical strains: the internal segments of recent viruses differ from them (H5N1 Goose/Guangdong/1/96 internal genes
  are 10-15 % divergent from 2.3.4.4b cattle viruses; PB2, PB1, PA, NP, M and NS of recent H1N1/H3N2 are also distant), so low mapping and
  reference-biased consensus are expected away from HA/NA. Only the HA contig is demixed.
- The `.gff3` header comments name the eight records; nothing parses them.
