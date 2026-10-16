# charlesfoster/wasteline: Citations

If you use wasteline, please cite the tools that did the work in your analysis. Entries marked
`<!-- verify citation -->` are given as a repository URL or from memory and should be checked against the publisher's
page before you cite them in a paper.

## Pipeline lineage

- **wasteline** is a Nextflow port of the wastewater workflow of `covid-illumina-snakemake` (see
  [`docs/migration_from_snakemake.md`](docs/migration_from_snakemake.md)).

  > Foster C. covid-illumina-snakemake. <https://github.com/charlesfoster/covid-illumina-snakemake>

- Conversion of the ivar variants table to VCF (`bin/ivar_variants_to_vcf.py`) was adapted, via the original Snakemake
  pipeline, from the nf-core/viralrecon converter.

  > Patel H, et al. nf-core/viralrecon. Zenodo. doi: 10.5281/zenodo.3901628 <!-- verify citation -->

## Workflow framework

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

- [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/) (the framework, the pipeline template, the nf-core modules and nf-core/tools)

  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

- [nf-test](https://www.nf-test.com/) (pipeline and module testing)

  > Forer L, Schönherr S. Improving the reliability and quality of Nextflow pipelines with nf-test. bioRxiv. 2024. doi: 10.1101/2024.05.25.595877 <!-- verify citation -->

- [nf-schema](https://github.com/nextflow-io/nf-schema) (parameter and samplesheet validation)

  > <https://github.com/nextflow-io/nf-schema> <!-- verify citation -->

## Pipeline tools

- [fastp](https://github.com/OpenGene/fastp) (read trimming and QC)

  > Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018 Sep 1;34(17):i884-i890. doi: 10.1093/bioinformatics/bty560.

  > Chen S. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta. 2023;2(2):e107. doi: 10.1002/imt2.107 <!-- verify citation -->

- [BWA (`bwa mem`)](https://github.com/lh3/bwa) (read mapping)

  > Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997 [q-bio.GN]. 2013.

- [minibwa](https://github.com/lh3/minibwa) (read mapping, only with the experimental `--aligner minibwa`)

  > Li H, Homer N. Fast genomic read alignment with minibwa. arXiv:2606.15357. 2026. <!-- verify citation -->

- [SAMtools](https://pubmed.ncbi.nlm.nih.gov/19505943/) (sorting, filtering, indexing, flagstat; `mpileup` runs inside `freyja variants`)

  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. PubMed PMID: 19505943.

  > Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008.

- [BCFtools](https://github.com/samtools/bcftools) (`bcftools view`/`consensus` for the consensus genomes; `bcftools csq` when `--annotate_variants` is used)

  > Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008. doi: 10.1093/gigascience/giab008.

  > Danecek P, McCarthy SA. BCFtools/csq: haplotype-aware variant consequences. Bioinformatics. 2017 Jul 1;33(13):2037-2039. doi: 10.1093/bioinformatics/btx100 <!-- verify citation -->

- [iVar](https://github.com/andersen-lab/ivar) (primer trimming; variant calling inside `freyja variants`)

  > Grubaugh ND, Gangavarapu K, Quick J, Matteson NL, De Jesus JG, Main BJ, Tan AL, Paul LM, Brackney DE, Grewal S, Gurfield N, Van Rompay KKA, Isern S, Michael SF, Coffey LL, Loman NJ, Andersen KG. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biol. 2019 Jan 8;20(1):8. doi: 10.1186/s13059-018-1618-7.

- [Freyja](https://github.com/andersen-lab/Freyja) (`freyja variants`, `demix`, `aggregate`, `plot`, `update`, and `boot` with `--freyja_boot`)

  > Karthikeyan S, Levy JI, De Hoff P, et al. Wastewater sequencing reveals early cryptic SARS-CoV-2 variant transmission. Nature. 2022;609:101-108. doi: 10.1038/s41586-022-05049-6

- [Python](https://www.python.org/) (helper scripts in `bin/`, run in a Python container) and [Jinja2](https://jinja.palletsprojects.com/) (HTML report template)

  > Python Software Foundation. Python Language Reference. <https://www.python.org/> <!-- verify citation -->

  > Pallets. Jinja2. <https://github.com/pallets/jinja> <!-- verify citation -->

  Note: the report renderer (`RENDER_REPORT`) runs in a small image (Python + Jinja2) built from conda-forge packages with the
  Seqera Wave CLI and hosted in Seqera's community registry (`community.wave.seqera.io/library/python_jinja2`); it appears in
  `run_info.json` as a container. MultiQC is **not** used by this pipeline.

## Data sources

- **Freyja barcodes and lineage data** (`usher_barcodes.feather|csv`, `lineage_mutations.json`, `curated_lineages.json`,
  `lineages.yml`, fetched by `freyja update` or taken from the Freyja container). Derived from the UShER SARS-CoV-2
  mutation-annotated tree and the Pango lineage designations; the barcode version used by a run is recorded in
  `pipeline_info/run_info.json`.

  > Turakhia Y, Thornlow B, Hinrichs AS, De Maio N, Gozashti L, Lanfear R, Haussler D, Corbett-Detig R. Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic. Nat Genet. 2021 Jun;53(6):809-816. doi: 10.1038/s41588-021-00862-7.

  > Rambaut A, Holmes EC, O'Toole Á, Hill V, McCrone JT, Ruis C, du Plessis L, Pybus OG. A dynamic nomenclature proposal for SARS-CoV-2 lineages to assist genomic epidemiology. Nat Microbiol. 2020 Nov;5(11):1403-1407. doi: 10.1038/s41564-020-0770-5 <!-- verify citation -->

  > Gangavarapu K, Latif AA, Mullen JL, Alkuzweny M, Hufbauer E, Tsueng G, Haag E, Zeller M, Aceves CM, Zaiets K, Cano M, Zhou X, Qian Z, Sattler R, Matteson NL, Levy JI, Lee RTC, Freitas L, Maurer-Stroh S; GISAID Core and Curation Team; Suchard MA, Wu C, Su AI, Andersen KG, Hughes LD. Outbreak.info genomic reports: scalable and dynamic surveillance of SARS-CoV-2 variants and mutations. Nat Methods. 2023 Apr;20(4):512-522. doi: 10.1038/s41592-023-01769-3 <!-- verify citation -->

- **SARS-CoV-2 reference genome** (NC_045512.2, bundled in `assets/references/`)

  > Wu F, Zhao S, Yu B, Chen YM, Wang W, Song ZG, Hu Y, Tao ZW, Tian JH, Pei YY, Yuan ML, Zhang YL, Dai FH, Liu Y, Wang QM, Zheng JJ, Xu L, Holmes EC, Zhang YZ. A new coronavirus associated with human respiratory disease in China. Nature. 2020 Mar;579(7798):265-269. doi: 10.1038/s41586-020-2008-3.

- **Amplicon primer schemes** (bundled BED files in `assets/primer_schemes/`, copied from the original Snakemake pipeline).
  Cite the scheme you used.
  - ARTIC schemes (`articv3`, `articv4.1`, `articv5.3.2`): <https://github.com/artic-network/primer-schemes> <!-- verify citation -->
  - Midnight (`midnight`, `midnightv2`):

    > Freed NE, Vlková M, Faisal B, Silander OK. Rapid and inexpensive whole-genome sequencing of SARS-CoV-2 using 1200 bp tiled amplicons and Oxford Nanopore Rapid Barcoding. Biol Methods Protoc. 2020;5(1):bpaa014. doi: 10.1093/biomethods/bpaa014 <!-- verify citation -->

  - `eden`, `swift`, `qiagenv1`, `xgen`: use the vendor or authors' documentation for the scheme. <!-- verify citation: source of each bundled BED not recorded in this repository -->

## Original pipeline's workflow engine

The original pipeline was built with Snakemake; it is no longer used but is credited here.

> Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH, Sochat V, Forster J, Lee S, Twardziok SO, Kanitz A, Wilm A, Holtgrewe M, Rahmann S, Nahnsen S, Köster J. Sustainable data analysis with Snakemake. F1000Research. 2021;10:33. doi: 10.12688/f1000research.29032.2 <!-- verify citation -->

## Software packaging and containerisation

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Seqera Containers / Wave](https://seqera.io/containers/) (the fastp, bwa and samtools images are built there)

  > <https://seqera.io/containers/> <!-- verify citation: no paper known; check the project's preferred citation -->

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel D. Docker: lightweight Linux containers for consistent development and deployment. Linux Journal. 2014;2014(239):2.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

## Software versions

The exact version of every tool and container used by a run is recorded, not listed here:

- `pipeline_info/software_versions.yml`: tool versions collected from each process.
- `pipeline_info/run_info.json`: pipeline version and git commit, Nextflow version, container tags, Freyja version, all
  parameters, input checksums, and (for any run that demixes) the `freyja_data` block with the barcode version and file
  checksums.

Report these two files, or the versions they contain, alongside any results you publish. Container images are pinned to
exact tags in the module definitions under `modules/`.
