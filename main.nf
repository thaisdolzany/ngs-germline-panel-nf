#!/usr/bin/env nextflow
/*
 * ngs-germline-panel-nf
 *
 * Germline SNV/indel calling for targeted (gene panel) Illumina paired-end data:
 * read QC and trimming, alignment, duplicate marking, target coverage QC,
 * variant calling, annotation and, when a truth set is supplied, benchmarking.
 */

include { FASTQC            } from './modules/local/fastqc'
include { FASTP             } from './modules/local/fastp'
include { PREPARE_REFERENCE } from './modules/local/prepare_reference'
include { PREPARE_TARGETS   } from './modules/local/prepare_targets'
include { ALIGN_MARKDUP     } from './modules/local/align_markdup'
include { COVERAGE_QC       } from './modules/local/coverage_qc'
include { CALL_VARIANTS     } from './modules/local/call_variants'
include { ANNOTATE_VARIANTS } from './modules/local/annotate_variants'
include { BENCHMARK         } from './modules/local/benchmark'
include { MULTIQC           } from './modules/local/multiqc'

// Absolute paths and URLs are used as given; anything else is taken
// relative to the folder that holds the samplesheet.
def resolvePath(String p, sheet_dir) {
    def target = (p.startsWith('/') || p.contains('://')) ? p : "${sheet_dir}/${p}"
    return file(target, checkIfExists: true)
}

workflow {

    if (!params.input) { error "Missing --input (samplesheet CSV)" }
    if (!params.fasta) { error "Missing --fasta (reference genome)" }
    if (!params.targets) { error "Missing --targets (panel BED)" }

    log.info """
    ngs-germline-panel-nf
    ---------------------
    input     : ${params.input}
    fasta     : ${params.fasta}
    targets   : ${params.targets}
    dbsnp     : ${params.dbsnp ?: 'none'}
    min depth : ${params.min_depth}x
    outdir    : ${params.outdir}
    """.stripIndent()

    // Samplesheet: sample,fastq_1,fastq_2[,truth_vcf]
    // Relative paths are resolved against the samplesheet location.
    def sheet_dir = file(params.input).parent

    def ch_samples = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2) {
                error "Samplesheet rows need sample, fastq_1 and fastq_2: ${row}"
            }
            def meta  = [id: row.sample]
            def truth = row.truth_vcf ? resolvePath(row.truth_vcf, sheet_dir) : null
            [meta, [resolvePath(row.fastq_1, sheet_dir), resolvePath(row.fastq_2, sheet_dir)], truth]
        }

    def ch_reads = ch_samples.map { meta, reads, _truth -> [meta, reads] }
    def ch_truth = ch_samples.filter { _meta, _reads, truth -> truth != null }.map { meta, _reads, truth -> [meta, truth] }

    // Reference and targets
    def fasta = file(params.fasta, checkIfExists: true)
    PREPARE_REFERENCE(fasta)
    def ch_reference = PREPARE_REFERENCE.out.index.map { fai, bwa -> [fasta, fai, bwa] }

    PREPARE_TARGETS(file(params.targets, checkIfExists: true))
    def ch_targets    = PREPARE_TARGETS.out.bed
    def ch_targets_gz = PREPARE_TARGETS.out.bed_gz

    // Reads
    FASTQC(ch_reads)
    FASTP(ch_reads)

    // Alignment, duplicate marking, coverage
    ALIGN_MARKDUP(FASTP.out.reads, ch_reference)
    COVERAGE_QC(ALIGN_MARKDUP.out.bam, ch_targets, params.min_depth)

    // Variants
    CALL_VARIANTS(ALIGN_MARKDUP.out.bam, ch_reference, ch_targets, params.min_qual, params.min_variant_depth)

    // --dbsnp false (or none) turns rsID annotation off, e.g. to override the test profile
    def use_dbsnp = params.dbsnp && !(params.dbsnp.toString().toLowerCase() in ['false', 'none', 'null'])
    def dbsnp     = use_dbsnp ? file(params.dbsnp, checkIfExists: true) : []
    def dbsnp_tbi = use_dbsnp ? file(params.dbsnp_tbi ?: "${params.dbsnp}.tbi", checkIfExists: true) : []
    ANNOTATE_VARIANTS(CALL_VARIANTS.out.vcf, ch_targets_gz, dbsnp, dbsnp_tbi)

    // Benchmark against truth sets, when provided
    def ch_bench_in = CALL_VARIANTS.out.vcf
        .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
        .join(ch_truth.map { meta, truth -> [meta.id, truth] })
        .map { _id, meta, vcf, tbi, truth -> [meta, vcf, tbi, truth] }
    BENCHMARK(ch_bench_in, ch_reference, ch_targets)

    // Software versions, one entry per process
    channel.empty()
        .mix(
            PREPARE_REFERENCE.out.versions,
            PREPARE_TARGETS.out.versions,
            FASTQC.out.versions.first(),
            FASTP.out.versions.first(),
            ALIGN_MARKDUP.out.versions.first(),
            COVERAGE_QC.out.versions.first(),
            CALL_VARIANTS.out.versions.first(),
            ANNOTATE_VARIANTS.out.versions.first(),
            BENCHMARK.out.versions.first()
        )
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", newLine: false, sort: { f -> f.name })

    // Run-level tables
    COVERAGE_QC.out.genes
        .map { _meta, f -> f }
        .collectFile(name: 'coverage_genes_all_samples.tsv', storeDir: "${params.outdir}/summary", keepHeader: true, skip: 1, sort: { f -> f.name })
    ANNOTATE_VARIANTS.out.table
        .map { _meta, f -> f }
        .collectFile(name: 'variants_all_samples.tsv', storeDir: "${params.outdir}/summary", keepHeader: true, skip: 1, sort: { f -> f.name })
    BENCHMARK.out.tsv
        .map { _meta, f -> f }
        .collectFile(name: 'benchmark_all_samples.tsv', storeDir: "${params.outdir}/summary", keepHeader: true, skip: 1, sort: { f -> f.name })

    // MultiQC
    def ch_multiqc = channel.empty().mix(
        FASTQC.out.zip.map { _meta, f -> f },
        FASTP.out.json.map { _meta, f -> f },
        ALIGN_MARKDUP.out.flagstat.map { _meta, f -> f },
        ALIGN_MARKDUP.out.stats.map { _meta, f -> f },
        CALL_VARIANTS.out.stats.map { _meta, f -> f },
        COVERAGE_QC.out.summary.collectFile(name: 'panel_coverage_summary.tsv', keepHeader: true, skip: 1, sort: { f -> f.name }),
        BENCHMARK.out.summary.collectFile(name: 'benchmark_summary.tsv', keepHeader: true, skip: 1, sort: { f -> f.name })
    ).collect()

    MULTIQC(ch_multiqc, file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true))
}
