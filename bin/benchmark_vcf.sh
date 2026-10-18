#!/usr/bin/env bash
# Compare PASS calls against a truth VCF inside the panel targets.
#
# Usage: benchmark_vcf.sh <sample> <query.vcf.gz> <truth.vcf[.gz]> <targets.bed> <reference.fasta>
#
# Both VCFs are split to biallelic records and left-normalised against the
# reference, then matched on CHROM, POS, REF and ALT (bcftools isec -c none).
# This is an exact-allele comparison. It does not reconcile different
# representations of the same complex haplotype, which tools such as hap.py
# or vcfeval do.
#
# Writes:
#   <sample>.benchmark.tsv         TP/FP/FN, precision, recall, F1 by variant type
#   <sample>.benchmark_summary.tsv one row per type, sample-prefixed, for MultiQC
#   <sample>.discordant.tsv        every FP, FN and genotype mismatch
set -euo pipefail

sample=$1
query=$2
truth=$3
bed=$4
ref=$5

prep() {  # $1 input, $2 output prefix, $3 extra view args
    bcftools norm --fasta-ref "$ref" --multiallelics -any -Ou "$1" \
        | bcftools view --targets-file "$bed" $3 -Ou \
        | bcftools sort -Oz -o "$2.vcf.gz" 2>/dev/null
    bcftools index --tbi "$2.vcf.gz"
}

prep "$query" q "--apply-filters PASS"
prep "$truth" t ""

bcftools isec --collapse none -p isec q.vcf.gz t.vcf.gz

count() { bcftools view -H $2 "$1" | wc -l | tr -d ' '; }

{
    printf 'sample\ttype\tTP\tFP\tFN\tprecision\trecall\tF1\tgenotype_concordance\n'
    for type in SNV INDEL ALL; do
        case $type in
            SNV)   sel="--types snps" ;;
            INDEL) sel="--types indels" ;;
            ALL)   sel="" ;;
        esac
        tp=$(count isec/0002.vcf "$sel")
        fp=$(count isec/0000.vcf "$sel")
        fn=$(count isec/0001.vcf "$sel")
        # genotype agreement among TP sites
        gt_ok=$(join -t $'\t' \
            <(bcftools view $sel isec/0002.vcf | bcftools query -f '%CHROM:%POS:%REF:%ALT\t[%GT]\n' | sed 's/|/\//' | sort) \
            <(bcftools view $sel isec/0003.vcf | bcftools query -f '%CHROM:%POS:%REF:%ALT\t[%GT]\n' | sed 's/|/\//' | sort) \
            | awk -F'\t' '$2 == $3' | wc -l | tr -d ' ')
        awk -v s="$sample" -v t="$type" -v tp="$tp" -v fp="$fp" -v fn="$fn" -v g="$gt_ok" 'BEGIN {
            p = (tp + fp) ? tp / (tp + fp) : 0
            r = (tp + fn) ? tp / (tp + fn) : 0
            f = (p + r) ? 2 * p * r / (p + r) : 0
            gc = tp ? g / tp : 0
            printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n", s, t, tp, fp, fn, p, r, f, gc
        }'
    done
} > "${sample}.benchmark.tsv"

{
    printf 'Sample\tTP\tFP\tFN\tprecision\trecall\tF1\tgenotype_concordance\n'
    awk -F'\t' -v OFS='\t' 'NR > 1 { print $1 "_" $2, $3, $4, $5, $6, $7, $8, $9 }' "${sample}.benchmark.tsv"
} > "${sample}.benchmark_summary.tsv"

{
    printf 'sample\tclass\tchrom\tpos\tref\talt\tcall_gt\ttruth_gt\n'
    bcftools query -f "${sample}\tFP\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t.\n" isec/0000.vcf
    bcftools query -f "${sample}\tFN\t%CHROM\t%POS\t%REF\t%ALT\t.\t[%GT]\n" isec/0001.vcf
    join -t $'\t' \
        <(bcftools query -f '%CHROM:%POS:%REF:%ALT\t[%GT]\n' isec/0002.vcf | sed 's/|/\//' | sort) \
        <(bcftools query -f '%CHROM:%POS:%REF:%ALT\t[%GT]\n' isec/0003.vcf | sed 's/|/\//' | sort) \
        | awk -F'\t' -v OFS='\t' -v s="$sample" '$2 != $3 { split($1, k, ":"); print s, "GT_MISMATCH", k[1], k[2], k[3], k[4], $2, $3 }'
} > "${sample}.discordant.tsv"
