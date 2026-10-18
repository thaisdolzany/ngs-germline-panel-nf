#!/usr/bin/env bash
# Target-level coverage QC from `samtools depth -a -b <bed>` output.
#
# Usage: coverage_qc.sh <sample> <targets.bed> <depth.tsv> <min_depth>
#
# The BED must be sorted by start within each chromosome and non-overlapping
# (the pipeline sorts and merges it before calling this script).
#
# Writes:
#   <sample>.coverage_targets.tsv   one row per target interval
#   <sample>.coverage_genes.tsv     one row per gene (4th BED column)
#   <sample>.low_coverage.bed       targets with at least one base below <min_depth>
#   <sample>.coverage_summary.tsv   one-line summary (MultiQC general stats)
set -euo pipefail

sample=$1
bed=$2
depth=$3
min=$4

: > "${sample}.low_coverage.bed"

awk -v OFS='\t' -v sample="$sample" -v min="$min" \
    -v t_out="${sample}.coverage_targets.tsv" \
    -v g_out="${sample}.coverage_genes.tsv" \
    -v low_out="${sample}.low_coverage.bed" \
    -v mqc_out="${sample}.coverage_summary.tsv" '
FNR == NR {
    n++
    chrom[n] = $1; start[n] = $2 + 0; end[n] = $3 + 0; gene[n] = (NF >= 4 ? $4 : ".")
    c = $1; k = ++cn[c]; idx[c, k] = n
    next
}
{
    c = $1; pos = $2 + 0; d = $3 + 0
    if (!(c in cn)) next
    # binary search: last target on this chromosome with start < pos
    lo = 1; hi = cn[c]; hit = 0
    while (lo <= hi) {
        mid = int((lo + hi) / 2); i = idx[c, mid]
        if (start[i] < pos) { hit = i; lo = mid + 1 } else { hi = mid - 1 }
    }
    if (!hit || pos > end[hit]) next
    sum[hit] += d; seen[hit]++
    if (d >= min) ok[hit]++
    if (d >= 100) ge100[hit]++
    if (!(hit in mn) || d < mn[hit]) mn[hit] = d
}
END {
    print "sample", "chrom", "start", "end", "gene", "length", "mean_depth", "min_depth", "pct_ge_" min "x", "pct_ge_100x" > t_out
    for (i = 1; i <= n; i++) {
        len = end[i] - start[i]
        mean = (len ? sum[i] / len : 0)
        printf "%s\t%s\t%d\t%d\t%s\t%d\t%.1f\t%d\t%.2f\t%.2f\n", sample, chrom[i], start[i], end[i], gene[i], len, mean, (seen[i] == len ? mn[i] : 0), 100 * ok[i] / len, 100 * ge100[i] / len > t_out
        if (ok[i] < len) printf "%s\t%d\t%d\t%s\n", chrom[i], start[i], end[i], gene[i] > low_out
        g = gene[i]
        if (!(g in glen)) gorder[++ng] = g
        glen[g] += len; gsum[g] += sum[i]; gok[g] += ok[i]
        tot_len += len; tot_sum += sum[i]; tot_ok += ok[i]; tot_100 += ge100[i]
    }
    print "sample", "gene", "target_bp", "mean_depth", "pct_ge_" min "x" > g_out
    for (j = 1; j <= ng; j++) {
        g = gorder[j]
        printf "%s\t%s\t%d\t%.1f\t%.2f\n", sample, g, glen[g], gsum[g] / glen[g], 100 * gok[g] / glen[g] > g_out
    }
    print "Sample", "mean_target_depth", "pct_target_ge_min_depth", "pct_target_ge_100x" > mqc_out
    printf "%s\t%.1f\t%.2f\t%.2f\n", sample, tot_sum / tot_len, 100 * tot_ok / tot_len, 100 * tot_100 / tot_len > mqc_out
}' "$bed" "$depth"
