process COVERAGE_QC {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::samtools=1.24 bioconda::htslib=1.24"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path targets
    val min_depth

    output:
    tuple val(meta), path("${meta.id}.coverage_targets.tsv"), emit: targets
    tuple val(meta), path("${meta.id}.coverage_genes.tsv"),   emit: genes
    tuple val(meta), path("${meta.id}.low_coverage.bed"),     emit: low
    path "${meta.id}.coverage_summary.tsv",                   emit: summary
    path "versions.yml",                                      emit: versions

    script:
    """
    # same read and base filters as the caller: MAPQ >= 20, base quality >= 20,
    # duplicates and secondary alignments excluded (samtools default), and
    # overlapping mates counted once (-s) so depth reflects unique fragments
    samtools depth -a -s -Q 20 -q 20 -b ${targets} ${bam} > ${meta.id}.depth.tsv
    coverage_qc.sh ${meta.id} ${targets} ${meta.id}.depth.tsv ${min_depth}
    rm ${meta.id}.depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
