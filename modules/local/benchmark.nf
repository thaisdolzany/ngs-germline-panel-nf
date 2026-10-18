process BENCHMARK {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bcftools=1.23.1 bioconda::htslib=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(truth)
    tuple path(fasta), path(fai), path(bwa_index)
    path targets

    output:
    tuple val(meta), path("${meta.id}.benchmark.tsv"),  emit: tsv
    tuple val(meta), path("${meta.id}.discordant.tsv"), emit: discordant
    path "${meta.id}.benchmark_summary.tsv",            emit: summary
    path "versions.yml",                                emit: versions

    script:
    """
    benchmark_vcf.sh ${meta.id} ${vcf} ${truth} ${targets} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
