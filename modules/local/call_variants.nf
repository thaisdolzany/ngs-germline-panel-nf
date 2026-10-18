process CALL_VARIANTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bcftools=1.23.1 bioconda::htslib=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(bwa_index)
    path targets
    val min_qual
    val min_dp

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${meta.id}.bcftools_stats.txt"),                  emit: stats
    path "versions.yml",                                                     emit: versions

    script:
    // --max-idepth: bcftools skips indel calling where depth exceeds this value
    // (default 250). Panels routinely sit above that, so it is raised here.
    """
    bcftools mpileup \\
        --fasta-ref ${fasta} \\
        --regions-file ${targets} \\
        --annotate FORMAT/AD,FORMAT/DP \\
        --min-MQ 20 \\
        --min-BQ 20 \\
        --max-depth 5000 \\
        --max-idepth 5000 \\
        --threads ${task.cpus} \\
        -Ou ${bam} \\
        | bcftools call --multiallelic-caller --variants-only -Ou \\
        | bcftools norm --fasta-ref ${fasta} --multiallelics -any -Ou \\
        | bcftools filter --soft-filter LowQual --exclude 'QUAL<${min_qual} || FORMAT/DP<${min_dp}' \\
            -Oz -o ${meta.id}.vcf.gz

    bcftools index --tbi ${meta.id}.vcf.gz
    bcftools stats --fasta-ref ${fasta} ${meta.id}.vcf.gz > ${meta.id}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
