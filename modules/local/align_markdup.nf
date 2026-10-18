process ALIGN_MARKDUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.22.1 bioconda::htslib=1.22.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta), path(reads)
    tuple path(fasta), path(fai), path(bwa_index)

    output:
    tuple val(meta), path("${meta.id}.md.bam"), path("${meta.id}.md.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.markdup.txt"),                          emit: markdup
    tuple val(meta), path("${meta.id}.flagstat.txt"),                         emit: flagstat
    tuple val(meta), path("${meta.id}.samtools_stats.txt"),                   emit: stats
    path "versions.yml",                                                      emit: versions

    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:ILLUMINA"
    """
    INDEX=\$(find -L ${bwa_index}/ -name "*.amb" | sed 's/\\.amb\$//')

    bwa mem -t ${task.cpus} -R '${rg}' "\$INDEX" ${reads[0]} ${reads[1]} \\
        | samtools fixmate -m -u - - \\
        | samtools sort -u -@ ${task.cpus} -T tmp_sort - \\
        | samtools markdup -@ ${task.cpus} -f ${meta.id}.markdup.txt -T tmp_md - ${meta.id}.md.bam

    samtools index ${meta.id}.md.bam
    samtools flagstat ${meta.id}.md.bam > ${meta.id}.flagstat.txt
    samtools stats --reference ${fasta} ${meta.id}.md.bam > ${meta.id}.samtools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
