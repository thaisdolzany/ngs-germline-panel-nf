process FASTP {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::fastp=1.3.6"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d013aad5427d824afe472e6607ea47685ff0181f1fb09e52a179e0ec39e43e88/data' :
        'community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.trim.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"),          emit: json
    tuple val(meta), path("${meta.id}.fastp.html"),          emit: html
    path "versions.yml",                                     emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_R1.trim.fastq.gz \\
        --out2 ${meta.id}_R2.trim.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/^fastp //')
    END_VERSIONS
    """
}
