process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::fastqc=0.13.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/816bda6dd014e810753f642611232443ddeabb2b47f9681ad43f9044816f00ff/data' :
        'community.wave.seqera.io/library/fastqc:0.13.0--b342db7b694c9af4' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"),  emit: zip
    path "versions.yml",                    emit: versions

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/^FastQC v//')
    END_VERSIONS
    """
}
