process PREPARE_REFERENCE {
    tag "${fasta.simpleName}"
    label 'process_medium'

    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.22.1 bioconda::htslib=1.22.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    path fasta

    output:
    tuple path("${fasta}.fai"), path("bwa"), emit: index
    path "versions.yml",                     emit: versions

    script:
    """
    samtools faidx ${fasta}
    mkdir bwa
    bwa index -p bwa/${fasta.baseName} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
