process ANNOTATE_VARIANTS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.23.1 bioconda::htslib=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(targets_gz), path(targets_tbi)
    path dbsnp
    path dbsnp_tbi

    output:
    tuple val(meta), path("${meta.id}.annotated.vcf.gz"), path("${meta.id}.annotated.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${meta.id}.variants.tsv"),                                             emit: table
    path "versions.yml",                                                                          emit: versions

    script:
    def id_step = dbsnp ?
        "bcftools annotate --annotations ${dbsnp} --columns ID -Oz -o ${meta.id}.annotated.vcf.gz genes.vcf.gz" :
        "mv genes.vcf.gz ${meta.id}.annotated.vcf.gz"
    """
    echo '##INFO=<ID=GENE,Number=.,Type=String,Description="Panel gene(s) overlapping the variant">' > gene.hdr

    # gene symbol from the panel BED
    bcftools annotate \\
        --annotations ${targets_gz} \\
        --columns CHROM,FROM,TO,INFO/GENE \\
        --header-lines gene.hdr \\
        -Oz -o genes.vcf.gz \\
        ${vcf}
    bcftools index --tbi genes.vcf.gz

    # rsID from dbSNP, when provided (both files must be indexed for this step)
    ${id_step}
    bcftools index --tbi ${meta.id}.annotated.vcf.gz

    printf 'sample\\tchrom\\tpos\\tid\\tref\\talt\\tgene\\tqual\\tfilter\\tgenotype\\tdepth\\tallele_depths\\n' > ${meta.id}.variants.tsv
    bcftools query \\
        -f '${meta.id}\\t%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/GENE\\t%QUAL\\t%FILTER[\\t%GT\\t%DP\\t%AD]\\n' \\
        ${meta.id}.annotated.vcf.gz >> ${meta.id}.variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
