process PREPARE_TARGETS {
    tag "${bed.simpleName}"
    label 'process_single'

    conda "bioconda::bcftools=1.23.1 bioconda::htslib=1.23.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    path bed

    output:
    path "targets.bed",                                 emit: bed
    tuple path("targets.bed.gz"), path("targets.bed.gz.tbi"), emit: bed_gz
    path "versions.yml",                                emit: versions

    script:
    """
    # keep the first four columns, drop headers, sort and merge overlapping intervals
    grep -v -E '^(#|track|browser)' ${bed} \\
        | awk -v OFS='\\t' 'NF >= 3 { print \$1, \$2, \$3, (NF >= 4 ? \$4 : ".") }' \\
        | sort -k1,1 -k2,2n \\
        | awk -v OFS='\\t' '
            \$1 == c && \$2 <= e { if (\$3 > e) e = \$3; if (index("," g ",", "," \$4 ",") == 0) g = g "," \$4; next }
            { if (c != "") print c, s, e, g; c = \$1; s = \$2; e = \$3; g = \$4 }
            END { if (c != "") print c, s, e, g }' \\
        > targets.bed

    bgzip -c targets.bed > targets.bed.gz
    tabix -p bed targets.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$(tabix --version | head -1 | sed 's/^tabix (htslib) //')
    END_VERSIONS
    """
}
