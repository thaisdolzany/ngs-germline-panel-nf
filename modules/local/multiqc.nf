process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc=1.35"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data' :
        'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc' }"

    input:
    path qc_files, stageAs: 'qc/*'
    path config

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    multiqc --force --config ${config} --filename multiqc_report.html qc/
    """
}
