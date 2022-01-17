process REANNOTATE_REPEATS {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path unclassified_domain_table  // Renamed repeat consensus library
    path hmmscan_domain_table       // hmmscan output (*.tbl + *.pfamtbl = *.domtbl )

    output:
    // path "*.fasta", emit: fasta
    path "domain_table.tsv", emit: domain_table

    script:
    """
    cat $unclassified_domain_table $hmmscan_domain_table > domain_table.tsv
    """
}
