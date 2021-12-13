process REANNOTATE_REPEATS {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path repeat_library            // Renamed repeat consensus library
    path pfam_table                // Pfam output (*.pfamtbl)

    output:
    path "*", emit: annotation

    script:
    """
    echo something
    """
}