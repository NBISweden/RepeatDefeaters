process ADD_TREP_ANNOTATION {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path repeat_library            // Reannotated repeat consensus library
    path trep_blast_hits           // TREP blast results ( *.blastn.tsv )

    output:
    path "*", emit: annotated_hits

    script:
    """
    If any Unknown sequence in "*.renamed.fasta" generates positive hits then
    take the name of the hit with smallest evalue from trep.blastn.out
    First three letters of TREP sequences indicates the transposon superfamily,
    for example
    An Unknown sequence with positive hit to >DHH_Mpol_A_RND-1 would have
    a new annotation formatted as "* /DHH"

    """

}