process REDUNDANT_HITS {

    conda (params.enable_conda ? 'conda-forge::sed=4.7' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    blast_tsv

    output:
    path "*.tsv", emit: tsv

    script:
    """
    grep Unknown $blast_tsv | \
        awk ' NR >= 1 {
            \$8=(\$6)/(\$2)
        }
        1
        \$8 > 0.6
        \$1 != \$3
        !a[\$5\$6]++
        ' > self_comparison.tsv
    """
}