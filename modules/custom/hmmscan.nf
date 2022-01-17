process HMMSCAN {

    conda (params.enable_conda ? 'bioconda::pfam_scan==1.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pfam_scan:1.6--hdfd78af_4' :
        'quay.io/biocontainers/pfam_scan:1.6--hdfd78af_4' }"

    input:
    path fasta
    path hmm_db

    output:
    path '*.hmmscan.out', emit: hmmscan_out
    path '*.tbl'        , emit: tbl
    path '*.pfamtbl'    , emit: pfam_table
    path '*.domtbl'    , emit: domain_table
    path "versions.yml" , emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    # Stage local copy of HMM DB
    mkdir -p HMM_DB
    gzip -cdf $hmm_db > HMM_DB/\$( basename $hmm_db .gz )
    find HMM_DB -name "*.hmm" -exec hmmpress {} \\;

    hmmscan $args \\
        --cpu $task.cpus \\
        --pfamtblout ${prefix}.pfamtbl \\
        --tblout ${prefix}.tbl \\
        -o ${prefix}.hmmscan.out \\
        HMM_DB/*.hmm \\
        $fasta

    # Combine tbl and pfamtbl data to
    # unclassified_domain_table format
    # col field               file-index   join-output-col
    #   1 <seq id>            1.4          1
    #   2 <alignment start>   2.9          2
    #   3 <alignment end>     2.10         3
    #   4 <envelope start>    2.7          4
    #   5 <envelope end>      2.8          5
    #   6 <hmm acc>           "-"
    #   7 <hmm name>          1.2          6
    #   8 <type>              "-"
    #   9 <hmm start>         2.11         7
    #  10 <hmm end>           2.12         8
    #  11 <hmm length>        `calculate`
    #  12 <bit score>         1.10         9
    #  13 <E-value>           1.9          10
    #  14 <significance>      "1"
    #  15 <clan>              "-"

    join -1 1 -2 1 -o1.4,2.9,2.10,2.7,2.8,1.2,2.11,2.12,1.10,1.9 \
        <( grep -v -e '^[[:space:]]*$' -e '^#' "$hmmscan_table" | \
            awk '{ print \$1"-"\$8"-"\$9" "\$0 } ' | sort -k1,1 ) \
        <( grep -v -e '^[[:space:]]*$' -e '^#' "$hmmscan_pfamtbl" | \
            awk 'NF > 10 { print \$1"-"\$3"-"\$2" "\$0 } ' | sort -k1,1 ) | \
        awk '{ \$5=\$5 "\\t-\\t"; \$6=\$6 "\\t-\\t"; \$8=\$8 "\\t-\\t"; \$10=\$10 "\\t1\\t-"; print \$0 }' \
        > ${prefix}.domtbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmscan: \$( hmmscan -h | sed '2 !d;s/[^0-9]*\\(\\([0-9]\\.\\)\\{0,4\\}[0-9][^.]\\).*/\\1/' )
    END_VERSIONS
    """

}
