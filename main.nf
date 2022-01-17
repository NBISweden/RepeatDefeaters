#! /usr/bin/env nextflow

// Enable DSL2 syntax for Nextflow
nextflow.enable.dsl = 2

include { RENAME_SEQUENCES as RENAME_REPEAT_MODELER_SEQUENCES } from './modules/custom/rename_sequences'
include { PFAM_TRANSPOSIBLE_ELEMENT_SEARCH                    } from './modules/custom/pfam_transposible_element_search'
include { MAKE_BLAST_DB as BUILD_PROTEIN_REF_BLAST_DB ;
          MAKE_BLAST_DB as BUILD_TREP_BLAST_DB ;
          MAKE_BLAST_DB as BUILD_ANNOTATED_LIB_BLAST_DB       } from './modules/custom/make_blast_db'
include { BLASTX_AND_FILTER                                   } from './modules/custom/blastx_and_filter'
include { PFAM_SCAN                                           } from './modules/custom/pfam_scan'
include { HMMSCAN as CUSTOM_HMM_SCAN                          } from './modules/custom/hmmscan'
include { ANNOTATE_REPEATS                                    } from './modules/custom/annotate_repeats'
include { BLASTN as TREP_BLASTN ;
          BLASTN as RECIPROCAL_BLASTN                         } from './modules/custom/blastn'
include { ADD_TREP_ANNOTATION                                 } from './modules/custom/add_trep_annotation'
include { REANNOTATE_REPEATS                                  } from './modules/custom/reannotate_repeats'
include { REDUNDANT_HITS                                      } from './modules/custom/redundant_hits'

log.info("""
NBIS support 5861

 Annotation of unclassified TEs
 ===================================
""")

// Check a project allocation is given for running on Uppmax clusters.
if( workflow.profile == "uppmax" && !params.project ){
    exit 1, "Please provide a SNIC project number ( --project )!\n"
}

// The primary analysis workflow
workflow {

    main:
        // Step 1: Change Fasta headers
        RENAME_REPEAT_MODELER_SEQUENCES (
            file (
                params.repeat_modeler_fasta,
                checkIfExists: true
            ),
            params.species_short_name
        )
        versions_ch = Channel.empty()

        // Step 2: Get ID's of PFAM Proteins with TE domains
        protein_te_domain_list = Channel.empty()
        if ( params.pfam_proteins_with_te_domain_list ){
            protein_te_domain_list = file (
                params.pfam_proteins_with_te_domain_list,
                checkIfExists: true
            )
        } else {
            PFAM_TRANSPOSIBLE_ELEMENT_SEARCH (
                file (
                    params.pfam_a_db,
                    checkIfExists: true
                ),
                file (
                    params.transposon_keywords,
                    checkIfExists: true
                ),
                file (
                    params.transposon_blacklist,
                    checkIfExists: true
                )
            )
            protein_te_domain_list = PFAM_TRANSPOSIBLE_ELEMENT_SEARCH.out.te_domain_proteins
            versions_ch = PFAM_TRANSPOSIBLE_ELEMENT_SEARCH.out.versions.first()
        }

        // Step 3: Strand specific Blast search of Repeats against
        // a protein reference database
        BUILD_PROTEIN_REF_BLAST_DB (
            Channel.fromPath(
                params.protein_reference,
                checkIfExists: true
            )
            .splitFasta( by: 1000000000 )   // Use high number to keep file intact
            .collectFile(                   // Merge multiple fasta's into one
                name: 'protein_reference.fasta'
            )
        )
        BLASTX_AND_FILTER (
            RENAME_REPEAT_MODELER_SEQUENCES.out.fasta,
            BUILD_PROTEIN_REF_BLAST_DB.out.db,
            ['plus','minus']
        )

        // Step 4: Scan Repeats for PFAM domains
        PFAM_SCAN (
            BLASTX_AND_FILTER.out.fasta,
            Channel.fromPath (
                [ params.pfam_hmm_db, params.pfam_hmm_dat ],
                checkIfExists: true
            ).collect()
        )

        // Step 5: Reannoate Fasta headers to emphasize
        // repeats with single-stranded non-TE domains
        ANNOTATE_REPEATS (
            RENAME_REPEAT_MODELER_SEQUENCES.out.fasta,
            PFAM_SCAN.out.pfam_table.collect(),
            protein_te_domain_list,
            params.species_short_name
        )

        // Step 6: TREP Annotation
        BUILD_TREP_BLAST_DB (
            Channel.fromPath( params.trep_db, checkIfExists: true )
            .splitFasta( by: 1000000000 )   // Auto unzip and keep files fasta intact
            .collectFile(                 // Merge multiple fasta's into one
                name: 'trep_db.fasta'
            )
        )
        TREP_BLASTN (
            ANNOTATE_REPEATS.out.fasta,
            BUILD_TREP_BLAST_DB.out.db
        )
        ADD_TREP_ANNOTATION (
            ANNOTATE_REPEATS.out.fasta,
            TREP_BLASTN.out.tsv
        ) // take Smallest e-value and rename with triplet added

        // Step 7: Custom HMM search
        CUSTOM_HMM_SCAN (
            BLASTX_AND_FILTER.out.fasta,
            Channel.fromPath(
                params.custom_hmms,
                checkIfExists: true
            ).collect()
        )
        REANNOTATE_REPEATS (
            ANNOTATE_REPEATS.out.unclassified_with_te_domains,
            CUSTOM_HMM_SCAN.out.domain_table.collect()
        )

        // Step 8: Reciprocal blast
        // BUILD_ANNOTATED_LIB_BLAST_DB ( REANNOTATE_REPEATS.out.fasta )
        // RECIPROCAL_BLASTN (
        //     REANNOTATE_REPEATS.out.fasta,
        //     BUILD_ANNOTATED_LIB_BLAST_DB.out.db
        // )
        // REDUNDANT_HITS ( RECIPROCAL_BLASTN.out.tsv )
        /*
    # The goal is to further annotate "*.renamed.fasta"
    ## Nucleotide identity search against TREP
    wget http://botserv2.uzh.ch/kelldata/trep-db/downloads/trep-db_nr_Rel-19.fasta.gz
    gunzip trep-db_nr_Rel-19.fasta.gz
    makeblastdb -in trep-db_nr_Rel-19.fasta -dbtype nucl
    blastn -db trep-db_nr_Rel-19.fasta -num_threads 8  -query "*.renamed.fasta" -outfmt 6 > trep.blastn.out

    ### If any Unknown sequence in "*.renamed.fasta" generates positive hits then
    ### take the name of the hit with smallest evalue from trep.blastn.out
    ### First three letters of TREP sequences indicates the transposon superfamily,
    ### for example
    ### An Unknown sequence with positive hit to >DHH_Mpol_A_RND-1 would have
    ### a new annotation formatted as "* /DHH"

    ## Custom HMM search
    ### INPUT: *.{plus,minus}.predicted.fasta from 03_Blastx
    pfam_scan.pl -fasta INPUT -dir Libraries/HMMs/ -outfile $INPUT.pfamtbl -e_seq 0.001
    ### Concatenate two $INPUT.pfamtbl output files with
    ### *.Unclassified_consensus_TEs from "05_Reannotated_Repeat_modeler_sequences".
    ### Unknown consensus sequences get annotated based on
    ### their best hit in the concatenated table. HMM best hits inherit HMM names.
    ### Pfam best hits inherit Pfam domain names similar to 05_Reannotated_Repeat_modeler_sequences
    ### THIS IS THE FINAL ANNOTATED LIBRARY.

    Custom HMM search
        INPUT: *.{plus,minus}.predicted.fasta from 03_Blastx
        pfam_scan.pl -fasta INPUT -dir Libraries/HMMs/ -outfile $INPUT.pfamtbl -e_seq 0.001
        Concatenate two $INPUT.pfamtbl output files with
        *.Unclassified_consensus_TEs from "05_Reannotated_Repeat_modeler_sequences".
        Unknown consensus sequences get annotated based on
        their best hit in the concatenated table. HMM best hits inherit HMM names.
        Pfam best hits inherit Pfam domain names similar to 05_Reannotated_Repeat_modeler_sequences
        THIS IS THE FINAL ANNOTATED LIBRARY.

    Show redundancy in the annotated library using blastn
        make database
        makeblastdb -in FINAL-ANNOTATED-LIBRARY -dbtype nucl
        Define custom output format. For queries and subjects that are both Unknown
        Get total matched bases/query length, which indicates the coverage of an overlap
        Filter coverage using 0.6 as cut-off then dedup
        blastn -db FINAL-ANNOTATED-LIBRARY
        -query FINAL-ANNOTATED-LIBRARY
        -outfmt "6 qseqid qlen sseqid slen length nident evalue"|
        grep Unknown|awk 'NR>=1{$8=($6)/($2)}1'|
        awk '$8>0.6'|awk '$1!=$3'|awk '!a[$5$6]++' > self.comparison

    ## Show redundancy in the annotated library using blastn
    ### make database
    makeblastdb -in FINAL-ANNOTATED-LIBRARY -dbtype nucl
    ### Define custom output format. For queries and subjects that are both Unknown
    ### Get total matched bases/query length, which indicates the coverage of an overlap
    ### Filter coverage using 0.6 as cut-off then dedup
    blastn -db FINAL-ANNOTATED-LIBRARY
    -query FINAL-ANNOTATED-LIBRARY
    -outfmt "6 qseqid qlen sseqid slen length nident evalue"|
    grep Unknown|awk 'NR>=1{$8=($6)/($2)}1'|
    awk '$8>0.6'|awk '$1!=$3'|awk '!a[$5$6]++' > self.comparison


        */

        // Report: Record software versions
        versions_ch.mix(
            RENAME_REPEAT_MODELER_SEQUENCES.out.versions.first(),
            BUILD_PROTEIN_REF_BLAST_DB.out.versions.first(),
            BLASTX_AND_FILTER.out.versions.first(),
            PFAM_SCAN.out.versions.first(),
            ANNOTATE_REPEATS.out.versions.first()
        ).collectFile(
            name: "software_versions.yml",
            cache: false,
            sort: true,
            storeDir: "${params.results}/pipeline_info"
        )

}






