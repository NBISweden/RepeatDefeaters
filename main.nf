#! /usr/bin/env nextflow

// Enable DSL2 syntax for Nextflow
nextflow.enable.dsl = 2

// Print parameters to screen before running workflow.
log.info("""
NBIS support 5861

 Annotation of unclassified TEs
 ===================================
""")

// Check a project allocation is given for running on Uppmax clusters.
if(workflow.profile == "uppmax" && !params.project){
    exit 1, "Please provide a SNIC project number ( --project )!\n"
}

include { BLAST_BLASTX as BLAST_POSITIVE_STRAND } from './modules/blast_blastx'     , addParams(options:params.modules['blast_positive_strand'])
include { BLAST_BLASTX as BLAST_NEGATIVE_STRAND } from './modules/blast_blastx'     , addParams(options:params.modules['blast_negative_strand'])
include { FILTER_BLAST_XML                      } from './modules/blast_xml_filter' , addParams(options:[:])
include { PFAM_SCAN                             } from './modules/pfam_scan'        , addParams(options:params.modules['pfam_scan'])
include { FILTER_PFAM                           } from './modules/pfam_filter'      , addParams(options:[:])
include { ANNOTATION                            } from './modules/annotation'       , addParams(options:[:])
include { PFAM_TRANSPOSIBLE_ELEMENT_SEARCH      } from './modules/pfam_te_search'   , addParams(options:[:])
include { DIVSUM                                } from './modules/divsum'           , addParams(options:[:])

// The main workflow
workflow {

    main:
        // Get data
        Channel.fromPath(params.samples)
            .ifEmpty { exit 1, "Cannot find reads from ${params.samples}!\n" }
            .set { readpairs }

        BLAST_POSITIVE_STRAND()
        BLAST_NEGATIVE_STRAND()
        FILTER_BLAST_XML()
        PFAM_SCAN()
        PFAM_TRANSPOSIBLE_ELEMENT_SEARCH()
        ANNOTATION()
        DIVSUM()

}
