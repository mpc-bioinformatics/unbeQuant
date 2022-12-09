#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.sp_embl_file = "proteins.txt"  // Databasefile in SP-EMBL
params.outdir = "$PWD/results"  // Output-Directory of the FASTA-results

// Parameters for Protein-Graph Generation
params.features_in_graphs = "-ft None" // Features added to the Protein-Graph
params.peptide_limits = "--pep_miscleavages 2 --pep_min_pep_length 5 --pep_max_weight 6500" // Limits used for Exporting peptides

workflow {
    // Create Protein-Graphs ,statistics and sqlite_file
    sp_embl_file = Channel.fromPath(params.sp_embl_file)
    create_sqlite_fasta_database(sp_embl_file)
    
    // Export sqlite into fasta
    convert_pepsqlite_to_fasta(create_sqlite_fasta_database.out[0])
}

process create_sqlite_fasta_database {
    publishDir "${params.outdir}/", mode:'copy'

    input:
    file input_sp_embl

    output:
    file "peptides.db"
    file "proteins_statistics.csv"

    """
    protgraph -epepsqlite --pep_sqlite_database peptides.db -eo . ${params.features_in_graphs} ${params.peptide_limits} -cnp -cnpm -amw -o proteins_statistics.csv ${input_sp_embl}

    """
}

process convert_pepsqlite_to_fasta {
    publishDir "${params.outdir}/", mode:'copy'

    input:
    file database

    output:
    file "peptides.fasta"

    """
    protgraph_pepsqlite_to_fasta ${database} -o peptides.fasta
    """
}
