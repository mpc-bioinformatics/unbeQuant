#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.one_fasta = true  // Generate a fasta for the whole dataset OR a fasta per mzML
params.mgf_folder = "$PWD/MZMLs"  // Input-Directory of mzMLs, which should be used to generate FASTA-files
params.sp_embl_file = "proteins.txt"  // Databasefile in SP-EMBL
params.sp_embl_file_2 = "<empty>"  // Databasefile in SP-EMBL for the second Protein-Graph generation run
params.outdir = "$PWD/results"  // Output-Directory of the FASTA-results

// Parameters for csv query generation
params.max_precursor = 6500  // Upper limit of search query. Precursors higher then this value (in Dalton) will be ommitted
params.query_ppm = 5  // Set ppm in measured precursors in MS2

// Parameters for Protein-Graph Generation
params.elbpcsr = 32 // Number of Intervals per Protein-Graph-Node
params.fixed_modifications = "-fm 'C:57.021464'"  // Fixed modifications (directly added to CLI in ProtGraph)
params.variable_modifications = "-vm 'M:15.994915'"  // Variable modifications (directly added to CLI in ProtGraph)
params.features_in_graphs = "-ft ALL" // Features added to the Protein-Graph
params.protgraph_run_2 = false  // Flag, to run a second run with ProtGraph (to append the database with more Protein-Graphs)
params.features_in_graphs_2 = "-ft NONE" // Features added to the Protein-Graph (in the second Protein-Graph generation run)

// Parameters for Protein-Graph-Traversal
params.num_procs_traversal = -1  // Number of process used by the graph traversal to generate a fasta (CAUTION: This can be very resource intensive!)
params.use_floats = 1  // Bool wheather to use floats for mono weights (1 --> use floats, 0 --> do not use floats)
params.limit_varcount = -1  // Decide wheather to limit the varcount (Default -1 --> Use up to infinite many variants)

workflow {
    // Get all mzML files and prepare queries file
    mgfs = Channel.fromPath(params.mgf_folder + "/*.mgf")
    generate_query_csvs(mgfs)
    // Decide whether we want a single fasta
    // or a fasta per file
    if (params.one_fasta) {
        concat_query_csvs(generate_query_csvs.out.collect())
        results = optimize_query(concat_query_csvs.out)
    } else {
        results = optimize_query(generate_query_csvs.out)
    }

    // Create Protein-Graphs
    sp_embl_file = Channel.fromPath(params.sp_embl_file)
    sp_embl_file_2 = Channel.fromPath(params.sp_embl_file_2)
    create_protein_graphs(sp_embl_file, sp_embl_file_2)

    // Execute ProtGraphCpp (graph traversal) and generate a/multiple fasta files
    db_query = create_protein_graphs.out[0].combine(results)
    create_fasta(db_query)
    compact_fasta(create_fasta.out)
}

process generate_query_csvs {

    input:
    file input_mgf

    output:
    file "${input_mgf.baseName}.csv"

    """
    01_mgf_parse_queries.py ${input_mgf} ${input_mgf.baseName}.csv ${params.max_precursor} ${params.query_ppm}
    """
}

process concat_query_csvs {

    input:
    file input_mgf_csv

    output:
    file "concatenated_query_csv.csv"

    """
    cat ${input_mgf_csv} > concatenated_query_csv.csv
    """
}

process optimize_query {

    input:
    file input_csv

    output:
    file "${input_csv.baseName}_optimized.csv"

    """
    02_query_optimize.py ${input_csv} -o ${input_csv.baseName}_optimized.csv
    """
}

process create_protein_graphs {
    publishDir "${params.outdir}/", mode:'copy'

    input:
    file input_sp_embl
    file input_sp_embl_2

    output:
    file "database.bpcsr"
    file "proteins_statistics.csv"
    file "proteins_statistics_2.csv"

    """
    protgraph -elbpcsr -eo . -elbpcsr_pdbs ${params.elbpcsr} ${params.fixed_modifications} ${params.variable_modifications} ${params.features_in_graphs} -cnp -cnpvar -amw -o proteins_statistics.csv ${input_sp_embl}

    touch proteins_statistics_2.csv
    if ${params.protgraph_run_2}
    then
        protgraph -elbpcsr -eo . -elbpcsr_pdbs ${params.elbpcsr} ${params.fixed_modifications} ${params.variable_modifications} ${params.features_in_graphs_2} -cnp -cnpvar -amw -o proteins_statistics_2.csv ${input_sp_embl_2}
    fi
    """
}

process create_fasta {

    input:
    tuple file(database), file(csv_query)

    output:
    file "${csv_query.baseName}.fasta"

    """
    if [ ${params.use_floats} -eq 1 ]
    then
        protgraphtraversefloat ${database} ${csv_query} ${params.num_procs_traversal} ${csv_query.baseName}.fasta ${params.limit_varcount}
    else
        protgraphtraverseint   ${database} ${csv_query} ${params.num_procs_traversal} ${csv_query.baseName}.fasta ${params.limit_varcount}
    fi
    """
}

process compact_fasta {
    publishDir "${params.outdir}/", mode:'copy'

    input:
    file input_fasta

    output:
    file "peptides.fasta"

    """
    protgraph_compact_fasta ${input_fasta} -o peptides.fasta
    """
}
