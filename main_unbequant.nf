#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters (Input-Files)
params.main_fasta_file = "" // The FASTA file of the species to be searched (downloadable from UniProtKB)
params.main_raw_files_folder = ""  // The folder where the RAW/.d-files are located.
params.main_comet_params = ""  // The comet parameter file for search. NOTE: Here the digestion should be explicitly turned on (or set appropiately, depending on the input FASTA.)).

// Optional Parameters
// See each nextflow script.

// Output Parameters which are fixed for the folder structure.
// NOTE: these can be changed and also individually turned off. See for the corresponding nextflow scripts
params.main_outdir = "$PWD/results"
params.ctm_outdir =  "${params.main_outdir}/mgfs_mzmls"
params.idc_outdir =  "${params.main_outdir}/identifications"
params.sir_outdir =  "${params.main_outdir}/identifications_summarized"
params.qal_outdir =  "${params.main_outdir}/quantifications"
params.outdir =      "${params.main_outdir}/extraced_xics"

// Set Parameters, since we use a FASTA (from UniProt most probably) which was not generated with ProtGraph
params.sir_identification_from_protgraph = false
params.sir_remove_variable_modifications = true
params.sir_count_same_protein_as_unique = true


// Import Workflows
PROJECT_DIR = workflow.projectDir
include {convert_to_mgf} from PROJECT_DIR + '/include/ProGFASTAGen/convert_to_mgf.nf'
include {convert_to_mzml} from PROJECT_DIR + '/include/ProGFASTAGen/convert_to_mzml.nf'
include {identification_via_comet} from PROJECT_DIR + '/include/ProGFASTAGen/identification_via_comet.nf'
include {summarize_ident_results} from PROJECT_DIR + '/include/ProGFASTAGen/summarize_ident_results.nf'
include {quantify_and_align} from PROJECT_DIR + '/quantify_and_align.nf'


// Standalone MAIN Workflow
workflow {
	fasta_file = Channel.fromPath(params.main_fasta_file)
    raw_files = Channel.fromPath(params.main_raw_files_folder  + "/*.raw")
    d_files = Channel.fromPath(params.main_raw_files_folder  + "/*.d", type: "dir")
    comet_params = Channel.fromPath(params.main_comet_params)

    main_unbequant(
        fasta_file,
        raw_files,
        d_files,
        comet_params
    )
}

// Importable MAIN Workflow
workflow main_unbequant {
    take:
        fasta_file
        raw_files
        d_files
        comet_parameters_file
    main:
        // Generate MGF-Files
        convert_to_mgf(raw_files, d_files)
        convert_to_mzml(raw_files, d_files)

        // Search via Comet (+ Percolator if set)
        identification_via_comet(
            convert_to_mgf.out,
            fasta_file,
            comet_parameters_file
        )

        // Group multiple FDRs
        grouped_fdrs = identification_via_comet.out.groupTuple()
        final_grouped_results = grouped_fdrs.map { tuple("___" + it[0].toString() + "_fdr",  it [1]) }
        
        // For each FDR get the identification summaries
        summarize_ident_results(final_grouped_results)
        
        // Get quantification results
        quantify_and_align(
            raw_files.concat(d_files), 
            convert_to_mzml.out,
            final_grouped_results
        )
}
