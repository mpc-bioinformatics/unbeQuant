#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters (Input-Files)
params.main_fasta_file = "/workspaces/unbeQuant/files/UP000000589_10090.fasta" // The FASTA file of the species to be searched (downloadable from UniProtKB)
params.main_raw_files_folder = "/workspaces/unbeQuant/files/raws"  // The folder where the RAW/.d-files are located.
params.main_comet_params = "/workspaces/unbeQuant/example_configurations/comet_config.txt"  // The comet parameter file for search. NOTE: Here the digestion should be explicitly turned on (or set appropiately, depending on the input FASTA.)).

// Optional Parameters
// See each nextflow script for detailed parameter documentation.
//
// Feature Pairing Parameters (map_mzml_features):
//   --mmf_optimize_pairing: Use optimized KD-tree pairing (default: true) or basic pairing
//   --mmf_best_match_only: Keep only best match for each feature pair (default: false)
//   --mmf_match_cutoff: Minimum match score cutoff 0.0-1.0 (default: 0.0)
//   --mmf_mz_cutoff: Maximum m/z coordinate difference for filtering (uses original coordinates when set)
//   --mmf_rt_cutoff: Maximum RT coordinate difference for filtering (uses original coordinates when set)
//   --mmf_edges_cutoff: Maximum euclidean distance cutoff for filtering (uses scaled coordinates)
//   --mmf_normalize_coordinates: Enable/disable coordinate normalization before KD-tree
//                               (null=auto: OFF for mz/rt cutoff, ON for euclidean cutoff)
//   --mmf_distance_calc_before_scaling: Calculate distance before coordinate scaling (default: false)
//   --mmf_normalize_edge_distances: Normalize edge distances to [0, 1] range for visualization (default: false)
//   --mmf_skip_json_output: Skip JSON serialization for speed (default: false - JSON enabled by default)
//   --mmf_analyze_pep_idents: Perform detailed pep_ident matching analysis (default: true - enabled by default)
//   --mmf_generate_html_report: Generate HTML pairing report (default: true)
//   --mmf_build_network_graph: Build and analyze network graph composition (default: false)
//
// Example usages:
//   With coordinate-based filtering:
//     nextflow run main_unbequant.nf --mmf_mz_cutoff 0.1 --mmf_rt_cutoff 0.2
//
//   With euclidean distance filtering:
//     nextflow run main_unbequant.nf --mmf_edges_cutoff 0.5
//
//   With manual coordinate normalization override:
//     nextflow run main_unbequant.nf --mmf_normalize_coordinates false --mmf_edges_cutoff 0.5

// Output Parameters which are fixed for the folder structure.
// NOTE: these can be changed and also individually turned off. See for the corresponding nextflow scripts
params.main_outdir = "$PWD/results"
params.ctm_outdir =  "${params.main_outdir}/mgfs_mzmls"
params.idc_outdir =  "${params.main_outdir}/identifications"
params.sir_outdir =  "${params.main_outdir}/identifications_summarized"
params.qal_outdir =  "${params.main_outdir}/quantifications"
params.outdir =      "${params.main_outdir}/extraced_xics"
params.map_alignment_features_tsv_dir = "${params.qal_outdir}/features_with_annotated_identifications"
params.map_alignment_outdir = "${params.main_outdir}/heatmaps_with_features"

// Set Parameters, since we use a FASTA (from UniProt most probably) which was not generated with ProtGraph
params.sir_identification_from_protgraph = false
params.sir_remove_variable_modifications = true
params.sir_count_same_protein_as_unique = true


// Import Workflows
PROJECT_DIR = workflow.projectDir
include {convert_to_mgf} from PROJECT_DIR + '/ProGFASTAGen/convert_to_mgf.nf'
include {convert_to_mzml} from PROJECT_DIR + '/ProGFASTAGen/convert_to_mzml.nf'
include {identification_via_comet} from PROJECT_DIR + '/ProGFASTAGen/identification_via_comet.nf'
include {summarize_ident_results} from PROJECT_DIR + '/ProGFASTAGen/summarize_ident_results.nf'
include {quantify_and_align} from PROJECT_DIR + '/quantify_and_align.nf'
include {map_mzml_features} from PROJECT_DIR + '/map_mzml_features.nf'



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
        // Map mzML features with identification data
        mzml_files = convert_to_mzml.out.collect().flatten()
        tsv_files = quantify_and_align.out[1]
            .map { fdr, consensus_xml, feature_tsvs -> feature_tsvs }
            .collect()
            .flatten()
        
        // Capture subworkflow outputs to enable proper caching with -resume
        map_mzml_features(mzml_files, tsv_files)
}
