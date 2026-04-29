#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// Test Network Graph Building
// ============================================================================
// This workflow tests the network graph visualization process.
// Smaller test case to verify integration before running full pipeline.

// Test parameters for network graph
params.test_output_dir = "$PWD/results"
params.test_input_edges_json = "${params.test_output_dir}/feature_analysis/feature_data_lists/*_edges.json"

// Network graph parameters - TEST SETTINGS
params.mmf_build_network_graph = true  // Enable graph building
params.mmf_graph_edge_cutoff = null  // Use all edges
params.mmf_graph_mz_cutoff = null  // No m/z filtering
params.mmf_graph_rt_cutoff = null  // No RT filtering
params.mmf_graph_output_image = true  // Generate SVG visualization
params.mmf_graph_output_graphml = true  // Also export GraphML
params.mmf_graph_skip_analysis = false  // Run analysis
params.mmf_graph_test_fraction = 0.2  // Use 20% of edges for testing (faster)

// Output configuration
params.main_outdir = "$PWD/results"

// Import the map_mzml_features workflow
PROJECT_DIR = workflow.projectDir
include {map_mzml_features} from PROJECT_DIR + '/map_mzml_features.nf'
include {build_network_graph_visualization} from PROJECT_DIR + '/map_mzml_features.nf'

// ============================================================================
// STANDALONE TEST WORKFLOW
// ============================================================================

workflow {
    println """
    ╔════════════════════════════════════════════════════════════════════╗
    ║         Network Graph Building - Test Workflow                     ║
    ║                                                                    ║
    ║  This test workflow validates the network graph integration      ║
    ║  using existing feature pairing results.                         ║
    ╠════════════════════════════════════════════════════════════════════╣
    ║  Configuration:                                                   ║
    ║    Input edges JSON:  ${params.test_input_edges_json}
    ║    Output directory: ${params.test_output_dir}
    ║    Test fraction:     20% (${params.mmf_graph_test_fraction})
    ║    Output formats:    SVG + GraphML
    ║    Analysis:          ${params.mmf_graph_skip_analysis ? 'SKIPPED' : 'ENABLED'}
    ╚════════════════════════════════════════════════════════════════════╝
    """.stripIndent()
    
    // Check if input files exist
    edges_files = Channel.fromPath(params.test_input_edges_json)
      .map { it.toAbsolutePath().toString() }
    edges_files = edges_files.ifEmpty {
        error """
        ERROR: No edge JSON files found.

        Expected path: ${params.test_input_edges_json}

        Make sure to run the full workflow first to generate edge data:
          nextflow main_unbequant.nf

        Or provide existing edge files with:
          nextflow test_network_graph.nf --test_input_edges_json "path/to/*_edges.json"
        """.stripIndent()
    }

    graph_results = edges_files | build_network_graph_visualization
    
    // Output final message
    graph_results | view { result ->
        println """
        ╔════════════════════════════════════════════════════════════════════╗
        ║                    Test Completed Successfully!                    ║
        ╚════════════════════════════════════════════════════════════════════╝
        """
    }
}
