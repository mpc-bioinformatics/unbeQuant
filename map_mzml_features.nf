#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// Map mzML Features - Subworkflow
// ============================================================================
// This workflow maps mzML spectra data with TSV feature information to create
// paired feature datasets across multiple files.

// Configuration parameters
params.mmf_generate_heatmap = false  // Generate heatmap images from mzML files
params.mmf_optimize_pairing = true  // Use optimized or basic KD-tree based pairing
params.mmf_best_match_only = true  // Keep only best match for each feature pair
params.mmf_match_cutoff = 0.0  // Minimum match score cutoff (0.0-1.0), discontinued
params.mmf_mz_cutoff = 0.2  // Maximum m/z coordinate difference for filtering (original space)
params.mmf_rt_cutoff = 100  // Maximum RT coordinate difference for filtering (original space)
params.mmf_edges_cutoff = null  // Maximum euclidean distance cutoff for filtering (scaled space)
params.mmf_normalize_coordinates = null  // Enable/disable coordinate normalization (null=auto)
params.mmf_postpair_normalize_coordinates = false  // Recalculate edge distances using normalized coords after pairing
params.mmf_distance_calc_before_scaling = false  // Calculate distance before coordinate scaling
params.mmf_normalize_edge_distances = false  // Normalize edge distances to [0, 1] range (default: disabled)
params.mmf_skip_json_output = false  // Skip JSON serialization for speed (default: disabled - JSON output enabled)
params.mmf_analyze_pep_idents = true  // Perform detailed pep_ident analysis (default: enabled)
params.mmf_round_up_to = 2  // Number of decimal places to round m/z values
params.mmf_log_scale = true  // Apply logarithmic scaling to intensities
params.mmf_scale_colors = true  // Scale colors based on min/max intensities
params.mmf_feature_mode = "CoM"  // Options: "rectangle", "CoM" (Center of Mass)
params.mmf_feature_center_diagnostic = true  // Generate diagnostic plots
params.mmf_intensity_method = "top3_sum"  // Intensity calculation method
params.mmf_compression_level = null  // HDF5 compression: null=none (fastest), 1=light, 4=balanced, 9=max
params.mmf_row_batch_size = 10  // Rows to process per batch (higher=faster but more RAM)
params.mmf_generate_html_report = true  // Generate HTML pairing report
params.mmf_build_network_graph = true  // Build and analyze network graph composition

// Network Graph Parameters (build_network_graph)
params.mmf_graph_edge_cutoff = null  // Maximum euclidean distance for graph edges
params.mmf_graph_mz_cutoff = null  // Maximum m/z difference for graph (overrides edge_cutoff if set)
params.mmf_graph_rt_cutoff = null  // Maximum RT difference for graph (overrides edge_cutoff if set)
params.mmf_graph_output_image = true  // Generate visualization image (SVG - fastest)
params.mmf_graph_output_graphml = true  // Export as GraphML format
params.mmf_graph_skip_analysis = false  // Skip detailed graph analysis
params.mmf_graph_test_fraction = null  // For testing: use fraction of edges (e.g., 0.1 for 10%)

// Output directories
params.mmf_outdir = "${params.main_outdir ?: "$PWD/results"}/feature_analysis"
params.mmf_heatmap_dir = "${params.mmf_outdir}/heatmaps"
params.mmf_features_dir = "${params.mmf_outdir}/feature_data_lists"
params.mmf_analysis_dir = "${params.mmf_outdir}/feature_analysis_plots"

// ============================================================================
// MEMORY CONFIGURATION NOTES
// ============================================================================
// Process memory requirements can be customized in nextflow.config:
//
//   process {
//       withName: 'map_mzml_features:create_heatmap_image' { memory = '16G' }
//       withName: 'map_mzml_features:pair_features' { memory = '12G' }
//   }
//
// Memory defaults (can be overridden):
//   - process_mzml_file: 16G
//   - create_heatmap_image: 20G (adjust based on image size and available RAM)
//   - extract_feature_data: 16G
//   - pair_features: 16G
//   - generate_pairing_report: 8G
//
// If running out of memory, try:
//   1. Disable heatmap generation: --mmf_generate_heatmap false
//   2. Reduce memory per process in nextflow.config
//   3. Process files one at a time instead of in parallel

// ============================================================================
// PROCESSES
// ============================================================================

process process_mzml_file {
    """
    Process a single mzML file to extract spectral data.
    """
    tag "${mzml.baseName}"
    container "luxii/unbequant:latest"
    memory "8G"

    input:
    file(mzml)

    output:
    tuple val(mzml.baseName), file("${mzml.baseName}_spectrum_data.pkl")

    script:
    """
    ${workflow.projectDir}/bin/process_mzml_file.py \
        --mzml ${mzml} \
        --output_pickle ${mzml.baseName}_spectrum_data.pkl \
        --round_up_to ${params.mmf_round_up_to}
    """
}

process create_heatmap_image {
    """
    Create heatmap image from processed mzML spectral data.
    Uses optimized HDF5 processing for memory efficiency, outputs PNG for visualization.
    """
    tag "${basename}"
    container "luxii/unbequant:latest"
    memory "8G"
    publishDir "${params.mmf_heatmap_dir}", mode: 'symlink'

    input:
    tuple file(spectrum_pkl), val(basename)

    output:
    file("${basename}_heatmap.png")

    script:
    def compression_arg = params.mmf_compression_level != null ? "--compression_level ${params.mmf_compression_level}" : ""
    """
    ${workflow.projectDir}/bin/create_heatmap_image_hdf5.py \
        --spectrum_pkl ${spectrum_pkl} \
        --output_png ${basename}_heatmap.png \
        --log_scale ${params.mmf_log_scale} \
        --scale_colors ${params.mmf_scale_colors} \
        --row_batch_size ${params.mmf_row_batch_size} \
        ${compression_arg}
    """
}

process extract_feature_data {
    """
    Extract feature data from TSV file.
    Intermediate output - stored in work directory.
    """
    tag "${basename}"
    container "luxii/unbequant:latest"
    memory "8G"

    input:
    tuple file(tsv), val(basename)

    output:
    tuple val(basename), path("${basename}_feature_data.pkl"), path("${basename}_feature_data.json")

    shell:
    '''
    !{workflow.projectDir}/bin/extract_feature_data.py \
        --tsv !{tsv} \
        --output_pkl !{basename}_feature_data.pkl \
        --output_json !{basename}_feature_data.json \
        --round_up_to !{params.mmf_round_up_to} \
        --feature_mode !{params.mmf_feature_mode} \
        --generate_diagnostic !{params.mmf_feature_center_diagnostic}
    '''
}

process pair_features {
    """
    Pair matching features across multiple files using KD-tree or basic matching.
    """
    tag "feature_pairing"
    container "luxii/unbequant:latest"
    memory "12G"
    publishDir "${params.mmf_features_dir}", mode: 'symlink'

    input:
    file(feature_files)

    output:
    tuple path("paired_features.pkl"), path("paired_features.json", optional: true), 
          path("paired_features_edges.pkl"), path("paired_features_edges.json", optional: true)

    script:
    def use_basic_flag = params.mmf_optimize_pairing ? "" : "--use-basic"
    def best_match_flag = params.mmf_best_match_only ? "--best_match_only" : ""
    def mz_cutoff_arg = params.mmf_mz_cutoff != null ? "--mz_cutoff ${params.mmf_mz_cutoff}" : ""
    def rt_cutoff_arg = params.mmf_rt_cutoff != null ? "--rt_cutoff ${params.mmf_rt_cutoff}" : ""
    def edges_cutoff_arg = params.mmf_edges_cutoff != null ? "--edges_cutoff ${params.mmf_edges_cutoff}" : ""
    def normalize_coords_arg = params.mmf_normalize_coordinates != null ? "--normalize_coordinates ${params.mmf_normalize_coordinates}" : ""
    def postpair_normalize_coords_flag = params.mmf_postpair_normalize_coordinates ? "--postpair_normalize_coordinates" : ""
    def distance_before_scaling_flag = params.mmf_distance_calc_before_scaling ? "--distance_calc_before_scaling" : ""
    def normalize_edge_distances_flag = params.mmf_normalize_edge_distances ? "--normalize_edge_distances" : ""
    def skip_json_output_flag = params.mmf_skip_json_output ? "--skip_json_output" : ""
    def analyze_pep_idents_flag = params.mmf_analyze_pep_idents ? "--analyze_pep_idents" : ""
    // Build the --input_files argument list from all files in the input
    def input_files_arg = feature_files.collect { "${it}" }.join(" ")
    
    """
    echo "DEBUG: Input feature files to process:"
    for f in ${input_files_arg}; do
        [ -f "\$f" ] && ls -lh "\$f" || echo "WARNING: File not found: \$f"
    done | head -3
    
    ${workflow.projectDir}/bin/pair_features.py \
        --input_files ${input_files_arg} \
        --output_pkl paired_features.pkl \
        --output_json paired_features.json \
        --output_edges_pkl paired_features_edges.pkl \
        --match_cutoff ${params.mmf_match_cutoff} \
        ${use_basic_flag} \
        ${best_match_flag} \
        ${mz_cutoff_arg} \
        ${rt_cutoff_arg} \
        ${edges_cutoff_arg} \
        ${normalize_coords_arg} \
        ${postpair_normalize_coords_flag} \
        ${distance_before_scaling_flag} \
        ${normalize_edge_distances_flag} \
        ${skip_json_output_flag} \
        ${analyze_pep_idents_flag}
    
    echo "DEBUG: Output files created:"
    ls -lah paired_features.pkl paired_features_edges.pkl 2>&1
    [ -f paired_features.json ] && echo "  ✓ paired_features.json" || echo "  ✗ paired_features.json (optional, skipped)"
    [ -f paired_features_edges.json ] && echo "  ✓ paired_features_edges.json" || echo "  ✗ paired_features_edges.json (optional, skipped)"
    
    """
}

process feature_pairing_report {
    """
    Generate comprehensive summary report of feature pairing results.
    Outputs text, CSV, and optional HTML formats with statistics and visualizations.
    """
    tag "pairing_report"
    container "luxii/unbequant:latest"
    memory "8G"
    publishDir "${params.mmf_features_dir}", mode: 'symlink'

    input:
    tuple path(paired_pkl), path(paired_json), path(edges_pkl), path(edges_json)

    output:
    tuple path("pairing_summary.txt"), path("pairing_statistics.csv"), 
          path("pairing_report.html", optional: true)

    script:
    def html_arg = params.mmf_generate_html_report ? "--output_html pairing_report.html" : ""
    
    """
    ${workflow.projectDir}/bin/generate_pairing_report.py \
        --paired_json ${paired_json} \
        --edges_json ${edges_json} \
        --output_summary pairing_summary.txt \
        --output_stats pairing_statistics.csv \
        ${html_arg}
    """
}

process build_network_graph_visualization {
    """
    Build and visualize network graph from paired feature edges using igraph.
    Analyze graph composition and generate network visualizations.
    """
    tag "network_graph"
    memory "8G"
    publishDir "${params.mmf_features_dir}", mode: 'copy'

    input:
    val(edges_json)

    output:
    tuple path("network_graph*.svg", optional: true), 
        path("network_graph*.graphml", optional: true),
        path("graph_analysis*.json", optional: true),
        path("graph_composition*.json", optional: true),
        path("degree_distribution*.png", optional: true)

    script:
    def edge_cutoff_arg = params.mmf_graph_edge_cutoff != null ? "--edge_cutoff ${params.mmf_graph_edge_cutoff}" : ""
    def mz_cutoff_arg = params.mmf_graph_mz_cutoff != null ? "--mz_cutoff ${params.mmf_graph_mz_cutoff}" : ""
    def rt_cutoff_arg = params.mmf_graph_rt_cutoff != null ? "--rt_cutoff ${params.mmf_graph_rt_cutoff}" : ""
    def output_image_arg = params.mmf_graph_output_image ? "--output_image network_graph.svg" : ""
    def output_graphml_arg = params.mmf_graph_output_graphml ? "--output_graphml network_graph.graphml" : ""
    def output_analysis_arg = "--output_analysis graph_analysis.json"
    def output_composition_arg = "--output_composition graph_composition.json"
    def output_histogram_arg = "--output_histogram degree_distribution.png"
    def skip_analysis_flag = params.mmf_graph_skip_analysis ? "--skip-analysis" : ""
    def test_fraction_arg = params.mmf_graph_test_fraction != null ? "--fraction ${params.mmf_graph_test_fraction}" : ""
    
    """
    python3 ${workflow.projectDir}/bin/build_network_graph.py \
        --input_pkl ${edges_json} \
        ${output_image_arg} \
        ${output_graphml_arg} \
        ${output_analysis_arg} \
        ${output_composition_arg} \
        ${output_histogram_arg} \
        ${edge_cutoff_arg} \
        ${mz_cutoff_arg} \
        ${rt_cutoff_arg} \
        ${skip_analysis_flag} \
        ${test_fraction_arg}
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow map_mzml_features {
    take:
        mzml_files      // Channel: mzML files
        tsv_files       // Channel: TSV files (tuples with mzML basename for matching)

    main:
        // Match mzML and TSV files
        mzml_named = mzml_files.map { mzml -> 
            tuple(mzml.baseName, mzml) 
        }
        
        tsv_named = tsv_files.map { tsv -> 
            tuple(tsv.baseName, tsv) 
        }
        
        // Join on basename to create file pairs
        file_pairs = mzml_named.join(tsv_named, failOnMismatch: false)
        
        // Step 1: Generate heatmaps if enabled (optional, for visualization only)
        if (params.mmf_generate_heatmap) {
            processed_mzml = process_mzml_file(file_pairs.map { name, mzml, tsv -> mzml })
            heatmap_input = processed_mzml
                .map { basename, pkl -> tuple(pkl, basename) }
            heatmap_output = create_heatmap_image(heatmap_input)
        }
        
        // Step 2: Extract feature data from all TSV files
        extraction_input = file_pairs
            .map { name, mzml, tsv -> tuple(tsv, name) }
        feature_data = extract_feature_data(extraction_input)
        
        // Step 3: Collect all feature data files
        all_feature_files = feature_data
            .map { name, pkl, json -> pkl }
            .collect()
        
        // Step 4: Pair features across files (always run, will handle single file)
        paired_output = pair_features(all_feature_files)
        
        // Step 5: Generate report with edges and composition data
        report_input = paired_output.map { pkl, json, edges_pkl, edges_json -> 
            tuple(pkl, json, edges_pkl, edges_json) 
        }
        pairing_report = feature_pairing_report(report_input)
        
        // Step 6: Build network graph visualization if enabled
        graph_output = null
        if (params.mmf_build_network_graph) {
            graph_edges_json = paired_output.map { pkl, json, e_pkl, e_json -> e_json.toAbsolutePath().toString() }
            graph_output = build_network_graph_visualization(graph_edges_json)
        }
        
        emit:
            feature_data = feature_data
            paired_features = paired_output.map { pkl, json, e_pkl, e_json -> pkl }
            paired_json = paired_output.map { pkl, json, e_pkl, e_json -> json }
            edges_pkl = paired_output.map { pkl, json, e_pkl, e_json -> e_pkl }.flatten()
            edges_json = paired_output.map { pkl, json, e_pkl, e_json -> e_json }.flatten()
            pairing_summary = pairing_report.map { summary, stats, html -> summary }
            pairing_stats = pairing_report.map { summary, stats, html -> stats }
            pairing_html = pairing_report.map { summary, stats, html -> html }
            graph_svg = graph_output ? graph_output.map { svg, graphml, analysis, comp -> svg } : Channel.empty()
            graph_graphml = graph_output ? graph_output.map { svg, graphml, analysis, comp -> graphml } : Channel.empty()
            graph_analysis = graph_output ? graph_output.map { svg, graphml, analysis, comp -> analysis } : Channel.empty()
            graph_composition = graph_output ? graph_output.map { svg, graphml, analysis, comp -> comp } : Channel.empty()
}
