#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// Map mzML Features - Subworkflow
// ============================================================================
// This workflow maps mzML spectra data with TSV feature information to create
// paired feature datasets across multiple files.

// Configuration parameters
params.mmf_generate_heatmap = true  // Generate heatmap images from mzML files
params.mmf_optimize_pairing = true  // Use optimized or basic KD-tree based pairing
params.mmf_best_match_only = true  // Keep only best match for each feature pair
params.mmf_match_cutoff = 0.0  // Minimum match score cutoff (0.0-1.0), discontinued
params.mmf_mz_cutoff = 0.2  // Maximum m/z coordinate difference for filtering (original space)
params.mmf_rt_cutoff = 100  // Maximum RT coordinate difference for filtering (original space)
params.mmf_edges_cutoff = null  // Maximum euclidean distance cutoff for filtering (scaled space)
params.mmf_normalize_coordinates = null  // Enable/disable coordinate normalization (null=auto)
params.mmf_postpair_normalize_coordinates = false  // Rescale axes to balance ranges and add rescaled distance fields (original preserved). Only has an effekt on the distances and the clustering. The Pairing is determined by mmf_normalize_coordinates and for distinct mz and rt cuttoffs disabled by default, since both cutoff points are used to apply a already "scaled" cutoff different in x and y axis
params.mmf_distance_calc_before_scaling = false  // Calculate distance before coordinate scaling
params.mmf_normalize_edge_distances = false  // Normalize edge distances to [0, 1] range (default: disabled)
params.mmf_skip_json_output = false  // Skip JSON serialization for speed (default: disabled - JSON output enabled)
params.mmf_analyze_pep_idents = true  // Perform detailed pep_ident analysis (default: enabled)

//Heatmap parameters
params.mmf_round_up_to = 2  // Number of decimal places to round m/z values (nur für heatmaps)
params.mmf_log_scale = true  // Apply logarithmic scaling to intensities
params.mmf_scale_colors = true  // Scale colors based on min/max intensities
params.mmf_feature_mode = "CoM"  // Options: "rectangle", "CoM" (Center of Mass)
params.mmf_feature_center_diagnostic = false  // Generate diagnostic plots
params.mmf_intensity_method = "top3_sum"  // Intensity calculation method
params.mmf_compression_level = null  // HDF5 compression: null=none (fastest), 1=light, 4=balanced, 9=max
params.mmf_row_batch_size = 10  // Rows to process per batch (higher=faster but more RAM)
params.mmf_heatmap_features_only = true  // Save only mapped feature points without heatmap background (creates much smaller file)
params.mmf_heatmap_save_report = true  // Generate HTML report with heatmap and statistics
// Trimming parameters for feature extraction (optional, can help remove noisy edges)
params.mmf_rt_start_trim = 900  // Seconds to trim from the start of retention time (0 = no trimming)
params.mmf_rt_end_trim = 6000  // Maximum RT value to keep (keep only features with RT <= this value, 0 = no limit)
//Pairing Parameters
params.mmf_generate_html_report = true  // Generate HTML pairing report

// Network Graph Parameters (build_network_graph)
params.mmf_build_network_graph = true  // Build and analyze network graph composition
params.mmf_graph_edge_cutoff = null  // Maximum euclidean distance for graph edges (carefull when cutoff is already set in pairing step)
params.mmf_graph_mz_cutoff = null  // Maximum m/z difference for graph (overrides edge_cutoff if set) - careful when cutoff is already set in pairing step
params.mmf_graph_rt_cutoff = null  // Maximum RT difference for graph (overrides edge_cutoff if set) - careful when cutoff is already set in pairing step
params.mmf_graph_output_image = true  // Generate visualization image (SVG - fastest)
params.mmf_graph_output_graphml = true  // Export as GraphML format
params.mmf_graph_output_gml = null  // Export as GML file path (optional)
params.mmf_graph_output_edgelist = null  // Export as edge list file path (optional)
params.mmf_graph_output_analysis = null  // Export as analysis JSON file path (optional)
params.mmf_graph_output_composition = null  // Export as graph composition JSON file path (optional)
params.mmf_graph_output_histogram = null  // Export as degree distribution histogram file path (optional)
params.mmf_graph_skip_analysis = false  // Skip detailed graph analysis
params.mmf_graph_skip_visualization = true  // Skip graph visualization rendering
params.mmf_graph_test_fraction = 1  // For testing: use fraction of edges (e.g., 0.1 for 10%)
params.mmf_graph_layout_engine = 'sfdp'  // Graphviz layout engine: dot, sfdp, fdp, neato, twopi, circo
params.mmf_graph_layout_use_weights = true  // Use edge distances to influence layout (neato/fdp/sfdp only)
params.mmf_graph_use_graphviz_clusters = true  // Group clusters into Graphviz subgraphs
params.mmf_graph_use_coordinate_layout = false  // Use original feature coordinates as initial positions (neato/fdp only)
params.mmf_graph_feature_data_jsons = []  // List of feature data JSON files for coordinate-based layout
// Clustering Parameters
params.mmf_enable_clustering = true  // Enable community detection clustering for HTML report
params.mmf_clustering_method = 'edge_betweenness'  // Clustering method: louvain, walktrap, label_propagation, edge_betweenness
params.mmf_clustering_use_weights = true  // Use edge distances as weights in clustering
params.mmf_clustering_weight_mode = 'inverse'  // Weight mode: inverse (1/(1+distance)) or distance
params.mmf_output_clusters_json = true  // Export as JSON file with cluster information (optional)
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
    In features-only mode, overlays feature boxes and centers on a white background.
    """
    tag "${basename}"
    container "luxii/unbequant:latest"
    memory "8G"
    publishDir "${params.mmf_heatmap_dir}", mode: 'symlink'

    input:
    tuple file(spectrum_pkl), val(basename), file(feature_json)

    output:
    path("${basename}_heatmap.png", optional: true)
    path("${basename}_heatmap_report.html", optional: true)

    script:
    def compression_arg = params.mmf_compression_level != null ? "--compression_level ${params.mmf_compression_level}" : ""
    def features_only_arg = params.mmf_heatmap_features_only ? "--features_only" : ""
    def save_report_arg = params.mmf_heatmap_save_report ? "--save_report" : ""
    def feature_json_arg = params.mmf_heatmap_features_only && feature_json ? "--feature_data_json ${feature_json}" : ""
    """
    ${workflow.projectDir}/bin/create_heatmap_image_hdf5.py \
        --spectrum_pkl ${spectrum_pkl} \
        --output_png ${basename}_heatmap.png \
        --log_scale ${params.mmf_log_scale} \
        --scale_colors ${params.mmf_scale_colors} \
        --row_batch_size ${params.mmf_row_batch_size} \
        ${compression_arg} \
        ${features_only_arg} \
        ${save_report_arg} \
        ${feature_json_arg}
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
        --generate_diagnostic !{params.mmf_feature_center_diagnostic} \
        --rt_start_trim !{params.mmf_rt_start_trim} \
        --rt_end_trim !{params.mmf_rt_end_trim}
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
    tag "pairing_report"
    container "luxii/unbequant:latest"
    memory "8G"
    publishDir "${params.mmf_features_dir}", mode: 'symlink'

    input:
    tuple path(paired_pkl), path(paired_json), path(edges_pkl), path(edges_json), val(clusters_json)

    output:
    tuple path("pairing_summary.txt"), path("pairing_statistics.csv"), 
          path("pairing_report.html", optional: true)

    script:
    def html_arg = params.mmf_generate_html_report ? "--output_html pairing_report.html" : ""
    def clusters_arg = clusters_json ? "--clusters_json ${clusters_json}" : ""
    
    // Pairing parameters
    def mz_cutoff_arg = params.mmf_mz_cutoff != null ? "--mz_cutoff ${params.mmf_mz_cutoff}" : ""
    def rt_cutoff_arg = params.mmf_rt_cutoff != null ? "--rt_cutoff ${params.mmf_rt_cutoff}" : ""
    def edges_cutoff_arg = params.mmf_edges_cutoff != null ? "--edges_cutoff ${params.mmf_edges_cutoff}" : ""
    def best_match_only_arg = params.mmf_best_match_only ? "--best_match_only" : ""
    def optimize_pairing_arg = params.mmf_optimize_pairing ? "--optimize_pairing" : ""
    
    // Trimming parameters
    def rt_start_trim_arg = params.mmf_rt_start_trim ? "--rt_start_trim ${params.mmf_rt_start_trim}" : ""
    def rt_end_trim_arg = params.mmf_rt_end_trim ? "--rt_end_trim ${params.mmf_rt_end_trim}" : ""
    
    // Network graph parameters
    def graph_edge_cutoff_arg = params.mmf_graph_edge_cutoff != null ? "--graph_edge_cutoff ${params.mmf_graph_edge_cutoff}" : ""
    def graph_mz_cutoff_arg = params.mmf_graph_mz_cutoff != null ? "--graph_mz_cutoff ${params.mmf_graph_mz_cutoff}" : ""
    def graph_rt_cutoff_arg = params.mmf_graph_rt_cutoff != null ? "--graph_rt_cutoff ${params.mmf_graph_rt_cutoff}" : ""
    def graph_test_fraction_arg = params.mmf_graph_test_fraction != 1 ? "--graph_test_fraction ${params.mmf_graph_test_fraction}" : ""
    def graph_layout_engine_arg = params.mmf_graph_layout_engine ? "--graph_layout_engine ${params.mmf_graph_layout_engine}" : ""
    
    // Clustering parameters
    def enable_clustering_arg = params.mmf_enable_clustering ? "--enable_clustering" : ""
    def clustering_method_arg = params.mmf_enable_clustering ? "--clustering_method ${params.mmf_clustering_method}" : ""
    def clustering_use_weights_arg = params.mmf_clustering_use_weights ? "--clustering_use_weights" : ""
    def clustering_weight_mode_arg = params.mmf_clustering_weight_mode ? "--clustering_weight_mode ${params.mmf_clustering_weight_mode}" : ""
    
    // Post-pairing normalization
    def postpair_normalize_coords_arg = params.mmf_postpair_normalize_coordinates ? "--postpair_normalize_coordinates" : ""
    
    // Build input files list for tracking
    def input_files_arg = "--input_files ${paired_json} ${edges_json}"
    
    """
    ${workflow.projectDir}/bin/generate_pairing_report.py \
        --paired_json ${paired_json} \
        --edges_json ${edges_json} \
        --output_summary pairing_summary.txt \
        --output_stats pairing_statistics.csv \
        ${html_arg} \
        ${clusters_arg} \
        ${mz_cutoff_arg} \
        ${rt_cutoff_arg} \
        ${edges_cutoff_arg} \
        ${best_match_only_arg} \
        ${optimize_pairing_arg} \
        ${rt_start_trim_arg} \
        ${rt_end_trim_arg} \
        ${graph_edge_cutoff_arg} \
        ${graph_mz_cutoff_arg} \
        ${graph_rt_cutoff_arg} \
        ${graph_test_fraction_arg} \
        ${graph_layout_engine_arg} \
        ${enable_clustering_arg} \
        ${clustering_method_arg} \
        ${clustering_use_weights_arg} \
        ${clustering_weight_mode_arg} \
        ${postpair_normalize_coords_arg} \
        ${input_files_arg}
    """
}

process build_network_graph_visualization {
    tag "network_graph"
    memory "8G"
    publishDir "${params.mmf_features_dir}", mode: 'copy'

    input:
    val(edges_json)
    path feature_data_jsons, stageAs: 'feature_*.json'

    output:
    tuple path("network_graph*.svg", optional: true), 
        path("network_graph*.graphml", optional: true),
        path("network_graph*.gml", optional: true),
        path("network_graph*.edgelist", optional: true),
        path("graph_analysis*.json", optional: true),
        path("graph_composition*.json", optional: true),
        path("*clusters*.json", optional: true),
        path("degree_distribution*.png", optional: true)

    script:
    def edge_cutoff_arg = params.mmf_graph_edge_cutoff != null ? "--edge_cutoff ${params.mmf_graph_edge_cutoff}" : ""
    def mz_cutoff_arg = params.mmf_graph_mz_cutoff != null ? "--mz_cutoff ${params.mmf_graph_mz_cutoff}" : ""
    def rt_cutoff_arg = params.mmf_graph_rt_cutoff != null ? "--rt_cutoff ${params.mmf_graph_rt_cutoff}" : ""
    def output_image_arg = params.mmf_graph_output_image ? "--output_image network_graph.svg" : ""
    def output_graphml_arg = params.mmf_graph_output_graphml ? "--output_graphml network_graph.graphml" : ""
    def output_gml_arg = params.mmf_graph_output_graphml ? "--output_gml network_graph.gml" : ""
    def output_edgelist_arg = params.mmf_graph_output_edgelist ? "--output_edgelist network_graph.edgelist" : ""
    def output_analysis_arg = params.mmf_graph_output_analysis ? "--output_analysis ${params.mmf_graph_output_analysis}" : "--output_analysis graph_analysis.json"
    def output_composition_arg = params.mmf_graph_output_composition ? "--output_composition ${params.mmf_graph_output_composition}" : "--output_composition graph_composition.json"
    def output_histogram_arg = params.mmf_graph_output_histogram ? "--output_histogram ${params.mmf_graph_output_histogram}" : "--output_histogram degree_distribution.png"
    def skip_analysis_flag = params.mmf_graph_skip_analysis ? "--skip-analysis" : ""
    def skip_visualization_flag = params.mmf_graph_skip_visualization ? "--skip-visualization" : ""
    def test_fraction_arg = params.mmf_graph_test_fraction != null ? "--fraction ${params.mmf_graph_test_fraction}" : ""
    def enable_clustering_arg = params.mmf_enable_clustering ? "--enable-clustering" : ""
    def clustering_method_arg = params.mmf_enable_clustering ? "--clustering-method ${params.mmf_clustering_method}" : ""
    def clustering_weights_arg = params.mmf_clustering_use_weights ? "--clustering-use-weights" : "--clustering-no-weights"
    def clustering_weight_mode_arg = params.mmf_clustering_weight_mode ? "--clustering-weight-mode ${params.mmf_clustering_weight_mode}" : ""
    def output_clusters_arg = params.mmf_output_clusters_json ? "--output-clusters-json clusters.json" : ""
    def layout_engine_arg = params.mmf_graph_layout_engine ? "--layout_engine ${params.mmf_graph_layout_engine}" : ""
    def layout_use_weights_arg = params.mmf_graph_layout_use_weights ? "--layout_use_weights" : ""
    def graphviz_clusters_arg = params.mmf_graph_use_graphviz_clusters ? "--use-graphviz-clusters" : ""
    def coordinate_layout_arg = params.mmf_graph_use_coordinate_layout ? "--use-coordinate-layout --feature-data-jsons feature_*.json" : ""
    
    """
    python3 ${workflow.projectDir}/bin/build_network_graph.py \
        --input_pkl ${edges_json} \
        ${output_image_arg} \
        ${output_graphml_arg} \
        ${output_gml_arg} \
        ${output_edgelist_arg} \
        ${output_analysis_arg} \
        ${output_composition_arg} \
        ${output_histogram_arg} \
        ${edge_cutoff_arg} \
        ${mz_cutoff_arg} \
        ${rt_cutoff_arg} \
        ${skip_analysis_flag} \
        ${skip_visualization_flag} \
        ${test_fraction_arg} \
        ${layout_engine_arg} \
        ${layout_use_weights_arg} \
        ${graphviz_clusters_arg} \
        ${coordinate_layout_arg} \
        ${enable_clustering_arg} \
        ${clustering_method_arg} \
        ${clustering_weights_arg} \
        ${clustering_weight_mode_arg} \
        ${output_clusters_arg}
    
    # Create placeholder files if outputs were skipped/not generated
    # SVG placeholder
    if ! ls network_graph*.svg 1> /dev/null 2>&1; then
        touch network_graph_placeholder.svg
    fi
    
    # GraphML placeholder
    if ! ls network_graph*.graphml 1> /dev/null 2>&1; then
        echo '<?xml version="1.0" encoding="UTF-8"?><graphml></graphml>' > network_graph_placeholder.graphml
    fi
    
    # GML placeholder
    if ! ls network_graph*.gml 1> /dev/null 2>&1; then
        echo 'graph []' > network_graph_placeholder.gml
    fi
    
    # Edgelist placeholder
    if ! ls network_graph*.edgelist 1> /dev/null 2>&1; then
        touch network_graph_placeholder.edgelist
    fi
    
    # Analysis JSON placeholder
    if ! ls graph_analysis*.json 1> /dev/null 2>&1; then
        echo '{}' > graph_analysis_placeholder.json
    fi
    
    # Composition JSON placeholder
    if ! ls graph_composition*.json 1> /dev/null 2>&1; then
        echo '{}' > graph_composition_placeholder.json
    fi
    
    # Clusters JSON placeholder
    if ! ls *clusters*.json 1> /dev/null 2>&1; then
        echo '{}' > clusters_placeholder.json
    fi
    
    # Histogram PNG placeholder
    if ! ls degree_distribution*.png 1> /dev/null 2>&1; then
        touch degree_distribution_placeholder.png
    fi
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
        
        // Step 1: Process mzML files and extract feature data first
        processed_mzml = process_mzml_file(file_pairs.map { name, mzml, tsv -> mzml })
        
        // Step 2: Extract feature data from all TSV files
        extraction_input = file_pairs
            .map { name, mzml, tsv -> tuple(tsv, name) }
        feature_data = extract_feature_data(extraction_input)
        
        // Step 3: Generate heatmaps if enabled (join spectrum data with feature data)
        if (params.mmf_generate_heatmap) {
            heatmap_input = processed_mzml
                .map { basename, pkl -> tuple(basename, pkl) }
                .join(feature_data, by: 0)
                .map { basename, pkl, _feature_pkl, feature_json -> 
                    tuple(pkl, basename, feature_json) 
                }
            heatmap_output = create_heatmap_image(heatmap_input)
        }
        
        // Step 4: Collect all feature data files
        all_feature_files = feature_data
            .map { name, pkl, json -> pkl }
            .collect()
        
        // Step 5: Pair features across files (always run, will handle single file)
        paired_output = pair_features(all_feature_files)
        
        // Step 6: Build network graph visualization if enabled (build before report so we can include clustering data)
        clusters_json_file = null
        if (params.mmf_build_network_graph) {
            graph_edges_json = paired_output.map { pkl, json, e_pkl, e_json -> e_json.toAbsolutePath().toString() }
            feature_data_jsons = feature_data.map { name, pkl, json -> json }.collect()
            graph_output = build_network_graph_visualization(graph_edges_json, feature_data_jsons)
            // Output tuple: [svg, graphml, gml, edgelist, analysis, composition, clusters, histogram]
            clusters_json_file = graph_output.map { it[6] }.flatten().first()
        }
        
        // Step 6: Generate report with paired data and cluster data
        if (clusters_json_file) {
            report_input = paired_output
                .combine(clusters_json_file)
                .map { pkl, json, edges_pkl, edges_json, clusters_json ->
                    tuple(pkl, json, edges_pkl, edges_json, clusters_json)
                }
        } else {
            report_input = paired_output
                .map { pkl, json, edges_pkl, edges_json ->
                    tuple(pkl, json, edges_pkl, edges_json, null)
                }
        }
        pairing_report = feature_pairing_report(report_input)
        
        emit:
            feature_data = feature_data
            paired_features = paired_output.map { pkl, json, e_pkl, e_json -> pkl }
            paired_json = paired_output.map { pkl, json, e_pkl, e_json -> json }
            edges_pkl = paired_output.map { pkl, json, e_pkl, e_json -> e_pkl }.flatten()
            edges_json = paired_output.map { pkl, json, e_pkl, e_json -> e_json }.flatten()
            pairing_summary = pairing_report.map { summary, stats, html -> summary }
            pairing_stats = pairing_report.map { summary, stats, html -> stats }
            pairing_html = pairing_report.map { summary, stats, html -> html }
            graph_svg = graph_output ? graph_output.map { it[0] } : Channel.empty()
            graph_graphml = graph_output ? graph_output.map { it[1] } : Channel.empty()
            graph_gml = graph_output ? graph_output.map { it[2] } : Channel.empty()
            graph_edgelist = graph_output ? graph_output.map { it[3] } : Channel.empty()
            graph_analysis = graph_output ? graph_output.map { it[4] } : Channel.empty()
            graph_composition = graph_output ? graph_output.map { it[5] } : Channel.empty()
            clusters_json = graph_output ? graph_output.map { it[6] } : Channel.empty()
            graph_histogram = graph_output ? graph_output.map { it[7] } : Channel.empty()
}
