#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// Map mzML Features - Subworkflow
// ============================================================================
// This workflow maps mzML spectra data with TSV feature information to create
// paired feature datasets across multiple files.

// Configuration parameters
params.main_outdir = "$PWD/results"  // Main output directory (default: PWD/results)
params.mmf_generate_heatmap = false  // Generate heatmap images from mzML files
// Note: Pairing method is now hardcoded to KD-tree based optimization in pair_features.py and pair_features_dedup_edges.py
// Note: best_match_only is now hardcoded to true in pair_features.py and build_network_graph.py
params.mmf_match_cutoff = 0.0  // Minimum match score cutoff (0.0-1.0), discontinued
params.mmf_mz_cutoff = 0.1  // Maximum m/z coordinate difference for filtering (original space)
params.mmf_rt_cutoff = 100  // Maximum RT coordinate difference for filtering (original space)
params.mmf_edges_cutoff = null  // Maximum euclidean distance cutoff for filtering (scaled space)
params.mmf_normalize_coordinates = true  // Enable/disable coordinate normalization (null=auto)
params.mmf_postpair_normalize_coordinates = false  // Rescale axes to balance ranges and add rescaled distance fields (original preserved). Only has an effekt on the distances and the clustering. The Pairing is determined by mmf_normalize_coordinates and for distinct mz and rt cuttoffs disabled by default, since both cutoff points are used to apply a already "scaled" cutoff different in x and y axis
params.mmf_distance_calc_before_scaling = false  // Calculate distance before coordinate scaling
params.mmf_normalize_edge_distances = false  // Normalize edge distances to [0, 1] range (default: disabled)
params.mmf_skip_json_output = false  // Skip JSON serialization for speed (default: disabled - JSON output enabled)
params.mmf_analyze_pep_idents = true  // Perform detailed pep_ident analysis (default: enabled)
params.mmf_disable_multiprocessing = false  // Disable multiprocessing (default: enabled)

//Heatmap parameters
params.mmf_round_up_to = 2  // Number of decimal places to round m/z values (nur für heatmaps)
params.mmf_log_scale = true  // Apply logarithmic scaling to intensities
params.mmf_scale_colors = true  // Scale colors based on min/max intensities
params.mmf_feature_mode = "first_mean"  // Feature center calculation mode: "CoM" (center of mass, all isotopes), "first_mean" (first isotope only), "single_peak" (strongest peak)
params.mmf_feature_center_diagnostic = false  // Generate diagnostic plots
params.mmf_intensity_method = "top3_sum"  // Intensity calculation method
params.mmf_compression_level = null  // HDF5 compression: null=none (fastest), 1=light, 4=balanced, 9=max
params.mmf_row_batch_size = 10  // Rows to process per batch (higher=faster but more RAM)
params.mmf_heatmap_features_only = true  // Save only mapped feature points without heatmap background (creates much smaller file)
params.mmf_heatmap_save_report = true  // Generate HTML report with heatmap and statistics
// Trimming parameters for feature extraction (optional, can help remove noisy edges)
params.mmf_rt_start_trim = 0  // Seconds to trim from the start of retention time (0 = no trimming)
params.mmf_rt_end_trim = 0  // Maximum RT value to keep (keep only features with RT <= this value, 0 = no limit)
// RT Correction Parameters (apply corrections to RT values before pairing)
params.mmf_compare_rt_alignment = false  // Enable RT alignment comparison and generate correction models
params.mmf_rt_alignment_method = "aligntree"  // Fitting method for RT alignment: 'polynomial', 'spline', 'rbf', 'piecewise', 'aligntree', or 'loess' (default: 'polynomial')
params.mmf_rt_correction_mode = 'pre-kdtree'  // RT correction mode: 'none', 'pre-kdtree', or 'per-distance' (default: 'pre-kdtree' to apply corrections before KD-tree matching)
params.mmf_rt_correction_json = null  // Path to JSON file with RT correction models (auto-generated from compare_rt_alignment if enabled)
params.mmf_pairing_rt_source = 'y_center'  // RT source for pairing: 'y_center', 'y_center_geo', 'rt_start', 'rt_end' (default: 'y_center' - the calculated RT center used for KD-tree matching)
params.mmf_mz_rt_weight_ratio = 0.005338  // Weight ratio for RT relative to m/z in distance calculations (default: 1.0 for equal weighting). Values > 1.0 emphasize RT, < 1.0 emphasize m/z
params.mmf_rt_alignment_outlier_threshold = 0.95  // Percentile threshold for outlier exclusion (0.0-1.0). 0.95 = keep central 95% of data, exclude outer 5% (default: 0.95). Set to 0 to disable filtering.
params.mmf_rt_alignment_one_to_one_only = true  // Only include pep_idents that appear exactly once in both files
params.mmf_rt_alignment_mz_source = 'start'  // M/Z source for alignment: 'geo_center', 'center', or 'start'
params.mmf_rt_alignment_rt_source = 'start'  // RT source for alignment: 'geo_center', 'center', or 'start'
params.mmf_rt_alignment_polynomial_degree = 8  // Polynomial degree for fitting RT correction curves (default: 8)
params.mmf_rt_alignment_spline_degree = 3  // Spline degree for cubic splines (default: 3, used with spline method), The spline degree controls the order of polynomial for each spline piece: 1 = Linear (straight line segments) - very rigid, 2 = Quadratic - more flexible, 3 = Cubic (default) - smooth, flexible curves, 4+ = Higher order - increasingly wiggly
params.mmf_rt_alignment_spline_smoothing = 50  // Smoothing factor for splines (default: 100, used with spline method)
params.mmf_rt_alignment_loess_fraction = 0.04  // LOESS smoothing fraction - proportion of data for local regression (default: 0.3, range 0.01-1.0), Lower values (e.g., 0.1) = tighter fit to data (less smooth), Higher values (e.g., 0.5) = more smoothing
params.mmf_rt_alignment_loess_iterations = 100  // LOESS iterations for robustness (default: 3)
params.mmf_rt_alignment_decimal_points = 40  // Decimal places for polynomial coefficients in equations (default: 35)
params.mmf_trafoxmls_dir = "results/quantifications/features_with_annotated_identifications"  // Directory containing trafoXML files for AlignTree method (OpenMS RT transformations, optional)

//Pairing Parameters
params.mmf_generate_html_report = true  // Generate HTML pairing report

// Network Graph Parameters (build_network_graph)
params.mmf_build_network_graph = true  // Build and analyze network graph composition
params.mmf_graph_edge_cutoff = null  // Maximum euclidean distance for graph edges (carefull when cutoff is already set in pairing step)
params.mmf_graph_mz_cutoff = null  // Maximum m/z difference for graph (overrides edge_cutoff if set) - careful when cutoff is already set in pairing step
params.mmf_graph_rt_cutoff = null  // Maximum RT difference for graph (overrides edge_cutoff if set) - careful when cutoff is already set in pairing step
params.mmf_graph_skip_analysis = false  // Skip detailed graph analysis
params.mmf_graph_feature_data_jsons = []  // List of feature data JSON files for coordinate-based layout
// Clustering Parameters
params.mmf_enable_clustering = true  // Enable community detection clustering for HTML report
params.mmf_clustering_method = 'leiden'  // Clustering method: louvain, walktrap, label_propagation, edge_betweenness, leiden (default: leiden)
params.mmf_clustering_use_weights = false  // Use edge distances as weights in clustering
params.mmf_clustering_weight_mode = 'inverse'  // Weight mode: inverse (1/(1+distance)) or distance
params.mmf_clustering_resolution_parameter = 0.98  // Resolution parameter for leiden clustering (default: 0.5, lower=more clusters, higher=fewer clusters)
params.mmf_clustering_objective_function = 'CPM'  // Objective function for leiden clustering (default: CPM, options: CPM or modularity)
params.mmf_output_clusters_json = true  // Export as JSON file with cluster information (optional)
// Resolution Optimization Parameters
params.mmf_auto_select_resolution = true  // Automatically select optimal resolution parameter
params.mmf_resolution_optimization_method = 'gradient-descent'  // Optimization method: 'gradient-descent' (adaptive, default) or 'grid-search' (exhaustive)
params.mmf_resolution_initial = 1.4  // Initial resolution for gradient descent (default: 0.1)
params.mmf_resolution_min = 0.79  // Minimum resolution for grid search (default: 0.01)
params.mmf_resolution_max = 0.79 // Maximum resolution for grid search (default: 2.0)
  // Maximum resolution for grid search (default: 2.0)
params.mmf_resolution_num_points = 1  // Number of points for grid search (default: 20)
params.mmf_resolution_metric = 'avg_cluster_size'  // Optimization metric: 'combined' (cluster size + delete penalty), 'modularity', or 'avg_cluster_size' (maximize overall average cluster size)
params.mmf_resolution_lambda = 100  // Lambda parameter controlling weight of vertices_deleted penalty (default: 0.1, only for 'combined' metric)
params.mmf_resolution_max_iterations = 5  // Maximum number of iterations for gradient descent (default: 15)
params.mmf_resolution_tolerance = 1e-2  // Convergence tolerance for resolution parameter (default: 1e-5, lower=more precise but slower)
// Filtering Parameters
params.mmf_filter_duplicate_file_vertices = true  // Enable filtering to remove duplicate file vertices
params.mmf_filter_method = 'modularity'  // Filter method: simplified-modularity or modularity
params.mmf_mega_complexity_threshold = 5000  // Complexity threshold for heuristic vs exhaustive duplicate filtering (default: 5000, set <=0 to disable heuristic, increase for slight quality increase at significant cost of runtime and resource usage)
params.mmf_save_deleted_vertices = true  // Enable intelligent recovery of deleted vertices with weighted cluster selection
// ConsensusXML Conversion Parameters
params.mmf_convert_clusters_to_consensusxml = true  // Convert recovered clusters (after recovery) to OpenMS ConsensusXML format
// Reproducibility
params.mmf_random_seed = 42  // Random seed for reproducible results (default: null = non-deterministic)
// Output directories
params.mmf_outdir = "${params.main_outdir ?: "$PWD/results"}/feature_analysis"
params.mmf_heatmap_dir = "${params.mmf_outdir}/heatmaps"
params.mmf_features_dir = "${params.mmf_outdir}/feature_data_lists"
params.mmf_analysis_dir = "${params.mmf_outdir}/feature_analysis_plots"
params.mmf_comparison_dir = "${params.mmf_outdir}/rt_alignment_comparison"
// Cluster Consensus Processing Parameters
params.mmf_process_clusters_consensus = true  // Process filtered/clusters consensusXML file
params.mmf_clusters_consensus_file = null  // Path to clusters consensusXML file
params.mmf_original_consensus_file = null  // Path to original consensusXML for comparison (optional)
params.mmf_clusters_outdir = "${params.mmf_outdir}/clusters_analysis"  // Output directory for cluster analysis

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
//   - extract_feature_data: 40G
//   - pair_features: 40G
//   - generate_pairing_report: 40G
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
    label "map_mzml_features"
    cache 'lenient'

    input:
    file(mzml)

    output:
    tuple val(mzml.baseName), path("${mzml.baseName}_spectrum_data.pkl")

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
    label "map_mzml_features_heatmap"
    cache 'lenient'
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
    label "map_mzml_features"
    cache 'lenient'

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
    label "map_mzml_features_pairing"
    cache 'lenient'
    publishDir "${params.mmf_features_dir}", mode: 'symlink'

    input:
    file(feature_files)
    val(rt_correction_json)

    output:
    tuple path("paired_features_edges.pkl"), path("paired_features_edges.json", optional: true),
          path("paired_features_metadata.json"), path("paired_features*_ident_lookup.json", optional: true),
          path("paired_features_unpaired.json", optional: true)

    script:
    def mz_cutoff_arg = params.mmf_mz_cutoff != null ? "--mz_cutoff ${params.mmf_mz_cutoff}" : ""
    def rt_cutoff_arg = params.mmf_rt_cutoff != null ? "--rt_cutoff ${params.mmf_rt_cutoff}" : ""
    def edges_cutoff_arg = params.mmf_edges_cutoff != null ? "--edges_cutoff ${params.mmf_edges_cutoff}" : ""
    def normalize_coords_arg = params.mmf_normalize_coordinates != null ? "--normalize_coordinates ${params.mmf_normalize_coordinates}" : ""
    def postpair_normalize_coords_flag = params.mmf_postpair_normalize_coordinates ? "--postpair_normalize_coordinates" : ""
    def distance_before_scaling_flag = params.mmf_distance_calc_before_scaling ? "--distance_calc_before_scaling" : ""
    def normalize_edge_distances_flag = params.mmf_normalize_edge_distances ? "--normalize_edge_distances" : ""
    def skip_json_output_flag = params.mmf_skip_json_output ? "--skip_json_output" : ""
    def analyze_pep_idents_flag = params.mmf_analyze_pep_idents ? "--analyze_pep_idents" : ""
    def disable_multiprocessing_flag = params.mmf_disable_multiprocessing ? "--disable_multiprocessing" : ""
    def mz_rt_weight_ratio_arg = params.mmf_mz_rt_weight_ratio != 1.0 ? "--mz_rt_weight_ratio ${params.mmf_mz_rt_weight_ratio}" : ""
    
    // Convert trafoXML directory to absolute path (handle both relative and absolute paths)
    // Always pass it to the script, even if not using aligntree (script will ignore it)
    def trafoxmls_dir_abs = null
    def trafoxmls_dir_arg = ""
    if (params.mmf_trafoxmls_dir) {
        def trafoxmls_path = file(params.mmf_trafoxmls_dir)
        if (trafoxmls_path.isAbsolute()) {
            trafoxmls_dir_abs = trafoxmls_path.toString()
        } else {
            trafoxmls_dir_abs = file("${workflow.projectDir}/${params.mmf_trafoxmls_dir}").toAbsolutePath().toString()
        }
        trafoxmls_dir_arg = "--trafoxmls_dir ${trafoxmls_dir_abs}"
    }
    
    // RT correction handling:
    // - Disable ALL RT alignment (including aligntree) if mmf_compare_rt_alignment = false
    // - For aligntree: trafoXML transformations are applied before pairing (trafoxmls_dir is all we need)
    //   No JSON-based corrections needed
    // - For other methods: JSON corrections from compare_rt_alignment are applied
    def rt_method = params.mmf_compare_rt_alignment ? params.mmf_rt_alignment_method : 'none'
    
    def rt_correction_args = ""
    if (rt_method != 'aligntree' && 
        params.mmf_compare_rt_alignment && params.mmf_rt_correction_mode != 'none' && 
        rt_correction_json && rt_correction_json.toString().trim()) {
        // Other methods: use JSON corrections from compare_rt_alignment
        rt_correction_args = "--rt_correction_json ${rt_correction_json} --rt_correction_mode ${params.mmf_rt_correction_mode} --rt_source ${params.mmf_pairing_rt_source}"
    }
    
    // Build the --input_files argument list from all files in the input
    def input_files_arg = feature_files.collect { "${it}" }.join(" ")
    
    """
    echo "DEBUG: RT Correction Configuration:"
    echo "  - mmf_compare_rt_alignment: ${params.mmf_compare_rt_alignment}"
    echo "  - mmf_rt_correction_mode: ${params.mmf_rt_correction_mode}"
    echo "  - RT JSON value: '${rt_correction_json}'"
    echo "  - RT Correction Args: '${rt_correction_args}'"
    echo "  - Configured rt_alignment_method: '${params.mmf_rt_alignment_method}'"
    echo "  - Actual rt_alignment_method passed to script: '${rt_method}'"
    
    echo "DEBUG: Input feature files to process:"
    for f in ${input_files_arg}; do
        [ -f "\$f" ] && ls -lh "\$f" || echo "WARNING: File not found: \$f"
    done | head -3
    
    echo "DEBUG: TrafoXML Argument:"
    echo "  - trafoxmls_dir_arg: '${trafoxmls_dir_arg}'"
    echo "  - params.mmf_trafoxmls_dir: '${params.mmf_trafoxmls_dir}'"
    echo "  - params.mmf_rt_alignment_method: '${params.mmf_rt_alignment_method}'"
    echo "  - Actual method being used: '${rt_method}'"
    
    ${workflow.projectDir}/bin/pair_features_1.py \
        --input_files ${input_files_arg} \
        --output_pkl paired_features.pkl \
        --output_json paired_features.json \
        --output_edges_pkl paired_features_edges.pkl \
        --output_unpaired_json paired_features_unpaired.json \
        --match_cutoff ${params.mmf_match_cutoff} \
        ${mz_cutoff_arg} \
        ${rt_cutoff_arg} \
        ${edges_cutoff_arg} \
        ${normalize_coords_arg} \
        ${postpair_normalize_coords_flag} \
        ${distance_before_scaling_flag} \
        ${normalize_edge_distances_flag} \
        ${skip_json_output_flag} \
        ${analyze_pep_idents_flag} \
        ${disable_multiprocessing_flag} \
        ${mz_rt_weight_ratio_arg} \
        --rt_alignment_method ${rt_method} \
        ${trafoxmls_dir_arg} \
        ${rt_correction_args}
    
    echo "DEBUG: Output files created:"
    ls -lah paired_features_edges.pkl 2>&1 || echo "✗ paired_features_edges.pkl (expected, required)"
    [ -f paired_features_edges.json ] && echo "  ✓ paired_features_edges.json" || echo "  ✗ paired_features_edges.json (optional, skipped)"
    [ -f paired_features_ident_lookup.json ] && echo "  ✓ paired_features_ident_lookup.json" || echo "  ✗ paired_features_ident_lookup.json (optional, skipped)"
    [ -f paired_features_unpaired.json ] && echo "  ✓ paired_features_unpaired.json" || echo "  ✗ paired_features_unpaired.json (optional, skipped)"
    
    echo "DEBUG: Component building removed - paired_features.pkl/json no longer generated (intentional optimization)"
    
    """
}

process feature_pairing_report {
    tag "pairing_report"
    label "map_mzml_features"
    cache false
    publishDir "${params.mmf_features_dir}", mode: 'symlink'

    input:
    tuple val(paired_json), path(edges_pkl), path(metadata_json), path(edges_json), val(clusters_json), val(filtered_clusters_json), val(deleted_vertices_json), val(resolution_iterations_json), val(graph_composition_json), val(ident_lookup_json)

    output:
    tuple path("pairing_summary.txt"), path("pairing_statistics.csv"), 
          path("pairing_report.html", optional: true)

    script:
    def html_arg = params.mmf_generate_html_report ? "--output_html pairing_report.html" : ""
    def clusters_arg = clusters_json ? "--clusters_json ${clusters_json}" : ""
    def filtered_clusters_arg = filtered_clusters_json ? "--filtered_clusters_json ${filtered_clusters_json}" : ""
    def deleted_vertices_arg = deleted_vertices_json ? "--deleted_vertices_json ${deleted_vertices_json}" : ""
    def resolution_iterations_arg = resolution_iterations_json ? "--resolution_iterations_json ${resolution_iterations_json}" : ""
    def graph_composition_arg = graph_composition_json ? "--graph_composition_json ${graph_composition_json}" : ""
    def ident_lookup_arg = ident_lookup_json ? "--ident_lookup_json ${ident_lookup_json}" : ""
    
    // Pairing parameters
    def mz_cutoff_arg = params.mmf_mz_cutoff != null ? "--mz_cutoff ${params.mmf_mz_cutoff}" : ""
    def rt_cutoff_arg = params.mmf_rt_cutoff != null ? "--rt_cutoff ${params.mmf_rt_cutoff}" : ""
    def edges_cutoff_arg = params.mmf_edges_cutoff != null ? "--edges_cutoff ${params.mmf_edges_cutoff}" : ""
    
    // Trimming parameters
    def rt_start_trim_arg = params.mmf_rt_start_trim ? "--rt_start_trim ${params.mmf_rt_start_trim}" : ""
    def rt_end_trim_arg = params.mmf_rt_end_trim ? "--rt_end_trim ${params.mmf_rt_end_trim}" : ""
    // RT Correction parameters
    def compare_rt_alignment_arg = params.mmf_compare_rt_alignment ? "--compare_rt_alignment" : ""
    def rt_correction_mode_arg = params.mmf_rt_correction_mode ? "--rt_correction_mode ${params.mmf_rt_correction_mode}" : ""
    def rt_alignment_method_arg = "--rt_alignment_method ${params.mmf_rt_alignment_method}"
    def mz_rt_weight_ratio_arg = params.mmf_mz_rt_weight_ratio != 1.0 ? "--mz_rt_weight_ratio ${params.mmf_mz_rt_weight_ratio}" : ""
    def rt_alignment_polynomial_degree_arg = "--rt_alignment_polynomial_degree ${params.mmf_rt_alignment_polynomial_degree}"
    def rt_alignment_spline_degree_arg = "--rt_alignment_spline_degree ${params.mmf_rt_alignment_spline_degree}"
    def rt_alignment_spline_smoothing_arg = "--rt_alignment_spline_smoothing ${params.mmf_rt_alignment_spline_smoothing}"
    def rt_alignment_loess_fraction_arg = "--rt_alignment_loess_fraction ${params.mmf_rt_alignment_loess_fraction}"
    def rt_alignment_outlier_threshold_arg =  "--rt_alignment_outlier_threshold ${params.mmf_rt_alignment_outlier_threshold}"
    def rt_alignment_loess_iterations_arg = "--rt_alignment_loess_iterations ${params.mmf_rt_alignment_loess_iterations}"
    def rt_alignment_decimal_points_arg = "--rt_alignment_decimal_points ${params.mmf_rt_alignment_decimal_points}"
    // Feature extraction parameters
    def feature_mode_arg = params.mmf_feature_mode ? "--feature_mode ${params.mmf_feature_mode}" : ""
    // Network graph parameters
    def graph_edge_cutoff_arg = params.mmf_graph_edge_cutoff != null ? "--graph_edge_cutoff ${params.mmf_graph_edge_cutoff}" : ""
    def graph_mz_cutoff_arg = params.mmf_graph_mz_cutoff != null ? "--graph_mz_cutoff ${params.mmf_graph_mz_cutoff}" : ""
    def graph_rt_cutoff_arg = params.mmf_graph_rt_cutoff != null ? "--graph_rt_cutoff ${params.mmf_graph_rt_cutoff}" : ""
    
    // Clustering parameters
    def enable_clustering_arg = params.mmf_enable_clustering ? "--enable_clustering" : ""
    def clustering_method_arg = params.mmf_enable_clustering ? "--clustering_method ${params.mmf_clustering_method}" : ""
    def clustering_use_weights_arg = params.mmf_clustering_use_weights ? "--clustering_use_weights" : ""
    def clustering_weight_mode_arg = params.mmf_clustering_weight_mode ? "--clustering_weight_mode ${params.mmf_clustering_weight_mode}" : ""
    def clustering_resolution_parameter_arg = params.mmf_clustering_resolution_parameter != null ? "--clustering_resolution_parameter ${params.mmf_clustering_resolution_parameter}" : ""
    def clustering_objective_function_arg = params.mmf_clustering_objective_function ? "--clustering_objective_function ${params.mmf_clustering_objective_function}" : ""
    def mega_complexity_threshold_arg = "--mega_complexity_threshold ${params.mmf_mega_complexity_threshold}"
    
    // Post-pairing normalization
    def postpair_normalize_coords_arg = params.mmf_postpair_normalize_coordinates ? "--postpair_normalize_coordinates" : ""
    def input_files_arg = "--input_files ${edges_json}"
    
    """
    ${workflow.projectDir}/bin/generate_pairing_report_v2.py \
        --edges_json ${edges_json} \
        --metadata_json ${metadata_json} \
        --output_summary pairing_summary.txt \
        --output_stats pairing_statistics.csv \
        ${html_arg} \
        ${clusters_arg} \
        ${filtered_clusters_arg} \
        ${deleted_vertices_arg} \
        ${resolution_iterations_arg} \
        ${graph_composition_arg} \
        ${ident_lookup_arg} \
        ${mz_cutoff_arg} \
        ${rt_cutoff_arg} \
        ${edges_cutoff_arg} \
        ${rt_start_trim_arg} \
        ${rt_end_trim_arg} \
        ${compare_rt_alignment_arg} \
        ${rt_correction_mode_arg} \
        ${rt_alignment_method_arg} \
        ${mz_rt_weight_ratio_arg} \
        ${rt_alignment_polynomial_degree_arg} \
        ${rt_alignment_spline_degree_arg} \
        ${rt_alignment_spline_smoothing_arg} \
        ${rt_alignment_loess_fraction_arg} \
        ${rt_alignment_outlier_threshold_arg} \
        ${rt_alignment_loess_iterations_arg} \
        ${rt_alignment_decimal_points_arg} \
        ${feature_mode_arg} \
        ${graph_edge_cutoff_arg} \
        ${graph_mz_cutoff_arg} \
        ${graph_rt_cutoff_arg} \
        ${enable_clustering_arg} \
        ${clustering_method_arg} \
        ${clustering_use_weights_arg} \
        ${clustering_weight_mode_arg} \
        ${clustering_resolution_parameter_arg} \
        ${clustering_objective_function_arg} \
        ${mega_complexity_threshold_arg} \
        ${postpair_normalize_coords_arg} \
        ${input_files_arg}
    """
}
process convert_clusters_to_consensusxml {
    """
    Convert filtered cluster data from JSON to OpenMS ConsensusXML format.
    Creates a ConsensusXML file compatible with OpenMS tools.
    Also includes deleted vertices (as new clusters in original components)
    and unpaired vertices (as new components).
    """
    tag "clusters_to_consensusxml"
    label "map_mzml_features"
    cache 'lenient'
    publishDir "${params.mmf_features_dir}", mode: 'copy'

    input:
    path clusters_json
    path deleted_vertices_json, stageAs: 'deleted_vertices.json'
    path unpaired_vertices_json, stageAs: 'unpaired_vertices.json'

    output:
    path("*.consensusXML")

    script:
    def deleted_arg = deleted_vertices_json ? "--deleted_vertices_json ${deleted_vertices_json}" : ""
    def unpaired_arg = unpaired_vertices_json ? "--unpaired_vertices_json ${unpaired_vertices_json}" : ""
    """
    ${workflow.projectDir}/bin/clusters_json_to_consensusxml.py \
        --input_json ${clusters_json} \
        --output_xml ${clusters_json.baseName}.consensusXML \
        ${deleted_arg} \
        ${unpaired_arg} \
        --mmf_mz_cutoff ${params.mmf_mz_cutoff} \
        --mmf_rt_cutoff ${params.mmf_rt_cutoff} \
        ${params.mmf_edges_cutoff != null ? "--mmf_edges_cutoff ${params.mmf_edges_cutoff}" : ""} \
        ${params.mmf_normalize_coordinates != null ? "--mmf_normalize_coordinates ${params.mmf_normalize_coordinates}" : ""} \
        --mmf_postpair_normalize_coordinates ${params.mmf_postpair_normalize_coordinates} \
        --mmf_distance_calc_before_scaling ${params.mmf_distance_calc_before_scaling} \
        --mmf_normalize_edge_distances ${params.mmf_normalize_edge_distances} \
        --mmf_round_up_to ${params.mmf_round_up_to} \
        --mmf_log_scale ${params.mmf_log_scale} \
        --mmf_scale_colors ${params.mmf_scale_colors} \
        --mmf_feature_mode ${params.mmf_feature_mode} \
        --mmf_build_network_graph ${params.mmf_build_network_graph} \
        ${params.mmf_graph_edge_cutoff != null ? "--mmf_graph_edge_cutoff ${params.mmf_graph_edge_cutoff}" : ""} \
        ${params.mmf_graph_mz_cutoff != null ? "--mmf_graph_mz_cutoff ${params.mmf_graph_mz_cutoff}" : ""} \
        ${params.mmf_graph_rt_cutoff != null ? "--mmf_graph_rt_cutoff ${params.mmf_graph_rt_cutoff}" : ""} \
        --mmf_enable_clustering ${params.mmf_enable_clustering} \
        --mmf_clustering_method ${params.mmf_clustering_method} \
        --mmf_clustering_use_weights ${params.mmf_clustering_use_weights} \
        --mmf_clustering_weight_mode ${params.mmf_clustering_weight_mode} \
        --mmf_clustering_resolution_parameter ${params.mmf_clustering_resolution_parameter} \
        --mmf_clustering_objective_function ${params.mmf_clustering_objective_function} \
        --mmf_auto_select_resolution ${params.mmf_auto_select_resolution} \
        --mmf_resolution_optimization_method ${params.mmf_resolution_optimization_method} \
        --mmf_resolution_min ${params.mmf_resolution_min} \
        --mmf_resolution_max ${params.mmf_resolution_max} \
        --mmf_resolution_num_points ${params.mmf_resolution_num_points} \
        --mmf_resolution_metric ${params.mmf_resolution_metric} \
        --mmf_resolution_lambda ${params.mmf_resolution_lambda} \
        --mmf_resolution_max_iterations ${params.mmf_resolution_max_iterations} \
        --mmf_filter_duplicate_file_vertices ${params.mmf_filter_duplicate_file_vertices} \
        --mmf_filter_method ${params.mmf_filter_method} \
        --mmf_save_deleted_vertices ${params.mmf_save_deleted_vertices} \
        --mmf_random_seed ${params.mmf_random_seed}
    """
}

process build_network_graph_visualization {
    tag "network_graph"
    label "map_mzml_features_graph"
    cache 'lenient'
    publishDir "${params.mmf_features_dir}", mode: 'copy'

    input:
    val(edges_json)
    path feature_data_jsons, stageAs: 'feature_*.json'
    val(ident_lookup_json)
    val(unpaired_json)

    output:
    tuple path("graph_analysis*.json", optional: true),
        path("graph_composition*.json", optional: true),
        path("*clusters*.json", optional: true),
        path("*clusters_*filtered.json", optional: true),
        path("deleted_vertices*.json", optional: true),
        path("resolution_iterations*.json", optional: true)

    script:
    def edge_cutoff_arg = params.mmf_graph_edge_cutoff != null ? "--edge_cutoff ${params.mmf_graph_edge_cutoff}" : ""
    def mz_cutoff_arg = params.mmf_graph_mz_cutoff != null ? "--mz_cutoff ${params.mmf_graph_mz_cutoff}" : ""
    def rt_cutoff_arg = params.mmf_graph_rt_cutoff != null ? "--rt_cutoff ${params.mmf_graph_rt_cutoff}" : ""
    def skip_analysis_flag = params.mmf_graph_skip_analysis ? "--skip-analysis" : ""
    def enable_clustering_arg = params.mmf_enable_clustering ? "--enable-clustering" : ""
    def clustering_method_arg = params.mmf_enable_clustering ? "--clustering-method ${params.mmf_clustering_method}" : ""
    def clustering_weights_arg = params.mmf_clustering_use_weights ? "--clustering-use-weights" : "--clustering-no-weights"
    def clustering_weight_mode_arg = params.mmf_clustering_weight_mode ? "--clustering-weight-mode ${params.mmf_clustering_weight_mode}" : ""
    def clustering_resolution_parameter_arg = params.mmf_clustering_resolution_parameter != null ? "--resolution-parameter ${params.mmf_clustering_resolution_parameter}" : ""
    def clustering_objective_function_arg = params.mmf_clustering_objective_function ? "--objective-function ${params.mmf_clustering_objective_function}" : ""
    def output_clusters_arg = params.mmf_output_clusters_json ? "--output-clusters-json clusters.json" : ""
    def output_deleted_vertices_arg = params.mmf_filter_duplicate_file_vertices ? "--output-deleted-vertices deleted_vertices.json" : ""
    def filter_duplicate_vertices_arg = params.mmf_filter_duplicate_file_vertices ? "--filter-duplicate-file-vertices" : ""
    def filter_method_arg = params.mmf_filter_duplicate_file_vertices ? "--filter-method ${params.mmf_filter_method}" : ""
    def save_deleted_vertices_arg = params.mmf_save_deleted_vertices ? "--save-deleted-vertices" : ""
    def mega_complexity_threshold_arg = "--mega_complexity_threshold ${params.mmf_mega_complexity_threshold}"
    def auto_select_resolution_arg = params.mmf_auto_select_resolution && params.mmf_clustering_method == 'leiden' ? "--auto-select-resolution" : ""
    def resolution_optimization_method_arg = params.mmf_auto_select_resolution ? "--resolution-optimization-method ${params.mmf_resolution_optimization_method}" : ""
    def resolution_initial_arg = params.mmf_auto_select_resolution ? "--resolution-initial ${params.mmf_resolution_initial}" : ""
    def resolution_min_arg = params.mmf_auto_select_resolution && params.mmf_resolution_optimization_method == 'grid-search' ? "--resolution-min ${params.mmf_resolution_min}" : ""
    def resolution_max_arg = params.mmf_auto_select_resolution && params.mmf_resolution_optimization_method == 'grid-search' ? "--resolution-max ${params.mmf_resolution_max}" : ""
    def resolution_num_points_arg = params.mmf_auto_select_resolution && params.mmf_resolution_optimization_method == 'grid-search' ? "--resolution-num-points ${params.mmf_resolution_num_points}" : ""
    def resolution_metric_arg = params.mmf_auto_select_resolution ? "--resolution-metric ${params.mmf_resolution_metric}" : ""
    def resolution_lambda_arg = params.mmf_auto_select_resolution ? "--resolution-lambda ${params.mmf_resolution_lambda}" : ""
    def resolution_max_iterations_arg = params.mmf_auto_select_resolution && params.mmf_resolution_optimization_method == 'gradient-descent' ? "--resolution-max-iterations ${params.mmf_resolution_max_iterations}" : ""
    def resolution_tolerance_arg = params.mmf_auto_select_resolution ? "--resolution-tolerance ${params.mmf_resolution_tolerance}" : ""
    def random_seed_arg = params.mmf_random_seed != null ? "--random-seed ${params.mmf_random_seed}" : ""
    def output_resolution_iterations_arg = params.mmf_auto_select_resolution ? "--output-resolution-iterations resolution_iterations.json" : ""
    def ident_lookup_arg = ident_lookup_json ? "--ident-lookup ${ident_lookup_json}" : ""
    def unpaired_json_arg = unpaired_json ? "--unpaired-json ${unpaired_json}" : ""
    def disable_multiprocessing_flag = params.mmf_disable_multiprocessing ? "--disable-multiprocessing" : ""
    
    """
    python3 ${workflow.projectDir}/bin/build_network_graph_multithread_v2.py \
        --input_pkl ${edges_json} \
        --output_analysis graph_analysis.json \
        --output_composition graph_composition.json \
        ${edge_cutoff_arg} \
        ${mz_cutoff_arg} \
        ${rt_cutoff_arg} \
        ${skip_analysis_flag} \
        ${enable_clustering_arg} \
        ${clustering_method_arg} \
        ${clustering_weights_arg} \
        ${clustering_weight_mode_arg} \
        ${clustering_resolution_parameter_arg} \
        ${clustering_objective_function_arg} \
        ${output_clusters_arg} \
        ${output_deleted_vertices_arg} \
        ${filter_duplicate_vertices_arg} \
        ${filter_method_arg} \
        ${save_deleted_vertices_arg} \
        ${mega_complexity_threshold_arg} \
        ${auto_select_resolution_arg} \
        ${resolution_optimization_method_arg} \
        ${resolution_initial_arg} \
        ${resolution_min_arg} \
        ${resolution_max_arg} \
        ${resolution_num_points_arg} \
        ${resolution_metric_arg} \
        ${resolution_lambda_arg} \
        ${resolution_max_iterations_arg} \
        ${resolution_tolerance_arg} \
        ${random_seed_arg} \
        ${output_resolution_iterations_arg} \
        ${ident_lookup_arg} \
        ${unpaired_json_arg} \
        ${disable_multiprocessing_flag}
    
    # Create placeholder files if outputs were skipped/not generated
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
    
    # Filtered clusters JSON placeholder
    if ! ls *clusters_*filtered.json 1> /dev/null 2>&1; then
        echo '{}' > clusters_filtered.json
    fi
    
    # Deleted vertices JSON placeholder
    if ! ls deleted_vertices*.json 1> /dev/null 2>&1; then
        echo '{}' > deleted_vertices.json
    fi
    
    # Resolution iterations JSON placeholder
    if ! ls resolution_iterations*.json 1> /dev/null 2>&1; then
        echo '{}' > resolution_iterations_placeholder.json
    fi
    """
}

process compare_rt_alignment {
    tag "rt_alignment_comparison"
    label "map_mzml_features"
    cache 'lenient'
    publishDir "${params.mmf_outdir}/rt_alignment_comparison", mode: 'copy'

    input:
    path(feature_jsons)

    output:
    tuple path("rt_alignment_report.html"), path("comparison_summary.txt"), 
          path("fitted_corrections.json")

    script:
    // Count input files
    def files_count = feature_jsons instanceof Collection ? feature_jsons.size() : 1
    
    """
    echo "Comparing RT alignment across ${files_count} feature files..."
    
    # Stage files with expected pattern for the script
    i=0
    for f in ${feature_jsons}; do
        ln -s "\$f" features_\$i.json 2>/dev/null || cp "\$f" features_\$i.json
        i=\$((i + 1))
    done
    
    if [ \$(ls -1 features_*.json 2>/dev/null | wc -l) -gt 1 ]; then
        echo "Running RT alignment comparison on \$(ls -1 features_*.json | wc -l) files"
        
        ${workflow.projectDir}/bin/compare_retention_time_alignment.py \
            --input_dir . \
            --input_pattern "features_*.json" \
            --output_html rt_alignment_report.html \
            --output_json fitted_corrections.json \
            --fitting_method ${params.mmf_rt_alignment_method} \
            --outlier_threshold ${params.mmf_rt_alignment_outlier_threshold} \
            --mz_source ${params.mmf_rt_alignment_mz_source} \
            --rt_source ${params.mmf_rt_alignment_rt_source} \
            --polynomial_degree ${params.mmf_rt_alignment_polynomial_degree} \
            --spline_degree ${params.mmf_rt_alignment_spline_degree} \
            --spline_smoothing ${params.mmf_rt_alignment_spline_smoothing} \
            --loess_fraction ${params.mmf_rt_alignment_loess_fraction} \
            --loess_iterations ${params.mmf_rt_alignment_loess_iterations} \
            --decimal_points ${params.mmf_rt_alignment_decimal_points} \
            ${params.mmf_rt_alignment_one_to_one_only ? '--one_to_one_only' : ''}
    else
        echo "Only one or zero files found - skipping RT alignment comparison"
        echo "RT Alignment Comparison: Skipped (Need at least 2 files)" > comparison_summary.txt
    fi
    
    # Ensure files exist
    [ -f comparison_summary.txt ] || echo "Comparison completed" > comparison_summary.txt
    [ -f fitted_corrections.json ] || echo '{"fitting_method":"${params.mmf_rt_alignment_method}","corrections":{},"note":"No RT corrections calculated"}' > fitted_corrections.json
    [ -f rt_alignment_report.html ] || echo '<!DOCTYPE html><html><body><p>RT Alignment: Skipped (Need at least 2 files)</p></body></html>' > rt_alignment_report.html
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
        
        // Step 1: Process mzML files ONLY if heatmap generation is enabled
        if (params.mmf_generate_heatmap) {
            processed_mzml = process_mzml_file(file_pairs.map { name, mzml, tsv -> mzml })
        } else {
            processed_mzml = Channel.empty()
        }
        
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
        
        // Step 4a: Optional RT alignment comparison (generate RT correction models if enabled)
        // Skip comparison for aligntree method since it uses pre-calculated trafoXML transformations
        rt_alignment_output = Channel.empty()
        if (params.mmf_compare_rt_alignment && params.mmf_rt_alignment_method != 'aligntree') {
            feature_jsons_for_comparison = feature_data
                .map { name, pkl, json -> json }
                .collect()
            rt_alignment_output = compare_rt_alignment(feature_jsons_for_comparison)
            rt_json_ch = rt_alignment_output.map { report, summary, json -> json }
        } else {
            // For aligntree or disabled compare_rt_alignment: create placeholder outputs
            if (params.mmf_rt_alignment_method == 'aligntree') {
                rt_json_ch = Channel.value("")
            } else {
                rt_json_ch = Channel.value("")
            }
        }
        
        // Step 5: Pair features across files (always run, will handle single file)
        paired_output = pair_features(all_feature_files, rt_json_ch)
        
        // Step 6: Build network graph visualization if enabled (build before report so we can include clustering data)
        graph_output = Channel.empty()
        consensus_output = Channel.empty()
        if (params.mmf_build_network_graph) {
            graph_edges_json = paired_output.map { e_pkl, e_json, metadata_json, ident_lookup, unpaired -> e_json.toAbsolutePath().toString() }
            feature_data_jsons = feature_data.map { name, pkl, json -> json }.collect()
            ident_lookup_json = paired_output.map { e_pkl, e_json, metadata_json, ident_lookup, unpaired -> ident_lookup ? ident_lookup.toAbsolutePath().toString() : "" }
            unpaired_json_channel = paired_output.map { e_pkl, e_json, metadata_json, ident_lookup, unpaired -> unpaired ? unpaired.toAbsolutePath().toString() : "" }
            graph_output = build_network_graph_visualization(graph_edges_json, feature_data_jsons, ident_lookup_json, unpaired_json_channel)
            // Output tuple: [analysis, composition, clusters, filtered_clusters, deleted_vertices, resolution_iterations]
            
            // Step 7a: Convert recovered clusters to ConsensusXML if enabled
            if (params.mmf_convert_clusters_to_consensusxml) {
                // Extract recovered clusters, deleted vertices, and unpaired vertices from graph and paired outputs
                // Prefer recovered clusters (after recovery) over filtered (pre-recovery)
                // Use multiMap to create separate channels for each process input
                consensusxml_input = graph_output.combine(paired_output).multiMap {
                    analysis, composition, clusters, filtered_file, deleted, resolution_iterations, e_pkl, e_json, metadata_json, ident_lookup, unpaired ->
                    // Extract recovered clusters from the collection
                    // clusters collection contains both filtered and recovered variants
                    def clusters_list = (clusters instanceof Collection) ? clusters.toList() : [clusters]
                    // Filter out nulls from clusters_list
                    clusters_list = clusters_list.findAll { it != null }
                    def recovered_file = clusters_list.find { it.toString().contains('_recovered') }
                    // Fallback to first clusters file if no _recovered variant found
                    def clusters_final = recovered_file ?: (clusters_list ? clusters_list.first() : file("${workflow.projectDir}/bin/.empty_clusters.json"))
                    
                    clusters: clusters_final
                    deleted: deleted ? ((deleted instanceof Collection) ? deleted.first() : deleted) : file("${workflow.projectDir}/bin/.empty_deleted.json")
                    unpaired: unpaired ? ((unpaired instanceof Collection) ? unpaired.first() : unpaired) : file("${workflow.projectDir}/bin/.empty_unpaired.json")
                }
                
                consensus_output = convert_clusters_to_consensusxml(consensusxml_input.clusters, consensusxml_input.deleted, consensusxml_input.unpaired)
            } else {
                consensus_output = Channel.empty()
            }

            // Step 7: Generate report with paired data and cluster data (graph enabled)
            report_input = paired_output
                .combine(
                    graph_output.map { 
                        analysis, composition, clusters, filtered, deleted, resolution_iterations ->
                        // Handle file collections from wildcard patterns
                        // Prefer recovered clusters file if it exists (contains file linkage scores)
                        def c_list = (clusters instanceof Collection) ? clusters.toList() : [clusters]
                        def c_file = c_list.find { it.toString().contains('_recovered') } ?: c_list.first()
                        
                        def f_file = (filtered instanceof Collection) ? filtered.first() : filtered
                        def d_file = (deleted instanceof Collection) ? deleted.first() : deleted
                        def r_file = (resolution_iterations instanceof Collection) ? resolution_iterations.first() : resolution_iterations
                        // Prefer the base graph_composition.json over fraction-specific variants
                        def comp_list = (composition instanceof Collection) ? composition.toList() : [composition]
                        def comp_file = comp_list.find { it != null && !it.toString().contains('_fraction') } ?: comp_list.find { it != null }
                        
                        tuple(
                            c_file?.toString() ?: null,
                            f_file?.toString() ?: null,
                            d_file?.toString() ?: null,
                            r_file?.toString() ?: null,
                            comp_file?.toString() ?: null
                        )
                    }
                )
                .map { edges_pkl, edges_json, metadata_json, ident_lookup, unpaired, clusters_json, filtered_clusters, deleted_verts, resolution_iters, composition_json ->
                    // paired_json is now optional/absent - pass null
                    tuple(null, edges_pkl, metadata_json, edges_json, clusters_json, filtered_clusters, deleted_verts, resolution_iters, composition_json, ident_lookup)
                }
        } else {
            // Step 7: Generate report with paired data only (graph disabled)
            report_input = paired_output
                .map { edges_pkl, edges_json, metadata_json, ident_lookup, unpaired ->
                    // paired_json is now optional/absent - pass null for clustering data
                    tuple(null, edges_pkl, metadata_json, edges_json, null, null, null, null, null, ident_lookup)
                }
        }
        pairing_report = feature_pairing_report(report_input)
        
        emit:
            feature_data = feature_data
            edges_pkl = paired_output.map { e_pkl, e_json, metadata_json, ident_lookup, unpaired -> e_pkl }.flatten()
            edges_json = paired_output.map { e_pkl, e_json, metadata_json, ident_lookup, unpaired -> e_json }.flatten()
            pairing_summary = pairing_report.map { summary, stats, html -> summary }
            pairing_stats = pairing_report.map { summary, stats, html -> stats }
            pairing_html = pairing_report.map { summary, stats, html -> html }
            rt_alignment_report = rt_alignment_output.map { report, summary, json -> report }.ifEmpty(Channel.empty())
            rt_alignment_summary = rt_alignment_output.map { report, summary, json -> summary }.ifEmpty(Channel.empty())
            rt_alignment_json = rt_alignment_output.map { report, summary, json -> json }.ifEmpty(Channel.empty())
            graph_analysis = graph_output.map { it[0] }.ifEmpty(Channel.empty())
            graph_composition = graph_output.map { it[1] }.ifEmpty(Channel.empty())
            clusters_json = graph_output.map { it[2] }.ifEmpty(Channel.empty())
            filtered_clusters_json = graph_output.map { it[3] }.ifEmpty(Channel.empty())
            deleted_vertices_json = graph_output.map { it[4] }.ifEmpty(Channel.empty())
            resolution_iterations_json = graph_output.map { it[5] }.ifEmpty(Channel.empty())
            consensus_xml = consensus_output.ifEmpty(Channel.empty())
}


// ============================================================================
// Process: Process Clusters ConsensusXML
// ============================================================================
process process_clusters_consensusxml {
    publishDir "${params.mmf_clusters_outdir}", mode: 'copy'
    label "map_mzml_features"
    time '1h'
    cache 'lenient'

    input:
    path clusters_consensus
    val feature_tsvs_dir
    path original_consensus

    output:
    tuple path("clusters_full.tsv"), path("clusters_reduced.tsv"), path("clusters_minimal.tsv")
    path("*.xlsx")
    path("*.html")
    path("summary.json")

    script:
    def original_arg = original_consensus.name != 'NO_FILE' ? "-original_consensus ${original_consensus}" : ""
    """
    python ${workflow.projectDir}/bin/process_clusters_consensusxml.py \
        -consensus ${clusters_consensus} \
        -feature_tsvs_dir ${feature_tsvs_dir} \
        -out_dir . \
        -verbose \
        ${original_arg}
    """
}

// Standalone workflow for direct execution
workflow {
    // Create channels from input parameters or globs and reference absolute paths
    mzml_pattern = params.mzml_files ?: "${workflow.launchDir}/results/mgfs_mzmls/*.mzML"
    tsv_pattern = params.tsv_files ?: "${workflow.launchDir}/results/quantifications/features_with_annotated_identifications/*.tsv"
    
    // Create channels from file paths
    // Do NOT collect here - let the subworkflow handle channel operations
    mzml_files_ch = Channel.fromPath(mzml_pattern)
    tsv_files_ch = Channel.fromPath(tsv_pattern)
    
    // Call the map_mzml_features workflow
    map_mzml_features(mzml_files_ch, tsv_files_ch)
    
    // Optionally process clusters consensus if enabled and file provided
    if (params.mmf_process_clusters_consensus && params.mmf_clusters_consensus_file) {
        clusters_file = Channel.fromPath(params.mmf_clusters_consensus_file)
        feature_dir = params.qal_outdir ? "${params.qal_outdir}/features_with_annotated_identifications" : "${workflow.launchDir}/results/quantifications/features_with_annotated_identifications"
        original_file = params.mmf_original_consensus_file ? Channel.fromPath(params.mmf_original_consensus_file) : Channel.fromPath("${workflow.launchDir}/results/quantifications/features_with_annotated_identifications/placeholder")
        
        process_clusters_consensusxml(clusters_file, feature_dir, original_file)
    }
}
