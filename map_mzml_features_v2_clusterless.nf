#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// Map mzML Features - Clusterless Iterative Component Building
// ============================================================================
// This workflow uses iterative cycle-based component building without clustering.
//
// Key differences from clustered version:
// - No leiden/louvain clustering
// - No resolution optimization
// - Iterative cycles: pair → filter → accept/delete → repeat
// - Simpler, faster, more predictable results

// ============================================================================
// PARAMETERS
// ============================================================================

// Input files — must be provided via --tsv_files "path/to/*.tsv"
params.tsv_files = null

// Main directories
params.main_outdir = "$PWD/results"
params.mmf_outdir = "${params.main_outdir}/feature_analysis_clusterless"

// Feature extraction parameters
params.mmf_round_up_to = 2  // Decimal places for m/z rounding (only for heatmaps, no impact on calculation)
params.mmf_feature_mode = "first_mean"  // Feature center calculation mode
params.mmf_rt_start_trim = 0  // Seconds to trim from start of RT
params.mmf_rt_end_trim = 0  // Maximum RT value to keep
params.mmf_intensity_method = "top3_sum"  // Intensity calculation method

// Feature pairing parameters
params.mmf_mz_cutoff = null  // Maximum m/z coordinate difference (mutually exclusive with mmf_ppm)
params.mmf_ppm = 10       // m/z cutoff as parts per million: cutoff = feature_mz * ppm / 10000000 (mutually exclusive with mmf_mz_cutoff)
params.mmf_rt_cutoff = 100  // Maximum RT coordinate difference
params.mmf_edges_cutoff = null  // Maximum euclidean distance cutoff (optional)
params.mmf_normalize_coordinates = true  // Enable/disable coordinate normalization (null=auto)
params.mmf_pairing_rt_source = 'y_center'  // RT field used for KD-tree pairing: 'y_center', 'y_center_geo', 'rt_start', 'rt_end'
params.mmf_mz_rt_weight_ratio = 1  // Weight ratio for RT relative to m/z in distance calculations (< 1.0 emphasizes m/z)

// RT alignment parameters
params.mmf_compare_rt_alignment = true  // Enable RT alignment
params.mmf_rt_alignment_method = "loess"  // RT alignment method: 'loess', 'polynomial', 'spline', 'rbf', 'piecewise', 'aligntree'
params.mmf_trafoxmls_dir = "results/quantifications/features_with_annotated_identifications"  // Directory containing trafoXML files (required when mmf_rt_alignment_method = 'aligntree')
params.mmf_rt_correction_mode = 'pre-kdtree'  // When to apply RT corrections
params.mmf_rt_alignment_loess_fraction = 0.04  // LOESS smoothing fraction
params.mmf_rt_alignment_outlier_threshold = 0.95  // Outlier exclusion threshold
params.mmf_rt_alignment_one_to_one_only = true  // Only include pep_idents that appear exactly once
params.mmf_rt_alignment_mz_source = 'start'  // M/Z source for alignment
params.mmf_rt_alignment_rt_source = 'start'  // RT source for alignment
params.mmf_rt_alignment_polynomial_degree = 8  // Polynomial degree for fitting
params.mmf_rt_alignment_spline_degree = 3  // Spline degree for cubic splines
params.mmf_rt_alignment_spline_smoothing = 50  // Smoothing factor for splines
params.mmf_rt_alignment_loess_iterations = 100  // LOESS iterations for robustness
params.mmf_rt_alignment_decimal_points = 40  // Decimal places for polynomial coefficients
params.mmf_comparison_dir = "${params.mmf_outdir}/rt_alignment_comparison"

// Filtering parameters
params.mmf_filter_duplicate_file_vertices = true  // Enable duplicate filtering
params.mmf_filter_method = 'modularity'  // Filter method, either 'simplified-modularity' or 'modularity'
params.mmf_mega_complexity_threshold = 5000  // Complexity threshold for heuristic

// Cycle parameters
params.mmf_max_cycles = 500  // Maximum number of iterative cycles
params.mmf_disable_multiprocessing = true  // Enable parallel processing for filtering

// Output parameters
params.mmf_generate_html_report = true  // Generate HTML report
params.mmf_convert_to_consensusxml = true  // Convert to ConsensusXML format

// Mass calculation parameters (for report pep_ident masses)
params.mmf_nterm_mod = "none"       // N-terminal modification: none | acetyl | biotin | fluorescein
params.mmf_cys_mod   = "carbamidomethyl"       // Cysteine modification:   none | oxidized | carbamidomethyl
params.mmf_cterm_mod = "none"       // C-terminal modification: none | amide (-CONH2 instead of -COOH)
// add oxidation of aminoacid M -> in M[o] in dataset
// ============================================================================
// PROCESSES
// ============================================================================

process extract_feature_data {
    """
    Extract feature data from TSV file.
    """
    tag "${basename}"
    container "luxii/unbequant:latest"
    memory "8G"
    cache 'lenient'

    input:
    tuple path(tsv), val(basename)

    output:
    tuple val(basename), path("${basename}_feature_data.pkl"), path("${basename}_feature_data.json")

    script:
    """
    python3 /workspaces/unbeQuant/bin/extract_feature_data.py \\
        --tsv ${tsv} \\
        --output_pkl ${basename}_feature_data.pkl \\
        --output_json ${basename}_feature_data.json \\
        --round_up_to ${params.mmf_round_up_to} \\
        --feature_mode ${params.mmf_feature_mode} \\
        --generate_diagnostic false \\
        --rt_start_trim ${params.mmf_rt_start_trim} \\
        --rt_end_trim ${params.mmf_rt_end_trim}
    """
}

process compare_retention_time_alignment {
    """
    Compare RT alignment between file pairs using identified peptides.
    Generates correction models for RT alignment.
    """
    tag "RT Alignment Comparison"
    container "luxii/unbequant:latest"
    memory "4G"
    publishDir "${params.mmf_comparison_dir}", mode: 'copy', overwrite: true
    
    input:
    path(feature_jsons)
    
    output:
    tuple path("rt_alignment_report.html"), path("comparison_summary.txt"), 
          path("fitted_corrections.json"), emit: rt_outputs
    path("fitted_corrections.json"), emit: rt_correction_json
    
    script:
    // Count input files
    def files_count = feature_jsons instanceof Collection ? feature_jsons.size() : 1
    // aligntree uses OpenMS TrafoXML-based alignment internally; skip the pre-fitting step
    def is_aligntree = params.mmf_rt_alignment_method == 'aligntree'

    """
    echo "Comparing RT alignment across ${files_count} feature files..."

    if ${is_aligntree}; then
        echo "RT alignment method is 'aligntree' - skipping pre-fitting (handled via TrafoXML internally)"
        echo "RT Alignment Comparison: Skipped (aligntree uses internal TrafoXML corrections)" > comparison_summary.txt
    else
        # Stage files with expected pattern for the script
        i=0
        for f in ${feature_jsons}; do
            ln -s "\$f" features_\$i.json 2>/dev/null || cp "\$f" features_\$i.json
            i=\$((i + 1))
        done

        if [ \$(ls -1 features_*.json 2>/dev/null | wc -l) -gt 1 ]; then
            echo "Running RT alignment comparison on \$(ls -1 features_*.json | wc -l) files"

            ${workflow.projectDir}/bin/compare_retention_time_alignment.py \\
                --input_dir . \\
                --input_pattern "features_*.json" \\
                --output_html rt_alignment_report.html \\
                --output_json fitted_corrections.json \\
                --fitting_method ${params.mmf_rt_alignment_method} \\
                --outlier_threshold ${params.mmf_rt_alignment_outlier_threshold} \\
                --mz_source ${params.mmf_rt_alignment_mz_source} \\
                --rt_source ${params.mmf_rt_alignment_rt_source} \\
                --polynomial_degree ${params.mmf_rt_alignment_polynomial_degree} \\
                --spline_degree ${params.mmf_rt_alignment_spline_degree} \\
                --spline_smoothing ${params.mmf_rt_alignment_spline_smoothing} \\
                --loess_fraction ${params.mmf_rt_alignment_loess_fraction} \\
                --loess_iterations ${params.mmf_rt_alignment_loess_iterations} \\
                --decimal_points ${params.mmf_rt_alignment_decimal_points} \\
                ${params.mmf_rt_alignment_one_to_one_only ? '--one_to_one_only' : ''}
        else
            echo "Only one or zero files found - skipping RT alignment comparison"
            echo "RT Alignment Comparison: Skipped (Need at least 2 files)" > comparison_summary.txt
        fi
    fi
    
    # Ensure files exist
    [ -f comparison_summary.txt ] || echo "Comparison completed" > comparison_summary.txt
    [ -f fitted_corrections.json ] || echo '{\"fitting_method\":\"${params.mmf_rt_alignment_method}\",\"corrections\":{},\"note\":\"No RT corrections calculated\"}' > fitted_corrections.json
    [ -f rt_alignment_report.html ] || echo '<!DOCTYPE html><html><body><p>RT Alignment: Skipped (Need at least 2 files)</p></body></html>' > rt_alignment_report.html
    """
}

process orchestrate_clusterless_cycles {
    """
    Main process: Orchestrates iterative cycle-based component building.
    Replaces separate pair_features + build_network_graph calls.
    """
    tag "Clusterless Cycles"
    container "luxii/unbequant:latest"
    memory "40G"
    publishDir "${params.mmf_outdir}", mode: 'copy', overwrite: true
    
    input:
    path(feature_pkls)
    path(rt_correction_json)
    
    output:
    path("final_components_all_cycles.json"), emit: components_json
    path("cycle_history.json"), emit: cycle_history
    path("edges_cycle*.pkl"), emit: edges_pkls, optional: true
    path("unpaired_cycle*.json"), emit: unpaired_jsons, optional: true
    path("accepted_components_cycle*.json"), emit: accepted_jsons, optional: true
    path("deleted_vertices_cycle*.json"), emit: deleted_jsons, optional: true
    
    script:
    def mz_rt_args = params.mmf_edges_cutoff ?
        "--edges_cutoff ${params.mmf_edges_cutoff}" :
        params.mmf_ppm != null ?
            "--ppm ${params.mmf_ppm} --rt_cutoff ${params.mmf_rt_cutoff}" :
            "--mz_cutoff ${params.mmf_mz_cutoff} --rt_cutoff ${params.mmf_rt_cutoff}"
    
    def filter_args = "--filter_method ${params.mmf_filter_method} --mega_complexity_threshold ${params.mmf_mega_complexity_threshold}"

    def mp_arg = params.mmf_disable_multiprocessing ? "--disable_multiprocessing" : ""

    def normalize_coords_arg = params.mmf_normalize_coordinates != null ? "--normalize_coordinates ${params.mmf_normalize_coordinates}" : ""

    def mz_rt_weight_ratio_arg = params.mmf_mz_rt_weight_ratio != 1.0 ? "--mz_rt_weight_ratio ${params.mmf_mz_rt_weight_ratio}" : ""

    def is_aligntree = params.mmf_rt_alignment_method == 'aligntree'
    // Resolve trafoxmls_dir to an absolute path so it's valid inside the Nextflow work dir
    def trafoxmls_arg = ""
    if (is_aligntree && params.mmf_trafoxmls_dir) {
        def trafoxmls_path = file(params.mmf_trafoxmls_dir)
        def trafoxmls_abs = trafoxmls_path.isAbsolute()
            ? trafoxmls_path.toString()
            : file("${workflow.projectDir}/${params.mmf_trafoxmls_dir}").toAbsolutePath().toString()
        trafoxmls_arg = "--trafoxmls_dir ${trafoxmls_abs}"
    }
    // For aligntree, rt_correction_mode must be 'none' (corrections applied via trafoXML, not JSON)
    def rt_correction_mode_arg = is_aligntree ? "none" : params.mmf_rt_correction_mode

    """
    python3 /workspaces/unbeQuant/bin/orchestrate_clusterless_cycles.py \\
        --input_files ${feature_pkls} \\
        --rt_correction_json ${rt_correction_json} \\
        --rt_correction_mode ${rt_correction_mode_arg} \\
        --rt_alignment_method ${params.mmf_rt_alignment_method} \\
        ${trafoxmls_arg} \\
        --rt_source ${params.mmf_pairing_rt_source} \\
        ${mz_rt_args} \\
        ${filter_args} \\
        ${mp_arg} \\
        ${normalize_coords_arg} \\
        ${mz_rt_weight_ratio_arg} \\
        --max_cycles ${params.mmf_max_cycles} \\
        --output_dir .
    """
}

process generate_clusterless_report {
    """
    Generate HTML and text summary reports for clusterless workflow.
    """
    tag "Generate Report"
    container "luxii/unbequant:latest"
    memory "8G"
    publishDir "${params.mmf_outdir}", mode: 'copy', overwrite: true

    input:
    path(components_json)
    path(cycle_history)
    path(tsv_files)

    output:
    path("clusterless_summary.txt"), emit: summary
    path("clusterless_report.html"), emit: html, optional: true

    script:
    def html_arg = params.mmf_generate_html_report ? "--output_html clusterless_report.html" : ""
    def run_params_json = groovy.json.JsonOutput.toJson([
        rt_alignment_method: params.mmf_rt_alignment_method,
        ppm                : params.mmf_ppm,
        mz_cutoff          : params.mmf_mz_cutoff,
        rt_cutoff          : params.mmf_rt_cutoff,
        mz_rt_weight_ratio : params.mmf_mz_rt_weight_ratio,
        nterm_mod          : params.mmf_nterm_mod,
        cys_mod            : params.mmf_cys_mod,
        cterm_mod          : params.mmf_cterm_mod,
    ])

    """
    echo '${run_params_json}' > run_parameters.json

    python3 /workspaces/unbeQuant/bin/generate_pairing_report_v2_clusterless.py \\
        --components_json ${components_json} \\
        --cycle_history_json ${cycle_history} \\
        --output_summary clusterless_summary.txt \\
        --nterm_mod ${params.mmf_nterm_mod} \\
        --cys_mod ${params.mmf_cys_mod} \\
        --cterm_mod ${params.mmf_cterm_mod} \\
        --parameters_json run_parameters.json \\
        --tsv_files ${tsv_files} \\
        ${html_arg}
    """
}

process convert_to_consensusxml {
    """
    Convert components to OpenMS ConsensusXML format.
    """
    tag "ConsensusXML Conversion"
    container "luxii/unbequant:latest"
    memory "8G"
    publishDir "${params.mmf_outdir}", mode: 'copy', overwrite: true
    
    input:
    path(components_json)
    
    output:
    path("components.consensusXML"), emit: consensusxml
    
    script:
    def mz_cutoff_arg = params.mmf_mz_cutoff != null ? "--mmf_mz_cutoff ${params.mmf_mz_cutoff}" : ""
    """
    python3 /workspaces/unbeQuant/bin/clusters_json_to_consensusxml_clusterless.py \\
        --input_json ${components_json} \\
        --output_xml components.consensusXML \\
        ${mz_cutoff_arg} \\
        --mmf_rt_cutoff ${params.mmf_rt_cutoff} \\
        --mmf_filter_duplicate_file_vertices ${params.mmf_filter_duplicate_file_vertices} \\
        --mmf_filter_method ${params.mmf_filter_method}
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow map_mzml_features_clusterless {
    take:
    tsv_files  // Channel of TSV files with feature data
    
    main:
    // Step 0: Extract feature data from TSV files
    tsv_with_basename = tsv_files.map { file -> 
        tuple(file, file.baseName)
    }
    feature_data = extract_feature_data(tsv_with_basename)
    
    // Collect feature pickle files and feature JSONs for next steps
    feature_pkls = feature_data.map { basename, pkl, json -> pkl }.collect()
    feature_jsons = feature_data.map { basename, pkl, json -> json }.collect()
    
    // Step 1: RT alignment comparison (if enabled)
    if (params.mmf_compare_rt_alignment) {
        rt_correction_result = compare_retention_time_alignment(feature_jsons)
        rt_correction_json = rt_correction_result.rt_correction_json
    } else {
        // Create empty RT correction file
        rt_correction_json = Channel.value("")
    }
    
    // Step 2: Orchestrate iterative clusterless cycles
    orchestrate_results = orchestrate_clusterless_cycles(
        feature_pkls,
        rt_correction_json
    )
    
    // Step 3: Generate reports
    if (params.mmf_generate_html_report) {
        generate_clusterless_report(
            orchestrate_results.components_json,
            orchestrate_results.cycle_history,
            tsv_files.collect()
        )
    }
    
    // Step 4: Convert to ConsensusXML (if enabled)
    if (params.mmf_convert_to_consensusxml) {
        convert_to_consensusxml(orchestrate_results.components_json)
    }
    
    emit:
    components_json = orchestrate_results.components_json
    cycle_history = orchestrate_results.cycle_history
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {
    if (!params.tsv_files) {
        error "No input files specified. Provide TSV files via --tsv_files \"path/to/*.tsv\""
    }

    tsv_ch = Channel
        .fromPath(params.tsv_files)
        .ifEmpty { error "No TSV files found matching: ${params.tsv_files}" }

    map_mzml_features_clusterless(tsv_ch)
}

// ============================================================================
// WORKFLOW COMPLETION HANDLER
// ============================================================================

workflow.onComplete {
    println ""
    println "=========================================="
    println "Clusterless Component Building Complete!"
    println "=========================================="
    println "Success:      ${workflow.success}"
    println "Exit status:  ${workflow.exitStatus}"
    println "Duration:     ${workflow.duration}"
    println "Output dir:   ${params.mmf_outdir}"
    println "=========================================="
    println ""
}
