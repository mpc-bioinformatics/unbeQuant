#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// Map mzML Features - Subworkflow
// ============================================================================
// This workflow maps mzML spectra data with TSV feature information to create
// paired feature datasets across multiple files.

// Configuration parameters
params.mmf_generate_heatmap = false  // Generate heatmap images from mzML files
params.mmf_optimize_pairing = true  // Use optimized KD-tree based pairing
params.mmf_match_cutoff = 0.98  // Minimum match score cutoff (0.0-1.0)
params.mmf_round_up_to = 2  // Number of decimal places to round m/z values
params.mmf_log_scale = true  // Apply logarithmic scaling to intensities
params.mmf_scale_colors = true  // Scale colors based on min/max intensities
params.mmf_feature_mode = "CoM"  // Options: "rectangle", "CoM" (Center of Mass)
params.mmf_feature_center_diagnostic = true  // Generate diagnostic plots
params.mmf_intensity_method = "top3_sum"  // Intensity calculation method
params.mmf_compression_level = null  // HDF5 compression: null=none (fastest), 1=light, 4=balanced, 9=max
params.mmf_row_batch_size = 10  // Rows to process per batch (higher=faster but more RAM)

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
    """
    tag "${basename}"
    container "luxii/unbequant:latest"
    memory "8G"
    publishDir "${params.mmf_features_dir}", mode: 'symlink'

    input:
    tuple file(tsv), val(basename)

    output:
    tuple val(basename), file("${basename}_feature_data.pkl"), file("${basename}_feature_data.json")

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
    Pair matching features across multiple files.
    """
    tag "feature_pairing"
    container "luxii/unbequant:latest"
    memory "8G"
    publishDir "${params.mmf_features_dir}", mode: 'symlink'

    input:
    file(feature_files)

    output:
    tuple file("paired_features.pkl"), file("paired_features.json")

    script:
    def optimize_flag = params.mmf_optimize_pairing ? "--optimize" : ""
    """
    ${workflow.projectDir}/bin/pair_features.py \
        --input_dir . \
        --output_pkl paired_features.pkl \
        --output_json paired_features.json \
        --match_cutoff ${params.mmf_match_cutoff} \
        ${optimize_flag}
    """
}

process feature_pairing_report {
    """
    Generate summary report of feature pairing results.
    """
    tag "pairing_report"
    container "luxii/unbequant:latest"
    publishDir "${params.mmf_features_dir}", mode: 'symlink'

    input:
    tuple file(paired_pkl), file(paired_json)

    output:
    tuple file("pairing_summary.txt"), file("pairing_statistics.csv")

    shell:
    '''
    !{workflow.projectDir}/bin/generate_pairing_report.py \
        --paired_json !{paired_json} \
        --output_summary pairing_summary.txt \
        --output_stats pairing_statistics.csv
    '''
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
        
        // Step 5: Generate report
        pairing_report = feature_pairing_report(paired_output)
        
        emit:
            feature_data = feature_data
            paired_features = paired_output.map { pkl, json -> pkl }
            paired_json = paired_output.map { pkl, json -> json }
            pairing_summary = pairing_report.map { summary, stats -> summary }
            pairing_stats = pairing_report.map { summary, stats -> stats }
}
