#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Standalone Nextflow workflow for processing cluster-filtered consensusXML files
 * Uses existing processes from quantify_and_align.nf to generate final reports
 */

// Required Parameters
params.pcc_clusters_consensus_xml = "results/feature_analysis/feature_data_lists/clusters_fraction1.0_filtered.consensusXML"  // Cluster-recovered consensusXML (includes recovered single-vertex clusters)
params.pcc_original_consensus_xml = "results/quantifications/features_with_annotated_identifications/consensus________0.01_fdr.consensusXML"  // Original consensus for reference
params.pcc_feature_tsvs_dir = "results/quantifications/features_with_annotated_identifications"  // Directory containing feature TSV files
params.pcc_outdir = "results/feature_analysis/clusters_final_report"  // Output directory
params.pcc_minlh = 7  // Minimum number of MS1 scans (needed for visualizations)

// Fixed FDR label for cluster processing
params.pcc_fdr_label = "filtered_clusters"

// Validate inputs
if (!file(params.pcc_clusters_consensus_xml).exists()) {
    error("ERROR: Cluster consensusXML file not found: ${params.pcc_clusters_consensus_xml}")
}

if (!file(params.pcc_original_consensus_xml).exists()) {
    error("ERROR: Original consensusXML file not found: ${params.pcc_original_consensus_xml}")
}

// Standalone MAIN Workflow
workflow {
    // Load input files
    clusters_consensus = file(params.pcc_clusters_consensus_xml)
    original_consensus = file(params.pcc_original_consensus_xml)
    
    // Load feature TSV files and collect into list
    feature_tsvs_list = Channel.fromPath(params.pcc_feature_tsvs_dir + "/*.tsv", checkIfExists: true)
        .toList()
    
    // Create tuple: [fdr_label, clusters_consensus, original_consensus, [tsvs]]
    consensus_input = feature_tsvs_list.map { tsvs ->
        tuple(params.pcc_fdr_label, clusters_consensus, original_consensus, tsvs)
    }
    
    // Execute cluster processing workflow
    process_clusters_consensus(consensus_input)
}

// Importable Workflow
workflow process_clusters_consensus {
    take:
        input_data  // tuple(fdr_label, clusters_consensus, original_consensus, [feature_tsvs])
    
    main:
        // Map to extract only fdr_label, clusters_consensus, and feature_tsvs for the process
        // (original_consensus is kept for future reference/comparison)
        input_data.map { fdr, clusters_cons, orig_cons, tsvs ->
            tuple(fdr, clusters_cons, tsvs)
        }.set { consensus_tsvs_input }
        
        // Generate TSV tables from cluster consensusXML and feature data
        generate_feature_ident_intesity_table(consensus_tsvs_input)
        
        // Generate Excel reports from the TSV tables
        generate_xlsx_reports_from_tables(generate_feature_ident_intesity_table.out)
        
        // Generate visualization plots from TSV data
        generate_plots_from_tables(generate_feature_ident_intesity_table.out)
    
    emit:
        tsv_outputs = generate_feature_ident_intesity_table.out
        xlsx_outputs = generate_xlsx_reports_from_tables.out
        plot_outputs = generate_plots_from_tables.out
}


// ============================================================================
// PROCESSES (adapted from quantify_and_align.nf for cluster processing)
// ============================================================================

process generate_feature_ident_intesity_table {
    publishDir "${params.pcc_outdir}", mode:'copy'
    container "luxii/unbequant:latest"
    stageInMode "copy"

    input:
    tuple val(fdr), file(consensus), file(tsvs)

    output:
    tuple val(fdr), file("clusters_raw_quantification_with_identifications.tsv"), file("clusters_quantification_with_identifications_reduced.tsv"), file("clusters_quantification_with_identifications_only_intensities_and_ids.tsv")

    """
    CONCAT_TSVS=""
    for file in ${tsvs}
    do
        CONCAT_TSVS+="\${file},"
    done
    CONCAT_TSVS=\$(echo \$CONCAT_TSVS | rev | cut -c2- | rev)

    PYTHONUNBUFFERED=1 python3 \$(which consensus_and_features_to_tsv.py) \
        -featurexmls_tsvs  \$CONCAT_TSVS \
        -consensus ${consensus} \
        -out_tsv clusters_raw_quantification_with_identifications.tsv \
        -out_tsv_reduced clusters_quantification_with_identifications_reduced.tsv \
        -out_tsv_minimal clusters_quantification_with_identifications_only_intensities_and_ids.tsv
    """
}


process generate_xlsx_reports_from_tables {
    publishDir "${params.pcc_outdir}", mode:'copy'
    container "luxii/unbequant:latest"
    stageInMode "copy"

    input:
    tuple val(fdr), file(full_file), file(reduced_file), file(minimal_file)

    output:
    file("*.xlsx")

    """
    PYTHONUNBUFFERED=1 python3 \$(which write_xlsx_report.py) -i ${full_file} -o ${full_file.baseName}.report.xlsx
    PYTHONUNBUFFERED=1 python3 \$(which write_xlsx_report.py) -i ${reduced_file} -o ${reduced_file.baseName}.report.xlsx
    PYTHONUNBUFFERED=1 python3 \$(which write_xlsx_report.py) -i ${minimal_file} -o ${minimal_file.baseName}.report.xlsx
    """
}


process generate_plots_from_tables {
    publishDir "${params.pcc_outdir}", mode:'copy'
    container "luxii/unbequant:latest"
    stageInMode "copy"

    input:
    tuple val(fdr), file(full_file), file(reduced_file), file(minimal_file)

    output:
    file("*.html")

    """
    PYTHONUNBUFFERED=1 python3 \$(which visualize_consensus_features_infos.py) \
        -minlh_parameter ${params.pcc_minlh} \
        -tsv_file ${full_file}
    """
}
