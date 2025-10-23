#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.qal_spectra_files = "$PWD/raws"  // .RAW/.d-files from which the XICS are extracted.
params.qal_mzmls = "$PWD/mzmls"  // Input mzML-files, which are used to generate the features. NOTE: They need to be named like the raw files, but with the .mzML extension (e.g. "sample1.raw" -> "sample1.mzML").
params.qal_idents = "$PWD/raws/tsvs"  // Folder containing Identifications in TSV-format. Default: They need to be called with the suffix: "*qvalue_no_decoys_fdr_0.0[15].tsv". Change the parameter below, if they are called differently. The columns: containing the columns: "charge", "plain_peptide", "used_score", "retention_time", "exp_mass_to_charge", "fasta_id", "fasta_desc" need to be present.
params.qal_idents_blob_filter = "*qvalue_no_decoys_fdr_0.0[15].tsv"  // Blob filter for the "params.qal_idents" folder (identification files). NOTE: "fdr_0.0[15]" needs to be present in the file name, since this workflow uses this to return fdr-specific results.
params.qal_outdir = "$PWD/results"  // Output-Directory of the quantification results (splitted by the fdr thresholds).

// Optional Parameters
params.qal_charge_low = 1  // Minimum charge of a feature to be considered (biosaur2 -cmin parameter)
params.qal_charge_high = 8  // Maximum charge of a feature to be considered (biosaur2 -cmax parameter)
params.qal_ppm_tolerance = 5  // Tolerance for the biosaur2 to use to map identifications on it, as well as to generate features (biosaur2 -htol and -itol parameters). The ppm tolerance should be set to the same value as it has been used for the search engine for identification of MS2 spectra.
params.qal_minlh = 7  // Minimum number of MS1 scans to be considered for a feature. Check out biosaur2 documnentation to set the correct value. (biosaur2 -minlh parameter)
params.qal_additional_biosaur_parameters  = ""  // Additional parameters for biosaur2. Check out the biosaur2 documentation if you want to set something specific.
params.qal_additional_feature_linker_params = ""  // Additional parameters for the OpenMS FeatureLinkerUnlabeled step. E.G.: "-algorithm:distance_RT:max_difference 300.0" in case of not a good alignment. Check out the OpenMS documentation. 
params.qal_mini = 1  // Minimum intensity for biosaur2 to consider a peak for feature finding. (biosaur2 -mini parameter).
params.qal_rt_enlarge_factor = 0.5  // A factor value to enlarge the RT-window for matching MS2 with features. This factor allow to match MS2 spectra to features, allowing a RT-error to the boundaries of a feature. Formula which enlarges the feature in the RT-dimension: "(MaxRT-MinRT)*enlarge_factor".
params.qal_protgraph_was_used = false  // A Flag which is needed for the output to know which parsing mode and which column of "fasta_id" and "fasta_desc" needs to be taken. If you searched prior with ProtGraph, set this to true. 
params.qal_intensity_method = "sum"  // Method to calculate the quantitative value of the individual XICs. Only available options are "top3_sum", "top3_mean", "maximum" and "sum".. See the python script for more details.
params.qal_cutoff = "t0"  // Minimum quantitative value of the XICs to be considered for the quantification (this is the cutoff, which is applied after retrieving the intensity via the method described in qal_intensity_method) Use tX or qX to set the cutoff (see Python script for more details). Default to no intensity cutoff (or cutoff at intensity smalle or equal to 0).

// Include the XIC-Extractor for Bruker and Thermo
PROJECT_DIR = workflow.projectDir
include {retrieve_xics_from_raw_spectra} from PROJECT_DIR + '/include/xic-extractor/main.nf'

// Standalone MAIN Workflow
workflow {
    raw_files = Channel.fromPath(params.qal_spectra_files  + "/*.raw")
    d_files = Channel.fromPath(params.qal_spectra_files  + "/*.d", type: "dir")
    spectra_files = raw_files.concat(d_files)
    mzmls = Channel.fromPath(params.qal_mzmls + "/*.mzML")
    identifications = Channel.fromPath(params.qal_idents + "/" + params.qal_idents_blob_filter)
    // Mapping and grouping identifications with fdrs
    identifications_tuple = identifications
        .map { file -> tuple(
            file.baseName.substring(file.baseName.indexOf("fdr_") + 4, file.baseName.indexOf("fdr_") + 8), 
            file
        ) }.groupTuple()

    // Execute workflow
    quantify_and_align(
        spectra_files,
        mzmls,
        identifications_tuple
    )
}

// Importable Workflow
workflow quantify_and_align {
    take: 
        // Takes raw/.d files, the corresponding converted mzmls and the corresponding identifications 
        // (in:  tuple(fdr, [list of identifications])) and maps them to features with quantitative 
        // values.
        spectra_files
        mzmls
        identifications_tuple
    main: 
        // Create featureXML from mzML
        create_feature_xml(mzmls)

        // Generate file_identifier and match features with identifications (on multiple fdrs)
        spectra_files_tuple = spectra_files.map { file -> tuple(file.baseName, file) }
        featurexmls_tuple = create_feature_xml.out.map { file -> tuple(file.baseName, file) }
        
        // Get all the single fdrs
        in_identifications_tuple = identifications_tuple.transpose().map { it -> tuple(it[1].baseName.split("_____")[0], it[0], it[1]) }
        
        // Match with identifications using file_identifier (it[0]) and fdr (it[1])
        in_featurexmls_tuple = in_identifications_tuple.map { it -> it[1] } .unique().combine(featurexmls_tuple).map { it -> tuple(it[1], it[0], it[2]) }
        matched_features_with_idents_tuple = in_identifications_tuple.join(in_featurexmls_tuple, by: [0,1])
        mzmls_tuple = mzmls.map { file -> tuple(file.baseName, file) }
        
        // Do the actual matching between features and identifications
        match_feature_with_idents(matched_features_with_idents_tuple)

        // Get all the single fdrs
        in_spectra_files_tuple = in_identifications_tuple.map { it -> it[1] } .unique().combine(spectra_files_tuple).map { it -> tuple(it[1], it[0], it[2]) }
        
        // Match with identifications using file_identifier (it[0]) and fdr (it[1])
        matched_ident_features_with_raws = match_feature_with_idents.out[0].join(in_spectra_files_tuple, by: [0,1])

        // Generate queries for XIC-Extraction from .raw and .d files 
        generate_queries_from_featurexmls(matched_ident_features_with_raws)
        
        // Do the actual XIC Retrieval.
        extract_xics_channel = generate_queries_from_featurexmls.out.map {
            it -> tuple(it[4], it[2])
        }
        retrieve_xics_from_raw_spectra(extract_xics_channel)
        
        // Convert the hdf5 extracted XICs to tsv and add there the additional information from the featureXMLs.
        xics_and_remaining_data = retrieve_xics_from_raw_spectra.out.join(generate_queries_from_featurexmls.out, by: [0])
        extracted_xics_from_hdf5_to_tsv(
            xics_and_remaining_data.map {it -> tuple(it[0], it[2], it[1], it[3]) }
        )

        // Apply map alignment and consensus generation (OpenMS) on the generated (and with MS2 and identification annotated) features.
        identified_features_by_fdr = match_feature_with_idents.out[0].map { it -> tuple(it[1], it[2]) }.groupTuple()
        map_alignment_and_consensus_generation(identified_features_by_fdr)

        // Preperation to generate the final output tables and visualizations.
        fdr_and_feature_tsvs = extracted_xics_from_hdf5_to_tsv.out.map { it -> tuple(it[1], it[2]) }.groupTuple()
        consensus_with_feature_tsvs = map_alignment_and_consensus_generation.out[0].join(fdr_and_feature_tsvs, by: 0)
        
        // RT transformations of each input file
        visualize_RT_transoformations(map_alignment_and_consensus_generation.out[1])

        // Final output table (in tsv and xlsx format) containing all "consensus features".
        generate_feature_ident_intesity_table(consensus_with_feature_tsvs)
        generate_xlsx_reports_from_tables(generate_feature_ident_intesity_table.out)
        
        // Further plots (currently only missingness count per input file)
        generate_plots_from_tables(generate_feature_ident_intesity_table.out)

    emit:
        generate_feature_ident_intesity_table.out[0]  // A tuple containing: [FDR, raw_quantification_full, raw_quantification_reduced, raw_quantification_minimal]
        consensus_with_feature_tsvs  // Consensus with features: [FDR, [consensusXML, [file1_features.tsv, file2_features.tsv, ...]]]
        match_feature_with_idents.out[0]  // [FDR, [file1.featureXML, file2.featureXML, ...]]  (Features with identifications and MS2 annotated)
}


process create_feature_xml {
    stageInMode "copy"
    container "luxii/unbequant:latest"
    memory "18G"

    input:
    file mzml

    output:
    file("${mzml.baseName}.featureXML")

    """
    biosaur2 ${params.qal_additional_biosaur_parameters} -minlh ${params.qal_minlh} -mini ${params.qal_mini} -htol ${params.qal_ppm_tolerance} -itol ${params.qal_ppm_tolerance} -iuse -1 -cmin ${params.qal_charge_low} -cmax ${params.qal_charge_high} -write_hills -write_extra_details ${mzml}

    convert_biosaur2_to_featurexml.py --mzml ${mzml} --feature_tsv ${mzml.baseName}.features.tsv --rt_enlarge_factor ${params.qal_rt_enlarge_factor} --hills_tsv ${mzml.baseName}.hills.tsv --mz_tolerance ${params.qal_ppm_tolerance} --output_featurexml ${mzml.baseName}.featureXML
    """
}


process match_feature_with_idents {
    container "luxii/unbequant:latest"
    publishDir "${params.qal_outdir}/features_with_annotated_identifications", mode:'copy', pattern:'*____with_identifications.featureXML'
    publishDir "${params.qal_outdir}/visualizations___${fdr}/cutoff_plots", mode:'copy', pattern:'*____cutoff_plot.html'

    input:
    tuple val(file_identifier), val(fdr), file(ident_tsv), file(featurexml)

    output:
    tuple val(file_identifier), val(fdr), file("${featurexml.baseName}_____${fdr}_____with_identifications.featureXML")
    file("${featurexml.baseName}_____cutoff_plot.html")

    """
    map_features_with_idents.py -cutoff ${params.qal_cutoff} -featureXML ${featurexml} -use_protgraph ${params.qal_protgraph_was_used} -tsv_file ${ident_tsv} -out_featurexml ${featurexml.baseName}_____${fdr}_____with_identifications.featureXML -out_plot_cutoff ${featurexml.baseName}_____cutoff_plot.html
    """
}

process generate_queries_from_featurexmls {
    publishDir "${params.qal_outdir}/extracted_xics", mode:'copy', pattern: '*-queries.csv'
    container "luxii/unbequant:latest"

    input:
    tuple val(file_identifier), val(fdr), file(feature_with_idents), file(raw)

    output:
    tuple val(file_identifier), val(fdr), file("${raw.baseName}-queries.csv"), file(feature_with_idents), file(raw)

    """
    features_to_xic_extractor_table.py -featurexml ${feature_with_idents} -out_csv ${raw.baseName}-queries.csv
    """
}

process extracted_xics_from_hdf5_to_tsv {
    publishDir "${params.qal_outdir}/features_with_annotated_identifications", mode:'copy'
    container "luxii/unbequant:latest"


    input:
    tuple val(file_identifier), val(fdr), file(hdf5), file(original_query_file)

    output:
    tuple val(file_identifier), val(fdr), file("${hdf5.baseName}.tsv")

    """
    extract_hdf5_to_tsv.py -method ${params.qal_intensity_method} -hdf5_xic_file ${hdf5} -xic_query_file ${original_query_file} -out_tsv ${hdf5.baseName}.tsv
    """
}


process map_alignment_and_consensus_generation {
    publishDir "${params.qal_outdir}/features_with_annotated_identifications", mode:'copy'
    container "luxii/unbequant:latest"
    cpus 8

    input:
    tuple val(fdr), file(features)

    output:
    tuple val(fdr), file("consensus_____${fdr}.consensusXML")
    tuple val(fdr), file("*.trafoXML")

    """
    # Ensure same input order after rerunning by sorting the arguments
    SORTED_FEATURES=\$(echo ${features} | xargs -n1 | sort | xargs)
    NEW_FEATURES=()
    NEW_FEATURES_TRAFO=()
    for file in \${SORTED_FEATURES[@]}
    do
        NEW_FEATURES+=("\$(basename -- "\$file")_____aligned.featureXML")
        NEW_FEATURES_TRAFO+=("\$(basename -- "\$file")_____aligned.trafoXML")
    done

    # MapAlignerTreeGuided
    \$(get_cur_bin_dir.sh)/openms/usr/bin/MapAlignerTreeGuided -in \${SORTED_FEATURES[@]} -out \${NEW_FEATURES[@]} -trafo_out \${NEW_FEATURES_TRAFO[@]}

    # Consensus Generation
    \$(get_cur_bin_dir.sh)/openms/usr/bin/FeatureLinkerUnlabeled -algorithm:distance_MZ:max_difference ${params.qal_ppm_tolerance} -algorithm:distance_MZ:unit "ppm" ${params.qal_additional_feature_linker_params} -in \${NEW_FEATURES[@]} -out consensus_____${fdr}.consensusXML
    """
}

process visualize_RT_transoformations {
    publishDir "${params.qal_outdir}/visualizations___${fdr}", mode:'copy'
    container "luxii/unbequant:latest"

    input:
    tuple val(fdr), file(trafo_xmls)

    output:
    file("*.png")
    file("*.html")
    file("*.json")

    """
    CONCAT_TRAFOS=""
    for file in $trafo_xmls
    do
        CONCAT_TRAFOS+="\$file,"
    done
    CONCAT_TRAFOS=\$(echo \$CONCAT_TRAFOS | rev | cut -c2- | rev)

    # Limit time, since kaleido does not stop (Bug: https://github.com/plotly/Kaleido/issues/134 )
    PYTHONUNBUFFERED=1 timeout 3m visualize_RT_alignment.py -trafo_xmls \$CONCAT_TRAFOS  | true

    """
}



process generate_feature_ident_intesity_table {
    publishDir "${params.qal_outdir}/final_report___${fdr}", mode:'copy'
    container "luxii/unbequant:latest"

    input:
    tuple val(fdr), file(consensus), file(tsvs)

    output:
    tuple val(fdr), file("raw_quantification_with_identifications.tsv"), file("quantification_with_identifications_reduced.tsv"), file("quantification_with_identifications_only_intensities_and_ids.tsv")

    """
    CONCAT_TSVS=""
    for file in $tsvs
    do
        CONCAT_TSVS+="\$file,"
    done
    CONCAT_TSVS=\$(echo \$CONCAT_TSVS | rev | cut -c2- | rev)

    PYTHONUNBUFFERED=1 consensus_and_features_to_tsv.py -featurexmls_tsvs  \$CONCAT_TSVS -consensus $consensus -out_tsv raw_quantification_with_identifications.tsv -out_tsv_reduced quantification_with_identifications_reduced.tsv -out_tsv_minimal quantification_with_identifications_only_intensities_and_ids.tsv
    """
}



process generate_xlsx_reports_from_tables {
    publishDir "${params.qal_outdir}/final_report___${fdr}", mode:'copy'
    container "luxii/unbequant:latest"

    input:
    tuple val(fdr), file(full_file), file(reduced_file), file(minimal_file)

    output:
    file("*.xlsx")

    """
    PYTHONUNBUFFERED=1 write_xlsx_report.py -i ${full_file} -o ${full_file.baseName}.report.xlsx
    PYTHONUNBUFFERED=1 write_xlsx_report.py -i ${reduced_file} -o ${reduced_file.baseName}.report.xlsx
    PYTHONUNBUFFERED=1 write_xlsx_report.py -i ${minimal_file} -o ${minimal_file.baseName}.report.xlsx
    """
}

process generate_plots_from_tables {
    publishDir "${params.qal_outdir}/visualizations___${fdr}", mode:'copy'
    container "luxii/unbequant:latest"

    input:
    tuple val(fdr), file(full_file), file(reduced_file), file(minimal_file)

    output:
    file("*.html")

    """
    PYTHONUNBUFFERED=1 visualize_consensus_features_infos.py -minlh_parameter ${params.qal_minlh} -tsv_file ${full_file}
    """
}

