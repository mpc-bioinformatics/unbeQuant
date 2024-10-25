#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.qal_spectra_files = "$PWD/raws"  // .RAW/.d-files
params.qal_mzmls = "$PWD/mzmls"  // mzML-files
params.qal_idents = "$PWD/raws/tsvs"  // Folder containing Identifications in TSV-format
params.qal_idents_blob_filter = "*qvalue_no_decoys_fdr_0.0[15].tsv"  // Should be TSV-files of Identification (already FDR-filtered), containing the columns: "charge", "plain_peptide", "used_score", "retention_time", "exp_mass_to_charge", "fasta_id", "fasta_desc"
params.qal_outdir = "$PWD/results"  // Output-Directory of the quantification results (splitted by the fdr)


// Parameters for Feature Detection
// params.resulolution_featurefinder = "-algortihm:mass_trace:mz_tolerance 0.02 -algorithm:isotopic_pattern:mz_tolerance 0.04"  // Parameters for Low Resolution Machines (E.G.: Q-TOF)
params.qal_resolution_featurefinder = "-algorithm:mass_trace:mz_tolerance 0.004 -algorithm:isotopic_pattern:mz_tolerance 0.005"   // Parameters for High Resolution Machines (E.G.: LTQ-OrbiTrap)
// REMOVE THE OTHTER TEST PARAMETERS!
// params.qal_resolution_featurefinder = "-algorithm:mass_trace:mz_tolerance 0.005 -algorithm:isotopic_pattern:mz_tolerance 0.005 -algorithm:intensity:bins 15"   // TESTING!
// params.qal_resolution_featurefinder = "-algorithm:mass_trace:mz_tolerance 0.005 -algorithm:isotopic_pattern:mz_tolerance 0.005 -algorithm:intensity:bins 15 -algorithm:feature:max_rt_span 5"   // TESTING!
// params.qal_resolution_featurefinder = "-algorithm:mass_trace:mz_tolerance 0.008 -algorithm:isotopic_pattern:mz_tolerance 0.010 -algorithm:intensity:bins 15 -algorithm:feature:max_rt_span 7.5 -algorithm:feature:min_isotope_fit 0.7  -algorithm:seed:min_score 0.7 -algorithm:mass_trace:max_missing 3 -algorithm:mass_trace:slope_bound 0.1 "   // TESTING!
params.qal_considered_charges_low = "2"  // Charges for the feature finder to use to extract features.
params.qal_considered_charges_high = "7"  // Charges for the feature finder to use to extract features.
params.qal_protgraph_was_used = false  // A Flag which is needed for the output to know which parsing mode and which column of "fasta_id" and "fasta_desc" needs to be taken
params.qal_limit_num_of_parallel_feature_finders = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)

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
    identifications_tuple = identifications
        .map { file -> tuple(
            file.baseName.substring(file.baseName.indexOf("fdr_") + 4, file.baseName.indexOf("fdr_") + 8), 
            file
        ) }.groupTuple()

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

        //// Generate file_identifier and match features with identifications (on multiple fdrs)
        spectra_files_tuple = spectra_files.map { file -> tuple(file.baseName, file) }
        featurexmls_tuple = create_feature_xml.out.map { file -> tuple(file.baseName, file) }
        // Get all the single fdrs
        in_identifications_tuple = identifications_tuple.transpose().map { it -> tuple(it[1].baseName.split("_____")[0], it[0], it[1]) }
        // Match with identifications using file_identifier (it[0]) and fdr (it[1])
        in_featurexmls_tuple = in_identifications_tuple.map { it -> it[1] } .unique().combine(featurexmls_tuple).map { it -> tuple(it[1], it[0], it[2]) }
        matched_features_with_idents_tuple = in_identifications_tuple.join(in_featurexmls_tuple, by: [0,1])
        // Do the actual matching
        match_feature_with_idents(matched_features_with_idents_tuple)

        //// Retrieve the quant absolute values via the TRFP XIC
        // Get all the single fdrs
        in_spectra_files_tuple = in_identifications_tuple.map { it -> it[1] } .unique().combine(spectra_files_tuple).map { it -> tuple(it[1], it[0], it[2]) }
        // Match with identifications using file_identifier (it[0]) and fdr (it[1])
        matched_ident_features_with_raws = match_feature_with_idents.out.join(in_spectra_files_tuple, by: [0,1])
        
        //// Retrieve XICs in hdf5 format
        // First, generate queries
        generate_queries_from_featurexmls(matched_ident_features_with_raws)
        // Then, extract via xic_extractor
        extract_xics_channel = generate_queries_from_featurexmls.out.map {
            it -> tuple(it[4], it[2])
        }
        retrieve_xics_from_raw_spectra(extract_xics_channel)
        // Finally generate the resulting tsv file containing all the data
        xics_and_remaining_data = retrieve_xics_from_raw_spectra.out.join(generate_queries_from_featurexmls.out, by: [0])

        extracted_xics_from_hdf5_to_tsv(
            xics_and_remaining_data.map {it -> tuple(it[0], it[2], it[1], it[3]) }
        )

        //// Apply the MapAligner and Consensus_generator
        identified_features_by_fdr = match_feature_with_idents.out.map { it -> tuple(it[1], it[2]) }.groupTuple()
        map_alignment_and_consensus_generation(identified_features_by_fdr)


        //// Generate the final statistics and visualizations
        fdr_and_feature_tsvs = extracted_xics_from_hdf5_to_tsv.out.map { it -> tuple(it[1], it[2]) }.groupTuple()
        consensus_with_feature_tsvs = map_alignment_and_consensus_generation.out[0].join(fdr_and_feature_tsvs, by: 0)
        visualize_RT_transoformations(map_alignment_and_consensus_generation.out[1])

        generate_feature_ident_intesity_table(consensus_with_feature_tsvs)

    emit:
        generate_feature_ident_intesity_table.out[0]
        consensus_with_feature_tsvs
        match_feature_with_idents
}


process create_feature_xml {
    maxForks params.qal_limit_num_of_parallel_feature_finders
    stageInMode "copy"
    container "luxii/unbequant:latest"

    input:
    file mzml

    output:
    file("${mzml.baseName}.featureXML")

    """
    # \$(get_cur_bin_dir.sh)/openms/usr/bin/NoiseFilterGaussian -in ${mzml} -out ${mzml.baseName}_filtered.mzML
    # \$(get_cur_bin_dir.sh)/openms/usr/bin/FeatureFinderIsotopeWavelet -in ${mzml.baseName}_filtered.mzML -out ${mzml.baseName}.featureXML -algorithm:hr_data

    # \$(get_cur_bin_dir.sh)/openms/usr/bin/FileFilter -in ${mzml} -out ${mzml.baseName}_filtered.mzML
    # \$(get_cur_bin_dir.sh)/openms/usr/bin/FeatureFinderMultiplex -in ${mzml.baseName}_filtered.mzML -out ${mzml.baseName}.featureXML

    \$(get_cur_bin_dir.sh)/openms/usr/bin/FeatureFinderCentroided -in ${mzml} -out ${mzml.baseName}.featureXML -algorithm:isotopic_pattern:charge_low ${params.qal_considered_charges_low} -algorithm:isotopic_pattern:charge_high ${params.qal_considered_charges_high} ${params.qal_resolution_featurefinder}
    # We do not use multiplex, it seems to be broken. Mem usage is way over 40 GB per RAW file following by a "std::bad_alloc"
    """
}


process match_feature_with_idents {
    container "luxii/unbequant:latest"

    input:
    tuple val(file_identifier), val(fdr), file(ident_tsv), file(featurexml)

    output:
    tuple val(file_identifier), val(fdr), file("${featurexml.baseName}_____${fdr}_____with_identifications.featureXML")

    """
    convert_ident_to_idXML.py -use_protgraph ${params.qal_protgraph_was_used} -tsv_file ${ident_tsv} -o ${ident_tsv.baseName}.idXML
    \$(get_cur_bin_dir.sh)/openms/usr/bin/IDMapper -id ${ident_tsv.baseName}.idXML -mz_reference precursor -in ${featurexml} -out ${featurexml.baseName}_____${fdr}_____with_identifications.featureXML
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
    extract_hdf5_to_tsv.py -hdf5_xic_file ${hdf5} -xic_query_file ${original_query_file} -out_tsv ${hdf5.baseName}.tsv
    """
}


process map_alignment_and_consensus_generation {
    publishDir "${params.qal_outdir}/features_with_annotated_identifications", mode:'copy'
    container "luxii/unbequant:latest"

    input:
    tuple val(fdr), file(features)

    output:
    tuple val(fdr), file("consensus_____${fdr}.consensusXML")
    tuple val(fdr), file("*.trafoXML")

    """
    NEW_FEATURES=()
    NEW_FEATURES_TRAFO=()
    for file in ${features}
    do
        NEW_FEATURES+=("\$(basename -- "\$file")_____aligned.featureXML")
        NEW_FEATURES_TRAFO+=("\$(basename -- "\$file")_____aligned.trafoXML")
    done

    # MapAlignerTreeGuided
    \$(get_cur_bin_dir.sh)/openms/usr/bin/MapAlignerTreeGuided -in ${features} -out \${NEW_FEATURES[@]} -trafo_out \${NEW_FEATURES_TRAFO[@]}

    # Consensus Generation
    \$(get_cur_bin_dir.sh)/openms/usr/bin/FeatureLinkerUnlabeled -in \${NEW_FEATURES[@]} -out consensus_____${fdr}.consensusXML
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
    publishDir "${params.qal_outdir}/statistics___${fdr}", mode:'copy'
    container "luxii/unbequant:latest"

    input:
    tuple val(fdr), file(consensus), file(tsvs)

    output:
    file("raw_quantification_with_identifications.tsv")
    file("quantification_with_identifications_reduced.tsv")
    file("quantification_with_identifications_only_intensities_and_ids.tsv")

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