#!/bin/bash
# Complete workflow script, from converting the Thermo-RAW-files up to the qvalue-calculation on the identifications.
# This script purposely divides each step to a single nextflow execution, to provide reusable snippets, like identification, 
# fasta-with-feature-generation or raw file conversion. 

# THIS IS A TEMPLATE SCRIPT! PLEASE SET THE PARAMETERS ACCORDINGLY!

set -e

### PARAMETERS
### REQUIRED 
FOLDER="Example_Folder"
FOLDER_COMET_PROTEINS_CONFIG="${FOLDER}/<Path to comets config, for proteins-databases, with digestion turned on>"
FOLDER_COMET_PEPTIDES_CONFIG="${FOLDER}/<Path to comets config, for peptide-databases, with digestion turned off>"
PROTEINS_TEXT="${FOLDER}/<Path to the SP-EMBL-file of the Database>"
PROTEINS_FASTA="${FOLDER}/<Path to the FASTA-file of the same DATABASE as above>"

# Global-FASTA-Parameters
GLOBAL_FASTA_FEATURES_IN_GRAPH=" -ft NONE"
GLOBAL_FASTA_PEPTIDE_LIMITS="--pep_miscleavages 2 --pep_min_pep_length 5 --pep_max_weight 4500"

# Local-FASTA-Parameters
MZML_SPECIFIC_FASTA_PROTEINS_TEXT="${FOLDER}/<could be the same file as in PROTEINS_TEXT>"
MZML_SPECIFIC_FASTA_PROTEINS_FEATURES=" -ft SIGNAL -ft INIT_MET -ft VARIANT -ft CONFLICT -ft VAR_SEQ -ft PEPTIDE -ft PROPEP -ft CHAIN "
MZML_SPECIFIC_FASTA_SECOND_RUN="<true or false, depending if a second run should be done to include further proteins>"
MZML_SPECIFIC_FASTA_PROTEINS2_TEXT="${FOLDER}/<could be another SP-EMBL-file, but named different to MZML_SPECIFIC_FASTA_PROTEINS_TEXT>"
MZML_SPECIFIC_FASTA_PROTEINS2_FEATURES=" -ft SIGNAL -ft INIT_MET -ft CONFLICT -ft VAR_SEQ -ft PEPTIDE -ft PROPEP -ft CHAIN"
MZML_SPECIFIC_FASTA_GENERATION_FIXED_MODIFICATIONS=" -fm 'C:57.021464' "
MZML_SPECIFIC_FASTA_GENERATION_VARIABLE_MODIFICATIONS=" -vm 'M:15.994915' "
MZML_SPECIFIC_FASTA_GENERATION_MAX_PRECURSOR="4500"
MZML_SPECIFIC_FASTA_GENERATION_MAX_PROCESSES="64"
MZML_SPECIFIC_FASTA_GENERATION_MAX_VARIANTS="5"

# Identification-Parameters
IDENTIFICATION_FDR="0.01"
IDENTIFICATION_NUM_PARALLEL_SEARCHES="8"


# Generate profile MS1-mzml-data (MS2 peak picked to save storage)
nextflow run -with-report "${FOLDER}_convert_raw_to_mzml_report.html" -with-timeline "${FOLDER}_convert_raw_to_mzml_timeline.html" \
    convert_raw_to_mzml.nf \
    --thermo_raws "${FOLDER}" \
    # --additional_params " --noPeakPicking=1 " \
    --outdir "${FOLDER}/mzmls"







# ### WORKFLOWS (seperated)
# # Create MZMLs
# nextflow run -with-report "${FOLDER}_convert_raw_to_mgf_report.html" -with-timeline "${FOLDER}_convert_raw_to_mgf_timeline.html" \
#     convert_raw_to_mgf.nf \
#     --thermo_raws "${FOLDER}" \
#     --outdir "${FOLDER}/mgfs"

# # Create FASTA (Global)
# nextflow run -with-report "${FOLDER}_create_global_fasta_report.html" -with-timeline "${FOLDER}_create_global_fasta_timeline.html" \
#     create_global_fasta.nf \
#     --sp_embl_file "${PROTEINS_TEXT}" \
#     --outdir "${FOLDER}/global_fasta" \
#     --features_in_graphs "$GLOBAL_FASTA_FEATURES_IN_GRAPH" \
#     --peptide_limits "$GLOBAL_FASTA_PEPTIDE_LIMITS"

# # Create FASTA (MZML-Specific)
# nextflow run -with-report "${FOLDER}_create_fasta_report.html" -with-timeline "${FOLDER}_create_fasta_timeline.html" \
#     create_fasta.nf \
#     --mgf_folder "${FOLDER}/mgfs" \
#     --sp_embl_file "$MZML_SPECIFIC_FASTA_PROTEINS_TEXT"  \
#     --features_in_graphs "$MZML_SPECIFIC_FASTA_PROTEINS_FEATURES" \
#     --protgraph_run_2 "$MZML_SPECIFIC_FASTA_SECOND_RUN" \
#     --sp_embl_file_2 "$MZML_SPECIFIC_FASTA_PROTEINS2_TEXT" \
#     --features_in_graphs_2 "$MZML_SPECIFIC_FASTA_PROTEINS2_FEATURES" \
#     --fixed_modifications "$MZML_SPECIFIC_FASTA_GENERATION_FIXED_MODIFICATIONS" \
#     --variable_modifications "$MZML_SPECIFIC_FASTA_GENERATION_VARIABLE_MODIFICATIONS" \
#     --max_precursor "$MZML_SPECIFIC_FASTA_GENERATION_MAX_PRECURSOR" \
#     --num_procs_traversal "$MZML_SPECIFIC_FASTA_GENERATION_MAX_PROCESSES" \
#     --limit_varcount "$MZML_SPECIFIC_FASTA_GENERATION_MAX_VARIANTS" \
#     --outdir "${FOLDER}/ms2_specific_fasta"

# # Identification via global fasta
# nextflow run -with-report "${FOLDER}_ident_global_fasta_report.html" -with-timeline "${FOLDER}_ident_global_fasta_timeline.html" \
#     identification_via_comet.nf \
#     --mgf_folder "${FOLDER}/mgfs" \
#     --fasta_file "${FOLDER}/global_fasta/peptides.fasta" \
#     --fasta_not_from_protgraph "false" \
#     --outdir "${FOLDER}/global_fasta_idents" \
#     --tda "1" \
#     --fdr "$IDENTIFICATION_FDR" \
#     --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
#     --search_parameter_file "$FOLDER_COMET_PEPTIDES_CONFIG" \
#     --use_percolator 0

# # Identification via ms2_specific fasta
# nextflow run -with-report "${FOLDER}_ident_ms2_specific_fasta_report.html" -with-timeline "${FOLDER}_ident_ms2_specific_fasta_timeline.html" \
#     identification_via_comet.nf \
#     --mgf_folder "${FOLDER}/mgfs" \
#     --fasta_file "${FOLDER}/ms2_specific_fasta/peptides.fasta" \
#     --fasta_not_from_protgraph "false" \
#     --outdir "${FOLDER}/ms2_specific_fasta_idents" \
#     --tda "1" \
#     --fdr "$IDENTIFICATION_FDR" \
#     --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
#     --search_parameter_file "$FOLDER_COMET_PEPTIDES_CONFIG" \
#     --use_percolator 0

# # Identification via protein_fasta
# nextflow run -with-report "${FOLDER}_ident_protein_fasta_report.html" -with-timeline "${FOLDER}_ident_protein_fasta_timeline.html" \
#     identification_via_comet.nf \
#     --mgf_folder "${FOLDER}/mgfs" \
#     --fasta_file "$PROTEINS_FASTA" \
#     --fasta_not_from_protgraph "true" \
#     --outdir "${FOLDER}/protein_fasta_idents" \
#     --tda 1 \
#     --fdr "$IDENTIFICATION_FDR" \
#     --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
#     --search_parameter_file "$FOLDER_COMET_PROTEINS_CONFIG" \
#     --use_percolator 0

# # Identification via global fasta
# nextflow run -with-report "${FOLDER}_ident_global_fasta_with_percolator_report.html" -with-timeline "${FOLDER}_ident_global_fasta_with_percolator_timeline.html" \
#     identification_via_comet.nf \
#     --mgf_folder "${FOLDER}/mgfs" \
#     --fasta_file "${FOLDER}/global_fasta/peptides.fasta" \
#     --fasta_not_from_protgraph "false" \
#     --outdir "${FOLDER}/global_fasta_with_percolator_idents" \
#     --tda "1" \
#     --fdr "$IDENTIFICATION_FDR" \
#     --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
#     --search_parameter_file "$FOLDER_COMET_PEPTIDES_CONFIG" \
#     --use_percolator 1

# # Identification via ms2_specific fasta
# nextflow run -with-report "${FOLDER}_ident_ms2_specific_fasta_with_percolator_report.html" -with-timeline "${FOLDER}_ident_ms2_specific_fasta_with_percolator_timeline.html" \
#     identification_via_comet.nf \
#     --mgf_folder "${FOLDER}/mgfs" \
#     --fasta_file "${FOLDER}/ms2_specific_fasta/peptides.fasta" \
#     --fasta_not_from_protgraph "false" \
#     --outdir "${FOLDER}/ms2_specific_fasta_with_percolator_idents" \
#     --tda "1" \
#     --fdr "$IDENTIFICATION_FDR" \
#     --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
#     --search_parameter_file "$FOLDER_COMET_PEPTIDES_CONFIG" \
#     --use_percolator 1

# # Identification via protein_fasta
# nextflow run -with-report "${FOLDER}_ident_protein_fasta_with_percolator_report.html" -with-timeline "${FOLDER}_ident_protein_fasta_with_percolator_timeline.html" \
#     identification_via_comet.nf \
#     --mgf_folder "${FOLDER}/mgfs" \
#     --fasta_file "$PROTEINS_FASTA" \
#     --fasta_not_from_protgraph "true" \
#     --outdir "${FOLDER}/protein_fasta_with_percolator_idents" \
#     --tda 1 \
#     --fdr "$IDENTIFICATION_FDR" \
#     --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
#     --search_parameter_file "$FOLDER_COMET_PROTEINS_CONFIG" \
#     --use_percolator 1

# ### Export to excel to get a quick overview of found results
# nextflow run export_psms_to_table.nf \
#     --identification_folder "${FOLDER}/global_fasta_idents" \
#     --outdir "${FOLDER}/global_fasta_idents" \
#     --fdr "$IDENTIFICATION_FDR"

# nextflow run export_psms_to_table.nf \
#     --identification_folder "${FOLDER}/ms2_specific_fasta_idents" \
#     --outdir "${FOLDER}/ms2_specific_fasta_idents" \
#     --fdr "$IDENTIFICATION_FDR"

# nextflow run export_psms_to_table.nf \
#     --identification_folder "${FOLDER}/protein_fasta_idents" \
#     --outdir "${FOLDER}/protein_fasta_idents" \
#     --fdr "$IDENTIFICATION_FDR"

# nextflow run export_psms_to_table.nf \
#     --identification_folder "${FOLDER}/global_fasta_with_percolator_idents" \
#     --outdir "${FOLDER}/global_fasta_with_percolator_idents" \
#     --fdr "$IDENTIFICATION_FDR"

# nextflow run export_psms_to_table.nf \
#     --identification_folder "${FOLDER}/ms2_specific_fasta_with_percolator_idents" \
#     --outdir "${FOLDER}/ms2_specific_fasta_with_percolator_idents" \
#     --fdr "$IDENTIFICATION_FDR"

# nextflow run export_psms_to_table.nf \
#     --identification_folder "${FOLDER}/protein_fasta_with_percolator_idents" \
#     --outdir "${FOLDER}/protein_fasta_with_percolator_idents" \
#     --fdr "$IDENTIFICATION_FDR"
