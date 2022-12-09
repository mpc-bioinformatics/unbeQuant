#!/bin/bash
# Complete workflow script, from converting the Thermo-RAW-files up to the qvalue-calculation on the identifications.
# This script purposely divides each step to a single nextflow execution, to provide reusable snippets, like identification, 
# fasta-with-feature-generation or raw file conversion. 

# THIS IS A TEMPLATE SCRIPT! PLEASE SET THE PARAMETERS ACCORDINGLY!

set -e

### PARAMETERS
### REQUIRED 
FOLDER="raws1"
FOLDER_COMET_PEPTIDES_CONFIG="${FOLDER}/<Path to comets config, for peptide-databases, with digestion turned off>"

# Local-FASTA-Parameters
MZML_SPECIFIC_FASTA_PROTEINS_TEXT="fasta_generation/f1/170508_SIHUMI_n8_species_cRAP_concatenated_target_reheadered.txt"
MZML_SPECIFIC_FASTA_PROTEINS_FEATURES=" -ft NONE "
MZML_SPECIFIC_FASTA_SECOND_RUN="false"
MZML_SPECIFIC_FASTA_PROTEINS2_TEXT="fasta_generation/f1/170508_SIHUMI_n8_species_cRAP_concatenated_target_reheadered.fasta"
MZML_SPECIFIC_FASTA_PROTEINS2_FEATURES=" -ft NONE"
MZML_SPECIFIC_FASTA_GENERATION_FIXED_MODIFICATIONS=" -fm 'C:57.021464' "
MZML_SPECIFIC_FASTA_GENERATION_VARIABLE_MODIFICATIONS=" -vm 'M:15.994915' "
MZML_SPECIFIC_FASTA_GENERATION_MAX_PRECURSOR="8000"
MZML_SPECIFIC_FASTA_GENERATION_MAX_PROCESSES="32"
MZML_SPECIFIC_FASTA_GENERATION_MAX_VARIANTS="-1"

# Identification-Parameters
IDENTIFICATION_FDR="0.01"
IDENTIFICATION_NUM_PARALLEL_SEARCHES="4"




# Create FASTA (MZML-Specific)
nextflow run -with-report "${FOLDER}_create_fasta_report.html" -with-timeline "${FOLDER}_create_fasta_timeline.html" \
    create_fasta.nf \
    --mgf_folder "${FOLDER}/mgfs" \
    --sp_embl_file "$MZML_SPECIFIC_FASTA_PROTEINS_TEXT"  \
    --features_in_graphs "$MZML_SPECIFIC_FASTA_PROTEINS_FEATURES" \
    --protgraph_run_2 "$MZML_SPECIFIC_FASTA_SECOND_RUN" \
    --sp_embl_file_2 "$MZML_SPECIFIC_FASTA_PROTEINS2_TEXT" \
    --features_in_graphs_2 "$MZML_SPECIFIC_FASTA_PROTEINS2_FEATURES" \
    --fixed_modifications "$MZML_SPECIFIC_FASTA_GENERATION_FIXED_MODIFICATIONS" \
    --variable_modifications "$MZML_SPECIFIC_FASTA_GENERATION_VARIABLE_MODIFICATIONS" \
    --max_precursor "$MZML_SPECIFIC_FASTA_GENERATION_MAX_PRECURSOR" \
    --num_procs_traversal "$MZML_SPECIFIC_FASTA_GENERATION_MAX_PROCESSES" \
    --limit_varcount "$MZML_SPECIFIC_FASTA_GENERATION_MAX_VARIANTS" \
    --outdir "${FOLDER}/ms2_specific_fasta" \
    --query_ppm "10"

# Identification via ms2_specific fasta
nextflow run -with-report "${FOLDER}_ident_ms2_specific_fasta_report.html" -with-timeline "${FOLDER}_ident_ms2_specific_fasta_timeline.html" \
    identification_via_comet.nf \
    --mgf_folder "${FOLDER}/mgfs" \
    --fasta_file "${FOLDER}/ms2_specific_fasta/peptides.fasta" \
    --fasta_not_from_protgraph "false" \
    --outdir "${FOLDER}/ms2_specific_fasta_idents" \
    --tda "1" \
    --fdr "$IDENTIFICATION_FDR" \
    --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
    --search_parameter_file "$FOLDER_COMET_PEPTIDES_CONFIG" \
    --use_percolator 0

# Identification via ms2_specific fasta
nextflow run -with-report "${FOLDER}_ident_ms2_specific_fasta_with_percolator_report.html" -with-timeline "${FOLDER}_ident_ms2_specific_fasta_with_percolator_timeline.html" \
    identification_via_comet.nf \
    --mgf_folder "${FOLDER}/mgfs" \
    --fasta_file "${FOLDER}/ms2_specific_fasta/peptides.fasta" \
    --fasta_not_from_protgraph "false" \
    --outdir "${FOLDER}/ms2_specific_fasta_with_percolator_idents" \
    --tda "1" \
    --fdr "$IDENTIFICATION_FDR" \
    --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
    --search_parameter_file "$FOLDER_COMET_PEPTIDES_CONFIG" \
    --use_percolator 1

nextflow run export_psms_to_table.nf \
    --identification_folder "${FOLDER}/ms2_specific_fasta_idents" \
    --outdir "${FOLDER}/ms2_specific_fasta_idents" \
    --fdr "$IDENTIFICATION_FDR"


nextflow run export_psms_to_table.nf \
    --identification_folder "${FOLDER}/ms2_specific_fasta_with_percolator_idents" \
    --outdir "${FOLDER}/ms2_specific_fasta_with_percolator_idents" \
    --fdr "$IDENTIFICATION_FDR"



FOLDER="raws2"
# Local-FASTA-Parameters
MZML_SPECIFIC_FASTA_PROTEINS_TEXT="fasta_generation/f2/gut_1_igc__concatenated_target_reheadered.txt"
MZML_SPECIFIC_FASTA_PROTEINS_FEATURES=" -ft NONE "
MZML_SPECIFIC_FASTA_SECOND_RUN="false"
MZML_SPECIFIC_FASTA_PROTEINS2_TEXT="fasta_generation/f2/gut_1_igc__concatenated_target_reheadered.fasta"
MZML_SPECIFIC_FASTA_PROTEINS2_FEATURES=" -ft NONE"
MZML_SPECIFIC_FASTA_GENERATION_FIXED_MODIFICATIONS=" -fm 'C:57.021464' "
MZML_SPECIFIC_FASTA_GENERATION_VARIABLE_MODIFICATIONS=" -vm 'M:15.994915' "
MZML_SPECIFIC_FASTA_GENERATION_MAX_PRECURSOR="8000"
MZML_SPECIFIC_FASTA_GENERATION_MAX_PROCESSES="32"
MZML_SPECIFIC_FASTA_GENERATION_MAX_VARIANTS="-1"


# Identification via ms2_specific fasta
nextflow run -with-report "${FOLDER}_ident_ms2_specific_fasta_report.html" -with-timeline "${FOLDER}_ident_ms2_specific_fasta_timeline.html" \
    identification_via_comet.nf \
    --mgf_folder "${FOLDER}/mgfs" \
    --fasta_file "${FOLDER}/ms2_specific_fasta/peptides.fasta" \
    --fasta_not_from_protgraph "false" \
    --outdir "${FOLDER}/ms2_specific_fasta_idents" \
    --tda "1" \
    --fdr "$IDENTIFICATION_FDR" \
    --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
    --search_parameter_file "$FOLDER_COMET_PEPTIDES_CONFIG" \
    --use_percolator 0 \
    --query_ppm "10"

# Identification via ms2_specific fasta
nextflow run -with-report "${FOLDER}_ident_ms2_specific_fasta_with_percolator_report.html" -with-timeline "${FOLDER}_ident_ms2_specific_fasta_with_percolator_timeline.html" \
    identification_via_comet.nf \
    --mgf_folder "${FOLDER}/mgfs" \
    --fasta_file "${FOLDER}/ms2_specific_fasta/peptides.fasta" \
    --fasta_not_from_protgraph "false" \
    --outdir "${FOLDER}/ms2_specific_fasta_with_percolator_idents" \
    --tda "1" \
    --fdr "$IDENTIFICATION_FDR" \
    --num_parallel_searches "$IDENTIFICATION_NUM_PARALLEL_SEARCHES" \
    --search_parameter_file "$FOLDER_COMET_PEPTIDES_CONFIG" \
    --use_percolator 1

nextflow run export_psms_to_table.nf \
    --identification_folder "${FOLDER}/ms2_specific_fasta_idents" \
    --outdir "${FOLDER}/ms2_specific_fasta_idents" \
    --fdr "$IDENTIFICATION_FDR"


nextflow run export_psms_to_table.nf \
    --identification_folder "${FOLDER}/ms2_specific_fasta_with_percolator_idents" \
    --outdir "${FOLDER}/ms2_specific_fasta_with_percolator_idents" \
    --fdr "$IDENTIFICATION_FDR"

