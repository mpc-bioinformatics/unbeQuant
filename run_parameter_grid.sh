#!/bin/bash

# ============================================================================
# Parameter Grid Runner for map_mzml_features.nf
# ============================================================================
# Tests multiple combinations of:
#   - RT alignment methods: aligntree, loess, none (disabled)
#   - m/z-RT weight ratios: 1.0 (unweighted), 0.001 (weighted)
#
# Creates 6 total runs with isolated output directories
# ============================================================================

set -e

# Configuration
CONFIG="nextflow_standalone.config"
MZML_FILES="results/mgfs_mzmls/*.mzML"
TSV_FILES="results/quantifications/features_with_annotated_identifications/*.tsv"
BASE_OUTDIR="results/parameter_grid"

# Define parameter grids
declare -a RT_METHODS=("aligntree" "loess" "none")
declare -a WEIGHT_RATIOS=(1.0 0.001)

# Create base output directory and log file
mkdir -p "$BASE_OUTDIR"
RESULTS_FILE="$BASE_OUTDIR/parameter_grid_results.txt"

# Initialize results log
echo "======================================================================" > "$RESULTS_FILE"
echo "Parameter Grid Execution Log" >> "$RESULTS_FILE"
echo "======================================================================" >> "$RESULTS_FILE"
echo "Started: $(date)" >> "$RESULTS_FILE"
echo "Work directory: $(pwd)" >> "$RESULTS_FILE"
echo "======================================================================" >> "$RESULTS_FILE"
echo "" >> "$RESULTS_FILE"

# track total runs
TOTAL_RUNS=0
SUCCESSFUL_RUNS=0
FAILED_RUNS=0

# Loop through all parameter combinations
for method in "${RT_METHODS[@]}"; do
    for weight_ratio in "${WEIGHT_RATIOS[@]}"; do
        TOTAL_RUNS=$((TOTAL_RUNS + 1))
        
        # Create descriptive run name
        if [ "$method" = "none" ]; then
            method_label="no_alignment"
        else
            method_label="$method"
        fi
        
        # Convert weight ratio to label (remove leading 0.)
        if [ "$weight_ratio" = "1.0" ] || [ "$weight_ratio" = "1" ]; then
            weight_label="unweighted"
        else
            weight_label="weighted_${weight_ratio}"
        fi
        
        run_name="${method_label}_${weight_label}"
        run_dir="$BASE_OUTDIR/$run_name"
        
        echo ""
        echo "═══════════════════════════════════════════════════════════════"
        echo "Run ${TOTAL_RUNS}/6: $run_name"
        echo "Output directory: $run_dir"
        echo "═══════════════════════════════════════════════════════════════"
        
        # Construct nextflow command
        cmd="nextflow run map_mzml_features.nf"
        cmd="$cmd -c $CONFIG"
        cmd="$cmd --mzml_files '$MZML_FILES'"
        cmd="$cmd --tsv_files '$TSV_FILES'"
        cmd="$cmd --main_outdir '$run_dir'"
        
        # Add RT alignment method (or disable it)
        if [ "$method" = "none" ]; then
            cmd="$cmd --mmf_compare_rt_alignment false"
        else
            cmd="$cmd --mmf_rt_alignment_method $method"
        fi
        
        # Add m/z-RT weight ratio
        cmd="$cmd --mmf_mz_rt_weight_ratio $weight_ratio"
        
        # Log the command
        echo "Command: $cmd" | tee -a "$RESULTS_FILE"
        echo "" | tee -a "$RESULTS_FILE"
        
        # Execute the workflow
        if eval "$cmd"; then
            SUCCESSFUL_RUNS=$((SUCCESSFUL_RUNS + 1))
            status="✓ SUCCESS"
            echo "$run_name: $status" >> "$RESULTS_FILE"
            echo "$status" 
        else
            FAILED_RUNS=$((FAILED_RUNS + 1))
            status="✗ FAILED"
            echo "$run_name: $status" >> "$RESULTS_FILE"
            echo "$status"
        fi
        
        echo "" | tee -a "$RESULTS_FILE"
        echo "---" >> "$RESULTS_FILE"
        
    done
done

# Print summary
echo ""
echo "======================================================================" | tee -a "$RESULTS_FILE"
echo "SUMMARY" | tee -a "$RESULTS_FILE"
echo "======================================================================" | tee -a "$RESULTS_FILE"
echo "Total runs: $TOTAL_RUNS" | tee -a "$RESULTS_FILE"
echo "Successful: $SUCCESSFUL_RUNS" | tee -a "$RESULTS_FILE"
echo "Failed: $FAILED_RUNS" | tee -a "$RESULTS_FILE"
echo "Completed: $(date)" | tee -a "$RESULTS_FILE"
echo "======================================================================" | tee -a "$RESULTS_FILE"
echo ""
echo "Output structure:"
echo "  $BASE_OUTDIR/"
echo "  ├── aligntree_unweighted/       (aligntree, weight=1.0)"
echo "  ├── aligntree_weighted_0.001/   (aligntree, weight=0.001)"
echo "  ├── loess_unweighted/           (loess, weight=1.0)"
echo "  ├── loess_weighted_0.001/       (loess, weight=0.001)"
echo "  ├── no_alignment_unweighted/    (disabled, weight=1.0)"
echo "  ├── no_alignment_weighted_0.001/(disabled, weight=0.001)"
echo "  └── parameter_grid_results.txt"
echo ""
echo "Full results saved to: $RESULTS_FILE"
