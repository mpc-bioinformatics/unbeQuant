#!/bin/bash

# ============================================================================
# Parameter Grid Runner for Complete UnbeQuant Pipeline - FILTER MODE TEST
# ============================================================================
# Derived from run_parameter_grid_feature_modes.sh, but the grid dimension is the
# component-reduction strategy instead of the feature centre mode.
#
# Runs the ANCHOR line: map_mzml_features_v3_anchor.nf + bin_anchor/*.py.
# The original bin/ and map_mzml_features_v2_clusterless.nf are not touched, so
# this grid can run alongside the normal pipeline.
#
# Strategies under test (--mmf_filter_method):
#   duplicate_only      - keep one vertex per file, evict nothing else. The
#                         historical behaviour, here as the control.
#   anchor_connectivity - pick the vertex covering the most distinct files inside
#                         its tolerance, keep only what fits around it, then
#                         resolve remaining duplicates by most-retained-connections
#                         (ties: shortest total distance).
#   anchor_nearest      - same anchor, but keep each file's vertex closest to it.
#   window              - solve the pairwise constraint EXACTLY. "every pair inside the
#                         cutoffs" == "the set's bounding box fits the cutoffs", so this
#                         enumerates (m/z, RT) windows of the component, takes the one
#                         covering the most files, and picks the tightest one-per-file
#                         selection inside it. Uses the FULL cutoff where BFS needs HALF,
#                         so it buys back coverage at the same zero-violation guarantee:
#                         cycle 3 gave 83.1% file coverage vs BFS 2.5:50's 78.1%, both at
#                         0.00% violating pairs. Orphans are negligible (3 of 15947 kept).
#   bfs                 - grow outward from a start vertex along EDGES instead of by
#                         coordinate tolerance, claiming at most one vertex per file
#                         (path nodes included). Bounded by accumulated ppm/RT budgets.
#                         Every vertex is tried as start; best coverage then shortest
#                         total distance wins.
#
# BFS budgets (--mmf_bfs_ppm_budget / --mmf_bfs_rt_budget), and why half the cutoff:
#   The budget caps the ACCUMULATED offset along the path from the start. The triangle
#   inequality means the direct offset never exceeds it, so a budget of B puts every
#   member within B of the start and every PAIR within 2B. At HALF the pairing cutoff
#   that makes the finished component pairwise-conforming - measured on cycle 3 at
#   exactly ZERO violating pairs, versus ~14% unbudgeted and 21.5% for the anchor.
#   The cost is coverage: 78% of file slots per cycle against the anchor's 98%.
#   Whether the cycle loop recovers that is what a full run has to answer.
#
# Anchor tolerance semantics (--mmf_anchor_radius_mode), anchor methods only:
#   full  - radius = the cutoff. Members are within the cutoff OF THE ANCHOR, so
#           two of them can sit on opposite sides and be up to 2x apart. Measured:
#           worst span exactly 1.99x, zero components beyond 2x.
#   prune - full radius, then drop members that clash pairwise. Reached 100%
#           pairwise conformance in the prototype.
#   half  - radius = cutoff/2, so the triangle inequality guarantees every pair
#           conforms. Strictest, evicts the most.
#
# WHY THIS GRID EXISTS
#   Connected components make membership transitive while the cutoff is pairwise:
#   A-B and B-C inside the cutoff does not keep A-C inside it, and nothing
#   downstream re-checks. A duplicate_only run produced a component spanning
#   7181 s of RT - 47.9x the cutoff - emitted as one consensus feature. The anchor
#   methods bound that. This grid measures what the bound costs.
#
# ---------------------------------------------------------------------------
# IMPORTANT - why -resume is OFF by default:
#   The .nf processes invoke the helper scripts by ABSOLUTE path, e.g.
#     python3 /workspaces/unbeQuant/bin_anchor/extract_feature_data.py ...
#   Nextflow hashes the .nf script and the params, but NOT the contents of those
#   bin_anchor/*.py files. So if you edit one and then use -resume, Nextflow can
#   silently reuse results computed with the OLD code. Since the filter methods
#   live in exactly those files, resuming would defeat the purpose of this grid.
#
#   THE ONE SAFE EXCEPTION: extraction and RT alignment do not depend on the
#   filter method at all, so enabling -resume across runs of THIS grid (with
#   bin_anchor/ unchanged) correctly reuses them and only re-runs the cycle loop.
#   That roughly halves total runtime. Set NF_RESUME="-resume" for that, but only
#   when you have not edited bin_anchor/ since the first run.
# ---------------------------------------------------------------------------
# ============================================================================

set -e

# ---- Optional start delay --------------------------------------------------
# Wait this many hours before starting the first run (0 = start immediately).
# Fractional values are allowed, e.g. 2.5. Can also be passed as the first
# argument, which overrides the value set here:
#   ./run_parameter_grid_filter_modes.sh 6      # start in 6 hours
# All configuration is validated BEFORE the delay, so a misconfiguration fails
# immediately instead of after the wait.
START_DELAY_HOURS="${1:-0}"

if ! [[ "$START_DELAY_HOURS" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "ERROR: start delay must be a non-negative number of hours, got: '$START_DELAY_HOURS'"
    echo "Usage: $0 [hours]     e.g. $0 6     or     $0 2.5"
    exit 1
fi

# Configuration
CONFIG="nextflow_standalone.config"
MZML_FILES="results/mgfs_mzmls/*.mzML"
TSV_FILES="results/quantifications/features_with_annotated_identifications/*.tsv"
BASE_OUTDIR="results/parameter_grid/Gold_Standard/Clusterless/FilterModes"

# The anchor workflow. Publishes to <run_dir>/feature_analysis_clusterless_anchor/.
WORKFLOW="map_mzml_features_v3_anchor.nf"

# Shared work directory across grid runs. Keeping one lets -resume reuse the
# extraction and RT alignment stages, which are identical for every filter mode.
NF_WORKDIR="$PWD/work_filtergrid"

# Set to "-resume" to enable Nextflow resume. See the warning above before doing so.
NF_RESUME=""

# ---- Which pipeline steps to run -------------------------------------------
#   RUN_STEP1M - component metrics (avg linkage, pep_ident conflicts, cutoff
#                conformance). Cheap, and it is the actual point of this grid.
#   RUN_STEP2  - process_clusters_consensusxml.nf
#   RUN_STEP3  - ProteoBench conversion (consumes step 2 output, so needs RUN_STEP2)
RUN_STEP1M=true
RUN_STEP2=true
RUN_STEP3=false

# ---- Pairwise within-component distance analysis ---------------------------
# Runs inside step 2 (process_clusters_consensusxml.nf) on the components JSON.
# Charts every pair inside each component against the pairing cutoffs, which is
# what shows how far transitive membership drifts past a pairwise cutoff - the
# same drift STEP 1M summarizes as pairwise_pct / worst_rt_x, but resolved by
# component size instead of collapsed to one number.
# Needs RUN_STEP2=true, since it is a process of that workflow.
RUN_PAIRWISE=true
# light | dark | both
PAIRWISE_THEME="both"

# ---- Step 2 sub-options: normalization + differential analysis -------------
RUN_NORMALIZATION=true
RUN_TTEST=true
SAMPLE_CONDITIONS="files/sample_conditions.tsv"
# loess | median | quantile | lts | nonorm
NORMALIZATION_METHOD="loess"
# Comma-separated pairs, e.g. "C1:C2,C1:C5". Empty = all pairwise combinations.
TTEST_COMPARISONS=""
R_CONTAINER="unbequant:r"

# ---- Preflight -------------------------------------------------------------
if [ ! -f "$WORKFLOW" ]; then
    echo "ERROR: workflow not found: $WORKFLOW"
    exit 1
fi
if [ ! -d "bin_anchor" ]; then
    echo "ERROR: bin_anchor/ not found. This grid runs the anchor file line."
    exit 1
fi
if [ "$RUN_NORMALIZATION" = "true" ] && [ ! -f "$SAMPLE_CONDITIONS" ]; then
    echo "ERROR: RUN_NORMALIZATION=true but sample conditions file not found: $SAMPLE_CONDITIONS"
    exit 1
fi
if [ "$RUN_NORMALIZATION" = "true" ] && ! docker image inspect "$R_CONTAINER" >/dev/null 2>&1; then
    echo "ERROR: RUN_NORMALIZATION=true but the container '$R_CONTAINER' is not available locally."
    echo "       Build it with: docker build -f Dockerfile.unbequant_r -t $R_CONTAINER ."
    exit 1
fi
if [ "$RUN_STEP3" = "true" ] && [ "$RUN_STEP2" != "true" ]; then
    echo "WARNING: RUN_STEP3=true requires RUN_STEP2=true (step 3 consumes step 2 output)."
    echo "         Disabling step 3."
    RUN_STEP3=false
fi
if [ "$RUN_PAIRWISE" = "true" ] && [ "$RUN_STEP2" != "true" ]; then
    echo "WARNING: RUN_PAIRWISE=true requires RUN_STEP2=true (it is a process of"
    echo "         process_clusters_consensusxml.nf). Disabling the pairwise analysis."
    RUN_PAIRWISE=false
fi
if [ "$RUN_PAIRWISE" = "true" ] && [ ! -f "bin/analyze_pairwise_components.py" ]; then
    echo "ERROR: RUN_PAIRWISE=true but bin/analyze_pairwise_components.py not found."
    exit 1
fi

TOTAL_STEPS=1
[ "$RUN_STEP1M" = "true" ] && TOTAL_STEPS=$((TOTAL_STEPS + 1))
[ "$RUN_STEP2" = "true" ] && TOTAL_STEPS=$((TOTAL_STEPS + 1))
[ "$RUN_STEP3" = "true" ] && TOTAL_STEPS=$((TOTAL_STEPS + 1))

# ---- Parameter grids -------------------------------------------------------
# Primary dimension under test.
declare -a FILTER_METHODS=(
    #"duplicate_only"
    #"anchor_connectivity"
    #"anchor_nearest"
    "window"
    "bfs"
    
)

# Secondary dimension, anchor methods only. duplicate_only ignores it, so it is
# run ONCE regardless of how many radius modes are listed (see the guard below).
# Each extra radius mode adds one full run per anchor method - roughly 1h45m each.
declare -a ANCHOR_RADIUS_MODES=(
    #"full"
    #"prune"
    "half"
)

# Anchor RT tolerance, seconds. Wider than mmf_rt_cutoff on purpose: the exported
# y_center is PRE-alignment, so a coordinate-based RT test carries the alignment
# offset (measured at up to ~30s per edge).
declare -a ANCHOR_RT_CUTOFFS=(
    100
    #150
)

# ---- BFS dimensions (filter_method = bfs only) ------------------------------
# Depth is a weak knob - 1 -> 2 moved linkage by 0.002 on cycle 3 - so only depth 1
# is enabled by default. Uncomment 2 to confirm that on the full loop.
declare -a BFS_DEPTHS=(
    #1
    #2
    #5
    20
)
# "ppm:rt" accumulated budget pairs. 2.5:50 is half the ppm5/rt100 pairing cutoff and
# gave zero violating pairs; 5:100 is the full cutoff, which bounds pairs at 2x and
# keeps more coverage.
declare -a BFS_BUDGETS=(
    #"2.5:50"
    "5:100"
    #"5:150"
)

# ---- Window dimensions (filter_method = window only) ------------------------
# "ppm:rt" window sizes, applied SET-WISE (mz_max-mz_min, rt_max-rt_min). Unlike the
# BFS budgets these take the FULL pairing cutoff, because the window checks the
# pairwise constraint directly instead of bounding it via the triangle inequality -
# which is what buys back coverage at the same zero-violation guarantee.
declare -a WINDOW_SIZES=(
    "5:100"
    #"2.5:50"
    #"10:150"
)

# Held fixed so the dimensions above are the only things varying.
declare -a FEATURE_MODES=("first_mean")
declare -a COORD_SOURCES=("center" )
declare -a MAX_MATCH_CANDIDATES=(1)
declare -a RT_METHODS=("loess" "aligntree")
declare -a WEIGHT_RATIOS=(0.001 1 0.1 0.0001 0.00001)
declare -a PPM_CUTOFFS=(
    #2.5
    5
)
declare -a RT_CUTOFFS=(
    #50
    100
)

# Planned runs: anchor methods take the full radius cross-product, duplicate_only
# takes exactly one.
N_ANCHOR_METHODS=0
N_BFS_METHODS=0
N_WINDOW_METHODS=0
N_PLAIN_METHODS=0
for m in "${FILTER_METHODS[@]}"; do
    case "$m" in
        anchor_*) N_ANCHOR_METHODS=$((N_ANCHOR_METHODS + 1)) ;;
        bfs)      N_BFS_METHODS=$((N_BFS_METHODS + 1)) ;;
        window)   N_WINDOW_METHODS=$((N_WINDOW_METHODS + 1)) ;;
        *)        N_PLAIN_METHODS=$((N_PLAIN_METHODS + 1)) ;;
    esac
done
FIXED_COMBOS=$(( ${#FEATURE_MODES[@]} * ${#COORD_SOURCES[@]} * ${#MAX_MATCH_CANDIDATES[@]} \
                 * ${#RT_METHODS[@]} * ${#WEIGHT_RATIOS[@]} \
                 * ${#PPM_CUTOFFS[@]} * ${#RT_CUTOFFS[@]} ))
TOTAL_PLANNED=$(( FIXED_COMBOS * ( N_PLAIN_METHODS \
                 + N_ANCHOR_METHODS * ${#ANCHOR_RADIUS_MODES[@]} * ${#ANCHOR_RT_CUTOFFS[@]} \
                 + N_BFS_METHODS * ${#BFS_DEPTHS[@]} * ${#BFS_BUDGETS[@]} \
                 + N_WINDOW_METHODS * ${#WINDOW_SIZES[@]} ) ))

# ---- Apply the optional start delay ----------------------------------------
DELAY_SECS=$(awk "BEGIN{printf \"%d\", ($START_DELAY_HOURS) * 3600}")
if [ "$DELAY_SECS" -gt 0 ]; then
    echo "======================================================================"
    echo "Delaying start by ${START_DELAY_HOURS} h"
    echo "  now:      $(date)"
    echo "  starting: $(date -d "+${DELAY_SECS} seconds")"
    echo "  (${TOTAL_PLANNED} run(s) queued; Ctrl-C to cancel)"
    echo "======================================================================"
    sleep "$DELAY_SECS"
    echo "Delay elapsed at $(date) — starting."
fi

mkdir -p "$BASE_OUTDIR"
RESULTS_FILE="$BASE_OUTDIR/filter_modes_results.txt"
METRICS_FILE="$BASE_OUTDIR/filter_modes_metrics.tsv"
# One row per run, collected from each run's published pairwise_summary.json.
PAIRWISE_FILE="$BASE_OUTDIR/pairwise_metrics.tsv"

echo "======================================================================" > "$RESULTS_FILE"
echo "Filter Mode Parameter Grid Execution Log" >> "$RESULTS_FILE"
echo "======================================================================" >> "$RESULTS_FILE"
echo "Started: $(date)" >> "$RESULTS_FILE"
echo "Work directory: $(pwd)" >> "$RESULTS_FILE"
echo "Workflow: $WORKFLOW (anchor line, bin_anchor/)" >> "$RESULTS_FILE"
echo "Nextflow work dir: $NF_WORKDIR" >> "$RESULTS_FILE"
echo "Filter methods under test: ${FILTER_METHODS[*]}" >> "$RESULTS_FILE"
echo "Anchor radius modes: ${ANCHOR_RADIUS_MODES[*]}" >> "$RESULTS_FILE"
echo "Anchor RT cutoffs: ${ANCHOR_RT_CUTOFFS[*]}" >> "$RESULTS_FILE"
echo "BFS depths: ${BFS_DEPTHS[*]}" >> "$RESULTS_FILE"
echo "BFS budgets (ppm:rt): ${BFS_BUDGETS[*]}" >> "$RESULTS_FILE"
echo "Window sizes (ppm:rt): ${WINDOW_SIZES[*]}" >> "$RESULTS_FILE"
echo "Planned runs: $TOTAL_PLANNED" >> "$RESULTS_FILE"
echo "Nextflow resume: ${NF_RESUME:-<disabled>}" >> "$RESULTS_FILE"
echo "Steps enabled: 1$([ "$RUN_STEP1M" = "true" ] && echo ",metrics")$([ "$RUN_STEP2" = "true" ] && echo ",2")$([ "$RUN_PAIRWISE" = "true" ] && echo ",pairwise")$([ "$RUN_STEP3" = "true" ] && echo ",3")" >> "$RESULTS_FILE"
echo "======================================================================" >> "$RESULTS_FILE"
echo "" >> "$RESULTS_FILE"

# Metrics table header - one row per run, written by the metrics step
if [ "$RUN_STEP1M" = "true" ] && [ ! -f "$METRICS_FILE" ]; then
    printf 'run\tfilter_method\tradius\tanchor_rt\tcycles\tcomponents\tavg_linkage\tavg_linkage_ge2\tconflict_pct\tminority_pct\tpairwise_pct\tstar_pct\tworst_ppm_x\tworst_rt_x\tbeyond_2x\n' > "$METRICS_FILE"
fi

TOTAL_RUNS=0
SUCCESSFUL_RUNS=0
FAILED_RUNS=0

for filter_method in "${FILTER_METHODS[@]}"; do

# duplicate_only ignores the anchor settings, so collapse those loops to a single
# pass instead of running identical pipelines under different names.
# The two inner dimension lists carry whatever the active method actually varies:
# radius/RT for the anchor methods, depth/budget for BFS, and a single dummy pass for
# duplicate_only, which ignores both and would otherwise run identical pipelines under
# different names.
if [[ "$filter_method" == anchor_* ]]; then
    radius_list=("${ANCHOR_RADIUS_MODES[@]}")
    anchor_rt_list=("${ANCHOR_RT_CUTOFFS[@]}")
elif [[ "$filter_method" == bfs ]]; then
    radius_list=("${BFS_DEPTHS[@]}")
    anchor_rt_list=("${BFS_BUDGETS[@]}")
elif [[ "$filter_method" == window ]]; then
    radius_list=("n/a")
    anchor_rt_list=("${WINDOW_SIZES[@]}")
else
    radius_list=("n/a")
    anchor_rt_list=("n/a")
fi

for radius_mode in "${radius_list[@]}"; do
for anchor_rt in "${anchor_rt_list[@]}"; do
for feature_mode in "${FEATURE_MODES[@]}"; do
for coord_source in "${COORD_SOURCES[@]}"; do
for max_candidates in "${MAX_MATCH_CANDIDATES[@]}"; do
for method in "${RT_METHODS[@]}"; do
    for weight_ratio in "${WEIGHT_RATIOS[@]}"; do

        if [ "$weight_ratio" = "1.0" ] || [ "$weight_ratio" = "1" ]; then
            weight_label="unweighted"
        else
            weight_label="weighted_${weight_ratio}"
        fi

        for ppm in "${PPM_CUTOFFS[@]}"; do
        for rt_cutoff in "${RT_CUTOFFS[@]}"; do
        TOTAL_RUNS=$((TOTAL_RUNS + 1))
        RUN_START=$(date +%s)

        if [ "$max_candidates" = "1" ]; then cand_label="k1_nofix"; else cand_label="k${max_candidates}_fix"; fi
        if [[ "$filter_method" == anchor_* ]]; then
            filter_label="${filter_method}_${radius_mode}_art${anchor_rt}"
        elif [[ "$filter_method" == bfs ]]; then
            # anchor_rt carries "ppm:rt" here; ':' and '.' are awkward in paths
            filter_label="bfs_d${radius_mode}_b$(echo "$anchor_rt" | tr ':.' 'p_')"
        elif [[ "$filter_method" == window ]]; then
            filter_label="window_w$(echo "$anchor_rt" | tr ':.' 'p_')"
        else
            filter_label="${filter_method}"
        fi
        # filter label first so runs sort together per strategy
        run_name="filter_${filter_label}_mode_${feature_mode}_src${coord_source}_${cand_label}_${method}_${weight_label}_ppm${ppm}_rt${rt_cutoff}"
        run_dir="$BASE_OUTDIR/$run_name"

        echo ""
        echo "═══════════════════════════════════════════════════════════════"
        echo "Run ${TOTAL_RUNS}/${TOTAL_PLANNED}: $run_name"
        echo "Filter method:     $filter_method"
        if [[ "$filter_method" == window ]]; then
            echo "Window size:       ${anchor_rt%%:*} ppm x ${anchor_rt##*:} s (set-wise)"
        elif [[ "$filter_method" == bfs ]]; then
            echo "BFS max depth:     $radius_mode"
            echo "BFS budgets:       ${anchor_rt%%:*} ppm / ${anchor_rt##*:} s accumulated"
        else
            echo "Anchor radius:     $radius_mode"
            echo "Anchor RT cutoff:  $anchor_rt"
        fi
        echo "Output directory:  $run_dir"
        echo "═══════════════════════════════════════════════════════════════"

        clusters_outdir="$run_dir/clusters_final_report"
        proteobench_outdir="$run_dir/ProteoBench"

        # ========== STEP 1: the anchor workflow ==========
        echo ""
        echo "--- Step 1/${TOTAL_STEPS}: Running $WORKFLOW ---"

        cmd="nextflow run $WORKFLOW"
        cmd="$cmd -c $CONFIG"
        cmd="$cmd -w '$NF_WORKDIR'"
        [ -n "$NF_RESUME" ] && cmd="$cmd $NF_RESUME"
        cmd="$cmd --mzml_files '$MZML_FILES'"
        cmd="$cmd --tsv_files '$TSV_FILES'"
        cmd="$cmd --main_outdir '$run_dir'"
        cmd="$cmd --mmf_feature_mode $feature_mode"
        cmd="$cmd --mmf_rt_alignment_method $method"
        cmd="$cmd --mmf_mz_rt_weight_ratio $weight_ratio"
        cmd="$cmd --mmf_ppm $ppm"
        cmd="$cmd --mmf_rt_cutoff $rt_cutoff"
        cmd="$cmd --mmf_rt_alignment_mz_source $coord_source"
        cmd="$cmd --mmf_rt_alignment_rt_source $coord_source"
        cmd="$cmd --mmf_pairing_max_match_candidates $max_candidates"
        # The dimension under test
        cmd="$cmd --mmf_filter_method $filter_method"
        if [[ "$filter_method" == anchor_* ]]; then
            cmd="$cmd --mmf_anchor_radius_mode $radius_mode"
            cmd="$cmd --mmf_anchor_rt_cutoff $anchor_rt"
        elif [[ "$filter_method" == window ]]; then
            cmd="$cmd --mmf_window_ppm ${anchor_rt%%:*}"
            cmd="$cmd --mmf_window_rt ${anchor_rt##*:}"
        elif [[ "$filter_method" == bfs ]]; then
            cmd="$cmd --mmf_bfs_max_depth $radius_mode"
            cmd="$cmd --mmf_bfs_ppm_budget ${anchor_rt%%:*}"
            cmd="$cmd --mmf_bfs_rt_budget ${anchor_rt##*:}"
        fi

        echo "Command: $cmd" | tee -a "$RESULTS_FILE"
        echo "" | tee -a "$RESULTS_FILE"

        if ! eval "$cmd"; then
            FAILED_RUNS=$((FAILED_RUNS + 1))
            status="✗ FAILED at step 1 ($WORKFLOW)"
            echo "$run_name: $status" >> "$RESULTS_FILE"
            echo "$status"
            echo "" | tee -a "$RESULTS_FILE"
            echo "---" >> "$RESULTS_FILE"
            continue
        fi

        echo "Step 1 completed successfully"

        # The anchor workflow publishes here, NOT to feature_analysis_clusterless/
        analysis_dir="$run_dir/feature_analysis_clusterless_anchor"
        final_components="$analysis_dir/final_components_all_cycles.json"

        # ========== STEP 1M: component metrics ==========
        # The numbers this grid exists to compare: linkage size (want it LARGE),
        # pep_ident conflicts (want them LOW), and how far components drift beyond
        # the pairing cutoff.
        if [ "$RUN_STEP1M" = "true" ]; then
            echo ""
            echo "--- Step metrics/${TOTAL_STEPS}: component metrics ---"
            if [ ! -f "$final_components" ]; then
                echo "WARNING: $final_components not found, skipping metrics" | tee -a "$RESULTS_FILE"
            else
                python3 bin_anchor/summarize_filter_metrics.py \
                    --components "$final_components" \
                    --ppm "$ppm" \
                    --rt "$(if [[ "$filter_method" == anchor_* ]]; then echo "$anchor_rt"; \
                             else echo "$rt_cutoff"; fi)" \
                    --run "$run_name" \
                    --filter_method "$filter_method" \
                    --radius "$radius_mode" \
                    --anchor_rt "$anchor_rt" \
                    --tsv_append "$METRICS_FILE" | tee -a "$RESULTS_FILE" || \
                    echo "WARNING: metrics step failed (non-fatal)" | tee -a "$RESULTS_FILE"
            fi
        fi

        if [ "$RUN_STEP2" = "true" ]; then
            # ========== STEP 2: process_clusters_consensusxml.nf ==========
            echo ""
            echo "--- Step 2/${TOTAL_STEPS}: Running process_clusters_consensusxml.nf ---"

            clusters_consensus_xml="$analysis_dir/components.consensusXML"
            original_consensus_xml="results/quantifications/features_with_annotated_identifications/consensus________0.01_fdr.consensusXML"
            feature_tsvs_dir="results/quantifications/features_with_annotated_identifications"

            if [ ! -f "$clusters_consensus_xml" ]; then
                FAILED_RUNS=$((FAILED_RUNS + 1))
                status="✗ FAILED at step 2 (components.consensusXML not found)"
                echo "$run_name: $status" >> "$RESULTS_FILE"
                echo "$status"
                echo "Expected file: $clusters_consensus_xml" | tee -a "$RESULTS_FILE"
                echo "" | tee -a "$RESULTS_FILE"
                echo "---" >> "$RESULTS_FILE"
                continue
            fi

            cmd2="nextflow run process_clusters_consensusxml.nf"
            cmd2="$cmd2 --pcc_clusters_consensus_xml '$clusters_consensus_xml'"
            cmd2="$cmd2 --pcc_original_consensus_xml '$original_consensus_xml'"
            cmd2="$cmd2 --pcc_feature_tsvs_dir '$feature_tsvs_dir'"
            cmd2="$cmd2 --pcc_outdir '$clusters_outdir'"
            cmd2="$cmd2 --pcc_normalize $RUN_NORMALIZATION"
            cmd2="$cmd2 --pcc_run_ttest $RUN_TTEST"
            if [ "$RUN_NORMALIZATION" = "true" ]; then
                cmd2="$cmd2 --pcc_sample_conditions '$(realpath "$SAMPLE_CONDITIONS")'"
                cmd2="$cmd2 --pcc_normalization_method $NORMALIZATION_METHOD"
                cmd2="$cmd2 --pcc_r_container '$R_CONTAINER'"
            fi
            if [ "$RUN_TTEST" = "true" ] && [ -n "$TTEST_COMPARISONS" ]; then
                cmd2="$cmd2 --pcc_ttest_comparisons '$TTEST_COMPARISONS'"
            fi
            # Pairwise analysis: hand step 2 the components JSON from step 1.
            # That file is produced by the mapping workflow, not by step 2, so
            # without this the analysis stays switched off.
            if [ "$RUN_PAIRWISE" = "true" ] && [ -f "$final_components" ]; then
                cmd2="$cmd2 --pcc_components_json '$(realpath "$final_components")'"
                cmd2="$cmd2 --pcc_pairwise_theme '$PAIRWISE_THEME'"
                # Chart against the cutoffs this run actually used.
                cmd2="$cmd2 --pcc_pairwise_ppm $ppm"
                cmd2="$cmd2 --pcc_pairwise_rt $rt_cutoff"
            elif [ "$RUN_PAIRWISE" = "true" ]; then
                echo "WARNING: $final_components not found, skipping pairwise analysis" \
                    | tee -a "$RESULTS_FILE"
            fi

            echo "Command: $cmd2" | tee -a "$RESULTS_FILE"
            echo "" | tee -a "$RESULTS_FILE"

            if ! eval "$cmd2"; then
                FAILED_RUNS=$((FAILED_RUNS + 1))
                status="✗ FAILED at step 2 (process_clusters_consensusxml.nf)"
                echo "$run_name: $status" >> "$RESULTS_FILE"
                echo "$status"
                echo "" | tee -a "$RESULTS_FILE"
                echo "---" >> "$RESULTS_FILE"
                continue
            fi

            echo "Step 2 completed successfully"

            # Collect this run's pairwise numbers into the grid-level table.
            # The process writes only inside its work dir (so it stays correct
            # under Docker and -resume); aggregating the published summary is
            # the caller's job.
            if [ "$RUN_PAIRWISE" = "true" ]; then
                pairwise_summary="$clusters_outdir/pairwise_components/pairwise_summary.json"
                if [ -f "$pairwise_summary" ]; then
                    python3 - "$pairwise_summary" "$PAIRWISE_FILE" "$run_name" <<'PYEOF' \
                        | tee -a "$RESULTS_FILE"
import json, os, sys
summary_path, tsv_path, run = sys.argv[1:4]
with open(summary_path) as fh:
    m = json.load(fh)
cols = ["run", "components", "pairs", "ppm_cutoff", "rt_cutoff",
        "ppm_exceed", "ppm_exceed_pct", "rt_exceed", "rt_exceed_pct",
        "mz_median", "mz_p99", "mz_max", "rt_median", "rt_p99", "rt_max",
        "size2_ppm_exceed", "size2_rt_exceed"]
m["run"] = run
new = not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0
with open(tsv_path, "a") as fh:
    if new:
        fh.write("\t".join(cols) + "\n")
    fh.write("\t".join(
        f"{m[c]:.4f}" if isinstance(m[c], float) else str(m[c]) for c in cols) + "\n")
print(f"  pairwise: {m['ppm_exceed']:,} pairs ({m['ppm_exceed_pct']:.4f}%) over "
      f"{m['ppm_cutoff']:g} ppm, {m['rt_exceed']:,} ({m['rt_exceed_pct']:.4f}%) over "
      f"{m['rt_cutoff']:g} s")
PYEOF
                else
                    echo "WARNING: $pairwise_summary not found (pairwise analysis did not run)" \
                        | tee -a "$RESULTS_FILE"
                fi
            fi
        fi  # RUN_STEP2

        if [ "$RUN_STEP3" = "true" ]; then
            # ========== STEP 3: ProteoBench conversion ==========
            echo ""
            echo "--- Step 3/${TOTAL_STEPS}: Converting to ProteoBench format ---"

            input_tsv="$clusters_outdir/clusters_raw_quantification_with_identifications.tsv"

            if [ ! -f "$input_tsv" ]; then
                FAILED_RUNS=$((FAILED_RUNS + 1))
                status="✗ FAILED at step 3 (input TSV not found)"
                echo "$run_name: $status" >> "$RESULTS_FILE"
                echo "$status"
                echo "Expected file: $input_tsv" | tee -a "$RESULTS_FILE"
                echo "" | tee -a "$RESULTS_FILE"
                echo "---" >> "$RESULTS_FILE"
                continue
            fi

            mkdir -p "$proteobench_outdir"

            cmd3="python3 convert_unbequant_to_proteobench_custom.py"
            cmd3="$cmd3 -i '$input_tsv'"
            cmd3="$cmd3 -f files/ProteoBenchFASTA_MixedSpecies_HYE.fasta"
            cmd3="$cmd3 -o '$proteobench_outdir/ProteoBench.tsv'"

            echo "Command: $cmd3" | tee -a "$RESULTS_FILE"
            echo "" | tee -a "$RESULTS_FILE"

            if ! eval "$cmd3"; then
                FAILED_RUNS=$((FAILED_RUNS + 1))
                status="✗ FAILED at step 3 (ProteoBench conversion)"
                echo "$run_name: $status" >> "$RESULTS_FILE"
                echo "$status"
                echo "" | tee -a "$RESULTS_FILE"
                echo "---" >> "$RESULTS_FILE"
                continue
            fi

            echo "Step 3 completed successfully"
        fi  # RUN_STEP3

        RUN_ELAPSED=$(( $(date +%s) - RUN_START ))
        SUCCESSFUL_RUNS=$((SUCCESSFUL_RUNS + 1))
        status="✓ SUCCESS (${TOTAL_STEPS} step(s) completed in ${RUN_ELAPSED}s)"
        echo "$run_name: $status" >> "$RESULTS_FILE"
        echo "$status"
        echo "Report: $analysis_dir/clusterless_report.html" | tee -a "$RESULTS_FILE"

        echo "" | tee -a "$RESULTS_FILE"
        echo "---" >> "$RESULTS_FILE"

        done  # rt_cutoff
        done  # ppm
    done      # weight_ratio
done          # method
done          # max_candidates
done          # coord_source
done          # feature_mode
done          # anchor_rt
done          # radius_mode
done          # filter_method

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
echo "Full results saved to: $RESULTS_FILE"

if [ "$RUN_STEP1M" = "true" ] && [ -f "$METRICS_FILE" ]; then
    echo ""
    echo "--- Filter mode comparison ---"
    column -t -s $'\t' "$METRICS_FILE" 2>/dev/null || cat "$METRICS_FILE"
    echo ""
    echo "Metrics table: $METRICS_FILE"
fi

if [ "$RUN_PAIRWISE" = "true" ] && [ -f "$PAIRWISE_FILE" ]; then
    echo ""
    echo "--- Pairwise within-component distances ---"
    column -t -s $'\t' "$PAIRWISE_FILE" 2>/dev/null || cat "$PAIRWISE_FILE"
    echo ""
    echo "Pairwise table: $PAIRWISE_FILE"
    echo "Per-run figures: <run_dir>/clusters_final_report/pairwise_components/"
fi

echo ""
echo "--- Generating grid comparison summary ---"
python3 bin/summarize_grid_results.py --results_dir "$BASE_OUTDIR" || \
    echo "(grid summary skipped - it expects the feature-mode layout)"
