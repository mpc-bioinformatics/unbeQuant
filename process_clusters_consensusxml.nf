#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Standalone Nextflow workflow for processing cluster-filtered consensusXML files
 * Uses existing processes from quantify_and_align.nf to generate final reports
 */

// Required Parameters
params.pcc_clusters_consensus_xml = "results/feature_analysis/feature_data_lists/clusters_recovered.consensusXML"  // Cluster-recovered consensusXML (includes recovered single-vertex clusters)
params.pcc_original_consensus_xml = "results/quantifications/features_with_annotated_identifications/consensus________0.01_fdr.consensusXML"  // Original consensus for reference
params.pcc_feature_tsvs_dir = "results/quantifications/features_with_annotated_identifications"  // Directory containing feature TSV files
params.pcc_outdir = "results/feature_analysis/clusters_final_report"  // Output directory

// Start-from-existing-tables bypass. When set to a directory that already holds
// the report tables (*raw_quantification_with_identifications.tsv, *_reduced.tsv,
// *_only_intensities_and_ids.tsv), the consensusXML -> table generation, the xlsx
// reports and the plots are all SKIPPED, and only normalization + t-test +
// spike-in QC run on those existing tables. Useful for re-running the downstream
// analysis on the original unbequant output without regenerating the (huge)
// tables. null = normal mode (generate tables from the consensusXML).
params.pcc_precomputed_report_dir = null
params.pcc_minlh = 7  // Minimum number of MS1 scans (needed for visualizations)

// Fixed FDR label for cluster processing
params.pcc_fdr_label = "filtered_clusters"

// ---------------------------------------------------------------------------
// Pairwise within-component distance analysis
// ---------------------------------------------------------------------------
// Measures EVERY pair inside each component (a 4-member component contributes
// 6 pairs) and charts those distances against the pairing cutoffs. Components
// are connected components, so membership is transitive while the cutoffs are
// pairwise: A-B and B-C both passing does not keep A-C inside the cutoff, and
// nothing downstream re-checks. This quantifies that drift.
//
// The components JSON is produced upstream by map_mzml_features_v*.nf, not by
// this workflow, so it is an optional extra input: null = skip the analysis.
params.pcc_components_json = null
// Cutoffs default to the ones recorded in the run's own metadata. Set these
// only to test against different values than the run used.
params.pcc_pairwise_ppm = null
params.pcc_pairwise_rt = null
// light | dark | both
params.pcc_pairwise_theme = "both"
// Per-run numbers are published as pairwise_summary.json next to the figures.
// Aggregating those across a grid is the caller's job: this process writes only
// inside its work directory, so it stays correct under Docker and -resume.

// ---------------------------------------------------------------------------
// Normalization (ProtStatsWF::workflow_QC) + differential analysis
// ---------------------------------------------------------------------------
params.pcc_normalize = true  // Run intensity normalization
// Maps each measurement file to an experimental condition. Required for
// normalization: workflow_QC derives sample groups from the intensity column
// name prefix, so the columns must be renamed <CONDITION>_<REPLICATE>.
params.pcc_sample_conditions = "$projectDir/files/sample_conditions.tsv"
// Container carrying R + ProtStatsWF (built from Dockerfile.unbequant_r)
params.pcc_r_container = "unbequant:r"
// workflow_QC normalization method: loess | median | quantile | lts | nonorm
params.pcc_normalization_method = "loess"
// MA-plots are a QC diagnostic and scale as n(n-1)/2 pairs. For 25 samples that
// is 300 pages (~11 min, ~1 GB), so they are off by default for routine runs.
params.pcc_generate_maplots = false

params.pcc_run_ttest = true  // Run t-test / volcano analysis on normalized data
// Condition pairs to compare, e.g. ["C1:C2", "C1:C5"]. Empty = all pairwise
// combinations found in the sample conditions file.
params.pcc_ttest_comparisons = []
// Fraction of missing values tolerated per group (0 - 1). The underlying script
// defaults to 0, which would discard any feature not seen in every replicate.
params.pcc_ttest_missingness = 0.3

// ---------------------------------------------------------------------------
// Spike-in intensity QC (needs normalization)
// ---------------------------------------------------------------------------
params.pcc_spike_in_qc = true  // Build the per-condition spike-in intensity QC report
// FASTA of the spike-in proteins (one boxplot section is produced per condition,
// one boxplot per spike-in protein found here).
params.pcc_spikeins_fasta = "$projectDir/files/spikeins.txt"

// ---------------------------------------------------------------------------
// HELPERS
// ---------------------------------------------------------------------------

// Params supplied on the command line arrive as Strings, and a non-empty String
// is truthy in Groovy - so `--pcc_normalize false` would otherwise enable it.
def asBool(value) {
    if (value == null) return false
    if (value instanceof Boolean) return value
    return value.toString().toLowerCase() in ["true", "yes", "1", "on"]
}

// Read the distinct conditions from the sample conditions file, preserving the
// order in which they first appear.
def readConditions(condition_file) {
    def conditions = []
    file(condition_file).eachLine { raw ->
        def line = raw.trim()
        if (!line || line.startsWith("#")) return
        def parts = line.split("\t")
        if (parts.size() < 2) parts = line.split(",")
        if (parts.size() < 2) return
        def fname = parts[0].trim()
        def cond  = parts[1].trim()
        if (fname.toLowerCase() in ["filename", "file", "sample"]) return
        if (cond && !conditions.contains(cond)) conditions << cond
    }
    return conditions
}

// Comparisons to run: either the explicit list, or every pairwise combination.
def resolveComparisons(condition_file) {
    def requested = params.pcc_ttest_comparisons
    if (requested) {
        // From the command line this is a String ("C1:C2,C1:C5"); from a config
        // or params file it is already a List.
        if (requested instanceof CharSequence) {
            return requested.toString()
                            .split(/[,\s]+/)
                            .findAll { it }
                            .collect { it.toString() }
        }
        return (requested as List).collect { it.toString() }
    }
    def conditions = readConditions(condition_file)
    def pairs = []
    conditions.eachWithIndex { first, idx ->
        conditions.drop(idx + 1).each { second ->
            pairs << "${first}:${second}".toString()
        }
    }
    return pairs
}

// Locate exactly one report table in a directory by filename suffix. Uses plain
// java.io.File so it resolves at DAG-construction time (before any process runs).
def pickTable(dirPath, suffix) {
    def matches = new File(dirPath.toString()).listFiles()
                      ?.findAll { it.name.endsWith(suffix) } ?: []
    if (matches.size() != 1) {
        error("Expected exactly one '*${suffix}' in ${dirPath}, found: ${matches*.name}")
    }
    return file(matches[0].toString())
}

// Validate inputs. The consensusXML is only needed when we actually generate the
// tables; in start-from-existing-tables mode it is irrelevant.
if (!params.pcc_precomputed_report_dir) {
    if (!file(params.pcc_clusters_consensus_xml).exists()) {
        error("ERROR: Cluster consensusXML file not found: ${params.pcc_clusters_consensus_xml}")
    }
    if (!file(params.pcc_original_consensus_xml).exists()) {
        error("ERROR: Original consensusXML file not found: ${params.pcc_original_consensus_xml}")
    }
}

// Standalone MAIN Workflow
workflow {
    if (params.pcc_precomputed_report_dir) {
        // Start-from-existing-tables mode: no consensusXML / feature TSVs needed.
        // The downstream (normalization + t-test + spike-in QC) reads the tables
        // directly, so pass placeholders for the unused consensus/tsv slots.
        consensus_input = Channel.of(
            tuple(params.pcc_fdr_label, file("${projectDir}/nextflow.config"),
                  file("${projectDir}/nextflow.config"), [])
        )
    } else {
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
    }

    // Execute cluster processing workflow
    process_clusters_consensus(consensus_input)
}

// Importable Workflow
workflow process_clusters_consensus {
    take:
        input_data  // tuple(fdr_label, clusters_consensus, original_consensus, [feature_tsvs])
    
    main:
        // The report tables feeding the downstream analysis. Either loaded from an
        // existing directory (start-from-existing-tables mode) or generated from
        // the cluster consensusXML. Tuple layout matches
        // generate_feature_ident_intesity_table.out: (fdr, full, reduced, minimal).
        tables_ch  = Channel.empty()
        xlsx_ch    = Channel.empty()
        plot_ch    = Channel.empty()

        if (params.pcc_precomputed_report_dir) {
            def dir = params.pcc_precomputed_report_dir
            log.info "Start-from-existing-tables mode: reading report tables from ${dir}"
            // normalize_quantification_table (the only downstream table consumer)
            // reads ONLY the reduced table; the full (multi-GB) and minimal tables
            // feed only the xlsx/plots, which are skipped here. Stage tiny distinct
            // placeholders for those two slots so the huge full table is not copied.
            def full_placeholder    = file("${projectDir}/nextflow.config")
            def minimal_placeholder = file("${projectDir}/nextflow_standalone.config")
            tables_ch = Channel.of(tuple(
                params.pcc_fdr_label,
                full_placeholder,
                pickTable(dir, "quantification_with_identifications_reduced.tsv"),
                minimal_placeholder,
            ))
            // xlsx reports and plots already exist for these tables - skip them.
        } else {
            // Map to extract only fdr_label, clusters_consensus, and feature_tsvs for the process
            // (original_consensus is kept for future reference/comparison)
            input_data.map { fdr, clusters_cons, orig_cons, tsvs ->
                tuple(fdr, clusters_cons, tsvs)
            }.set { consensus_tsvs_input }

            // Generate TSV tables from cluster consensusXML and feature data
            generate_feature_ident_intesity_table(consensus_tsvs_input)
            tables_ch = generate_feature_ident_intesity_table.out

            // Generate Excel reports from the TSV tables
            generate_xlsx_reports_from_tables(tables_ch)
            xlsx_ch = generate_xlsx_reports_from_tables.out

            // Generate visualization plots from TSV data
            generate_plots_from_tables(tables_ch)
            plot_ch = generate_plots_from_tables.out
        }

        // ── Pairwise within-component distance analysis ──────────────────────
        // Independent of the table/normalization chain: it reads the components
        // JSON from the upstream mapping workflow, not the consensusXML.
        pairwise_ch = Channel.empty()
        if (params.pcc_components_json) {
            components_json = file(params.pcc_components_json, checkIfExists: true)
            pairwise_component_analysis(components_json)
            pairwise_ch = pairwise_component_analysis.out.figures
        }

        // ── Normalization + differential analysis ────────────────────────────
        normalized_ch    = Channel.empty()
        ttest_ch         = Channel.empty()
        spike_in_qc_ch   = Channel.empty()

        if (asBool(params.pcc_normalize)) {
            conditions_file = file(params.pcc_sample_conditions, checkIfExists: true)

            normalize_quantification_table(
                tables_ch,
                conditions_file
            )
            normalized_ch = normalize_quantification_table.out.normalized

            if (asBool(params.pcc_run_ttest)) {
                comparisons = resolveComparisons(params.pcc_sample_conditions)
                if (comparisons) {
                    log.info "Running t-tests for comparisons: ${comparisons.join(', ')}"
                    ttest_input = normalized_ch.combine(Channel.fromList(comparisons))
                    ttest_analysis(ttest_input)
                    ttest_ch = ttest_analysis.out
                } else {
                    log.warn "No condition pairs resolved — skipping t-test analysis"
                }
            }

            // Spike-in intensity QC: per-condition boxplots + relative-mean table.
            // Consumes the same normalized table as the t-test.
            if (asBool(params.pcc_spike_in_qc)) {
                spikeins_fasta = file(params.pcc_spikeins_fasta, checkIfExists: true)
                spike_in_intensity_qc(normalized_ch, spikeins_fasta)
                spike_in_qc_ch = spike_in_intensity_qc.out
            }
        }

    emit:
        tsv_outputs = tables_ch
        xlsx_outputs = xlsx_ch
        plot_outputs = plot_ch
        normalized_outputs = normalized_ch
        ttest_outputs = ttest_ch
        spike_in_qc_outputs = spike_in_qc_ch
        pairwise_outputs = pairwise_ch
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


// Normalize the quantification intensities with ProtStatsWF::workflow_QC.
//
// Runs in three steps:
//   1. prepare_qc_input.py          rename intensity columns to <CONDITION>_<NN>
//                                   (workflow_QC derives groups from the name)
//   2. normalize_intensities_qc.R   workflow_QC -> D_norm_ID.tsv (log2, normalized)
//   3. merge_normalized_intensities.py
//                                   put the normalized values back into the full
//                                   unbeQuant table
//
// The published result is the MERGED table: the complete quantification table
// with normalized intensities. D_norm_ID.tsv (intensities only) is kept as a
// consumable intermediate.
process normalize_quantification_table {
    tag "Normalize (${params.pcc_normalization_method})"
    publishDir "${params.pcc_outdir}/normalization", mode:'copy'
    container "${params.pcc_r_container}"
    memory "16G"
    stageInMode "copy"

    input:
    tuple val(fdr), file(full_file), file(reduced_file), file(minimal_file)
    file(sample_conditions)

    output:
    tuple val(fdr), path("clusters_quantification_normalized.tsv"), path("clusters_quantification_normalized_group_columns.json"), emit: normalized
    path("D_norm_ID.tsv"), emit: normalized_intensities
    path("qc_input.tsv"), emit: qc_input
    path("qc_input_column_map.json"), emit: column_map
    path("*.png"), optional: true, emit: qc_plots
    path("*.csv"), optional: true, emit: qc_tables
    // MA-plots are written as a single multi-page PDF (only when
    // pcc_generate_maplots is enabled), so they need their own pattern.
    path("*.pdf"), optional: true, emit: ma_plots

    script:
    def ma_plots = asBool(params.pcc_generate_maplots) ? "TRUE" : "FALSE"
    """
    # 1. Rename intensity columns to the <CONDITION>_<REPLICATE> form workflow_QC needs
    python3 \$(which prepare_qc_input.py) \
        --input_tsv ${reduced_file} \
        --sample_conditions ${sample_conditions} \
        --output_tsv qc_input.tsv \
        --output_column_map qc_input_column_map.json

    # The intensity block is contiguous; read its range from the column map so
    # this adapts to any number of measurement files.
    INT_START=\$(python3 -c "import json;print(json.load(open('qc_input_column_map.json'))['intensity_column_start'])")
    INT_END=\$(python3 -c "import json;print(json.load(open('qc_input_column_map.json'))['intensity_column_end'])")

    # 2. Normalize (log2 + ${params.pcc_normalization_method})
    Rscript \$(which normalize_intensities_qc.R) \
        --data_path=qc_input.tsv \
        --output_path=. \
        --intensity_col_start=\$INT_START \
        --intensity_col_end=\$INT_END \
        --normalization_method=${params.pcc_normalization_method} \
        --generate_MAplots=${ma_plots}

    # 3. Merge normalized intensities back into the full table
    python3 \$(which merge_normalized_intensities.py) \
        --reduced_tsv ${reduced_file} \
        --normalized_tsv D_norm_ID.tsv \
        --column_map qc_input_column_map.json \
        --output_tsv clusters_quantification_normalized.tsv \
        --output_group_columns clusters_quantification_normalized_group_columns.json
    """
}


// Two-group t-test + volcano plot on the normalized quantification table.
//
// The comparison is given as "<conditionA>:<conditionB>"; the actual column
// headers for each group are expanded from the group-columns JSON written by
// merge_normalized_intensities.py.
process ttest_analysis {
    tag "t-test ${comparison}"
    publishDir "${params.pcc_outdir}/ttest", mode:'copy'
    // Same image as the normalization step, which carries plotly/numpy; the
    // script is invoked from the live workspace bin/ (projectDir is mounted) so
    // edits take effect without a container rebuild.
    container "${params.pcc_r_container}"
    memory "8G"
    stageInMode "copy"

    input:
    tuple val(fdr), path(normalized_tsv), path(group_columns), val(comparison)

    output:
    tuple val(fdr), val(comparison), path("ttest_*.tsv"), path("volcano_*.html")

    script:
    def (groupA, groupB) = comparison.tokenize(':')
    def label = "${groupA}_vs_${groupB}"
    """
    # Expand each condition into its -gA/-gB column-header arguments
    GA=\$(python3 -c "import json;print(' '.join('-gA '+c for c in json.load(open('${group_columns}'))['${groupA}']))")
    GB=\$(python3 -c "import json;print(' '.join('-gB '+c for c in json.load(open('${group_columns}'))['${groupB}']))")

    if [ -z "\$GA" ] || [ -z "\$GB" ]; then
        echo "ERROR: no columns found for '${groupA}' or '${groupB}' in ${group_columns}" >&2
        exit 1
    fi

    # -lip: the merged table carries log2 values, so fold changes are computed
    # by un-logging before averaging.
    python3 ${projectDir}/bin/unbequant_ttest_and_plotting.py \
        -in ${normalized_tsv} \
        \$GA \$GB \
        -p ${params.pcc_ttest_missingness} \
        -lip 1 \
        -out ttest_${label}.tsv \
        -volcano volcano_${label}.html
    """
}


// Per-condition spike-in intensity QC: one boxplot section per condition (one
// box per file, one panel per spike-in protein) plus a relative-mean table
// across conditions. Reads the normalized (log2) quantification table.
process spike_in_intensity_qc {
    tag "spike-in intensity QC"
    publishDir "${params.pcc_outdir}/spike_in_qc", mode: 'copy'
    // luxii/unbequant:latest carries plotly/numpy; the script is invoked from the
    // live workspace bin/ (projectDir is mounted) so edits take effect without a
    // container rebuild.
    container "luxii/unbequant:latest"
    memory "8G"
    stageInMode "copy"

    input:
    tuple val(fdr), path(normalized_tsv), path(group_columns)
    path spikeins_fasta

    output:
    tuple val(fdr), path("spike_in_intensity_qc.html"), path("spike_in_relative_means.tsv")

    script:
    """
    python3 ${projectDir}/bin/spike_in_intensity_qc.py \\
        --input_tsv ${normalized_tsv} \\
        --spikeins_fasta ${spikeins_fasta} \\
        --group_columns ${group_columns} \\
        --output_html spike_in_intensity_qc.html \\
        --output_table spike_in_relative_means.tsv \\
        --log2_input 1
    """
}


// Pairwise within-component distance charts + summary.
//
// Reads final_components_all_cycles.json (from the upstream mapping workflow)
// and writes, per theme: the |Delta m/z| and |Delta RT| histograms, and the two
// cutoff-conformance figures broken down by component size. Cutoffs come from
// the run's own metadata unless overridden, so a 5 ppm run is never charted
// against a 10 ppm line.
//
// luxii/unbequant:latest already carries matplotlib + numpy; the script is
// invoked from the live workspace bin/ (projectDir is mounted) so edits take
// effect without a container rebuild, matching spike_in_intensity_qc above.
process pairwise_component_analysis {
    tag "pairwise component distances"
    publishDir "${params.pcc_outdir}/pairwise_components", mode: 'copy'
    container "luxii/unbequant:latest"
    memory "8G"
    // No stageInMode 'copy' here (unlike the table processes): the components
    // JSON runs to ~700 MB and is only read, so the default symlink avoids a
    // full copy per run.

    input:
    path components_json

    output:
    path("*.png"), emit: figures
    path("pairwise_summary.json"), emit: summary

    script:
    def ppm_arg = params.pcc_pairwise_ppm ? "--ppm ${params.pcc_pairwise_ppm}" : ""
    def rt_arg  = params.pcc_pairwise_rt  ? "--rt ${params.pcc_pairwise_rt}"   : ""
    """
    PYTHONUNBUFFERED=1 python3 ${projectDir}/bin/analyze_pairwise_components.py \\
        --components_json ${components_json} \\
        --outdir . \\
        --theme ${params.pcc_pairwise_theme} \\
        --run '${params.pcc_fdr_label}' \\
        ${ppm_arg} ${rt_arg}
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
