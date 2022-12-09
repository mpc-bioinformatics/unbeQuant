#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.identification_folder = "$PWD/results"  // Input-Directory of the tsvs statistics folder. E.G.: "ms2_specific_fasta_idents"
params.outdir = "$PWD/results" // Output-Directory of the tsvs statistics folder. E.G.: "ms2_specific_fasta_idents"
params.fdr = "0.01" // The same FDR-Value used as in all other workflows (used to cut the filenames only)


workflow{
    // Execute for each group of files individually
    unique()
    shared()
    unique_ft()
    shared_ft()
    summary()
}


workflow unique {
    unique = Channel.fromPath(params.identification_folder + "/tsvs_psms_statistics/*ident_count_unique.tsv")
    create_summary_tsv(unique.toList(), Channel.of("unique_psms.tsv"), unique.count())
    summarize_psms_to_excel(create_summary_tsv.out, Channel.of("unique_psms_heatmap.xlsx"))
}


workflow shared {
    shared = Channel.fromPath(params.identification_folder + "/tsvs_psms_statistics/*ident_count_shared.tsv")
    create_summary_tsv(shared.toList(), Channel.of("shared_psms.tsv"), shared.count())
    summarize_psms_to_excel(create_summary_tsv.out, Channel.of("shared_psms_heatmap.xlsx"))
}


workflow unique_ft{
    unique_ft = Channel.fromPath(params.identification_folder + "/tsvs_psms_statistics/*ident_count_unique_with_features.tsv")
    create_summary_tsv(unique_ft.toList(), Channel.of("unique_feature_psms.tsv"), unique_ft.count())
    summarize_psms_to_excel(create_summary_tsv.out, Channel.of("unique_feature_psms_heatmap.xlsx"))
}


workflow shared_ft{
    shared_ft = Channel.fromPath(params.identification_folder + "/tsvs_psms_statistics/*ident_count_shared_with_only_features.tsv")
    create_summary_tsv(shared_ft.toList(), Channel.of("shared_only_feature_psms.tsv"), shared_ft.count())
    summarize_psms_to_excel(create_summary_tsv.out, Channel.of("shared_only_feature_psms_heatmap.xlsx"))
}


workflow summary {
    summaries = Channel.fromPath(params.identification_folder + "/tsvs_psms_statistics/*summary.tsv")
    create_summary_tsv(summaries.toList(), Channel.of("psms_summary.tsv"), summaries.count())
}


process create_summary_tsv {
    publishDir "${params.outdir}", mode:'copy'
    debug true

    input:
    file(tsv_file)
    val(file_name)
    val(num_files)

    output: 
    file "${params.fdr}_fdr_${file_name}"

    when:
        num_files != 0

    """
    { for filename in ${tsv_file}; do head -n 1 \$filename | while read line; do echo -e "source_file\t\${line}"; done; break; done;
    for file in ${tsv_file}
    do
    filename=\${file%_${params.fdr}fdr*}
    tail -n +2 \$file | while read line; do echo -e "\$filename\t\${line}"; done
    done
    } > ${params.fdr}_fdr_${file_name}
    """  
}


process summarize_psms_to_excel {
    publishDir "${params.outdir}", mode:'copy'

    input:
    file(tsv_file)
    val(output_tsv_file)

    output: 
    file "${params.fdr}_fdr_${output_tsv_file}"

    """
    export_excel_from_found_psms.py ${tsv_file} ${params.fdr}_fdr_${output_tsv_file}
    """  
}
