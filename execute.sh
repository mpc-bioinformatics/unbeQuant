nextflow run main_workflow_protein_fasta.nf \
    -resume \
    --main_fasta_file ../mus_musculus_UP000000589_2024_02_13.fasta \
    --main_raw_files_folder ../raws \
    --main_comet_params ../Goldstandard_search_params.txt\
    --ctm_num_procs_conversion 8 \
    --idc_num_parallel_threads_per_search 6 \
    --qal_limit_num_of_parallel_feature_finders 6
