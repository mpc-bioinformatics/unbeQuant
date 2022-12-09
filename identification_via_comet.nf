#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.mgf_folder = "$PWD/MGFs"  // Input-Directory of MGFs, which should be used for identification
params.fasta_file = "peptides.fasta"  // Database (FASTA-file) used for identification WITHOUT DECOYS!
params.fasta_not_from_protgraph = false // Flag if fasta is from protgraph (peptides-fasta) or not
params.statistics_remove_varmods = true // Flag if peptides with variable PTMs are counted as unique
params.statistics_unique_proteins = true // Flag if proteins are counted as unique (in case of multiple sequences in the same protein)

// Results output dir
params.outdir = "$PWD/results"  // Output-Directory of the FASTA-results

// Optional, but should be set!
params.search_parameter_file = "$PWD/example_configurations/comet_config.txt"  //Search Parameters for Comet
params.tda = 1 // 0 --> No Target-Decoy appraoch | 1 --> Target-Decoy appraoch
params.use_percolator = 0 // 0 --> Do not use, 1 --> Use percolator instead of the coment PSMSs results TODO Default should be 0
params.fdr = 0.05 // If TDA is set to 1, also return the pruned tsvs, whit up to fdr. (originals are also kept)
params.num_parallel_searches = Runtime.runtime.availableProcessors()

workflow {
    // Get all mzML files which should be identified
    mgfs = Channel.fromPath(params.mgf_folder + "/*.mgf")

    // Get FASTA-file
    fasta_file = Channel.fromPath(params.fasta_file)

    // Get Modification Parameters file
    modifications_file = Channel.fromPath(params.search_parameter_file)

    // Combined channel replicated the indexed fasta for each MZML to be reused
    combined_channel = fasta_file
        .combine(modifications_file)
        .combine(mgfs)
        
    // Start search
    comet_search_mgf(combined_channel)

    // if tda, then also return the pruned results
    if (params.tda == 1) {
        if (params.use_percolator == 1) {
            execute_percolator(comet_search_mgf.out[2])
            combined_txt_channel = execute_percolator.out
                .combine(fasta_file)
        } else {
            combined_txt_channel = comet_search_mgf.out[1]
                .combine(fasta_file)
        }

        comet_txt_to_tsv_with_fdr(combined_txt_channel)
        count_and_seperate_psms(comet_txt_to_tsv_with_fdr.out[1], comet_txt_to_tsv_with_fdr.out[2])
    }
}

process comet_search_mgf {
    maxForks params.num_parallel_searches
    publishDir "${params.outdir}/idents", mode:'copy'

    input:
    tuple file(input_fasta), file(mod_file), file(mgf_file)

    output: 
    file "${mgf_file.baseName}.mzid"
    file "${mgf_file.baseName}.txt"
    file "${mgf_file.baseName}.pin"

    """
    sed 's/^decoy_search.*/decoy_search = ${params.tda} /' ${mod_file} > ${mod_file.baseName}_new.txt
    sed -i 's/^output_mzidentmlfile.*/output_mzidentmlfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^output_txtfile.*/output_txtfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^output_percolatorfile.*/output_percolatorfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^decoy_prefix.*/decoy_prefix = DECOY_/' ${mod_file.baseName}_new.txt


    comet.linux_v2022.01.2.exe -P${mod_file.baseName}_new.txt -D${input_fasta} ${mgf_file}
    """  
}


process execute_percolator {
    maxForks params.num_parallel_searches
    publishDir "${params.outdir}/idents", mode:'copy'

    input:
    file(comet_pin)

    output: 
    file "${comet_pin.baseName}.tsv"

    """
    percolator3.06 --trainFDR ${params.fdr} --testFDR ${params.fdr} -P DECOY_ --only-psms --decoy-xml-output -X ${comet_pin.baseName}.xml ${comet_pin}
    percout_to_tsv.py ${comet_pin.baseName}.xml ${comet_pin.baseName}.tsv
    """  
}

process comet_txt_to_tsv_with_fdr {
    publishDir "${params.outdir}/tsvs", mode:'copy'

    input:
    tuple file(txt), file(fasta)

    output:
    file "${txt.baseName}*.tsv"
    file "${txt.baseName}*${params.fdr}.tsv"
    file txt

    """
    comet_txt_to_tsv_prune.py ${fasta} ${txt} ${txt.baseName} ${params.fdr} ${params.use_percolator}
    """
}

process count_and_seperate_psms {
    publishDir "${params.outdir}/tsvs_psms_statistics", mode:'copy'
    
    input:
    file(fdr_filtered_tsv)
    file(txt_file)

    output:
    file "${txt_file.baseName}_${params.fdr}fdr_ident*.tsv"

    """
    count_and_create_fdr_statistics.py ${params.fasta_not_from_protgraph} ${fdr_filtered_tsv} ${txt_file.baseName}_${params.fdr}fdr_ident ${params.statistics_remove_varmods} ${params.statistics_unique_proteins}
    """
}
