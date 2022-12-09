#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.mzml_folder = "$PWD/MZMLs"  // Input-Directory of mzMLs, which should be used for identification
params.fasta_file = "peptides.fasta"  // Database (FASTA-file) used for identification WITHOUT DECOYS!

// Results output dir
params.outdir = "$PWD/results"  // Output-Directory of the FASTA-results

// Optional, but should be set!
params.search_parameter_file = "$PWD/example_configurations/msgfplus_config.txt"  //Search Parameters for MSGFPlus
params.tda = 1 // 0 --> No Target-Decoy appraoch | 1 --> Target-Decoy appraoch
params.fdr = 0.05 // If TDA is set to 1, also return the pruned tsvs, whit up to fdr. (originals are also kept)
params.num_parallel_searches = Runtime.runtime.availableProcessors()
params.decoy_prefix = "XXX"

// Parameters for the JVM
params.jvm_params = "Xmx3500M"

workflow {
    // Get all mzML files which should be identified
    mzmls = Channel.fromPath(params.mzml_folder + "/*.mzML")

    // Get FASTA-file
    fasta_file = Channel.fromPath(params.fasta_file)

    // Get Modification Parameters file
    modifications_file = Channel.fromPath(params.search_parameter_file)

    // Build indexed fasta for MSGFPLUS
    msgfplus_buildsa(fasta_file)

    // Combined channel replicated the indexed fasta for each MZML to be reused
    combined_channel = fasta_file
        .combine(modifications_file)
        .combine(mzmls)
        .combine(msgfplus_buildsa.out.toList())
        
    // Start search
    msgfplus_search_mzml(combined_channel)

    // Convert to user readable tsv format
    mzid_to_tsv(msgfplus_search_mzml.out)

    // if tda, then also return the pruned results
    if (params.tda == 1) {
        tsv_pruned_via_fdr(mzid_to_tsv.out)
    }
}

process msgfplus_buildsa {

    input:
    file input_fasta

    output:
    file "${input_fasta.baseName}*"
    // file "${input_fasta.baseName}*.fasta"
    // file "${input_fasta.baseName}*.canno"
    // file "${input_fasta.baseName}*.cnlcp"
    // file "${input_fasta.baseName}*.csarr"
    // file "${input_fasta.baseName}*.cseq"    

    """
    buildsa_msgfplus.sh ${params.jvm_params} -d ${input_fasta} -tda ${params.tda} -decoy ${params.decoy_prefix}
    """
}

process msgfplus_search_mzml {
    maxForks params.num_parallel_searches
    stageInMode "copy"
    publishDir "${params.outdir}/mzid", mode:'copy'

    input:
    tuple file(input_fasta), file(mod_file), file(mzml_file), file(remainder)

    output: 
    file "${mzml_file.baseName}.mzid"


    """
    touch -m ${input_fasta.baseName}*

    run_msgfplus.sh ${params.jvm_params} -d ${input_fasta} -s ${mzml_file} -tda ${params.tda} -decoy ${params.decoy_prefix} -conf ${mod_file} -o ${mzml_file.baseName}.mzid

    # Since we copy the input files, we later remove them manually to save memory
    rm ${input_fasta.baseName}*
    rm ${mzml_file}
    rm ${mod_file}
    """  
}

process mzid_to_tsv {
    publishDir "${params.outdir}/tsvs", mode:'copy'

    input:
    file mzid

    output:
    file "${mzid.baseName}.tsv"

    """
    convert_msgfplus.sh -mzid ${mzid} -geneid \\'.*\\' -tsv ${mzid.baseName}.tsv
    """
}


process tsv_pruned_via_fdr {
    publishDir "${params.outdir}/tsvs_pruned", mode:'copy'

    input:
    file tsv

    output:
    file "${tsv.baseName}_pruned.tsv"

    """
    tsv_qvalue_prune.py ${tsv} ${tsv.baseName}_pruned.tsv ${params.fdr}
    """
}
