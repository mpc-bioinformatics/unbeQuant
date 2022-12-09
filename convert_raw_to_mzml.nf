#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.thermo_raws = "$PWD/raws"  // Databasefile in SP-EMBL
params.outdir = "$PWD/results"  // Output-Directory of the FASTA-results

params.additional_params = ""
params.num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)

workflow {
    // Convert the file to mzML
    rawfiles = Channel.fromPath(params.thermo_raws + "/*.raw")
    convert_to_raw_via_thermorawfileparser(rawfiles)
}

process convert_to_raw_via_thermorawfileparser {
    maxForks params.num_procs_conversion
    stageInMode "copy"

    publishDir "${params.outdir}/", mode:'copy'

    input:
    file raw

    output:
    file "${raw.baseName}.mzML"

    """
    run_thermorawfileparser.sh ${params.additional_params} --format=1 --output_file=${raw.baseName}.mzML --input=${raw} 
    rm ${raw}
    """
}
