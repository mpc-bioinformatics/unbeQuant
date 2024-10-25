# UnbeQuant

A worklfow which allows quantification of identified and unidentified feature points across multiple runs. We use a basic identification workflow, followed by a feature detection algorithm and use the results to annotated quantitative data points with identifiacation. These points are used as anchor points to align the retention times across multiple runs and summarizes the results in a large `tsv`-table.

## Executing this workflow

The workflow is containerized via docker. Please follow the [installation guide](https://docs.docker.com/engine/install/ubuntu/) for docker to have it installed in your system, if not already. The workflow can be run as is. Examples of some calls can be found further below.

## Executing this workflow locally (Development)

All dependencies can be downloaded via `compile_and_setup_dependencies.sh` for a local execution. Make sure you are in a Python-Environment,
to prevent errors or dependency conflicts with the operating systems package manager.

The docker can be created via: 

> docker build -t unbeqonet:local . -f docker/Dockerfile


# Example Call

``` shell
# PXD028605 Protein FASTA
nextflow run main_workflow_protein_fasta.nf \
    --main_fasta_file < Path to the FASTA files >.fasta \
    --main_raw_files_folder < Path to the .raw/.d of Thermo/Bruker files which should be analyzed > \
    --main_comet_params < Path to the comet configuration parameters text >

# Run unbeQuant together with ProGFASTAGen (ProtGraph)
nextflow run main_workflow_precursor_specific_fasta.nf \
    --main_sp_embl_file < Path to the SP-EMBL file >.txt \
    --main_raw_files_folder < Path to the .raw/.d of Thermo/Bruker files which should be analyzed > \
    --main_comet_params < Path to the comet configuration parameters text, here with digestion turned off >
```

For a detailed list of possible parameters, check the scripts directly.
