# UnbeQuant

Qauntifying the unknown on MS1 level. UnbeQuant allows the quantification of measured ions without identification annotations in a DDA setting for mass spectrometry proteomics data from Bruker or Thermo mass spectrometers. To achieve this it uses identification results and sets same identifications across runs as anchors to align multiple runs, providing a mixture of the following: Identified ions with quantitative values, only some identified ions with quantitative values as well as ions without any identification information.

The realized workflows are implemented with the [Nextflow DSL 2](https://www.nextflow.io/docs/latest/module.html) language utilizing the importing capabilties and reusing parts of existing workflows from [ProtGraph/ProGFASTAGen](https://github.com/mpc-bioinformatics/ProGFASTAGen) and [XIC-Extractor](https://github.com/mpc-bioinformatics/xic-extractor). The provided workflow can be imported and further used. These workflows utilize [biosaur2](https://github.com/markmipt/biosaur2) for feature detection, [OpenMS](https://openms.de/) for retention time alignment and consensus generation and many python scripts for intermediate steps and data conversion.

If you are interested in running **unbeQuant**, checkout **Prerequisites** of how to setup unbeQuant. This repository provides two workflow variants and the quantification, containing the actual modules. A brief description on how thes can be called individually can be found in **Workflow Scripts**.

![unbeQuant_logo](resources/logo_full.png)

## Prerequisites

## Cloning this repository

This project uses git submodules to include other git projects into this one. To clone all project at once use the following:

> git clone --recursive https://github.com/mpc-bioinformatics/unbeQuant


Alternatively, use the following commands:

```shell
git clone https://github.com/mpc-bioinformatics/unbeQuant
cd unbeQuant
git submodule init
git submodule update
```

## Executing unbeQuant in Docker (locally)

The workflow is containerized completely via docker. Please follow the [installation guide](https://docs.docker.com/engine/install/ubuntu/) for docker to have it installed in your system, if not already. The workflow can be run as is. Examples of some calls can be found further below.

**NOTE**: If the docker image for unbeQuant is not available anymore, a local docker-container can be build with all needed dependencies for the workflow. if desired. We provide a Dockerfile in the `docker`-folder. To build it, simply execute the following while inside the root folder:

> docker build -t luxii/unbequant:latest . -f docker/Dockerfile

## Executing on Linux (locally without docker)

If desiered, all dependencies can be installed on a linux computer (tested on Ubuntu 22.04 and ArchLinux) and the workflow can be run locally. the following packages need to be installed via apt:

```text
build-essential
wget
curl
unzip
cmake
python3-pip
mono-complete
python-is-python3
libqt5network5
libqt5sql5
git
```

Remaining dependencies can be downloaded via `compile_and_setup_dependencies.sh` for a local execution. Make this file excutable (`chmod +x`) and run it once, but make sure you are in a Python-Environment (e.g. pyenv or conda). This will set up the python dependencies and download additional software. In essence, the `bin` folder is set up.

If the scripts exits without errors, only the `nextflow.config` setting needs to be adapted (deactivating docker) and the provided workflows can be executed without docker.


## Workflow Scripts

TBD

Further below are some example calls 

If you are interested for details of each possible parameter which can be set, check out the scripts in all workflows directly, as all parameters contain a brief description.

# Example Calls

``` shell
# Run UnbeQuant
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
