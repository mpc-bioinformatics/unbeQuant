![unbeQuant_logo](resources/logo_full_25.png)

---

Qauntifying the unknown on MS1 level. UnbeQuant allows the quantification of measured ions without identification annotations in a DDA setting for mass spectrometry proteomics data from Bruker or Thermo mass spectrometers. To achieve this it uses identification results and sets same identifications across runs as anchors to align multiple runs, providing a mixture of the following: Identified ions with quantitative values, only some identified ions with quantitative values as well as ions without any identification information.

The realized workflows are implemented with the [Nextflow DSL 2](https://www.nextflow.io/docs/latest/module.html) language utilizing the importing capabilties and reusing parts of existing workflows from [ProtGraph/ProGFASTAGen](https://github.com/mpc-bioinformatics/ProGFASTAGen) and [XIC-Extractor](https://github.com/mpc-bioinformatics/xic-extractor). The provided workflow can be imported and further used. These workflows utilize [biosaur2](https://github.com/markmipt/biosaur2) for feature detection, [OpenMS](https://openms.de/) for retention time alignment and consensus generation and many python scripts for intermediate steps and data conversion.

If you are interested in **unbeQuant**, checkout **materials_and_posters** to learn more about it or see **Prerequisites** of how to setup unbeQuant locally. This repository provides two workflow variants containing the quantification. A brief description on how these can be called individually can be found in **Workflow Scripts**.


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

## Results Structure

UnbeQuant organizes all its results into a single output folder. You can define the location of this folder using `--main_outdir`. The results are structured by the following:

* Conversion and Identification results from [ProtGraph/ProGFASTAGen](https://github.com/mpc-bioinformatics/ProGFASTAGen)
* XIC-Extraction results from [XIC-Extractor](https://github.com/mpc-bioinformatics/xic-extractor) as well as
* UnbeQuant results in the folder `quantification` (This contains un-/identified features and linked features between multiple runs)


Within the `quantification` folder, you will find the followins sub folders:

| Folder                                  | Contents            | Description                                                                                                            |
|-----------------------------------------|---------------------|------------------------------------------------------------------------------------------------------------------------|
| extracted_xics                          | Extracted XICs      | Contains the Extracted Ion Chromatograms generated from the features found by the feature finder (*.tsv).              |
| features_with_annotated_identifications | Mapped Feature Data | Contains various files after the mapping of features with identification results:                                      |
|                                         |                     | - *.featurexml (with identifications) and *.trafoXML (for alignment).                                                  |
|                                         |                     | - *.tsv files containing the features, their boundaries, associated identifications (if any), and the extracted XICs.  |
| visualizations                          | HTML-Plots          | Contains *.htmls-plots which were generated during execution of UnbeQuant                                              |
| final_report                            | Final Result Tables | Contains the final quantitative results of un-/identified features across multiple runs in tabular format (*.csv and *.xlsx)             |
|                                         |                     | - NOTE: Reduced versions of these tables (omitting specific columns) are also available.                               |
|                                         |                     | - CAUTION: *.xlsx files may truncate cells or fail to display all rows. Use the *.csv files if possible                |

### Visualization

UnbeQuant generates, as of now, during and at the end following plots in the `visualizations`-folder: 

* Retention-Time-Alignment Plots, which show how the data was transformed to overlap with each other. X-axis describes the original retention time, while the y-axis desribes the transformed retention time.
* Cutoff Plots, showcasing how many features would be filtered per file The plot shows in an ascending order, the quantitative value (y-axis) for all features per file. The x-axis shows the percentage of how many features are remaining (quantile, similar as to a boxplot). If setting the `--qal_cutoff` parameter you can set a threshold for each axis globally.
* MINLH-parameter estimated plots, which may be used to tune the `--qal_minlh` for the feature finder [biosaur2](https://github.com/markmipt/biosaur2). These plots should be used as an aid, as they estimate the ratios of traces with and without MS2-spectra or identification annotation while also showing the total number of ions and how many potentially get lost with a higher minlh parameter. Additionally a trade off point may be present, where as many MS2-spectra (or identified traces) are present as traces without MS2-spectra (or identifications). The user can choose to select a higher minlh parameter as the trade off resulting in a higher ratio of traces with MS2-spectra (or identifications) or vice versa, depending on the use-case.

## Workflow Exectuion

### Exeution of single script

UnbeQuant only adds a single additonal nextflow script `quantify_and_align.nf`, which could be used as a standalone script. In this step, it is expected to have already converted data as well as identified data. You can call it with the following:

```txt
nextflow run quantify_and_align.nf \
    --qal_spectra_files < Folder containing RAW-files > \
    --qal_mzmls < Folder containing mzmls-files, must be named as in RAW-files > \
    --qal_idents < Folder with identifications files in tsv format, which need  to be named as the raw_files with a specfic suffix: "*qvalue_no_decoys_fdr_0.0[15].tsv" > \
    --qal_outdir < Output-Folder, where the quantification results will be saved >
```

### Whole Workflow execution

It is easier to run UnbeQuant with the `main_workflow_*` scripts, which are already available in this repository and wrap around `quantify_and_align.nf`.  We provide the following scripts:

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
