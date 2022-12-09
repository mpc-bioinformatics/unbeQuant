# TODO

The LFQ is currently not implemented in Nextflow

Current Roadmap: 

0. ~~ThermoRawFileParser (generate profile mzML of MS1)~~
1. OpenMS Feature Detection (FeatureFinderIsotopeWaveloet)
2. Alignstein (https://github.com/grzsko/Alignstein)
3. XIC Extraction (for AUR, ThermoRawFileParser via JSON)
4. Mapping Identifications with Features (custom Python script)
5. Generate readable and usable XLSX and/or tsv (which can be later used for statistical analysis)

For 2.

```text
/home/luxdo/Desktop/spass/openms/usr/bin/FeatureFinderIsotopeWavelet -in Example_Folder/mzmls/QEXHF19830std.mzML -out QEXHF19830std.mzML.featureXML -algorithm:hr_data
```

# Readme

In this project we summarized the workflows as well as the software for XXX.

## FASTA-Generation

We define the two workflow to generate FASTA-files using the precursor masses from mzML (MS2).

Before executing those, make sure to execute `setup_befeore_nextflow.sh` prior to compile the cpp graph traversal and to install ProtGraph from PyPI.

Execute the nextflow workflow in the `pipenv` environment (please use `pipenv shell` before running `nextflow run ...`)

## Identification Requirements MSGFPlus

Make sure that for the identification with MSGFPlus, the `java` executable is available.

The exectuable `mono` is also required to convert the mzid into tsv-tables (human readable formats)

## Converting RAW file into mzMLs

Please use `nextflow run convert_raw_to_mzml.nf -with-docker alpine`
