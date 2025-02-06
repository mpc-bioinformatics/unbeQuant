#!/bin/env python

import sys
import argparse
import os

import pyopenms
import pandas as pd


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-tsv_file", help="The tsv-Identification file, containing the needed columns")
    parser.add_argument("-featureXML", help="The featureXML file, containing the user param about the MS2 scans")
    parser.add_argument("-use_protgraph", help="Flag to either include fasta_ids (not protgraph) or fasta_descs (protgraph)")
    parser.add_argument("-out_featurexml", help="The Output feature XML file with annotated identifications")

    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    if args.use_protgraph == "true":
        column_ident = "fasta_desc"
    else:
        column_ident = "fasta_acc"


    # Load features without idents
    features = pyopenms.FeatureMap()
    fh = pyopenms.FeatureXMLFile()
    fh.load(args.featureXML, features)


    # Load the identifications
    run_name = args.tsv_file.split(os.sep)[-1]
    # Taken and modified from https://github.com/OpenMS/OpenMS/issues/4417
    psms = pd.read_csv(args.tsv_file, sep='\t', header=0)
    peptide_ids = []
    for _, psm in psms.iterrows():
        peptide_id = pyopenms.PeptideIdentification()
        peptide_id.setRT(psm['retention_time'])
        peptide_id.setMZ(psm['exp_mass_to_charge'])
        peptide_id.setScoreType('q-value')
        peptide_id.setHigherScoreBetter(False) # This was the suggested change from the comment
        # peptide_id.setIdentifier(psm['spectra_ref'])
        peptide_id.setIdentifier(run_name)
        peptide_id.setMetaValue("proteins", psm[column_ident])
        peptide_hit = pyopenms.PeptideHit()
        peptide_hit.setScore(psm['qvalue'])
        peptide_hit.setRank(1)
        peptide_hit.setCharge(psm['charge'])
        peptide_hit.setSequence(pyopenms.AASequence.fromString(psm['plain_peptide']))
        peptide_id.setHits([peptide_hit])
        peptide_ids.append(peptide_id)


    protein_id = pyopenms.ProteinIdentification()
    protein_id.setIdentifier(run_name)


    # Write the featureXML with identifications
    ident_features = pyopenms.FeatureMap()

    ident_features.setUniqueId(features.getUniqueId())  # Needed for the ConsensusMap generation via OpenMs
    ident_features.setMetaValue("spectra_data", features.getMetaValue("spectra_data"))  # Needed for the ConsensusMap generation
    ident_fh = pyopenms.FeatureXMLFile()
    ident_features.setProteinIdentifications([protein_id])

    ident_idcs = []
    for f in features:
        # Get Scan Idx and check if it has an identification
        scans = f.getMetaValue("unbeQuant_MS2_Scan_Map")        
        idcs = psms.index[psms["scan"].isin(scans)].tolist()

        # Add the Identification to the features and keep track of add rows
        f.setPeptideIdentifications([peptide_ids[idx] for idx in idcs])
        f.setMetaValue("unbeQuant_MS2_Scan_Map", scans)
        ident_idcs.extend(idcs)

        # Append feature to the output
        ident_features.push_back(f)

    # Add unassigned identifications
    unassigned_idcs = set(range(len(psms))).difference(set(ident_idcs))
    ident_features.setUnassignedPeptideIdentifications([peptide_ids[idx] for idx in unassigned_idcs])
    ident_fh.store(args.out_featurexml, ident_features)

