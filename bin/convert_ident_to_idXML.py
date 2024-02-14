#!/bin/env python

import sys
import argparse
import os

import pyopenms
import pandas as pd


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-tsv_file", help="The tsv-Identification file, containing the needed columns")
    parser.add_argument("-use_protgraph", help="Flag to either include fasta_ids (not protgraph) or fasta_descs (protgraph)")
    parser.add_argument("-out_idxml", help="The Output idxml")

    return parser.parse_args()


if __name__ == "__main__":
    args = argparse_setup()

    if args.use_protgraph == "true":
        column_ident = "fasta_desc"
    else:
        column_ident = "fasta_acc"

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
    pyopenms.IdXMLFile().store(args.out_idxml, [protein_id], peptide_ids)
