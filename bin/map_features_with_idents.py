#!/bin/env python

import argparse
import os

import numpy as np
import pyopenms
import pandas as pd
import plotly.express as ex


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-tsv_file", help="The tsv-Identification file, containing the needed columns")
    parser.add_argument("-featureXML", help="The featureXML file, containing the user param about the MS2 scans")
    parser.add_argument("-cutoff", help="Cutoff value used to remove low abundant features. CAUTION: This accounts for the previous set intensity calculation method which happens acroos all traces for a feature., Use 'tX' to remove all features with intensity <= X. Alternatively you can use 'qX' to remove the lowest X percent of intensities. E.G.: 't1000' or 'q0.05'.", default="t0")
    parser.add_argument("-use_protgraph", help="Flag to either include fasta_ids (not protgraph) or fasta_descs (protgraph)")
    parser.add_argument("-out_featurexml", help="The Output feature XML file with annotated identifications")
    parser.add_argument("-out_plot_cutoff", help="The Output directory for the plotsfeature XML file with annotated identifications", default="feature_cutoff_plot.html")

    return parser.parse_args()

if __name__ == "__main__":
    # Map feature with identifications (and remove features, corresponding to the cutoff)
    # the features are mapped with their corresponding identification, using the already added MS2 scan information in previous a step
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
    raw_peptide_col = "peptide_seq" if "peptide_seq" in psms else "modified_peptide"
    peptide_ids = []
    for _, psm in psms.iterrows():
        peptide_id = pyopenms.PeptideIdentification()
        peptide_id.setRT(psm['retention_time'])
        peptide_id.setMZ(psm['exp_mass_to_charge'])
        peptide_id.setScoreType('q-value')
        peptide_id.setHigherScoreBetter(False) # This was the suggested change from the comment
        peptide_id.setIdentifier(run_name)
        peptide_id.setMetaValue("proteins", psm[column_ident])
        peptide_hit = pyopenms.PeptideHit()
        peptide_hit.setScore(psm['qvalue'])
        peptide_hit.setRank(1)
        peptide_hit.setCharge(psm['charge'])
        peptide_hit.setSequence(pyopenms.AASequence.fromString(psm['plain_peptide']))
        peptide_hit.setMetaValue("raw_peptide", psm[raw_peptide_col])
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

    # Get intensity cutoff
    feature_intensities = pd.Series([x.getIntensity() for x in features])
    if args.cutoff.startswith("t"):
        cutoff = float(args.cutoff[1:])
        q_cutoff = sum(feature_intensities <= cutoff)/len(feature_intensities)
    elif args.cutoff.startswith("q"):
        q_cutoff = float(args.cutoff[1:])
        if q_cutoff > 1 or q_cutoff < 0:
            raise ValueError("Quantile cutoff must be between 0 and 1")
        cutoff = feature_intensities.quantile(q_cutoff)
    else:
        raise ValueError("Cutoff value not correctly formatted. Please use 't' for intensity cutoff or 'q' for quantile cutoff.")

    # Plot the intensities and the cutoff
    q_pts = np.linspace(0,1, num=1000)
    i_pts = feature_intensities.quantile(q=q_pts)
    fig = ex.scatter(x=q_pts, y=i_pts, log_y=False, title=f"Intensities of Features broken into q-quantiles ({args.featureXML})", labels={"x":"Quantile", "y":"Intensity"})
    fig.add_vline(x=q_cutoff, line_color="red", line_width=1, line_dash="dash", annotation_text=f"q_cutoff = {q_cutoff}", annotation_position="bottom right")
    fig.add_hline(y=cutoff, line_color="red", line_width=1, line_dash="dash", annotation_text=f"i_cutoff = {cutoff}", annotation_position="top left", annotation_y=np.log10(cutoff))
    fig.add_shape(type="rect", x0=0, y0=min(feature_intensities), x1=q_cutoff, y1=cutoff, line=dict( color="Gray", width=2,), fillcolor="Gray", opacity=0.25)
    fig.update_yaxes(type="log")
    fig.write_html(args.out_plot_cutoff)  

    # Only keep features with intensity above the cutoff
    ident_idcs = []
    for f in features:
        # Append feature to the output, only if min_intensity is met
        if f.getIntensity() >= float(cutoff):
            # Get Scan Idx and check if it has an identification
            scans = f.getMetaValue("unbeQuant_MS2_Scan_Map")        
            idcs = psms.index[psms["scan"].isin(scans)].tolist()

            # Add the Identification to the features and keep track of add rows
            f.setPeptideIdentifications([peptide_ids[idx] for idx in idcs])
            f.setMetaValue("unbeQuant_MS2_Scan_Map", scans)
            ident_idcs.extend(idcs)

            # We only add the feature if the intensity is above the threshold
            ident_features.push_back(f)

    # Add unassigned identifications
    unassigned_idcs = set(range(len(psms))).difference(set(ident_idcs))
    ident_features.setUnassignedPeptideIdentifications([peptide_ids[idx] for idx in unassigned_idcs])

    # Save final feature XML file with identifications
    ident_fh.store(args.out_featurexml, ident_features)
