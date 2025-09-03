#!/bin/env python

import argparse
import sys
import csv
csv.field_size_limit(sys.maxsize)

import pyopenms


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-featurexml", help="Input-featureXML-File")
    parser.add_argument("-out_csv", help="Output-csv-file which can be used with the xic-extractor")

    return parser.parse_args()

if __name__ == "__main__":
    # Preparing the xic-extractor query file from a featureXML file
    args = parse_args()

    # First load features
    fmap = pyopenms.FeatureMap()
    pyopenms.FeatureXMLFile().load(args.featurexml, fmap)
    
    # Create Feature-Queries for xic extraction
    queries = []
    features_mapping = []

    # For each feature
    for feature in fmap:
        scans = feature.getMetaValue("unbeQuant_MS2_Scan_Map")

        pep_idents = []
        raw_pep_idents = []
        prot_idents = []
        for pep_id in feature.getPeptideIdentifications():
            seq = pep_id.getHits()[0].getSequence().toString()
            raw_pep = pep_id.getHits()[0].getMetaValue("raw_peptide")
            prots = pep_id.getMetaValue("proteins")
            

            # If a peptide is already in idents, we do not append, otherwise we append.
            # E.G.: PEPa -> ProtA, PEPb -> ProtA, we want [PEPa, PEPb] and [ProtA, ProtA]
            # E.G.: PEPa -> ProtA, PEPa -> ProtA, we want only [PEPa] and [ProtA]
            if seq not in pep_idents:
                pep_idents.append(seq)
                raw_pep_idents.append(raw_pep)
                prot_idents.append(prots)

        # Prepare query for XIC-Extractor and collect other information from the feature (isotope wise)
        for hull in feature.getConvexHulls():
            # Query for each isotope
            bb = hull.getBoundingBox()
            rt_start, mz_start = bb.minPosition()
            rt_end, mz_end = bb.maxPosition()

            d = dict()
            d["mz_start"] = mz_start
            d["mz_end"] = mz_end
            d["rt_start"] = rt_start / 60 # in minutes
            d["rt_end"] = rt_end / 60 # in minutes
            queries.append(d)

            # Append the feature mapping with all the information for a single query (for a single isotope)
            features_mapping.append(("f_" + str(feature.getUniqueId()), mz_start, mz_end, rt_start, rt_end, feature.getCharge(), pep_idents, raw_pep_idents, prot_idents, scans))


    # Write query file for xic-extractor
    with open(args.out_csv, "w") as out_csv_file:
        csv_out = csv.writer(out_csv_file)

        csv_out.writerow([
            "identifier", "mz_start", "mz_end", "rt_start", "rt_end", "mz", "ppm", "ms_level", ## Needed for the xic_extractor
            "charge", "pep_idents", "raw_pep_idents", "prot_idents", "ms2_scans"  ## Additional information to be caried over to the next step
        ])

        for fm, q in zip(features_mapping, queries):
            csv_out.writerow([
                fm[0], q["mz_start"], q["mz_end"], q["rt_start"], q["rt_end"], None, None, "ms",
                fm[5], fm[6], fm[7], fm[8], fm[9]
            ])
