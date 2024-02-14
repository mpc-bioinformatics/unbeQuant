#!/bin/env python

import sys
import argparse
import os
import csv
import json
import subprocess

import pyopenms


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-featurexml", help="Input-featureXML-File")
    parser.add_argument("-rawfile", help="Input-RAW-file, from which the XICs should be extracted")
    parser.add_argument("-trfp_executable", help="Executable to ThermoRawFileParser")
    parser.add_argument("-out_tsv", help="Output-tsv-file, summarizing features")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # First load features
    fmap = pyopenms.FeatureMap()
    pyopenms.FeatureXMLFile().load(args.featurexml, fmap)
    

    # Create Feature-Queries for TRFP
    queries = []
    features_mapping = []

    for feature in fmap:
        # For each feature

        pep_idents = []
        prot_idents = []
        for pep_id in feature.getPeptideIdentifications():
            seq = pep_id.getHits()[0].getSequence().toString()
            prots = pep_id.getMetaValue("proteins")

            # If a peptide is already in idents, we do not append, otherwise we append.
            # E.G.: PEPa -> ProtA, PEPb -> ProtA, we want [PEPa, PEPb] and [ProtA, ProtA]
            # E.G.: PEPa -> ProtA, PEPa -> ProtA, we want only [PEPa] and [ProtA]
            if seq not in pep_idents:
                pep_idents.append(seq)
                prot_idents.append(prots)

        # Prepare query for TRFP and collect other information from the feature
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

            # TODO Identification is missing
            features_mapping.append(("f_" + str(feature.getUniqueId()), mz_start, mz_end, rt_start, rt_end, feature.getCharge(), pep_idents, prot_idents))


    # Write query file for TRFP
    with open('queries.json', 'w') as json_out:
        json.dump(queries, json_out)


    # Retrieve XICs, which are saved in the order as the queries (we can itereate them, via zip)
    subprocess.call([
        "mono",
        args.trfp_executable,  # Executable
        "xic",  # use XIC
        "-i", args.rawfile, # INPUT-RAW-File
        "-j", "queries.json",  # INPUT-Queries
        "-o",  os.getcwd() # OUTPUT- Folder
    ])

    # Parse output and prepare for saving as tsv
    with open(os.getcwd() + os.sep + ".".join(args.rawfile.split(os.sep)[-1].split(".")[:-1]) + ".json", "r") as in_trfp:
        xics = json.load(in_trfp)

        # Below it is broken down which columns we have
        # We itereate at results and other saved information in parallel
        prev_id = None
        tsv_rows = []
        for fm, xic in zip(features_mapping, xics["Content"]):
            if fm[0] ==  prev_id:
                # We append to the last entry
                tsv_rows[-1][1] += sum(xic["Intensities"])
                tsv_rows[-1][5].append(fm[1]) # Append mz_start
                tsv_rows[-1][6].append(fm[2]) # Append mz_end
                tsv_rows[-1][7].append(fm[3]) # Append rt_start
                tsv_rows[-1][8].append(fm[4]) # Append rt_end
                tsv_rows[-1][9].append(xic["RetentionTimes"]) # Append retentiontime
                tsv_rows[-1][10].append(xic["Intensities"]) # Append intensity

            else:
                prev_id = fm[0]
                # We generate a new entry:
                tsv_rows.append([
                    fm[0], # openms_fid
                    sum(xic["Intensities"]), # sum_intensity
                    fm[5], # charge
                    # These two are in lists in case that there are many (different) idents for a feature. (can be iterated via a zip-loop)
                    fm[6], # pep_idents (if any)
                    fm[7], # prot_idents
                    # These are in lists (and list in lists for each isotope)
                    # These are parallel lists and can be itereated in parallel (e.g. with a zip-loop)
                    [fm[1]], # mz_starts
                    [fm[2]], # mz_ends
                    [fm[3]], # rt_starts in seconds
                    [fm[4]], # rt_ends in seconds
                    [xic["RetentionTimes"]], # retentiontimes
                    [xic["Intensities"]] # intensities
                ])   

    # Write final table
    with open(args.out_tsv, "w") as out_file:
        csv_out = csv.writer(out_file, delimiter="\t")

        # Write header
        csv_out.writerow([
            "openms_fid",
            "sum_intensity",
            "charge",
            "l_pep_ident",
            "l_prot_ident",
            "l_mz_start",
            "l_mz_end",
            "l_rt_start",
            "l_rt_end",
            "l_retention_times",
            "l_intensities"
        ])

        # Write data
        csv_out.writerows(tsv_rows)


