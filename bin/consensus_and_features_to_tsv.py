#!/bin/env python

import argparse
import os
import csv
from ast import literal_eval

import tqdm
import pandas as pd
import pyopenms


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-featurexmls_tsvs", help="Input-featureXML-TSV-Files, comma-seperated and converted by 'features_to_tsv'")
    parser.add_argument("-consensus", help="Input")
    parser.add_argument("-out_tsv", help="Output-tsv-file, summarizing all the results across all input_files")
    parser.add_argument("-out_tsv_reduced", help="Output-tsv-file, summarizing a smaller result set in the table")
    parser.add_argument("-out_tsv_minimal", help="Output-tsv-file, summarizing a minimal result set in the table")

    return parser.parse_args()


if __name__ == "__main__":
    # Loads the consensus file and all the feature tsv files.
    # This file then generates the final report tsv files (full, reduced and minimal)
    args = parse_args()

    # Loading all feature tsv files
    dict_of_single_features = dict()
    dict_of_single_features_sanity_check = dict()
    for ftsv in args.featurexmls_tsvs.split(","):
        filename = ".".join(ftsv.split(os.sep)[-1].split(".")[:-1])
        print("Loading features of file: {}".format(filename))
        dict_of_single_features[filename] = pd.read_csv(ftsv, sep='\t',
            converters={
                "l_pep_ident": literal_eval,
                "l_prot_ident": literal_eval,
                "l_mz_start": literal_eval,
                "l_mz_end": literal_eval,
                "l_rt_start": literal_eval,
                "l_rt_end": literal_eval,
                "l_retention_times": literal_eval,
                "l_mass_to_charges": literal_eval,
                "l_intensities": literal_eval,
                "l_ms2_scans": literal_eval
                }
        )
        dict_of_single_features_sanity_check[filename] = [False]*len(dict_of_single_features[filename])

    # First Conensus Map (for each consensus we will write a row in the final report)
    print("Loading consensus...")
    cmap = pyopenms.ConsensusMap()
    pyopenms.ConsensusXMLFile().load(args.consensus, cmap)

    # Get number mapping:
    map_list =  cmap.getColumnHeaders()
    for key,val in map_list.items():
        map_list[key] = ".".join(val.filename.split(".")[:-1])

    # Set header order for  full, reduced and minimal (and all) output
    header_per_file = [
        "intensity",
        "l_pep_ident",
        "l_prot_ident",
        "openms_fid",
        "charge",
        "l_ms2_scans",
        "l_mz_start",
        "l_mz_end",
        "l_rt_start",
        "l_rt_end",
        "l_retention_times",
        "l_mass_to_charges",
        "l_intensities"
    ]

    header_per_file_reduced = [
        "intensity",
        "l_pep_ident",
        "l_prot_ident",
        "openms_fid",
        "charge",
        "l_ms2_scans",
    ]

    header_per_minimal = [
        "intensity",
        "openms_fid",
        "charge",
        "l_ms2_scans",
    ]

    global_headers = [
        "first_iso_global_min_mz",
        "first_iso_global_max_mz",
        "first_iso_global_min_rt",
        "first_iso_global_max_rt",
        "feature_global_min_mz",
        "feature_global_max_mz",
        "feature_global_min_rt",
        "feature_global_max_rt"
    ]
    # Set File-Order
    filenames = sorted(dict_of_single_features.keys())

    # Create final results_file
    with open(args.out_tsv, "w") as out_file, open(args.out_tsv_reduced, "w") as out_file_reduced, open(args.out_tsv_minimal, "w") as out_file_minimal:
        out_csv = csv.writer(out_file, delimiter="\t")
        out_csv_reduced = csv.writer(out_file_reduced, delimiter="\t")
        out_csv_minimal = csv.writer(out_file_minimal, delimiter="\t")

        # For each feature, we get all infos from all files and headers in the above mentioned order
        # First get header
        row_header = ["openms_ceid"]
        for h in header_per_file:
            for f in filenames:
                row_header.append(
                    f + "_____" + h
                )
        for gh in global_headers:
            row_header.append(gh)
        out_csv.writerow(row_header)

        row_header = ["openms_ceid"]
        for h in header_per_file_reduced:
            for f in filenames:
                row_header.append(
                    f + "_____" + h
                )
        for gh in global_headers:
            row_header.append(gh)
        out_csv_reduced.writerow(row_header)

        row_header = ["openms_ceid"]
        for h in header_per_minimal:
            for f in filenames:
                row_header.append(
                    f + "_____" + h
                )
        for gh in global_headers:
            row_header.append(gh)
        out_csv_minimal.writerow(row_header)

        print("Writing consensus table")
        # Now get the actual data
        for ce in tqdm.tqdm(cmap, total=len([None for _ in cmap]), unit="consensus_element(s)"):
            entry = ["e_" + str(ce.getUniqueId())]
            entry_reduced = ["e_" + str(ce.getUniqueId())]
            entry_minimal = ["e_" + str(ce.getUniqueId())]
            feature_list = ce.getFeatureList()

            files_involved = [map_list[fl.getMapIndex()] for fl in feature_list]
            feature_in_files_involved = [fl.getUniqueId() for fl in feature_list]

            # Preload rows
            file_rows = []
            for fi, fi_fid in zip(files_involved, feature_in_files_involved):
                filter_row = dict_of_single_features[fi]["openms_fid"] == "f_" + str(fi_fid)
                file_rows.append(dict_of_single_features[fi][filter_row])

                # SanityCheck if Feature was also in ConsensusMap
                dict_of_single_features_sanity_check[fi][file_rows[-1].index[0]] = True

            mzs, mze, rts, rte = [], [], [], []
            for f in filenames:
                if f in files_involved:
                    fidx = files_involved.index(f)
                    if len(file_rows[fidx]["l_mz_start"].values[0]) != 0:
                        mzs.append(float(file_rows[fidx]["l_mz_start"].values[0][0]))
                    if len(file_rows[fidx]["l_mz_end"].values[0]) != 0:
                        mze.append(float(file_rows[fidx]["l_mz_end"].values[0][0]))
                    if len(file_rows[fidx]["l_rt_start"].values[0]) != 0:
                        rts.append(float(file_rows[fidx]["l_rt_start"].values[0][0]))
                    if len(file_rows[fidx]["l_rt_end"].values[0]) != 0:
                        rte.append(float(file_rows[fidx]["l_rt_end"].values[0][0]))

            g_mzs, g_mze, g_rts, g_rte = [], [], [], []
            for f in filenames:
                if f in files_involved:
                    fidx = files_involved.index(f)
                    if len(file_rows[fidx]["l_mz_start"].values[0]) != 0:
                        g_mzs.extend([float(x) for x in file_rows[fidx]["l_mz_start"].values[0]])
                    if len(file_rows[fidx]["l_mz_end"].values[0]) != 0:
                        g_mze.extend([float(x) for x in file_rows[fidx]["l_mz_end"].values[0]])
                    if len(file_rows[fidx]["l_rt_start"].values[0]) != 0:
                        g_rts.extend([float(x) for x in file_rows[fidx]["l_rt_start"].values[0]])
                    if len(file_rows[fidx]["l_rt_end"].values[0]) != 0:
                        g_rte.extend([float(x) for x in file_rows[fidx]["l_rt_end"].values[0]])
            global_info = [
                min(mzs), max(mze), min(rts), max(rte),
                min(g_mzs), max(g_mze), min(g_rts), max(g_rte)
            ]

            # Now Add entry information and write to the tsv files
            for h in header_per_file:
                for f in filenames:
                    if f in files_involved:  # If data, append it
                        fidx = files_involved.index(f)
                        entry.append(
                            file_rows[fidx][h].values[0]
                        )
                    else:  # If no data, append None
                        entry.append(
                            None
                        )
            out_csv.writerow(entry + global_info)

            for h in header_per_file_reduced:
                for f in filenames:
                    if f in files_involved:  # If data, append it
                        fidx = files_involved.index(f)
                        entry_reduced.append(
                            file_rows[fidx][h].values[0]
                        )
                    else:  # If no data, append None
                        entry_reduced.append(
                            None
                        )
            out_csv_reduced.writerow(entry_reduced + global_info)

            for h in header_per_minimal:
                for f in filenames:
                    if f in files_involved:  # If data, append it
                        fidx = files_involved.index(f)
                        entry_minimal.append(
                            file_rows[fidx][h].values[0]
                        )
                    else:  # If no data, append None
                        entry_minimal.append(
                            None
                        )
            out_csv_minimal.writerow(entry_minimal + global_info)

    # Sanity-Check, check if all features were used by generating the final consensus table
    for key, val in dict_of_single_features_sanity_check.items():
        if len(val) != sum(val):
            print("WARNING: Some features were not described by consensus for {}".format(key))
            indices = [f_idx for f_idx, f in enumerate(dict_of_single_features_sanity_check[key]) if f == True]
            print("Not included feature: {}".format(
                str(list(dict_of_single_features[key]["openms_fid"].iloc[indices]))
            ))
