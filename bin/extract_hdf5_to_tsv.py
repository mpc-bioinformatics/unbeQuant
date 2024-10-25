#!/bin/env python

import sys
import argparse
import os
import csv
import json
import subprocess

import h5py

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-hdf5_xic_file", help="Input-hdf5 file")
    parser.add_argument("-xic_query_file", help="Input query file, which has been used to generate the xic")
    parser.add_argument("-out_tsv", help="Output-tsv-file, summarizing features")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Read feature mapping (query file)
    feature_mapping = []
    with open(args.xic_query_file, "r") as in_query:
        in_csv = csv.reader(in_query)
        feature_mapping_headers = next(in_csv)

        for l in in_csv:
            feature_mapping.append(l)

    # Open hdf5 file
    h5 = h5py.File(args.hdf5_xic_file, "r")

    fm_identifier = feature_mapping_headers.index("identifier")
    fm_charge = feature_mapping_headers.index("charge")
    fm_pep_idents = feature_mapping_headers.index("pep_idents")
    fm_prot_idents = feature_mapping_headers.index("prot_idents")
    fm_mz_start = feature_mapping_headers.index("mz_start")
    fm_mz_end = feature_mapping_headers.index("mz_end")
    fm_rt_start = feature_mapping_headers.index("rt_start")
    fm_rt_end = feature_mapping_headers.index("rt_end")


    prev_id = None
    tsv_rows = []
    for idx, fm in enumerate(feature_mapping):
        # XIC-Extractor returns the xics in same order as in the query file, we can simply iterate through both at once
        # We itereate at results and other saved information in parallel
        if fm[fm_identifier] ==  prev_id:
            # We append to the last entry
            tsv_rows[-1][1] += h5["intensities"][idx].sum()
            tsv_rows[-1][5].append(fm[fm_mz_start]) # Append mz_start
            tsv_rows[-1][6].append(fm[fm_mz_end]) # Append mz_end
            tsv_rows[-1][7].append(float(fm[fm_rt_start]) * 60) # Append rt_start
            tsv_rows[-1][8].append(float(fm[fm_rt_end]) * 60) # Append rt_end
            tsv_rows[-1][9].append(list(h5["intensities"][idx])) # Append retentiontime
            tsv_rows[-1][10].append(list(h5["retention_times"][idx])) # Append intensity

        else:
            prev_id = fm[fm_identifier]
            # We generate a new entry:
            tsv_rows.append([
                fm[fm_identifier], # openms_fid
                h5["intensities"][idx].sum(), # sum_intensity
                fm[fm_charge], # charge
                # These two are in lists in case that there are many (different) idents for a feature. (can be iterated via a zip-loop)
                fm[fm_pep_idents], # pep_idents (if any)
                fm[fm_prot_idents], # prot_idents (if any)
                # These are in lists (and list in lists for each isotope)
                # These are parallel lists and can be itereated in parallel (e.g. with a zip-loop)
                [fm[fm_mz_start]], # mz_starts
                [fm[fm_mz_end]], # mz_ends
                [float(fm[fm_rt_start]) * 60], # rt_starts in seconds
                [float(fm[fm_rt_end]) * 60], # rt_ends in seconds
                [list(h5["retention_times"][idx])], # retentiontimes
                [list(h5["intensities"][idx])] # intensities
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
