#!/bin/env python

import sys
import csv
csv.field_size_limit(sys.maxsize)
import re 

use_fasta_id = sys.argv[1].lower() == "true" # Flag to use Fasta id for protein-fasta or fasta_desc for protgraph generated peptide-fasta
input_statistics = sys.argv[2] # Result-file from Comet with already added information from the fasta (containing the fasta_id and fasta_desc column)
write_base_name = sys.argv[3] # Base_name how the files should be written

remove_varmods = sys.argv[4].lower() == "true" # Remove Variable Modifications (Count a peptide with and without variable modification as unique)
unique_proteins = sys.argv[5].lower() == "true" # Count PSMs, which occur multiple times in the same proteins as unique 

unique_psms = []
shared_psms = []

# Only calculated when fasta_desc is used.
unique_feature_psms = []
feature_psms = []

strings_considered_features = ["CONFLICT", "SIGNAL", "INIT_MET", "PROPEP", "PEPTIDE", "MUTAGEN", "VARIANT"]

# Regex to retrieve single proteins from header
regex_pg_header = r"[A-Z0-9-_]+?\(.*?\)"
regex_get_varmods = r"VARMOD\[.*?\](?:,|)"

def fits_in_list_fasta_id(line):
    matches = line.split(",")
    if unique_proteins:
        # Count unique proteins instead (in case of same sequence in 1 proteins)
        matches = set(matches)

    is_unique = len(matches) == 1
    is_shared = len(matches) > 1
    return is_unique, is_shared, False, False

def fits_in_list_fasta_desc(line):
    # Parse every match
    line = line.replace(" ", "")
    matches_regex = re.finditer(regex_pg_header, line, re.MULTILINE)
    matches = [x.group() for x in matches_regex]
    if remove_varmods:  
        # Special Case, since FASTA-files CANNOT encode modifications. 
        # We consider the VARMODs (since FIXMODS are globally present in the description)
        # from the same protein as unique. 
        # E.G.: >pg|ID_XXXX|P68871(42:60,mssclvg:0,),P68871(42:60,mssclvg:0,VARMOD[56:56,M:15.994915])
        # would be considered unique, since only a variable modification was applied.
        # However, for other modification it could be interesting to consider them as unique (ubiquitination, phosporylation)
        matches = [re.sub(regex_get_varmods, "", m).replace(",", "") for m in matches]
        matches = set(matches)

    if unique_proteins:
        # Count unique proteins instead (in case of same sequence in 1 proteins)
        matches = set([m[:m.find("(")] for m in matches])

    # Get shared and unique information
    is_unique = len(matches) == 1
    is_shared = len(matches) > 1

    # Get shared unique information about features
    feature_check = [any([y in x for y in strings_considered_features]) for x in matches]
    is_feature_unique = (len(matches) == 1) and feature_check[0]
    is_feature_shared = len(matches) > 1 and all(feature_check)

    return is_unique, is_shared, is_feature_unique, is_feature_shared


if __name__ == "__main__":
    with open(input_statistics, "r") as input_csv_file:
        csv_in = csv.reader(input_csv_file, delimiter = "\t")

        # Read header
        header = next(csv_in)

        # Set index for reading and set parsing method
        fasta_id_desc_index = header.index("fasta_id") if use_fasta_id else header.index("fasta_desc")
        classify_psm = fits_in_list_fasta_id if use_fasta_id else fits_in_list_fasta_desc

        for l in csv_in:
            unique, shared, feature_unique, feature_shared = classify_psm(l[fasta_id_desc_index])

            if unique:
                unique_psms.append(l)
            if shared:
                shared_psms.append(l)
            if feature_unique:
                unique_feature_psms.append(l)
            if feature_shared:
                feature_psms.append(l)                                                

    with open(write_base_name + "_summary.tsv", "w") as summary_out:
        csv_out = csv.writer(summary_out, delimiter="\t")
        csv_out.writerow(["count_idents", "count_unique", "count_shared", "count_unique_with_features", "count_shared_with_only_features"])
        csv_out.writerow([len(unique_psms) + len(shared_psms), len(unique_psms), len(shared_psms), len(unique_feature_psms), len(feature_psms)])


    for psms_list, postfix in zip(
        [unique_psms, shared_psms, unique_feature_psms, feature_psms], 
        ["count_unique", "count_shared", "count_unique_with_features", "count_shared_with_only_features"]
    ):
        if len(psms_list) != 0:
            with open(write_base_name + "_" + postfix + ".tsv", "w") as csv_out_file:
                csv_out = csv.writer(csv_out_file, delimiter="\t")
                csv_out.writerow(header)
                csv_out.writerows(psms_list)
