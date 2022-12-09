#!/bin/env python

import sys
import csv
csv.field_size_limit(sys.maxsize)
from collections import defaultdict


INPUT_POUT  = sys.argv[1]
OUTPUT_TSV  = sys.argv[2]

if __name__ == "__main__":

    # Retrieve all psms from percolator and save as dict
    psms = []
    with open(INPUT_POUT, "r") as in_xml:
        in_psms = False
        retrieved_info = defaultdict(lambda: list())
        for l in in_xml:
            l = l.strip()
            if l.startswith("<psms>"):
                # We have read all available psms
                in_psms = True
                continue

            if in_psms and l.startswith("</psm>"):
                # We have read all information for a psm, save!
                psms.append(retrieved_info)
                retrieved_info = defaultdict(lambda: list())

            if in_psms and l.startswith("</psms>"):
                # We have read all available psms
                break

            if in_psms and l.startswith("<psm "):
                # We have read all available psms
                retrieved_info["psm_id"].append(l[l.index("\"")+1:l.index("\" p:decoy")])
                continue


            if in_psms and l.startswith("<peptide_seq"):
                # We have read all available psms
                # We also remove any modifications from the sequence
                peptide = l[l.index("seq=\"")+5:-3]
                while "[" in peptide:
                    peptide = peptide[0:peptide.index("[")] + peptide[peptide.index("]")+1:]
                retrieved_info["peptide_seq"].append(peptide)
                continue


            if in_psms and l.startswith("<"):
                # We have read an attribute not in the tag
                key = l[l.index("<")+1:l.index(">")]
                retrieved_info[key].append(
                    l[len(key)+2:-(len(key)+3)]
                )
                continue

    # Export dict as tsv_value file and rename header
    header = list(psms[0].keys())
    psms_lists = []
    for p in psms:
        new_entry = []
        for h in header:
            new_entry.append(",".join(p[h]))
        psms_lists.append(new_entry)
    

    header[header.index("peptide_seq")] = "plain_peptide"
    header[header.index("protein_id")] = "protein"

    with open(OUTPUT_TSV, "w") as out_file:
        csv_out = csv.writer(out_file, delimiter="\t")
        csv_out.writerows([header, *psms_lists])
