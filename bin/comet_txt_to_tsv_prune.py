#!/bin/env python

import sys
import csv
csv.field_size_limit(sys.maxsize)
import re

fasta = sys.argv[1]
input_txt = sys.argv[2]
output_tsv_base = sys.argv[3]
limit_fdr =  float(sys.argv[4])
percolator_file =  int(sys.argv[5]) # 1 for True


def get_header_from_fasta(protein_ids):
    """ Memory-efficient solution of getting the FASTA-Headers """
    fasta_d = dict()
    with open(fasta, "r") as in_fasta:
        for l in in_fasta:
            if l.startswith(">"):
                line_split = l.split("|", 2)
                if line_split[1] in protein_ids:
                    fasta_d[line_split[1]] = line_split[-1]
    return fasta_d


fasta_header_pattern = r"(.*?[pg|tr|sp|lcl|ref])\|([a-zA-Z0-9\-\_]+)\|"

if __name__ == "__main__":
    # Open all files

    with open(input_txt, "r") as in_file:

        if percolator_file == 0:
            next(in_file) # Skip the comet header
            csv_in = csv.reader(in_file, delimiter="\t")

            header_idx = next(csv_in)
            score_idx = header_idx.index("xcorr")
            protein_idx = header_idx.index("protein")
        else: 
            csv_in = csv.reader(in_file, delimiter="\t")

            header_idx = next(csv_in)
            score_idx = header_idx.index("svm_score")
            protein_idx = header_idx.index("protein")

        entries = [x for x in csv_in]
        entries = sorted(entries, key=lambda x: float(x[score_idx]), reverse=True)
        
        # Calculate decoys
        num_hits = 0
        num_decs = 0

        additional_info_entries = []

        # Get fasta_d
        protein_ids = []
        for entry in entries:
            matches = re.finditer(fasta_header_pattern, entry[protein_idx], re.MULTILINE)
            for matchNum, match in enumerate(matches, start=1):
                protein_ids.append(match.groups()[1])
        fasta_d = get_header_from_fasta(set(protein_ids))

        for entry in entries:
            proteins_hits = entry[protein_idx]
            is_decoy = []
            protein_id = []

            matches = re.finditer(fasta_header_pattern, entry[protein_idx], re.MULTILINE)
            for matchNum, match in enumerate(matches, start=1):
                is_decoy.append(True if "DECOY_" in match.groups()[0] else False)
                protein_id.append(match.groups()[1])

            if all(is_decoy):
                num_decs += 1
            else:
                num_hits += 1 


            additional_info_entries.append(
                [
                    num_decs / (num_decs + num_hits),
                    ",".join([y for x,y in zip(is_decoy, protein_id)]), # Concatenate IDs
                    ",".join([fasta_d[y] for x,y in zip(is_decoy, protein_id)]) # Concatenate Descriptions
                ]
            )
            
        qvalue = additional_info_entries[-1][0]
        for add_in in additional_info_entries[::-1]:
            if add_in[0] > qvalue:
                add_in[0] = qvalue
            else:
                qvalue = add_in[0]

        # Write all with decoys
        with open(output_tsv_base + "_qvalue.tsv", "w") as out_file:
            csv_out = csv.writer(out_file, delimiter="\t")
            csv_out.writerow(header_idx + ["qvalue", "fasta_id", "fasta_desc"])
            for a, b in zip(entries, additional_info_entries):
                csv_out.writerow(
                        a + b
                )

        # Write all without decoys
        with open(output_tsv_base + "_qvalue_no_decoys.tsv", "w") as out_file:
            csv_out = csv.writer(out_file, delimiter="\t")
            csv_out.writerow(header_idx + ["qvalue", "fasta_id", "fasta_desc"])
            for a, b in zip(entries, additional_info_entries):
                if not a[protein_idx].startswith("DECOY_"):
                    csv_out.writerow(
                        a + b
                    )

        # Write all without decoys up to fdr
        with open(output_tsv_base + "_qvalue_no_decoys_fdr_" + str(limit_fdr) + ".tsv", "w") as out_file:
            csv_out = csv.writer(out_file, delimiter="\t")
            csv_out.writerow(header_idx + ["qvalue", "fasta_id", "fasta_desc"])
            for a, b in zip(entries, additional_info_entries):
                if not a[protein_idx].startswith("DECOY_") and b[0] < limit_fdr:
                    csv_out.writerow(
                        a + b
                    )
