#!/bin/env python

import sys
import csv
csv.field_size_limit(sys.maxsize)
import xlsxwriter


# Generate a more comprehendable overview of found PSMs in a Dataset
input_file = sys.argv[1]
output_file = sys.argv[2]

if __name__ == "__main__":
    # Parse identification results from comet (containing the extra fasta header AND the Source_file_name)
    with open(input_file, "r") as in_file:
        csv_in = csv.reader(in_file, delimiter="\t")

        header = next(csv_in)

        source_file_idx = header.index("source_file")
        peptide_idx = header.index("plain_peptide")
        fasta_desc_idx = header.index("fasta_desc")

        source_files = set()

        # Map peptide to source_file, to count
        peptide_to_souce_file_dict = dict()
        for l in csv_in:
            source_files.add(l[source_file_idx])
            if l[peptide_idx] not in peptide_to_souce_file_dict:
                peptide_to_souce_file_dict[l[peptide_idx]] = dict()

            if l[source_file_idx] not in peptide_to_souce_file_dict[l[peptide_idx]]:
                peptide_to_souce_file_dict[l[peptide_idx]][l[source_file_idx]] = [[], 0]

            peptide_to_souce_file_dict[l[peptide_idx]][l[source_file_idx]][0].append(l)
            peptide_to_souce_file_dict[l[peptide_idx]][l[source_file_idx]][1] += 1

        pass
        source_files = sorted(list(source_files))
        # Summarize and write table
        summarized = [
            (
                sum([x[1] for x in source_count.values()]),  # Number of PSM matches
                len(source_count.keys()),  # Number of matches of RAW-files
                pep,  # Matched peptide
                next(iter(source_count.values()))[0][0][fasta_desc_idx],  # FASTA-Description fitting to peptide
                [source_count[x][1] if x in source_count else 0 for x in source_files]
            ) 
            for pep, source_count in peptide_to_souce_file_dict.items()
        ]


        # Generate a Excel Sheet, containing the counts and found peptides
        flatted_summarized = sorted([list(x[0:4]) + x[4] for x in summarized], key=lambda x: x[0], reverse=True )
        new_header = [
                "#PSMs",
                "#Source_files",
                "plan_peptide",
                "fasta_desc",
                *source_files
        ]

        workbook = xlsxwriter.Workbook(output_file)
        worksheet1 = workbook.add_worksheet()
        worksheet1.freeze_panes(1, 0)
        for j in range(len(new_header)):
            worksheet1.write(0, j, new_header[j])

        for i in range(len(flatted_summarized)):
            for j in range(len(new_header)):
                worksheet1.write(i+1,j, flatted_summarized[i][j])

        worksheet1.conditional_format(1, 4, len(flatted_summarized)+1, len(new_header), {'type': '3_color_scale'})





        workbook.close()
