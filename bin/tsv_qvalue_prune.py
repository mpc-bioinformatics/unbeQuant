#!/bin/env python

import sys
import csv
csv.field_size_limit(sys.maxsize)

input_tsv = sys.argv[1]
output_tsv = sys.argv[2]
limit_fdr =  float(sys.argv[3])


if __name__ == "__main__":
    # Open all files
    with open(output_tsv, "w") as out_file, open(input_tsv, "r") as in_file:
        csv_out = csv.writer(out_file, delimiter="\t")
        csv_in = csv.reader(in_file, delimiter="\t")

        # Read Header, get qval index and write header to new line
        header = next(csv_in)
        idx_qval = header.index("QValue")
        csv_out.writerow(header)

        # Write as long as we do not go over the limit
        # MSGFPlus returns the qvalue already in order!
        for line in csv_in:
            if float(line[idx_qval]) > limit_fdr:
                pass  # Ignore line and do not write
            else:
                csv_out.writerow(line)
