#!/bin/env python

import argparse
import csv
import os
import sys

csv.field_size_limit(sys.maxsize)


def _check_if_file_exists(s: str):
    """ checks if a file exists. If not: raise Exception """
    # TODO copied from prot_graph.py
    if os.path.isfile(s):
        return s
    else:
        raise Exception("File '{}' does not exists".format(s))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Graph-Generator for Proteins/Peptides and Exporter to various formats"
    )

    # Statistics file
    parser.add_argument(
        "input_csv", type=_check_if_file_exists, nargs=1,
        help="Input_CSV"
    )

    # Number of entries in csv
    parser.add_argument(
        "--output_csv", "-o", type=str, default="output.csv",
        help="Output csv file"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Parameters
    in_csv_file = args.input_csv[0]
    out_csv_file = args.output_csv

    # Open al files and sort them accordingly
    with open(out_csv_file, "w") as out_file, open(in_csv_file, "r") as in_file:
        # Initialize CSV writer
        csv_out = csv.writer(out_file)
        csv_in = csv.reader(in_file)

        queries = [[float(y) for y in x] for x in csv_in]
        queries = sorted(queries)
        
        num_queries = len(queries)
        index = 0
        while index < len(queries) - 1:
            if queries[index][1] > queries[index+1][0]:  # Overlap
                queries[index][1] = max(queries[index][1], queries[index+1][1])
                del queries[index+1]
            else:
                index += 1

        print("#Queries reduced to: {}%".format(len(queries)*100/num_queries))

        for q in queries:
            csv_out.writerow(q)

