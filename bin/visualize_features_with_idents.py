#!/bin/env python

import sys
import csv
csv.field_size_limit(sys.maxsize)
import os
import matplotlib
import pyopenms
import argparse
import math

import plotly
import plotly.express as ex
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-feature_tsv", help="The feature_tsv to plot rectangles")
    parser.add_argument("-ident_tsv", help="The identification tsv, containing the retention time and mz values to plot points")
    return parser.parse_args()


if __name__ == "__main__":
    # Parse Arguments
    args = parse_args()

    args.feature_tsv = "/home/luxii/Desktop/DELME/Goldstandard_Metapipeline_without_SPikeIns/source/results/quantifications/features_with_annotated_identifications/C2_FLI04008_____fdr____0.01_fdr_____with_identifications.tsv"
    args.ident_tsv = "/home/luxii/Desktop/DELME/Goldstandard_Metapipeline_without_SPikeIns/source/results/identifications/C2_FLI04008_____perc_fdr_0.01_____qvalue_no_decoys_fdr_0.01.tsv"


    with open(args.feature_tsv, "r") as features_in, open(args.ident_tsv, "r") as idents_in:
        features_in_csv = csv.reader(features_in, delimiter="\t")
        idents_in_csv = csv.reader(idents_in, delimiter="\t")

        print("help mee")
        header_features = next(features_in_csv)

