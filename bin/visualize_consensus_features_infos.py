#!/bin/env python

import argparse
import os
import csv

from ast import literal_eval
import plotly.graph_objects as go
import pandas as pd


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-tsv_file", help="The final tsv file UnbeQuant produces (with all columns)")
    parser.add_argument("-output_folder_plots", help="The featureXML file, containing the user param about the MS2 scans", default=".")

    return parser.parse_args()

if __name__ == "__main__":
    # Parse Arguments
    args = argparse_setup()

    # Check if the output folder exists
    if not os.path.exists(args.output_folder_plots):
        os.makedirs(args.output_folder_plots)


    ### PLOT MISSING VALUES
    # Read only needed columns from the tsv file
    with open(args.tsv_file, "r") as in_file:
        csv_in = csv.reader(in_file, delimiter="\t")
        header = next(csv_in)
        header_dict = {x: header.index(x) for x in header}
        columns_to_be_read = [x for x in header_dict.keys() if x.endswith("_____intensity") or x.endswith("_____l_pep_ident")]
        filenames = [x[:-14] for x in header_dict.keys() if x.endswith("_____intensity")]
    df = pd.read_csv(args.tsv_file, sep="\t", usecols=columns_to_be_read)

    # Check if identification is there (-1 --> No Entry, 0 No identification, 1 Identification)
    for file in filenames:
        df[file + "_____l_pep_ident"] = df[file + "_____l_pep_ident"].apply(lambda x: len(literal_eval(x)) if not pd.isna(x) else -1)
        df.loc[df[file + "_____l_pep_ident"] > 1, file + "_____l_pep_ident"] = 1


    # Get Number of missing values per row
    df_intens = df[[x for x in df.columns if x.endswith("_____intensity")]]
    missing_values = df_intens.isna().sum(axis=1)

    # Count and do a bar plot of the missing values
    missing_values_count = missing_values.value_counts().sort_index()
    fig = go.Figure()
    fig.add_trace(go.Bar(x=missing_values_count.index, y=missing_values_count.values))
    fig.update_layout(
        title="Missing Values per C.Feature",
        xaxis_title="Number of Missing Values per C.Feature",
        yaxis_title="Count",
    )
    fig.update_layout(updatemenus=[
                dict(
                    buttons=list([
                        dict(label="Log", 
                            method="relayout", 
                            args=[{"yaxis.type": "log"}]),

                        dict(label="Linear",  
                            method="relayout", 
                            args=[{"yaxis.type": "linear"}]),
                    ]),
                )])
    fig.update_yaxes(type="log")
    fig.write_html(os.path.join(args.output_folder_plots, "missing_values_per_CFeature.html"))


    # Get Number of missing values per row for identified and unidentified 
    df["identified"] = df[[x for x in df.columns if x.endswith("_____l_pep_ident")]].replace(-1, 0).sum(axis=1) == 1
    missing_values_ident = df[df["identified"] == True].isna().sum(axis=1)
    missing_values_unident = df[df["identified"] == False].isna().sum(axis=1)

    # Create dataframe for the bar plot
    missing_values_ident_count = missing_values_ident.value_counts().sort_index()
    missing_values_unident_count = missing_values_unident.value_counts().sort_index()
    bar_df = pd.DataFrame(
        {
            "count": pd.concat([missing_values_ident_count, missing_values_unident_count]),
            "x": list(missing_values_ident_count.index) + list(missing_values_unident_count.index),  
            "color": len(missing_values_ident_count) * ["Identified"] + len(missing_values_unident_count) * ["Unidentified"],
        }
    )

    # Color the stacked bars of missing values by either identifications or unidentifications
    import plotly.express as ex
    fig = ex.bar(bar_df, x="x", y="count", color="color")
    fig.update_layout(
        title="Missing Values per C.Feature",
        xaxis_title="Number of Missing Values per C.Feature",
        yaxis_title="Count",
    )
    fig.update_layout(updatemenus=[
                dict(
                    buttons=list([
                        dict(label="Log", 
                            method="relayout", 
                            args=[{"yaxis.type": "log"}]),

                        dict(label="Linear",  
                            method="relayout", 
                            args=[{"yaxis.type": "linear"}]),
                    ]),
                )])
    fig.update_yaxes(type="log")
    fig.write_html(os.path.join(args.output_folder_plots, "missing_values_per_CFeature_by_identification.html"))



