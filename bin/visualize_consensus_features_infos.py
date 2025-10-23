#!/bin/env python

import argparse
import os
import sys
import csv
csv.field_size_limit(sys.maxsize)

from ast import literal_eval

import plotly.graph_objects as go
import pandas as pd
import numpy as np


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-tsv_file", help="The final tsv file UnbeQuant produces (with all columns)")
    parser.add_argument("-output_folder_plots", help="The featureXML file, containing the user param about the MS2 scans", default=".")
    parser.add_argument("-minlh_parameter", help="The minlh parameter used for the quantification", type=int, default=12)
    parser.add_argument("-minlh_up_to", help="The maximum minlh value to be considered for the estimation", type=int, default=100)

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

    # Count the missing values and do a bar plot
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


    ### PLOT minlh ESTIMATIONS


    with open(args.tsv_file, 'r') as f:
        # Read TSV and get relevant columns
        reader = csv.reader(f, delimiter='\t')
        all_headers = next(reader)
        rts_idcs = [idx for idx, h in enumerate(all_headers) if h.endswith("_____l_retention_times")]
        files = [all_headers[rt_idx].replace("_____l_retention_times", "") for rt_idx in rts_idcs]
        ms2s_idcs = [all_headers.index(f + "_____l_ms2_scans") for f in files]
        raw_peps_idcs = [all_headers.index(f + "_____l_raw_pep_ident") for f in files]

        # Get data for minlh estimation
        all_count = np.zeros(args.minlh_up_to, dtype=int)
        all_count_first_iso = np.zeros(args.minlh_up_to, dtype=int)
        ms2s_count = np.zeros(args.minlh_up_to, dtype=int)
        ms2s_count_first_iso = np.zeros(args.minlh_up_to, dtype=int)
        no_ms2s_count = np.zeros(args.minlh_up_to, dtype=int)
        no_ms2s_count_first_iso = np.zeros(args.minlh_up_to, dtype=int)
        ident_count = np.zeros(args.minlh_up_to, dtype=int)
        ident_count_first_iso = np.zeros(args.minlh_up_to, dtype=int)
        unident_count = np.zeros(args.minlh_up_to, dtype=int)
        unident_count_first_iso = np.zeros(args.minlh_up_to, dtype=int)
        count_traces = 0
        count_features = 0
        for row in reader:
            for rt_idx, ms2s_idx, raw_pep_idx in zip(rts_idcs, ms2s_idcs, raw_peps_idcs):
                # Get RTs
                if row[rt_idx]:
                    rts = literal_eval(row[rt_idx])
                    count_features += 1
                    # Check if it contains an MS2
                    is_with_ms2 = False
                    if row[ms2s_idx]:
                        
                        ms2s_set = set(literal_eval(row[ms2s_idx]))
                        is_with_ms2 = len(ms2s_set) > 0
                    # Check if it is identified:
                    is_ident = False
                    if row[raw_pep_idx]:
                        raw_pep_set = set(literal_eval(row[raw_pep_idx]))
                        is_ident = len(raw_pep_set) > 0

                    # For all traces
                    all_count_first_iso[:len(rts[0])+1] += 1
                    for rt_list in rts:
                        count_traces += 1
                        all_count[:len(rt_list)+1] += 1

                    # For MS2s
                    if is_with_ms2:
                        ms2s_count_first_iso[:len(rts[0])+1] += 1
                        for rt_list in rts:
                            ms2s_count[:len(rt_list)+1] += 1
                    else:
                        no_ms2s_count_first_iso[:len(rts[0])+1] += 1
                        for rt_list in rts:
                            no_ms2s_count[:len(rt_list)+1] += 1
                    
                    # For identification status
                    if is_ident:
                        ident_count_first_iso[:len(rts[0])+1] += 1
                        for rt_list in rts:
                            ident_count[:len(rt_list)+1] += 1
                    else:
                        unident_count_first_iso[:len(rts[0])+1] += 1
                        for rt_list in rts:
                            unident_count[:len(rt_list)+1] += 1

        # Plot MS2s
        df = pd.DataFrame({
            "All Traces": all_count[:100],
            "Traces with MS2": ms2s_count[:100],
            "Traces without MS2": no_ms2s_count[:100],
            "index": range(len(all_count))[:100],
        })
        melt_df = df.melt(id_vars=["index"])

        fig = px.scatter(melt_df, x="index", y="value", color="variable",
                        title=f"Estimated 'minlh' for > {args.minlh_parameter} split by MS2 status. Number of found traces: {count_traces}. Number of found features: {count_features}",
                        labels={"index": "minlh-Parameter", "value": "#Traces", "variable": "MS2 Status"}
        )
        fig.add_vline(x=args.minlh_parameter, line_dash="dash", line_color="red", annotation_position="top right",
            annotation_text=f"minlh-parameter: {args.minlh_parameter} ",
            annotation_font_color="red"
        )
        to_pt = np.argmax(ms2s_count > no_ms2s_count)
        if to_pt > 0:
            fig.add_vline(x=to_pt, line_dash="dash", line_color="blue",
                annotation_text=f"Trade-Off at {to_pt} ", annotation_position="bottom right",
                annotation_font_color="blue"
            )
        fig.write_html(os.path.join(args.output_folder_plots, "unbequant_minlh_estimations_ms2s.html"))

        # Plot Identifications
        df = pd.DataFrame({
            "All Traces": all_count[:100],
            "Identified Traces": ident_count[:100],
            "Unidentified Traces": unident_count[:100],
            "index": range(len(all_count))[:100],
        })
        melt_df = df.melt(id_vars=["index"])

        fig = px.scatter(melt_df, x="index", y="value", color="variable",
                        title=f"Estimated 'minlh' for > {args.minlh_parameter} split by Identification status. Number of found traces: {count_traces}. Number of found features: {count_features}",
                        labels={"index": "minlh-Parameter", "value": "#Traces", "variable": "Identification Status"}
        )
        fig.add_vline(x=args.minlh_parameter, line_dash="dash", line_color="red", annotation_position="top right",
            annotation_text=f"minlh-parameter: {args.minlh_parameter} ",
            annotation_font_color="red"
        )
        to_pt = np.argmax(ident_count > unident_count)
        if to_pt > 0:
            fig.add_vline(x=to_pt, line_dash="dash", line_color="blue",
                annotation_text=f"Trade-Off at {to_pt} ", annotation_position="bottom right",
                annotation_font_color="blue"
            )
        fig.write_html(os.path.join(args.output_folder_plots, "unbequant_minlh_estimations_identifications.html"))
