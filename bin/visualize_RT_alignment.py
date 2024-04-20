#!/bin/env python

import sys
import csv
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
    parser.add_argument("-trafo_xmls", help="The trafo_xmls, seperated by a ',' which should be plotted ")
    parser.add_argument("-out_file_postfix", help="Postfix for the file_output", default="")
    return parser.parse_args()


if __name__ == "__main__":
    # Parse Arguments
    args = parse_args()

    trafo_info = dict()
    for trafo in args.trafo_xmls.split(","):
        # Load Transformation Information
        trdesc = pyopenms.TransformationDescription()
        pyopenms.TransformationXMLFile().load(trafo, trdesc, False)

        # Save Information
        data_points = [(x.first, x.second) for x in trdesc.getDataPoints()]
        trafo_info[os.path.basename(trafo).split("_____", 1)[0]] = (
            [x[0] for x in data_points], # From RT
            [x[1] for x in data_points]  # To RT
        )

    # Transform to pandas frame
    df = pd.DataFrame([(x,z1, z2) for x,y in trafo_info.items() for z1, z2 in zip(y[0],y[1])], columns=["File", "Original RT", "Transformed RT"])

    # Plot Information in one plot
    fig = ex.scatter(df, x="Original RT", y="Transformed RT", color="File")
    plotly.offline.plot(fig, filename="RT_transformation_all_measurements_" + args.out_file_postfix + ".html", auto_open=False)
    fig.write_image("RT_transformation_all_measurements_" + args.out_file_postfix + ".png")

    # Plot Information in single plots
    fig = make_subplots(cols=5, rows=math.ceil((len(trafo_info)) /5 ), subplot_titles=list(trafo_info.keys()))
    row_idx = 1
    col_idx = 1
    for file, entries in trafo_info.items():
        fig.add_trace(
            go.Scatter(
                x=entries[0],
                y=entries[1],
                mode='markers'
            ), 
            row=row_idx, col=col_idx
        )
        if col_idx == 5:
            row_idx+=1
            col_idx=1
        else:
            col_idx+=1

    fig.update_layout(
        title='title',
        autosize=True,
        height=400 * math.ceil((len(trafo_info)) /5 )
    )
    fig.update_xaxes(rangeslider=dict(visible=False))
    plotly.offline.plot(fig, filename="RT_transformation_all_measurements_single_plots_" + args.out_file_postfix + ".html", auto_open=False)
    plotly.io.write_json(fig, "RT_transformation_all_measurements_single_plots_" + args.out_file_postfix + ".json")
    fig.write_image("RT_transformation_all_measurements_single_plots_" + args.out_file_postfix + ".png")
