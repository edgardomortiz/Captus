#!/usr/bin/env python3
"""
Copyright 2020 Gentaro Shigita (gentaro.shigita@tum.de)
https://github.com/edgardomortiz/Captus

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""


import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from . import settings_assembly
from .misc import dim, elapsed_time, red


def build_qc_report(out_dir, qc_extras_dir):
    start = time.time()


    ### Stats on reads/bases ###
    df = pd.read_table(Path(qc_extras_dir, settings_assembly.QC_FILES["REBA"]),
                       usecols=[*range(0,3),4,5,7,8,10,11])

    df["reads_input_%"] = 100
    df["bases_input_%"] = 100
    df["reads_passed_round1_%"] = df["reads_passed_round1"] / df["reads_input"] * 100
    df["bases_passed_round1_%"] = df["bases_passed_round1"] / df["bases_input"] * 100
    df["reads_passed_round2_%"] = df["reads_passed_round2"] / df["reads_input"] * 100
    df["bases_passed_round2_%"] = df["bases_passed_round2"] / df["bases_input"] * 100
    df["reads_passed_cleaning_%"] = df["reads_passed_cleaning"] / df["reads_input"] * 100
    df["bases_passed_cleaning_%"] = df["bases_passed_cleaning"] / df["bases_input"] * 100

    sample_list = df["sample"].unique()
    var_suffix_list = ["_input", "_passed_round1", "_passed_round2", "_passed_cleaning"]
    legend_list = ["Input", "Round1", "Round2", "Cleaned"]
    colors = ["#CC79A7", "#E69F00", "#009E73", "#56B4E9"]

    fig1 = make_subplots(
        cols=2,
        shared_yaxes=True,
        horizontal_spacing=0.02,
        subplot_titles=["Reads", "Bases"],
    )
    for var, name, color in zip(var_suffix_list, legend_list, colors):
        # For reads
        fig1.append_trace(
            go.Bar(
                x=df["reads" + var],
                y=df["sample"],
                orientation="h",
                marker_color=color,
                marker_line_color="rgb(8,8,8)",
                marker_line_width=0.25,
                legendgroup=name,
                name=name,
                customdata=df,
                hovertemplate="<b>%{y}</b><br>" +
                              "Count: %{x:,.0f} reads",
            ),
            row=1,
            col=1
        )
        # For bases
        fig1.append_trace(
            go.Bar(
                x=df["bases" + var],
                y=df["sample"],
                orientation="h",
                marker_color=color,
                marker_line_color="rgb(8,8,8)",
                marker_line_width=0.25,
                legendgroup=name,
                name=name,
                customdata=df,
                hovertemplate="<b>%{y}</b><br>" +
                              "Count: %{x:,.0f} bases",
                showlegend=False,
            ),
            row=1,
            col=2
        )

    buttons = [
        dict(label="Count",
             method="update",
             args=[
                 dict(
                     x=[df["reads_input"],
                        df["bases_input"],
                        df["reads_passed_round1"],
                        df["bases_passed_round1"],
                        df["reads_passed_round2"],
                        df["bases_passed_round2"],
                        df["reads_passed_cleaning"],
                        df["bases_passed_cleaning"]],
                     hovertemplate=["<b>%{y}</b><br>Count: %{x:,.0f} reads",
                                    "<b>%{y}</b><br>Count: %{x:,.0f} bases"] * 4,
                 ),
                 dict(
                     xaxis=dict(
                         title="Count",
                         domain=[0, 0.49],
                         ticks="outside",
                         gridcolor="rgb(64,64,64)",
                     ),
                     xaxis2=dict(
                         title="Count",
                         domain=[0.51, 1],
                         ticks="outside",
                         gridcolor="rgb(64,64,64)",
                     ),
                 ),
             ],
        ),
        dict(label="Percentage",
             method="update",
             args=[
                 dict(
                     x=[df["reads_input_%"],
                        df["bases_input_%"],
                        df["reads_passed_round1_%"],
                        df["bases_passed_round1_%"],
                        df["reads_passed_round2_%"],
                        df["bases_passed_round2_%"],
                        df["reads_passed_cleaning_%"],
                        df["bases_passed_cleaning_%"]],
                     hovertemplate=["<b>%{y}</b><br>Proportion: %{x:.2f}%"] * 8,
                 ),
                 dict(
                     xaxis=dict(
                         title="Proportion (%)",
                         range=[0, 100],
                         domain=[0, 0.49],
                         ticks="outside",
                         gridcolor="rgb(64,64,64)",
                     ),
                     xaxis2=dict(
                         title="Proportion (%)",
                         range=[0,100],
                         domain=[0.51, 1],
                         ticks="outside",
                         gridcolor="rgb(64,64,64)",
                     ),
                 ),
             ],
        ),
    ]
    updatemenus = [
        dict(
            buttons=buttons,
            type="buttons",
            direction="right",
            pad={"t": 10, "b": 10},
            showactive=True,
            x=1,
            xanchor="right",
            y=1,
            yanchor="bottom",
        ),
    ]

    fig1.update_layout(
        font_family="Arial",
        title="<b>Stats on Reads/Bases</b>",
        yaxis=dict(title="Sample"),
        barmode="overlay",
        bargap=0,
        bargroupgap=0.1,
        legend_tracegroupgap=0,
        height=180 + 15 * len(sample_list),
        plot_bgcolor="rgb(8,8,8)",
        updatemenus=updatemenus,
    )
    fig1.update_xaxes(
        title="Count",
        ticks="outside",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig1.update_yaxes(
        type="category",
        categoryorder="category descending",
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(8,8,8)",
    )


    ### Per Base Quality ###
    df = pd.read_table(Path(qc_extras_dir, settings_assembly.QC_FILES["PBSQ"]))
    df["stage"] = df["stage"].str.capitalize()
    var_list = [
        "mean",
        "percentile_90",
        "lower_quartile",
        "median",
        "upper_quartile",
        "percentile_10",
    ]
    var_lab_list = [
        "Mean",
        "90th Percentile",
        "75th Percentile",
        "50th Percentile",
        "25th Percentile",
        "10th Percentile",
    ]

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig2 = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig2.append_trace(
            go.Heatmap(
                x=df_R1["base"],
                y=[df_R1["sample_name"], df_R1["stage"]],
                z=df_R1["mean"],
                customdata=df_R1,
                hovertemplate="<br>".join([
                    "<b>%{y}</b>",
                    "Position: %{x} bp",
                    "Mean: %{customdata[4]:.0f}",
                    "90<sub>th</sub> percentile: %{customdata[9]}",
                    "75<sub>th</sub> percentile: %{customdata[7]}",
                    "50<sub>th</sub> percentile: %{customdata[5]}",
                    "25<sub>th</sub> percentile: %{customdata[6]}",
                    "10<sub>th</sub> percentile: %{customdata[8]}<extra></extra>",
                ]),
                coloraxis="coloraxis",
                hoverongaps=False,
                ygap=0.75,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig2.append_trace(
            go.Heatmap(
                x=df_R2["base"],
                y=[df_R2["sample_name"], df_R2["stage"]],
                z=df_R2["mean"],
                customdata=df_R2,
                hovertemplate="<br>".join([
                    "<b>%{y}</b>",
                    "Position: %{x} bp",
                    "Mean: %{customdata[4]:.0f}",
                    "90<sub>th</sub> percentile: %{customdata[9]}",
                    "75<sub>th</sub> percentile: %{customdata[7]}",
                    "50<sub>th</sub> percentile: %{customdata[5]}",
                    "25<sub>th</sub> percentile: %{customdata[6]}",
                    "10<sub>th</sub> percentile: %{customdata[8]}<extra></extra>",
                ]),
                coloraxis="coloraxis",
                hoverongaps=False,
                ygap=0.75,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig2 = go.Figure()
        fig2.add_trace(
            go.Heatmap(
                x=df["base"],
                y=[df["sample_name"], df["stage"]],
                z=df["mean"],
                customdata=df,
                hovertemplate="<br>".join([
                    "<b>%{y}</b>",
                    "Position: %{x} bp",
                    "Mean: %{customdata[4]:.0f}",
                    "90<sub>th</sub> percentile: %{customdata[9]}",
                    "75<sub>th</sub> percentile: %{customdata[7]}",
                    "50<sub>th</sub> percentile: %{customdata[5]}",
                    "25<sub>th</sub> percentile: %{customdata[6]}",
                    "10<sub>th</sub> percentile: %{customdata[8]}<extra></extra>",
                ]),
                coloraxis="coloraxis",
                hoverongaps=False,
                ygap=0.75,
            ),
        )

    # Draw boundaries between samples
    y = 1.5
    while y < (len(sample_list) - 1) * 2:
        fig2.add_hline(y=y, line_width=2, line_color="rgb(8,8,8)")
        y = y + 2

    buttons = []
    for var, lab in zip(var_list, var_lab_list):
        if "R2" in df["read"].to_list():
            z=[df_R1[var], df_R2[var]]
        else:
            z=[df[var]]
        button = dict(
            label=lab,
            method="restyle",
            args=[dict(z=z)]
        )
        buttons.append(button)

    updatemenus = [
        dict(
            buttons=buttons,
            type="dropdown",
            direction="down",
            pad=dict(t=10, b=10),
            x=1,
            xanchor="right",
            y=1,
            yanchor="bottom",
        ),
    ]

    fig2.update_layout(
        font_family="Arial",
        title_text="<b>Per Base Quality</b>",
        plot_bgcolor="rgb(8,8,8)",
        yaxis=dict(title="Sample - Stage"),
        coloraxis=dict(
            colorscale="Spectral",
            cmin=2,
            cmax=41,
            colorbar=dict(
                title="PHRED<br>Score",
                lenmode="pixels",
                len=200,
                outlinecolor="rgb(8,8,8)",
                outlinewidth=1,
                ticks="outside",
                yanchor="top" if len(sample_list) > 7 else "middle",
                y=1 if len(sample_list) > 7 else 0.5,
            ),
        ),
        height=180 + 30 * len(sample_list),
        updatemenus=updatemenus,
    )
    fig2.update_xaxes(
        title="Position (bp)",
        ticks="outside",
        tickson="labels",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig2.update_yaxes(
        showdividers=False,
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(64,64,64)",
        autorange="reversed",
    )


    ### Per Read Quality ###
    df = pd.read_table(Path(qc_extras_dir, settings_assembly.QC_FILES["PSQS"]))
    df["stage"] = df["stage"].str.capitalize()
    df_pivot = df.pivot(
        index=["sample_name", "read", "stage"], columns="quality", values="count"
    )
    col = 0
    while df_pivot.iloc[:,col].isnull().any() == True:
        df_pivot.iloc[:,col].fillna(0, inplace=True)
        col += 1
    df = df_pivot.reset_index().melt(
        id_vars=["sample_name", "read", "stage"], value_name="count").sort_values(
            by=["sample_name", "stage"], ascending=[True, False]
    )
    df_grouped = df.groupby(["sample_name", "read", "stage"], as_index=False)["count"].sum()
    df_merged = pd.merge(df, df_grouped, on=["sample_name", "read", "stage"], how="outer")
    df_merged["freq"] = df_merged["count_x"] / df_merged["count_y"] * 100
    df = df_merged

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig3 = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig3.add_trace(
            go.Heatmap(
                x=df_R1["quality"],
                y=[df_R1["sample_name"], df_R1["stage"]],
                z=df_R1["freq"],
                coloraxis="coloraxis",
                customdata=df_R1,
                hovertemplate="<b>%{y}</b><br>" +
                              "Mean PHRED Score: %{x}<br>" +
                              "Proportion: %{z:.2f}%<br>" +
                              "Count: %{customdata[4]:,.0f} reads<extra></extra>",
                hoverongaps=False,
                ygap=0.75,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig3.add_trace(
            go.Heatmap(
                x=df_R2["quality"],
                y=[df_R2["sample_name"], df_R2["stage"]],
                z=df_R2["freq"],
                coloraxis="coloraxis",
                customdata=df_R2,
                hovertemplate="<b>%{y}</b><br>" +
                              "Mean PHRED Score: %{x}<br>" +
                              "Proportion: %{z:.2f}%<br>" +
                              "Count: %{customdata[4]:,.0f} reads<extra></extra>",
                hoverongaps=False,
                ygap=0.75,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig3 = go.Figure()
        fig3.add_trace(
            go.Heatmap(
                x=df["quality"],
                y=[df["sample_name"], df["stage"]],
                z=df["freq"],
                coloraxis="coloraxis",
                customdata=df,
                hovertemplate="<b>%{y}</b><br>" +
                              "Mean PHRED Score: %{x}<br>" +
                              "Proportion: %{z:.2f}%<br>" +
                              "Count: %{customdata[4]:,.0f} reads<extra></extra>",
                ygap=0.75,
            )
        )

    # Draw boundaries between samples
    y = 1.5
    while y < (len(sample_list) - 1) * 2:
        fig3.add_hline(y=y, line_width=2, line_color="rgb(8,8,8)")
        y = y + 2

    fig3.update_layout(
        font_family="Arial",
        title_text="<b>Per Read Quality</b>",
        yaxis=dict(title="Sample - Stage"),
        coloraxis=dict(
            colorscale="Spectral_r",
            colorbar=dict(
                title="Proportion",
                lenmode="pixels",
                len=200,
                outlinecolor="rgb(8,8,8)",
                outlinewidth=1,
                ticks="outside",
                ticksuffix="%",
                yanchor="top" if len(sample_list) > 7 else "middle",
                y=1 if len(sample_list) > 7 else 0.5,
            ),
        ),
        height=180 + 30 * len(sample_list),
        plot_bgcolor="rgb(8,8,8)",
    )
    fig3.update_xaxes(
        title="PHRED Score",
        ticks="outside",
        matches="x",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
        range=[df["quality"].min() - 0.5, df["quality"].max() + 0.5],
    )
    fig3.update_yaxes(
        autorange="reversed",
        showdividers=False,
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(64,64,64)",
    )


    ### Read Length Distribution ###
    df = pd.read_table(Path(qc_extras_dir, settings_assembly.QC_FILES["SLEN"]))
    df = df[df["stage"] == "after"]
    df_pivot = df.pivot(
        index=["sample_name", "read"], columns="length", values="count"
    )
    col = 0
    while df_pivot.iloc[:,col].isnull().any() == True:
        df_pivot.iloc[:,col].fillna(0, inplace=True)
        col = col + 1
    df = df_pivot.reset_index().melt(
        id_vars=["sample_name", "read"], value_name="count").sort_values(
            by="sample_name", ascending=True
    )
    df_grouped = df.groupby(["sample_name", "read"], as_index=False)["count"].sum()
    df_merged = pd.merge(df, df_grouped, on=["sample_name", "read"], how="outer")
    df_merged["freq"] = df_merged["count_x"] / df_merged["count_y"] * 100
    df = df_merged

    colorscale = [
        [0,    "#5E4FA2"],
        [0.01, "#3683BB"],
        [0.05, "#5DB7A9"],
        [0.1,  "#98D6A4"],
        [0.2,  "#D1EC9C"],
        [0.3,  "#F4FAAD"],
        [0.4,  "#FFF1A7"],
        [0.5,  "#FECE7C"],
        [0.6,  "#FB9C59"],
        [0.7,  "#EE6445"],
        [0.8,  "#D0384E"],
        [1,    "#9E0142"],
    ]

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig4 = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig4.add_trace(
            go.Heatmap(
                x=df_R1["length"],
                y=df_R1["sample_name"],
                z=df_R1["freq"],
                coloraxis="coloraxis",
                customdata=df_R1,
                hovertemplate="<b>%{y}</b><br>" +
                              "Length: %{x} bp<br>" +
                              "Proportion: %{z:.2f}%<br>" +
                              "Count: %{customdata[3]:,.0f} reads<extra></extra>",
                hoverongaps=False,
                ygap=2,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig4.add_trace(
            go.Heatmap(
                x=df_R2["length"],
                y=df_R2["sample_name"],
                z=df_R2["freq"],
                coloraxis="coloraxis",
                customdata=df_R2,
                hovertemplate="<b>%{y}</b><br>" +
                              "Length: %{x} bp<br>" +
                              "Proportion: %{z:.2f}%<br>" +
                              "Count: %{customdata[3]:,.0f} reads<extra></extra>",
                hoverongaps=False,
                ygap=2,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig4 = go.Figure()
        fig4.add_trace(
            go.Heatmap(
                x=df["length"],
                y=df["sample_name"],
                z=df["freq"],
                coloraxis="coloraxis",
                customdata=df,
                hovertemplate="<b>%{y}</b><br>" +
                              "Length: %{x} bp<br>" +
                              "Proportion: %{z:.2f}%<br>" +
                              "Count: %{customdata[3]:,.0f} reads<extra></extra>",
                ygap=2,
            ),
        )

    fig4.update_layout(
        font_family="Arial",
        title_text="<b>Read Length Distribution</b>",
        yaxis=dict(title="Sample"),
        coloraxis=dict(
            colorscale=colorscale,
            colorbar=dict(
                title="Proportion",
                lenmode="pixels",
                len=200,
                outlinecolor="rgb(8,8,8)",
                outlinewidth=1,
                ticks="outside",
                ticksuffix="%",
                yanchor="top" if len(sample_list) > 7 else "middle",
                y=1 if len(sample_list) > 7 else 0.5,
            ),
        ),
        height=180 + 15 * len(sample_list),
        plot_bgcolor="rgb(8,8,8)",
    )
    fig4.update_xaxes(
        title="Read Length (bp)",
        matches="x",
        ticks="outside",
        tickson="labels",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig4.update_yaxes(
        autorange="reversed",
        showdividers=False,
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(64,64,64)",
    )


    ### Per Base Nucleotide Content ###
    df = pd.read_table(Path(qc_extras_dir, settings_assembly.QC_FILES["PBSC"]))
    df["stage"] = df["stage"].str.capitalize()
    df_melt = df.melt(
        id_vars=["sample_name", "read", "stage", "base"],
        value_vars=["A", "T", "G", "C"],
        var_name="nuc",
        value_name="pct",
    )
    df_melt["stage-nuc"] = df_melt["stage"] + " - " + df_melt["nuc"]
    df = df_melt.sort_values(by=["sample_name", "stage", "base"], ascending=[True, False, True])

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig5 = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig5.add_trace(
            go.Heatmap(
                x=df_R1["base"],
                y=[df_R1["sample_name"], df_R1["stage-nuc"]],
                z=df_R1["pct"],
                name="Read 1",
                coloraxis="coloraxis",
                hovertemplate="<b>%{y}</b><br>" +
                              "Position: %{x} bp<br>" +
                              "Proportion: %{z:.2f}%<extra></extra>",
                hoverongaps=False,
                ygap=0.25,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig5.add_trace(
            go.Heatmap(
                x=df_R2["base"],
                y=[df_R2["sample_name"], df_R2["stage-nuc"]],
                z=df_R2["pct"],
                name="Read 2",
                coloraxis="coloraxis",
                hovertemplate="<b>%{y}</b><br>" +
                              "Position: %{x} bp<br>" +
                              "Proportion: %{z:.2f}%<extra></extra>",
                hoverongaps=False,
                ygap=0.25,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig5 = go.Figure()
        fig5.add_trace(
            go.Heatmap(
                x=df["gc_content"],
                y=[df["sample_name"], df["stage"]],
                z=df["freq"],
                hovertemplate="<b>%{y}</b><br>" +
                              "Position: %{x} bp<br>" +
                              "Proportion: %{z:.2f}%<extra></extra>",
                ygap=0.25,
            ),
        )

    # Draw boundaries between stages
    y = 3.5
    while y < (len(sample_list)) * 8:
        fig5.add_hline(y=y, line_width=0.75, line_color="rgb(8,8,8)")
        y = y + 8

    # Draw boundaries between samples
    y = 7.5
    while y < (len(sample_list) - 1) * 8:
        fig5.add_hline(y=y, line_width=2, line_color="rgb(8,8,8)")
        y = y + 8

    fig5.update_layout(
        font_family="Arial",
        title_text="<b>Per Base Nucleotide Content</b>",
        yaxis=dict(title="Sample - Stage - Nucleotide"),
        coloraxis=dict(
            colorscale="Spectral_r",
            colorbar=dict(
                title="Proportion",
                lenmode="pixels",
                len=200,
                outlinecolor="rgb(8,8,8)",
                outlinewidth=1,
                ticks="outside",
                ticksuffix="%",
                yanchor="top" if len(sample_list) > 1 else "middle",
                y=1 if len(sample_list) > 1 else 0.5,
            ),
        ),
        height=270 + 120 * len(sample_list),
        plot_bgcolor="rgb(8,8,8)",
    )
    fig5.update_xaxes(
        title="Position (bp)",
        linecolor="black",
        ticks="outside",
        matches="x",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig5.update_yaxes(
        autorange="reversed",
        ticks="outside",
        dtick=1,
        tickson="labels",
        showdividers=False,
        gridcolor="rgb(64,64,64)",
    )


    ### Per Read GC Content ###
    df = pd.read_table(Path(qc_extras_dir, settings_assembly.QC_FILES["PSGC"]))
    df["stage"] = df["stage"].str.capitalize()
    df_grouped = df.groupby(["sample_name", "read", "stage"], as_index=False)["count"].sum()
    df_merged = pd.merge(df, df_grouped, on=["sample_name", "read", "stage"], how="outer")
    df_merged["freq"] = df_merged["count_x"] / df_merged["count_y"] * 100
    df = df_merged.drop(columns=["count_x", "count_y"])

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig6 = make_subplots(
            cols=2,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig6.add_trace(
            go.Heatmap(
                x=df_R1["gc_content"],
                y=[df_R1["sample_name"], df_R1["stage"]],
                z=df_R1["freq"],
                coloraxis="coloraxis",
                hovertemplate="<b>%{y}</b><br>" +
                              "GC Content: %{x}%<br>" +
                              "Proportion: %{z:.2f}%<extra></extra>",
                ygap=0.75,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig6.add_trace(
            go.Heatmap(
                x=df_R2["gc_content"],
                y=[df_R2["sample_name"], df_R2["stage"]],
                z=df_R2["freq"],
                coloraxis="coloraxis",
                hovertemplate="<b>%{y}</b><br>" +
                              "GC Content: %{x}%<br>" +
                              "Proportion: %{z:.2f}%<extra></extra>",
                ygap=0.75,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig6 = go.Figure()
        fig6.add_trace(
            go.Heatmap(
                x=df["gc_content"],
                y=[df["sample_name"], df["stage"]],
                z=df["freq"],
                hovertemplate="<b>%{y}</b><br>" +
                              "GC Content: %{x}%<br>" +
                              "Proportion: %{z:.2f}%<extra></extra>",
                ygap=0.75,
            ),
        )

    # Draw boundaries between samples
    y = 1.5
    while y < (len(sample_list) - 1) * 2:
        fig6.add_hline(y=y, line_width=2, line_color="rgb(8,8,8)")
        y = y + 2

    fig6.update_layout(
        font_family="Arial",
        title_text="<b>Per Read GC Content</b>",
        yaxis=dict(title="Sample - Stage"),
        coloraxis=dict(
            colorscale="Spectral_r",
            colorbar=dict(
                title="Proportion",
                lenmode="pixels",
                len=200,
                outlinecolor="rgb(8,8,8)",
                outlinewidth=1,
                ticks="outside",
                ticksuffix="%",
                yanchor="top" if len(sample_list) > 7 else "middle",
                y=1 if len(sample_list) > 7 else 0.5,
            ),
        ),
        height=180 + 30 * len(sample_list),
        plot_bgcolor="rgb(8,8,8)",
    )
    fig6.update_xaxes(
        title="GC Content (%)",
        linecolor="black",
        ticks="outside",
        matches="x",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig6.update_yaxes(
        autorange="reversed",
        linecolor="black",
        ticks="outside",
        dtick=1,
        tickson="labels",
        showdividers=False,
        gridcolor="rgb(64,64,64)",
    )


    ### Sequence Duplication Level ###
    df = pd.read_table(Path(qc_extras_dir, settings_assembly.QC_FILES["SDUP"]),
                       usecols=[*range(0,4), 5])
    df["stage"] = df["stage"].str.capitalize()
    sample_list = df["sample_name"].unique()
    dup_lev_list = df["duplication_level"].unique()
    colors = {
        "1":     "#5e4fa2",
        "2":     "#4175b4",
        "3":     "#439bb5",
        "4":     "#66c2a5",
        "5":     "#94d4a4",
        "6":     "#bfe5a0",
        "7":     "#e6f598",
        "8":     "#f7fcb2",
        "9":     "#fff5ae",
        ">10":   "#fee08b",
        ">50":   "#fdbf6f",
        ">100":  "#fa9857",
        ">500":  "#f46d43",
        ">1k":   "#df4e4b",
        ">5k":   "#c32a4b",
        ">10k+": "#9e0142",
    }
    # For paired-end
    if "R2" in df["read"].to_list():
        fig7 = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        for dup_lev in dup_lev_list:
            # Read 1
            data = df_R1[df_R1["duplication_level"] == dup_lev]
            fig7.add_trace(
                go.Bar(
                    x=data["percentage_of_total"],
                    y=[data["sample_name"], data["stage"]],
                    name=dup_lev,
                    meta=[dup_lev],
                    legendgroup=dup_lev,
                    marker_color=colors[dup_lev],
                    marker_line_color="rgb(8,8,8)",
                    marker_line_width=0.25,
                    orientation="h",
                    hovertemplate="<b>%{y}</b><br>" +
                                  "Duplication Level: %{meta[0]}<br>" +
                                  "Percentage: %{x:.2f}%<extra></extra>",
                ),
                row=1,
                col=1,
            )
            # Read 2
            data = df_R2[df_R2["duplication_level"] == dup_lev]
            fig7.add_trace(
                go.Bar(
                    x=data["percentage_of_total"],
                    y=[data["sample_name"], data["stage"]],
                    name=dup_lev,
                    meta=[dup_lev],
                    legendgroup=dup_lev,
                    marker_color=colors[dup_lev],
                    marker_line_color="rgb(8,8,8)",
                    marker_line_width=0.25,
                    showlegend=False,
                    orientation="h",
                    hovertemplate="<b>%{y}</b><br>" +
                                  "Duplication Level: %{meta[0]}<br>" +
                                  "Percentage: %{x:.2f}%<extra></extra>",
                ),
                row=1,
                col=2,
            )

    # For single-end
    else:
        fig7 = go.Figure()
        for dup_lev in dup_lev_list:
            fig7.add_trace(
                go.Bar(
                    x=df["percentage_of_total"],
                    y=[df["sample_name"], df["stage"]],
                    name=dup_lev,
                    meta=[dup_lev],
                    legendgroup=dup_lev,
                    marker_color=colors[dup_lev],
                    marker_line_color="rgb(8,8,8)",
                    marker_line_width=0.25,
                    orientation="h",
                    hovertemplate="<b>%{y}</b><br>" +
                                  "Duplication Level: %{meta[0]}<br>" +
                                  "Percentage: %{x:.2f}%<extra></extra>",
                ),
            )

    # Draw boundaries between samples
    y = 1.5
    while y < (len(sample_list) - 1) * 2:
        fig7.add_hline(y=y, line_width=2, line_color="rgb(8,8,8)")
        y = y + 2

    fig7.update_layout(
        font_family="Arial",
        title_text="<b>Sequence Duplication Level</b>",
        yaxis=dict(title="Sample - Stage"),
        barmode="stack",
        bargap=0,
        bargroupgap=0.0375,
        height=180 + 30 * len(sample_list),
        legend=dict(title="Duplication<br>Level", tracegroupgap=0),
        plot_bgcolor="rgb(8,8,8)",
    )
    fig7.update_xaxes(
        title="Proportion (%)",
        ticks="outside",
        matches="x",
        range=[0, 100],
        zeroline=False,
        gridcolor="rgb(64,64,64)",
    )
    fig7.update_yaxes(
        autorange="reversed",
        showdividers=False,
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(64,64,64)",
    )


    ### Adapter Content ###
    df = pd.read_table(Path(qc_extras_dir, settings_assembly.QC_FILES["ADCO"]))
    df["stage"] = df["stage"].str.capitalize()
    df["adapter_content"] = df.iloc[:,4:].sum(axis=1)

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig8 = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig8.add_trace(
            go.Heatmap(
                x=df_R1["position"],
                y=[df_R1["sample_name"], df_R1["stage"]],
                z=df_R1["adapter_content"],
                coloraxis="coloraxis",
                name="Read 1",
                customdata=df_R1,
                hovertemplate="<b>%{y}</b><br>" +
                              "Position: %{x} bp<br>" +
                              "Illumina universal adaptor: %{customdata[4]:.2f}%<br>" +
                              "Illumina small RNA 3' adaptor: %{customdata[5]:.2f}%<br>" +
                              "Illumina small RNA 5' adaptor: %{customdata[6]:.2f}%<br>" +
                              "Nextera transposase sequence: %{customdata[7]:.2f}%<br>" +
                              "SOLID small RNA adaptor: %{customdata[8]:.2f}%<br>" +
                              "Total adapter content: %{z:.2f}%<extra></extra>",
                hoverongaps=False,
                ygap=0.75,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig8.add_trace(
            go.Heatmap(
                x=df_R2["position"],
                y=[df_R2["sample_name"], df_R2["stage"]],
                z=df_R2["adapter_content"],
                coloraxis="coloraxis",
                name="Read 2",
                customdata=df_R2,
                hovertemplate="<b>%{y}</b><br>" +
                              "Position: %{x} bp<br>" +
                              "Illumina universal adaptor: %{customdata[4]:.2f}%<br>" +
                              "Illumina small RNA 3' adaptor: %{customdata[5]:.2f}%<br>" +
                              "Illumina small RNA 5' adaptor: %{customdata[6]:.2f}%<br>" +
                              "Nextera transposase sequence: %{customdata[7]:.2f}%<br>" +
                              "SOLID small RNA adaptor: %{customdata[8]:.2f}%<br>" +
                              "Total adapter content: %{z:.2f}%<extra></extra>",
                hoverongaps=False,
                ygap=0.75,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig8 = go.Figure()
        fig8.add_trace(
            go.Heatmap(
                x=df["position"],
                y=[df["sample_name"], df["stage"]],
                z=df["adapter_content"],
                customdata=df_R2,
                hovertemplate="<b>%{y}</b><br>" +
                              "Position: %{x} bp<br>" +
                              "Illumina universal adaptor: %{customdata[4]:.2f}%<br>" +
                              "Illumina small RNA 3' adaptor: %{customdata[5]:.2f}%<br>" +
                              "Illumina small RNA 5' adaptor: %{customdata[6]:.2f}%<br>" +
                              "Nextera transposase sequence: %{customdata[7]:.2f}%<br>" +
                              "SOLID small RNA adaptor: %{customdata[8]:.2f}%<br>" +
                              "Total adapter content: %{z:.2f}%<extra></extra>",
                hoverongaps=False,
                ygap=0.75,
            )
        )

    # Draw boundaries between samples
    y = 1.5
    while y < (len(sample_list) - 1) * 2:
        fig8.add_hline(y=y, line_width=2, line_color="rgb(8,8,8)")
        y = y + 2

    fig8.update_layout(
        font_family="Arial",
        title_text="<b>Adapter Content</b>",
        yaxis=dict(title="Sample - Stage"),
        coloraxis=dict(
            colorscale="Spectral_r",
            colorbar=dict(
                title="Proportion",
                lenmode="pixels",
                len=200,
                outlinecolor="rgb(8,8,8)",
                outlinewidth=1,
                ticks="outside",
                ticksuffix="%",
                yanchor="top" if len(sample_list) > 7 else "middle",
                y=1 if len(sample_list) > 7 else 0.5,
            )
        ),
        height=180 + 30 * len(sample_list),
        plot_bgcolor="rgb(8,8,8)",
    )
    fig8.update_xaxes(
        title="Position (bp)",
        ticks="outside",
        matches="x",
        zeroline=False,
        gridcolor="rgb(64,64,64)",
    )
    fig8.update_yaxes(
        autorange="reversed",
        ticks="outside",
        dtick=1,
        tickson="labels",
        showdividers=False,
        gridcolor="rgb(64,64,64)",
    )

    # Export all plots into one html file
    qc_html_report = Path(out_dir, "captus-assembly_clean.report.html")
    with open(qc_html_report, "w") as f:
        f.write(fig1.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(fig2.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(fig3.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(fig4.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(fig5.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(fig6.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(fig7.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(fig8.to_html(full_html=False, include_plotlyjs="cdn"))
    if qc_html_report.exists() and qc_html_report.is_file():
        qc_html_msg = dim(f"Report generated in {elapsed_time(time.time() - start)}")
    else:
        qc_html_msg = red(f"Report not generated, verify your Python environment")

    return qc_html_report, qc_html_msg


def build_assembly_report(out_dir, asm_stats_tsv):
    start = time.time()

    # Load datatable
    df = pd.read_table(asm_stats_tsv)

    # Variables available as drop-down menu
    var_list = [
        "n_contigs",
        "total_length",
        "N50",
        "avg_length",
        "median_length",
        "longest_contig",
        "shortest_contig",
        [
            "pct_contigs_>=_1kbp",
            "pct_contigs_>=_2kbp",
            "pct_contigs_>=_5kbp",
            "pct_contigs_>=_10kbp",
        ],
        [
            "pct_length_>=_1kbp",
            "pct_length_>=_2kbp",
            "pct_length_>=_5kbp",
            "pct_length_>=_10kbp",
        ],
        "GC_content",
        "avg_depth",
        [
            "pct_contigs_>=_1x",
            "pct_contigs_>=_2x",
            "pct_contigs_>=_5x",
            "pct_contigs_>=_10x",
        ],
    ]

    # Button labels
    button_lab_list = [
        "Number of Contigs",
        "Total Length (bp)",
        "N50 (bp)",
        "Average Length (bp)",
        "Median Length (bp)",
        "Maximum Length (bp)",
        "Minimum Length (bp)",
        "Contig Breakdown by Length (%)",
        "Length Breakdown by Contig Length (%)",
        "GC Content (%)",
        "Average Depth (x)",
        "Contigs Breakdown by Depth (%)",
    ]

    # X axis labels
    xlab_list = [
        "Number of Contigs",
        "Total Length (bp)",
        "N50 (bp)",
        "Average Length (bp)",
        "Median Length (bp)",
        "Maximum Length (bp)",
        "Minimum Length (bp)",
        "Contigs (%)",
        "Length (%)",
        "GC Content (%)",
        "Average Depth (x)",
        "Contigs (%)",
    ]

    colors = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]

    # Create figure
    fig = go.Figure()

    for i in range(4):
        if i == 0:
            visible = True
        else:
            visible = False
        fig.add_trace(
            go.Bar(
                x=df[var_list[0]],
                y=df["sample"],
                orientation="h",
                visible=visible,
                marker_color=colors[i],
                hovertemplate="Sample: %{y}<br>" +
                              xlab_list[0] + ": %{x}<extra></extra>",
            ),
        )

    # Dropdown setting
    buttons = []
    for j in range(len(var_list)):
        if type(var_list[j]) == list:
            x = [
                df[var_list[j][0]],
                df[var_list[j][1]],
                df[var_list[j][2]],
                df[var_list[j][3]],
            ]
            name = [re.sub(".*_>=_", "≥ ", name) for name in var_list[j]]
            visible = [True] * 4
        else:
            x = [df[var_list[j]], None, None, None]
            name = []
            visible = [True, False, False, False]

        hovertemplate = [
            "Sample: %{y}<br>" +
            xlab_list[j] + ": %{x}<extra></extra>"
        ]
        button = dict(
            label=button_lab_list[j],
            method="update",
            args=[
                dict(
                    x=x,
                    name=name,
                    visible=visible,
                    hovertemplate=hovertemplate,
                ),
                dict(
                    xaxis=dict(
                        title=xlab_list[j],
                        showgrid=True,
                    ),
                ),
            ],
        )

        buttons.append(button)

    updatemenus = [
        dict(
            buttons=buttons,
            type = "dropdown",
            direction="down",
            pad={"t": 10, "b": 10},
            showactive=True,
            x=1,
            xanchor="right",
            y=1,
            yanchor="bottom",
        ),
    ]

    # Layout setting
    fig.update_layout(
        plot_bgcolor="rgb(8,8,8)",
        font_family="Arial",
        title="<b>Captus-assembly: Assemble (Assembly Report)</b>",
        xaxis=dict(
            title=xlab_list[0],
            showgrid=True,
            gridcolor="rgb(64,64,64)",
        ),
        yaxis=dict(
            title="Sample",
            type="category",
            autorange="reversed",
            gridcolor="rgb(64,64,64)",
        ),
        barmode="overlay",
        updatemenus=updatemenus,
    )

    # Save plot in HTML
    asm_html_report = Path(out_dir, "captus-assembly_assemble.report.html")
    fig.write_html(asm_html_report)
    if asm_html_report.exists() and asm_html_report.is_file():
        asm_html_msg = dim(f"Report generated in {elapsed_time(time.time() - start)}")
    else:
        asm_html_msg = red(f"Report not generated, verify your Python environment")

    return asm_html_report, asm_html_msg


def build_extraction_report(out_dir, ext_stats_tsv):
    start = time.time()

    # Load datatable
    df = pd.read_table(ext_stats_tsv)

    # Preprocess
    df_best = (
        df.loc[
            df.groupby(["sample_name", "marker_type", "locus"])["lwscore"].idxmax(),:
        ].reset_index(drop=True)
    )
    df_best["hit"] = (
        df.groupby(["sample_name", "marker_type", "locus"], as_index=False).count()["hit"]
    )

    # Define variables
    marker_type = df_best["marker_type"].unique()
    if len(marker_type) > 1:
        marker_type = np.insert(marker_type, 0, "ALL")
    var_list = ["pct_recovered", "pct_identity", "hit", "score", "lwscore"]
    var_lab_list = ["Recovered Length (%)", "Identity (%)", "Hit", "Score", "Length-weighted Score"]
    cmap = "Spectral_r"

    # Create figure
    fig = go.Figure()

    # Make heatmap for each marker type
    buttons1, buttons2 = [], []
    for i in range(len(marker_type)):
        if marker_type[i] == "ALL":
            data = df_best
            data["marker_type - locus"] = data["marker_type"] + " - " + data["locus"].astype(str)
            d = data.pivot(
                index="sample_name",
                columns="marker_type - locus",
                values=var_list[0],
            )
        else:
            data = df_best[df_best["marker_type"] == marker_type[i]]
            d = data.pivot(
                index="sample_name",
                columns="locus",
                values=var_list[0],
            )

        heatmap = go.Heatmap(
            x=d.columns,
            y=d.index,
            z=d.values,
            colorscale=cmap,
            colorbar=dict(
                outlinecolor="rgb(8,8,8)",
                outlinewidth=1,
            ),
            hovertemplate="Sample: %{y}<br>" +
                          "Locus: %{x}<br>" +
                          var_lab_list[0] + ": %{z}<extra></extra>",
            visible=False,
            xgap=0.5,
            ygap=0.5,
        )

        if i == 0:
            heatmap.visible = True
        fig.add_trace(heatmap)

        # Add dropdown for marker types
        visible = [False] * len(marker_type)
        visible[i] = True
        button = dict(
            label=marker_type[i],
            method="update",
            args=[{"visible": visible}],
        )
        buttons1.append(button)

    # Add dropdown for variables
    for j in range(len(var_list)):
        z = []
        for i in range(len(marker_type)):
            if marker_type[i] == "ALL":
                data = df_best
                d = data.pivot(
                    index="sample_name",
                    columns="marker_type - locus",
                    values=var_list[j],
                )
            else:
                data = df_best[df_best["marker_type"] == marker_type[i]]
                d = data.pivot(
                    index="sample_name",
                    columns="locus",
                    values=var_list[j],
                )
            z.append(d.values)
        hovertemplate = [
            "Sample: %{y}<br>" +
            "Locus: %{x}<br>" +
            var_lab_list[j] + ": %{z}<extra></extra>"
        ] * len(marker_type)
        button = dict(
            label=var_lab_list[j],
            method="restyle",
            args=[{"z": z, "hovertemplate": hovertemplate}],
        )
        buttons2.append(button)

    updatemenus = [
        dict(
            buttons=buttons1,
            type = "dropdown",
            direction="down",
            pad={"t": 10, "b": 10},
            showactive=True,
            x=0.5,
            xanchor="right",
            y=1,
            yanchor="bottom",
        ),
        dict(
            buttons=buttons2,
            type = "dropdown",
            direction="down",
            pad={"t": 10, "b": 10},
            showactive=True,
            x=1,
            xanchor="right",
            y=1,
            yanchor="bottom",
        ),
    ]

    # Layout setting
    fig.update_layout(
        plot_bgcolor="rgb(8,8,8)",
        font_family="Arial",
        title="<b>Captus-assembly: Extract (Marker Recovery Report)</b>",
        xaxis=dict(
            title="Locus",
            type="category",
            gridcolor="rgb(64,64,64)",
        ),
        yaxis=dict(
            title="Sample",
            type="category",
            autorange="reversed",
            gridcolor="rgb(64,64,64)",
        ),
        annotations=[
            dict(
                text="<b>Marker:</b>",
                x=0.5,
                xref="paper",
                xanchor="right",
                xshift=-70,
                y=1,
                yref="paper",
                yanchor="top",
                yshift=36,
                align="right",
                showarrow=False,
            ),
            dict(
                text="<b>Variable:</b>",
                x=1,
                xref="paper",
                xanchor="right",
                xshift=-170,
                y=1,
                yref="paper",
                yanchor="top",
                yshift=36,
                align="right",
                showarrow=False,
            ),
        ],
        updatemenus=updatemenus,
    )

    # Save plot in HTML
    ext_html_report = Path(out_dir, "captus-assembly_extract.report.html")
    fig.write_html(ext_html_report)
    if ext_html_report.exists() and ext_html_report.is_file():
        ext_html_msg = dim(f"Report generated in {elapsed_time(time.time() - start)}")
    else:
        ext_html_msg = red(f"Report not generated, verify your Python environment")

    return ext_html_report, ext_html_msg

