#!/usr/bin/env python3
"""
Copyright 2020-2023 Gentaro Shigita (gentaro.shigita@tum.de)
https://github.com/edgardomortiz/Captus

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""


import itertools
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
from plotly.subplots import make_subplots

from . import settings
from .misc import dim, elapsed_time, red


def normalize(l):
    if len(l) > 0:
        l_min = min(l)
        l_max = max(l)
        return [(i - l_min) / (l_max - l_min) for i in l] if l_max > l_min else 0.5
    else:
        return 0.5

def build_qc_report(out_dir, qc_extras_dir):
    start = time.time()

    ### Summary table ###
    df1 = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["REBA"]))
    df2 = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["PSQS"]))
    df3 = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["SLEN"]))
    df4 = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["PSGC"]))
    df5 = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["ADCO"]))

    # Output reads/bases%
    df1["reads_passed_cleaning_%"] = df1["reads_passed_cleaning"] / df1["reads_input"] * 100
    df1["bases_passed_cleaning_%"] = df1["bases_passed_cleaning"] / df1["bases_input"] * 100
    # Mean read length%
    df3["length * count"] = df3["length"] * df3["count"]
    avg_read_len = (
        df3.groupby(["sample_name", "stage"])["length * count"].sum()
        / df3.groupby(["sample_name", "stage"])["count"].sum()
    )
    read_len_pct = avg_read_len.loc[:, "after"] / avg_read_len.loc[:, "before"] * 100
    # Q20, Q30 reads%
    q20 = (df2[df2["quality"] >= 20].groupby(["sample_name", "stage"])["count"].sum()
        / df2.groupby(["sample_name", "stage"])["count"].sum() * 100
    )
    q30 = (df2[df2["quality"] >= 30].groupby(["sample_name", "stage"])["count"].sum()
        / df2.groupby(["sample_name", "stage"])["count"].sum() * 100
    )
    # GC%
    df4["gc_content * count"] = df4["gc_content"] * df4["count"]
    gc = (df4.groupby(["sample_name", "stage"])["gc_content * count"].sum()
        / df4.groupby(["sample_name", "stage"])["count"].sum()
    )
    # Adapter%
    df5 = df5[df5["stage"] == "before"]
    df5["total_adapter_content%"] = df5.iloc[:,4:].sum(axis=1)
    df5 = df5.groupby(["sample_name", "read"])["total_adapter_content%"].max().reset_index()
    df5 = df5.groupby(["sample_name"])["total_adapter_content%"].mean().reset_index()

    df = pd.DataFrame({
        "Sample": df1["sample"],
        "Input Reads": df1["reads_input"],
        "Input Bases": df1["bases_input"],
        "Output Reads": df1["reads_passed_cleaning"],
        "Output Reads%": df1["reads_passed_cleaning_%"],
        "Output Bases": df1["bases_passed_cleaning"],
        "Output Bases%": df1["bases_passed_cleaning_%"],
        "Mean Read Length%": read_len_pct.reset_index()[0],
        "≥Q20 Reads%": q20.xs("after", level="stage").reset_index()["count"],
        "≥Q30 Reads%": q30.xs("after", level="stage").reset_index()["count"],
        "GC%": gc.xs("after", level="stage").reset_index()[0],
        "Adapter%": df5["total_adapter_content%"],
    })

    sample_list = df["Sample"].unique()

    colorscale = [
        [0, "#F0D9E6"],
        [0.5, "#F5F5F5"],
        [1, "#BDE3D8"]
    ]

    fill_color = []
    for col in df.columns:
        if col == "Sample":
            fill_color.append("#F5F5F5")
        else:
            fill_color.append(
                sample_colorscale(
                    colorscale,
                    normalize(df[col].fillna(0)),
                )
            )

    figs = []
    fig = go.Figure()
    fig.add_trace(go.Table(
        header=dict(values=["<b>" + col + "</b>" for col in df.columns]),
        cells=dict(
            values=[df[col] for col in df.columns],
            fill_color=fill_color,
            format=[None, ",", ",", ",", ".2f", ",", ".2f"],
            align=["left", "right"],
            height=21,
        )
    ))

    buttons = []
    for col in df.columns:
        df.sort_values(
            by=col,
            ascending=True if col == "Sample" else False,
            inplace=True,
        )
        fill_color = []
        for col2 in df.columns:
            if col2 == "Sample":
                fill_color.append("#F5F5F5")
            else:
                fill_color.append(
                    sample_colorscale(
                        colorscale,
                        normalize(df[col2].fillna(0)),
                    )
                )
        button = dict(
            label=col,
            method="restyle",
            args=[
                dict(
                    cells=dict(
                        values=[df[col] for col in df.columns],
                        fill=dict(
                            color=fill_color,
                        ),
                        format=[None, ",", ",", ",", ".2f", ",", ".2f"],
                        align=["left", "right"],
                        height=21,
                )
                )
            ],
        )
        buttons.append(button)

    updatemenus = [dict(
        buttons=buttons,
        type="dropdown",
        direction="down",
        pad=dict(t=10, b=10),
        x=1,
        xanchor="right",
        y=1,
        yanchor="bottom",
    )]

    annotations=[dict(
        text="<b>Sort by:</b>",
        x=1,
        xref="paper",
        xanchor="right",
        xshift=-160,
        y=1,
        yref="paper",
        yanchor="top",
        yshift=36,
        align="right",
        showarrow=False
    )]

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>Captus-assembly: Clean (Quality Control Report)<br>"
            "1. Summary Table</b><br>"
            "<sup>(Source: 03_qc_extras/"
            + ", ".join([
                settings.QC_FILES["REBA"],
                settings.QC_FILES["PSQS"],
                settings.QC_FILES["SLEN"],
                settings.QC_FILES["PSGC"],
                settings.QC_FILES["ADCO"],
            ])
            + ")</sup>"
        ),
        height=230 + 21 * len(sample_list) if len(sample_list) < 31 else None,
        updatemenus=updatemenus,
        annotations=annotations,
    )
    figs.append(fig)

    ### Stats on reads/bases ###
    df = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["REBA"]),
                       usecols=[*range(0,3),4,5,7,8,10,11])

    df["reads_input_%"] = 100
    df["bases_input_%"] = 100
    df["reads_passed_round1_%"] = df["reads_passed_round1"] / df["reads_input"] * 100
    df["bases_passed_round1_%"] = df["bases_passed_round1"] / df["bases_input"] * 100
    df["reads_passed_round2_%"] = df["reads_passed_round2"] / df["reads_input"] * 100
    df["bases_passed_round2_%"] = df["bases_passed_round2"] / df["bases_input"] * 100
    df["reads_passed_cleaning_%"] = df["reads_passed_cleaning"] / df["reads_input"] * 100
    df["bases_passed_cleaning_%"] = df["bases_passed_cleaning"] / df["bases_input"] * 100

    var_suffix_list = ["_input", "_passed_round1", "_passed_round2", "_passed_cleaning"]
    legend_list = ["Input", "Round1", "Round2", "Cleaned"]
    colors = ["#CC79A7", "#E69F00", "#009E73", "#56B4E9"]

    fig = make_subplots(
        cols=2,
        # shared_yaxes=True,
        horizontal_spacing=0.02,
        subplot_titles=["Reads", "Bases"],
    )
    df.sort_values(
        by="bases_passed_cleaning",
        inplace=True,
    )
    for var, name, color in zip(var_suffix_list, legend_list, colors):
        # For reads
        fig.append_trace(
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
                              "Count: <b>%{x:,} reads</b>",
            ),
            row=1,
            col=1
        )
        # For bases
        fig.append_trace(
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
                              "Count: <b>%{x:,} bases</b>",
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
                     y=[df["sample"]],
                     hovertemplate=["<b>%{y}</b><br>Count: <b>%{x:,} reads</b>",
                                    "<b>%{y}</b><br>Count: <b>%{x:,} bases</b>"] * 4,
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
                 )
             ]
        )
    ]
    df.sort_values(
        by="bases_passed_cleaning_%",
        inplace=True,
    )
    buttons.append(
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
                     y=[df["sample"]],
                     hovertemplate=["<b>%{y}</b><br>Proportion: <b>%{x:.2f}%</b>"] * 8,
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
                 )
             ]
        )
    )
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
        )
    ]

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>2. Stats on Reads/Bases</b><br>"
            "<sup>(Source: 03_qc_extras/"
            + settings.QC_FILES["REBA"]
            + ")</sup>"
        ),
        yaxis=dict(title="Sample"),
        yaxis2=dict(showticklabels=False),
        hoverlabel=dict(
            font_color="rgb(64,64,64)",
            bordercolor="rgb(64,64,64)",
        ),
        barmode="overlay",
        bargap=0,
        bargroupgap=0.05,
        legend_tracegroupgap=0,
        height=180 + 15 * len(sample_list),
        plot_bgcolor="rgb(8,8,8)",
        updatemenus=updatemenus,
    )
    fig.update_xaxes(
        title="Count",
        ticks="outside",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig.update_yaxes(
        type="category",
        # categoryorder="category descending",
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(8,8,8)",
    )
    figs.append(fig)

    ### Per Base Quality ###
    df = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["PBSQ"]))
    df["stage"] = df["stage"].str.capitalize()

    # Covert Phred64 to Phred33
    if df["percentile_90"].max() > 42:
        phred64_sample_list = df.query("percentile_90 > 42")["sample_name"].unique()
        phred64_index = df[(df["sample_name"].isin(phred64_sample_list)) & (df["stage"] == "Before")].index
        df.iloc[phred64_index,4:] = df.iloc[phred64_index,4:] - 31

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

    df_pivot = df.pivot_table(
        index=["sample_name", "stage"],
        columns=["read", "base"],
        values=var_list
    )
    df = df_pivot.stack(
        level=["read", "base"],
        dropna=False,
    ).reset_index()
    df.sort_values(
            by=["sample_name", "stage"],
            ascending=[True, False],
            inplace=True,
    )

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig.append_trace(
            go.Heatmap(
                x=df_R1["base"],
                y=[df_R1["sample_name"], df_R1["stage"]],
                z=df_R1["mean"],
                customdata=df_R1,
                hovertemplate="<br>".join([
                    "<b>%{y}</b>",
                    "Position: <b>%{x} bp</b>",
                    "Mean: <b>%{customdata[5]:.0f}</b>",
                    "90<sub>th</sub> percentile: <b>%{customdata[8]}</b>",
                    "75<sub>th</sub> percentile: <b>%{customdata[4]}</b>",
                    "50<sub>th</sub> percentile: <b>%{customdata[6]}</b>",
                    "25<sub>th</sub> percentile: <b>%{customdata[9]}</b>",
                    "10<sub>th</sub> percentile: <b>%{customdata[7]}</b><extra></extra>",
                ]),
                coloraxis="coloraxis",
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig.append_trace(
            go.Heatmap(
                x=df_R2["base"],
                y=[df_R2["sample_name"], df_R2["stage"]],
                z=df_R2["mean"],
                customdata=df_R2,
                hovertemplate="<br>".join([
                    "<b>%{y}</b>",
                    "Position: <b>%{x} bp</b>",
                    "Mean: <b>%{customdata[5]:.0f}</b>",
                    "90<sub>th</sub> percentile: <b>%{customdata[8]}</b>",
                    "75<sub>th</sub> percentile: <b>%{customdata[4]}</b>",
                    "50<sub>th</sub> percentile: <b>%{customdata[6]}</b>",
                    "25<sub>th</sub> percentile: <b>%{customdata[9]}</b>",
                    "10<sub>th</sub> percentile: <b>%{customdata[7]}</b><extra></extra>",
                ]),
                coloraxis="coloraxis",
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                x=df["base"],
                y=[df["sample_name"], df["stage"]],
                z=df["mean"],
                customdata=df,
                hovertemplate="<br>".join([
                    "<b>%{y}</b>",
                    "Position: <b>%{x} bp</b>",
                    "Mean: <b>%{customdata[5]:.0f}</b>",
                    "90<sub>th</sub> percentile: <b>%{customdata[8]}</b>",
                    "75<sub>th</sub> percentile: <b>%{customdata[4]}</b>",
                    "50<sub>th</sub> percentile: <b>%{customdata[6]}</b>",
                    "25<sub>th</sub> percentile: <b>%{customdata[9]}</b>",
                    "10<sub>th</sub> percentile: <b>%{customdata[7]}</b><extra></extra>",
                ]),
                coloraxis="coloraxis",
                hoverongaps=False,
                ygap=1,
            )
        )

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
        )
    ]

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>3. Per Base Quality</b><br>"
            "<sup>(Source: 03_qc_extras/"
            + settings.QC_FILES["PBSQ"]
            + ")</sup>"
        ),
        plot_bgcolor="rgb(8,8,8)",
        yaxis=dict(title="Sample - Stage"),
        coloraxis=dict(
            colorscale="Spectral",
            cmin=0,
            cmax=42,
            colorbar=dict(
                title="Phred<br>Score",
                lenmode="pixels",
                len=200,
                outlinecolor="rgb(8,8,8)",
                outlinewidth=1,
                ticks="outside",
                yanchor="top" if len(sample_list) > 7 else "middle",
                y=1 if len(sample_list) > 7 else 0.5,
            )
        ),
        height=180 + 30 * len(sample_list),
        updatemenus=updatemenus,
    )
    fig.update_xaxes(
        title="Position (bp)",
        ticks="outside",
        tickson="labels",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig.update_yaxes(
        showdividers=False,
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(64,64,64)",
        autorange="reversed",
    )
    figs.append(fig)

    ### Per Read Quality ###
    df = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["PSQS"]))
    df["stage"] = df["stage"].str.capitalize()

    # Convert Phred64 to Phred33
    if df["quality"].max() > 42:
        phred64_sample_list = df.query("quality > 42")["sample_name"].unique()
        phred64_index = df[(df["sample_name"].isin(phred64_sample_list)) & (df["stage"] == "Before")].index
        df.iloc[phred64_index,3] = df.iloc[phred64_index,3] - 31

    df_pivot = df.pivot(
        index=["sample_name", "read", "stage"], columns="quality", values="count"
    )
    col = 0
    while df_pivot.iloc[:,col].isnull().any() == True:
        df_pivot.iloc[:,col].fillna(0, inplace=True)
        col += 1
    df = df_pivot.reset_index().melt(
        id_vars=["sample_name", "read", "stage"],
        value_name="count",
    )
    df_grouped = df.groupby(["sample_name", "read", "stage"], as_index=False)["count"].sum()
    df_merged = pd.merge(df, df_grouped, on=["sample_name", "read", "stage"], how="outer")
    df_merged["freq"] = df_merged["count_x"] / df_merged["count_y"] * 100
    df_pivot = df_merged.pivot_table(
        index=["sample_name", "stage"],
        columns=["read", "quality"],
    )
    df = df_pivot.stack(
        level=["read", "quality"],
        dropna=False,
    ).reset_index()
    df.sort_values(
            by=["sample_name", "stage"],
            ascending=[True, False],
            inplace=True,
    )

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig.add_trace(
            go.Heatmap(
                x=df_R1["quality"],
                y=[df_R1["sample_name"], df_R1["stage"]],
                z=df_R1["freq"],
                coloraxis="coloraxis",
                customdata=df_R1,
                hovertemplate="<b>%{y}</b><br>" +
                              "Mean Phred Score: <b>%{x}</b><br>" +
                              "Proportion: <b>%{z:.2f}%</b><br>" +
                              "Count: <b>%{customdata[4]:,.0f} reads</b><extra></extra>",
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig.add_trace(
            go.Heatmap(
                x=df_R2["quality"],
                y=[df_R2["sample_name"], df_R2["stage"]],
                z=df_R2["freq"],
                coloraxis="coloraxis",
                customdata=df_R2,
                hovertemplate="<b>%{y}</b><br>" +
                              "Mean Phred Score: <b>%{x}</b><br>" +
                              "Proportion: <b>%{z:.2f}%</b><br>" +
                              "Count: <b>%{customdata[4]:,.0f} reads</b><extra></extra>",
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                x=df["quality"],
                y=[df["sample_name"], df["stage"]],
                z=df["freq"],
                coloraxis="coloraxis",
                customdata=df,
                hovertemplate="<b>%{y}</b><br>" +
                              "Mean Phred Score: <b>%{x}</b><br>" +
                              "Proportion: <b>%{z:.2f}%</b><br>" +
                              "Count: <b>%{customdata[4]:,.0f} reads</b><extra></extra>",
                hoverongaps=False,
                ygap=1,
            )
        )

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>4. Per Read Quality</b><br>"
            "<sup>(Source: 03_qc_extras/"
            + settings.QC_FILES["PSQS"]
            + ")</sup>"
        ),
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
    fig.update_xaxes(
        title="Phred Score",
        ticks="outside",
        matches="x",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
        range=[df["quality"].min() - 0.5, df["quality"].max() + 0.5],
    )
    fig.update_yaxes(
        autorange="reversed",
        showdividers=False,
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(64,64,64)",
    )
    figs.append(fig)

    ### Read Length Distribution ###
    df = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["SLEN"]))
    stage_list = df["stage"].unique()
    read_list = df["read"].unique()
    length_range = range(
        df["length"].min(),
        df["length"].max() + 1,
    )
    df_skeleton = pd.DataFrame(
        list(
            itertools.product(
                sample_list,
                stage_list,
                read_list,
                length_range,
            )
        ),
        columns=[
            "sample_name",
            "stage",
            "read",
            "length",
        ]
    )
    df = pd.merge(
        df_skeleton,
        df,
        on=[
            "sample_name",
            "stage",
            "read",
            "length",
        ],
        how="left"
    )
    df_grouped = df.groupby(
        [
            "sample_name",
            "stage",
            "read"
        ],
        as_index=False
    ).agg(total = ("count", np.sum))
    df_merged = pd.merge(df, df_grouped, on=["sample_name", "stage", "read"], how="outer")
    df_merged["freq"] = df_merged["count"] / df_merged["total"] * 100
    df = df_merged
    df["stage"] = df["stage"].str.capitalize()

    colorscale = [
        [0,     "#5E4FA2"],
        [0.002, "#3683BB"],
        [0.006, "#5DB7A9"],
        [0.01,  "#98D6A4"],
        [0.05,  "#D1EC9C"],
        [0.09,  "#F4FAAD"],
        [0.2,   "#FFF1A7"],
        [0.4,   "#FECE7C"],
        [0.6,   "#FB9C59"],
        [0.8,   "#EE6445"],
        [1.0,   "#D0384E"],
        # [1,    "#9E0142"],
    ]

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig.add_trace(
            go.Heatmap(
                x=df_R1["length"],
                y=[df_R1["sample_name"], df_R1["stage"]],
                z=df_R1["freq"],
                coloraxis="coloraxis",
                customdata=df_R1,
                hovertemplate="<b>%{y}</b><br>" +
                              "Length: <b>%{x} bp</b><br>" +
                              "Proportion: <b>%{z:.2f}%</b><br>" +
                              "Count: <b>%{customdata[4]:,.0f} reads</b><extra></extra>",
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig.add_trace(
            go.Heatmap(
                x=df_R2["length"],
                y=[df_R2["sample_name"], df_R2["stage"]],
                z=df_R2["freq"],
                coloraxis="coloraxis",
                customdata=df_R2,
                hovertemplate="<b>%{y}</b><br>" +
                              "Length: <b>%{x} bp</b><br>" +
                              "Proportion: <b>%{z:.2f}%</b><br>" +
                              "Count: <b>%{customdata[4]:,.0f} reads</b><extra></extra>",
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                x=df["length"],
                y=[df["sample_name"], df["stage"]],
                z=df["freq"],
                coloraxis="coloraxis",
                customdata=df,
                hovertemplate="<b>%{y}</b><br>" +
                              "Length: <b>%{x} bp</b><br>" +
                              "Proportion: <b>%{z:.2f}%</b><br>" +
                              "Count: <b>%{customdata[4]:,.0f} reads</b><extra></extra>",
                hoverongaps=False,
                ygap=1,
            )
        )

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>5. Read Length Distribution</b><br>"
            "<sup>(Source: 03_qc_extras/"
            + settings.QC_FILES["SLEN"]
            + ")</sup>"
        ),
        yaxis=dict(title="Sample - Stage"),
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
            )
        ),
        height=180 + 30 * len(sample_list),
        plot_bgcolor="rgb(8,8,8)",
    )
    fig.update_xaxes(
        title="Read Length (bp)",
        matches="x",
        ticks="outside",
        tickson="labels",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig.update_yaxes(
        autorange="reversed",
        showdividers=False,
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(64,64,64)",
    )
    figs.append(fig)

    ### Per Base Nucleotide Content ###
    df = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["PBSC"]))
    df["stage"] = df["stage"].str.capitalize()

    color_dict = {
        "A": "#DD2000",
        "T": "#36B4F5",
        "C": "#25FE5B",
        "G": "#FBFF00",
    }
    hovertemplate = "<br>".join(["<b>%{y}</b>",
        "Position: <b>%{x} bp</b>",
        "A: <b>%{customdata[5]:.2f}%</b>",
        "T: <b>%{customdata[6]:.2f}%</b>",
        "G: <b>%{customdata[4]:.2f}%</b>",
        "C: <b>%{customdata[7]:.2f}%</b><extra></extra>",
    ])

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        for N in ["A", "T", "G", "C"]:
            # Read 1
            fig.add_trace(
                go.Heatmap(
                    x=df_R1["base"],
                    y=[df_R1["sample_name"], df_R1["stage"]],
                    z=df_R1[N],
                    zmin=0,
                    zmax=100,
                    colorscale=[
                        [0, "rgba(255,255,255,0)"],
                        [1, color_dict[N]],
                    ],
                    showscale=False,
                    hoverinfo="skip" if N == "C" else None,
                    customdata=df_R1 if N == "C" else None,
                    hovertemplate=hovertemplate,
                    hoverongaps=False,
                    ygap=1,
                ),
                row=1,
                col=1,
            )
            # Read 2
            fig.add_trace(
                go.Heatmap(
                    x=df_R2["base"],
                    y=[df_R2["sample_name"], df_R2["stage"]],
                    z=df_R2[N],
                    zmin=0,
                    zmax=100,
                    colorscale=[
                        [0, "rgba(255,255,255,0)"],
                        [1, color_dict[N]],
                    ],
                    showscale=False,
                    hoverinfo="skip" if N == "C" else None,
                    customdata=df_R2 if N == "C" else None,
                    hovertemplate=hovertemplate,
                    hoverongaps=False,
                    ygap=1,
                ),
                row=1,
                col=2,
            )

    # For single-end
    else:
        for N in ["A", "T", "G", "C"]:
            fig = go.Figure()
            fig.add_trace(
                go.Heatmap(
                    x=df["base"],
                    y=[df["sample_name"], df["stage"]],
                    z=df[N],
                    zmin=0,
                    zmax=100,
                    colorscale=[
                        [0, "rgba(255,255,255,0)"],
                        [1, color_dict[N]],
                    ],
                    showscale=False,
                    hoverinfo="skip" if N == "C" else None,
                    customdata=df if N == "C" else None,
                    hovertemplate=hovertemplate,
                    hoverongaps=False,
                    ygap=1,
                )
            )

    for N in ["A", "T", "G", "C"]:
        fig.add_trace(
            go.Bar(
                x=[None],
                y=[None],
                name=N,
                marker_color=color_dict[N],
                marker_line_color="rgb(8,8,8)",
                marker_line_width=0.25,
            )
        )

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>6. Per Base Nucleotide Content</b><br>"
            "<sup>(Source: 03_qc_extras/"
            + settings.QC_FILES["PBSC"]
            + ")</sup>"
        ),
        yaxis=dict(title="Sample - Stage"),
        legend=dict(
            title_text="Nucleotide",
        ),
        height=180 + 30 * len(sample_list),
        plot_bgcolor="rgb(8,8,8)",
    )
    fig.update_xaxes(
        title="Position (bp)",
        linecolor="black",
        ticks="outside",
        matches="x",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig.update_yaxes(
        autorange="reversed",
        ticks="outside",
        dtick=1,
        tickson="labels",
        showdividers=False,
        gridcolor="rgb(64,64,64)",
    )
    figs.append(fig)

    ### Per Read GC Content ###
    df = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["PSGC"]))
    df["stage"] = df["stage"].str.capitalize()
    df_grouped = df.groupby(["sample_name", "read", "stage"], as_index=False)["count"].sum()
    df_merged = pd.merge(df, df_grouped, on=["sample_name", "read", "stage"], how="outer")
    df_merged["freq"] = df_merged["count_x"] / df_merged["count_y"] * 100
    df = df_merged.drop(columns=["count_x", "count_y"])
    df_pivot = df.pivot(
        index=["sample_name", "stage"],
        columns=["read", "gc_content"],
        values="freq",
    ).reset_index()
    df = df_pivot.melt(
        id_vars=["sample_name", "stage"],
        value_name="freq",
    )
    df.sort_values(
        by=["sample_name", "stage"],
        ascending=[True, False],
        inplace=True,
    )
    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig = make_subplots(
            cols=2,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig.add_trace(
            go.Heatmap(
                x=df_R1["gc_content"],
                y=[df_R1["sample_name"], df_R1["stage"]],
                z=df_R1["freq"],
                coloraxis="coloraxis",
                hovertemplate="<b>%{y}</b><br>" +
                              "GC Content: <b>%{x}%</b><br>" +
                              "Proportion: <b>%{z:.2f}%</b><extra></extra>",
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig.add_trace(
            go.Heatmap(
                x=df_R2["gc_content"],
                y=[df_R2["sample_name"], df_R2["stage"]],
                z=df_R2["freq"],
                coloraxis="coloraxis",
                hovertemplate="<b>%{y}</b><br>" +
                              "GC Content: <b>%{x}%</b><br>" +
                              "Proportion: <b>%{z:.2f}%</b><extra></extra>",
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                x=df["gc_content"],
                y=[df["sample_name"], df["stage"]],
                z=df["freq"],
                coloraxis="coloraxis",
                hovertemplate="<b>%{y}</b><br>" +
                              "GC Content: <b>%{x}%</b><br>" +
                              "Proportion: <b>%{z:.2f}%</b><extra></extra>",
                hoverongaps=False,
                ygap=1,
            )
        )

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>7. Per Read GC Content</b><br>"
            "<sup>(Source: 03_qc_extras/"
            + settings.QC_FILES["PSGC"]
            + ")</sup>"
        ),
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
    fig.update_xaxes(
        title="GC Content (%)",
        linecolor="black",
        ticks="outside",
        matches="x",
        gridcolor="rgb(64,64,64)",
        zeroline=False,
    )
    fig.update_yaxes(
        autorange="reversed",
        linecolor="black",
        ticks="outside",
        dtick=1,
        tickson="labels",
        showdividers=False,
        gridcolor="rgb(64,64,64)",
    )
    figs.append(fig)

    ### Sequence Duplication Level ###
    df = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["SDUP"]),
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
        fig = make_subplots(
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
            fig.add_trace(
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
                                  "Duplication Level: <b>%{meta[0]}</b><br>" +
                                  "Percentage: <b>%{x:.2f}%</b><extra></extra>",
                ),
                row=1,
                col=1,
            )
            # Read 2
            data = df_R2[df_R2["duplication_level"] == dup_lev]
            fig.add_trace(
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
                                  "Duplication Level: <b>%{meta[0]}</b><br>" +
                                  "Percentage: <b>%{x:.2f}%</b><extra></extra>",
                ),
                row=1,
                col=2,
            )

    # For single-end
    else:
        fig = go.Figure()
        for dup_lev in dup_lev_list:
            data = df[df["duplication_level"] == dup_lev]
            fig.add_trace(
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
                                  "Duplication Level: <b>%{meta[0]}</b><br>" +
                                  "Percentage: <b>%{x:.2f}%</b><extra></extra>",
                )
            )

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>8. Sequence Duplication Level</b><br>"
            "<sup>(Source: 03_qc_extras/"
            + settings.QC_FILES["SDUP"]
            + ")</sup>"
        ),
        yaxis=dict(title="Sample - Stage"),
        barmode="stack",
        bargap=0,
        bargroupgap=0.05,
        height=180 + 30 * len(sample_list),
        legend=dict(title="Duplication<br>Level", tracegroupgap=0),
        plot_bgcolor="rgb(8,8,8)",
    )
    fig.update_xaxes(
        title="Proportion (%)",
        ticks="outside",
        matches="x",
        range=[0, 100],
        zeroline=False,
        gridcolor="rgb(64,64,64)",
    )
    fig.update_yaxes(
        autorange="reversed",
        showdividers=False,
        ticks="outside",
        dtick=1,
        tickson="labels",
        gridcolor="rgb(64,64,64)",
    )
    figs.append(fig)

    ### Adapter Content ###
    df = pd.read_table(Path(qc_extras_dir, settings.QC_FILES["ADCO"]))
    df["stage"] = df["stage"].str.capitalize()
    df["Total adapter content"] = df.iloc[:,4:].sum(axis=1)
    df_pivot = df.pivot(
        index=["sample_name", "stage"],
        columns=["read", "position"],
    )
    df = df_pivot.stack(
        level=["read", "position"],
        dropna=False,
    ).reset_index()
    df.sort_values(
        by=["sample_name", "stage"],
        ascending=[True, False],
        inplace=True,
    )
    hover_info_list = []
    col_num = 4
    for col in df.columns[4:]:
        hover_info_list.append(col + ": <b>%{customdata[" + str(col_num) + "]:.2f}%</b>")
        col_num += 1
    hover_info_list.insert(0, "<b>%{y}</b><br>Position: <b>%{x} bp</b>")
    hovertemplate = "<br>".join(hover_info_list) + "<extra></extra>"

    # For paired-end
    if "R2" in df["read"].to_list():
        df_R1 = df[df["read"] == "R1"]
        df_R2 = df[df["read"] == "R2"]
        fig = make_subplots(
            cols=2,
            shared_xaxes=True,
            shared_yaxes=True,
            horizontal_spacing=0.02,
            subplot_titles=["Read 1", "Read 2"],
        )
        # Read 1
        fig.add_trace(
            go.Heatmap(
                x=df_R1["position"],
                y=[df_R1["sample_name"], df_R1["stage"]],
                z=df_R1["Total adapter content"],
                coloraxis="coloraxis",
                name="Read 1",
                customdata=df_R1,
                hovertemplate=hovertemplate,
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=1,
        )
        # Read 2
        fig.add_trace(
            go.Heatmap(
                x=df_R2["position"],
                y=[df_R2["sample_name"], df_R2["stage"]],
                z=df_R2["Total adapter content"],
                coloraxis="coloraxis",
                name="Read 2",
                customdata=df_R2,
                hovertemplate=hovertemplate,
                hoverongaps=False,
                ygap=1,
            ),
            row=1,
            col=2,
        )

    # For single-end
    else:
        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                x=df["position"],
                y=[df["sample_name"], df["stage"]],
                z=df["Total adapter content"],
                coloraxis="coloraxis",
                customdata=df,
                hovertemplate=hovertemplate,
                hoverongaps=False,
                ygap=1,
            )
        )

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>9. Adapter Content</b><br>"
            "<sup>(Source: 03_qc_extras/"
            + settings.QC_FILES["ADCO"]
            + ")</sup>"
        ),
        yaxis=dict(title="Sample - Stage"),
        coloraxis=dict(
            colorscale="Spectral_r",
            cmin=0 if max(df["Total adapter content"]) < 10 else None,
            cmax=10 if max(df["Total adapter content"]) < 10 else None,
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
    fig.update_xaxes(
        title="Position (bp)",
        ticks="outside",
        matches="x",
        zeroline=False,
        gridcolor="rgb(64,64,64)",
    )
    fig.update_yaxes(
        autorange="reversed",
        ticks="outside",
        dtick=1,
        tickson="labels",
        showdividers=False,
        gridcolor="rgb(64,64,64)",
    )
    figs.append(fig)

    # Save plots in HTML
    config = dict(
        toImageButtonOptions=dict(
            format="svg",
        ),
        modeBarButtonsToAdd=[
            "v1hovermode",
            "hoverclosest",
            "hovercompare",
            "togglehover",
            "togglespikelines",
            "drawline",
            "drawopenpath",
            "drawclosedpath",
            "drawcircle",
            "drawrect",
            "eraseshape",
        ]
    )
    qc_html_report = Path(out_dir, "captus-assembly_clean.report.html")
    with open(qc_html_report, "w") as f:
        for fig in figs:
            f.write(
                fig.to_html(
                    full_html=False,
                    include_plotlyjs="cdn",
                    config=config
                )
            )
    if qc_html_report.exists() and qc_html_report.is_file():
        qc_html_msg = dim(f"Report generated in {elapsed_time(time.time() - start)}")
    else:
        qc_html_msg = red(f"Report not generated, verify your Python environment")

    return qc_html_report, qc_html_msg


def build_assembly_report(out_dir, asm_stats_tsv):
    start = time.time()

    ### Table ###
    # Load datatable
    df = pd.read_table(asm_stats_tsv)

    if df["filtered_n_contigs"].max() == 0:
        filtered = False
    else:
        filtered = True

    if filtered == True:
        df = df.reindex(
            columns=[
                "sample",
                "total_length",
                "filtered_total_length",
                "n_contigs",
                "filtered_n_contigs",
                "avg_length",
                "N50",
                "longest_contig",
                "shortest_contig",
                "pct_contigs_>=_1kbp",
                "GC_content",
                "filtered_gc_content",
                "avg_depth",
            ],
        )
        df.rename(
            columns={
                "sample": "Sample",
                "total_length": "Total Length Passed (bp)",
                "filtered_total_length": "Total Length Removed (bp)",
                "n_contigs": "Number of Contigs Passed",
                "filtered_n_contigs": "Number of Contigs Removed",
                "avg_length": "Mean Length (bp)",
                "N50": "Contig N50 (bp)",
                "longest_contig": "Longest Contig (bp)",
                "shortest_contig": "Shortest Contig (bp)",
                "pct_contigs_>=_1kbp": "≥1 kbp Contigs (%)",
                "GC_content": "GC Content Passed (%)",
                "filtered_gc_content": "GC Content Removed (%)",
                "avg_depth": "Mean Depth (x)",
            },
            inplace=True,
        )
        format = [
            None,
            ",",
            ",",
            ",",
            ",",
            ",",
            ",",
            ",",
            ",",
            ".2f",
        ]
    else:
        df = df.reindex(
            columns=[
                "sample",
                "total_length",
                "n_contigs",
                "avg_length",
                "N50",
                "longest_contig",
                "shortest_contig",
                "pct_contigs_>=_1kbp",
                "GC_content",
                "avg_depth",
            ],
        )
        df.rename(
            columns={
                "sample": "Sample",
                "total_length": "Total Length (bp)",
                "n_contigs": "Number of Contigs",
                "avg_length": "Mean Length (bp)",
                "N50": "Contig N50 (bp)",
                "longest_contig": "Longest Contig (bp)",
                "shortest_contig": "Shortest Contig (bp)",
                "pct_contigs_>=_1kbp": "≥1 kbp Contigs (%)",
                "GC_content": "GC Content(%)",
                "avg_depth": "Mean Depth (x)",
            },
            inplace=True,
        )
        format = [
            None,
            ",",
            ",",
            ",",
            ",",
            ",",
            ",",
            ".2f",
        ]

    sample_list = df["Sample"].unique()

    colorscale = [
        [0, "#F0D9E6"],
        [0.5, "#F5F5F5"],
        [1, "#BDE3D8"]
    ]

    fill_color = []
    for col in df.columns:
        if col == "Sample":
            fill_color.append("#F5F5F5")
        else:
            fill_color.append(
                sample_colorscale(
                    colorscale,
                    normalize(df[col])
                )
            )

    figs = []
    fig = go.Figure()
    fig.add_trace(
        go.Table(
            header=dict(values=["<b>" + col + "</b>" for col in df.columns]),
            cells=dict(
                values=[df[col] for col in df.columns],
                fill_color=fill_color,
                format=format,
                align=["left", "right"],
                height=21,
            )
        )
    )

    buttons = []
    for col in df.columns:
        df.sort_values(
            by=col,
            ascending=True if col == "Sample" else False,
            inplace=True,
        )
        fill_color = []
        for col2 in df.columns:
            if col2 == "Sample":
                fill_color.append("#F5F5F5")
            else:
                fill_color.append(
                    sample_colorscale(
                        colorscale,
                        normalize(df[col2])
                    )
                )
        button = dict(
            label=col,
            method="restyle",
            args=[
                dict(
                    cells=dict(
                        values=[df[col] for col in df.columns],
                        fill=dict(
                            color=fill_color,
                        ),
                        format=format,
                        align=["left", "right"],
                        height=21,
                    )
                )
            ],
        )
        buttons.append(button)

    updatemenus = [dict(
        buttons=buttons,
        type="dropdown",
        direction="down",
        pad=dict(t=10, b=10),
        x=1,
        xanchor="right",
        y=1,
        yanchor="bottom",
    )]

    annotations=[dict(
        text="<b>Sort by:</b>",
        x=1,
        xref="paper",
        xanchor="right",
        xshift=-200 if filtered == True else -155,
        y=1,
        yref="paper",
        yanchor="top",
        yshift=36,
        align="right",
        showarrow=False
    )]

    fig.update_layout(
        font_family="Arial",
        title_text=(
            "<b>Captus-assembly: Assemble (<i>De Novo</i> Assembly Report)<br>"
            "1. Summary Table</b><br>"
            "<sup>(Source: "
            + str(asm_stats_tsv.name)
            + ")</sup>"
        ),
        height=230 + 21 * len(sample_list) if len(sample_list) < 31 else None,
        updatemenus=updatemenus,
        annotations=annotations,
    )
    figs.append(fig)

    ### Bar plot ###
    # Load datatable
    df = pd.read_table(asm_stats_tsv)

    # Variables available as drop-down menu
    var_dict = {
        "Total Length (bp)": [
            "total_length",
            None,
            None,
            "filtered_total_length" if filtered == True else None,
        ],
        "Number of Contigs": [
            "n_contigs",
            None,
            None,
            "filtered_n_contigs" if filtered == True else None,
        ],
        "Mean Length (bp)": [
            "avg_length",
            None,
            None,
            "filtered_avg_length" if filtered == True else None,
        ],
        "Median Length (bp)": [
            "median_length",
            None,
            None,
            None,
        ],
        "Contig N50 (bp)": [
            "N50",
            None,
            None,
            None,
        ],
        "Longest Contig (bp)": [
            "longest_contig",
            None,
            None,
            None,
        ],
        "Shortest Contig (bp)": [
            "shortest_contig",
            None,
            None,
            None,
        ],
        "Contig Breakdown<br> by Length (%)": [
            "pct_contigs_>=_1kbp",
            "pct_contigs_>=_2kbp",
            "pct_contigs_>=_5kbp",
            "pct_contigs_>=_10kbp",
        ],
        "Length Breakdown<br> by Contig Length (%)": [
            "pct_length_>=_1kbp",
            "pct_length_>=_2kbp",
            "pct_length_>=_5kbp",
            "pct_length_>=_10kbp",
        ],
        "GC Content (%)": [
            "GC_content",
            None,
            None,
            "filtered_gc_content" if filtered == True else None,
        ],
        "Mean Depth (x)": [
            "avg_depth",
            None,
            None,
            None,
        ],
        "Contig Breakdown<br> by Depth (%)": [
            "pct_contigs_>=_1x",
            "pct_contigs_>=_2x",
            "pct_contigs_>=_5x",
            "pct_contigs_>=_10x",
        ],
    }

    hover_template_dict = {
        "Total Length (bp)": "<br>".join(
            [
                "Sample: <b>%{x}</b>",
                "Total length: <b>%{y:,} bp</b>",
                "<extra></extra>" if filtered == False else "",
            ]
        ),
        "Number of Contigs": "<br>".join(
            [
                "Sample: <b>%{x}</b>",
                "Number of contigs: <b>%{y:,}</b>",
                "<extra></extra>" if filtered == False else "",
            ]
        ),
        "Mean Length (bp)": "<br>".join(
            [
                "Sample: <b>%{x}</b>",
                "Mean contig length: <b>%{y:,} bp</b>",
                "<extra></extra>" if filtered == False else "",
            ]
        ),
        "Median Length (bp)": (
            "Sample: <b>%{x}</b><br>" +
            "Median contig length: <b>%{y:,} bp</b><extra></extra>"
        ),
        "Contig N50 (bp)": (
            "Sample: <b>%{x}</b><br>" +
            "Contig N50: <b>%{y:,} bp</b><extra></extra>"
        ),
        "Longest Contig (bp)": (
            "Sample: <b>%{x}</b><br>" +
            "Longest contig length: <b>%{y:,} bp</b><extra></extra>"
        ),
        "Shortest Contig (bp)": (
            "Sample: <b>%{x}</b><br>" +
            "Shortest contig length: <b>%{y:,} bp</b><extra></extra>"
        ),
        "Contig Breakdown<br> by Length (%)": (
            "Sample: <b>%{x}</b><br>" +
            "Proportion of contigs: <b>%{y:.2f}%</b>"
        ),
        "Length Breakdown<br> by Contig Length (%)": (
            "Sample: <b>%{x}</b><br>" +
            "Proportion of length: <b>%{y:.2f}%</b>"
        ),
        "GC Content (%)": "<br>".join(
            [
                "Sample: <b>%{x}</b>",
                "GC content: <b>%{y:.2f}%</b>",
                "<extra></extra>" if filtered == False else "",
            ]
        ),
        "Mean Depth (x)": (
            "Sample: <b>%{x}</b><br>" +
            "Mean depth: <b>%{y:.2f}x</b><extra></extra>"
        ),
        "Contig Breakdown<br> by Depth (%)": (
            "Sample: <b>%{x}</b><br>" +
            "Proportion of contigs: <b>%{y:.2f}%</b>"
        ),
    }

    trace_name_dict = {
        "n_contigs": "Passed",
        "filtered_n_contigs": "Removed",
        "total_length": "Passed",
        "filtered_total_length": "Removed",
        "avg_length": "Passed",
        "filtered_avg_length": "Removed",
        "GC_content": "Passed",
        "filtered_gc_content": "Removed",
    }

    colors = [
        "#56B4E9",
        "#009E73",
        "#E69F00",
        "#CC79A7"
    ]

    # Create figure
    fig = go.Figure()
    for i in range(len(list(var_dict.values())[0])):
        if var_dict[list(var_dict.keys())[0]][i] is None:
            fig.add_trace(
                go.Bar(
                    marker_color=colors[i],
                    marker_line_color="rgb(8,8,8)",
                )
            )
        else:
            fig.add_trace(
                go.Bar(
                    x=df["sample"],
                    y=df[var_dict[list(var_dict.keys())[0]][i]],
                    name=trace_name_dict[var_dict[list(var_dict.keys())[0]][i]],
                    marker_color=colors[i],
                    marker_line_color="rgb(8,8,8)",
                    hovertemplate=hover_template_dict[list(var_dict.keys())[0]],
                )
            )

    # Dropdown setting
    buttons = []
    for key, values in var_dict.items():
        x_list, y_list, name_list, hovertemplate_list = [], [], [], []
        for v in values:
            if v is None:
                x_list.append(None)
                y_list.append(None)
                name_list.append(None)
                hovertemplate_list.append(None)
            else:
                x_list.append(df["sample"])
                y_list.append(df[v])
                if v in list(trace_name_dict.keys()):
                    name_list.append(trace_name_dict[v])
                else:
                    name_list.append(re.sub(".*_>=_(\d+)", r"≥\1 ", v))
                hovertemplate_list.append(hover_template_dict[key])
        if key in [
            "Contig Breakdown<br> by Length (%)",
            "Length Breakdown<br> by Contig Length (%)",
            "Contig Breakdown<br> by Depth (%)"
        ]:
            barmode = "overlay"
        elif key in [
            "Mean Length (bp)",
            "GC Content (%)",
        ]:
            barmode = "group"
        else:
            barmode = "stack"
        button = dict(
            label=key,
            method="update",
            args=[
                dict(
                    x=x_list,
                    y=y_list,
                    name=name_list,
                    hovertemplate=hovertemplate_list,
                ),
                dict(
                    barmode=barmode,
                )
            ]
        )
        buttons.append(button)

    updatemenus = [
        dict(
            buttons=buttons,
            type="dropdown",
            direction="down",
            pad={"t": 10, "b": 10, "r": 50},
            showactive=True,
            x=0,
            xanchor="right",
            y=0.5,
            yanchor="middle"
        )
    ]
    # Layout setting
    fig.update_layout(
        plot_bgcolor="rgb(8,8,8)",
        font_family="Arial",
        title_text=(
            "<b>2. Visual Stats</b><br>"
            "<sup>(Source: "
            + str(asm_stats_tsv.name)
            + ")</sup>"
        ),
        yaxis=dict(
            gridcolor="rgb(64,64,64)",
            ticks="outside",
            zeroline=False,
        ),
        xaxis=dict(
            title="Sample",
            showgrid=False,
            ticks="outside",
        ),
        hoverlabel=dict(
            font_color="rgb(64,64,64)",
            bordercolor="rgb(64,64,64)",
        ),
        barmode="stack",
        legend=dict(
            traceorder="normal",
        ),
        # height=180 + 15 * len(sample_list),
        updatemenus=updatemenus,
    )
    figs.append(fig)

    # Save plot in HTML
    config = dict(
        toImageButtonOptions=dict(
            format="svg",
        ),
        modeBarButtonsToAdd=[
            "v1hovermode",
            "hoverclosest",
            "hovercompare",
            "togglehover",
            "togglespikelines",
            "drawline",
            "drawopenpath",
            "drawclosedpath",
            "drawcircle",
            "drawrect",
            "eraseshape",
        ]
    )
    asm_html_report = Path(out_dir, "captus-assembly_assemble.report.html")
    with open(asm_html_report, "w") as f:
        for fig in figs:
            f.write(
                fig.to_html(
                    full_html=False,
                    include_plotlyjs="cdn",
                    config=config
                )
            )
    if asm_html_report.exists() and asm_html_report.is_file():
        asm_html_msg = dim(f"Report generated in {elapsed_time(time.time() - start)}")
    else:
        asm_html_msg = red(f"Report not generated, verify your Python environment")

    return asm_html_report, asm_html_msg


def build_extraction_report(out_dir, ext_stats_tsv):
    start = time.time()

    # Load datatable
    df = pd.read_table(
        ext_stats_tsv,
        low_memory=False,
        usecols=[*range(0,22)],
    )
    df.sort_values(
        by=["sample_name", "marker_type", "locus"],
        inplace=True,
    )
    # Preprocess
    df_best = df[df["hit"] == 0].reset_index(drop=True).fillna("NA")
    df_best["hit"] = df.groupby(["sample_name", "marker_type", "locus"], as_index=False).count()["hit"]
    df_best["ref_len_unit"] = np.where(df_best["ref_type"] == "prot", "aa", "bp")
    df_best.loc[df_best["frameshifts"] == "NA", "n_frameshifts"] = 0
    df_best.loc[df_best["frameshifts"] != "NA", "n_frameshifts"] = df_best["frameshifts"].str.count(",") + 1
    # Define variables
    marker_type = df_best["marker_type"].unique()
    if len(marker_type) > 1:
        marker_type = np.insert(marker_type, 0, "ALL")
    var_list = [
        "pct_recovered",
        "pct_identity",
        "hit",
        "score",
        "wscore",
        "n_frameshifts",
        "hit_contigs",
        "hit_l50",
        "hit_l90",
        "hit_lg50",
        "hit_lg90",
    ]
    var_lab_list = [
        "Recovered Length (%)",
        "Identity (%)",
        "Total Hits (Copies)",
        "Score",
        "Weighted Score",
        "Number of Frameshifts",
        "Contigs in Best Hit",
        "Best Hit L50",
        "Best Hit L90",
        "Best Hit LG50",
        "Best Hit LG90",
    ]
    colorscale = [
        [0.0, "rgb(94,79,162)"],
        [0.1, "rgb(50,136,189)"],
        [0.2, "rgb(102,194,165)"],
        [0.3, "rgb(171,221,164)"],
        [0.4, "rgb(230,245,152)"],
        [0.5, "rgb(255,255,191)"],
        [0.6, "rgb(254,224,139)"],
        [0.7, "rgb(253,174,97)"],
        [0.8, "rgb(244,109,67)"],
        [0.9, "rgb(213,62,79)"],
        [1.0, "rgb(158,1,66)"],
        [1.0, "rgb(255, 64, 255)"],
        [1.0, "rgb(255, 64, 255)"],
    ]
    colorscale2 = [
        [0.0, "rgb(94,79,162)"],
        [0.1, "rgb(50,136,189)"],
        [0.2, "rgb(102,194,165)"],
        [0.3, "rgb(171,221,164)"],
        [0.4, "rgb(230,245,152)"],
        [0.5, "rgb(255,255,191)"],
        [0.6, "rgb(254,224,139)"],
        [0.7, "rgb(253,174,97)"],
        [0.8, "rgb(244,109,67)"],
        [0.9, "rgb(213,62,79)"],
        [1.0, "rgb(158,1,66)"],
    ]
    # Dropdown for sorting
    buttons2 = [
        dict(
            label="None",
            method="relayout",
            args=[{
                "xaxis.categoryorder": "category ascending",
                "yaxis.categoryorder": "category descending",
            }],
        ),
        dict(
            label="Mean X",
            method="relayout",
            args=[{
                "xaxis.categoryorder": "mean descending",
                "yaxis.categoryorder": "category descending",
            }],
        ),
        dict(
            label="Mean Y",
            method="relayout",
            args=[{
                "xaxis.categoryorder": "category ascending",
                "yaxis.categoryorder": "mean ascending",
            }],
        ),
        dict(
            label="Mean Both",
            method="relayout",
            args=[{
                "xaxis.categoryorder": "mean descending",
                "yaxis.categoryorder": "mean ascending",
            }],
        ),
        dict(
            label="Total X",
            method="relayout",
            args=[{
                "xaxis.categoryorder": "total descending",
                "yaxis.categoryorder": "category descending",
            }],
        ),
        dict(
            label="Total Y",
            method="relayout",
            args=[{
                "xaxis.categoryorder": "category ascending",
                "yaxis.categoryorder": "total ascending",
            }],
        ),
        dict(
            label="Total Both",
            method="relayout",
            args=[{
                "xaxis.categoryorder": "total descending",
                "yaxis.categoryorder": "total ascending",
            }],
        ),
    ]

    # Make heatmap for each marker type
    figs = []
    for i, marker in enumerate(marker_type):
        if marker == "ALL":
            data = df_best
            data["marker_type - locus"] = data["marker_type"] + " - " + data["locus"].astype(str)
            matrix_size = len(data["sample_name"].unique()) * len(data["marker_type - locus"].unique())
            bar1 = data.groupby("marker_type - locus", as_index=False)["sample_name"].count()
            bar2 = data.groupby("sample_name", as_index=False)["marker_type - locus"].count()
        else:
            data = df_best[df_best["marker_type"] == marker]
            matrix_size = len(data["sample_name"].unique()) * len(data["locus"].unique())
            bar1 = data.groupby("locus", as_index=False)["sample_name"].count()
            bar2 = data.groupby("sample_name", as_index=False)["locus"].count()
        if matrix_size > 500000:
            customdata = None
            hovertemplate = "<br>".join([
                "Sample: <b>%{y}</b>",
                "Locus: <b>%{x}</b>",
                "Recovered length: <b>%{z:.2f}%</b><extra></extra>",
            ])
        else:
            customdata = data
            hovertemplate = "<br>".join([
                "Sample: <b>%{customdata[0]}</b>",
                "Marker type: <b>%{customdata[1]}</b>",
                "Locus: <b>%{customdata[2]}</b>",
                "Ref name: <b>%{customdata[3]}</b>",
                "Ref coords: <b>%{customdata[4]}</b>",
                "Ref type: <b>%{customdata[5]}</b>",
                "Ref len matched: <b>%{customdata[6]:,.0f} %{customdata[22]}</b>",
                "Total hits (copies): <b>%{customdata[7]}</b>",
                "Recovered length: <b>%{customdata[8]:.2f}%</b>",
                "Identity: <b>%{customdata[9]:.2f}%</b>",
                "Score: <b>%{customdata[10]:.3f}</b>",
                "Weighted score: <b>%{customdata[11]:.3f}</b>",
                "Hit length: <b>%{customdata[12]:,.0f} bp</b>",
                "CDS length: <b>%{customdata[13]:,.0f} bp</b>",
                "Intron length: <b>%{customdata[14]:,.0f} bp</b>",
                "Flanking length: <b>%{customdata[15]:,.0f} bp</b>",
                "Number of frameshifts: <b>%{customdata[23]}</b>",
                "Position of frameshifts: <b>%{customdata[16]}</b>",
                "Contigs in best hit: <b>%{customdata[17]}</b>",
                "Best hit L50: <b>%{customdata[18]}</b>",
                "Best hit L90: <b>%{customdata[19]}</b>",
                "Best hit LG50: <b>%{customdata[20]}</b>",
                "Best hit LG90: <b>%{customdata[21]}</b><extra></extra>",
            ])
        hovertemplate_bar1 = "<br>".join([
            "Locus: <b>%{x}</b>",
            "Samples recovered: <b>%{y}</b><extra></extra>",
        ])
        hovertemplate_bar2 = "<br>".join([
            "Sample: <b>%{y}</b>",
            "Loci recovered: <b>%{x}</b><extra></extra>",
        ])

        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                x=data["marker_type - locus"] if marker == "ALL" else data["locus"],
                y=data["sample_name"],
                z=data[var_list[0]],
                zmin=data[var_list[0]].min(),
                zmax=data[var_list[0]].max() if data[var_list[0]].max() < 200 else 200,
                colorscale=colorscale2 if data[var_list[0]].max() < 200 else colorscale,
                colorbar=dict(
                    ticks="outside",
                    ticksuffix="%",
                    outlinecolor="rgb(8,8,8)",
                    outlinewidth=1,
                ),
                customdata=customdata,
                hovertemplate=hovertemplate,
                hoverongaps=False,
                xgap=0.5,
                ygap=0.5,
            )
        )
        fig.add_trace(
            go.Bar(
                x=bar1["marker_type - locus"] if marker == "ALL" else bar1["locus"],
                y=bar1["sample_name"],
                yaxis="y2",
                hovertemplate=hovertemplate_bar1,
                marker=dict(
                    color="#56B4E9",
                    line=dict(
                        color="rgb(8,8,8)",
                        width=0.5,
                    ),
                ),
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Bar(
                x=bar2["marker_type - locus"] if marker == "ALL" else bar2["locus"],
                y=bar2["sample_name"],
                xaxis="x2",
                hovertemplate=hovertemplate_bar2,
                marker=dict(
                    color="#56B4E9",
                    line=dict(
                        color="rgb(8,8,8)",
                        width=0.5,
                    ),
                ),
                orientation="h",
                showlegend=False,
            )
        )

        # Dropdown for variables
        buttons1 = []
        for j, var in enumerate(var_list):
            if var == "pct_recovered":
                zmax = data[var].max() if data[var].max() < 200 else 200
                cmap = [colorscale2] if data[var].max() < 200 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Recovered length: <b>%{z:.2f}%</b><extra></extra>",
                    ])
            elif var == "pct_identity":
                zmax = data[var].max() if data[var].max() < 100 else 100
                cmap = [colorscale2] if data[var].max() <= 100 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Identity: <b>%{z:.2f}%</b><extra></extra>",
                    ])
            elif var == "hit":
                zmax = data[var].max() if data[var].max() < 50 else 50
                cmap = [colorscale2] if data[var].max() < 50 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Total hits (copies): <b>%{z}</b><extra></extra>",
                    ])
            elif var == "score":
                zmax = data[var].max() if data[var].max() < 2 else 2
                cmap = [colorscale2] if data[var].max() < 2 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Score: <b>%{z:.3f}</b><extra></extra>",
                    ])
            elif var == "wscore":
                zmax = data[var].max() if data[var].max() < 2 else 2
                cmap = [colorscale2] if data[var].max() < 2 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Weighted score: <b>%{z:.3f}</b><extra></extra>",
                    ])
            elif var == "n_frameshifts":
                zmax = data[var].max() if data[var].max() < 10 else 10
                cmap = [colorscale2] if data[var].max() < 10 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Number of frameshifts: <b>%{z}</b><extra></extra>",
                    ])
            elif var == "hit_contigs":
                zmax = data[var].max() if data[var].max() < 10 else 10
                cmap = [colorscale2] if data[var].max() < 10 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Contigs in best hit: <b>%{z}</b><extra></extra>",
                    ])
            elif var == "hit_l50":
                zmax = data[var].max() if data[var].max() < 10 else 10
                cmap = [colorscale2] if data[var].max() < 10 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Best hit L50: <b>%{z}</b><extra></extra>",
                    ])
            elif var == "hit_l90":
                zmax = data[var].max() if data[var].max() < 10 else 10
                cmap = [colorscale2] if data[var].max() < 10 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Best hit L90: <b>%{z}</b><extra></extra>",
                    ])
            elif var == "hit_lg50":
                zmax = data[var].max() if data[var].max() < 10 else 10
                cmap = [colorscale2] if data[var].max() < 10 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Best hit LG50: <b>%{z}</b><extra></extra>",
                    ])
            elif var == "hit_lg90":
                zmax = data[var].max() if data[var].max() < 10 else 10
                cmap = [colorscale2] if data[var].max() < 10 else [colorscale]
                if matrix_size > 500000:
                    hovertemplate = "<br>".join([
                        "Sample: <b>%{y}</b>",
                        "Locus: <b>%{x}</b>",
                        "Best hit LG90: <b>%{z}</b><extra></extra>",
                    ])
            button = dict(
                label=var_lab_list[j],
                method="restyle",
                args=[
                    {
                        "z": [data[var]],
                        "zmin": data[var].min(),
                        "zmax": zmax,
                        "colorscale": cmap,
                        "colorbar.ticksuffix": None if j > 1 else "%",
                        "hovertemplate": [
                            hovertemplate,
                            hovertemplate_bar1,
                            hovertemplate_bar2,
                        ],
                    }
                ],
            )
            buttons1.append(button)

        updatemenus=[
            dict(
                buttons=buttons1,
                type="dropdown",
                direction="down",
                pad={"t": 10, "b": 10},
                showactive=True,
                x=0.475,
                xanchor="right",
                y=1,
                yanchor="bottom",
            ),
            dict(
                buttons=buttons2,
                type="dropdown",
                direction="down",
                pad={"t": 10, "b": 10},
                showactive=True,
                x=0.95,
                xanchor="right",
                y=1,
                yanchor="bottom",
            ),
        ]

        annotations = [
            dict(
                text="<b>Variable:</b>",
                x=0.475,
                xref="paper",
                xanchor="right",
                xshift=-165,
                y=1,
                yref="paper",
                yanchor="top",
                yshift=36,
                align="right",
                showarrow=False,
            ),
            dict(
                text="<b>Sort by Value:</b>",
                x=0.95,
                xref="paper",
                xanchor="right",
                xshift=-102,
                y=1,
                yref="paper",
                yanchor="top",
                yshift=36,
                align="right",
                showarrow=False,
            ),
        ]

        # Layout setting
        if i == 0:
            if len(marker_type) == 1:
                title = (
                        "<b>Captus-assembly: Extract (Marker Recovery Report)</b><br>"
                        "<sup>(Source: "
                        + str(ext_stats_tsv.name)
                        + ")</sup>"
                )
            else:
                title = (
                        "<b>Captus-assembly: Extract (Marker Recovery Report)<br>"
                        "1. Marker Type: "
                        + marker
                        + "</b><br><sup>(Source: "
                        + str(ext_stats_tsv.name)
                        + ")</sup>"
                )
        else:
            title = (
                "<b>"
                + str(i + 1)
                + ". Marker Type: "
                + marker
                + "</b><br><sup>(Source: "
                + str(ext_stats_tsv.name)
                + ")</sup>"
            )
        fig.update_layout(
            font_family="Arial",
            plot_bgcolor="rgb(8,8,8)",
            title=title,
            xaxis1=dict(
                title="Marker type - Locus" if marker == "ALL" else "Locus",
                type="category",
                categoryorder="category ascending",
                gridcolor="rgb(64,64,64)",
                ticks="outside",
                domain=[0, 0.95],
            ),
            yaxis1=dict(
                title="Sample",
                type="category",
                categoryorder="category descending",
                gridcolor="rgb(64,64,64)",
                ticks="outside",
                domain=[0, 0.925],
            ),
            xaxis2=dict(
                title="Loci",
                gridcolor="rgb(64,64,64)",
                ticks="outside",
                zeroline=True,
                domain=[0.95, 1],
            ),
            yaxis2=dict(
                title="Samples",
                gridcolor="rgb(64,64,64)",
                ticks="outside",
                zeroline=True,
                domain=[0.925, 1],
            ),
            hoverlabel=dict(
                font_color="#FFFFFF",
                bordercolor="#FFFFFF",
                bgcolor="#444444",
            ),
            bargap=0,
            annotations=annotations,
            updatemenus=updatemenus
        )
        figs.append(fig)

    # Save plot in HTML
    config = dict(
        scrollZoom=True,
        toImageButtonOptions=dict(
            format="svg",
        ),
        modeBarButtonsToAdd=[
            "v1hovermode",
            "hoverclosest",
            "hovercompare",
            "togglehover",
            "togglespikelines",
            "drawline",
            "drawopenpath",
            "drawclosedpath",
            "drawcircle",
            "drawrect",
            "eraseshape",
        ]
    )
    ext_html_report = Path(out_dir, "captus-assembly_extract.report.html")
    with open(ext_html_report, "w") as f:
        for fig in figs:
            f.write(
                fig.to_html(
                    full_html=False,
                    include_plotlyjs="cdn",
                    config=config
                )
            )
    if ext_html_report.exists() and ext_html_report.is_file():
        ext_html_msg = dim(f"Report generated in {elapsed_time(time.time() - start)}")
    else:
        ext_html_msg = red(f"Report not generated, verify your Python environment")

    return ext_html_report, ext_html_msg


def build_alignment_report(out_dir, aln_stats_tsv, sam_stats_tsv):
    start = time.time()

    df = pd.read_table(aln_stats_tsv)
    df["stage"] = df["trimmed"].astype("str").str.cat([df["paralog_filter"], df["with_refs"].astype("str")], sep="-")
    marker_type_list = df["marker_type"].unique()

    color_dict = {
        "AA": "#CC79A7",
        "NT": "#56B4E9",
        "GE": "#E69F00",
        "GF": "#009E73",
        "MA": "#CC79A7",
        "MF": "#009E73",
    }

    stage_dict = {
        "False-unfiltered-True": "02_untrimmed/ <br>01_unfiltered_w_refs -",
        "False-naive-True": "02_untrimmed/ <br>02_naive_w_refs -",
        "False-informed-True": "02_untrimmed/ <br>03_informed_w_refs -",
        "False-unfiltered-False": "02_untrimmed/ <br>04_unfiltered -",
        "False-naive-False": "02_untrimmed/ <br>05_naive -",
        "False-informed-False": "02_untrimmed/ <br>06_informed -",
        "True-unfiltered-True": "03_trimmed/ <br>01_unfiltered_w_refs -",
        "True-naive-True": "03_trimmed/ <br>02_naive_w_refs -",
        "True-informed-True": "03_trimmed/ <br>03_informed_w_refs -",
        "True-unfiltered-False": "03_trimmed/ <br>04_unfiltered -",
        "True-naive-False": "03_trimmed/ <br>05_naive -",
        "True-informed-False": "03_trimmed/ <br>06_informed -",
    }

    var_dict = {
        "Sequences": "seqs",
        "Samples": "samples",
        "Sequences Per Sample": "avg_copies",
        "Alignment Length": "sites",
        "Informative Sites": "informative",
        "Informativeness (%)": "informativeness",
        "Uninformative Sites": "uninformative",
        "Constant Sites": "constant",
        "Singleton Sites": "singleton",
        "Patterns": "patterns",
        "Mean Pairwise Identity (%)": "avg_pid",
        "Missingness (%)": "missingness",
        "GC Content (%)": "gc",
        "GC Content at<br>1st Codon Position (%)": "gc_codon_p1",
        "GC Content at<br>2nd Codon Position (%)": "gc_codon_p2",
        "GC Content at<br>3rd Codon Position (%)": "gc_codon_p3",
    }

    figs = []

    fig = go.Figure()
    for j, marker_type in enumerate(marker_type_list):
        data = df[df["marker_type"] == marker_type]
        fmt_list = list(reversed(data["format"].unique()))
        for i, fmt in enumerate(fmt_list):
            fig.add_trace(
                go.Violin(
                    x=data[data["format"] == fmt][list(var_dict.values())[3]],
                    y=data[data["format"] == fmt]["stage"].map(stage_dict),
                    name=fmt,
                    line=dict(
                        color=color_dict[fmt],
                    ),
                    spanmode="hard",
                    meanline_visible=True,
                    hoverinfo="x",
                    visible=True if j == 0 else False,
                    orientation="h",
                    side="negative",
                    width=1.9,
                    points=False,
                )
            )
        if j == 0:
            false_list = [[False] * len(fmt_list)]
        else:
            false_list += [[False] * len(fmt_list)]

    buttons = []
    for j, marker_type in enumerate(marker_type_list):
        data = df[df["marker_type"] == marker_type]
        fmt_list = data["format"].unique()
        visible = list(false_list)
        visible[j] = [True] * len(fmt_list)
        visible = sum(visible, [])
        button = dict(
            label=marker_type,
            method="restyle",
            args=[
                dict(
                    visible=visible,
                )
            ]
        )
        buttons.append(button)

    buttonsX = []
    for lab, var in var_dict.items():
        x_list, hoverinfo_list = [], []
        for j, marker_type in enumerate(marker_type_list):
            data = df[df["marker_type"] == marker_type]
            fmt_list = list(reversed(data["format"].unique()))
            for i, fmt in enumerate(fmt_list):
                if var == "gc" and fmt == "AA":
                    hoverinfo = "skip"
                elif "codon" in var and fmt in ["AA", "GE", "GF", "MA", "MF"]:
                    hoverinfo = "skip"
                else:
                    hoverinfo = "x"
                x = data[data["format"] == fmt][var]
                x_list.append(x)
                hoverinfo_list.append(hoverinfo)
        buttonX = dict(
            label=lab,
            method="restyle",
            args=[
                dict(
                    x=x_list,
                    hoverinfo=hoverinfo_list,
                ),
            ]
        )
        buttonsX.append(buttonX)

    updatemenus = [
        dict(
            buttons=buttonsX,
            type="dropdown",
            direction="up",
            pad={"t": 30, "b": 10},
            active=3,
            showactive=True,
            x=0.5,
            xanchor="center",
            y=0,
            yanchor="top"
        ),
    ]
    if len(marker_type_list) > 1:
        updatemenus.append(
            dict(
                buttons=buttons,
                type="dropdown",
                direction="down",
                pad={"t": 10, "b": 10},
                showactive=True,
                x=1,
                xanchor="right",
                y=1,
                yanchor="bottom"
            )
        )

        annotations = [
            dict(
                text="<b>Marker Type:</b>",
                x=1,
                xref="paper",
                xanchor="right",
                xshift=-75,
                y=1,
                yref="paper",
                yanchor="top",
                yshift=36,
                align="right",
                showarrow=False,
            ),
        ]

    fig.update_layout(
        font_family="Arial",
        plot_bgcolor="rgb(8,8,8)",
        title=(
            "<b>Captus-assembly: Align (Alignment Report)<br>"
            + "1. Stats Comparison at Each Processing Step"
            + "</b><br><sup>(Source: "
            + str(aln_stats_tsv.name)
            + ")</sup>"
        ),
        xaxis=dict(
            showgrid=True,
            gridcolor="rgb(64,64,64)",
            ticks="outside",
            rangemode="nonnegative",
            zeroline=False,
        ),
        yaxis=dict(
            title="Processing Step",
            showgrid=False,
            # ticks="outside",
            ticklabelposition="outside top",
            zeroline=False,
            autorange="reversed",
        ),
        hoverlabel=dict(
                font_color="rgb(64,64,64)",
                bordercolor="rgb(64,64,64)",
            ),
        legend=dict(
            title=dict(
                text="<b>Format</b>",
                side="top",
            ),
            traceorder="reversed+grouped",
        ),
        annotations=annotations if len(marker_type_list) > 1 else None,
        updatemenus=updatemenus,
    )
    figs.append(fig)

    stage_dict = {
        "True-informed-False":    "03_trimmed/<br>06_informed",
        "True-naive-False":       "03_trimmed/<br>05_naive",
        "True-unfiltered-False":  "03_trimmed/<br>04_unfiltered",
        "True-informed-True":     "03_trimmed/<br>03_informed_w_refs",
        "True-naive-True":        "03_trimmed/<br>02_naive_w_refs",
        "True-unfiltered-True":   "03_trimmed/<br>01_unfiltered_w_refs",
        "False-informed-False":   "02_untrimmed/<br>06_informed",
        "False-naive-False":      "02_untrimmed/<br>05_naive",
        "False-unfiltered-False": "02_untrimmed/<br>04_unfiltered",
        "False-informed-True":    "02_untrimmed/<br>03_informed_w_refs",
        "False-naive-True":       "02_untrimmed/<br>02_naive_w_refs",
        "False-unfiltered-True":  "02_untrimmed/<br>01_unfiltered_w_refs",
    }

    hovertemplate_aa = "<br>".join([
            "Locus: <b>%{customdata[6]}</b>",
            "Marker type: <b>%{customdata[4]}</b>",
            "Sequences: <b>%{customdata[7]:,.0f}</b>",
            "Samples: <b>%{customdata[8]:,.0f}</b>",
            "Sequences per sample: <b>%{customdata[9]:.2f}</b>",
            "Alignment length: <b>%{customdata[10]:,.0f} aa</b>",
            "Informative sites: <b>%{customdata[11]:,.0f}</b>",
            "Informativeness: <b>%{customdata[12]:.2f}%</b>",
            "Uninformative sites: <b>%{customdata[13]:,.0f}</b>",
            "Constant sites: <b>%{customdata[14]:,.0f}</b>",
            "Singleton sites: <b>%{customdata[15]:,.0f}</b>",
            "Patterns: <b>%{customdata[16]:,.0f}</b>",
            "Mean pairwise identity: <b>%{customdata[17]:.2f}%</b>",
            "Missingness: <b>%{customdata[18]:.2f}%</b>",
    ])
    hovertemplate_nt = "<br>".join([
            "Locus: <b>%{customdata[6]}</b>",
            "Marker type: <b>%{customdata[4]}</b>",
            "Sequences: <b>%{customdata[7]:,.0f}</b>",
            "Samples: <b>%{customdata[8]:,.0f}</b>",
            "Sequences per sample: <b>%{customdata[9]:.2f}</b>",
            "Alignment length: <b>%{customdata[10]:,.0f} bp</b>",
            "Informative sites: <b>%{customdata[11]:,.0f}</b>",
            "Informativeness: <b>%{customdata[12]:.2f}%</b>",
            "Uninformative sites: <b>%{customdata[13]:,.0f}</b>",
            "Constant sites: <b>%{customdata[14]:,.0f}</b>",
            "Singleton sites: <b>%{customdata[15]:,.0f}</b>",
            "Patterns: <b>%{customdata[16]:,.0f}</b>",
            "Mean pairwise identity: <b>%{customdata[17]:.2f}%</b>",
            "Missingness: <b>%{customdata[18]:.2f}%</b>",
            "GC content: <b>%{customdata[19]:.2f}%</b>",
            "GC content at 1st codon position: <b>%{customdata[20]:.2f}%</b>",
            "GC content at 2nd codon position: <b>%{customdata[21]:.2f}%</b>",
            "GC content at 3rd codon position: <b>%{customdata[22]:.2f}%</b>",
    ])
    hovertemplate_other = "<br>".join([
            "Locus: <b>%{customdata[6]}</b>",
            "Marker type: <b>%{customdata[4]}</b>",
            "Sequences: <b>%{customdata[7]:,.0f}</b>",
            "Samples: <b>%{customdata[8]:,.0f}</b>",
            "Sequences per sample: <b>%{customdata[9]:.2f}</b>",
            "Alignment length: <b>%{customdata[10]:,.0f} bp</b>",
            "Informative sites: <b>%{customdata[11]:,.0f}</b>",
            "Informativeness: <b>%{customdata[12]:.2f}%</b>",
            "Uninformative sites: <b>%{customdata[13]:,.0f}</b>",
            "Constant sites: <b>%{customdata[14]:,.0f}</b>",
            "Singleton sites: <b>%{customdata[15]:,.0f}</b>",
            "Patterns: <b>%{customdata[16]:,.0f}</b>",
            "Mean pairwise identity: <b>%{customdata[17]:.2f}%</b>",
            "Missingness: <b>%{customdata[18]:.2f}%</b>",
            "GC content: <b>%{customdata[19]:.2f}%</b>",
    ])

    for j, marker_type in enumerate(marker_type_list):
        data = df[df["marker_type"] == marker_type]
        stage_list = list(reversed(data["stage"].unique()))
        fmt_list = list(reversed(data["format"].unique()))
        fig = go.Figure()
        buttons_stage = []
        for i, stage in enumerate(stage_list):
            for fmt in fmt_list:
                d = data[(data["stage"] == stage) & (data["format"] == fmt)]
                if fmt == "AA":
                    hovertemplate = hovertemplate_aa
                elif fmt == "NT":
                    hovertemplate = hovertemplate_nt
                else:
                    hovertemplate = hovertemplate_other
                fig.add_trace(
                    go.Histogram2dContour(
                        x=d[list(var_dict.values())[3]],
                        y=d[list(var_dict.values())[4]],
                        contours_coloring="fill",
                        colorscale=[
                            [0, "rgba(8,8,8,0)"],
                            [1, color_dict[fmt]],
                        ],
                        opacity=0.5,
                        showscale=False,
                        line_width=0.1,
                        visible=True if i == 0 else False,
                        hoverinfo="skip",
                        legendgroup=fmt,
                    )
                )
                fig.add_trace(
                    go.Scatter(
                        x=d[list(var_dict.values())[3]],
                        y=d[list(var_dict.values())[4]],
                        name=fmt,
                        mode="markers",
                        visible=True if i == 0 else False,
                        legendgroup=fmt,
                        customdata=d,
                        hovertemplate=hovertemplate,
                        marker=dict(
                            size=7 if len(d) < 1000 else 5,
                            color=color_dict[fmt],
                            opacity=0.7,
                            line=dict(
                                # color="white",
                                width=1,
                            )
                        ),
                    )
                )
                fig.add_trace(
                    go.Histogram(
                        x=d[list(var_dict.values())[3]],
                        yaxis="y2",
                        name=fmt,
                        bingroup=stage + "_x",
                        hovertemplate="<br>".join([
                            "Bin: <b>%{x}</b>",
                            "Count: <b>%{y}</b>",
                        ]),
                        marker=dict(
                            color=color_dict[fmt],
                            opacity=0.5,
                            line=dict(
                                color="rgb(8,8,8)",
                                width=1,
                            ),
                        ),
                        visible=True if i == 0 else False,
                        legendgroup=fmt,
                        showlegend=False,
                    )
                )
                fig.add_trace(
                    go.Histogram(
                        y=d[list(var_dict.values())[4]],
                        xaxis="x2",
                        name=fmt,
                        bingroup=stage + "_y",
                        hovertemplate="<br>".join([
                            "Bin: <b>%{y}</b>",
                            "Count: <b>%{x}</b>",
                        ]),
                        marker=dict(
                            color=color_dict[fmt],
                            opacity=0.5,
                            line=dict(
                                color="rgb(8,8,8)",
                                width=1,
                            ),
                        ),
                        visible=True if i == 0 else False,
                        legendgroup=fmt,
                        showlegend=False,
                    )
                )
            visible = [[False] * len(fmt_list) * 4] * len(stage_list)
            visible[i] = [True] * len(fmt_list) * 4
            visible = sum(visible, [])
            button_stage = dict(
                label=stage_dict[stage],
                method="restyle",
                args=[
                    dict(
                        visible=visible
                    )
                ]
            )
            buttons_stage.append(button_stage)

        buttonsX, buttonsY = [], []
        for lab, var in var_dict.items():
            if marker_type in ["DNA", "CLR"] and "codon" in var:
                pass
            else:
                x_list, y_list = [], []
                for stage in stage_list:
                    for fmt in fmt_list:
                        d = data[(data["stage"] == stage) & (data["format"] == fmt)]
                        if var == "gc" and fmt == "AA":
                            x = [[], [], [], []]
                            y = [[], [], [], []]
                        elif "codon" in var and fmt in ["AA", "GE", "GF", "MA", "MF"]:
                            x = [[], [], [], []]
                            y = [[], [], [], []]
                        else:
                            x = [d[var], d[var], d[var], None]
                            y = [d[var], d[var], None, d[var]]
                        x_list += x
                        y_list += y
                buttonX = dict(
                    label=lab,
                    method="update",
                    args=[
                        dict(
                            x=x_list,
                        ),
                    ]
                )
                buttonY = dict(
                    label=lab,
                    method="update",
                    args=[
                        dict(
                            y=y_list,
                        ),
                    ]
                )
                buttonsX.append(buttonX)
                buttonsY.append(buttonY)

        updatemenus = [
            dict(
                buttons=buttons_stage,
                type="dropdown",
                direction="down",
                pad={"t": 10, "b": 10},
                showactive=True,
                x=0.95,
                xanchor="right",
                y=1,
                yanchor="bottom"
            ),
            dict(
                buttons=buttonsX,
                type="dropdown",
                direction="up",
                pad={"t": 30, "b": 10},
                active=3,
                showactive=True,
                x=0.475,
                xanchor="center",
                y=0,
                yanchor="top"
            ),
            dict(
                buttons=buttonsY,
                type="dropdown",
                direction="down",
                pad={"t": 10, "b": 10, "r": 40},
                active=4,
                showactive=True,
                x=0,
                xanchor="right",
                y=0.4625,
                yanchor="middle"
            ),
        ]

        annotations = [
            dict(
                text="<b>Processing Step:</b>",
                x=0.95,
                xref="paper",
                xanchor="right",
                xshift=-160,
                y=1,
                yref="paper",
                yanchor="top",
                yshift=36,
                align="right",
                showarrow=False,
            ),
        ]

        if j == 0:
            if len(marker_type_list) == 1:
                title = (
                    "<b>2. Bivariate Relationships and Distributions</b><br>"
                    + "<sup>(Source: "
                    + str(aln_stats_tsv.name)
                    + ")</sup>"
                )
            else:
                title = (
                    "<b>2. Bivariate Relationships and Distributions<br>"
                    + "Marker Type: "
                    + marker_type
                    + "</b><br><sup>(Source: "
                    + str(aln_stats_tsv.name)
                    + ")</sup>"
                )
        else:
            title = (
                "<b>Marker Type: "
                + marker_type
                + "</b><br><sup>(Source: "
                + str(aln_stats_tsv.name)
                + ")</sup>"
            )
        fig.update_layout(
            font_family="Arial",
            plot_bgcolor="rgb(8,8,8)",
            title=title,
            xaxis=dict(
                showgrid=True,
                gridcolor="rgb(64,64,64)",
                ticks="outside",
                domain=[0, 0.95],
            ),
            yaxis=dict(
                showgrid=True,
                gridcolor="rgb(64,64,64)",
                ticks="outside",
                domain=[0, 0.925],
            ),
            xaxis2=dict(
                title="Count",
                gridcolor="rgb(64,64,64)",
                ticks="outside",
                zeroline=True,
                domain=[0.95, 1],
            ),
            yaxis2=dict(
                title="Count",
                gridcolor="rgb(64,64,64)",
                ticks="outside",
                zeroline=True,
                domain=[0.925, 1],
            ),
            hoverlabel=dict(
                font_color="rgb(64,64,64)",
                bordercolor="rgb(64,64,64)",
            ),
            legend=dict(
                title=dict(
                    text="<b>Format</b>",
                    side="top",
                ),
                traceorder="reversed",
            ),
            barmode="overlay",
            annotations=annotations,
            updatemenus=updatemenus,
        )
        figs.append(fig)

    df = pd.read_table(sam_stats_tsv)
    df = df[~df["stage_marker_format"].str.contains("_w_refs")]
    df["gaps"] = df["len_gapped"] - df["len_ungapped"]
    df[["stage", "filter", "marker", "format"]] = df["stage_marker_format"].str.split(" / ", expand=True)
    for var in ["stage", "filter", "marker", "format"]:
        df[var] = df[var].str[3:]

    df = df.groupby(
        [
            "sample",
            "stage_marker_format",
            "stage",
            "filter",
            "marker",
            "format",
        ]
    ).agg(
        num_loci = ("locus", "nunique"),
        mean_len = ("len_ungapped", np.mean),
        total_len = ("len_ungapped", np.sum),
        mean_gaps = ("gaps", np.mean),
        total_gaps = ("gaps", np.sum),
        mean_ambig = ("ambigs", np.mean),
        mean_gc = ("gc", np.mean),
        mean_gc_codon1 = ("gc_codon_p1", np.mean),
        mean_gc_codon2 = ("gc_codon_p2", np.mean),
        mean_gc_codon3 = ("gc_codon_p3", np.mean),
        mean_copies = ("num_copies", np.mean),
    ).reset_index()

    var_dict = {
        "num_loci": "Number of Loci",
        "mean_len": "Mean Ungapped Length",
        "total_len": "Total Ungapped Length",
        "mean_gaps": "Mean Gaps",
        "total_gaps": "Total Gaps",
        "mean_ambig": "Mean Ambiguities",
        "mean_gc": "Mean GC Content (%)",
        "mean_gc_codon1": "Mean GC Content at<br>1st Codon Position (%)",
        "mean_gc_codon2": "Mean GC Content at<br>2nd Codon Position (%)",
        "mean_gc_codon3": "Mean GC Content at<br>3rd Codon Position (%)",
        "mean_copies": "Mean Copies",
    }

    marker_type_list = df["marker"].unique()

    color_dict = {
        "01_AA": "#CC79A7",
        "02_NT": "#56B4E9",
        "03_genes": "#E69F00",
        "04_genes_flanked": "#009E73",
        "01_matches": "#CC79A7",
        "02_matches_flanked": "#009E73",
    }

    shape_dict = {
        "04_unfiltered": "circle",
        "05_naive": "star-triangle-up",
        "06_informed": "star-square",
    }

    annotations = [
        dict(
            text="<b>Sort Samples by:</b>",
            x=1,
            xref="paper",
            xanchor="right",
            xshift=-80,
            y=1,
            yref="paper",
            yanchor="top",
            yshift=36,
            align="right",
            showarrow=False,
        ),
    ]
    buttons_sort = [
        dict(
            label="Name",
            method="relayout",
            args=[{
                "xaxis.categoryorder": "category ascending",
            }],
        ),
        dict(
            label="Value",
            method="relayout",
            args=[{
                "xaxis.categoryorder": "mean descending",
            }],
        ),
    ]
    hovertemplate_aa = "<br>".join([
        "Sample: <b>%{customdata[0]}</b>",
        "Trimming: <b>%{customdata[2]}</b>",
        "Paralog filter: <b>%{customdata[3]}</b>",
        "Marker type: <b>%{customdata[4]}</b>",
        "Format: <b>%{customdata[5]}</b>",
        "Number of loci: <b>%{customdata[6]}</b>",
        "Mean ungapped length: <b>%{customdata[7]:,.0f} aa</b>",
        "Total ungapped length: <b>%{customdata[8]:,.0f} aa</b>",
        "Mean gaps: <b>%{customdata[9]:,.0f} aa</b>",
        "Total gaps: <b>%{customdata[10]:,.0f} aa</b>",
        "Mean ambiguities: <b>%{customdata[11]:,.2f}</b>",
        "Mean copies: <b>%{customdata[16]:,.2f}</b><extra></extra>",
    ])
    hovertemplate_nt = "<br>".join([
        "Sample: <b>%{customdata[0]}</b>",
        "Trimming: <b>%{customdata[2]}</b>",
        "Paralog filter: <b>%{customdata[3]}</b>",
        "Marker type: <b>%{customdata[4]}</b>",
        "Format: <b>%{customdata[5]}</b>",
        "Number of loci: <b>%{customdata[6]}</b>",
        "Mean ungapped length: <b>%{customdata[7]:,.0f} bp</b>",
        "Total ungapped length: <b>%{customdata[8]:,.0f} bp</b>",
        "Mean gaps: <b>%{customdata[9]:,.0f} bp</b>",
        "Total gaps: <b>%{customdata[10]:,.0f} bp</b>",
        "Mean ambiguities: <b>%{customdata[11]:,.2f}</b>",
        "Mean GC content: <b>%{customdata[12]:.2f}%</b>",
        "Mean GC content (1st codon pos.): <b>%{customdata[13]:.2f}%</b>",
        "Mean GC content (2nd codon pos.): <b>%{customdata[14]:.2f}%</b>",
        "Mean GC content (3rd codon pos.): <b>%{customdata[15]:.2f}%</b>",
        "Mean copies: <b>%{customdata[16]:,.2f}</b><extra></extra>",
    ])
    hovertemplate_other = "<br>".join([
        "Sample: <b>%{customdata[0]}</b>",
        "Trimming: <b>%{customdata[2]}</b>",
        "Paralog filter: <b>%{customdata[3]}</b>",
        "Marker type: <b>%{customdata[4]}</b>",
        "Format: <b>%{customdata[5]}</b>",
        "Number of loci: <b>%{customdata[6]}</b>",
        "Mean ungapped length: <b>%{customdata[7]:,.0f} bp</b>",
        "Total ungapped length: <b>%{customdata[8]:,.0f} bp</b>",
        "Mean gaps: <b>%{customdata[9]:,.0f} bp</b>",
        "Total gaps: <b>%{customdata[10]:,.0f} bp</b>",
        "Mean ambiguities: <b>%{customdata[11]:,.2f}</b>",
        "Mean GC content: <b>%{customdata[12]:.2f}%</b>",
        "Mean copies: <b>%{customdata[16]:,.2f}</b><extra></extra>",
    ])

    for i, marker_type in enumerate(marker_type_list):
        fig = go.Figure()
        d = df[df["marker"] == marker_type]
        for stage_marker_format in d["stage_marker_format"].unique():
            stage, filter, marker, format = stage_marker_format.split(" / ")
            if format == "01_AA":
                hovertemplate = hovertemplate_aa
            elif format == "02_NT":
                hovertemplate = hovertemplate_nt
            else:
                hovertemplate = hovertemplate_other
            fig.add_trace(
                go.Scatter(
                    x=d[d["stage_marker_format"] == stage_marker_format]["sample"],
                    y=d[d["stage_marker_format"] == stage_marker_format][list(var_dict.keys())[0]],
                    mode="markers",
                    # visible=True if stage == "03_trimmed" else "legendonly",
                    name=stage_marker_format,
                    legendgroup=format,
                    # legendgrouptitle_text="<b>" + format + "</b>",
                    customdata=d[d["stage_marker_format"] == stage_marker_format],
                    hovertemplate=hovertemplate,
                    marker=dict(
                        size=7 if len(d[d["stage_marker_format"] == stage_marker_format]) < 500 else 5,
                        symbol=shape_dict[filter] if stage == "03_trimmed" else shape_dict[filter] + "-open-dot",
                        color=color_dict[format],
                        opacity=0.7 if stage == "03_trimmed" else 1,
                        line_width=1,
                    ),
                )
            )

        buttonsY = []
        for var in var_dict.keys():
            if marker_type in ["misc_DNA", "clusters"] and "codon" in var:
                pass
            else:
                y_list = []
                for stage in d["stage_marker_format"].unique():
                    if var == "mean_gc" and "01_AA" in stage:
                        y = []
                    elif "codon" in var and "02_NT" not in stage:
                        y = []
                    else:
                        y = d[d["stage_marker_format"] == stage][var]
                    y_list.append(y)
                buttonY = dict(
                    label=var_dict[var],
                    method="update",
                    args=[
                        dict(
                            y=y_list,
                        ),
                    ]
                )
                buttonsY.append(buttonY)

        updatemenus = [
            dict(
                buttons=buttonsY,
                type="dropdown",
                direction="down",
                pad={"t": 10, "b": 10, "r": 40},
                showactive=True,
                x=0,
                xanchor="right",
                y=0.5,
                yanchor="middle",
            ),
            dict(
                buttons=buttons_sort,
                type="dropdown",
                direction="down",
                pad={"t": 10, "b": 10},
                showactive=True,
                x=1,
                xanchor="right",
                y=1,
                yanchor="bottom",
            ),
        ]

        fig.update_xaxes(
            title="Sample",
            ticks="outside",
            categoryorder="category ascending",
            gridcolor="rgb(64,64,64)",
        )
        fig.update_yaxes(
            ticks="outside",
            gridcolor="rgb(64,64,64)",
            zeroline=False,
        )
        if i == 0:
            if len(marker_type_list) == 1:
                title = (
                    "<b>3. Stats Per Sample</b><br>"
                    + "<sup>(Source: "
                    + str(sam_stats_tsv.name)
                    + ")</sup>"
                )
            else:
                title = (
                    "<b>3. Stats Per Sample<br>"
                    + "Marker Type: "
                    + marker_type
                    + "</b><br><sup>(Source: "
                    + str(sam_stats_tsv.name)
                    + ")</sup>"
                )
        else:
            title = (
                "<b>Marker Type: "
                + marker_type
                + "</b><br><sup>(Source: "
                + str(sam_stats_tsv.name)
                + ")</sup>"
            )
        fig.update_layout(
            title_text=title,
            font_family="Arial",
            plot_bgcolor="rgb(8,8,8)",
            legend=dict(
                title="<b>Trimming / Paralog Filter / Marker Type / Format",
                groupclick="toggleitem",
            ),
            hoverlabel=dict(
                font_color="rgb(64,64,64)",
                bordercolor="rgb(64,64,64)",
            ),
            annotations=annotations,
            updatemenus=updatemenus,
        )
        figs.append(fig)

    # Save plot in html
    aln_html_report = Path(out_dir, "captus-assembly_align.report.html")
    with open(aln_html_report, "w") as f:
        for i, fig in enumerate(figs):
            f.write(
                fig.to_html(
                    full_html=False,
                    include_plotlyjs="cdn",
                    config=dict(
                        scrollZoom=False if i == 0 else True,
                        toImageButtonOptions=dict(
                            format="svg",
                        ),
                        modeBarButtonsToAdd=[
                            "v1hovermode",
                            "hoverclosest",
                            "hovercompare",
                            "togglehover",
                            "togglespikelines",
                            "drawline",
                            "drawopenpath",
                            "drawclosedpath",
                            "drawcircle",
                            "drawrect",
                            "eraseshape",
                        ]
                    ),
                )
            )
    if aln_html_report.exists() and aln_html_report.is_file():
        aln_html_msg = dim(f"Report generated in {elapsed_time(time.time() - start)}")
    else:
        aln_html_msg = red(f"Report not generated, verify your Python environment")

    return aln_html_report, aln_html_msg

def build_design_report(out_dir, des_stats_tsv, step):
    start = time.time()

    df = pd.read_table(des_stats_tsv)
    df.drop(columns=["path"], inplace=True)
    var_dict = {
        "Length (bp)": "length",
        "GC content (%)": "gc_content",
        "Mean pairwise identity (%)": "avg_pid",
        "Informative sites": "informative_sites",
        "Informativeness (%)": "informativeness",
        "Missingness (%)": "missingness",
        "Sequences": "sequences",
        "Samples": "samples",
        "Focal species": "focal_species",
        "Outgroup species": "outgroup_species",
        "Add-on samples": "addon_samples",
        "Species": "species",
        "Genera": "genera",
        "CDS length (bp)": "cds_len",
        "Length of long exons retained (bp)": "len_long_exons_retained",
        "Length of short exons retained (bp)": "len_short_exons_retained",
        "Percentage of exons retained (%)": "perc_exons_retained",
        "Percentage of long exons retained (%)": "perc_long_exons_retained",
        "Percentage of short exons retained (%)": "perc_short_exons_retained",
    }
    hovertemplate = "<br>".join([
        "Locus: <b>%{customdata[0]}</b>",
        "Single-copy: <b>%{customdata[1]}</b>",
        "Length: <b>%{customdata[2]:,.0f} bp</b>",
        "GC content: <b>%{customdata[3]:,.2f}%</b>",
        "Mean pairwise identity: <b>%{customdata[4]:.2f}%</b>",
        "Informative sites: <b>%{customdata[5]:,.0f}</b>",
        "Informativeness: <b>%{customdata[6]:.2f}%</b>",
        "Missingness: <b>%{customdata[7]:,.2f}%</b>",
        "Sequences: <b>%{customdata[8]:,.0f}</b>",
        "Samples: <b>%{customdata[9]:,.0f}</b>",
        "Focal species: <b>%{customdata[10]:,.0f}</b>",
        "Outgroup species: <b>%{customdata[11]:,.0f}</b>",
        "Add-on samples: <b>%{customdata[12]:,.0f}</b>",
        "Species: <b>%{customdata[13]:,.0f}</b>",
        "Genera: <b>%{customdata[14]:,.0f}</b>",
        "CDS id: <b>%{customdata[15]}</b>",
        "CDS length: <b>%{customdata[16]:,.0f} bp</b>",
        "Length of long exons retained: <b>%{customdata[17]:,.0f} bp</b>",
        "Length of short exons retained: <b>%{customdata[18]:,.0f} bp</b>",
        "Percentage of exons retained: <b>%{customdata[19]:.2f}%</b>",
        "Percentage of long exons retained: <b>%{customdata[20]:.2f}%</b>",
        "Percentage of short exons retained: <b>%{customdata[21]:.2f}%</b><extra></extra>",
    ])
    figs = []
    fig = go.Figure()
    for single_copy in sorted(df["single_copy"].unique(), reverse=True):
        d = df.query('single_copy == @single_copy')
        x=d[list(var_dict.values())[0]]
        y=d[list(var_dict.values())[3]]
        color = "#56B4E9" if single_copy == True else "#CC79A7"
        fig.add_trace(
            go.Histogram2dContour(
                x=x,
                y=y,
                contours_coloring="fill",
                colorscale=[
                    [0, "rgba(8,8,8,0)"],
                    [1, color],
                ],
                opacity=0.5,
                showscale=False,
                line_width=0.1,
                hoverinfo="skip",
                legendgroup=str(single_copy),
            )
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name=str(single_copy),
                mode="markers",
                legendgroup=str(single_copy),
                customdata=d,
                hovertemplate=hovertemplate,
                marker=dict(
                    size=7 if len(df) < 1000 else 4,
                    color=color,
                    opacity=0.7,
                    line=dict(
                        width=1,
                    )
                ),
            )
        )
        fig.add_trace(
            go.Histogram(
                x=x,
                yaxis="y2",
                bingroup="x",
                hovertemplate="<br>".join([
                    "Bin: <b>%{x}</b>",
                    "Count: <b>%{y}</b><extra></extra>",
                ]),
                marker=dict(
                    color=color,
                    opacity=0.7,
                    line=dict(
                        color="rgb(8,8,8)",
                        width=1,
                    ),
                ),
                legendgroup=str(single_copy),
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Histogram(
                y=y,
                xaxis="x2",
                bingroup="y",
                hovertemplate="<br>".join([
                    "Bin: <b>%{y}</b>",
                    "Count: <b>%{x}</b><extra></extra>",
                ]),
                marker=dict(
                    color=color,
                    opacity=0.7,
                    line=dict(
                        color="rgb(8,8,8)",
                        width=1,
                    ),
                ),
                legendgroup=str(single_copy),
                showlegend=False,
            )
        )
    buttonsX, buttonsY = [], []
    for lab, var in var_dict.items():
        x_list, y_list = [], []
        for single_copy in sorted(df["single_copy"].unique(), reverse=True):
            d = df.query('single_copy == @single_copy')
            x = [d[var], d[var], d[var], None]
            y = [d[var], d[var], None, d[var]]
            x_list += x
            y_list += y
        buttonX = dict(
            label=lab,
            method="update",
            args=[
                dict(
                    x=x_list,
                ),
            ]
        )
        buttonY = dict(
            label=lab,
            method="update",
            args=[
                dict(
                    y=y_list,
                ),
            ]
        )
        buttonsX.append(buttonX)
        buttonsY.append(buttonY)

    updatemenus = [
            dict(
                buttons=buttonsX,
                type="dropdown",
                direction="up",
                pad={"t": 30, "b": 10},
                active=0,
                showactive=True,
                x=0.475,
                xanchor="center",
                y=0,
                yanchor="top"
            ),
            dict(
                buttons=buttonsY,
                type="dropdown",
                direction="down",
                pad={"t": 10, "b": 10, "r": 40},
                active=3,
                showactive=True,
                x=0,
                xanchor="right",
                y=0.4625,
                yanchor="middle"
            ),
        ]
    if step == "cluster":
        des_html_report = Path(out_dir, "captus-design_cluster.report.html")
        title = (
            "<b>Captus-design: Cluster (Alignment Report)</b><br>"
            + "<sup>(Source: "
            + str(des_stats_tsv.name)
            + ")</sup>"
        )
    elif step == "select":
        des_html_report = Path(out_dir, "captus-design_select.report.html")
        title = (
            "<b>Captus-design: Select (Alignment Report)</b><br>"
            + "<sup>(Source: "
            + str(des_stats_tsv.name)
            + ")</sup>"
        )
    fig.update_layout(
        font_family="Arial",
        plot_bgcolor="rgb(8,8,8)",
        title=title,
        xaxis=dict(
            showgrid=True,
            gridcolor="rgb(64,64,64)",
            ticks="outside",
            domain=[0, 0.95],
        ),
        yaxis=dict(
            showgrid=True,
            gridcolor="rgb(64,64,64)",
            ticks="outside",
            domain=[0, 0.925],
        ),
        xaxis2=dict(
            title="Count",
            gridcolor="rgb(64,64,64)",
            ticks="outside",
            zeroline=True,
            domain=[0.95, 1],
        ),
        yaxis2=dict(
            title="Count",
            gridcolor="rgb(64,64,64)",
            ticks="outside",
            zeroline=True,
            domain=[0.925, 1],
        ),
        hoverlabel=dict(
            font_color="rgb(64,64,64)",
            bordercolor="rgb(64,64,64)",
        ),
        legend=dict(
            title=dict(
                text="<b>Single-copy</b>",
                side="top",
            ),
            itemsizing="constant",
        ),
        barmode="overlay",
        updatemenus=updatemenus,
    )
    figs.append(fig)
    with open(des_html_report, "w") as f:
        for fig in figs:
            f.write(
                fig.to_html(
                    full_html=False,
                    include_plotlyjs="cdn",
                    config=dict(
                        scrollZoom=True,
                        toImageButtonOptions=dict(
                            format="svg",
                        ),
                        modeBarButtonsToAdd=[
                            "v1hovermode",
                            "hoverclosest",
                            "hovercompare",
                            "togglehover",
                            "togglespikelines",
                            "drawline",
                            "drawopenpath",
                            "drawclosedpath",
                            "drawcircle",
                            "drawrect",
                            "eraseshape",
                        ]
                    ),
                )
            )
    if des_html_report.exists() and des_html_report.is_file():
        des_html_msg = dim(f"Report generated in {elapsed_time(time.time() - start)}")
    else:
        des_html_msg = red(f"Report not generated, verify your Python environment")

    return des_html_report, des_html_msg
