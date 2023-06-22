#!/usr/bin/env python3

"""
File:         gene_cluster_plot.py
Created:      2020/04/27
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
from __future__ import print_function
from pathlib import Path
import math
import os

# Third party imports.
from colour import Color
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Gene Cluster plot"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class main():
    def __init__(self):
        self.input_dirs = ["cis_output", "alzheimer_output", "ms_output",
                           "schizophrenia_output", "alzheimer_trans_output",
                           "ms_trans_output", "schizophrenia_trans_output",
                           "depression_output", "parkinson_output ",
                           "depression_trans_output", "parkinson_trans_output"]
        self.colormap = {"Neuron": ["#e8e5e5", "#b38d84"],
                         "Oligodendrocyte": ["#dfdddb", "#5d9166"],
                         "EndothelialCell": ["#f5f0f0", "#f2a7a7"],
                         "Microglia": ["#eee8e6", "#e8c06f"],
                         "Macrophage": ["#eee8e6", "#e8c06f"],
                         "Astrocyte": ["#e7e8e4", "#9b7bb8"]}
        self.method = "CellMapNNLS_"
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'pca_plots')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        for dir in self.input_dirs:
            print("Working on dir: {}".format(dir))

            interaction_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/{}/covariates/interaction_table.txt.gz".format(
                dir)
            expression_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/{}/create_matrices/expression_table.txt.gz".format(
                dir)

            print("\tLoad the interaction z-scores")
            inter_df = pd.read_csv(interaction_path, sep="\t", header=0,
                                   index_col=0)
            print("\t\tshape: {}".format(inter_df.shape))

            print("\tSubsetting interactions of interest")
            inter_df = inter_df.loc[inter_df.index.str.startswith(self.method), :]
            print("\t\tshape: {}".format(inter_df.shape))

            print("\tLoad expression")
            expr_df = pd.read_csv(expression_path, sep="\t", header=0,
                                  index_col=0)
            print("\t\tshape: {}".format(expr_df.shape))

            print("\tPerforming PCA")
            projections, var_expl = self.pca(expr_df)

            print("\tPlotting")
            for x_pc in range(50):
                for y_pc in range(50):
                    if x_pc < y_pc:
                        for index, row in inter_df.iterrows():
                            celltype = index.replace(self.method, "")

                            print("\t\tPC{} vs PC{} [color: {}]".format(x_pc,
                                                                        y_pc,
                                                                        celltype))

                            hue = row.round(2)
                            hue.name = celltype

                            colormap = self.create_colormap(hue.min(),
                                                            hue.max(),
                                                            self.colormap[
                                                                celltype])

                            self.plot_single_pca(projections, var_expl, hue,
                                                 x_pc, y_pc, celltype, colormap,
                                                 dir, self.outdir)
            exit()

    @staticmethod
    def load_hits(path):
        data = []
        cols = []

        with open(path, 'r') as f:
            for i, line in enumerate(f):
                if len(line) == 1:
                    break
                line = line.strip("\n").split("\t")
                if i == 0:
                    cols = line
                else:
                    data.append(line)
        f.close()

        return pd.DataFrame(data, columns=cols)

    @staticmethod
    def create_colormap(min, max, color, precision=100):
        print("\t\tmin: {}\tmax: {}\tcolor: {}".format(min, max, color))
        color_min = 0
        if min < 0:
            color_min = min
        grey_range = math.ceil(abs(color_min * precision)) + 1
        color_range = math.ceil(max * precision) + 1

        palette = list(
            Color("#000000").range_to(Color("#FFFFFF"), grey_range)) + \
                  list(Color(color[0]).range_to(Color(color[1]), color_range))
        colors = [str(x).upper() for x in palette]
        values = [x / precision for x in range(-1 * grey_range, color_range)]
        color_map = {}
        for val, col in zip(values, colors):
            color_map[val] = col
        return color_map

    @staticmethod
    def pca(data):
        correlation_matrix = np.dot(data.T, data) / (data.shape[0] - 1)

        eigenvalues, eigenvectors = np.linalg.eig(correlation_matrix)
        order = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]

        ex_variance_pcnt = (eigenvalues / eigenvalues.sum()) * 100

        projections = np.dot(data, eigenvectors)
        return projections, ex_variance_pcnt

    @staticmethod
    def plot_single_pca(data, var_expl, hue, x_pc, y_pc, celltype, palette,
                        title, outdir):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")

        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        ax.axvline(0, ls='--', color="#000000", alpha=0.15, zorder=-2)
        ax.axhline(0, ls='--', color="#000000", alpha=0.15, zorder=-2)

        x = data[:, x_pc].copy()
        y = data[:, y_pc].copy()
        df = pd.DataFrame({"x": x.real, "y": y.real, "z": hue})
        negative = df.loc[df["z"] <= 0, :].copy()
        positive = df.loc[df["z"] > 0, :].copy()

        for df, zorder in zip([negative, positive], [-1, 1]):
            sns.scatterplot(x="x",
                            y="y",
                            hue="z",
                            data=df,
                            palette=palette,
                            s=200,
                            legend=False,
                            zorder=zorder,
                            ax=ax)
        ax.set_title(
            "{} - {}".format(title.split("_")[0].lower(), celltype.lower()),
            fontsize=18,
            fontweight='bold')
        ax.set_xlabel("PC{} [{:.2f}%]".format(x_pc, var_expl[x_pc].real),
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel("PC{} [{:.2f}%]".format(y_pc, var_expl[y_pc].real),
                      fontsize=14,
                      fontweight='bold')

        plt.tight_layout()
        plt.show()
        fig.savefig(os.path.join(outdir,
                                 "{}_PC{}_PC{}_{}.png".format(title, x_pc, y_pc,
                                                              celltype)))
        plt.close(fig)
        plt.clf()


if __name__ == '__main__':
    m = main()
    m.start()
