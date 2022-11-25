#!/usr/bin/env python3

"""
File:         covariate_cluster_plot.py
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
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Covariate Cluster plot"
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
        self.input_dir = "cis_output"
        self.input_file = "interaction_table"
        self.label_dict = {"Neuron": "Neuron",
                           "Oligodendrocyte": "Oligodendrocyte",
                           "EndothelialCell": "EndothelialCell",
                           "Microglia": "Microglia\nMacrophage",
                           "Macrophage": "Microglia\nMacrophage",
                           "Astrocyte": "Astrocyte",
                           "Comp": "PC",
                           "SEX": "Sex",
                           "before-cov-correction": "Context"}
        self.colormap = {"Neuron": "#b38d84",
                         "Oligodendrocyte": "#5d9166",
                         "EndothelialCell": "#f2a7a7",
                         "Microglia\nMacrophage": "#e8c06f",
                         "Astrocyte": "#9b7bb8",
                         "PC": "#808080",
                         "Sex": "#B22222",
                         "Context": "#000000"}
        self.labels = {"McKenzie": [2, ""], "CellMap": [0, "CellMap"]}
        self.markers = {"Other": ".", "McKenzie": "d", "PCA": "^", "NMF": "s", "NNLS": "P"}
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'pca_plots')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        interaction_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/{}/covariates/{}.txt.gz".format(
            self.input_dir, self.input_file)
        # cov_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/{}/create_cov_matrix/{}.txt.gz".format(
        #     self.input_dir, self.input_file)

        print("Load {}".format(self.input_file))
        data_df = pd.read_csv(interaction_path, sep="\t", header=0, index_col=0)
        print("\tshape: {}".format(data_df.shape))

        # print("Performing PCA")
        # projections, var_expl = self.pca(data_df)
        # print("\tshape: {}".format(projections.shape))
        # print("\tshape: {}".format(var_expl.shape))

        projections_file = os.path.join(self.outdir, "projections_{}_{}.npy".format(self.input_dir, self.input_file))
        var_expl_file = os.path.join(self.outdir, "var_expl_{}_{}.npy".format(self.input_dir, self.input_file))
        # np.save(projections_file, projections)
        # np.save(var_expl_file, var_expl)

        projections = np.load(projections_file)
        var_expl = np.load(var_expl_file)

        print("Creating groups")
        groups = []
        labels = []
        styles = []
        for index in data_df.index:
            group = "Other"
            for key, value in self.label_dict.items():
                if key in index:
                    group = value
                    break

            label = ""
            for key, (label_index, str_replace) in self.labels.items():
                if index.startswith(key):
                    label = index.split("_")[label_index].replace(str_replace, "")
                    break
            labels.append(label)

            style = "Other"
            for key, value in self.markers.items():
                if key in index:
                    style = key
                    break
            styles.append(style)

            groups.append(group)

        print("Plotting")
        for x_pc in range(5):
            for y_pc in range(5):
                if x_pc < y_pc:
                    self.plot_single_pca(projections, var_expl, groups, labels, styles,
                                         x_pc, y_pc, self.colormap, self.markers,
                                         self.input_file, self.outdir)

    @staticmethod
    def pca(data):
        data.dropna(axis=1, inplace=True)
        zscores = data.subtract(data.mean(axis=1), axis=0).divide(data.std(axis=1), axis=0)

        correlation_matrix = np.dot(zscores.T, zscores) / (zscores.shape[0] - 1)

        eigenvalues, eigenvectors = np.linalg.eig(correlation_matrix)
        order = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]

        ex_variance_pcnt = (eigenvalues / eigenvalues.sum()) * 100

        projections = np.dot(data, eigenvectors)
        return projections, ex_variance_pcnt

    @staticmethod
    def label_point(x, y, val, ax):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            ax.text(point['x'] + .02, point['y'], str(point['val']))

    @staticmethod
    def plot_single_pca(data, var_expl, hue, labels, styles,  x_pc, y_pc,
                        palette, markers, file_suffix, outdir):
        # if x_pc != 0 or y_pc != 4:
        #     return
        sns.set(rc={'figure.figsize': (8.4, 6.3)})
        sns.set_style("ticks")

        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        ax.axvline(0, ls='--', color="#000000", alpha=0.15, zorder=-2)
        ax.axhline(0, ls='--', color="#000000", alpha=0.15, zorder=-2)

        df = pd.DataFrame({'x': data[:, x_pc].copy(), 'y': data[:, y_pc].copy(),
                           'Color': hue, 'Style': styles, 'label': labels})

        sns.scatterplot(x='x',
                        y='y',
                        hue='Color',
                        style='Style',
                        data=df,
                        palette=palette,
                        markers=markers,
                        s=250,
                        ax=ax)

        # for i, point in df.iterrows():
        #     ax.text(point['x'].real + .02, point['y'].real, point['label'])

        ax.set_title("",
                     fontsize=18,
                     fontweight='bold')
        ax.set_xlabel("PC{} [{:.2f}%]".format(x_pc + 1, var_expl[x_pc].real),
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel("PC{} [{:.2f}%]".format(y_pc + 1, var_expl[y_pc].real),
                      fontsize=14,
                      fontweight='bold')

        plt.legend(bbox_to_anchor=(1.05, 0.75), borderaxespad=0.)
        plt.setp(ax.get_legend().get_texts(), fontsize='10')
        plt.setp(ax.get_legend().get_title(), fontsize='12', fontweight='bold')

        plt.show()
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "PC{}_PC{}_{}.png".format(x_pc, y_pc, file_suffix)))
        plt.close(fig)
        plt.clf()


if __name__ == '__main__':
    m = main()
    m.start()
