#!/usr/bin/env python3

"""
File:         visualise_gene_celltype_correlations.py
Created:      2020/09/09
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
import argparse
import gzip
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Visualise Gene Cell Type Correlations"
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
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.df_path = getattr(arguments, 'data')
        self.extension = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "gene_cellcount%_corr")
        self.colormap = {
            "": "#FFFFFF",
            "Neuron": "#b38d84",
            "Oligodendrocyte": "#5d9166",
            "EndothelialCell": "#f2a7a7",
            "Microglia": "#e8c06f",
            "Macrophage": "#e8c06f",
            "Astrocyte": "#9b7bb8",
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the gene - cell fracttions "
                                 "correlation matrix.")
        parser.add_argument("-e",
                            "--extension",
                            type=str,
                            choices=["png", "pdf"],
                            default="png",
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        print("Loading gene - cell fractions correlation matrix.")
        df = pd.read_csv(self.df_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.df_path),
                                      df.shape))
        print(df)

        self.plot_histogram(df)
        self.plot_clustermap(df)

    def plot_histogram(self, df):
        dfm = df.copy().melt()

        sns.set(style="ticks", color_codes=True)
        g = sns.FacetGrid(dfm, col='variable', sharex=True, sharey=True,
                          palette=self.colormap)
        g.map(sns.distplot, 'value')
        g.map(self.vertical_mean_line, 'value')
        g.set_titles('{col_name}')
        plt.tight_layout()
        g.savefig(os.path.join(self.outdir, "correlation_distributions.{}".format(self.extension)))
        plt.close()

    def plot_clustermap(self, df):
        df_tmp = df.copy().T
        row_colors = [self.colormap[ct] for ct in df_tmp.index]

        sns.set(color_codes=True)
        g = sns.clustermap(df_tmp, center=0, cmap="RdBu_r",
                           row_colors=row_colors,
                           yticklabels=True, xticklabels=False,
                           dendrogram_ratio=(.1, .1),
                           figsize=(12, 9))
        plt.setp(g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=10))
        g.fig.subplots_adjust(bottom=0.05, top=0.7)
        plt.tight_layout()
        g.savefig(os.path.join(self.outdir, "correlation_clustermap.{}".format(self.extension)))
        plt.close()

    @staticmethod
    def vertical_mean_line(x, **kwargs):
        plt.axvline(x.mean(), ls="--", c="black")


if __name__ == '__main__':
    m = main()
    m.start()
