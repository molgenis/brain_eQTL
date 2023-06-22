#!/usr/bin/env python3

"""
File:         plot_pca_outliers.py
Created:      2021/09/28
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
import time
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Plot PCA Outliers"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
Syntax: 
./plot_pca_outliers.py -i /groups/umcg-biogen/tmp01/output/2021-FreezeThree/2021-02-18-splicing/2021-07-22-cortex-rmats-data/7-removeOutlierSamples/output
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.input_path = getattr(arguments, 'input')
        self.cut_off = 3

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
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
                            help="show program's version number and exit.")
        parser.add_argument("-i",
                            "--input",
                            type=str,
                            required=True,
                            help="The path to the input directory.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        for fpath in glob.glob(os.path.join(self.input_path, "*.txt")):
            name = os.path.basename(fpath).replace(".txt", "")
            print("Working on '{}'".format(name))

            print("\tLoading file.")
            df = self.load_file(fpath, header=0, index_col=0)

            if "PC1" not in df.columns or "PC2" not in df.columns:
                print("Error, PC1 and/or PC2 not in input data.")
                exit()

            print("\tSubsetting columns")
            df = df.loc[:, ["PC1", "PC2"]]

            print("\tCalculating z-scores")
            zscores_df = self.perform_zscore_transform(df=df)
            zscores_df.columns = ["{} z-score".format(x) for x in df.columns]

            print("\tMerge.")
            df = df.merge(zscores_df, left_index=True, right_index=True)
            print(df)

            print("\tAdding color.")
            df["outlier"] = False
            df.loc[(df["PC1 z-score"].abs() > self.cut_off) | (df["PC2 z-score"].abs() > self.cut_off), "outlier"] = True

            print("\tPlotting")
            self.plot(df=df,
                      x="PC1 z-score",
                      y="PC2 z-score",
                      hue="outlier",
                      palette={True: "#0072B2", False: "#808080"},
                      line_pos=self.cut_off,
                      title=name,
                      name=name)

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def perform_zscore_transform(df):
        return df.subtract(df.mean(axis=0), axis=1).divide(df.std(axis=0), axis=1)

    def plot(self, df, x="x", y="y", hue=None, palette=None, line_pos=None,
             xlabel=None, ylabel=None, title="", name="PCA_plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, gridspec_kw={"width_ratios": [0.9, 0.1]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        s=100,
                        linewidth=0,
                        legend=None,
                        palette=palette,
                        ax=ax1)

        if line_pos is not None:
            if df[x].min() < (-1 * line_pos):
                ax1.axvline(-1 * line_pos, ls='--', color="#000000", alpha=0.15, zorder=-1)

            if df[x].max() > line_pos:
                ax1.axvline(line_pos, ls='--', color="#000000", alpha=0.15, zorder=-1)

            if df[y].min() < (-1 * line_pos):
                ax1.axhline(-1 * line_pos, ls='--', color="#000000", alpha=0.15, zorder=-1)

            if df[y].max() > line_pos:
                ax1.axhline(line_pos, ls='--', color="#000000", alpha=0.15, zorder=-1)

        ax1.annotate(
            'N = {:,}'.format(df.shape[0]),
            xy=(0.03, 0.94),
            xycoords=ax1.transAxes,
            color="#000000",
            fontsize=18,
            fontweight='bold')
        ax1.annotate(
            'N = {:,}'.format(df.loc[df[hue] == True, :].shape[0]),
            xy=(0.03, 0.9),
            xycoords=ax1.transAxes,
            color="#0072B2",
            fontsize=18,
            fontweight='bold')

        ax1.set_title(title,
                      fontsize=20,
                      fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_xlabel(xlabel,
                       fontsize=14,
                       fontweight='bold')

        fig.savefig(os.path.join(self.outdir, "{}_PCA_OutlierPlot.png".format(name)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Input path: {}".format(self.input_path))
        print("  > Output directory {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
