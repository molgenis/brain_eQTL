#!/usr/bin/env python3

"""
File:         plot_pca.py
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
__program__ = "Plot PCA"
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
./plot_pca.py -i /groups/umcg-biogen/tmp01/output/2021-FreezeThree/2021-02-18-splicing/2021-07-22-cortex-rmats-data/10-removeSplicingPcs/pc-output -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.input_path = getattr(arguments, 'input')
        self.std_path = getattr(arguments, 'sample_to_dataset')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
           "GTE-EUR-AMPAD-MAYO-V2": "#9c9fa0",
           "GTE-EUR-CMC_HBCC_set2": "#0877b4",
           "GTE-EUR-GTEx": "#0fa67d",
           "GTE-EUR-AMPAD-ROSMAP-V2": "#6950a1",
           "GTE-EUR-BrainGVEX-V2": "#48b2e5",
           "GTE-EUR-TargetALS": "#d5c77a",
           "GTE-EUR-AMPAD-MSBB-V2": "#5cc5bf",
           "GTE-EUR-NABEC-H610": "#6d743a",
           "GTE-EUR-LIBD_1M": "#e49d26",
           "GTE-EUR-ENA": "#d46727",
           "GTE-EUR-LIBD_h650": "#e49d26",
           "GTE-EUR-GVEX": "#000000",
           "GTE-EUR-NABEC-H550": "#6d743a",
           "GTE-EUR-CMC_HBCC_set3": "#0877b4",
           "GTE-EUR-UCLA_ASD": "#f36d2a",
           "GTE-EUR-CMC": "#eae453",
        }

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
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")


        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading sample to dataset")
        std_df = self.load_file(self.std_path, header=0, index_col=None)
        std_dict = dict(zip(std_df.iloc[:, 0], std_df.iloc[:, 1]))

        for fpath in glob.glob(os.path.join(self.input_path, "*.txt")):
            name = os.path.basename(fpath).replace(".txt", "")
            print("Working on '{}'".format(name))

            print("\tLoading file.")
            df = self.load_file(fpath, header=0, index_col=0)

            print("\tAdding color.")
            overlap = set(df.index.values).intersection(set(std_df.iloc[:, 0].values))
            if len(overlap) != df.shape[0]:
                print("Error, some samples do not have a dataset.")
                exit()
            df["dataset"] = df.index.map(std_dict)

            print("\tPlotting")
            self.plot(df=df,
                      columns=df.columns[:5],
                      hue="dataset",
                      palette=self.palette,
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

    def plot(self, df, columns, hue, palette, name):
        ncols = len(columns)
        nrows = len(columns)

        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='col',
                                 sharey='row',
                                 figsize=(10 * ncols, 10 * nrows))
        sns.set(color_codes=True)
        sns.set_style("ticks")

        for i, y_col in enumerate(columns):
            for j, x_col in enumerate(columns):
                print(i, j)
                ax = axes[i, j]
                if i == 0 and j == (ncols - 1):
                    ax.set_axis_off()
                    if hue is not None and palette is not None:
                        groups_present = df[hue].unique()
                        handles = []
                        for key, value in self.palette.items():
                            if key in groups_present:
                                handles.append(mpatches.Patch(color=value, label=key))
                        ax.legend(handles=handles, loc=4, fontsize=25)

                elif i < j:
                    ax.set_axis_off()
                    continue
                elif i == j:
                    ax.set_axis_off()

                    ax.annotate(y_col,
                                xy=(0.5, 0.5),
                                ha='center',
                                xycoords=ax.transAxes,
                                color="#000000",
                                fontsize=40,
                                fontweight='bold')
                else:
                    sns.despine(fig=fig, ax=ax)

                    sns.scatterplot(x=x_col,
                                    y=y_col,
                                    hue=hue,
                                    data=df,
                                    s=100,
                                    palette=palette,
                                    linewidth=0,
                                    legend=False,
                                    ax=ax)

                    ax.set_ylabel("",
                                  fontsize=20,
                                  fontweight='bold')
                    ax.set_xlabel("",
                                  fontsize=20,
                                  fontweight='bold')

        fig.savefig(os.path.join(self.outdir, "{}_PCA_ColorByCohort.png".format(name)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Input path: {}".format(self.input_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Output directory {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
