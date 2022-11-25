#!/usr/bin/env python3

"""
File:         visualise_single_cell_expression.py
Created:      2020/10/27
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
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Visualise Single Cell Expression"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.input_dir_path = getattr(arguments, 'input_dir')
        self.interest = getattr(arguments, 'interest')
        self.extension = getattr(arguments, 'extension')

        self.colormap = {
            "In": "#56B4E9",
            "Ex": "#0072B2",
            "Oli": "#009E73",
            "End": "#CC79A7",
            "Mic": "#E69F00",
            "Ast": "#D55E00",
            "Per": "#808080",
            "Opc": "#F0E442"
        }

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'sc_expression')

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
        parser.add_argument("-dir",
                            "--input_dir",
                            type=str,
                            required=True,
                            help="The input directory path.")
        parser.add_argument("-i",
                            "--interest",
                            nargs="*",
                            type=str,
                            default=[],
                            help="The indices of interst.")
        parser.add_argument("-e",
                            "--extension",
                            type=str,
                            choices=["png", "pdf"],
                            default="png",
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        print("### STEP1 ###")
        print("Loading data.")
        data_frames = self.load_data()

        print("### STEP2 ###")
        print("Calculating expression info.")
        expr_info = self.calc_expr_info(data_frames)

        print("### STEP3 ###")
        print("Subsetting data.")
        for interest in self.interest:
            gene_sc_df = self.subset(data_frames, interest)
            print(gene_sc_df)

            gene_expr_info = self.subset(expr_info, interest)
            print(gene_expr_info)

            self.plot(df=gene_sc_df.T.melt(),
                      name=interest,
                      ylabel="normalized expression",
                      xlabel="cell type")

    def load_data(self):
        sc_data = {}

        for filepath in glob.glob(os.path.join(self.input_dir_path + "*_expression.tsv")):
            celltype = os.path.basename(filepath).split("_")[0]

            print("\tLoading {} expression matrix.".format(celltype))

            sc_data[celltype] = pd.read_csv(filepath,
                                            sep="\t",
                                            header=0,
                                            index_col=0)

        return sc_data

    @staticmethod
    def calc_expr_info(dataframes):
        expr_info = {}

        for key, df in dataframes.items():
            print("\tAnalyzing {}".format(key))
            new_indicies = []
            new_data = []

            for index, row in df.iterrows():
                new_indicies.append(index)
                new_data.append([np.min(row), np.max(row), np.mean(row), np.median(row), np.std(row), np.sum(row)])

            info_df = pd.DataFrame(new_data,
                                   index=new_indicies,
                                   columns=["min", "max", "mean", "median", "sd", "sum"])
            info_df["rank"] = info_df.loc[:, "sum"].rank(ascending=False)
            info_df["total"] = info_df.shape[0] + 1
            info_df["pcnt"] = info_df["rank"] / info_df["total"]

            expr_info[key] = info_df

        return expr_info

    @staticmethod
    def subset(data_frames, index):
        print("\tSubsetting {} data.".format(index))

        combined_df = None

        for celltype, df in data_frames.items():
            df = df.T
            if index in df.columns:
                subset = df[[index]].copy()
                subset.columns = [celltype]
                if combined_df is None:
                    combined_df = subset
                else:
                    combined_df = combined_df.merge(subset,
                                                    left_index=True,
                                                    right_index=True)

        return combined_df.T

    def plot(self, df, xlabel="", ylabel="", name=""):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.boxplot(x="variable",
                    y="value",
                    data=df,
                    palette=self.colormap,
                    ax=ax)
        fig.suptitle(name, fontsize=20, fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_boxplot.{}".format(name, self.extension)))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
