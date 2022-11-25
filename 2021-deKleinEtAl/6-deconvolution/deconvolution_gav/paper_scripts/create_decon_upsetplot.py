#!/usr/bin/env python3

"""
File:         create_decon_upsetplot.py
Created:      2020/11/24
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
import itertools
import argparse
import os

# Third party imports.
import pandas as pd
import matplotlib
from statsmodels.stats import multitest
matplotlib.use('Agg')
import upsetplot as up
import matplotlib.pyplot as plt

# Local application imports.

"""
Syntax:
./create_decon_upsetplot.py -d ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -e pdf
"""


# Metadata
__program__ = "Create Decon UpsetPlot"
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
        self.decon_path = getattr(arguments, 'decon')
        self.alpha = getattr(arguments, 'alpha')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'create_decon_upsetplot')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix")
        parser.add_argument("-a",
                            "--alpha",
                            type=float,
                            required=False,
                            default=0.05,
                            help="The significance cut-off. Default: 0.05.")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()
        print("Loading data.")
        decon_df = self.load_file(self.decon_path, nrows=None)

        print("Calculating FDR.")
        _, decon_fdr_df = self.bh_correct(decon_df)

        print("Preprocessing data.")
        data = self.parse_df(decon_fdr_df, self.alpha)
        counts = self.count(data)
        counts = counts[counts > 0]
        print(counts)

        print("Creating plot.")
        up.plot(counts, sort_by='cardinality', show_counts=True)
        for extension in self.extensions:
            plt.savefig(os.path.join(self.outdir, "eQTL_upsetplot.{}".format(extension)))
        plt.close()

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def bh_correct(pvalue_df):
        df = pvalue_df.copy()
        pval_data = []
        fdr_data = []
        indices = []
        for col in df.columns:
            if col.endswith("_pvalue"):
                pval_data.append(df.loc[:, col])
                fdr_data.append(multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
                indices.append(col.replace("_pvalue", ""))
        pval_df = pd.DataFrame(pval_data, index=indices, columns=df.index)
        fdr_df = pd.DataFrame(fdr_data, index=indices, columns=df.index)

        return pval_df.T, fdr_df.T

    @staticmethod
    def parse_df(df, alpha):
        tmp_df = df.copy()
        tmp_df.reset_index(drop=False, inplace=True)
        tmp_df = tmp_df.melt(id_vars="index", var_name="CellType", value_name="FDR")
        tmp_df = tmp_df.loc[tmp_df["FDR"] < alpha, :]
        # print(tmp_df)
        # genes = set()
        # for index in tmp_df["index"]:
        #     genes.add(index.split("_")[0])
        # print(genes)
        # print(len(genes))
        # exit()

        data = {}
        for ct in list(tmp_df["CellType"].unique()):
            data[ct] = set(tmp_df.loc[tmp_df["CellType"] == ct, "index"])
        return data

    @staticmethod
    def count(input_data):
        combinations = []
        cols = list(input_data.keys())
        for i in range(1, len(cols) + 1):
            combinations.extend(list(itertools.combinations(cols, i)))

        abbreviations = {"CellMapNNLS_Neuron": "neuro",
                         "CellMapNNLS_Oligodendrocyte": "oligo",
                         "CellMapNNLS_EndothelialCell": "endo",
                         "CellMapNNLS_Macrophage": "macro",
                         "CellMapNNLS_Astrocyte": "astro"}
        abbr_cols = []
        for col in cols:
            if col in abbreviations.keys():
                abbr_cols.append(abbreviations[col])
            else:
                abbr_cols.append(col)

        indices = []
        data = []
        for combination in combinations:
            index = []
            for col in cols:
                if col in combination:
                    index.append(True)
                else:
                    index.append(False)

            background = set()
            for key in cols:
                if key not in combination:
                    work_set = input_data[key].copy()
                    background.update(work_set)

            overlap = None
            for key in combination:
                work_set = input_data[key].copy()
                if overlap is None:
                    overlap = work_set
                else:
                    overlap = overlap.intersection(work_set)

            duplicate_set = overlap.intersection(background)
            length = len(overlap) - len(duplicate_set)

            indices.append(index)
            data.append(length)

        s = pd.Series(data,
                      index=pd.MultiIndex.from_tuples(indices, names=abbr_cols))
        s.name = "value"
        return s

    def print_arguments(self):
        print("Arguments:")
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Alpha: {}".format(self.alpha))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
