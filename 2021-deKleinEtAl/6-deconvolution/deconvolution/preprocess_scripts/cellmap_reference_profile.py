#!/usr/bin/env python3

"""
File:         cellmap_reference_profile.py
Created:      2021/07/06
Last Changed: 2021/10/13
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
import sys
import os

# Third party imports.
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Cellmap Reference Profile"
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
./cellmap_reference_profile.py -d ../data/CellMap_brain_CNS7_log2CPM.txt -s ../data/CNS7_signature_genes.txt -o CellMap_brain_CNS7_avgCPM.txt
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.signature_genes_path = getattr(arguments, 'signature_genes')
        self.outfile = getattr(arguments, 'outfile')

        self.plotdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.plotdir):
            os.makedirs(self.plotdir)

        self.cell_type_dict = {
            "Oligodendrocytes": "Oligodendrocyte",
            "Endothelial": "EndothelialCell",
            "Pericyte": "Pericytes",
            "Astrocytes": "Astrocyte"
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
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the Cellmap reference profile")
        parser.add_argument("-s",
                            "--signature_genes",
                            type=str,
                            required=True,
                            help="The path to the signature genes list")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the output file")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")
        df = self.load_file(self.data_path, header=0, index_col=0)
        signature_genes_df = self.load_file(self.signature_genes_path, header=0, index_col=None)

        print("Selecting signature genes")
        overlap = set(signature_genes_df["gene"]).intersection(set(df.index))
        print("\toverlap = {} [{}/{}]".format(len(overlap), len(overlap), signature_genes_df.shape[0]))
        df = df.loc[signature_genes_df["gene"], :]

        print("Preprocessing header")
        info = []
        for col in df.columns:
            info.append(col.split("|"))
        info_df = pd.DataFrame(info, columns=["cell type", "dataset", "num", "iter"])
        print(info_df)

        print("Step 1: inversing log2 transform.")
        df = np.power(2, df)

        print("Step 2: Plotting input matrix")
        self.plot_clustermap(df=df, filename=self.outfile.split(".txt")[0] + "_full")

        print("Step 3: take average per cell type")
        avg_df = None
        for cell_type in info_df["cell type"].unique():
            mask = (info_df["cell type"] == cell_type).to_numpy()
            ct_avg_expr = df.loc[:, mask].mean(axis=1).to_frame()
            adj_cell_type = cell_type
            if adj_cell_type in self.cell_type_dict:
                adj_cell_type = self.cell_type_dict[cell_type]
            ct_avg_expr.columns = adj_cell_type

            if avg_df is None:
                avg_df = ct_avg_expr
            else:
                avg_df = pd.concat([avg_df, ct_avg_expr], axis=1)
        print(avg_df)

        print("Step 4: Plotting average matrix")
        self.plot_clustermap(df=avg_df, filename=self.outfile.split(".txt")[0] + "_average")

        print("Step 5: prepare output matrix")
        avg_df.insert(0, "GeneSymbol", avg_df.index)

        print("Saving file")
        self.save_file(avg_df,
                       outpath=os.path.join(os.path.dirname(self.data_path), self.outfile),
                       index=False)

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def plot_clustermap(self, df, filename):
        sys.setrecursionlimit(100000)

        sns.set(color_codes=True)
        g = sns.clustermap(df, cmap="Blues",
                           row_cluster=True, col_cluster=False,
                           yticklabels=False, xticklabels=False,
                           figsize=(12, 9))

        outpath = os.path.join(self.plotdir, "{}_clustermap.png".format(filename))
        plt.tight_layout()
        g.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(outpath))

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Data path: {}".format(self.data_path))
        print("  > Signature genes path: {}".format(self.signature_genes_path))
        print("  > Output filename: {}".format(self.outfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
