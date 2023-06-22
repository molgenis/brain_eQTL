#!/usr/bin/env python3

"""
File:         filter_interaction_eqtls.py
Created:      2020/09/25
Last Changed: 2022/02/10
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
import os

# Third party imports.
import pandas as pd
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Filter Interaction eQTLs"
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
        self.decon_path = getattr(arguments, 'decon')
        self.gene_info_path = getattr(arguments, 'gene_info')
        self.alpha = getattr(arguments, 'alpha')
        self.interest = getattr(arguments, 'interest')
        outfile = getattr(arguments, 'outfile')

        # Set variables.
        outdir = str(Path(__file__).parent.parent)
        self.outfile = None
        if outfile is not None:
            self.outfile = os.path.join(outdir, "{}_{}FDRFiltered.txt.gz".format(outfile, self.alpha))

    def create_argument_parser(self):
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
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix.")
        parser.add_argument("-g",
                            "--gene_info",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the gene info matrix. "
                                 "Default: None.")
        parser.add_argument("-a",
                            "--alpha",
                            type=float,
                            required=False,
                            default=0.05,
                            help="The significance cut-off. Default: 0.05.")
        parser.add_argument("-i",
                            "--interest",
                            nargs="+",
                            type=str,
                            required=False,
                            default=None,
                            help="The HGNCSymbols to print. Default: none.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            default=None,
                            help="The name of the output file (no extension)."
                                 " Default: None.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        decon_df = self.load_file(self.decon_path)

        gene_dict = None
        if self.gene_info_path is not None:
            gene_info_df = self.load_file(self.gene_info_path, index_col=None)
            gene_dict = dict(zip(gene_info_df["ArrayAddress"], gene_info_df["Symbol"]))

        decon_fdr_df = self.bh_correct(decon_df)
        filtered_decon_df, pre_shape, post_shape = self.filter(decon_fdr_df,
                                                               gene_dict,
                                                               self.alpha)

        self.print_data(filtered_decon_df,
                        pre_shape,
                        post_shape,
                        self.alpha,
                        self.interest)

        if self.outfile is not None:
            self.save_file(filtered_decon_df,
                           self.outfile)

    @staticmethod
    def bh_correct(pvalue_df):
        df = pvalue_df.copy()
        data = []
        indices = []
        for col in df.columns:
            if col.endswith("_pvalue"):
                data.append(multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
                indices.append(col.replace("_pvalue", ""))
        fdr_df = pd.DataFrame(data, index=indices, columns=df.index)

        return fdr_df.T

    @staticmethod
    def filter(df, trans_dict, alpha):
        print(df)
        df.reset_index(drop=False, inplace=True)
        df = df.melt(id_vars="index", var_name="CellType", value_name="FDR")
        df[['ProbeName', 'SNPName']] = df["index"].str.split("_", n=1, expand=True)
        print(df)

        if trans_dict is None:
            df["HGNCSymbol"] = None
        else:
            df["HGNCSymbol"] = df["ProbeName"].map(trans_dict)

        df = df[['ProbeName', 'HGNCSymbol', 'SNPName', 'CellType', 'FDR']]
        print(df)

        pre_shape = df.shape

        df = df.loc[df["FDR"] <= alpha, :]
        df = df.sort_values(by="FDR", ascending=True)
        df.reset_index(drop=True, inplace=True)
        print(df)

        post_shape = df.shape

        return df, pre_shape, post_shape

    @staticmethod
    def print_data(df, pre_shape, post_shape, alpha, interest):
        print("")
        print("{}/{} tests where significant (FDR <= {}) [{:.2f}%]".format(post_shape[0], pre_shape[0], alpha, (100 / pre_shape[0]) * post_shape[0]))
        print("")

        print("Top 10 hits:")
        print(df.iloc[0:10, :])
        print("")

        counts = df["CellType"].value_counts()
        print("CellType info:")
        for index, value in counts.iteritems():
            print("\t{}\tN: {}".format(index, value))
        print("\t----------------------")
        print("\tTotal\tN: {}".format(counts.sum()))
        print("")

        if interest is not None:
            for gene in interest:
                print("HGNCSymbol: {}".format(gene))
                print(df.loc[df["HGNCSymbol"] == gene, :])
                print("")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print(
            "\tSaved dataframe: {} "
            "with shape: {}".format(os.path.basename(outpath),
                                    df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Alpha: {}".format(self.alpha))
        print("  > Interest: {}".format(self.interest))
        print("  > Output file: {}".format(self.outfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
