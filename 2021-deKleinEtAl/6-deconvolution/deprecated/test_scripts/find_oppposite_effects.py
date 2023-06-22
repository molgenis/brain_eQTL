#!/usr/bin/env python3

"""
File:         find_opposite_effects.py
Created:      2020/10/22
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
import argparse
import os

# Third party imports.
import pandas as pd
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Find Opposite Effects"
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
        self.decon1_path = getattr(arguments, 'decon1')
        self.name1 = getattr(arguments, 'name1')
        self.decon2_path = getattr(arguments, 'decon2')
        self.name2 = getattr(arguments, 'name2')
        self.alpha = getattr(arguments, 'alpha')

        # Check input.
        if not self.validate():
            exit()

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
        parser.add_argument("-d1",
                            "--decon1",
                            type=str,
                            required=True,
                            help="The path to the first deconvolution matrix")
        parser.add_argument("-n1",
                            "--name1",
                            type=str,
                            default="one",
                            help="The name for the first deconvolution matrix")
        parser.add_argument("-d2",
                            "--decon2",
                            type=str,
                            required=True,
                            help="The path to the second deconvolution matrix")
        parser.add_argument("-n2",
                            "--name2",
                            type=str,
                            default="two",
                            help="The name for the second deconvolution matrix")
        parser.add_argument("-a",
                            "--alpha",
                            type=float,
                            default=0.05,
                            help="The signifiance cutoff")

        return parser.parse_args()

    def validate(self):
        for filepath in [self.decon1_path, self.decon2_path]:
            if not os.path.exists(filepath) or not os.path.isfile(filepath):
                return False
        return True

    def start(self):
        self.print_arguments()

        print("### Loading data ###")
        d1_df, ct1, d2_df, ct2 = self.load_data()

        print("### Find overlapping celltypes ###")
        overlap_ct = ct1.intersection(ct2)

        print("### Add FDR ###")
        d1_df = self.add_bh_fdr(d1_df)
        d2_df = self.add_bh_fdr(d2_df)

        print("### Loop over each celltype ###")
        data = []
        for celltype in overlap_ct:
            print("Celltype: {}".format(celltype))

            # filter.
            subset1 = d1_df.loc[d1_df["{}_FDR".format(celltype)] <= self.alpha, [x for x in d1_df.columns if celltype in x]].copy()
            subset2 = d2_df.loc[d2_df["{}_FDR".format(celltype)] <= self.alpha, [x for x in d2_df.columns if celltype in x]].copy()
            print("\tN-signif. {}: {}".format(self.name1, subset1.shape[0]))
            print("\tN-signif. {}: {}".format(self.name2, subset2.shape[0]))

            subset1, subset2 = self.filter(subset1, subset2)

            beta_cols = [x for x in subset1.columns if "Beta" in x]
            n = 0
            for i in range(subset1.shape[0]):
                for col in beta_cols:
                    beta1 = subset1.loc[subset1.index[i], col]
                    beta2 = subset2.loc[subset2.index[i], col]
                    if beta1 * beta2 < 0:
                        n += 1
                        index = subset1.index[i].split(":")
                        data.append([index[0],
                                     ":".join(index[1:]),
                                     celltype,
                                     beta1,
                                     beta2,
                                     subset1.loc[subset1.index[i], "{}_FDR".format(celltype)],
                                     subset2.loc[subset2.index[i], "{}_FDR".format(celltype)]])
            print("Opposite betas: {}".format(n))

        opposite_df = pd.DataFrame(data,
                                   columns=["ProbeName",
                                            "SNPName",
                                            "CellType",
                                            "{}_Beta".format(self.name1),
                                            "{}_Beta".format(self.name2),
                                            "{}_FDR".format(self.name1),
                                            "{}_FDR".format(self.name2)])
        opposite_df.sort_values(by="{}_FDR".format(self.name1), inplace=True)
        print(opposite_df)

    @staticmethod
    def add_bh_fdr(df):
        for col in df.columns:
            if col.endswith("_pvalue"):
                df.loc[:, col.replace("_pvalue", "_FDR")] = multitest.multipletests(df.loc[:, col], method='fdr_bh')[1]

        return df

    def load_data(self):
        print("Loading decon matrix 1.")
        d1_df = pd.read_csv(self.decon1_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon1_path),
                                      d1_df.shape))
        ct1 = set(["_".join(x.split("_")[:-1]) for x in d1_df.columns if x.endswith("_pvalue")])

        print("Loading decon matrix 2.")
        d2_df = pd.read_csv(self.decon2_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon2_path),
                                      d2_df.shape))
        d2_df.columns = [x.replace("CellMapNNLS_", "") for x in d2_df.columns]
        ct2 = set(["_".join(x.split("_")[:-1]) for x in d2_df.columns if x.endswith("_pvalue")])

        return d1_df, ct1, d2_df, ct2

    @staticmethod
    def filter(d1_df, d2_df):
        overlap = set(d1_df.index).intersection(set(d2_df.index))
        print("\t{} of eQTLs are overlapping".format(len(overlap)))

        return d1_df.loc[overlap, :], d2_df.loc[overlap, :]

    def print_arguments(self):
        print("Arguments:")
        print("  > Deconvolution 1 path: {}".format(self.decon1_path))
        print("  > Deconvolution 1 name: {}".format(self.name1))
        print("  > Deconvolution 2 path: {}".format(self.decon2_path))
        print("  > Deconvolution 2 name: {}".format(self.name2))
        print("  > Alpha: {}".format(self.alpha))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
