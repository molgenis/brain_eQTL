#!/usr/bin/env python3

"""
File:         rank_gene_expression.py
Created:      2020/10/28
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

# Local application imports.

# Metadata
__program__ = "Rank Gene Expresion"
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
        self.expr_inpath = getattr(arguments, 'expression')
        self.filter_path = getattr(arguments, 'filter')
        self.interest = getattr(arguments, 'interest')
        self.print_interval = 500

        # Validate input.
        if not self.validate():
            exit()

        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Prepare an output directory.
        outdir = os.path.join(current_dir, 'output')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # Define the output file.
        suffix = ""
        if self.filter_path is not None:
            suffix = "_filtered"
        self.outpath = os.path.join(outdir, os.path.basename(self.expr_inpath).replace(".txt.gz", "_ranks{}.txt.gz".format(suffix)))

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
        parser.add_argument("-e",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the decon matrix.")
        parser.add_argument("-f",
                            "--filter",
                            type=str,
                            default=None,
                            help="The path to the filter matrix. " \
                                 "Default: None.")
        parser.add_argument("-i",
                            "--interest",
                            nargs="*",
                            type=str,
                            default=None,
                            help="The indices of interest. Default: None.")

        return parser.parse_args()

    def validate(self):
        if not (os.path.exists(self.expr_inpath) and os.path.isfile(self.expr_inpath)):
            print("File {} does not exist".format(self.expr_inpath))
            return False

        # Check if correct extension.
        if not self.expr_inpath.endswith(".txt.gz"):
            print("Matrix input must be in .txt.gz format")
            return False

        return True

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        filter = None
        if self.filter_path is not None:
            print("\tLoading filter data frame.")
            filter_df = pd.read_csv(self.filter_path,
                                    sep="\t",
                                    header=None,
                                    index_col=None)

            filter = set(filter_df.iloc[:, 0].to_list())
            print("Filter length: {}".format(len(filter)))

        print("### Step2 ###")
        df = self.work(filter)

        print("### Step3 ###")
        print("Sort")
        df.sort_values(by="sum", ascending=False, inplace=True)
        df["rank"] = df.loc[:, "sum"].rank(ascending=False)

        print("Top genes:")
        print(df.head(5).round(4))

        print("Bottom genes:")
        print(df.tail(5).round(4))

        print("### Step4 ###")
        if self.interest is not None:
            for interest in self.interest:
                if interest in df["index"]:
                    print(df.loc[interest, :].round(4))
        else:
            print("No intersts.")

        print("### Step5 ###")
        print("Safe")
        df.to_csv(self.outpath,
                  sep="\t",
                  index=True,
                  header=True,
                  compression='gzip')

    def work(self, filter_list):
        new_data = []

        print("Start parsing file.")
        with gzip.open(self.expr_inpath, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % self.print_interval == 0):
                    print("\tprocessed {} lines".format(i))

                if i > 0:
                    splitted_line = line.decode().strip('\n').split('\t')
                    index = splitted_line[0]

                    if filter_list is None or index in filter_list:
                        data = np.array(splitted_line[1:], dtype=float)
                        new_data.append([index, np.min(data), np.max(data), np.mean(data), np.median(data), np.std(data), np.sum(data)])
        f.close()

        return pd.DataFrame(new_data, columns=["index", "min", "max", "mean", "median", "sd", "sum"])

    def print_arguments(self):
        print("Arguments:")
        print("  > Expression path: {}".format(self.expr_inpath))
        print("  > Filter path: {}".format(self.filter_path))
        print("  > Interest: {}".format(self.interest))
        print("  > Outpath {}".format(self.outpath))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
