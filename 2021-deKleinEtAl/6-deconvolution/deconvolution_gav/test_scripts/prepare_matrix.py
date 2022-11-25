#!/usr/bin/env python3

"""
File:         prepare_matrix.py
Created:      2020/10/20
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
__program__ = "Prepare Matrix"
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
        self.matrix_inpath = getattr(arguments, 'matrix')
        self.sample_dict_inpath = getattr(arguments, 'sample_dict')
        self.print_interval = 500
        self.write_interval = 1000

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
        self.outpath = os.path.join(outdir, os.path.basename(self.matrix_inpath))

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
        parser.add_argument("-m",
                            "--matrix",
                            type=str,
                            required=True,
                            help="The path to the matrix")
        parser.add_argument("-s",
                            "--sample_dict",
                            type=str,
                            required=True,
                            help="A matrix with two columns for translating "
                                 "sample identifiers (left = key, "
                                 "right = value).")

        return parser.parse_args()

    def validate(self):
        for filepath in [self.matrix_inpath, self.sample_dict_inpath]:
            if filepath is not None and not (os.path.exists(filepath) and os.path.isfile(filepath)):
                print("File {} does not exist".format(filepath))
                return False

        # Check if correct extension.
        if not self.matrix_inpath.endswith(".txt.gz"):
            print("Matrix input must be in .txt.gz format")
            return False

        return True

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        matrix_df, sample_order = self.load_data()

        print("### Step2 ###")
        order = self.check_overlap(list(matrix_df.columns.to_list()), sample_order)

        print("### Step3 ###")
        new_df = self.work(matrix_df, order, sample_order)

        print("### Step4 ###")
        self.save(new_df)

    def load_data(self):
        print("Loading matrix.")
        matrix_df = pd.read_csv(self.matrix_inpath,
                                sep="\t",
                                header=0,
                                index_col=0)
        print("\tLoaded dataframe: {} "
                    "with shape: {}".format(os.path.basename(self.matrix_inpath),
                                            matrix_df.shape))

        print("Loading sample_dict.")
        sample_dict_df = pd.read_csv(self.sample_dict_inpath,
                                     sep="\t",
                                     header=None)
        print("\tLoaded dataframe: {} "
                    "with shape: {}".format(os.path.basename(self.sample_dict_inpath),
                                            sample_dict_df.shape))

        sample_dict = dict(zip(sample_dict_df.iloc[:, 0],
                               sample_dict_df.iloc[:, 1]))
        sample_order = list(sample_dict_df.iloc[:, 1])

        matrix_df.columns = [sample_dict[x] if x in sample_dict else x for x in matrix_df.columns]

        return matrix_df, sample_order

    @staticmethod
    def check_overlap(matrix_df_columns, sample_order):
        order = []
        for sample in sample_order:
            if sample in matrix_df_columns:
                order.append(matrix_df_columns.index(sample))

        print("\t{}/{} [{:.2f}%] of columns are overlapping".format(len(order), len(matrix_df_columns), (100/len(matrix_df_columns)*len(order))))

        return order

    def work(self, df, order, sample_order):
        print("Filtering data.")
        subset = df.iloc[:, order]
        if not np.array_equal(subset.columns.to_numpy(), np.array(sample_order)):
            print("Error! Something went wrong in the sample order.")
            exit()
        return subset

    @staticmethod
    def write_buffer(filename, buffer):
        with gzip.open(filename, 'wb') as f:
            for line in buffer:
                f.write(line.encode())
        f.close()

    def save(self, new_df):
        new_df.to_csv(self.outpath,
                      sep="\t",
                      index=True,
                      header=True,
                      compression='gzip')

    def print_arguments(self):
        print("Arguments:")
        print("  > Matrix file: {}".format(self.matrix_inpath))
        print("  > Sample dict file: {}".format(self.sample_dict_inpath))
        print("  > Outpath {}".format(self.outpath))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
