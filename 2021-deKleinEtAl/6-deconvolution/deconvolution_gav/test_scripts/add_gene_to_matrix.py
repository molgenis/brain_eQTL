#!/usr/bin/env python3

"""
File:         add_gene_to_matrix.py
Created:      2020/10/20
Last Changed: 2020/10/28
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
import gzip
import os

# Third party imports.
import pandas as pd
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Add Gene to Matrix"
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
        self.expression_inpath = getattr(arguments, 'expression')
        self.genes = getattr(arguments, 'genes')
        self.sample_dict_inpath = getattr(arguments, 'sample_dict')
        self.print_interval = 500

        # Validate input.
        if not self.validate():
            exit()

        # Define the output file.
        self.outpath = self.matrix_inpath.replace(".txt.gz", "_with{}Genes.txt.gz".format(len(self.genes)))

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
        parser.add_argument("-m",
                            "--matrix",
                            type=str,
                            required=True,
                            help="The path to the matrix.")
        parser.add_argument("-e",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
        parser.add_argument("-g",
                            "--genes",
                            nargs="*",
                            type=str,
                            required=True,
                            help="The name(s) of the gene to add.")
        parser.add_argument("-s",
                            "--sample_dict",
                            type=str,
                            default=None,
                            required=False,
                            help="A matrix with two columns for translating "
                                 "sample identifiers (left = key, "
                                 "right = value. Default: None.")

        return parser.parse_args()

    def validate(self):
        for filepath in [self.matrix_inpath, self.expression_inpath]:
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
        matrix_df, sample_dict = self.load_data()
        print(matrix_df)
        print(sample_dict)

        print("### Step2 ###")
        new_matrix_df = self.work(sample_dict)
        print(new_matrix_df)

        print("### Step2 ###")
        print("Combine and safe.")
        combined_df = matrix_df.T.merge(new_matrix_df.T, left_index=True, right_index=True).T
        print(combined_df)
        combined_df.to_csv(self.outpath,
                           sep="\t",
                           index=True,
                           header=True,
                           compression='gzip')

    def load_data(self):
        print("Loading matrix.")
        matrix_df = pd.read_csv(self.matrix_inpath,
                                sep="\t",
                                header=0,
                                index_col=0)

        sample_dict = {}
        if self.sample_dict_inpath is not None:
            print("Loading sample_dict.")
            sample_dict_df = pd.read_csv(self.sample_dict_inpath,
                                         sep="\t",
                                         header=0)
            print("\tLoaded dataframe: {} "
                  "with shape: {}".format(os.path.basename(self.sample_dict_inpath),
                                          sample_dict_df.shape))

            sample_dict = dict(zip(sample_dict_df.iloc[:, 0],
                                   sample_dict_df.iloc[:, 1]))

        return matrix_df, sample_dict

    def work(self, trans_dict):
        new_columns = []
        new_indices = []
        new_data = []

        search_list = set(self.genes)
        print("Searching for {}".format(", ".join(search_list)))
        with gzip.open(self.expression_inpath, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % self.print_interval == 0):
                    print("\tprocessed {} lines".format(i))
                if len(search_list) == 0:
                    break

                splitted_line = line.decode().strip('\n').split('\t')
                index = splitted_line[0]
                data = splitted_line[1:]

                if i == 0:
                    new_columns = [trans_dict[x] if x in trans_dict else x for x in data]
                else:
                    if index in search_list:
                        print("\t\tFound {}".format(index))
                        new_indices.append(index)
                        new_data.append(data)

                        search_list.remove(index)
        f.close()

        return pd.DataFrame(new_data, index=new_indices, columns=new_columns)

    def print_arguments(self):
        print("Arguments:")
        print("  > Matrix file: {}".format(self.matrix_inpath))
        print("  > Expression file: {}".format(self.expression_inpath))
        print("  > Sample dict file: {}".format(self.sample_dict_inpath))
        print("  > Outpath {}".format(self.outpath))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
