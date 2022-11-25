#!/usr/bin/env python3

"""
File:         find_top_interacting_covariate.py
Created:      2020/10/26
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
__program__ = "Find Top Interacting Covariate"
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
        self.decon_inpath = getattr(arguments, 'decon')
        self.gene_info_inpath = getattr(arguments, 'gene_info')
        self.outfile_prefix = getattr(arguments, 'outfile_prefix')
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
        self.outpath = os.path.join(outdir, self.outfile_prefix + "_" + os.path.basename(self.decon_inpath).replace(".txt.gz", "_chi2sum.txt.gz"))

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
                            help="The path to the decon matrix.")
        parser.add_argument("-gi",
                            "--gene_info",
                            type=str,
                            default=None,
                            help="The path to the gene info matrix.")
        parser.add_argument("-p",
                            "--outfile_prefix",
                            type=str,
                            default="",
                            help="The outputfile prefix. Default: ''.")

        return parser.parse_args()

    def validate(self):
        if not (os.path.exists(self.decon_inpath) and os.path.isfile(self.decon_inpath)):
            print("File {} does not exist".format(self.decon_inpath))
            return False

        # Check if correct extension.
        if not self.decon_inpath.endswith(".txt.gz"):
            print("Matrix input must be in .txt.gz format")
            return False

        return True

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        gene_info_df = self.load_data()

        print("### Step2 ###")
        chi2sum_df = self.work()
        print(chi2sum_df)

        print("### Step3 ###")
        print("Combine")
        combined_df = chi2sum_df.merge(gene_info_df, left_on="index", right_on="ArrayAddress")
        combined_df.sort_values(by="chi2sum", ascending=False, inplace=True)
        print(combined_df)

        print("Top covariates:")
        print(combined_df.head(5))

        print("Bottom covariates:")
        print(combined_df.tail(5))

        print("### Step4 ###")
        print("Safe")
        combined_df.to_csv(self.outpath,
                           sep="\t",
                           index=True,
                           header=True,
                           compression='gzip')

    def load_data(self):
        gene_info_df = None
        if self.gene_info_inpath is not None:
            print("Loading gene info matrix.")
            gene_info_df = pd.read_csv(self.gene_info_inpath,
                                       sep="\t",
                                       header=0)
            print("\tLoaded dataframe: {} "
                  "with shape: {}".format(os.path.basename(self.gene_info_inpath),
                                          gene_info_df.shape))

        return gene_info_df

    def work(self):
        new_data = []
        sign_steps = [0.01, 0.02, 0.03, 0.04, 0.05]

        print("Start parsing file.")
        with gzip.open(self.decon_inpath, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % self.print_interval == 0):
                    print("\tprocessed {} lines".format(i))

                if i > 0:
                    splitted_line = line.decode().strip('\n').split('\t')
                    index = splitted_line[0]
                    data = np.array(splitted_line[1:], dtype=float)
                    counts = []
                    for step in sign_steps:
                        counts.append(np.sum(data < step))
                    new_data.append([index] + counts + [np.sum(np.power(data, 2))])

        f.close()

        return pd.DataFrame(new_data, columns=["index"] + ["p<{}".format(x) for x in sign_steps] + ["chi2sum"])

    def print_arguments(self):
        print("Arguments:")
        print("  > Decon file: {}".format(self.decon_inpath))
        print("  > Gene info file: {}".format(self.gene_info_inpath))
        print("  > Outpath {}".format(self.outpath))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
