#!/usr/bin/env python3

"""
File:         combine_star_output.py
Created:      2021/12/01
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
import re

# Third party imports.
import numpy as np
import pandas as pd


# Local application imports.

# Metadata
__program__ = "Combine STAR output"
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
# Calculon
./combine_star_output.py -d /groups/umcg-biogen/tmp04/input/rawdata/psychEncode/pipelines/no_patch_chromosomes/ -o 2021-12-01-MetaBrain
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        self.outdir = os.path.join(Path().resolve(), 'combine_star_output', outfolder)
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
                            help="show program's version number and exit.")
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the data matrix.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")

        datasets = [os.path.basename(fpath) for fpath in glob.glob(os.path.join(self.data_path, "*")) if os.path.isdir(fpath)]
        print("\tdatasets: {}".format(", ".join(datasets)))

        print("Loading data.")
        for dataset in datasets:
            unstranded_counts_df_list = []
            first_strand_counts_df_list = []
            second_strand_counts_df_list = []
            samples = []

            fpaths = glob.glob(os.path.join(self.data_path, dataset, "results", "star", "*.ReadsPerGene.out.tab.gz"))
            for fpath in fpaths:
                sample_id = os.path.basename(fpath).replace(".ReadsPerGene.out.tab.gz", "")

                df = self.load_file(fpath, header=None, index_col=0)
                unstranded_counts_df_list.append(df.iloc[4:, [0]])
                first_strand_counts_df_list.append(df.iloc[4:, [1]])
                second_strand_counts_df_list.append(df.iloc[4:, [2]])
                samples.append(sample_id)

            unstranded_counts_df = pd.concat(unstranded_counts_df_list, axis=1)
            unstranded_counts_df.columns = samples
            unstranded_counts_df.index.name = None
            print(unstranded_counts_df)
            self.save_file(df=unstranded_counts_df, outpath=os.path.join(self.outdir, "{}.unstrandedCounts.txt.gz".format(dataset)))

            first_strand_counts_df = pd.concat(first_strand_counts_df_list, axis=1)
            first_strand_counts_df.columns = samples
            first_strand_counts_df.index.name = None
            print(first_strand_counts_df)
            self.save_file(df=first_strand_counts_df, outpath=os.path.join(self.outdir, "{}.firstStrandCounts.txt.gz".format(dataset)))

            second_strand_counts_df = pd.concat(second_strand_counts_df_list, axis=1)
            second_strand_counts_df.columns = samples
            second_strand_counts_df.index.name = None
            print(second_strand_counts_df)
            self.save_file(df=second_strand_counts_df, outpath=os.path.join(self.outdir, "{}.secondStandCounts.txt.gz".format(dataset)))

            del unstranded_counts_df, first_strand_counts_df, second_strand_counts_df, unstranded_counts_df_list, first_strand_counts_df_list, second_strand_counts_df_list

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

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
        print("  > Gene info path: {}".format(self.data_path))
        print("  > Output directory {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
