#!/usr/bin/env python3

"""
File:         prepare_eqtl_file.py
Created:      2021/12/07
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
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Prepare eQTL file"
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
./prepare_eqtl_file.py -h

./prepare_eqtl_file.py -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/80PCs-1000perm-IterationsMerged.txt -p 80PCs-1000perm-

"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.n_iterations = getattr(arguments, 'n_iterations')
        self.prefix = getattr(arguments, 'prefix')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'prepare_eqtl_file')
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
        parser.add_argument("-eq",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the replication eqtl matrix.")
        parser.add_argument("-ni",
                            "--n_iterations",
                            type=int,
                            default=4,
                            help="The number of eQTL iterations to include. "
                                 "Default: 4.")
        parser.add_argument("-p",
                            "--prefix",
                            type=str,
                            required=True,
                            help="Prefix for the output file.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        eqtl_df = self.load_file(self.eqtl_path, header=0, index_col=None)
        print(eqtl_df)

        print("Preprocessing data")
        # Translating column names.
        trans_dict = {
            "SNP": "SNPName",
            "SNPPos": "SNPChrPos",
            "Gene": "ProbeName",
            "GeneChr": "ProbeChr",
            "GenePos": "ProbeCenterChrPos",
            "SNPAlleles": "SNPType",
            "SNPEffectAllele": "AlleleAssessed",
            "QTLRank": "Iteration"}
        eqtl_df.columns = [trans_dict[x] if x in trans_dict else x for x in eqtl_df.columns]

        # Removing rows with a too high iteration.
        eqtl_df = eqtl_df.loc[eqtl_df["Iteration"] <= self.n_iterations, :]
        print(eqtl_df["Iteration"].value_counts())

        # Sorting file.
        eqtl_df.sort_values(by=["Iteration", "MetaP"], inplace=True)
        print(eqtl_df)

        print("Parsing eQTL file.")
        mask = np.zeros(eqtl_df.shape[0], dtype=bool)
        found_genes = {i: set() for i in range(1, self.n_iterations+1)}
        for i, (_, row) in enumerate(eqtl_df.iterrows()):
            if (i == 0) or (i % 5000 == 0):
                print("\tprocessed {} lines".format(i))

            if row["ProbeName"] not in found_genes[row["Iteration"]]:
                mask[i] = True
                found_genes[row["Iteration"]].add(row["ProbeName"])
        print(np.size(mask) - np.sum(mask))

        top_eqtl_df = eqtl_df.loc[mask, :]
        print(top_eqtl_df["Iteration"].value_counts())

        print("Saving file.")
        print(top_eqtl_df)
        self.save_file(df=top_eqtl_df, outpath=os.path.join(self.outdir, "{}eQTLProbesFDR0.05-ProbeLevel.txt.gz".format(self.prefix)), index=False)

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
        print("  > eQTL: {}".format(self.eqtl_path))
        print("  > Prefix: {}".format(self.prefix))
        print("  > N-iterations: {}".format(self.n_iterations))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
