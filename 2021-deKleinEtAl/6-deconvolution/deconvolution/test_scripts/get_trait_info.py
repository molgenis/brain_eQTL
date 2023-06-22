#!/usr/bin/env python3

"""
File:         get_trait_info.py
Created:      2020/11/17
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
import argparse
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Get Trait Info"
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
        # Default arguments.
        self.snp_to_gwasid_path = "/groups/umcg-biogen/tmp03/annotation/ieugwas/2020-05-03-allTopAssociations-wgwascatalog-wALS-MetaBrain2dot1IDs.txt.gz"
        self.gwasid_to_trait_path = "/groups/umcg-biogen/tmp03/annotation/ieugwas/2020-05-03-gwaslist-wgwascatalog-wALS.txt.gz"

        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.snps = getattr(arguments, 'snps')

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
        parser.add_argument("-s",
                            "--snps",
                            nargs="*",
                            type=str,
                            required=True,
                            help="The SNPs identifiers to look up.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading GWAS ID - Trait data frame.")
        gwas_trait_df = pd.read_csv(self.gwasid_to_trait_path, sep="\t", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.gwasid_to_trait_path),
                                      gwas_trait_df.shape))

        print("### Step2 ###")
        print("Loading SNP - GWAS ID data frame.")
        snp_gwas_df = pd.read_csv(self.snp_to_gwasid_path, sep="\t", header=0,
                                  index_col=False, low_memory=False)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.snp_to_gwasid_path),
                                      snp_gwas_df.shape))

        print("### Step2 ###")
        for snp in self.snps:
            print("Searching SNP: '{}'".format(snp))
            snp_gwas_subset = snp_gwas_df.loc[snp_gwas_df["RsID"] == snp, :]

            if snp_gwas_subset.shape[0] == 0:
                print("\tNo GWAS ID found.")
                continue

            trait_df = None
            for _, row in snp_gwas_subset.iterrows():
                gwas_trait_subset = gwas_trait_df.loc[gwas_trait_df["ID"] == row["ID"], :]
                if gwas_trait_subset.shape[0] == 0:
                    continue

                data = gwas_trait_subset.merge(row.to_frame().T, left_on="ID", right_on=["ID"])

                if trait_df is None:
                    trait_df = data
                else:
                    trait_df = pd.concat([trait_df, data], axis=0)

            if trait_df.shape[0] == 0:
                print("\tNo traits found.")
                continue

            print(trait_df)
            print("")


    def print_arguments(self):
        print("Arguments:")
        print("  > SNPs path: {}".format(self.snps))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
