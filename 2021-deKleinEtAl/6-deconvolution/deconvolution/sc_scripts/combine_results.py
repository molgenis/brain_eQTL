#!/usr/bin/env python3

"""
File:         combine_results.py
Created:      2020/11/04
Last Changed: 2020/12/01
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
__program__ = "Combine Results"
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
        self.eqtl_type = getattr(arguments, 'eqtl_type')

        # Set the other arguments.
        self.infolder = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/{}_100Perm".format(self.eqtl_type)
        self.cell_types = ["AST", "END", "EX", "IN", "MIC", "OLI", "OPC", "PER"]
        self.result_files = ["eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                             "eQTLSNPsFDR0.05-ProbeLevel.txt.gz",
                             "eQTLs.txt.gz",
                             "eQTLsFDR-ProbeLevel.txt.gz",
                             "eQTLsFDR0.05-ProbeLevel.txt.gz"]
        self.outdir = os.path.join(self.infolder, "COMBINED")
        self.interest = ["ENSG00000140873.16", "ENSG00000019186.10", "ENSG00000184293.7"]

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
                            help="show program's version number and exit")
        parser.add_argument("-e",
                            "--eqtl_type",
                            type=str,
                            choices=["cis", "trans"],
                            required=True,
                            help="The type of eQTLs to combine.")

        return parser.parse_args()

    def start(self):
        for filename in self.result_files:
            print("Combining {}".format(filename))

            combined_df = None
            for cell_type in self.cell_types:

                inpath = os.path.join(self.infolder, cell_type, filename)
                df = pd.read_csv(inpath, sep="\t")
                if df.shape[0] > 0:
                    print("\tCell type {} data has shape: {}".format(cell_type, df.shape))
                else:
                    print("\tCell type {} has no results".format(cell_type))
                    continue

                df["CellType"] = cell_type

                if combined_df is None:
                    combined_df = df
                else:
                    combined_df = pd.concat([combined_df, df],
                                            axis=0,
                                            ignore_index=True)

            if combined_df is not None:
                combined_df.loc[:, "ProbeName"] = combined_df.loc[:, "ProbeName"].astype(str)

                if self.interest:
                    print("Searching for results of interest:")
                    interest_frames = {x: None for x in self.interest}
                    for _, row in combined_df.iterrows():
                        for interest in self.interest:
                            if str(row["ProbeName"]) == interest:
                                if interest_frames[interest] is None:
                                    interest_frames[interest] = row.to_frame()
                                else:
                                    interest_frames[interest] = interest_frames[interest].merge(row.to_frame(),
                                                                                                left_index=True,
                                                                                                right_index=True)

                    for interest, interest_df in interest_frames.items():
                        if interest_df is not None:
                            interest_df = interest_df.T
                            interest_df.sort_values(by=["FDR", "PValue"], ascending=True, inplace=True)
                            print(interest_df)
                        else:
                            print("\t{} not found.".format(interest))

                combined_df.sort_values(by="FDR", ascending=True, inplace=True)
                combined_df.to_csv(os.path.join(self.outdir, filename),
                                   sep="\t",
                                   index=False,
                                   header=True,
                                   compression='gzip')


if __name__ == '__main__':
    m = main()
    m.start()
