#!/usr/bin/env python3

"""
File:         get_gene_expression.py
Created:      2021/05/19
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

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Get Gene Expression"
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
./export_gene_expression.py -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/expression_noSCTransform

"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.indir = getattr(arguments, 'input_directory')
        self.gene_info_path = "/groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz"

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "get_gene_expression")
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
        parser.add_argument("-i",
                            "--input_directory",
                            type=str,
                            required=True,
                            help="The path to the folder containing the "
                                 "expression matrices")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Load gene info.")
        gene_info_df = self.load_file(self.gene_info_path)
        gene_dict = dict(zip(gene_info_df["ArrayAddress"], gene_info_df["Symbol"]))
        del gene_info_df

        data = {}
        order = []
        for filepath in glob.glob(os.path.join(self.indir, "*")):
            cell_type = os.path.basename(filepath).split("_")[0].upper()
            order.append(cell_type)

            df = self.load_file(filepath)
            data[cell_type] = df.mean(axis=1)

        merged_df = pd.DataFrame(data)
        merged_df.index.name = "ENSG"
        merged_df.reset_index(inplace=True, drop=False)
        merged_df["HGNC"] = merged_df["ENSG"].map(gene_dict)
        merged_df = merged_df.loc[:, ["ENSG", "HGNC"] + order]
        print(merged_df)
        exit()

        # % van cellen dat een gen tot expressie brengt
        # per individu wat de gemiddelde expressie is per gen gegeven dat het genormaliseerd is

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
