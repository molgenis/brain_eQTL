#!/usr/bin/env python3

"""
File:         get_sn_avg_gene_expr.py
Created:      2020/11/30
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
__program__ = "Get SN Average Gene Expression"
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
./get_sn_avg_gene_expr.py -ex /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/ -cc /groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/2020-10-22-MVochteloo-Copy/cell_counts.txt -g /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/multiple_sclerosis_patsopoulos_harm_jan_enrichtments_exHla.txt
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.expr_path = getattr(arguments, 'expression')
        self.suffix = getattr(arguments, 'suffix')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.genes_path = getattr(arguments, 'genes')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "get_sn_avg_gene_expr")
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
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the folder containing the "
                                 "expression matrices")
        parser.add_argument("-s",
                            "--suffix",
                            type=str,
                            required=False,
                            default="_expression.txt",
                            help="The expressio matrix file suffix. Default:"
                                 " <cell_type>_expression.txt")
        parser.add_argument("-cc",
                            "--cell_counts",
                            type=str,
                            required=True,
                            help="The path to the cell counts matrix")
        parser.add_argument("-g",
                            "--genes",
                            type=str,
                            required=True,
                            help="The path to a list of genes")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data")
        cc_df = self.load_file(self.cc_path)
        cc_df.index = [x.upper() for x in cc_df.index]
        print(cc_df)

        genes = list(self.load_file(self.genes_path, header=None, index_col=None).iloc[:, 0].values)
        print(genes)

        print("### Step2 ###")
        print("Filter data")

        for filepath in glob.glob(os.path.join(self.expr_path, "*" + self.suffix)):
            cell_type = os.path.basename(filepath).replace(self.suffix, "").upper()
            if cell_type in cc_df.index:
                ct_expr_df = self.load_file(filepath)
                ct_expr_df.index = [x.split(".")[0] for x in ct_expr_df.index]

                row_overlap = set(ct_expr_df.index).intersection(set(genes))
                col_overlap = set(ct_expr_df.columns).intersection(set(cc_df.columns))
                ct_expr_df = ct_expr_df.loc[row_overlap, col_overlap]
                print("\t Cell type: {} has {} / {} genes and {} / {} samples.".format(cell_type, len(row_overlap), len(genes), len(col_overlap), len(cc_df.columns)))

                row_missing = set(genes).symmetric_difference(set(ct_expr_df.index))
                if len(row_missing) > 0:
                    missing_df = pd.DataFrame(np.nan, index=row_missing, columns=ct_expr_df.columns)
                    ct_expr_df = pd.concat([ct_expr_df, missing_df])
                    ct_expr_df = ct_expr_df.loc[genes, :]

                col_missing = set(cc_df.columns).symmetric_difference(set(ct_expr_df.columns))
                if len(col_missing) > 0:
                    missing_df = pd.DataFrame(np.nan, index=genes, columns=col_missing)
                    ct_expr_df = ct_expr_df.merge(missing_df, left_index=True, right_index=True)
                    ct_expr_df = ct_expr_df.loc[:, cc_df.columns]

                print(ct_expr_df.shape)




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
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Suffix: {}".format(self.suffix))
        print("  > Genes path: {}".format(self.genes_path))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
