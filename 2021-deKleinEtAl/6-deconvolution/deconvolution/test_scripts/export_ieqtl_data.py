#!/usr/bin/env python3

"""
File:         export_ieqtl_data.py
Created:      2021/02/26
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
from colour import Color
import argparse
import os

# Third party imports.
import pandas as pd
import seaborn as sns
import matplotlib
from statsmodels.stats import multitest
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Export ieQTL Data"
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
./export_ieqtl_data.py -eq ../../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ge ../../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -al ../../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/create_matrices/genotype_alleles.txt.gz -ex ../../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/create_matrices/expression_table.txt.gz -cc ../../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/perform_deconvolution/deconvolution_table.txt.gz -i CYP24A1 CLECL1 -n 600
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.geno_path = getattr(arguments, 'genotype')
        self.alleles_path = getattr(arguments, 'alleles')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cellcount%')
        self.interest = getattr(arguments, 'interest')
        self.nrows = getattr(arguments, 'nrows')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'export_ieqtl_data')

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
        parser.add_argument("-eq",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the eqtl matrix")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix")
        parser.add_argument("-al",
                            "--alleles",
                            type=str,
                            required=True,
                            help="The path to the alleles matrix")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix")
        parser.add_argument("-cc",
                            "--cellcount%",
                            type=str,
                            required=True,
                            help="The path to the cell count % matrix")
        parser.add_argument("-i",
                            "--interest",
                            nargs="+",
                            type=str,
                            required=False,
                            default=None,
                            help="The HGNCSymbols to plot. Default: none.")
        parser.add_argument("-n",
                            "--nrows",
                            type=int,
                            required=False,
                            default=None,
                            help="Cap the number of runs to load. "
                                 "Default: None.")

        return parser.parse_args()

    def start(self):
        print("Loading data.")
        eqtl_df = self.load_file(self.eqtl_path, index_col=None, nrows=self.nrows)
        geno_df = self.load_file(self.geno_path, nrows=self.nrows)
        alleles_df = self.load_file(self.alleles_path, nrows=self.nrows)
        expr_df = self.load_file(self.expr_path, nrows=self.nrows)
        cc_df = self.load_file(self.cc_path)

        # pre-process cc_df.
        cc_df.columns = [x.split("_")[1].lower() for x in cc_df.columns]

        print("Iterating over data.")
        for i, (index, row) in enumerate(eqtl_df.iterrows()):
            if self.interest is not None and row["HGNCName"] not in self.interest:
                continue

            # prepare output directory.
            outdir = os.path.join(self.outdir, row["HGNCName"])
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            info_s = row.copy()

            # Get the genotype / expression data.
            genotype = geno_df.iloc[i, :].T.to_frame()
            if genotype.columns != [info_s["SNPName"]]:
                print("\t\tGenotype file not in identical order as eQTL file.")
                exit()
            expression = expr_df.iloc[i, :].T.to_frame()
            if expression.columns != [info_s["ProbeName"]]:
                print("\t\tExpression file not in identical order as eQTL file.")
                exit()
            data = genotype.merge(expression, left_index=True, right_index=True)
            data.columns = ["genotype", "expression"]
            data = data.merge(cc_df, left_index=True, right_index=True)

            # Saving data.
            self.save_file(data, os.path.join(outdir, "data.txt.gz"), sep=",")

            # Get the allele data.
            alleles = alleles_df.iloc[i, :]
            if alleles.name != info_s["SNPName"]:
                print("\t\tAlleles file not in identical order as eQTL file.")
                exit()
            info_s["0_allele"] = alleles["Alleles"].split("/")[0]
            info_s["2_allele"] = alleles["Alleles"].split("/")[1]

            # Count the genotypes.
            counts = data["genotype"].round(0).value_counts()
            for x in [-1.0, 0.0, 1.0, 2.0]:
                if x in counts:
                    info_s.loc[x] = counts.loc[x]
                else:
                    info_s.loc[x] = 0

            # Saving info.
            self.save_file(info_s, os.path.join(outdir, "info.txt.gz"), sep=",")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
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
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Alleles path: {}".format(self.alleles_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell count % path: {}".format(self.cc_path))
        print("  > Interest: {}".format(self.interest))
        print("  > Nrows: {}".format(self.nrows))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
