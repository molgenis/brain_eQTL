#!/usr/bin/env python3

"""
File:         correlate_genes_with_cell_fractions.py
Created:      2022/05/03
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
import math
import json
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Correlate Genes with Cell Fractions"
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
./correlate_genes_with_cell_fractions.py -h

./correlate_genes_with_cell_fractions.py \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/expression_table.txt.gz \
    -gi /groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz \
    -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -o 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected \
    -sn no_ENA_0_PCs
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.expr_path = getattr(arguments, 'expression')
        self.gene_info_path = getattr(arguments, 'gene_info')
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.out_filename = getattr(arguments, 'outfile')
        self.sheet_name = getattr(arguments, 'sheet_name')

        # Set variables.
        base_dir = str(os.path.dirname(os.path.abspath(__file__)))
        self.outdir = os.path.join(base_dir, 'correlate_genes_with_cell_fractions')
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
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
        parser.add_argument("-gi",
                            "--gene_info",
                            type=str,
                            required=True,
                            help="The path to the gene info matrix.")
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fracttions matrix")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the sample-dataset link matrix.")
        parser.add_argument("-avge",
                            "--average_gene_expression",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the average gene expression "
                                 "matrix.")
        parser.add_argument("-p",
                            "--palette",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to a json file with the"
                                 "dataset to color combinations.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            default="output",
                            help="The name of the outfile. Default: output.")
        parser.add_argument("-sn",
                            "--sheet_name",
                            type=str,
                            default="Sheet1",
                            help="The name of the excel sheet. "
                                 "Default: Sheet1.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # Load data.
        print("Loading data.")
        expr_df = self.load_file(self.expr_path, header=0, index_col=0)
        expr_df = expr_df.groupby(expr_df.index).first()
        gene_info_df = self.load_file(self.gene_info_path, header=0, index_col=None)
        gene_dict = dict(zip(gene_info_df["ArrayAddress"], gene_info_df["Symbol"]))
        del gene_info_df
        cf_df = self.load_file(self.cf_path, header=0, index_col=0)
        std_df = self.load_file(self.std_path, header=0, index_col=None)

        print("Pre-processing data.")
        if cf_df.shape[0] > cf_df.shape[1]:
            cf_df = cf_df.T

        # Make sure order is the same.
        samples = std_df.iloc[:, 0].tolist()
        expr_df = expr_df.loc[:, samples]
        cf_df = cf_df.loc[:, samples]

        print("Correlating.")
        corr_df, pvalue_df = self.correlate(df1=expr_df,
                                            df2=cf_df)
        corr_df.insert(0, "ProbeName", corr_df.index)
        corr_df.insert(1, 'HGNCName', corr_df["ProbeName"].map(gene_dict))
        print(corr_df)

        print("Saving output")
        self.save_file(df=corr_df,
                       outpath=os.path.join(self.outdir, "{}.txt.gz".format(self.out_filename)),
                       index=False)
        self.save_file(df=corr_df,
                       outpath=os.path.join(self.outdir, "{}.xlsx".format(self.out_filename)),
                       index=False,
                       sheet_name=self.sheet_name.replace("_", " "))

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
    def correlate(df1, df2):
        all_coefficients = []
        all_pvalues = []
        for _, row1 in df1.iterrows():
            coefficients = []
            pvalues = []
            for _, row2 in df2.iterrows():
                coef, p = stats.spearmanr(row1, row2)
                coefficients.append(coef)
                pvalues.append(p)
            all_coefficients.append(coefficients)
            all_pvalues.append(pvalues)

        corr_coef_df = pd.DataFrame(all_coefficients,
                                    index=df1.index,
                                    columns=df2.index)

        pvalue_df = pd.DataFrame(all_pvalues,
                                 index=df1.index,
                                 columns=df2.index)

        return corr_coef_df, pvalue_df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t", na_rep="NA",
                  sheet_name="Sheet1"):
        if outpath.endswith('xlsx'):
            df.to_excel(outpath,
                        sheet_name=sheet_name,
                        na_rep=na_rep,
                        header=header,
                        index=index)
        else:
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
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Gene info: {}".format(self.gene_info_path))
        print("  > Cell fraction path: {}".format(self.cf_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Output filename: {}".format(self.out_filename))
        print("  > Sheet name: {}".format(self.sheet_name))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
