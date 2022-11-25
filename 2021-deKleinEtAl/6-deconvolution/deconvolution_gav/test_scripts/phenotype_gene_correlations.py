#!/usr/bin/env python3

"""
File:         phenotype_gene_correlations.py
Created:      2021/02/02
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
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Phenotype Gene Correlations"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
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
./phenotype_gene_correlations.py -ph /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-03-09.brain.phenotypes.txt -ge ../data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.pkl -c RIN PMI
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.pheno_path = getattr(arguments, 'phenotype')
        self.ge_path = getattr(arguments, 'gene_expression')
        self.columns = getattr(arguments, 'columns')
        self.expr_id = getattr(arguments, 'expr_id')
        self.outfile_prefix = getattr(arguments, 'outfile_prefix')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "phenotype_gene_correlations")

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
        parser.add_argument("-ph",
                            "--phenotype",
                            type=str,
                            required=True,
                            help="The path to the phenotype matrix")
        parser.add_argument("-ge",
                            "--gene_expression",
                            type=str,
                            required=True,
                            help="The path to the gene expression matrix")
        parser.add_argument("-eid",
                            "--expr_id",
                            type=str,
                            required=False,
                            default="rnaseq_id",
                            help="The ID column for matchging the expression"
                                 " to the phenotype. Default: rnaseq_id.")
        parser.add_argument("-c",
                            "--columns",
                            nargs="*",
                            type=str,
                            required=False,
                            default=None,
                            help="Which columns to test for.")
        parser.add_argument("-op",
                            "--outfile_prefix",
                            type=str,
                            required=False,
                            default="",
                            help="A prefix for the output files.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading phenotype matrix.")
        pheno_df = self.load_file(self.pheno_path, index_col=None, low_memory=False)
        pheno_df.set_index(self.expr_id, inplace=True)

        print("Loading gene expression matrix.")
        ge_df = self.load_file(self.ge_path).T

        print("Pre-process")
        filter = []
        for column in pheno_df.columns:
            for col in self.columns:
                if column.startswith(col):
                    filter.append(column)
        pheno_df = pheno_df[filter]
        pheno_df.dropna(axis=0, inplace=True, how="all")

        sample_overlap = set(ge_df.index).intersection(set(pheno_df.index))
        print("\tSample overlap: {}".format(len(sample_overlap)))
        pheno_df = pheno_df.loc[sample_overlap, :]
        ge_df = ge_df.loc[sample_overlap, :]

        print(pheno_df)
        print(ge_df)

        print("Correlate data.")
        corr_coefficients = []
        pvalues = []
        for phenotype, pheno_data in pheno_df.T.iterrows():
            print("\tProcessing '{}'".format(phenotype))
            ct_coefficients = []
            ct_pvalues = []
            for i, (gene, expr_data) in enumerate(ge_df.T.iterrows()):
                # data = pd.concat([pheno_data, expr_data], axis=1)
                # data.dropna(axis=0, inplace=True)
                coef, p = stats.spearmanr(pheno_data, expr_data, nan_policy='omit')
                ct_coefficients.append(coef)
                ct_pvalues.append(p)
                # del data
            corr_coefficients.append(ct_coefficients)
            pvalues.append(ct_pvalues)

        corr_coef_df = pd.DataFrame(corr_coefficients,
                                    index=pheno_df.columns,
                                    columns=ge_df.columns)
        print(corr_coef_df)
        pvalue_df = pd.DataFrame(pvalues,
                                 index=pheno_df.columns,
                                 columns=ge_df.columns)
        print(pvalue_df)

        # Save.
        corr_coef_df.to_csv(os.path.join(self.outdir, '{}coefficients.txt.gz'.format(self.outfile_prefix)),
                            sep="\t", index=True, header=True, compression="gzip")
        pvalue_df.to_csv(os.path.join(self.outdir, '{}pvalues.txt.gz'.format(self.outfile_prefix)),
                         sep="\t", index=True, header=True, compression="gzip")

        #corr_coef_df = pd.read_csv(os.path.join(self.outdir, '{}coefficients.txt.gz'.format(self.outfile_prefix)), sep="\t", index_col=0, header=0)
        for index, values in corr_coef_df.iterrows():
            self.plot_distribution(data=values, name=index)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  low_memory=True):
        if path.endswith(".pkl"):
            df = pd.read_pickle(path)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def plot_distribution(self, data, name):
        fig, ax = plt.subplots()
        sns.set_style("darkgrid", {"axes.facecolor": ".9"})
        g = sns.distplot(data)
        g.set_title('Spearman Corr.Coef. Distribution')
        g.set_ylabel('Frequency')
        g.set_xlabel('Correlation Coefficient')
        fig.savefig(os.path.join(self.outdir, "{}_corr_coef_distribution.png".format(name)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Phenotype path: {}".format(self.pheno_path))
        print("  > Gene expression path: {}".format(self.ge_path))
        print("  > Expression ID: {}".format(self.expr_id))
        print("  > Columns: {}".format(self.columns))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
