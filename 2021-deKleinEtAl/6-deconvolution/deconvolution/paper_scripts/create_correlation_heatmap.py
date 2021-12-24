#!/usr/bin/env python3

"""
File:         create_correlation_heatmap.py
Created:      2021/12/08
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
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Create Correlation Matrix"
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
./create_correlation_heatmap.py -h

./create_correlation_heatmap.py -rd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table.txt.gz -rn cell fractions -o CellFrations_PsychENCODE_NoDev

./create_correlation_heatmap.py -rd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table.txt.gz -rn new_profile -cd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/partial_deconvolution/OLD_PROFILE_CORTEX_EUR_TMM_LOG2/IHC_0CPM_LOG2_FILTERED_CC/deconvolution.txt.gz -cn old_profile -o CellFrations_CellMapOriginal_vs_PsychENCODENoDev

./create_correlation_heatmap.py -rd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table_NoInhibitory.txt.gz -rn cell fractions -o CellFrations_PsychENCODE_NoDevNoInhibitory
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.row_data_path = getattr(arguments, 'row_data')
        self.row_name = " ".join(getattr(arguments, 'row_name'))
        self.col_data_path = getattr(arguments, 'col_data')
        self.col_name = " ".join(getattr(arguments, 'col_name'))
        self.method = getattr(arguments, 'method')
        self.out_filename = getattr(arguments, 'outfile')

        # Set variables.
        self.outdir = os.path.join(Path().resolve(), 'create_correlation_heatmap')
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
        parser.add_argument("-rd",
                            "--row_data",
                            type=str,
                            required=True,
                            help="The path to the data matrix.")
        parser.add_argument("-rn",
                            "--row_name",
                            nargs="*",
                            type=str,
                            required=False,
                            default="",
                            help="The name of -r / --row_data.")
        parser.add_argument("-cd",
                            "--col_data",
                            type=str,
                            required=False,
                            help="The path to the data matrix.")
        parser.add_argument("-cn",
                            "--col_name",
                            nargs="*",
                            type=str,
                            required=False,
                            default="",
                            help="The name of -c / --col_data.")
        parser.add_argument("-a",
                            "--axis",
                            type=int,
                            default=0,
                            choices=[0, 1],
                            help="The axis that denotes the samples. "
                                 "Default: 0")
        parser.add_argument("-m",
                            "--method",
                            type=str,
                            choices=["Pearson", "Spearman"],
                            default="Spearman",
                            help="The correlation method. Default: Spearman.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the outfile.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading row data.")
        row_df = self.load_file(self.row_data_path, header=0, index_col=0)

        col_df = row_df
        triangle = True
        if self.col_data_path is not None:
            print("Loading column data.")
            col_df = self.load_file(self.col_data_path, header=0, index_col=0)
            triangle = False

        if row_df.shape[1] > row_df.shape[0]:
            row_df = row_df.T

        if col_df.shape[1] > col_df.shape[0]:
            col_df = col_df.T

        print("Removing columns without variance.")
        row_df = row_df.loc[:, row_df.std(axis=0) != 0]
        col_df = col_df.loc[:, col_df.std(axis=0) != 0]

        print("Getting overlap.")
        overlap = set(row_df.index).intersection(set(col_df.index))
        print("\tN = {}".format(len(overlap)))
        if len(overlap) == 0:
            print("No data overlapping.")
            exit()
        row_df = row_df.loc[overlap, :]
        col_df = col_df.loc[overlap, :]

        print("\tCorrelating.")
        corr_df, pvalue_df = self.correlate(index_df=row_df, columns_df=col_df, triangle=triangle)
        print(corr_df)
        print(pvalue_df)

        print("\tMasking non-significant")
        signif_df = self.mask_non_significant(df=corr_df, pvalue_df=pvalue_df)

        print("\tPlotting.")
        self.plot_heatmap(df=signif_df, annot_df=corr_df.round(2), xlabel=self.col_name, ylabel=self.row_name, appendix="correlations")

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def correlate(self, index_df, columns_df, triangle=False):
        index_df_n_nan_values = index_df.shape[0] - index_df.isna().sum(axis=0)
        column_df_n_nan_values = columns_df.shape[0] - columns_df.isna().sum(axis=0)

        index_df_colnames = ["{} [N={:,}]".format(colname, index_df_n_nan_values[colname]) for colname in index_df.columns]
        column_df_colnames = ["{} [N={:,}]".format(colname, column_df_n_nan_values[colname]) for colname in columns_df.columns]

        corr_df = pd.DataFrame(np.nan, index=index_df_colnames, columns=column_df_colnames)
        pvalue_df = pd.DataFrame(np.nan, index=index_df_colnames, columns=column_df_colnames)

        for i, (index_column, index_colname) in enumerate(zip(index_df.columns, index_df_colnames)):
            for j, (column_column, column_colname) in enumerate(zip(columns_df.columns, column_df_colnames)):
                if triangle and i < j:
                    continue
                corr_data = pd.concat([index_df[index_column], columns_df[column_column]], axis=1)
                corr_data.dropna(inplace=True)

                coef = np.nan
                pvalue = np.nan
                if np.min(corr_data.std(axis=0)) > 0:
                    if self.method == "Pearson":
                        coef, pvalue = stats.pearsonr(corr_data.iloc[:, 1], corr_data.iloc[:, 0])
                    elif self.method == "Spearman":
                        coef, pvalue = stats.spearmanr(corr_data.iloc[:, 1], corr_data.iloc[:, 0])

                corr_df.loc[index_colname, column_colname] = coef
                pvalue_df.loc[index_colname, column_colname] = pvalue

        return corr_df, pvalue_df

    @staticmethod
    def mask_non_significant(df, pvalue_df, a=0.05):
        signif_df = df.copy()
        for i in range(signif_df.shape[0]):
            for j in range(signif_df.shape[1]):
                if np.isnan(pvalue_df.iloc[i, j]) or pvalue_df.iloc[i, j] >= a:
                    signif_df.iloc[i, j] = 0

        return signif_df

    def plot_heatmap(self, df, annot_df, xlabel="", ylabel="", appendix="",
                     vmin=-1, vmax=1):
        cmap = sns.diverging_palette(246, 24, as_cmap=True)

        fig, axes = plt.subplots(nrows=2,
                                 ncols=2,
                                 figsize=(1 * df.shape[1] + 10, 1 * df.shape[0] + 10),
                                 gridspec_kw={"width_ratios": [0.2, 0.8],
                                              "height_ratios": [0.8, 0.2]})
        sns.set(color_codes=True)

        annot_df.fillna("", inplace=True)

        row_index = 0
        col_index = 0
        for _ in range(4):
            ax = axes[row_index, col_index]
            if row_index == 0 and col_index == 1:
                sns.heatmap(df, cmap=cmap, vmin=vmin, vmax=vmax, center=0,
                            square=True, annot=annot_df, fmt='',
                            cbar=False, annot_kws={"size": 14, "color": "#000000"},
                            ax=ax)

                plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
                plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=90))

                ax.set_xlabel(xlabel, fontsize=14)
                ax.xaxis.set_label_position('top')

                ax.set_ylabel(ylabel, fontsize=14)
                ax.yaxis.set_label_position('right')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > 1:
                col_index = 0
                row_index += 1

        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_corr_heatmap_{}_{}.png".format(self.out_filename, self.method, appendix)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Row data path: {}".format(self.row_data_path))
        print("  > Row name: {}".format(self.row_name))
        print("  > Col data path: {}".format(self.col_data_path))
        print("  > Col name: {}".format(self.col_name))
        print("  > Correlation method: {}".format(self.method))
        print("  > Output filename: {}".format(self.out_filename))
        print("  > Outpath {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
