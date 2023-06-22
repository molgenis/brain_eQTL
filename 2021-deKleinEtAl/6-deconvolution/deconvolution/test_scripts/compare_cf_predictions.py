#!/usr/bin/env python3

"""
File:         comnpare_cf_predictions.py
Created:      2021/11/22
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
import math
import os

# Third party imports.
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Compare Cell Fraction Predictions"
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
./compare_cf_predictions.py -cf1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/partial_deconvolution/PSYCHENCODE_PROFILE_CORTEX_EUR_TMM_LOG2/IHC_0CPM_LOG2/deconvolution_raw.txt.gz -transpose1 -n1 MetaBrain -cf2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/Cell_fractions_Raw.txt.gz -n2 PsychENCODE_MyPredictions -o MetaBrain_vs_PsychENCODE_MyPredictions_CellFraction_comparison

./compare_cf_predictions.py -cf1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/partial_deconvolution/PSYCHENCODE_PROFILE_CORTEX_EUR_TPM_LOG2/NoENA_IHC_0CPM_LOG2/deconvolution_raw.txt.gz -transpose1 -n1 MetaBrain_TPM_NoENA -cf2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/Cell_fractions_Raw.txt.gz -n2 PsychENCODE_MyPredictions -o MetaBrain_TPM_NoENA_vs_PsychENCODE_MyPredictions_CellFraction_comparison

./compare_cf_predictions.py -cf1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table.txt.gz -transpose1 -n1 TPMLog2 -cf2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-19-CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -transpose2 -n2 TMMLog2CovariatesRemovedOLS -o MetaBrain_CortexEUR_cis_TPM_Log2_CovariatesRemovedOrNot

./compare_cf_predictions.py -cf1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table.txt.gz -transpose1 -n1 TPMLog2 -cf2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-19-CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -transpose2 -n2 TMMLog2CovariatesRemovedOLS -o MetaBrain_CortexEUR_cis_TPM_Log2_CovariatesRemovedOrNot

./compare_cf_predictions.py -cf1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-19-CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -transpose1 -n1 NegativeToZero -cf2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-20-CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -transpose2 -n2 ShiftedPositive -o MetaBrain_CortexEUR_cis_NegativeToZero_vs_ShiftedPositive

./compare_cf_predictions.py -cf1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-20-CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -transpose1 -n1 WithIntercept -cf2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-20-CortexEUR-cis-NoIntercept/perform_deconvolution/deconvolution_table.txt.gz -transpose2 -n2 WithoutIntercept -o MetaBrain_CortexEUR_cis_ShiftedPositive_WithOrWithoutIntercept

./compare_cf_predictions.py -cf1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-RAMCorrected/perform_deconvolution/deconvolution_table.txt.gz -transpose1 -n1 NoDatasetCorrection -cf2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz -transpose2 -n2 DatasetCorrection -o 2022-01-21-CortexEUR-cis-NegativeToZero-RAMCorrected-DatasetCorrectedOrNot

"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cell_fractions1_path = getattr(arguments, 'cell_fractions1')
        self.transpose1 = getattr(arguments, 'transpose1')
        self.name1 = getattr(arguments, 'name1')
        self.cell_fractions2_path = getattr(arguments, 'cell_fractions2')
        self.transpose2 = getattr(arguments, 'transpose2')
        self.name2 = getattr(arguments, 'name2')
        self.outfile = getattr(arguments, 'outfile')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.colormap = {
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "OPC": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Pericytes": "#808080",
            "OtherNeuron": "#2690ce"
        }

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-cf1",
                            "--cell_fractions1",
                            type=str,
                            required=True,
                            help="The path to the first cell fraction matrix")
        parser.add_argument("-transpose1",
                            action='store_true',
                            help="Transpose 1.")
        parser.add_argument("-n1",
                            "--name1",
                            type=str,
                            required=False,
                            default="x",
                            help="The name for the first profile matrix")
        parser.add_argument("-cf2",
                            "--cell_fractions2",
                            type=str,
                            required=True,
                            help="The path to the second cell fraction matrix")
        parser.add_argument("-transpose2",
                            action='store_true',
                            help="Transpose 2.")
        parser.add_argument("-n2",
                            "--name2",
                            type=str,
                            required=False,
                            default="y",
                            help="The name for the second profile matrix")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=False,
                            default="deconvolution",
                            help="The name of the plot file.")

        return parser.parse_args()

    def start(self):
        print("Loading profile 1.")
        cf1_df = pd.read_csv(self.cell_fractions1_path, sep="\t", header=0, index_col=0)
        if self.transpose1:
            cf1_df = cf1_df.T
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.cell_fractions1_path),
                                      cf1_df.shape))
        print(cf1_df)

        print("Loading profile 2.")
        cf2_df = pd.read_csv(self.cell_fractions2_path, sep="\t", header=0, index_col=0)
        if self.transpose2:
            cf2_df = cf2_df.T
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.cell_fractions2_path),
                                      cf2_df.shape))
        print(cf2_df)

        print("Subset overlap in samples and cell types")
        col_overlap = set(cf1_df.columns).intersection(set(cf2_df.columns))
        row_overlap = set(cf1_df.index).intersection(set(cf2_df.index))
        print("\tN col-overlap: {}".format(len(col_overlap)))
        print("\tN row-overlap: {}".format(len(row_overlap)))
        cf1_df = cf1_df.loc[row_overlap, col_overlap].copy()
        cf2_df = cf2_df.loc[row_overlap, col_overlap].copy()

        print("Correlating.")
        corr_df, pvalue_df = self.correlate(index_df=cf1_df.T, columns_df=cf2_df.T)
        signif_df = self.mask_non_significant(df=corr_df, pvalue_df=pvalue_df)
        self.plot_heatmap(df=signif_df,
                          annot_df=corr_df.round(2),
                          xlabel=self.name2,
                          ylabel=self.name1)

        print("Merge")
        cf1_df = cf1_df.T
        cf1_df.reset_index(drop=False, inplace=True)
        cf1_dfm = cf1_df.melt(id_vars="index")
        cf1_dfm.columns = ["sample", "cell type", "x"]

        cf2_df = cf2_df.T
        cf2_df.reset_index(drop=False, inplace=True)
        cf2_dfm = cf2_df.melt(id_vars="index")
        cf2_dfm.columns = ["sample", "cell type", "y"]

        cf_dfm = cf1_dfm.merge(cf2_dfm, on=["sample", "cell type"])

        print("Plotting.")
        self.plot_regplot(
            df=cf_dfm,
            group_column="cell type",
            xlabel=self.name1,
            ylabel=self.name2,
            palette=self.colormap)

    @staticmethod
    def correlate(index_df, columns_df, triangle=False):
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
                    coef, pvalue = stats.pearsonr(corr_data.iloc[:, 1], corr_data.iloc[:, 0])

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

    def plot_heatmap(self, df, annot_df, xlabel="", ylabel="", vmin=-1, vmax=1):
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
        fig.savefig(os.path.join(self.outdir, "{}_corr_heatmap_Pearson.png".format(self.outfile)))
        plt.close()

    def plot_regplot(self, df, group_column, x="x", y="y", xlabel="",
                     ylabel="", palette=None, filename=""):

        if df.shape[0] <= 2:
            return

        group_counts = list(zip(*np.unique(df[group_column].to_numpy(), return_counts=True)))
        group_counts.sort(key=lambda x: -x[1])
        groups = [x[0] for x in group_counts]
        groups.sort()

        nplots = len(groups)
        ncols = math.ceil(np.sqrt(nplots))
        nrows = math.ceil(nplots / ncols)

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 figsize=(12 * ncols, 12 * nrows))
        sns.set(color_codes=True)

        row_index = 0
        col_index = 0
        for i in range(ncols * nrows):
            if nrows == 1 and ncols == 1:
                ax = axes
            elif nrows == 1 and ncols > 1:
                ax = axes[col_index]
            elif nrows > 1 and ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            if i < nplots:
                plot_df = df.loc[df[group_column] == groups[i], :]

                sns.despine(fig=fig, ax=ax)

                lower_quadrant = plot_df.loc[(plot_df[x] < 0) & (plot_df[y] < 0), :]
                upper_quadrant = plot_df.loc[(plot_df[x] > 0) & (plot_df[y] > 0), :]
                concordance = (100 / plot_df.shape[0]) * (lower_quadrant.shape[0] + upper_quadrant.shape[0])

                coef, _ = stats.spearmanr(plot_df[y], plot_df[x])

                accent_color = "#b22222"
                if palette is not None:
                    accent_color = palette[groups[i]]

                sns.regplot(x=x, y=y, data=plot_df, ci=None,
                            scatter_kws={'facecolors': "#000000",
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": accent_color,
                                      'linewidth': 5},
                            ax=ax)

                ax.annotate(
                    'N = {}'.format(plot_df.shape[0]),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color=accent_color,
                    alpha=1,
                    fontsize=18,
                    fontweight='bold')
                ax.annotate(
                    'r = {:.2f}'.format(coef),
                    xy=(0.03, 0.90),
                    xycoords=ax.transAxes,
                    color=accent_color,
                    alpha=1,
                    fontsize=18,
                    fontweight='bold')
                ax.annotate(
                    'concordance = {:.0f}%'.format(concordance),
                    xy=(0.03, 0.86),
                    xycoords=ax.transAxes,
                    color=accent_color,
                    alpha=1,
                    fontsize=18,
                    fontweight='bold')

                ax.axhline(0, ls='--', color="#000000", zorder=-1)
                ax.axvline(0, ls='--', color="#000000", zorder=-1)

                tmp_xlabel = ""
                if row_index == (nrows - 1):
                    tmp_xlabel = xlabel
                ax.set_xlabel(tmp_xlabel,
                              fontsize=20,
                              fontweight='bold')
                tmp_ylabel = ""
                if col_index == 0:
                    tmp_ylabel = ylabel
                ax.set_ylabel(tmp_ylabel,
                              fontsize=20,
                              fontweight='bold')

                ax.set_title(groups[i],
                             fontsize=25,
                             fontweight='bold')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        outpath = os.path.join(self.outdir, "{}_regplot.png".format(self.outfile))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved: {}".format(outpath))


if __name__ == '__main__':
    m = main()
    m.start()
