#!/usr/bin/env python3

"""
File:         visualise_decon_eqtl.py
Created:      2021/12/20
Last Changed: 2022/02/10
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
from scipy import stats
from statsmodels.stats import multitest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Visualise Decon-eQTL"
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
./visualise_decon_eqtl.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron -o 2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron

./visualise_decon_eqtl.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDevNoInhibitoryCT -o 2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDevNoInhibitoryCT
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_dir = getattr(arguments, 'data')
        self.output_filename = getattr(arguments, 'output')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "OtherNeuron": "#2690ce",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Astrocyte": "#D55E00"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the input data.")
        parser.add_argument("-o",
                            "--output",
                            type=str,
                            default="PlotPerColumn_ColorByCohort",
                            help="The name of the output file.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")
        decon_df = self.load_file(os.path.join(self.data_dir, "deconvolutionResults.txt.gz"), header=0, index_col=0)
        columns = [x for x in decon_df.columns if "pvalue" in x]
        decon_pval_df = decon_df[columns]
        decon_pval_df.columns = [col.replace("_pvalue", "") for col in decon_pval_df.columns]

        print("Plotting")
        self.plot_distributions(df=decon_pval_df, name="pvalue_")

        print("Calculating FDR.")
        decon_fdr_df = self.bh_correct(decon_pval_df)
        print(decon_fdr_df)
        print(decon_fdr_df.loc[decon_fdr_df["EndothelialCell"] <= 0.05, :])
        exit()

        print("Plotting")
        self.plot_distributions(df=decon_fdr_df, name="FDR_")

        print("Calculating max p-value with FDR <0.05")
        for col in decon_pval_df.columns:
            compare_df = decon_pval_df[[col]].merge(decon_fdr_df[[col]], left_index=True, right_index=True)
            compare_df.columns = ["p-value", "FDR"]
            signif_compare_df = compare_df.loc[compare_df["FDR"] <= 0.05, :]
            print("\t{}: {:.2e}".format(col, signif_compare_df["p-value"].max()))

        # print("Calculating the avg missingess per ieQTL")
        geno_stats_df = self.load_file(os.path.join(self.data_dir, "geno_stats.txt.gz"), header=0, index_col=0)
        geno_stats_copy_df = geno_stats_df.copy()
        geno_stats_copy_df = geno_stats_copy_df.loc[geno_stats_copy_df["mask"] == 1, :]
        geno_stats_copy_df.index = decon_fdr_df.index
        print(geno_stats_copy_df)
        ct_geno_stats_min = []
        ct_geno_stats_max = []
        ct_geno_stats_mean = []
        for col in decon_fdr_df.columns:
            mask = (decon_fdr_df[col] <= 0.05).to_numpy(dtype=bool)
            ct_geno_stats_min.append(geno_stats_copy_df.loc[mask, :].min())
            ct_geno_stats_max.append(geno_stats_copy_df.loc[mask, :].max())
            ct_geno_stats_mean.append(geno_stats_copy_df.loc[mask, :].mean())
        for data in [ct_geno_stats_min, ct_geno_stats_max, ct_geno_stats_mean]:
            ct_geno_stats_df = pd.concat(data, axis=1)
            ct_geno_stats_df.columns = decon_fdr_df.columns
            print(ct_geno_stats_df)

        print("Plotting")
        for col in decon_fdr_df.columns:
            plot_df = decon_fdr_df[[col]].merge(geno_stats_copy_df[["NaN"]], left_index=True, right_index=True)
            plot_df.columns = ["x", "y"]
            self.single_regplot(df=plot_df,
                                xlabel="FDR",
                                ylabel="#NaN",
                                title=col,
                                filename="{}_{}_FDR_vs_NaN".format(self.output_filename, col))

        print("Plotting eQTL FDR vs ieQTL FDR")
        eqtl_df = self.load_file("/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/combine_eqtlprobes/eQTLprobes_combined.txt.gz", header=0, index_col=None)
        eqtl_df = eqtl_df.loc[(geno_stats_df["mask"] == 1).to_numpy(dtype=bool), :]
        eqtl_df.index = eqtl_df["ProbeName"] + "_" + eqtl_df["SNPName"]
        print(eqtl_df)
        print(decon_fdr_df)
        for col in decon_fdr_df.columns:
            plot_df = eqtl_df[["qval"]].merge(decon_fdr_df[[col]], left_index=True, right_index=True)
            plot_df.columns = ["x", "y"]
            plot_df = plot_df.loc[plot_df["y"] <= 0.05, :]
            plot_df = np.log10(plot_df) * -1
            self.single_regplot(df=plot_df,
                                xlabel="-log10 eQTL q-value",
                                ylabel="-log10 ieQTL FDR",
                                title=col,
                                filename="{}_{}_eQTL_qval_vs_ieQTL_FDR".format(self.output_filename, col))

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def bh_correct(pvalue_df):
        df = pvalue_df.copy()
        fdr_data = []
        for col in df.columns:
            fdr_data.append(multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
        fdr_df = pd.DataFrame(fdr_data, index=df.columns, columns=df.index)

        return fdr_df.T

    def plot_distributions(self, df, name=""):
        cols = list(df.columns)
        cols.sort()

        ngroups = len(cols)
        ncols = int(np.ceil(np.sqrt(ngroups)))
        nrows = int(np.ceil(ngroups / ncols))

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharey="all",
                                 sharex="all",
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

            if i < ngroups:
                sns.despine(fig=fig, ax=ax)

                label = cols[i]
                color = self.palette[label.split("_")[0]]
                x = df.loc[:, label].to_numpy()

                sns.despine(fig=fig, ax=ax)

                sns.kdeplot(x, shade=True, color=color, ax=ax, cut=0, zorder=-1)
                ax.axvline(x.mean(), ls='--', color="#808080", zorder=-1)

                ax.annotate(
                    'N = {:,}'.format(np.size(x)),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=15,
                    fontweight='bold')
                ax.annotate(
                    'Mean = {:.2f}'.format(np.mean(x)),
                    xy=(0.03, 0.9),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=15,
                    fontweight='bold')
                ax.annotate(
                    'sd = {:.2f}'.format(np.std(x)),
                    xy=(0.03, 0.86),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=15,
                    fontweight='bold')

                ax.set_xlabel("",
                              fontsize=14,
                              fontweight='bold')
                ax.set_ylabel("density",
                              fontsize=14,
                              fontweight='bold')
                ax.set_title(label,
                             color=color,
                             fontsize=18,
                             fontweight='bold')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        outpath = os.path.join(self.outdir, "{}_{}distribution.png".format(self.output_filename, name))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def single_regplot(self, df, x="x", y="y", hue=None, palette=None,
                       xlabel=None, ylabel=None, title="", filename="plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
                                       gridspec_kw={"width_ratios": [0.8, 0.2]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        # Set annotation.
        pearson_coef, _ = stats.pearsonr(df[y], df[x])
        ax1.annotate(
            'total N = {:,}'.format(df.shape[0]),
            xy=(0.03, 0.94),
            xycoords=ax1.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold')
        ax1.annotate(
            'total r = {:.2f}'.format(pearson_coef),
            xy=(0.03, 0.90),
            xycoords=ax1.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold')

        group_column = hue
        if hue is None:
            df["hue"] = "#000000"
            group_column = "hue"

        group_corr_coef = {}
        group_sizes = {}
        for i, hue_group in enumerate(df[group_column].unique()):
            subset = df.loc[df[group_column] == hue_group, :]
            if subset.shape[0] < 2:
                continue

            facecolors = "#000000"
            color = "#b22222"
            if palette is not None:
                facecolors = palette[hue_group]
                color = facecolors

            sns.regplot(x=x, y=y, data=subset, ci=None,
                        scatter_kws={'facecolors': facecolors,
                                     'linewidth': 0},
                        line_kws={"color": color},
                        ax=ax1)

            if hue is not None:
                subset_pearson_coef, _ = stats.pearsonr(subset[y], subset[x])
                group_corr_coef[hue_group] = subset_pearson_coef
                group_sizes[hue_group] = subset.shape[0]

        if hue is not None:
            handles = []
            for hue_group in df[group_column].unique():
                if hue_group in palette:
                    n = "0"
                    if hue_group in group_sizes:
                        n = "{:,}".format(group_sizes[hue_group])
                    r = "NA"
                    if hue_group in group_corr_coef:
                        n = "{:.2f}".format(group_corr_coef[hue_group])
                    handles.append(mpatches.Patch(color=palette[hue_group],
                                                  label="{} [n={}; r={}]".format(hue_group, n, r)))
            ax2.legend(handles=handles, loc="center", fontsize=8)

        ax1.set_xlabel(xlabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_title(title,
                      fontsize=18,
                      fontweight='bold')

        # Change margins.
        xlim = ax1.get_xlim()
        ylim = ax1.get_ylim()

        xmargin = (xlim[1] - xlim[0]) * 0.05
        ymargin = (ylim[1] - ylim[0]) * 0.05

        new_xlim = (xlim[0] - xmargin, xlim[1] + xmargin)
        new_ylim = (ylim[0] - ymargin, ylim[1] + ymargin)

        ax1.set_xlim(new_xlim[0], new_xlim[1])
        ax1.set_ylim(new_ylim[0], new_ylim[1])

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Data directory: {}".format(self.data_dir))
        print("  > Output filename: {}".format(self.output_filename))
        print("  > Output directory {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
