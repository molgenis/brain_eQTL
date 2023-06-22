#!/usr/bin/env python3

"""
File:         compare_fim_vs_decon_eqtl.py
Created:      2021/12/21
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
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
from statsmodels.stats import multitest
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.


# Metadata
__program__ = "Compare Fast Interaction Mapper vs Decon-eQTL"
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
./compare_fim_vs_decon_eqtl.py -h

./compare_fim_vs_decon_eqtl.py -fim /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/fast_interaction_mapper/2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron-NoDatasetCorrection-NoFN -de /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron -o 2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.fim_indir = getattr(arguments, 'fim_indir')
        self.de_indir = getattr(arguments, 'de_indir')
        self.out_filename = getattr(arguments, 'outfile')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.point_colormap = {
            "no signif": "#808080",
            "x signif": "#0072B2",
            "y signif": "#D55E00",
            "both signif": "#009E73"
        }

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
        parser.add_argument("-fim",
                            "--fim_indir",
                            type=str,
                            required=True,
                            help="The path to the Fast Interaction Mapper"
                                 " interaction directory.")
        parser.add_argument("-de",
                            "--de_indir",
                            type=str,
                            required=True,
                            help="The path to the Decon-eQTL interaction "
                                 "directory.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            default="output",
                            help="The name of the outfile. Default: output.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data")
        fim_df = self.load_fim_results(indir=self.fim_indir)
        de_df = self.load_de_results(indir=self.de_indir)
        print(fim_df)
        print(de_df)

        print("### Step2 ###")
        print("Overlap")
        row_overlap = set(fim_df.index).intersection(set(de_df.index))
        col_overlap = set(fim_df.columns).intersection(set(de_df.columns))
        fim_df = fim_df.loc[row_overlap, col_overlap]
        de_df = de_df.loc[row_overlap, col_overlap]
        print(fim_df)
        print(de_df)

        cols = list(col_overlap)
        cols.sort()

        print("### Step3 ###")
        print("Plot")
        self.plot_multiple_regplot(df1=fim_df,
                                   df2=de_df,
                                   xlabel="Standard ieQTL",
                                   ylabel="Decon-eQTL",
                                   cols=cols,
                                   filename=self.out_filename)

    @staticmethod
    def load_fim_results(indir):
        fpaths = glob.glob(os.path.join(indir, "*.txt.gz"))
        fdr_data = []
        for fpath in fpaths:
            if os.path.basename(fpath) in ["call_rate.txt.gz", "genotype_stats.txt.gz"]:
                continue
            df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
            df.index = df["gene"] + "_" + df["snp"]
            fdr_df = df[["ieQTL FDR"]]
            fdr_df.columns = [os.path.basename(fpath).split(".")[0]]
            fdr_data.append(fdr_df)
        fdr_df = pd.concat(fdr_data, axis=1)
        return fdr_df

    @staticmethod
    def load_de_results(indir):
        df = pd.read_csv(os.path.join(indir, "deconvolutionResults.txt.gz"), sep="\t", header=0, index_col=0)

        fdr_data = []
        indices = []
        for col in df.columns:
            if col.endswith("_pvalue"):
                fdr_data.append(multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
                indices.append(col.replace("_pvalue", ""))
        fdr_df = pd.DataFrame(fdr_data, index=indices, columns=df.index)

        return fdr_df.T

    def plot_multiple_regplot(self, df1, df2, cols, xlabel="", ylabel="",
                              filename="plot"):
        # Calculate number of rows / columns.
        ngroups = len(cols)
        ncols = int(np.ceil(np.sqrt(ngroups)))
        nrows = int(np.ceil(ngroups / ncols))

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharey="none",
                                 sharex="none",
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
                df = df1[[cols[i]]].merge(df2[[cols[i]]], left_index=True, right_index=True)
                df.columns = ["x", "y"]

                # Add color
                df["hue"] = self.point_colormap["no signif"]
                df.loc[(df["x"] <= 0.05) & (df["y"] > 0.05), "hue"] = self.point_colormap["x signif"]
                df.loc[(df["x"] > 0.05) & (df["y"] <= 0.05), "hue"] = self.point_colormap["y signif"]
                df.loc[(df["x"] <= 0.05) & (df["y"] <= 0.05), "hue"] = self.point_colormap["both signif"]

                counts = df["hue"].value_counts()
                for value in self.point_colormap.values():
                    if value not in counts.index:
                        counts[value] = 0

                df["x"] = np.log10(df["x"]) * -1
                df["y"] = np.log10(df["y"]) * -1

                coef, _ = stats.spearmanr(df["y"], df["x"])

                # Plot.
                g = sns.regplot(x="x",
                                y="y",
                                data=df,
                                scatter_kws={'facecolors': df["hue"],
                                             'linewidth': 0,
                                             'alpha': 0.5},
                                line_kws={"color": "#000000"},
                                ax=ax
                                )
                # Add the text.
                ax.annotate(
                    'r = {:.2f}'.format(coef),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color="#404040",
                    fontsize=18,
                    fontweight='bold')
                ax.annotate(
                    'total N = {:,}'.format(df.shape[0]),
                    xy=(0.03, 0.9),
                    xycoords=ax.transAxes,
                    color="#404040",
                    fontsize=18,
                    fontweight='bold')
                ax.annotate(
                    'N = {:,}'.format(counts[self.point_colormap["both signif"]]),
                    xy=(0.03, 0.86),
                    xycoords=ax.transAxes,
                    color=self.point_colormap["both signif"],
                    fontsize=18,
                    fontweight='bold')
                ax.annotate(
                    'N = {:,}'.format(counts[self.point_colormap["x signif"]]),
                    xy=(0.03, 0.82),
                    xycoords=ax.transAxes,
                    color=self.point_colormap["x signif"],
                    fontsize=18,
                    fontweight='bold')
                ax.annotate(
                    'N = {:,}'.format(counts[self.point_colormap["y signif"]]),
                    xy=(0.03, 0.78),
                    xycoords=ax.transAxes,
                    color=self.point_colormap["y signif"],
                    fontsize=18,
                    fontweight='bold')
                ax.annotate(
                    'N = {:,}'.format(counts[self.point_colormap["no signif"]]),
                    xy=(0.03, 0.74),
                    xycoords=ax.transAxes,
                    color=self.point_colormap["no signif"],
                    fontsize=18,
                    fontweight='bold')

                xlabel_str = ""
                if row_index == (nrows - 1):
                    xlabel_str = "-log10 p-value " + xlabel
                ax.set_xlabel(xlabel_str,
                              fontsize=14,
                              fontweight='bold')

                ylabel_str = ""
                if col_index == 0:
                    ylabel_str = "-log10 p-value " + ylabel
                ax.set_ylabel(ylabel_str,
                              fontsize=14,
                              fontweight='bold')

                ax.set_title(label,
                             color="#000000",
                             fontsize=18,
                             fontweight='bold')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        outpath = os.path.join(self.outdir, "{}_FastInteractionMapper_vs_DeconeQTL_regplot.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Fast Interaction Mapper input directory: {}".format(self.fim_indir))
        print("  > Decon-eQTL input directory: {}".format(self.de_indir))
        print("  > Output directory {}".format(self.outdir))
        print("  > Output filename: {}".format(self.out_filename))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
