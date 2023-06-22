#!/usr/bin/env python3

"""
File:         compare_eqtl_results.py
Created:      2022/04/01
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
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Local application imports.

# Metadata
__program__ = "Compare eQTL Results"
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
./compare_eqtl_results.py \
    -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -da AlleleAssessed \
    -dv OverallZScore \
    -ds FDR \
    -dn meta-analysis \
    -rd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/eqtl_mapper/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/eQTLResults.txt.gz \
    -ra AltAllele \
    -rv tvalue-genotype \
    -rs FDR \
    -rn bulk \
    -o 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected
    
./compare_eqtl_results.py \
    -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -da AlleleAssessed \
    -dv OverallZScore \
    -ds FDR \
    -dn meta-analysis \
    -rd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/eqtl_mapper/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-DeconExpressionMatrix/eQTLResults.txt.gz \
    -ra AltAllele \
    -rv tvalue-genotype \
    -rs FDR \
    -rn bulk \
    -o 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-DeconExpressionMatrix
    
### Harm-Jan Files ###

./compare_eqtl_results.py \
    -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/eqtl_mapper/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-DeconExpressionMatrix/eQTLResults.txt.gz \
    -da AltAllele \
    -dv tvalue-genotype \
    -ds FDR \
    -dn my_files \
    -rd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/eqtl_mapper/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-DeconExpressionMatrix-HJFiles/eQTLResults.txt.gz \
    -ra AltAllele \
    -rv tvalue-genotype \
    -rs FDR \
    -rn hj_files \
    -o 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-DeconExpressionMatrix-Myfiles_vs_HJFiles
    
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.disc_data_path = getattr(arguments, 'discovery_data')
        self.disc_aa = getattr(arguments, 'discovery_allele_assessed')
        self.disc_value = getattr(arguments, 'discovery_value')
        self.disc_signif = getattr(arguments, 'discovery_signif')
        self.disc_name = getattr(arguments, 'discovery_name')
        self.repl_data_path = getattr(arguments, 'replication_data')
        self.repl_aa = getattr(arguments, 'replication_allele_assessed')
        self.repl_value = getattr(arguments, 'replication_value')
        self.repl_signif = getattr(arguments, 'replication_signif')
        self.repl_name = getattr(arguments, 'replication_name')
        self.outname = getattr(arguments, 'outname')

        # Set variables.
        self.outdir = os.path.join(Path().resolve(), 'compare_eqtl_results')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
            "no signif": "#808080",
            "x signif": "#0072B2",
            "y signif": "#D55E00",
            "both signif": "#009E73"
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
        parser.add_argument("-dd",
                            "--discovery_data",
                            type=str,
                            required=True,
                            help="The path to the discovery data file")
        parser.add_argument("-da",
                            "--discovery_allele_assessed",
                            type=str,
                            required=False,
                            default="x",
                            help="The discovery allele assessed column name")
        parser.add_argument("-dv",
                            "--discovery_value",
                            type=str,
                            required=False,
                            default="x",
                            help="The discovery data value column name")
        parser.add_argument("-ds",
                            "--discovery_signif",
                            type=str,
                            required=False,
                            default="x",
                            help="The discovery data significance column name")
        parser.add_argument("-dn",
                            "--discovery_name",
                            type=str,
                            required=True,
                            help="The name for the discovery results")
        parser.add_argument("-rd",
                            "--replication_data",
                            type=str,
                            required=True,
                            help="The path to the replication data file")
        parser.add_argument("-ra",
                            "--replication_allele_assessed",
                            type=str,
                            required=False,
                            default="x",
                            help="The replication allele assessed column name")
        parser.add_argument("-rv",
                            "--replication_value",
                            type=str,
                            required=False,
                            default="x",
                            help="The replication data value column name")
        parser.add_argument("-rs",
                            "--replication_signif",
                            type=str,
                            required=False,
                            default="x",
                            help="The replication data significance column name")
        parser.add_argument("-rn",
                            "--replication_name",
                            type=str,
                            required=True,
                            help="The name for the replication results")
        parser.add_argument("-o",
                            "--outname",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output file.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        disc_df = self.load_file(self.disc_data_path, header=0, index_col=None)
        repl_df = self.load_file(self.repl_data_path, header=0, index_col=None)
        print(disc_df)
        print(repl_df)

        print("Preprocess data`")
        disc_df.index = disc_df["ProbeName"] + "_" + disc_df["SNPName"]
        repl_df.index = repl_df["ProbeName"] + "_" + repl_df["SNPName"]
        df = disc_df[[self.disc_value, self.disc_aa, self.disc_signif]].merge(repl_df[[self.repl_value, self.repl_aa, self.repl_signif]], left_index=True, right_index=True)
        df.columns = ["x", "x_aa", "x_signif", "y", "y_aa", "y_signif"]
        print(df)
        del disc_df, repl_df

        print("Fllipping effects")
        df["y"] = df["y"] * (df["x_aa"] == df["y_aa"]).map({True: 1, False: -1})
        print(df)

        print("Adding hue")
        df["hue"] = self.palette["no signif"]
        df.loc[(df["x_signif"] <= 0.05) & (df["y_signif"] > 0.05), "hue"] = self.palette["x signif"]
        df.loc[(df["x_signif"] > 0.05) & (df["y_signif"] <= 0.05), "hue"] = self.palette["y signif"]
        df.loc[(df["x_signif"] <= 0.05) & (df["y_signif"] <= 0.05), "hue"] = self.palette["both signif"]

        print("Plotting data")
        self.plot(df=df,
                  hue="hue",
                  palette=self.palette,
                  xlabel=self.disc_name,
                  ylabel=self.repl_name,
                  title="eQTL comparison",
                  filename="{}_eqtl_comparison".format(self.outname)
                  )

    @staticmethod
    def load_file(path, sep="\t", header=None, index_col=None, nrows=None):
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
              "with shape: {}".format(os.path.basename(outpath), df.shape))

    def plot(self, df, x="x", y="y", hue=None, palette=None,
                       xlabel=None, ylabel=None, title="", filename="plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        pearson_coef, _ = stats.pearsonr(df[y], df[x])

        lower_quadrant = df.loc[(df[x] < 0) & (df[y] < 0), :]
        upper_quadrant = df.loc[(df[x] > 0) & (df[y] > 0), :]
        concordance = (100 / df.shape[0]) * (lower_quadrant.shape[0] + upper_quadrant.shape[0])

        # Add the text.
        ax.annotate(
            'r = {:.2f}'.format(pearson_coef),
            xy=(0.03, 0.9),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=18,
            fontweight='bold')
        ax.annotate(
            'concordance = {:.0f}%'.format(concordance),
            xy=(0.03, 0.85),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=18,
            fontweight='bold')
        ax.annotate(
            'total N = {:,}'.format(df.shape[0]),
            xy=(0.03, 0.8),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=18,
            fontweight='bold')
        if hue is not None:
            counts = df["hue"].value_counts()
            for value in palette.values():
                if value not in counts.index:
                    counts[value] = 0

            ax.annotate(
                'N = {:,}'.format(counts[palette["both signif"]]),
                xy=(0.03, 0.75),
                xycoords=ax.transAxes,
                color=palette["both signif"],
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {:,}'.format(counts[palette["x signif"]]),
                xy=(0.03, 0.7),
                xycoords=ax.transAxes,
                color=palette["x signif"],
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {:,}'.format(counts[palette["y signif"]]),
                xy=(0.03, 0.65),
                xycoords=ax.transAxes,
                color=palette["y signif"],
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {:,}'.format(counts[palette["no signif"]]),
                xy=(0.03, 0.60),
                xycoords=ax.transAxes,
                color=palette["no signif"],
                fontsize=18,
                fontweight='bold')

        facecolors = "#000000"
        linecolor = "#b22222"
        if hue is not None:
            facecolors = df["hue"]
            linecolor = "#000000"

        sns.regplot(x=x, y=y, data=df,
                    scatter_kws={'facecolors': facecolors,
                                 'linewidth': 0,
                                 'alpha': 0.75},
                    line_kws={"color": linecolor,
                              'linewidth': 5},
                    ax=ax)

        ax.set_xlabel(xlabel,
                      fontsize=20,
                      fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=20,
                      fontweight='bold')
        ax.set_title(title,
                     fontsize=25,
                     fontweight='bold')

        # Change margins.
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        xmargin = (xlim[1] - xlim[0]) * 0.05
        ymargin = (ylim[1] - ylim[0]) * 0.05

        new_xlim = (xlim[0] - xmargin, xlim[1] + xmargin)
        new_ylim = (ylim[0] - ymargin, ylim[1] + ymargin)

        ax.set_xlim(new_xlim[0], new_xlim[1])
        ax.set_ylim(new_ylim[0], new_ylim[1])

        # Add diagonal.
        ax.axhline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)
        ax.axvline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)
        ax.axline((0, 0), slope=1, ls='--', color="#000000", alpha=0.15, zorder=-1)

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def plot_boxplot(self, df_m, x="x", y="value", hue=None, palette=None,
                     coefs_per_dataset=None, title="", xlabel="", ylabel="", filename=""):
        # Get the order.
        x_order = df_m[x].unique().tolist()
        x_order.sort()

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
                                       gridspec_kw={"width_ratios": [0.8, 0.2]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        sns.despine(fig=fig, ax=ax1)

        sns.violinplot(x=x,
                       y=y,
                       hue=hue,
                       data=df_m,
                       order=x_order,
                       palette=palette,
                       cut=0,
                       dodge=False,
                       ax=ax1)

        plt.setp(ax1.collections, alpha=.75)

        sns.boxplot(x=x,
                    y=y,
                    hue=hue,
                    data=df_m,
                    order=x_order,
                    whis=np.inf,
                    color="white",
                    dodge=False,
                    ax=ax1)

        ax1.set_title(title,
                      fontsize=20,
                      fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=20,
                       fontweight='bold')
        ax1.set_xlabel(xlabel,
                       fontsize=20,
                       fontweight='bold')

        if ax1.get_legend() is not None:
            ax1.get_legend().remove()

        if hue is not None:
            handles = []
            for x_group in x_order:
                if x_group in palette:
                    coef = np.nan
                    if x_group in coefs_per_dataset.index:
                        coef = coefs_per_dataset[x_group]
                    handles.append(mpatches.Patch(color=palette[x_group],
                                                  label="{} [r={:.2f}]".format(x_group, coef)))
            ax2.legend(handles=handles, loc="center")

        fig.savefig(os.path.join(self.outdir, "{}_correlations_PerDataset.png".format(filename)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Discovery:")
        print("    > Path: {}".format(self.disc_data_path))
        print("    > Allele assessed: {}".format(self.disc_aa))
        print("    > Value: {}".format(self.disc_value))
        print("    > Significance: {}".format(self.disc_signif))
        print("    > Name: {}".format(self.disc_name))
        print("  > Replication:")
        print("    > Path: {}".format(self.repl_data_path))
        print("    > Repl assessed: {}".format(self.repl_aa))
        print("    > Value: {}".format(self.repl_value))
        print("    > Significance: {}".format(self.repl_signif))
        print("    > Name: {}".format(self.repl_name))
        print("  > Output name: {}".format(self.outname))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
