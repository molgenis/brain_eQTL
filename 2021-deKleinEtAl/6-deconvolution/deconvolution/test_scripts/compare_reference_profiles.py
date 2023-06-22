#!/usr/bin/env python3

"""
File:         compare_reference_profiles.py
Created:      2021/09/01
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
__program__ = "Compare Reference Profiles"
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
./compare_reference_profiles.py -p1 ../data/CellMap_brain_celltype_avgCPM.txt -n1 OLD -p2 ../data/CellMap_brain_CNS7_avgCPM.txt -n2 NEW -o OLD_vs_NEW

./compare_reference_profiles.py -p1 ../data/CellMap_brain_CNS7_avgCPM.txt -n1 METABRAIN -p2 ../data/PyschENCODE_average_signature_df.txt.gz -n2 PsychENCODE -o METABRAIN_VS_PsycheENCODE

./compare_reference_profiles.py -p1 ../data/CellMap_brain_CNS7_avgCPM.txt -n1 METABRAIN -p2 ../data/PsychENCODE_ref_profile_unique_avg.txt.gz -n2 PsychENCODE -o METABRAIN_VS_PsycheENCODE
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.profile1_path = getattr(arguments, 'profile1')
        self.profile1_name = getattr(arguments, 'name1')
        self.profile2_path = getattr(arguments, 'profile2')
        self.profile2_name = getattr(arguments, 'name2')
        self.outfile = getattr(arguments, 'outfile')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.accent_colormap = {
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Average": "#E8E8E8"
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
        parser.add_argument("-p1",
                            "--profile1",
                            type=str,
                            required=True,
                            help="The path to the first profile matrix")
        parser.add_argument("-n1",
                            "--name1",
                            type=str,
                            required=False,
                            default="x",
                            help="The name for the first profile matrix")
        parser.add_argument("-p2",
                            "--profile2",
                            type=str,
                            required=True,
                            help="The path to the second profile matrix")
        parser.add_argument("-n2",
                            "--name2",
                            type=str,
                            required=False,
                            default="y",
                            help="The name for the second profile matrix")
        parser.add_argument("-log10",
                            action='store_true',
                            help="Log10 transform the profile values."
                                 " Default: False.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=False,
                            default="deconvolution",
                            help="The name of the plot file.")

        return parser.parse_args()

    def start(self):
        print("Loading profile 1.")
        profile1_df = pd.read_csv(self.profile1_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.profile1_path),
                                      profile1_df.shape))
        # profile1_df.columns = [x.split("_")[1] for x in profile1_df.columns]
        profile1_df = np.log2(profile1_df + 1)
        print(profile1_df)

        print("Loading profile 2.")
        profile2_df = pd.read_csv(self.profile2_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.profile2_path),
                                      profile2_df.shape))
        profile2_df = np.log2(profile2_df + 1)

        # PsychENCODE stuff.
        gene_trans_df = pd.read_csv("../data/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.PsychENCODEGenesExtended.txt.gz", sep="\t", index_col=None)
        print(gene_trans_df)
        gene_trans_dict = dict(zip(gene_trans_df["Symbol"], [x.split(".")[0] for x in gene_trans_df["ArrayAddress"]]))
        profile1_df.index = [gene_trans_dict[hgnc] if hgnc in gene_trans_dict else hgnc for hgnc in profile1_df.index]
        profile2_df.index = [gene_trans_dict[hgnc] if hgnc in gene_trans_dict else hgnc for hgnc in profile2_df.index]
        profile2_df = profile2_df.loc[:, ['Adult-Ex1', 'Adult-Ex2', 'Adult-Ex3', 'Adult-Ex4', 'Adult-Ex5', 'Adult-Ex6', 'Adult-Ex7', 'Adult-Ex8', 'Adult-In1', 'Adult-In2', 'Adult-In3', 'Adult-In4', 'Adult-In5', 'Adult-In6', 'Adult-In7', 'Adult-In8', 'Adult-Astrocytes', 'Adult-Endothelial', 'Dev-Quiescent', 'Dev-Replicating', 'Adult-Microglia', 'Adult-OtherNeuron', 'Adult-OPC', 'Adult-Oligo']]
        print(profile1_df)
        print(profile2_df)

        for profile_df, name in ((profile1_df, self.profile1_name), (profile2_df, self.profile2_name)):
            corr_df = self.correlate(index_df=profile_df,
                                     columns_df=profile_df,
                                     triangle=True)
            self.plot_heatmap(corr_df=corr_df, filename="{}_{}".format(self.outfile, name))

        print("Correlate profiles")
        row_overlap = set(profile1_df.index).intersection(set(profile2_df.index))
        print("\tN-overlap: {}".format(len(row_overlap)))
        overlap1_df = profile1_df.loc[row_overlap, :].copy()
        overlap2_df = profile2_df.loc[row_overlap, :].copy()
        # overlap1_df = overlap1_df.sort_index(axis=1)
        # overlap2_df = overlap2_df.sort_index(axis=1)
        # overlap1_df.columns = [x.split("_")[1] for x in overlap1_df.columns]
        # overlap2_df.columns = [x.split("_")[1] for x in overlap2_df.columns]
        corr_df = self.correlate(index_df=overlap1_df, columns_df=overlap2_df)
        print(corr_df)

        print("\tPlotting.")
        self.plot_heatmap(corr_df=corr_df, xlabel=self.profile2_name, ylabel=self.profile1_name, filename="{}_{}".format(self.outfile, "overlap"))
        exit()

        print("Compare cell types")
        for profile_df, name in ((profile1_df, self.profile1_name), (profile2_df, self.profile2_name)):
            for i, column1 in enumerate(profile_df.columns):
                column1_name = column1.split("_")[1]
                for j, column2 in enumerate(profile_df.columns):
                    column2_name = column2.split("_")[1]

                    if i < j:
                        self.plot_regplot(df=profile_df,
                                          x=column1,
                                          y=column2,
                                          xlabel=column1_name,
                                          ylabel=column2_name,
                                          title="{} profile".format(name),
                                          filename="{}_{}_{}_vs_{}".format(self.outfile, name, column1_name, column2_name))

    @staticmethod
    def correlate(index_df, columns_df, triangle=False):
        out_df = pd.DataFrame(np.nan, index=index_df.columns, columns=columns_df.columns)

        for i, index_column in enumerate(index_df.columns):
            for j, column_column in enumerate(columns_df.columns):
                if triangle and i <= j:
                    continue
                corr_data = pd.concat([index_df[index_column], columns_df[column_column]], axis=1)
                filtered_corr_data = corr_data.dropna()
                coef, _ = stats.pearsonr(filtered_corr_data.iloc[:, 1], filtered_corr_data.iloc[:, 0])

                out_df.loc[index_column, column_column] = coef

        return out_df

    def plot_heatmap(self, corr_df, xlabel="", ylabel="", filename="plot"):
        cmap = sns.diverging_palette(246, 24, as_cmap=True)

        fig, axes = plt.subplots(nrows=2,
                                 ncols=2,
                                 figsize=(1 * corr_df.shape[1] + 5, 1 * corr_df.shape[0] + 5),
                                 gridspec_kw={"width_ratios": [0.2, 0.8],
                                              "height_ratios": [0.8, 0.2]})
        sns.set(color_codes=True)

        row_index = 0
        col_index = 0
        for _ in range(4):
            ax = axes[row_index, col_index]
            if row_index == 0 and col_index == 1:

                sns.heatmap(corr_df, cmap=cmap, vmin=-1, vmax=1, center=0,
                            square=True, annot=corr_df.round(2), fmt='',
                            cbar=False, annot_kws={"size": 16, "color": "#000000"},
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

        # plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_corr_heatmap.png".format(filename)))
        plt.close()

    def plot_regplot(self, df, x="x", y="y", xlabel="", ylabel="", title="",
                     filename="plot"):

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        # Plot.
        sns.regplot(x=x, y=y, data=df, ci=95,
                    scatter_kws={'facecolors': "#000000",
                                 'linewidth': 0,
                                 'alpha': 0.75},
                    line_kws={"color": "#b22222"},
                    ax=ax
                    )

        # Regression.
        coef, _ = stats.spearmanr(df[y], df[x])

        # Add the text.
        ax.annotate(
            'r = {:.2f}'.format(coef),
            xy=(0.03, 0.9),
            xycoords=ax.transAxes,
            color="#b22222",
            alpha=0.75,
            fontsize=14,
            fontweight='bold')

        ax.set_title(title,
                     fontsize=18,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        fig.savefig(os.path.join(self.outdir, "{}_regplot.png".format(filename)))
        print("\tSaved '{}'".format(filename))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
