#!/usr/bin/env python3

"""
File:         compare_bulk_expression.py
Created:      2021/11/25
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
__program__ = "Compare Bulk Expression"
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
./compare_bulk_expression.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/translated/DER-02_PEC_Gene_expression_matrix_TPM_MGFiltered_LOG2.txt.gz -n1 PsychENCODE -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/partial_deconvolution/PSYCHENCODE_PROFILE_CORTEX_EUR_TPM_LOG2/expression.txt.gz -n2 MetaBrain -o PsychENCODE_vs_MetaBrain

./compare_bulk_expression.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/translated/DER-02_PEC_Gene_expression_matrix_TPM.txt.gz -n1 PsychENCODE -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrain_Cortex/data/geneCounts.TPM.MergedExonLength.txt.gz -n2 MetaBrain -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex.txt.gz

./compare_bulk_expression.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/translated/DER-02_PEC_Gene_expression_matrix_TPM.txt.gz -n1 PsychENCODE -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrain_Cortex/data/geneCounts.TPM.TotalGeneLength.txt.gz -n2 MetaBrain -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex.txt.gz

./compare_bulk_expression.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/translated/DER-02_PEC_Gene_expression_matrix_TPM.txt.gz -n1 PsychENCODE -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrain_Unstranded/data/geneCounts.TPM.MergedExonLength.txt.gz -n2 MetaBrain -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_all.txt.gz -o Unstranded_MergedExonLength

./compare_bulk_expression.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/translated/DER-02_PEC_Gene_expression_matrix_TPM.txt.gz -n1 PsychENCODE -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrain_FirstStrand/data/geneCounts.TPM.MergedExonLength.txt.gz -n2 MetaBrain -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_all.txt.gz -o FirstStrand_MergedExonLength

./compare_bulk_expression.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/translated/DER-02_PEC_Gene_expression_matrix_TPM.txt.gz -n1 PsychENCODE -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrain_SecondStrand/data/geneCounts.TPM.MergedExonLength.txt.gz -n2 MetaBrain -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_all.txt.gz -o SecondStrand_MergedExonLength

./compare_bulk_expression.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/translated/DER-02_PEC_Gene_expression_matrix_TPM.txt.gz -n1 PsychENCODE -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrain_SecondStrand/data/geneCounts.TPM.TotalGeneLength.txt.gz -n2 MetaBrain -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_all.txt.gz -o SecondStrand_TotalGeneLength
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data1_path = getattr(arguments, 'data1')
        self.name1 = getattr(arguments, 'name1')
        self.data2_path = getattr(arguments, 'data2')
        self.name2 = getattr(arguments, 'name2')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        self.outdir = os.path.join(Path().resolve(), 'compare_bulk_expression', self.outfolder)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
            "AMPAD-MAYO-V2": "#9C9FA0",
            "CMC_HBCC_set2": "#0877B4",
            "GTEx": "#0FA67D",
            "AMPAD-ROSMAP-V2": "#6950A1",
            "BrainGVEX-V2": "#48B2E5",
            "TargetALS": "#D5C77A",
            "AMPAD-MSBB-V2": "#5CC5BF",
            "NABEC-H610": "#6D743A",
            "LIBD_1M": "#808080",
            "LIBD_5M": "#808080",
            "ENA": "#D46727",
            "LIBD_h650": "#808080",
            "GVEX": "#48B2E5",
            "NABEC-H550": "#6D743A",
            "CMC_HBCC_set3": "#0877B4",
            "UCLA_ASD": "#F36D2A",
            "CMC": "#EAE453",
            "CMC_HBCC_set1": "#0877B4",
            "Braineac": "#E49D26",
            "Bipseq_1M": "#000000",
            "Bipseq_h650": "#000000",
            "Brainseq": "#C778A6",
            "EUR": "#0fa67d",
            "AFR": "#0877b4",
            "EAS": "#d46727",
            "AMR": "#e49d26",
            "SAS": "#6950A1"
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
        parser.add_argument("-d1",
                            "--data1",
                            type=str,
                            required=True,
                            help="The path to the first cell fraction matrix")
        parser.add_argument("-n1",
                            "--name1",
                            type=str,
                            required=False,
                            default="x",
                            help="The name for the first profile matrix")
        parser.add_argument("-d2",
                            "--data2",
                            type=str,
                            required=True,
                            help="The path to the second cell fraction matrix")
        parser.add_argument("-n2",
                            "--name2",
                            type=str,
                            required=False,
                            default="y",
                            help="The name for the second profile matrix")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=False,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # print("Loading data")
        # std = self.load_file(self.std_path, header=0, index_col=None)
        # std.index = std["rnaseq_id"]
        # df1 = self.load_file(self.data1_path, header=0, index_col=0)
        # df2 = self.load_file(self.data2_path, header=0, index_col=0)
        # df2.index = [x.split(".")[0] for x in df2.index]
        #
        # print(std)
        # print(df1)
        # print(df2)
        #
        # print("Removing duplicates")
        # df1 = df1.groupby(df1.index).first()
        # df2 = df2.groupby(df2.index).first()
        #
        # print("Subset sample overlap")
        # col_overlap = set(df1.columns).intersection(set(df2.columns))
        # print("\tSample-overlap: {}".format(len(col_overlap)))
        # df1 = df1.loc[:, col_overlap]
        # df2 = df2.loc[:, col_overlap]
        # std = std.loc[col_overlap, :]
        #
        # print("Removing zero variance probes / samples")
        # df1 = df1.loc[df1.std(axis=1) != 0, df1.std(axis=0) != 0]
        # df2 = df2.loc[df2.std(axis=1) != 0, df2.std(axis=0) != 0]
        #
        # print("Subset gene overlap")
        # row_overlap = set(df1.index).intersection(set(df2.index))
        # print("\tGene-overlap: {}".format(len(row_overlap)))
        # df1 = df1.loc[row_overlap, :]
        # df2 = df2.loc[row_overlap, :]
        # genes = list(df1.index)
        #
        # print(df1)
        # print(df2)
        #
        # datasets = list(std["dataset"].unique())
        #
        # print("Correlate")
        # pearson_coefs = np.empty((df1.shape[0], len(datasets)))
        # spearman_coefs = np.empty((df1.shape[0], len(datasets)))
        # for i, index in enumerate(genes):
        #     for j, dataset in enumerate(datasets):
        #         dataset_samples = list(std.loc[std["dataset"] == dataset, "rnaseq_id"].values)
        #
        #         pearson_coef = np.nan
        #         spearman_coef = np.nan
        #         if min(np.std(df1.loc[index, dataset_samples]), np.std(df2.loc[index, dataset_samples])) > 0:
        #             pearson_coef, _ = stats.pearsonr(df1.loc[index, dataset_samples], df2.loc[index, dataset_samples])
        #             spearman_coef, _ = stats.spearmanr(df1.loc[index, dataset_samples], df2.loc[index, dataset_samples])
        #
        #         pearson_coefs[i, j] = pearson_coef
        #         spearman_coefs[i, j] = spearman_coef
        #
        # pearson_coefs_df = pd.DataFrame(pearson_coefs, index=genes, columns=datasets)
        # self.save_file(df=pearson_coefs_df, outpath=os.path.join(self.outdir, "pearson_coefs.txt.gz"))
        pearson_coefs_df = self.load_file(os.path.join(self.outdir, "pearson_coefs.txt.gz"), header=0, index_col=0)
        print("Pearson:")
        pearson_coefs_per_dataset = pearson_coefs_df.mean(axis=0)
        print(pearson_coefs_per_dataset)
        print("")

        # spearman_coefs_df = pd.DataFrame(spearman_coefs, index=genes, columns=datasets)
        # self.save_file(df=spearman_coefs_df, outpath=os.path.join(self.outdir, "spearman_coefs.txt.gz"))
        # print("Spearman:")
        # print(spearman_coefs_df.mean(axis=0))
        # print("")

        print("Plotting")

        plot_df_m = pearson_coefs_df.melt()
        print(plot_df_m)
        self.plot_boxplot(df_m=plot_df_m,
                          x="variable",
                          hue="variable",
                          palette=self.palette,
                          title=self.outfolder,
                          ylabel="pearson r",
                          coefs_per_dataset=pearson_coefs_per_dataset,
                          filename="pearson")
        exit()

        # for i in range(5):
        #     plot_df = df1.iloc[[i], ].T.merge(df2.iloc[[i], :].T, left_index=True, right_index=True).merge(std.loc[:, ["dataset"]], left_index=True, right_index=True)
        #     plot_df.columns = ["x", "y", "hue"]
        #     plot_df.dropna(inplace=True)
        #
        #     self.single_regplot(df=plot_df,
        #                         hue="hue",
        #                         palette=self.palette,
        #                         xlabel=self.name1,
        #                         ylabel=self.name2,
        #                         filename=genes[i])

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

    def single_regplot(self, df, x="x", y="y", hue=None, palette=None,
                       xlabel=None, ylabel=None, title="", filename="plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
                                       gridspec_kw={"width_ratios": [0.9, 0.1]})
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
                    handles.append(mpatches.Patch(color=palette[hue_group],
                                                  label="{} [n={:,}; r={:.2f}]".format(hue_group, group_sizes[hue_group],group_corr_coef[hue_group])))
            ax2.legend(handles=handles, loc="center", fontsize=20)

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

        # Add diagonal.
        ax1.axline((0, 0), slope=1, ls='--', color="#000000", alpha=0.15, zorder=-1)

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
        print("  > Dataset 1:")
        print("    > Path: {}".format(self.data1_path))
        print("    > Name: {}".format(self.name1))
        print("  > Dataset 2:")
        print("    > Path: {}".format(self.data2_path))
        print("    > Name: {}".format(self.name2))
        print("  > Sample to dataset: {}".format(self.std_path))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
