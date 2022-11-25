#!/usr/bin/env python3

"""
File:         compare_bulk_and_sn_zscores_with_decon.py
Created:      2020/11/04
Last Changed: 2020/12/16
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
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from adjustText import adjust_text
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Compare Bulk and SN Z-scores with Decon"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.eqtl_type = getattr(arguments, 'type')
        self.extensions = getattr(arguments, 'extension')

        # Declare input files.
        self.bulk_eqtl_infile = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_{}/combine_eqtlprobes/eQTLprobes_combined.txt.gz".format(
            self.eqtl_type)
        self.bulk_alleles_infile = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_{}/create_matrices/genotype_alleles.txt.gz".format(
            self.eqtl_type)
        self.decon_infile = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/2020-11-20-decon-QTL/{}/cortex/decon_out/deconvolutionResults.csv".format(
            self.eqtl_type)
        self.sn_infolder = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/{}_100Perm/".format(
            self.eqtl_type)
        self.sn_filename = "eQTLsFDR-ProbeLevel.txt.gz"
        self.cell_types = [
            ("AST", "CellMapNNLS_Astrocyte", "Astrocyte", "#D55E00"),
            ("END", "CellMapNNLS_EndothelialCell", "Endothelial Cell",
             "#CC79A7"),
            ("EX", "CellMapNNLS_Neuron", "Ex. Neuron VS Neuron", "#0072B2"),
            ("IN", "CellMapNNLS_Neuron", "In. Neuron VS Neuron", "#0072B2"),
            ("MIC", "CellMapNNLS_Macrophage", "Microglia VS Macrophage",
             "#E69F00"),
            (
                "OLI", "CellMapNNLS_Oligodendrocyte", "Oligodendrocyte",
                "#009E73")]
        self.extensions = ["png", "pdf"]

        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'compare_bulk_and_sn_zscores_with_decon_{}'.format(
                                       self.eqtl_type))

        self.shared_xlim = {i: (0, 1) for i in range(len(self.cell_types))}
        self.shared_ylim = {i: (0, 1) for i in range(5)}

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        matplotlib.rcParams['pdf.fonttype'] = 42

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
        parser.add_argument("-type",
                            type=str,
                            required=False,
                            choices=["cis", "trans"],
                            default="cis",
                            help="The type of eQTLs to plot.")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        print("Loading data")

        bulk_eqtl_df = pd.read_csv(self.bulk_eqtl_infile,
                                   sep="\t",
                                   header=0,
                                   index_col=None)
        bulk_eqtl_df.reset_index(drop=False, inplace=True)
        print("\tBulk eQTL data frame: {}".format(bulk_eqtl_df.shape))
        print(bulk_eqtl_df)

        bulk_alleles_df = pd.read_csv(self.bulk_alleles_infile,
                                      sep="\t",
                                      header=0,
                                      index_col=None,
                                      nrows=bulk_eqtl_df.shape[0])
        bulk_alleles_df.columns = ["SNPName", "Alleles", "MinorAllele"]
        bulk_alleles_df["DeconAllele"] = [x.split("/")[1] for x in
                                          bulk_alleles_df["Alleles"]]
        bulk_alleles_df.reset_index(drop=False, inplace=True)
        print("\tBulk alleles data frame: {}".format(bulk_alleles_df.shape))
        print(bulk_alleles_df)

        bulk_df = bulk_eqtl_df.merge(bulk_alleles_df,
                                     left_on=["index", "SNPName"],
                                     right_on=["index", "SNPName"])
        print(bulk_df)

        decon_df = pd.read_csv(self.decon_infile,
                               sep="\t",
                               header=0,
                               index_col=0)
        print("\tDeconvolution results data frame: {}".format(decon_df.shape))
        for col in decon_df.columns:
            if col.endswith("_pvalue"):
                new_colname = col.replace("_pvalue", "_FDR")
                decon_df[new_colname] = \
                    multitest.multipletests(decon_df.loc[:, col],
                                            method='fdr_bh')[
                        1]
                print(col, decon_df[decon_df[new_colname] < 0.05].shape[0])

        print("Preprocessing deconvolution data frame")
        probe_names = []
        snp_names = []
        for index in decon_df.index:
            probe_names.append(index.split("_")[0])
            snp_names.append("_".join(index.split("_")[1:]))
        decon_df["ProbeName"] = probe_names
        decon_df["SNPName"] = snp_names
        decon_df.reset_index(drop=True, inplace=True)
        print(decon_df)

        print("Merging results")
        bulk_decon_df = bulk_df.merge(decon_df,
                                      left_on=["SNPName", "ProbeName"],
                                      right_on=["SNPName", "ProbeName"])
        del bulk_df
        bulk_decon_df.index = bulk_decon_df["SNPName"] + "_" + bulk_decon_df[
            "ProbeName"]
        print(bulk_decon_df)

        print("Visualizing")

        sns.set(rc={'figure.figsize': (len(self.cell_types) * 8, 5 * 6)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=5, ncols=len(self.cell_types),
                                 sharex='col', sharey='row')

        for col_index, (sn_ct, b_ct, title, color) in enumerate(
                self.cell_types):
            sn_df = pd.read_csv(
                os.path.join(self.sn_infolder, sn_ct, self.sn_filename),
                sep="\t",
                header=0,
                index_col=0)
            sn_df.index = sn_df["SNPName"] + "_" + sn_df["ProbeName"]
            print("\tSingle-nucleus {} eQTL data frame: {}".format(sn_ct,
                                                                   sn_df.shape))

            # test = sn_df[["OverallZScore", "FDR"]].copy()
            # test.sort_values(by="FDR", inplace=True)
            # print(test.iloc[0:50, :])
            # exit()

            merged_df = bulk_decon_df.merge(sn_df, left_index=True,
                                            right_index=True,
                                            suffixes=["_bulk", "_sn"])

            merged_df["eQTLFlipMask"] = (
                    merged_df["AlleleAssessed_sn"] == merged_df[
                "AlleleAssessed_bulk"]).replace({0: -1, 1: 1})
            merged_df["DeconFlipMask"] = (
                    merged_df["AlleleAssessed_sn"] == merged_df[
                "DeconAllele"]).replace({0: -1, 1: 1})

            beta_col = None
            for col in merged_df.columns:
                if col.startswith("Beta") and col.endswith(
                        "_{}:GT".format(b_ct)):
                    beta_col = col
                    break

            merged_df["OverallZScore_bulk_flipped"] = merged_df[
                                                          "OverallZScore_bulk"] * \
                                                      merged_df["eQTLFlipMask"]
            merged_df[
                "log_{}_flipped".format(beta_col)] = self.log_modulus_beta(
                merged_df[beta_col]) * merged_df["DeconFlipMask"]

            print("\tPrepare plot df.")
            plot_df = merged_df[["HGNCName_bulk", "OverallZScore_sn", "FDR_sn",
                                 "OverallZScore_bulk_flipped",
                                 "log_{}_flipped".format(beta_col),
                                 "{}_FDR".format(b_ct)]].copy()
            del merged_df
            plot_df["hue"] = "#808080"
            plot_df.loc[(plot_df["FDR_sn"] < 0.05) & (
                    plot_df["{}_FDR".format(b_ct)] < 0.05), "hue"] = color

            # plot_df.sort_values(by="FDR_sn", inplace=True)
            # print(plot_df.iloc[0:50, :])

            include_ylabel = False
            if col_index == 0:
                include_ylabel = True

            print("\tPlotting row 1.")
            xlim, ylim = self.plot(df=plot_df,
                      fig=fig,
                      ax=axes[0, col_index],
                      x="OverallZScore_sn",
                      y="OverallZScore_bulk_flipped",
                      xlabel="",
                      ylabel="cortex eQTL z-score",
                      title=title,
                      color=color,
                      include_ylabel=include_ylabel)
            self.update_limits(xlim, ylim, 0, col_index)

            print("\tPlotting row 2.")
            xlim, ylim = self.plot(df=plot_df.loc[plot_df["{}_FDR".format(b_ct)] < 0.05, :],
                      fig=fig,
                      ax=axes[1, col_index],
                      x="OverallZScore_sn",
                      y="OverallZScore_bulk_flipped",
                      xlabel="",
                      ylabel="cortex eQTL z-score",
                      title="",
                      color=color,
                      include_ylabel=include_ylabel)
            self.update_limits(xlim, ylim, 1, col_index)

            print("\tPlotting row 3.")
            xlim, ylim = self.plot(df=plot_df.loc[plot_df["{}_FDR".format(b_ct)] < 0.05, :],
                      fig=fig,
                      ax=axes[2, col_index],
                      x="OverallZScore_sn",
                      y="log_{}_flipped".format(beta_col),
                      xlabel="",
                      ylabel="log deconvolution beta",
                      title="",
                      color=color,
                      include_ylabel=include_ylabel)
            self.update_limits(xlim, ylim, 2, col_index)

            print("\tPlotting row 4.")
            xlim, ylim = self.plot(df=plot_df.loc[plot_df["FDR_sn"] < 0.05, :],
                      fig=fig,
                      ax=axes[3, col_index],
                      x="OverallZScore_sn",
                      y="OverallZScore_bulk_flipped",
                      facecolors="hue",
                      xlabel="",
                      ylabel="cortex eQTL z-score",
                      title="",
                      color=color,
                      include_ylabel=include_ylabel)
            self.update_limits(xlim, ylim, 3, col_index)

            print("\tPlotting row 5.")
            xlim, ylim = self.plot(df=plot_df.loc[(plot_df["FDR_sn"] < 0.05) & (plot_df["{}_FDR".format(b_ct)] < 0.05), :],
                      fig=fig,
                      ax=axes[4, col_index],
                      x="OverallZScore_sn",
                      y="log_{}_flipped".format(beta_col),
                      label="HGNCName_bulk",
                      xlabel="single-nucleus z-score",
                      ylabel="log deconvolution beta",
                      title="",
                      color=color,
                      ci=None,
                      include_ylabel=include_ylabel)
            self.update_limits(xlim, ylim, 4, col_index)

            test = plot_df.loc[(plot_df["FDR_sn"] < 0.05) & (
                    plot_df["{}_FDR".format(b_ct)] < 0.05), :]
            with pd.option_context('display.max_rows', None,
                                   'display.max_columns',
                                   None):
                print(test)

            print("")

        for (m, n), ax in np.ndenumerate(axes):
            (xmin, xmax) = self.shared_xlim[n]
            (ymin, ymax) = self.shared_ylim[m]

            xmargin = (xmax - xmin) * 0.05
            ymargin = (ymax - ymin) * 0.05

            ax.set_xlim(xmin - xmargin, xmax + xmargin)
            ax.set_ylim(ymin - ymargin, ymax + ymargin)

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir,
                                     "compare_bulk_and_sn_zscores_with_decon_{}.{}".format(
                                         self.eqtl_type, extension)))
        plt.close()

    @staticmethod
    def log_modulus_beta(series):
        s = series.copy()
        data = []
        for index, beta in s.T.iteritems():
            data.append(np.log(abs(beta) + 1) * np.sign(beta))
        new_df = pd.Series(data, index=s.index)

        return new_df

    @staticmethod
    def create_flisk_mask(df1, df1_key, df2, df2_key):
        if df1.shape[0] != df2.shape[0]:
            print("ERROR, unequal data frame sizes")
            exit()

        flip_mask = []
        count = 0
        for i in range(df1.shape[0]):
            if df1.iloc[i, :][df1_key] != df2.iloc[i, :][df2_key]:
                flip_mask.append(-1)
                count += 1
            else:
                flip_mask.append(1)
        print("\t\t{} vs {}:\t-1 = {}\t1 = {}".format(df1_key, df2_key, count,
                                                      len(flip_mask) - count))

        return flip_mask

    def plot(self, df, fig, ax, x="x", y="y", facecolors=None, label=None,
             xlabel="", ylabel="", title="", color="#000000", ci=95,
             include_ylabel=True):
        sns.despine(fig=fig, ax=ax)

        if not include_ylabel:
            ylabel = ""

        if facecolors is None:
            facecolors = "#808080"
        else:
            facecolors = df[facecolors]

        n = df.shape[0]
        coef = np.nan
        concordance = 0

        if n > 0:
            lower_quadrant = df.loc[(df[x] < 0) & (df[y] < 0), :]
            upper_quadrant = df.loc[(df[x] > 0) & (df[y] > 0), :]
            concordance = (100 / n) * (lower_quadrant.shape[0] + upper_quadrant.shape[0])

            if n > 1:
                # coef, p = stats.spearmanr(subset[x], subset[y])
                coef, p = stats.pearsonr(df[x], df[y])

            sns.regplot(x=x, y=y, data=df, ci=ci,
                        scatter_kws={'facecolors': facecolors,
                                     'edgecolors': "#808080"},
                        line_kws={"color": color},
                        ax=ax
                        )

            if label is not None:
                texts = []
                for i, point in df.iterrows():
                    texts.append(ax.text(point[x],
                                         point[y],
                                         str(point[label]),
                                         ha='center',
                                         va='center',
                                         color=color))
                adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#808080'))

        ax.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        ax.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)

        ax.text(0.5, 1.1, title,
                fontsize=18, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02, "Allelic concordance: {:.2f}% [N = {}]".format(concordance, n),
                fontsize=14, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        ax.legend(handles=[mpatches.Patch(color=color, label="r = {:.2f}".format(coef))], loc=4)

        return (df[x].min(), df[x].max()), (df[y].min(), df[y].max())

    def update_limits(self, xlim, ylim, row, col):
        row_ylim = self.shared_ylim[row]
        if ylim[0] < row_ylim[0]:
            row_ylim = (ylim[0], row_ylim[1])
        if ylim[1] > row_ylim[1]:
            row_ylim = (row_ylim[0], ylim[1])
        self.shared_ylim[row] = row_ylim

        col_xlim = self.shared_xlim[col]
        if xlim[0] < col_xlim[0]:
            col_xlim = (xlim[0], col_xlim[1])
        if xlim[1] > col_xlim[1]:
            col_xlim = (col_xlim[0], xlim[1])
        self.shared_xlim[col] = col_xlim


if __name__ == '__main__':
    m = main()
    m.start()
