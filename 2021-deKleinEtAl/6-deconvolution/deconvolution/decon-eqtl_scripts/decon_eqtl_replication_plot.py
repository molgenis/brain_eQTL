#!/usr/bin/env python3

"""
File:         afr_replication.py
Created:      2021/06/24
Last Changed: 2021/09/02
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
import numpy as np
import pandas as pd
from statsmodels.stats import multitest
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

"""
Syntax:
./afr_replication.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2021-06-24-decon-QTL/CortexEUR-cis-PrimaryeQTLs/deconvolutionResults.csv -dn EUR -r /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2021-06-24-decon-QTL/CortexAFR-cis-Replication-EUR/deconvolutionResults.csv -rn AFR

./afr_replication.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/CortexEUR-cis-PrimaryeQTLs/deconvolutionResults.txt.gz -dn EUR -r /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/CortexAFR-cis-Replication-EUR/deconvolutionResults.txt.gz -rn AFR

./afr_replication.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/CortexEUR-cis-WithPermutations-StaticInteractionShuffle/deconvolutionResults.txt.gz -dn EUR -r /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/CortexAFR-cis-Replication-EUR-HalfNormalised/deconvolutionResults.txt.gz -rn AFR

./afr_replication.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/CortexEUR-cis-NormalisedMAF5-LimitedConfigs-NewProfileNoPericytes/deconvolutionResults.txt.gz -dn EUR -r /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/CortexAFR-cis-EURReplication-NormalisedMAF5-LimitedConfigs-NewProfileNoPericytes/deconvolutionResults.txt.gz -rn AFR

./afr_replication.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/CortexEUR-cis-NormalisedMAF5-LimitedConfigs-OldProfile/deconvolutionResults.txt.gz -dn EUR -r /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/CortexAFR-cis-EURReplication-NormalisedMAF5-LimitedConfigs-OldProfile/deconvolutionResults.txt.gz -rn AFR

### 2021-12-23 ###

./afr_replication.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2021-12-22-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron/deconvolutionResults.txt.gz -dn EUR -r /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-ProbesWithZeroVarianceRemoved-InhibitorySummedWithOtherNeuron/deconvolutionResults.txt.gz -rn AFR


### 2022-02-08 ###

./afr_replication.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron/deconvolutionResults_EMP_FDR.txt.gz \
    -dn EUR \
    -r /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-02-09-CortexAFR-cis-replicationOfCortexEUR-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron/deconvolutionResults.txt.gz \
    -rn AFR
"""

# Metadata
__program__ = "AFR Replication Plot"
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
        self.discovery_path = getattr(arguments, 'discovery_path')
        self.discovery_name = getattr(arguments, 'discovery_name')
        self.replication_path = getattr(arguments, 'replication_path')
        self.replication_name = getattr(arguments, 'replication_name')

        # Set variables.
        outdir = os.path.join(str(Path(__file__).parent.parent), 'decon_eqtl_replication_plot')
        self.plot_outdir = os.path.join(outdir, 'plot')
        self.data_outdir = os.path.join(outdir, 'data')
        for tmp_outdir in [self.plot_outdir, self.data_outdir]:
            if not os.path.exists(tmp_outdir):
                os.makedirs(tmp_outdir)

        self.colormap = {
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Neuron": "#2690ce",
            "Oligodendrocyte": "#009E73",
            "OPC": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Pericytes": "#808080",
            "OtherNeuron": "#2690ce"
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
                            help="show program's version number and exit")
        parser.add_argument("-d",
                            "--discovery_path",
                            type=str,
                            required=True,
                            help="The path to the discovery deconvolution "
                                 "results matrix")
        parser.add_argument("-dn",
                            "--discovery_name",
                            type=str,
                            default="one",
                            help="The name for the discovery data")
        parser.add_argument("-r",
                            "--replication_path",
                            type=str,
                            required=True,
                            help="The path to the replication deconvolution "
                                 "results matrix")
        parser.add_argument("-rn",
                            "--replication_name",
                            type=str,
                            default="one",
                            help="The name for the replication data")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        discovery_df = self.load_file(self.discovery_path, header=0, index_col=0)
        discovery_cell_types = [x.replace("_pvalue", "") for x in discovery_df if x.endswith("_pvalue")]
        replication_df = self.load_file(self.replication_path, header=0, index_col=0)
        replication_cell_types = [x.replace("_pvalue", "") for x in discovery_df if x.endswith("_pvalue")]

        if len(set(discovery_cell_types).symmetric_difference(set(replication_cell_types))) != 0:
            print("Unequal cell types used.")
            exit()
        cell_types = discovery_cell_types
        del discovery_cell_types, replication_cell_types

        print("Transfer betas to log scale")
        discovery_df = self.log_transform(discovery_df)
        replication_df = self.log_transform(replication_df)
        print(replication_df)

        print("Pre-processing")
        discovery_beta_df_m = discovery_df.loc[:, [col for col in discovery_df.columns if "Beta" in col]].melt()
        discovery_beta_df_m = self.add_boxplot_columns(discovery_beta_df_m)
        discovery_beta_df_m["file"] = "discovery"
        replication_beta_df_m = discovery_df.loc[:, [col for col in replication_df.columns if "Beta" in col]].melt()
        replication_beta_df_m = self.add_boxplot_columns(replication_beta_df_m)
        replication_beta_df_m["file"] = "replication"
        beta_df_m = pd.concat([discovery_beta_df_m, replication_beta_df_m], axis=0)
        beta_df_m = beta_df_m.sort_values(by="value", ascending=False)
        self.plot_boxplot(df_m=beta_df_m,
                          x="type",
                          col="file",
                          hue="hue",
                          filename="")
        del discovery_beta_df_m, replication_beta_df_m

        # print("Calculating BH-FDR for discovery")
        # discovery_df = self.add_fdr_columns(discovery_df)
        # print(discovery_df)

        print("Merging")
        combined_df = None
        for cell_type in cell_types:
            print("\tProcessing cell type: {}".format(cell_type))

            # Find the interaction beta columns for this cell type.
            discovery_beta_column = None
            for col in discovery_df:
                if cell_type in col and col.endswith(":GT"):
                    discovery_beta_column = col
                    break

            replication_beta_column = None
            for col in replication_df:
                if cell_type in col and col.endswith(":GT"):
                    replication_beta_column = col
                    break

            discovery_subset_df = discovery_df[["{}_pvalue".format(cell_type), "{}_FDR".format(cell_type), discovery_beta_column]].copy()
            replication_subset_df = replication_df[["{}_pvalue".format(cell_type), replication_beta_column]].copy()

            ct_df = discovery_subset_df.merge(replication_subset_df, left_index=True, right_index=True)
            ct_df.columns = ["discovery p-value", "discovery FDR", "discovery beta", "replication p-value", "replication beta"]
            ct_df.insert(0, "cell type", cell_type)
            ct_df.reset_index(drop=False, inplace=True)

            # Adding the replication FDR.
            ct_df["replication FDR"] = np.nan
            mask = (ct_df["discovery FDR"] <= 0.05).to_numpy()
            print("\t  Discovery N-ieqtls: {}".format(np.sum(mask)))
            ct_df.loc[mask, "replication FDR"] = multitest.multipletests(ct_df.loc[mask, "replication p-value"], method='fdr_bh')[1]
            print("\t  Replication N-ieqtls: {}".format(ct_df.loc[ct_df["replication FDR"] <= 0.05, :].shape[0]))

            if combined_df is None:
                combined_df = ct_df
            else:
                combined_df = pd.concat([combined_df, ct_df], axis=0)

            del discovery_subset_df, replication_subset_df

        print("Saving")
        print(combined_df)
        self.save_file(df=combined_df, outpath=os.path.join(self.data_outdir, "replication_table.txt.gz"))
        opposite_df = combined_df.loc[(combined_df["discovery FDR"] <= 0.05) & (combined_df["replication FDR"] <= 0.05) & (((combined_df["discovery beta"] < 0) & (combined_df["replication beta"] > 0)) | ((combined_df["discovery beta"] > 0) & (combined_df["replication beta"] < 0))), :]
        print(opposite_df)
        self.save_file(df=opposite_df, outpath=os.path.join(self.data_outdir, "replication_opposite_efects_table.txt.gz"))

        print("Plotting")
        self.plot_scatterplot(
            df=combined_df,
            group_column="cell type",
            x="discovery beta",
            y="replication beta",
            xlabel="{} (log beta)".format(self.discovery_name),
            ylabel="{} (log beta)".format(self.replication_name),
            palette=self.colormap,
            filename="{}_vs_{}_all".format(self.discovery_name, self.replication_name)
        )
        self.plot_scatterplot(
            df=combined_df.loc[combined_df["discovery FDR"] <= 0.05, :],
            group_column="cell type",
            x="discovery beta",
            y="replication beta",
            xlabel="{} (log beta)".format(self.discovery_name),
            ylabel="{} (log beta)".format(self.replication_name),
            palette=self.colormap,
            filename="{}_signif_vs_{}".format(self.discovery_name,
                                              self.replication_name)
            )
        self.plot_scatterplot(
            df=combined_df.loc[combined_df["replication FDR"] <= 0.05, :],
            group_column="cell type",
            x="discovery beta",
            y="replication beta",
            xlabel="{} (log beta)".format(self.discovery_name),
            ylabel="{} (log beta)".format(self.replication_name),
            palette=self.colormap,
            filename="{}_vs_{}_signif".format(self.discovery_name,
                                              self.replication_name)
            )
        self.plot_scatterplot(
            df=combined_df.loc[(combined_df["discovery FDR"] <= 0.05) & (combined_df["replication FDR"] <= 0.05), :],
            group_column="cell type",
            x="discovery beta",
            y="replication beta",
            xlabel="{} (log beta)".format(self.discovery_name),
            ylabel="{} (log beta)".format(self.replication_name),
            palette=self.colormap,
            filename="{}_vs_{}_both_signif".format(self.discovery_name, self.replication_name)
        )

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
    def log_transform(df):
        for col in df.columns:
            if "Beta" in col:
                beta_a = df.loc[:, col].to_numpy()
                pos_beta_mask = beta_a > 0
                neg_beta_mask = beta_a < 0

                log_beta_a = np.zeros_like(beta_a)
                log_beta_a[pos_beta_mask] = np.log(beta_a[pos_beta_mask] + 1)
                log_beta_a[neg_beta_mask] = np.log(beta_a[neg_beta_mask] * -1 + 1) * -1

                df[col] = log_beta_a

        return df

    def add_boxplot_columns(self, df):
        df["hue"] = np.nan
        for key in self.colormap.keys():
            df.loc[df["variable"].str.contains(key), "hue"] = key

        df["type"] = "beta"
        df.loc[df["variable"].str.contains(":GT"), "type"] = "interaction beta"
        return df

    @staticmethod
    def is_outlier(s):
        lower_limit = s.mean() - (s.std() * 3)
        upper_limit = s.mean() + (s.std() * 3)
        return ~s.between(lower_limit, upper_limit)

    @staticmethod
    def add_fdr_columns(df):
        for col in df.columns:
            if col.endswith("_pvalue"):
                df[col.replace("_pvalue", "_FDR")] = multitest.multipletests(df.loc[:, col], method='fdr_bh')[1]

        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def plot_boxplot(self, df_m, x="x", y="value", col=None, hue=None,
                     palette=None, xlabel="", ylabel="", title="", filename=""):
        cols = [None]
        if col is not None:
            cols = df_m[col].unique().tolist()
            cols.sort()

        ncols = len(cols) + 1
        nrows = 1

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 figsize=(12 * ncols, 12 * nrows),
                                 gridspec_kw={"width_ratios": [0.99 / len(cols)] * len(cols) + [0.01]}
                                 )
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

            if i < len(cols):
                if col is not None:
                    subset = df_m.loc[df_m[col] == cols[i], :]
                else:
                    subset = df_m

                sns.despine(fig=fig, ax=ax)

                sns.violinplot(x=x,
                               y=y,
                               hue=hue,
                               data=subset,
                               palette=palette,
                               ax=ax)

                plt.setp(ax.collections, alpha=.75)

                sns.boxplot(x=x,
                            y=y,
                            hue=hue,
                            data=subset,
                            whis=np.inf,
                            color="white",
                            ax=ax)

                ax.set_title(title,
                             fontsize=20,
                             fontweight='bold')
                ax.set_ylabel(ylabel,
                              fontsize=20,
                              fontweight='bold')
                ax.set_xlabel(xlabel,
                              fontsize=20,
                              fontweight='bold')

                ax.get_legend().remove()

                ax.tick_params(axis='both', which='major', labelsize=14)
            else:
                ax.set_axis_off()

                if palette is not None:
                    handles = []
                    for label, color in palette.items():
                        handles.append(mpatches.Patch(color=color, label=label))
                    ax.legend(handles=handles, loc=4)

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        fig.savefig(os.path.join(self.plot_outdir, "decon_eqtl_betas{}.png".format( filename)))
        plt.close()

    def plot_scatterplot(self, df, group_column, x="x", y="y", xlabel="",
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

        outpath = os.path.join(self.plot_outdir, "decon_eqtl_replication_{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved: {}".format(outpath))

    def print_arguments(self):
        print("Arguments:")
        print("  > Discovery:")
        print("    > File path: {}".format(self.discovery_path))
        print("    > Name: {}".format(self.discovery_name))
        print("  > Replication:")
        print("    > File path: {}".format(self.replication_path))
        print("    > Name: {}".format(self.replication_name))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
