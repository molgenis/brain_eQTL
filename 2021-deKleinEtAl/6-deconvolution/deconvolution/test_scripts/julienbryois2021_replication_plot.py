#!/usr/bin/env python3

"""
File:         julienbryois2021_replication_plot.py
Created:      2022/01/12
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
./julienbryois2021_replication_plot.py -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2021-12-22-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron/deconvolutionResults.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_alleles.txt.gz -dn MetaBrain_Decon-eQTL

./julienbryois2021_replication_plot.py -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-19-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron-ExprMatrixCovariatesRemovedOLS-NegativeToZero/deconvolutionResults.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_alleles.txt.gz -dn MetaBrain_Decon-eQTL

./julienbryois2021_replication_plot.py -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-20-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron-ExprMatrixCovariatesRemovedOLS-ShiftedPositive/deconvolutionResults.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_alleles.txt.gz -dn MetaBrain_Decon-eQTL

./julienbryois2021_replication_plot.py -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2021-12-22-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron-GT4SD/deconvolutionResults.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz -dn MetaBrain_Decon-eQTL
"""

# Metadata
__program__ = "Julien Bryois et al. 2021 Replication Plot"
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
        self.discovery_data_path = getattr(arguments, 'discovery_data_path')
        self.discovery_alleles_path = getattr(arguments, 'discovery_alleles_path')
        self.discovery_name = getattr(arguments, 'discovery_name')
        self.replication_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/julienbryois2021"
        self.replication_folders = ["Astrocytes", "EndothelialCells", "ExcitatoryNeurons", "InhibitoryNeurons", "Microglia", "Oligodendrocytes", "Oligodendrocytes", "OPCsCOPs", "Pericytes"]
        self.replication_name = "JulienBryois2021"

        # Set variables.
        outdir = os.path.join(str(Path(__file__).parent.parent), 'julienbryois2021_replication_plot')
        self.outdir = os.path.join(outdir, 'plot')
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
        parser.add_argument("-dd",
                            "--discovery_data_path",
                            type=str,
                            required=True,
                            help="The path to the discovery deconvolution "
                                 "results matrix")
        parser.add_argument("-da",
                            "--discovery_alleles_path",
                            type=str,
                            required=True,
                            help="The path to the discovery alleles matrix")
        parser.add_argument("-dn",
                            "--discovery_name",
                            type=str,
                            default="one",
                            help="The name for the discovery data")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        discovery_df = self.load_file(self.discovery_data_path, header=0, index_col=None)
        discovery_df["SNPName"] = ["_".join(x.split("_")[1:]) for x in discovery_df.iloc[:, 0]]
        discovery_cell_types = [x.replace("_pvalue", "") for x in discovery_df if x.endswith("_pvalue")]

        discovery_alleles_df = self.load_file(self.discovery_alleles_path, header=0, index_col=None)
        discovery_alleles_df.columns = ["SNPName", "Alleles", "MinorAllele"]
        discovery_alleles_df["DeconAllele"] = [x.split("/")[1] for x in discovery_alleles_df["Alleles"]]
        snp_to_da_dict = dict(zip(discovery_alleles_df["SNPName"], discovery_alleles_df["DeconAllele"]))

        discovery_df["DeconAllele"] = discovery_df["SNPName"].map(snp_to_da_dict)
        discovery_df.index = ["{}_{}".format(x.split(".")[0], x.split(":")[2]) for x in discovery_df.iloc[:, 0]]
        discovery_df = discovery_df.iloc[:, 1:]
        discovery_df.reset_index(drop=False, inplace=True)

        discovery_dfm_list = []
        for discovery_ct in discovery_cell_types:
            discovery_fdr_column = "{}_FDR".format(discovery_ct)
            discovery_beta_column = None
            for col in discovery_df:
                if discovery_ct in col and col.endswith(":GT"):
                    discovery_beta_column = col
                    break

            discovery_dfm = discovery_df.loc[:, ["index", discovery_beta_column, discovery_fdr_column, "DeconAllele"]].copy()
            discovery_dfm.columns = ["index", "discovery beta", "discovery fdr", "discovery allele"]
            discovery_dfm["cell type"] = discovery_ct
            discovery_dfm_list.append(discovery_dfm)

        discovery_dfm = pd.concat(discovery_dfm_list, axis=0)
        print(discovery_dfm)

        cell_type_trans_dict = {"Astrocytes": "Astrocyte",
                                "EndothelialCells": "EndothelialCell",
                                "ExcitatoryNeurons": "Excitatory",
                                "InhibitoryNeurons": "Inhibitory",
                                "Microglia": "Microglia",
                                "Oligodendrocytes": "Oligodendrocyte",
                                "OPCsCOPs": "OPCs...COPs",
                                "Pericytes": "Pericytes"}

        replication_df = self.load_file(os.path.join(self.replication_path, "JulienBryois2021SummaryStats.txt.gz"), header=0, index_col=0, nrows=None)
        replication_cell_types = [x.replace(" p-value", "") for x in replication_df if x.endswith(" p-value")]
        replication_df.reset_index(drop=False, inplace=True)

        replication_dfm_list = []
        for replication_ct in replication_cell_types:
            replication_pvalue_column = "{} p-value".format(replication_ct)
            replication_beta_column = "{} beta".format(replication_ct)

            replication_dfm = replication_df.loc[:, ["index", replication_beta_column, replication_pvalue_column, "effect_allele"]].copy()
            replication_dfm.columns = ["index", "replication beta", "replication pvalue", "replication allele"]
            replication_dfm["cell type"] = cell_type_trans_dict[replication_ct]
            replication_dfm_list.append(replication_dfm)

        replication_dfm = pd.concat(replication_dfm_list, axis=0)
        print(replication_dfm)

        print("Combine data")
        df = discovery_dfm.merge(replication_dfm, on=["index", "cell type"])
        print(df)

        print("Flip beta's")
        df["flip"] = (df["discovery allele"] == df["replication allele"]).replace({False: -1, True: 1})
        df["replication beta flipped"] = df["replication beta"] * df["flip"]
        print(df)

        print("Transfer betas to log scale")
        df = self.log_transform(df)
        print(df)

        print("Plotting")
        self.plot_scatterplot(
            df=df,
            group_column="cell type",
            x="discovery beta",
            y="replication beta flipped",
            xlabel="{} (log beta)".format(self.discovery_name),
            ylabel="{} (log beta)".format(self.replication_name),
            palette=self.colormap,
            filename="{}_vs_{}_all".format(self.discovery_name, self.replication_name)
        )
        self.plot_scatterplot(
            df=df.loc[df["discovery fdr"] <= 0.05, :],
            group_column="cell type",
            x="discovery beta",
            y="replication beta flipped",
            xlabel="{} (log beta)".format(self.discovery_name),
            ylabel="{} (log beta)".format(self.replication_name),
            palette=self.colormap,
            filename="{}_signif_vs_{}".format(self.discovery_name,
                                              self.replication_name)
            )
        self.plot_scatterplot(
            df=df.loc[df["replication pvalue"] <= 0.05, :],
            group_column="cell type",
            x="discovery beta",
            y="replication beta flipped",
            xlabel="{} (log beta)".format(self.discovery_name),
            ylabel="{} (log beta)".format(self.replication_name),
            palette=self.colormap,
            filename="{}_vs_{}_signif".format(self.discovery_name,
                                              self.replication_name)
            )
        self.plot_scatterplot(
            df=df.loc[(df["discovery fdr"] <= 0.05) & (df["replication pvalue"] <= 0.05), :],
            group_column="cell type",
            x="discovery beta",
            y="replication beta flipped",
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
            if "beta" in col:
                beta_a = df.loc[:, col].to_numpy()
                pos_beta_mask = beta_a > 0
                neg_beta_mask = beta_a < 0

                log_beta_a = np.zeros_like(beta_a)
                log_beta_a[pos_beta_mask] = np.log(beta_a[pos_beta_mask] + 1)
                log_beta_a[neg_beta_mask] = np.log(beta_a[neg_beta_mask] * -1 + 1) * -1

                df[col] = log_beta_a

        return df

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
                plot_df = df.loc[df[group_column] == groups[i], :].copy()
                plot_df.dropna(inplace=True)

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

        outpath = os.path.join(self.outdir, "julienbryois2021_replication_{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved: {}".format(outpath))

    def print_arguments(self):
        print("Arguments:")
        print("  > Discovery:")
        print("    > Data path: {}".format(self.discovery_data_path))
        print("    > Alleles path: {}".format(self.discovery_alleles_path))
        print("    > Name: {}".format(self.discovery_name))
        print("  > Replication:")
        print("    > File path: {}".format(self.replication_path))
        print("    > Name: {}".format(self.replication_name))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
