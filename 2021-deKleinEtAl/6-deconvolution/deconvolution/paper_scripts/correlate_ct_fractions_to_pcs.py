#!/usr/bin/env python3

"""
File:         correlate_ct_fractions_to_pcs.py
Created:      2021/02/02
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
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Correlate Cell Type Fractions to PCs"
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
./correlate_ct_fractions_to_pcs.py \
    -cf ../../2020-03-12-deconvolution/partial_deconvolution/NEW_PROFILE_NOPERICYTES_CORTEX_EUR_TMM_LOG2/IHC_0CPM_LOG2_FILTERED_CC/deconvolution.txt.gz \
    -pc ../data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.100PrincipalComponents.pkl \
    -e pdf

./correlate_ct_fractions_to_pcs.py \
    -cf ../../2020-11-10-decon-optimizer/ftestAsPvalue/cycle0/optimized_cell_fractions.txt.gz \
    -pc ../data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.100PrincipalComponents.pkl \
    -e pdf

./correlate_ct_fractions_to_pcs.py \
    -cf ../../2020-11-10-decon-optimizer/ftestAsPvalue/cycle0/optimized_cell_fractions.txt.gz \
    -pc ../data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.100PrincipalComponents.pkl \
    -e pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.pc_path = getattr(arguments, 'principal_components')
        self.n_pc = getattr(arguments, 'n_principal_components')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'correlate_ct_fractions_to_pcs')
        self.palette = {
            "cortex": "#0072B2",
            "cerebellum": "#D55E00",
            "basalganglia": "#009E73",
            "spinalcord": "#56B4E9",
            "hypothalamus": "#E69F00",
            "hippocampus": "#F0E442",
            "amygdala": "#CC79A7",
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00"
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fracttions matrix")
        parser.add_argument("-pc",
                            "--principal_components",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the principal components.")
        parser.add_argument("-npc",
                            "--n_principal_components",
                            type=int,
                            required=False,
                            default=50,
                            help="The number of principal components to plot. "
                                 "Default = 50.")
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
        self.print_arguments()

        print("Loading data.")
        cf_df = self.load_file(self.cf_path)
        print(cf_df)
        pc_df = self.load_file(self.pc_path)
        print(pc_df)

        print("Plot PCA")
        self.plot_pca(df=pc_df)
        self.plot_cf(df=cf_df)

        print("Preprocessing data.")
        overlap_samples = set(cf_df.index).intersection(set(pc_df.index))
        print("\tOverlap: N = {}".format(len(overlap_samples)))
        cc_df = cf_df.loc[overlap_samples, :]
        pc_df = pc_df.loc[overlap_samples, :]

        print("Plotting correlation scatterplot.")
        self.plot_corr_scatterplots(df_row=cc_df, df_col=pc_df)

        if self.n_pc > pc_df.shape[1]:
            self.n_pc = pc_df.shape[1]
        pc_df = pc_df.iloc[:, 0:self.n_pc]

        print("Correlate data.")
        corr_coefficients = []
        pvalues = []
        for cell_type, cf_data in cc_df.T.iterrows():
            print("\tProcessing '{}'".format(cell_type))
            ct_coefficients = []
            ct_pvalues = []
            for i, (pc, pc_data) in enumerate(pc_df.T.iterrows()):
                coef, p = stats.spearmanr(cf_data, pc_data)
                ct_coefficients.append(coef)
                ct_pvalues.append(p)
            corr_coefficients.append(ct_coefficients)
            pvalues.append(ct_pvalues)
        corr_coef_df = pd.DataFrame(corr_coefficients, index=cc_df.columns, columns=pc_df.columns)
        print(corr_coef_df)
        pvalue_df = pd.DataFrame(pvalues, index=cc_df.columns, columns=pc_df.columns)
        print(pvalue_df)

        print("Saving data.")
        corr_coef_df.to_csv(os.path.join(self.outdir, ".".join(os.path.basename(self.pc_path).split(".")[:-1]) + ".CorrelationCoefficients.txt.gz"),
                            sep="\t", compression="gzip", header=True, index=True)
        pvalue_df.to_csv(os.path.join(self.outdir, ".".join(os.path.basename(self.pc_path).split(".")[:-1]) + ".CorrelationCoefficientPvalues.txt.gz"),
                         sep="\t", compression="gzip", header=True, index=True)

        print("Plotting data.")
        self.create_clustermap(corr_coef_df)


    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None, low_memory=True):
        if path.endswith(".pkl"):
            df = pd.read_pickle(path)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def create_clustermap(self, df):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           yticklabels=True, xticklabels=True,
                           row_cluster=True, col_cluster=False,
                           dendrogram_ratio=(.1, .1),
                           figsize=(12, 9))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=10))
        g.fig.subplots_adjust(bottom=0.05, top=0.7)
        plt.tight_layout()
        for extension in self.extensions:
            g.savefig(os.path.join(self.outdir, "cf_vs_pcs_clustermap.{}".format(extension)))
        plt.close()

    def plot_corr_scatterplots(self, df_row, df_col):
        max_cols = 5

        # Create combinations.
        nrows = df_row.shape[1]
        ncols = min(df_col.shape[1], max_cols)

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 figsize=(12 * ncols, 9 * nrows),
                                 sharex='col',
                                 sharey='row')

        for i, (y_col, y_data) in enumerate(df_row.T.iterrows()):
            for j, (x_col, x_data) in enumerate(df_col.T.iterrows()):
                if j >= max_cols:
                    break

                ax = axes[i, j]

                sns.despine(fig=fig, ax=ax)

                plot_df = pd.concat([x_data, y_data], axis=1)

                coef, _ = stats.spearmanr(plot_df[y_col], plot_df[x_col])

                sns.regplot(x=x_col, y=y_col, data=plot_df, ci=95,
                            scatter_kws={'facecolors': "#000000",
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": "#b22222",
                                      'linewidth': 5},
                            ax=ax)

                ax.annotate(
                    'r = {:.2f}'.format(coef),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color="#000000",
                    alpha=1,
                    fontsize=30,
                    fontweight='bold')

                xlabel = ""
                if i == (nrows - 1):
                    xlabel = x_col
                ax.set_xlabel(xlabel,
                              fontsize=30,
                              fontweight='bold')

                ylabel = ""
                if j == 0:
                    ylabel = y_col
                ax.set_ylabel(ylabel,
                              fontsize=30,
                              fontweight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "cf_vs_pcs_regression_plot.{}".format(extension)))
        plt.close()

    def plot_pca(self, df):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, gridspec_kw={"width_ratios": [0.9, 0.1]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        x = "PC1"
        y = "PC2"

        sns.scatterplot(x=x,
                        y=y,
                        data=df,
                        legend=False,
                        ax=ax1)

        ax1.set_xlabel(x,
                       fontsize=30,
                       fontweight='bold')
        ax1.set_ylabel(y,
                       fontsize=30,
                       fontweight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "pca_plot.{}".format(extension)))
        plt.close()

    def plot_cf(self, df):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, gridspec_kw={"width_ratios": [0.9, 0.1]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        x = "Neuron"
        y = "Astrocyte"

        coef, _ = stats.spearmanr(df[y], df[x])

        sns.regplot(x=x,
                    y=y,
                    data=df,
                    ax=ax1)

        ax1.set_xlabel(x,
                       fontsize=30,
                       fontweight='bold')
        ax1.set_ylabel(y,
                       fontsize=30,
                       fontweight='bold')

        ax1.annotate(
             'r = {:.2f}'.format(coef),
             xy=(0.03, 0.94),
             xycoords=ax1.transAxes,
             color="#000000",
             alpha=1,
             fontsize=30,
             fontweight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "neuron_vs_astro.{}".format(extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Cell fractions path: {}".format(self.cf_path))
        print("  > Principal components path: {}".format(self.pc_path))
        print("  > N Principal components: {}".format(self.n_pc))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
