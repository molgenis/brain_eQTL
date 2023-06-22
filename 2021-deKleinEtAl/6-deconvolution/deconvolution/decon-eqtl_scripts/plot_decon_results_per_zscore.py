#!/usr/bin/env python3

"""
File:         plot_decon_results_per_zscore.py
Created:      2022/01/24
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
import argparse
import math
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Plot Decon-eQTL Results per Z-score"
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
./plot_decon_results_per_zscore.py -of 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron

./plot_decon_results_per_zscore.py -of 2022-02-27-CortexEUR-cis-ForceNormalised-MAF5-NegativeToZero-DatasetAndRAMCorrected
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), "plot")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "OtherNeuron": '#2690ce',
            "Oligodendrocyte": "#009E73",
            "Astrocyte": "#D55E00",
            "Microglia": "#E69F00",
            "EndothelialCell": "#CC79A7"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading decon-eQTL data")
        data = []
        for zscore in [1, 2, 3, 4, 5, 0]:
            label = ">{}SD".format(zscore)
            if zscore == 0:
                label = "All"

            dataset_result_path = os.path.join("decon_eqtl", self.outfolder + "-GT{}SD".format(zscore), "deconvolutionResults.txt.gz")

            if os.path.exists(dataset_result_path):
                print("\tLoading results from z-score cut-off >'{}'".format(zscore))
                decon_df = self.load_file(dataset_result_path)
                decon_fdr_df = decon_df.loc[:, [x for x in decon_df.columns if x.endswith("_FDR")]]
                decon_fdr_df.columns = [x.split("_")[0] for x in decon_fdr_df.columns]

                for ct in decon_fdr_df.columns:
                    n_ieqtls = decon_fdr_df.loc[decon_fdr_df[ct] <= 0.05, :].shape[0]
                    data.append([label, ct, n_ieqtls])

        df = pd.DataFrame(data, columns=["variable", "group", "value"])
        print(df)

        print("Plotting")
        self.lineplot(df_m=df,
                      x="variable",
                      y="value",
                      hue="group",
                      palette=self.palette,
                      xlabel="z-score cut-off",
                      ylabel="#ieQTLs (FDR <0.05)",
                      filename=self.outfolder + "_lineplot",
                      outdir=self.outdir)

        sum_df = df.groupby(df["variable"]).sum()
        sum_df["variable"] = sum_df.index
        sum_df["group"] = ""
        print(sum_df)
        self.plot_barplot(
            df=sum_df,
            group_column="group",
            groups=[""],
            x="value",
            y="variable",
            palette={key: "#000000" for key in sum_df["variable"].unique()},
            xlabel="#ieQTLs (FDR<0.05)",
            ylabel="z-score cut-off",
            filename=self.outfolder + "_summed_lineplot",
            outdir=self.outdir
        )

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df


    @staticmethod
    def lineplot(df_m, x="x", y="y", hue=None, palette=None, title="",
                 xlabel="", ylabel="", filename="plot", outdir=None):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1,
                                       ncols=2,
                                       gridspec_kw={"width_ratios": [0.99, 0.01]})
        sns.despine(fig=fig, ax=ax1)

        g = sns.lineplot(data=df_m,
                         x=x,
                         y=y,
                         units=hue,
                         hue=hue,
                         palette=palette,
                         estimator=None,
                         legend=None,
                         sort=False,
                         ax=ax1)

        ax1.set_title(title,
                      fontsize=14,
                      fontweight='bold')
        ax1.set_xlabel(xlabel,
                       fontsize=10,
                       fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=10,
                       fontweight='bold')

        if palette is not None:
            handles = []
            for key, color in palette.items():
                if key in df_m[hue].values.tolist():
                    handles.append(mpatches.Patch(color=color, label=key))
            ax2.legend(handles=handles, loc="center")
        ax2.set_axis_off()

        plt.tight_layout()
        outpath = "{}.png".format(filename)
        if outdir is not None:
            outpath = os.path.join(outdir, outpath)
        fig.savefig(outpath)
        plt.close()

    def plot_barplot(self, df, group_column, groups, x="x", y="y", xlabel="",
                         ylabel="", palette=None, filename="", outdir=None):
        if df.shape[0] <= 2:
            return

        nplots = len(groups)
        ncols = math.ceil(np.sqrt(nplots))
        nrows = math.ceil(nplots / ncols)

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharey="all",
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

                g = sns.barplot(x=x,
                                y=y,
                                hue=y,
                                palette=palette,
                                dodge=False,
                                data=plot_df,
                                orient="h",
                                ax=ax)
                g.legend_.remove()

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

        plt.tight_layout()
        outpath = "{}.png".format(filename)
        if outdir is not None:
            outpath = os.path.join(outdir, outpath)
        fig.savefig(outpath)
        plt.close()
        print("\tSaved: {}".format(outpath))

    def print_arguments(self):
        print("Arguments:")
        print("  > Output folder: {}".format(self.outfolder))
        print("  > Output directory {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
