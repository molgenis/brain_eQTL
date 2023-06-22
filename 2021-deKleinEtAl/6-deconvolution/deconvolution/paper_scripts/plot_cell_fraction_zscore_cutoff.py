#!/usr/bin/env python3

"""
File:         plot_cell_fraction_zscore_cutoff.py
Created:      2022/01/26
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
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Plot Cell Fraction Z-score Cut-off"
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
./plot_cell_fraction_zscore_cutoff.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex_EUR.txt.gz -o 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected -e png pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.zscore = getattr(arguments, 'zscore')
        self.outfile = getattr(arguments, 'outfile')
        self.extensions = getattr(arguments, 'extensions')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "plot_cell_fraction_zscore_cutoff")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
            "ENA": "#D46727",
            "LIBD_h650": "#808080",
            "GVEX": "#48B2E5",
            "NABEC-H550": "#6D743A",
            "CMC_HBCC_set3": "#0877B4",
            "UCLA_ASD": "#F36D2A",
            "CMC": "#EAE453",
            "Outlier": "#000000"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__,
                                         add_help=False)

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
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=False,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-z",
                            "--zscore",
                            type=int,
                            required=False,
                            default=None,
                            help="Zscore cut-off")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the output file")
        parser.add_argument("-e",
                            "--extensions",
                            type=str,
                            nargs="+",
                            default=["png"],
                            choices=["eps", "pdf", "pgf", "png", "ps", "raw", "rgba", "svg", "svgz"],
                            help="The output file format(s), default: ['png']")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")
        cf_df = self.load_file(self.cf_path, header=0, index_col=0)
        std_df = self.load_file(self.std_path, header=0, index_col=None)
        std_df.index = std_df["rnaseq_id"]
        sa_df = std_df[["dataset"]].copy()
        del std_df

        # Get stats.
        celltypes = cf_df.columns.tolist()
        celltypes.sort()

        # Merge.
        plot_df = cf_df.merge(sa_df, left_index=True, right_index=True)
        del cf_df, sa_df

        print("Plotting data.")
        self.plot(df=plot_df,
                  columns=celltypes,
                  hue="dataset",
                  palette=self.palette,
                  )

        if self.zscore is None:
            exit()

        print("Calculting z-score.")
        z_score_columns = []
        for ct in celltypes:
            plot_df["{} z-score".format(ct)] = (plot_df[ct] - plot_df[ct].mean()) / plot_df[ct].std()
            z_score_columns.append("{} z-score".format(ct))
        print(plot_df)
        plot_df["outlier"] = "False"
        plot_df.loc[plot_df.loc[:, z_score_columns].abs().max(axis=1) > self.zscore, "outlier"] = "True"
        print(plot_df)
        outlier_df = plot_df.loc[plot_df["outlier"] == "True", :].copy()
        print(outlier_df["dataset"].value_counts())

        print("Plotting data.")
        self.plot(df=plot_df,
                  columns=celltypes,
                  hue="outlier",
                  palette={"True": "#b22222", "False": "#000000"},
                  zscore=self.zscore
                  )

    @staticmethod
    def load_file(path, sep="\t", header=None, index_col=None, nrows=None, low_memory=True):
        if path.endswith(".pkl"):
            df = pd.read_pickle(path)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def plot(self, df, columns, hue=None, palette=None, zscore=None):
        nplots = len(columns) + 1
        ncols = math.ceil(np.sqrt(nplots))
        nrows = math.ceil(nplots / ncols)

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 # sharex='col',
                                 # sharey='row',
                                 figsize=(12 * ncols, 12 * nrows))
        sns.set(color_codes=True)

        row_index = 0
        col_index = 0
        groups_present = set()
        for i in range(ncols * nrows):
            print(i)
            if nrows == 1:
                ax = axes[col_index]
            elif ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            if i < len(columns):
                sns.despine(fig=fig, ax=ax)

                pi = columns[i]

                # Merge.
                plot_df = df[[pi, "dataset", hue]].copy()
                plot_df.columns = ["y", "dataset", "hue"]
                plot_df.dropna(inplace=True)

                # set order.
                counter = 0
                for group in plot_df["dataset"].unique():
                    mask = plot_df["dataset"] == group
                    subset = plot_df.loc[mask, :]
                    plot_df.loc[mask, "x"] = subset["y"].argsort() + counter
                    counter += np.sum(mask)
                    groups_present.add(group)

                # Plot.
                self.plot_scatterplot(ax=ax,
                                      df=plot_df,
                                      palette=palette,
                                      title=pi)
            else:
                ax.set_axis_off()

                if palette is not None and i == (nplots - 1):
                    handles = []
                    for key, value in palette.items():
                        if key in groups_present:
                            handles.append(mpatches.Patch(color=value, label=key))
                    ax.legend(handles=handles, loc=4, fontsize=25)

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        appendix = ""
        if zscore is not None:
            appendix = "-{}SD".format(zscore)

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}_scatterplot{}.{}".format(self.outfile, appendix, extension)))
        plt.close()

    @staticmethod
    def plot_scatterplot(ax, df, x="x", y="y", hue="hue", palette=None,
                     xlabel="", ylabel="", title=""):

        # Plot.
        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        palette=palette,
                        linewidth=0,
                        legend=False,
                        ax=ax)

        ax.set_title(title,
                     fontsize=40,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=20,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=20,
                      fontweight='bold')

    def print_arguments(self):
        print("Arguments:")
        print("  > Cell fraction file: {}".format(self.cf_path))
        print("  > Sample-to-dataset file: {}".format(self.std_path))
        print("  > Z-score cut-off: {}".format(self.zscore))
        print("  > Outfile: {}".format(self.outfile))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
