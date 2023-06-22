#!/usr/bin/env python3

"""
File:         plot_cell_fraction_per_cohort2.py
Created:      2021/11/29
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
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Plot Cell Fraction per Cohort 2"
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


./plot_cell_fraction_per_cohort2.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex_EUR.txt.gz -h dataset -o 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected -e png pdf

./plot_cell_fraction_per_cohort2.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table_complete.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex_EUR.txt.gz -h dataset -o 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected-Complete -e png pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.hue = getattr(arguments, 'hue')
        self.outfile = getattr(arguments, 'outfile')
        self.extensions = getattr(arguments, 'extensions')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "cell_fraction_per_cohort2")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Color map.
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

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
        parser.add_argument("-h",
                            "--hue",
                            type=str,
                            default="dataset",
                            choices=["dataset", "ethnicity"],
                            help="")
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

        print("### Step1 ###")
        print("Loading data.")
        cf_df = self.load_file(self.cf_path, header=0, index_col=0)
        std_df = self.load_file(self.std_path, header=0, index_col=None)

        print("### Step2 ###")
        print("Pre-processing data frame")
        # Merge.
        plot_df = cf_df.merge(std_df, left_index=True, right_on="rnaseq_id", how="inner")

        # Get stats.
        celltypes = cf_df.columns.tolist()

        hue_sample_counts = list(zip(*np.unique(plot_df[self.hue], return_counts=True)))
        hue_sample_counts.sort(key=lambda x: -x[1])
        x_order = [csc[0] for csc in hue_sample_counts]
        x_order.sort()
        hue_sample_counts = {key:value for key, value in hue_sample_counts}

        averages = plot_df.loc[:, celltypes].mean(axis=0)
        averages.sort_values(inplace=True, ascending=False)
        col_order = [index for index, _ in averages.iteritems()]

        # Melt.
        plot_df_m = plot_df.melt(id_vars=[self.hue], value_vars=celltypes)
        del plot_df

        print("### Step3 ###")
        print("Plotting data.")
        self.plot_boxplot(df_m=plot_df_m,
                          x=self.hue,
                          col="variable",
                          palette=self.palette,
                          cut=0,
                          x_order=x_order,
                          col_order=col_order,
                          hue_counts=hue_sample_counts,
                          ylabel="cell fraction %")

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

    def plot_boxplot(self, df_m, x="x", y="value", col=None, hue=None,
                     palette=None, cut=None, x_order=None, col_order=None,
                     hue_counts=None, xlabel="", ylabel="", filename=""):
        cols = [None]
        if col is not None:
            cols = df_m[col].unique().tolist()
            cols.sort()

            if col_order is None:
                col_order = cols

        if x_order is None:
            x_order = df_m[x].unique().tolist()
            x_order.sort()

        ngroups = len(cols) + 1
        ncols = int(np.ceil(np.sqrt(ngroups)))
        nrows = int(np.ceil(ngroups / ncols))

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex="all",
                                 figsize=(12 * ncols, 12 * nrows)
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

            if i < len(col_order):
                title = ""
                if col is not None:
                    subset = df_m.loc[df_m[col] == col_order[i], :]
                    title = col_order[i]
                else:
                    subset = df_m

                sns.despine(fig=fig, ax=ax)

                sns.violinplot(x=x,
                               y=y,
                               hue=hue,
                               data=subset,
                               order=x_order,
                               palette=palette,
                               cut=cut,
                               ax=ax)

                plt.setp(ax.collections, alpha=.75)

                sns.boxplot(x=x,
                            y=y,
                            hue=hue,
                            data=subset,
                            order=x_order,
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

                if ax.get_legend() is not None:
                    ax.get_legend().remove()

                ax.tick_params(axis='both', which='major', labelsize=14)
                ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

                ax.axhline(subset[y].mean(), ls='--', lw=5, color="#000000", zorder=-2)
            else:
                ax.set_axis_off()

                if row_index == (nrows - 1) and col_index == (ncols - 1):
                    if palette is not None:
                        handles = []
                        for group in x_order:
                            label = group
                            if hue_counts is not None:
                                label = "{} [n={:,}]".format(group, hue_counts[group])
                            handles.append(mpatches.Patch(color=self.palette[group], label=label))
                        ax.legend(handles=handles, loc=4, fontsize=20)

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}_CFPer{}{}.{}".format(self.outfile, self.hue, filename, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Cell fraction file: {}".format(self.cf_path))
        print("  > Sample-to-dataset file: {}".format(self.std_path))
        print("  > Hue: {}".format(self.hue))
        print("  > Outfile: {}".format(self.outfile))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
