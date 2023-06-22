#!/usr/bin/env python3

"""
File:         plot_cell_fraction_per_disease_state.py
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
__program__ = "Plot Cell Fraction per disease state"
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
./plot_cell_fraction_per_disease_state.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex_EUR.txt.gz -p /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt -o 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected -e png pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.pheno_path = getattr(arguments, 'phenotype')
        self.outfile = getattr(arguments, 'outfile')
        self.extensions = getattr(arguments, 'extensions')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "cell_fraction_per_disease_state")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

        self.datasets = ["AMPAD-MAYO-V2", "AMPAD-MSBB-V2", "AMPAD-ROSMAP-V2"]

        self.palette = {
            'Alzheimer disease': '#000000',
            'Non-Neurological Control': '#D55E00'
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
        parser.add_argument("-p",
                            "--phenotype",
                            type=str,
                            required=False,
                            help="The path to the phenotype matrix.")
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
        pheno_df = self.load_file(self.pheno_path, header=0, index_col=None, low_memory=False)

        print("### Step2 ###")
        print("Pre-processing data frame")
        plot_df = cf_df.merge(std_df, left_index=True, right_on="rnaseq_id", how="inner")
        plot_df = plot_df.merge(pheno_df[["rnaseq_id", "reannotated_diangosis"]], on="rnaseq_id", how="inner")
        plot_df = plot_df.loc[(plot_df["dataset"].isin(self.datasets)) & (plot_df["reannotated_diangosis"].isin(["Alzheimer disease", "Non-Neurological Control"])), :]
        print(plot_df)

        # Get stats.
        celltypes = cf_df.columns.tolist()
        celltypes.sort()

        plot_dfm = plot_df.melt(id_vars=["dataset", "reannotated_diangosis"], value_vars=celltypes)
        print(plot_dfm)

        print("### Step3 ###")
        print("Plotting data.")
        self.plot_boxplot(df_m=plot_dfm,
                          row="dataset",
                          col="variable",
                          x="reannotated_diangosis",
                          y="value",
                          col_order=celltypes,
                          row_order=self.datasets,
                          palette=self.palette
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

    def plot_boxplot(self, df_m, row, col, x, y, col_order, row_order, palette):
        ncols = len(col_order)
        nrows = len(row_order)

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex="all",
                                 figsize=(5 * ncols, 12 * nrows)
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

            df = df_m.loc[(df_m[row] == row_order[row_index]) & (
                        df_m[col] == col_order[col_index]), :]

            sns.despine(fig=fig, ax=ax)

            sns.violinplot(x=x,
                           y=y,
                           hue=x,
                           data=df,
                           palette=palette,
                           cut=0,
                           dodge=False,
                           ax=ax)

            plt.setp(ax.collections, alpha=.75)

            sns.boxplot(x=x,
                        y=y,
                        hue=x,
                        data=df,
                        dodge=False,
                        whis=np.inf,
                        color="white",
                        ax=ax)

            tmp_title = ""
            if row_index == 0:
                tmp_title = col_order[col_index]
            ax.set_title(tmp_title,
                         fontsize=20,
                         fontweight='bold')

            tmp_ylabel = ""
            if col_index == 0:
                tmp_ylabel = row_order[row_index]
            ax.set_ylabel(tmp_ylabel,
                          fontsize=20,
                          fontweight='bold')

            ax.set_xlabel("",
                          fontsize=20,
                          fontweight='bold')

            if ax.get_legend() is not None:
                ax.get_legend().remove()

            ax.tick_params(axis='both', which='major', labelsize=14)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}_CFPerDiseaseState.{}".format(self.outfile, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Cell fraction file: {}".format(self.cf_path))
        print("  > Sample-to-dataset file: {}".format(self.std_path))
        print("  > Phenotype file: {}".format(self.pheno_path))
        print("  > Outfile: {}".format(self.outfile))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
