#!/usr/bin/env python3

"""
File:         visualise_deconvolution_covariates.py
Created:      2020/09/30
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
import scipy.stats as stats

# Local application imports.

# Metadata
__program__ = "Visualise Deconvolution Covariates"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.decon_path = getattr(arguments, 'decon')
        self.lt_path = getattr(arguments, 'lower_triangle')
        self.ut_path = getattr(arguments, 'upper_triangle')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = str(Path(__file__).parent.parent)
        self.abbreviations = {"Neuron": "neuro",
                              "Oligodendrocyte": "oligo",
                              "EndothelialCell": "endo",
                              "Microglia": "micro",
                              "Macrophage": "macro",
                              "Astrocyte": "astro"}

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
        parser.add_argument("-d",
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix")
        parser.add_argument("-lt",
                            "--lower_triangle",
                            type=str,
                            required=True,
                            help="The path to the samples for the "
                                 "lower triangle.")
        parser.add_argument("-ut",
                            "--upper_triangle",
                            type=str,
                            required=True,
                            help="The path to the samples for the"
                                 "upper triangle.")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default="png",
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        print("Loading deconvolution matrix.")
        decon_df = pd.read_csv(self.decon_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon_path),
                                      decon_df.shape))

        print("Loading lower triangle GTEx matrix.")
        lt_df = pd.read_csv(self.lt_path, sep="\t", header=0, index_col=0,
                            nrows=1)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.lt_path),
                                      lt_df.shape))
        lt_samples = list(lt_df.columns)

        print("Loading upper triangle GTEx matrix.")
        ut_df = pd.read_csv(self.ut_path, sep="\t", header=0, index_col=0,
                            nrows=1)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.ut_path),
                                      ut_df.shape))
        ut_samples = list(ut_df.columns)

        # Split the data.
        lower_df = decon_df.loc[lt_samples, :].copy()
        upper_df = decon_df.loc[ut_samples, :].copy()

        print("Lower triangle: N = {}".format(lower_df.shape[0]))
        print("Upper triangle: N = {}".format(upper_df.shape[0]))

        # Correlate and merge.
        data_df, annot_df = self.correlate(lower_df, upper_df)

        # Plot.
        self.plot(data_df, annot_df)

    @staticmethod
    def correlate(lower_df, upper_df):
        if not lower_df.columns.equals(upper_df.columns):
            raise ValueError("Dataframes are not identical.")

        corr_df = pd.DataFrame(np.nan, index=lower_df.columns, columns=lower_df.columns)
        pval_df = pd.DataFrame("", index=lower_df.columns, columns=lower_df.columns)
        mask_df = pd.DataFrame("", index=lower_df.columns, columns=lower_df.columns)
        for i, col1 in enumerate(lower_df.columns):
            for j, col2 in enumerate(lower_df.columns):
                if i < j:
                    x = upper_df.loc[:, col1]
                    y = upper_df.loc[:, col2]
                    mask = "upper"
                else:
                    x = lower_df.loc[:, col1]
                    y = lower_df.loc[:, col2]
                    mask = "lower"
                coef, p = stats.spearmanr(x, y)
                corr_df.loc[col1, col2] = coef
                pval_df.loc[col1, col2] = "{:.2e}".format(p)
                mask_df.loc[col1, col2] = mask

        print(mask_df)

        return corr_df, pval_df

    def plot(self, data_df, annot_df):
        print(data_df)
        print(annot_df)

        data_df.index = data_df.index.map(self.abbreviations)
        data_df.columns = data_df.columns.map(self.abbreviations)

        cmap = sns.diverging_palette(246, 24, as_cmap=True)

        sns.set(color_codes=True)
        g = sns.clustermap(data_df, cmap=cmap,
                           row_cluster=False, col_cluster=False,
                           yticklabels=True, xticklabels=True, square=True,
                           vmin=-1, vmax=1, annot=data_df.round(2),
                           fmt='',
                           annot_kws={"size": 16, "color": "#000000"},
                           figsize=(12, 12))
        plt.setp(
            g.ax_heatmap.set_yticklabels(
                g.ax_heatmap.get_ymajorticklabels(),
                fontsize=25, rotation=0))
        plt.setp(
            g.ax_heatmap.set_xticklabels(
                g.ax_heatmap.get_xmajorticklabels(),
                fontsize=25, rotation=45))
        # g.cax.set_visible(False)
        plt.tight_layout()
        for extension in self.extensions:
            g.savefig(os.path.join(self.outdir, "deconvolution_covariate_clustermap.{}".format(extension)))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
