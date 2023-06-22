#!/usr/bin/env python3

"""
File:         visualise_decon_eqtl_betas.py
Created:      2021/09/07
Last Changed: 2021/10/07
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
import glob
import os
import re

# Third party imports.
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Visualise Decon-eQTL Beta's"
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
./visualise_decon_eqtl_betas.py -id /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts -if CortexEUR-cis-NormalisedMAF5-SaveAllBetas

./visualise_decon_eqtl_betas.py -if CortexEUR-cis-NormalisedMAF5-AllConfigs

./visualise_decon_eqtl_betas.py -if CortexEUR-cis-NormalisedMAF5-LimitedConfigs
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        indir = getattr(arguments, 'indir')
        infolder = getattr(arguments, 'infolder')

        # Set variables.
        if indir is None:
            indir = str(Path(__file__).parent.parent)
        self.indir = os.path.join(indir, "decon_eqtl", infolder)
        self.outdir = os.path.join(self.indir, "plot")
        for dir in [self.indir, self.outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-id",
                            "--indir",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of the input path.")
        parser.add_argument("-if",
                            "--infolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the input folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Processing nominal alternative model beta's")
        print("\tLoading data")
        betas_alt_m = self.load_matrix(os.path.join(self.indir, "betas_alternative_model.npy"))

        n_zero_betas_alt_s = np.sum((betas_alt_m == 0), axis=1)
        n_zero_betas_alt_counts_df = pd.DataFrame(dict(zip(*np.unique(n_zero_betas_alt_s, return_counts=True))), index=["x"]).T
        n_zero_betas_alt_counts_df["y"] = n_zero_betas_alt_counts_df.index
        self.plot_barplot(df=n_zero_betas_alt_counts_df,
                          xlabel="# beta's == 0 per alternative model",
                          filename="n_zero_betas_alternative_model")

        self.plot(m=betas_alt_m,
                  xlabel="all beta's (alternative model)",
                  filename="all_betas_alternative_model")
        self.plot(m=betas_alt_m[:, :5],
                  xlabel="cell type beta's (alternative model)",
                  filename="cell_type_betas_alternative_model")
        self.plot(m=betas_alt_m[:, 5:],
                  xlabel="interaction beta's (alternative model)",
                  filename="interaction_betas_alternative_model")
        del betas_alt_m
        print("")

        ########################################################################

        print("Processing nominal null model beta's")
        print("\tLoading data")
        betas_null_m = self.load_matrix(os.path.join(self.indir, "betas_null_model.npy"))

        n_zero_betas_null_s = np.sum((betas_null_m == 0), axis=1)
        n_zero_betas_null_counts_df = pd.DataFrame(dict(zip(*np.unique(n_zero_betas_null_s, return_counts=True))), index=["x"]).T
        n_zero_betas_null_counts_df["y"] = n_zero_betas_null_counts_df.index
        self.plot_barplot(df=n_zero_betas_null_counts_df,
                          xlabel="# beta's == 0 per null model",
                          filename="n_zero_betas_null_model")

        self.plot(m=betas_null_m,
                  xlabel="all beta's (null model)",
                  filename="all_betas_null_model")
        self.plot(m=betas_null_m[:, :, :5],
                  xlabel="cell type beta's (null model)",
                  filename="cell_type_betas_null_model")
        self.plot(m=betas_null_m[:, :, 5:],
                  xlabel="interaction beta's (null model)",
                  filename="interaction_betas_null_model")
        del betas_null_m
        print("")

        ########################################################################

        print("Loading permutation alternative model beta's")
        print("\tLoading data")
        perm_betas_alt_list = []
        perm_betas_alt_inpaths = glob.glob(os.path.join(self.indir, "permutation_betas_alternative_model_*"))
        perm_betas_alt_inpaths.sort(key=self.natural_keys)
        for i, perm_betas_alt_inpath in enumerate(perm_betas_alt_inpaths):
            perm_betas_alt_list.append(self.load_matrix(perm_betas_alt_inpath))
            break
        perm_betas_alt_m = np.dstack(perm_betas_alt_list)
        del perm_betas_alt_list
        print("\tShape: {}".format(perm_betas_alt_m.shape))
        print("")

        self.plot(m=perm_betas_alt_m,
                  xlabel="all permutation beta's (alternative model)",
                  filename="all_permutation_betas_alternative_model")
        self.plot(m=perm_betas_alt_m[:, :, :, :5],
                  xlabel="cell type permutation beta's (alternative model)",
                  filename="cell_type_permutation_betas_alternative_model")
        self.plot(m=perm_betas_alt_m[:, :, :, 5:],
                  xlabel="interaction permutation beta's (alternative model)",
                  filename="interaction_permutation_betas_alternative_model")
        del perm_betas_alt_m
        print("")

    @staticmethod
    def load_matrix(inpath):
        m = np.load(inpath)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      m.shape))
        return m

    @staticmethod
    def natural_keys(text):
        return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', text)]

    def plot_barplot(self, df, x="x", y="y", xlabel="", ylabel="", title="",
                     filename=""):
        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(12, 12))

        sns.despine(fig=fig, ax=ax)

        # Plot.
        g = sns.barplot(x=x, y=y, color="#000000", dodge=False, data=df, orient="h", ax=ax)

        offset = df["x"].max() / 25
        for index, row in df.iterrows():
            g.text(row["x"] + offset, row["y"], "{:,.0f}".format(row["x"]), color='#b22222', ha="center")

        ax.set_title(title,
                     fontsize=22,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\t  Saved figure: {} ".format(os.path.basename(outpath)))

    def plot(self, m, xlabel, filename):
        print("\tPre-process")
        x = m.flatten()
        x = np.abs(x[~np.isnan(x)])

        print("\tPlotting")
        self.distplot(x=x,
                      xlabel=xlabel,
                      filename=filename)
        del x

    def distplot(self, x, xlabel="", title="", filename="distribution"):
        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(12, 12))

        sns.despine(fig=fig, ax=ax)

        sns.kdeplot(x, shade=True, color="#808080", ax=ax, cut=0, zorder=-1)
        ax.axvline(x.mean(), ls='--', color="#808080", zorder=-1)

        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel("density",
                      fontsize=14,
                      fontweight='bold')
        ax.set_title(title,
                     fontsize=18,
                     fontweight='bold')

        ax.annotate(
            'N = {:,}'.format(np.size(x)),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=15,
            fontweight='bold')

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\t  Saved figure: {} ".format(os.path.basename(outpath)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
