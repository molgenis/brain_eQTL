#!/usr/bin/env python3

"""
File:         compare_pvalues.py
Created:      2022/04/06
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
import argparse
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Compare Pvalues"
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
./compare_pvalues.py -h

./compare_pvalues.py \
    -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/interaction_mapper/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/Alzheimerdisease_InteractionResults.txt.gz \
    -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/interaction_mapper/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-AMPAD/Alzheimerdisease_InteractionResults.txt.gz \
    -n1 all_samples \
    -n2 AMPAD_samples \
    -t CortexEUR_and_AFR_noENA_trans_0PCs_ADStatusInteraction \
    -o 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-all_vs_AMPAD
    
./compare_pvalues.py \
    -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/interaction_mapper/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/Alzheimerdisease_InteractionResults.txt.gz \
    -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/interaction_mapper/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected-AMPAD/Alzheimerdisease_InteractionResults.txt.gz \
    -n1 all_samples \
    -n2 AMPAD_samples \
    -t CortexEUR_and_AFR_noENA_trans_100PCs_ADStatusInteraction \
    -o 2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected-all_vs_AMPAD
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data1_path = getattr(arguments, 'data1')
        self.name1 = getattr(arguments, 'name1')
        self.data2_path = getattr(arguments, 'data2')
        self.name2 = getattr(arguments, 'name2')
        self.title = getattr(arguments, 'title')
        self.out_filename = getattr(arguments, 'outfile')

        # Set variables.
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
            "no signif": "#808080",
            "x signif": "#0072B2",
            "y signif": "#D55E00",
            "both signif": "#009E73"
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
                            help="show program's version number and exit.")
        parser.add_argument("-d1",
                            "--data1",
                            type=str,
                            required=True,
                            help="The path to the first cell fraction matrix")
        parser.add_argument("-n1",
                            "--name1",
                            type=str,
                            required=False,
                            default="x",
                            help="The name for the first profile matrix")
        parser.add_argument("-d2",
                            "--data2",
                            type=str,
                            required=True,
                            help="The path to the second cell fraction matrix")
        parser.add_argument("-n2",
                            "--name2",
                            type=str,
                            required=False,
                            default="y",
                            help="The name for the second profile matrix")
        parser.add_argument("-t",
                            "--title",
                            type=str,
                            default="",
                            help="The title of the plot. Default: .")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the outfile.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        df1 = self.load_file(self.data1_path, header=0, index_col=None)
        df2 = self.load_file(self.data2_path, header=0, index_col=None)
        print(df1)
        print(df2)

        # Merge.
        df1.index = df1["SNPName"] + "_" + df1["ProbeName"]
        df2.index = df2["SNPName"] + "_" + df2["ProbeName"]
        df = df1.loc[:, ["p-value"]].merge(df2.loc[:, ["p-value"]], left_index=True, right_index=True)
        df.columns = ["x", "y"]

        # # Adding color.
        df["hue"] = "no signif"
        df.loc[(df["x"] <= 0.05) & (df["y"] > 0.05), "hue"] = "x signif"
        df.loc[(df["x"] > 0.05) & (df["y"] <= 0.05), "hue"] = "y signif"
        df.loc[(df["x"] <= 0.05) & (df["y"] <= 0.05), "hue"] = "both signif"

        # Log10 transform.
        print(df)
        df["x"] = np.log10(df["x"]) * -1
        df["y"] = np.log10(df["y"]) * -1
        print(df)

        self.plot(df=df,
                  hue="hue",
                  palette=self.palette,
                  xlabel="{} -log10(p-value)".format(self.name1.replace("_", " ")),
                  ylabel="{} -log10(p-value)".format(self.name2.replace("_", " ")),
                  title=self.title.replace("_", " "),
                  filename=self.out_filename,
                  outdir=self.outdir)

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
    def plot(df, x="x", y="y", hue=None, palette=None, xlabel="", ylabel="",
             title="", filename="plot", outdir=None):
        print(df)

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        coef, _ = stats.spearmanr(df[x], df[y])

        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        palette=palette,
                        linewidth=0,
                        legend=None,
                        ax=ax)

        ax.axhline(1.3010299956639813, ls='--', color="#000000", alpha=0.15, zorder=-1)
        ax.axvline(1.3010299956639813, ls='--', color="#000000", alpha=0.15, zorder=-1)
        ax.axline((0, 0), slope=1, ls='--', color="#000000", alpha=0.15, zorder=-1)

        # Add the text.
        ax.annotate(
            'r = {:.2f}'.format(coef),
            xy=(0.75, 0.94),
            xycoords=ax.transAxes,
            color="#404040",
            fontsize=14,
            fontweight='bold')
        ax.annotate(
            'total N = {:,}'.format(df.shape[0]),
            xy=(0.75, 0.9),
            xycoords=ax.transAxes,
            color="#404040",
            fontsize=14,
            fontweight='bold')
        if hue is not None:
            for i, group in enumerate(df[hue].unique()):
                ax.annotate(
                    'N = {:,}'.format(df.loc[df[hue] == group, :].shape[0]),
                    xy=(0.75, 0.86 - (i * 0.04)),
                    xycoords=ax.transAxes,
                    color=palette[group],
                    fontsize=14,
                    fontweight='bold')

        ax.set_title("",
                     fontsize=20,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        # Add the main title.
        fig.suptitle(title,
                     fontsize=20,
                     weight='bold')

        outpath = "{}.png".format(filename)
        if outdir is not None:
            outpath = os.path.join(outdir, outpath)
        fig.savefig(outpath)
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Dataset 1:")
        print("    > Path: {}".format(self.data1_path))
        print("    > Name: {}".format(self.name1))
        print("  > Dataset 2:")
        print("    > Path: {}".format(self.data2_path))
        print("    > Name: {}".format(self.name2))
        print("  > Title {}".format(self.title))
        print("  > Outpath {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
