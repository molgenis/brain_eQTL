#!/usr/bin/env python3

"""
File:         plot_cell_fraction_vs_IHC.py
Created:      2022/01/25
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
import json
import os

# Third party imports.
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Plot Cell Fraction vs IHC"
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

./plot_cell_fraction_vs_IHC.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz -o 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected -e png pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.outfile = getattr(arguments, 'outfile')
        self.extensions = getattr(arguments, 'extensions')

        # Set other input arguments.
        self.ihc_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz"

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "plot_cell_fraction_vs_IHC")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

        self.palette = {
            "Macrophage/Microglia": "#E69F00",
            "EndothelialCell": "#CC79A7",
            "Oligodendrocyte": "#009E73",
            "Neuron": "#0072B2",
            "Astrocyte": "#D55E00"
        }

        self.ct_match = {
            "Macrophage/Microglia": (["Microglia"], ["Macrophage"]),
            "EndothelialCell": (["EndothelialCell"], ["EndothelialCell"]),
            "Oligodendrocyte": (["Oligodendrocyte"], ["Oligodendrocyte"]),
            "Neuron": (["Excitatory", "Inhibitory", "OtherNeuron"], ["Neuron"]),
            "Astrocyte": (["Astrocyte"], ["Astrocyte"]),
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
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fracttions matrix")
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
        ihc_df = self.load_file(self.ihc_path, header=0, index_col=0)

        print(cf_df)
        print(ihc_df)

        print("Merging.")
        df_list = []
        for ct_label, (cf_columns, ihc_columns) in self.ct_match.items():
            df1 = cf_df[cf_columns].sum(axis=1).to_frame()
            df2 = ihc_df[ihc_columns].sum(axis=1).to_frame()
            df = df1.merge(df2, left_index=True, right_index=True)
            df.columns = ["x", "y"]
            df["hue"] = ct_label
            df_list.append(df)

        df = pd.concat(df_list, axis=0)
        df["x"] = df["x"] * 100
        df["y"] = df["y"] * 100
        print(df)

        print("Plotting.")
        self.single_regplot(df=df,
                            hue="hue",
                            palette=self.palette,
                            xlabel="Predicted cell type proportion",
                            ylabel="IHC cell type proportion",
                            filename=self.outfile + "_vs_IHC_regplot")

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def single_regplot(self, df, x="x", y="y", hue=None, palette=None,
                       xlabel=None, ylabel=None, title="", filename="plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
                                       gridspec_kw={"width_ratios": [0.8, 0.2]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        n = df.shape[0]
        if hue is not None:
            n = df.shape[0] / len(df[hue].unique())

        # Set annotation.
        pearson_coef, _ = stats.pearsonr(df[y], df[x])
        ax1.annotate(
            'N = {:,.0f}'.format(n),
            xy=(0.03, 0.94),
            xycoords=ax1.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold')
        ax1.annotate(
            'Pearson r = {:.2f}'.format(pearson_coef),
            xy=(0.03, 0.90),
            xycoords=ax1.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold')

        group_column = hue
        if hue is None:
            df["hue"] = "#000000"
            group_column = "hue"

        group_corr_coef = {}
        for i, hue_group in enumerate(df[group_column].unique()):
            subset = df.loc[df[group_column] == hue_group, :]
            if subset.shape[0] < 2:
                continue

            facecolors = "#000000"
            color = "#b22222"
            if palette is not None:
                facecolors = palette[hue_group]
                color = facecolors

            sns.regplot(x=x, y=y, data=subset,
                        scatter_kws={'facecolors': facecolors,
                                     's': 75,
                                     'linewidth': 0},
                        line_kws={"color": color},
                        ax=ax1)

            if hue is not None:
                subset_pearson_coef, _ = stats.pearsonr(subset[y], subset[x])
                group_corr_coef[hue_group] = subset_pearson_coef

        if hue is not None:
            handles = []
            for hue_group in df[group_column].unique():
                if hue_group in palette:
                    r = "NA"
                    if hue_group in group_corr_coef:
                        r = "{:.2f}".format(group_corr_coef[hue_group])
                    handles.append([mpatches.Patch(color=palette[hue_group], label="{} [r={}]".format(hue_group, r)), group_corr_coef[hue_group]])
            handles.sort(key=lambda x: -x[1])
            handles = [x[0] for x in handles]
            ax2.legend(handles=handles, loc="center", fontsize=8)

        ax1.set_xlabel(xlabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_title(title,
                      fontsize=18,
                      fontweight='bold')

        # Change margins.
        xlim = ax1.get_xlim()
        ylim = ax1.get_ylim()

        xmargin = (xlim[1] - xlim[0]) * 0.05
        ymargin = (ylim[1] - ylim[0]) * 0.05

        new_xlim = (xlim[0] - xmargin, xlim[1] + xmargin)
        new_ylim = (ylim[0] - ymargin, ylim[1] + ymargin)

        min_point = min(new_xlim[0], new_ylim[0])
        max_point = min(new_xlim[1], new_ylim[1])

        ax1.set_xlim(min_point, max_point)
        ax1.set_ylim(min_point, max_point)

        ax1.plot([min_point, max_point], [min_point, max_point], ls="--", c=".3")

        for extension in self.extensions:
            outpath = os.path.join(self.outdir, "{}.{}".format(filename, extension))
            fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(filename)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Cell fraction file: {}".format(self.cf_path))
        print("  > Outfile: {}".format(self.outfile))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
