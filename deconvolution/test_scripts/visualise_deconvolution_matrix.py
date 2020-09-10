#!/usr/bin/env python3

"""
File:         visualise_deconvolution_matrix.py
Created:      2020/09/04
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
import itertools
import argparse
import os

# Third party imports.
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

# Local application imports.

# Metadata
__program__ = "Visualise Deconvolution Matrix"
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
        self.info_path = getattr(arguments, 'info')
        self.sample_id = getattr(arguments, 'sample')
        self.x_id = getattr(arguments, 'x')
        self.row_id = getattr(arguments, 'row')
        self.extension = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = str(Path(__file__).parent.parent)
        self.colormap = {
            "": "#FFFFFF",
            "Neuron": "#b38d84",
            "Oligodendrocyte": "#5d9166",
            "EndothelialCell": "#f2a7a7",
            "Microglia": "#e8c06f",
            "Macrophage": "#e8c06f",
            "Astrocyte": "#9b7bb8"
        }

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
        parser.add_argument("-i",
                            "--info",
                            type=str,
                            required=False,
                            default="x",
                            help="The name for the sample information matrix")
        parser.add_argument("-sample",
                            type=str,
                            required=True,
                            help="The sample column.")
        parser.add_argument("-x",
                            type=str,
                            required=True,
                            help="The information column for the x-axis.")
        parser.add_argument("-row",
                            type=str,
                            required=True,
                            help="The information column for the rows.")
        parser.add_argument("-e",
                            "--extension",
                            type=str,
                            choices=["png", "pdf"],
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

        print("Loading information matrix.")
        info_df = pd.read_csv(self.info_path, sep="\t", header=0,
                              low_memory=False)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.info_path),
                                      info_df.shape))

        print("Preprocessing.")
        x_dict = dict(zip(info_df.loc[:, self.sample_id], info_df.loc[:, self.x_id]))
        row_dict = dict(zip(info_df.loc[:, self.sample_id], info_df.loc[:, self.row_id]))

        df = decon_df.reset_index().melt(id_vars=["index"])
        df[self.x_id] = df["index"].map(x_dict).astype(str)
        df[self.row_id] = df["index"].map(row_dict).astype(str)

        print("Visualizing.")

        order = list(df[self.x_id].unique())
        order.sort()

        rows = list(df[self.row_id].unique())
        rows.sort()

        cols = list(df["variable"].unique())
        cols.sort()

        sns.set(style="ticks")
        fig, axes = plt.subplots(ncols=len(cols),
                                 nrows=len(rows),
                                 figsize=(4*len(cols), 4*len(rows)),
                                 sharex="all")

        for i, row_variable in enumerate(rows):
            for j, col_variable in enumerate(cols):
                subset = df[(df[self.row_id] == row_variable) & (df["variable"] == col_variable)].copy()
                counts = subset[self.x_id].value_counts()

                ax = axes[i, j]
                print(i, j, row_variable, col_variable)
                sns.despine(fig=fig, ax=ax)
                g = sns.boxplot(x=self.x_id, y="value", data=subset,
                                order=order, ax=ax)

                medians = subset.groupby([self.x_id])['value'].median()
                sizes = subset.groupby([self.x_id])['value'].count()
                vertical_offset = subset['value'].median() * 0.01

                for xtick, label in enumerate(order):
                    if label in medians and label in sizes:
                        g.text(xtick, medians[label] + vertical_offset,
                               sizes[label],
                               horizontalalignment='center', size='xx-small',
                               color='lightgrey', weight='semibold')

                ylabel = ""
                if j == 0:
                    ylabel = row_variable

                title = ""
                if i == 0:
                    title = col_variable
                if i == (len(rows) - 1):
                    ax.set_xticks(range(len(order)))
                    ax.set_xticklabels(order,
                                       rotation=45,
                                       ha='right',
                                       fontsize=10,
                                       fontweight='bold')

                ax.set_title(title,
                             fontsize=11,
                             fontweight='bold',
                             color=self.colormap[title])
                ax.set_xlabel("",
                              fontsize=10,
                              fontweight='bold')
                ax.set_ylabel(ylabel,
                              fontsize=10,
                              fontweight='bold')

                box_pairs = itertools.combinations(counts.index, 2)
                if len(counts.index) >= 2:
                    add_stat_annotation(ax, data=subset, x=self.x_id, y="value",
                                        order=order, box_pairs=box_pairs,
                                        test='Mann-Whitney', text_format='star',
                                        loc='inside', fontsize='small',
                                        line_height=0.01, text_offset=0.1)

        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "x{}_row{}_catplot.{}".format(self.x_id.replace(" ", "_"), self.row_id.replace(" ", "_"), self.extension)))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
