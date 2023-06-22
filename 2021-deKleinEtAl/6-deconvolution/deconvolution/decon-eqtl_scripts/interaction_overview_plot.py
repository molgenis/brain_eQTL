#!/usr/bin/env python3

"""
File:         interaction_overview_plot.py
Created:      2021/10/25
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
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Interaction Overview Plot"
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
./interaction_overview_plot.py -d ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -calc_fdr
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.decon_path = getattr(arguments, 'decon')
        self.calc_fdr = getattr(arguments, 'calc_fdr')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'interaction_overview_plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00"
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
        parser.add_argument("-d",
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix")
        parser.add_argument("-calc_fdr",
                            action='store_true',
                            help="Treat input as p-values and convert them to "
                                 "FDR using BH multiple testing. "
                                 " default: False.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()
        print("Loading data.")
        decon_df = self.load_file(self.decon_path, nrows=None)
        decon_df.columns = [x.replace("CellMapNNLS_", "") for x in decon_df.columns]

        print("Calculating FDR.")
        decon_pval_df, decon_fdr_df = self.bh_correct(decon_df)

        variable = "p-value"
        df = decon_pval_df
        if self.calc_fdr:
            df = decon_fdr_df
            variable = "FDR"
        self.print_n_signif(df=df, variable=variable)

        print("Preprocessing data.")
        # Convert to booleans.
        df = (df <= 0.05).astype(int)

        # Split data.
        no_interaction_df = df.loc[df.sum(axis=1) == 0, :]
        single_interaction_df = df.loc[df.sum(axis=1) == 1, :]
        multiple_interactions_df = df.loc[df.sum(axis=1) > 1, :]

        n_ieqtls_unique_per_pic = [(name, value) for name, value in single_interaction_df.sum(axis=0).iteritems()]
        n_ieqtls_unique_per_pic.sort(key=lambda x: -x[1])
        cell_types = [x[0] for x in n_ieqtls_unique_per_pic]
        n_ieqtls = [x[1] for x in n_ieqtls_unique_per_pic]

        data = [no_interaction_df.shape[0]] + n_ieqtls + [multiple_interactions_df.shape[0]]
        labels = ["None"] + cell_types + ["Multiple"]
        colors = None
        if self.palette is not None:
            colors = ["#D3D3D3"] + [self.palette[ct] for ct in cell_types] + ["#808080"]
        explode = [0.] + [0.1 for _ in range(df.shape[1])] + [0.1]

        data_sum = np.sum(data)
        for i in range(len(data)):
            print("{} (N = {}) = {:.2f}%".format(labels[i], data[i], (100 / data_sum) * data[i]))

        total_n_ieqtls = df.shape[0] - no_interaction_df.shape[0]
        print("Total interactions (N = {}) = {:.2f}%".format(total_n_ieqtls, (100 / df.shape[0]) * total_n_ieqtls))

        print("Creating plot.")
        self.plot(data=data, labels=labels, explode=explode, colors=colors, extension="png")
        self.plot(data=data, labels=labels, explode=explode, colors=colors, extension="pdf", label_threshold=100)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def bh_correct(pvalue_df):
        df = pvalue_df.copy()
        pval_data = []
        fdr_data = []
        indices = []
        for col in df.columns:
            if col.endswith("_pvalue"):
                pval_data.append(df.loc[:, col])
                fdr_data.append(multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
                indices.append(col.replace("_pvalue", ""))
        pval_df = pd.DataFrame(pval_data, index=indices, columns=df.index)
        fdr_df = pd.DataFrame(fdr_data, index=indices, columns=df.index)

        return pval_df.T, fdr_df.T

    def print_n_signif(self, df, variable):
        m = df.to_numpy()
        colnames = df.columns.tolist()

        print("\nN-interaction ({} <= 0.05):".format(variable))
        n_hits_a = (m <= 0.05).sum(axis=0)
        n_hits_total = np.sum(n_hits_a)
        cov_length = np.max([len(x) for x in colnames])
        hits_length = np.max([len(str(x)) for x in n_hits_a] + [len(str(n_hits_total))])
        for n_hits, cell_type in zip(n_hits_a, colnames):
            print("\t{:{}s}  {:{}d}".format(cell_type, cov_length, n_hits, hits_length))
        print("\t{}".format("".join(["-"] * cov_length)))
        print("\t{:{}s}  {:{}d}".format("total", cov_length, n_hits_total, hits_length))

        print("", flush=True)

    def plot(self, data, labels, explode=None, colors=None, extension='png', label_threshold=0):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        plt.pie(data,
                autopct=lambda pct: self.autopct_func(pct, data, label_threshold),
                explode=explode,
                labels=labels,
                shadow=False,
                colors=colors,
                startangle=90,
                wedgeprops={'linewidth': 1})
        ax.axis('equal')

        fig.savefig(os.path.join(self.outdir, "interaction_piechart.{}".format(extension)))
        plt.close()

    @staticmethod
    def autopct_func(pct, allvalues, label_threshold):
        if pct >= label_threshold:
            absolute = int(pct / 100. * np.sum(allvalues))
            return "{:.1f}%\n(N = {:,.0f})".format(pct, absolute)
        else:
            return ""

    def print_arguments(self):
        print("Arguments:")
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Calc FDR: {}".format(self.calc_fdr))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
