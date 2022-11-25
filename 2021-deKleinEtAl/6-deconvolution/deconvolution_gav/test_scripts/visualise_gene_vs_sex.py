#!/usr/bin/env python3

"""
File:         visualise_gene_vs_sex.py
Created:      2020/10/30
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
import gzip
import argparse
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Visualise Gene VS Sex"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.expr_path = getattr(arguments, 'expression')
        self.pheno_path = getattr(arguments, 'phenotype')
        self.interest = getattr(arguments, 'interest')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'gene_vs_sex')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

        # Create color map.
        self.color_map = {"Male": "#56B4E9", "Female": "#CC79A7"}

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
        parser.add_argument("-expr",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
        parser.add_argument("-pheno",
                            "--phenotype",
                            type=str,
                            required=True,
                            help="The path to the phenotype matrix")
        parser.add_argument("-i",
                            "--interest",
                            nargs="+",
                            type=str,
                            required=True,
                            help="The expression indices to plot")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step 1 ###")
        print("Loading phenotype matrix.")
        sex_df = self.load_sex_df()
        print(sex_df)
        
        print("### Step 2 ###")
        print("Loading expression data.")
        _, expr_df = self.search_file(self.expr_path,
                                      interest=self.interest)
        print(expr_df)
        
        print("### Step 4 ###")
        print("Combine data frames.")
        df = expr_df.T.merge(sex_df,
                             left_index=True,
                             right_index=True)
        print("\tSHape: {}".format(df.shape))
        print(df)

        print("### Step 5 ###")
        print("Plot.")
        for interest in self.interest:
            subset = df[["SEX", interest]].copy()
            subset["X"] = 1
            subset.columns = ["hue", "value", "variable"]

            for extension in self.extensions:
                self.plot(subset,
                          title="{} vs SEX".format(interest),
                          ylabel="{} expression".format(interest),
                          extension=extension)

    def load_sex_df(self):
        pheno_df = self.load_file(self.pheno_path,
                                  header=0,
                                  index_col=4, 
                                  low_memory=False)

        # Combine the two gender columns, keep 'sex.by.expression' as main
        # gender ans use 'Gender' when no information is available.
        pheno_df = pheno_df.loc[:, ["Gender", "sex.by.expression"]]
        pheno_df.replace("no expression available", np.nan, inplace=True)
        pheno_df["SEX"] = pheno_df['sex.by.expression'].combine_first(
            pheno_df['Gender'])
        sex_df = pheno_df["SEX"].to_frame()
        del pheno_df
        sex_df = sex_df.replace({"SEX": {"M": "Male",
                                         "F": "Female",
                                         np.nan: "NA"}})
        
        return sex_df

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, low_memory=True):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df
    
    @staticmethod
    def search_file(path, sep="\t", header_index=0, index_col=0, nrows=None,
                    interest=None, print_interval=500):
        columns = []
        indices_int = []
        indices_str = []
        data = []

        if interest is None:
            return None
        else:
            search_list = set(interest)

        print("Searching file.")
        with gzip.open(path, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % print_interval == 0):
                    print("\tprocessed {} lines".format(i))
                if len(search_list) == 0:
                    break

                splitted_line = line.decode().strip('\n').split(sep)
                index = None
                if index_col is not None:
                    index = splitted_line[index_col]

                data_start_index = 0
                if index_col is not None:
                    data_start_index = index_col + 1

                content = splitted_line[data_start_index:]

                if header_index is not None and i == header_index:
                    columns = content
                else:
                    if interest is None or index in interest:
                        indices_int.append(i)
                        indices_str.append(index)
                        data.append(np.array(content, dtype=float))

                        search_list.remove(index)

                if nrows is not None and i > nrows:
                    break

        f.close()

        df = pd.DataFrame(data, columns=columns, index=indices_str)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return indices_int, df

    def plot(self, df, title="", xlabel="", ylabel="", extension="png"):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.violinplot(x="variable", y="value", hue="hue",
                       data=df, ax=ax, split=True, scale="count",
                       inner="quartile", palette=self.color_map)
        fig.suptitle(title)
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_vs_{}_violin.{}".format(title, ylabel, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Phenotype path: {}".format(self.pheno_path))
        print("  > Interest path: {}".format(self.interest))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
