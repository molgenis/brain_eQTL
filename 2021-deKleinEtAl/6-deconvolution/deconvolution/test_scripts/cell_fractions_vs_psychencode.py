#!/usr/bin/env python3

"""
File:         cell_fractions_vs_psychencode.py
Created:      2021/10/01
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
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Cell Fractions VS PsychENCODE"
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


"""
Syntax:
./cell_fractions_vs_psychencode.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/OLD/ContainsDuplicateSamples/CortexEUR-cis/perform_deconvolution/deconvolution_table_CNS7.txt.gz
"""

class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.input_cf_path = getattr(arguments, 'cell_fractions')
        self.phenotype_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-03-09.brain.phenotypes.txt"
        self.psychencode_cf_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/DER-24_Cell_fractions_Normalized.xlsx"

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.cell_type_dict = {
            "Astrocytes" : "Astrocyte",
            "Endothelial" : "EndothelialCell",
            "Ex" : "Excitatory",
            "In" : "Inhibitory",
            "Microglia" : "Microglia",
            "OPC" : "OPC",
            "Oligo" : "Oligodendrocyte",
            "OtherNeuron" : "OtherNeuron",
            "Quiescent" : "Quiescent",
            "Replicating" : "Replicating"
        }

        self.palette = {
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "OPC": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Pericytes": "#808080",
            "Microglia/Macrophage": "#E69F00",
            "Excitatory/Neuron": "#56B4E9",
            "Inhibitory/Neuron": "#0072B2",
            "Excitatory+Inhibitory/Neuron": "#BEBEBE",
            "-": "#000000"
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
                            help="show program's version number and exit")
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fracttions matrix")

        return parser.parse_args()

    def start(self):

        print("Step 1: loading data")
        input_cf_df = self.load_file(path=self.input_cf_path)
        psychencode_cf_df = self.load_file(path=self.psychencode_cf_path, sheet_name="Sheet1")
        phenotype_df = self.load_file(path=self.phenotype_path, header=0, index_col=None, low_memory=False)

        print("Step 2: pre-processing")
        # Summarize.
        psychencode_cf_df["cell type"] = [''.join([i for i in x if not i.isdigit()]) for x in [x.split("-")[1] for x in psychencode_cf_df.index]]
        psychencode_cf_df = psychencode_cf_df.groupby(["cell type"]).sum()
        psychencode_cf_df = psychencode_cf_df.T

        # Remove columns with no cell fraction.
        psychencode_cf_df = psychencode_cf_df.loc[:, psychencode_cf_df.sum(axis=0) != 0]

        # Align sample ID's.
        rownames = []
        for sample_id in psychencode_cf_df.index:
            rnaseq_id = sample_id
            rnaseq_id_s = phenotype_df.loc[phenotype_df['SampleFull'].str.contains(sample_id), 'rnaseq_id']
            if rnaseq_id_s.shape[0] == 1:
                rnaseq_id = rnaseq_id_s.values[0]
            rownames.append(rnaseq_id)
        psychencode_cf_df.index = rownames

        # Align cell types.
        input_cf_df.columns = [cell_type.split("_")[1] for cell_type in input_cf_df.columns]
        psychencode_cf_df.columns = [self.cell_type_dict[cell_type] for cell_type in psychencode_cf_df.columns]

        print(input_cf_df)
        print(psychencode_cf_df)

        print("Step 3: merging data")
        input_cf_df.reset_index(drop=False, inplace=True)
        psychencode_cf_df.reset_index(drop=False, inplace=True)
        input_cf_df_m = input_cf_df.melt(id_vars=["index"])
        psychencode_cf_df_m = psychencode_cf_df.melt(id_vars=["index"])

        comparison_df = input_cf_df_m.merge(psychencode_cf_df_m, on=["index", "variable"])
        comparison_df.columns = ["index", "cell type", "MetaBrain", "PsychENCODE"]
        print(comparison_df)

        print("Step 4: plotting")
        self.plot_regplot(df=comparison_df,
                          x="MetaBrain",
                          y="PsychENCODE",
                          hue="cell type",
                          palette=self.palette,
                          title="Cell fraction predictions",
                          name="MetaBrain_vs_PsychENCODE")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  skiprows=0, sheet_name=None, low_memory=True):
        if path.endswith(".xlsx"):
            df = pd.read_excel(path, header=header, index_col=index_col,
                         nrows=nrows, skiprows=skiprows, sheet_name=sheet_name)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, skiprows=skiprows,
                             low_memory=low_memory)

        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def plot_regplot(self, df, x="x", y="y", hue=None, palette=None,
                     xlabel=None, ylabel=None, title="", name="plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        total_coef, _ = stats.pearsonr(df[x], df[y])

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        if hue is None:
            df[hue] = "-"

        groups = df[hue].unique()

        coefs = {}
        for group in groups:
            subset = df.loc[df[hue] == group, :].copy()
            color = self.palette[group]

            coef, p = stats.pearsonr(subset[x], subset[y])
            coefs[group] = "r = {:.2f}".format(coef)

            sns.regplot(x=x, y=y, data=subset,
                        scatter_kws={'facecolors': color,
                                     'edgecolors': "#808080"},
                        line_kws={"color": color},
                        ax=ax
                        )

        ax.annotate(
            'r = {:.2f}'.format(total_coef),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=14,
            fontweight='bold')
        ax.annotate(
            'N = {:,.0f}'.format(df.shape[0] / len(groups)),
            xy=(0.03, 0.9),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=14,
            fontweight='bold')

        max_value = max(df[x].max(), df[y].max()) + 0.05

        ax.set_xlim(-0.01, max_value)
        ax.set_ylim(-0.01, max_value)
        ax.plot([-0.01, max_value], [-0.01, max_value], ls="--", c=".3")

        ax.set_title(title,
                     fontsize=16,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        if hue is not None and palette is not None:
            handles = []
            for key, color in self.palette.items():
                if key in groups and key in coefs.keys():
                    handles.append(mpatches.Patch(color=color, label="{} [{}]".format(key, coefs[key])))
            ax.legend(handles=handles, loc=4)

        fig.savefig(os.path.join(self.outdir, "{}_regression.png".format(name)))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
