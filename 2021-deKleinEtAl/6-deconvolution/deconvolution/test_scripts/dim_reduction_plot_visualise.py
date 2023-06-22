#!/usr/bin/env python3

"""
File:         dim_reduction_plot_visualise.py
Created:      2020/12/04
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
import time
import os

# Third party imports.
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Dimension Reduction Plot Visualise"
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
./dim_reduction_plot_visualise.py -dr /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/dim_reduction_plot/dim_reduction_df.txt.gz -m /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt -mid rnaseq_id

./dim_reduction_plot_visualise.py -dr /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/dim_reduction_plot/dim_reduction_df.txt.gz -m /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/partial_deconvolution/NEW_PROFILE_NOPERICYTES_CORTEX_EUR_TMM_LOG2/IHC_0CPM_LOG2_FILTERED_CC/deconvolution.txt.gz
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.dim_reduction_path = getattr(arguments, 'dim_reduction')
        self.meta_path = getattr(arguments, 'meta_table')
        self.meta_id = getattr(arguments, 'meta_id')
        self.color_id = getattr(arguments, 'color_id')
        self.extensions = getattr(arguments, 'extensions')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent), 'dim_reduction_plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set palettes.
        self.cohort_palette = {
            "AMP-AD": "#A9AAA9",
            "Braineac": "#E7A023",
            "Brainseq": "#CC79A9",
            "CMC": "#EEE643",
            "GTEx": "#1A9F74",
            "NABEC": "#767833",
            "TargetALS": "#DDCC78",
            "ENA": "#D56228",
            "PsychEncode": "#5AB4E5"
        }

        self.tissue_palette = {
            "cortex": "#0072B2",
            "basalganglia": "#009E73",
            "hippocampus": "#F0E442",
            "amygdala": "#CC79A7",
            "spinalcord": "#56B4E9",
            "cerebellum": "#D55E00",
            "hypothalamus": "#E69F00",
            "no_predicted_region_available": "#000000"
        }

        self.cell_type_palette = {
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
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
                            help="show program's version number and exit")
        parser.add_argument("-dr",
                            "--dim_reduction",
                            type=str,
                            required=True,
                            help="The path to the dim. reduction matrix.")
        parser.add_argument("-m",
                            "--meta_table",
                            type=str,
                            required=True,
                            help="The path to the meta data matrix.")
        parser.add_argument("-mid",
                            "--meta_id",
                            type=str,
                            required=False,
                            default=None,
                            help="The ID column in the -m / --meta_table matrix"
                                 "on which to merge. Default: index.")
        parser.add_argument("-c",
                            "--color_id",
                            nargs="*",
                            type=str,
                            required=False,
                            default=None,
                            help="The ID column in the -m / --meta_table matrix"
                                 "on which to color. Default: None.")
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

        print("Load the dim. reduction file")
        dim_red_df = self.load_file(self.dim_reduction_path)
        print(dim_red_df)

        meta_df = self.load_file(self.meta_path, index_col=None, low_memory=False)
        if self.meta_id is None:
            self.meta_id = meta_df.columns[0]
        meta_df.set_index(self.meta_id, inplace=True)
        meta_df.fillna("NA", inplace=True)
        print(meta_df)
        print("", flush=True)

        if self.color_id is None:
            self.color_id = list(meta_df.columns)

        print("Merging data")
        df = meta_df.merge(dim_red_df, left_index=True, right_index=True)
        del dim_red_df, meta_df

        for algorithm, x, y in [("PCA", "PC1", "PC2"), ("TSNE", "TSNE1", "TSNE2"), ("UMAP", "UMAP1", "UMAP2")]:
            if x in df.columns and y in df.columns:
                print("Plotting {}".format(algorithm), flush=True)
                for hue in self.color_id:
                    if len(df.loc[:, hue].unique()) < 0 or "/" in hue:
                        continue

                    outdir = os.path.join(self.outdir, algorithm)
                    if not os.path.exists(outdir):
                        os.makedirs(outdir)

                    self.plot(df=df,
                              x=x,
                              y=y,
                              hue=hue,
                              xlabel=x,
                              ylabel=y,
                              title="{} - {}".format(algorithm, hue),
                              outdir=outdir)
                print("", flush=True)

        print("Program finished", flush=True)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None, low_memory=True):
        if path.endswith(".pkl"):
            df = pd.read_pickle(path)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def plot(self, df, x="x", y="y", hue=None, xlabel="", ylabel="",
             title="", outdir=None):
        if outdir is None:
            outdir = self.outdir

        palette = None
        legend = True
        if hue == "MetaCohort":
            palette = self.cohort_palette
            legend = False
        elif hue == "predicted.brain.region":
            palette = self.tissue_palette
            legend = False

        if str(dict(df.dtypes)[hue]).startswith('float') or str(dict(df.dtypes)[hue]).startswith('int'):
            palette = sns.color_palette("light:{}".format(self.cell_type_palette[hue]), as_cmap=True)

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, gridspec_kw={"width_ratios": [0.9, 0.1]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        legend=legend,
                        palette=palette,
                        ax=ax1)

        ax1.set_title(title,
                      fontsize=20,
                      fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_xlabel(xlabel,
                       fontsize=14,
                       fontweight='bold')

        if legend:
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        else:
            handles = []
            for label, color in palette.items():
                handles.append(mpatches.Patch(color=color, label=label))
            ax2.legend(handles=handles, loc="center")

        plt.show()

        # Safe the plot.
        for extension in self.extensions:
            filename = "dimReductionPlot_{}_vs_{}_coloredBy{}.{}".format(x, y, hue, extension)
            print("\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(outdir, filename))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Dim Reduction path: {}".format(self.dim_reduction_path))
        print("  > Meta table path: {}".format(self.meta_path))
        print("  > Meta table ID: {}".format(self.meta_id))
        print("  > Meta table color ID: {}".format(self.color_id))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("", flush=True)


if __name__ == '__main__':
    m = main()
    m.start()
