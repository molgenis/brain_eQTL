#!/usr/bin/env python3

"""
File:         dim_reduction_plot.py
Created:      2020/12/03
Last Changed: 2020/12/04
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
__program__ = "Dimension Reduction Plot"
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
/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/dim_reduction_plot.py -ex /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.pkl -m /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt -mid rnaseq_id -a PCA TSNE

./dim_reduction_plot.py -ex /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.pkl -m /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt -mid rnaseq_id -a PCA TSNE
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.expr_path = getattr(arguments, 'expression')
        self.meta_path = getattr(arguments, 'meta_table')
        self.meta_id = getattr(arguments, 'meta_id')
        self.color_id = getattr(arguments, 'color_id')
        self.algorithm = getattr(arguments, 'algorithm')
        self.extensions = getattr(arguments, 'extension')

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
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
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
                            default=["MetaCohort", "predicted.brain.region"],
                            help="The ID column in the -m / --meta_table matrix"
                                 "on which to color. Default: [cohort,"
                                 " predicted.brain.region].")
        parser.add_argument("-a",
                            "--algorithm",
                            nargs="*",
                            type=str,
                            required=False,
                            choices=["PCA", "TSNE", "UMAP"],
                            default=["PCA", "TSNE", "UMAP"],
                            help="The type of dimensionality reduction plots to"
                                 "create. Default: PCA, TSNE, and UMAP.")
        parser.add_argument("-e",
                            "--extension",
                            type=str,
                            choices=["png", "pdf"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Load the expression file")
        expr_df = self.load_file(self.expr_path)
        print(expr_df)

        meta_df = self.load_file(self.meta_path, index_col=None, low_memory=False)
        if self.meta_id is None:
            self.meta_id = meta_df.columns[0]
        meta_df.set_index(self.meta_id, inplace=True)
        print(meta_df)
        print("", flush=True)

        dim_reduction_df = None
        if "PCA" in self.algorithm:
            print("Performing PCA")
            pca_projection, pca_expl_var = self.get_pca_components(expr_df.T, n=2)
            pca_plot_df = pca_projection.merge(meta_df, left_index=True, right_index=True)
            print(pca_plot_df[["PC1", "PC2"] + self.color_id])

            outdir = os.path.join(self.outdir, "PCA")
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            print("Plotting PCA", flush=True)
            for hue in self.color_id:
                self.plot(df=pca_plot_df,
                          x="PC1",
                          y="PC2",
                          hue=hue,
                          xlabel="PC1 [{:.2f}%]".format(pca_expl_var["PC1"]),
                          ylabel="PC2 [{:.2f}%]".format(pca_expl_var["PC2"]),
                          title="PCA - {}".format(hue),
                          outdir=outdir)
            print("", flush=True)

            if dim_reduction_df is None:
                dim_reduction_df = pca_projection
            else:
                dim_reduction_df = dim_reduction_df.merge(pca_projection, left_index=True, right_index=True)

        if "TSNE" in self.algorithm:
            print("Performing TSNE")
            tsne_embedding = self.get_tsne_components(expr_df.T, n=2)
            tsne_plot_df = tsne_embedding.merge(meta_df, left_index=True, right_index=True)
            print(tsne_plot_df[["TSNE1", "TSNE2"] + self.color_id])

            outdir = os.path.join(self.outdir, "TSNE")
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            print("Plotting TSNE", flush=True)
            for hue in self.color_id:
                self.plot(df=tsne_plot_df,
                          x="TSNE1",
                          y="TSNE2",
                          hue=hue,
                          xlabel="TSNE1",
                          ylabel="TSNE2",
                          title="TSNE - {}".format(hue),
                          outdir=outdir)
            print("")

            if dim_reduction_df is None:
                dim_reduction_df = tsne_embedding
            else:
                dim_reduction_df = dim_reduction_df.merge(tsne_embedding, left_index=True, right_index=True)

        if "UMAP" in self.algorithm:
            print("Performing UMAP")
            umap_embedding = self.get_umap_components(expr_df.T)
            umap_plot_df = umap_embedding.merge(meta_df, left_index=True,
                                                right_index=True)
            print(umap_plot_df[["UMAP1", "UMAP2"] + self.color_id])

            outdir = os.path.join(self.outdir, "UMAP")
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            print("Plotting UMAP", flush=True)
            for hue in self.color_id:
                self.plot(df=umap_plot_df,
                          x="UMAP1",
                          y="UMAP2",
                          hue=hue,
                          xlabel="UMAP1",
                          ylabel="UMAP2",
                          title="UMAP - {}".format(hue),
                          outdir=outdir)
            print("")

            if dim_reduction_df is None:
                dim_reduction_df = umap_embedding
            else:
                dim_reduction_df = dim_reduction_df.merge(umap_embedding, left_index=True, right_index=True)

        print("Saving dim. reduction data frame")
        dim_reduction_path = os.path.join(self.outdir, "dim_reduction_df.txt.gz")
        if os.path.exists(dim_reduction_path):
            print("\tOutput file 'dim_reduction_df.txt.gz' found.")
            prev_dim_reduction_df = self.load_file(dim_reduction_path)

            filter = []
            for col in prev_dim_reduction_df:
                if col not in dim_reduction_df:
                    filter.append(col)

            print("\tOutput file contains {} columns not present in current "
                  "data frame: '{}'.".format(len(filter), ", ".join(filter)))

            dim_reduction_df = dim_reduction_df.merge(prev_dim_reduction_df.loc[:, filter], left_index=True, right_index=True)

        dim_reduction_df.index.name = "-"
        dim_reduction_df.to_csv(dim_reduction_path, compression="gzip", sep="\t", header=True, index=True)
        print("\tSaved dataframe: dim_reduction_df.txt.gz "
              "with shape: {}".format(dim_reduction_df.shape))

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

    @staticmethod
    def get_pca_components(X, n=2):
        start_time = time.time()

        pca = PCA(n_components=n)
        pca.fit(X)
        projections = pd.DataFrame(pca.transform(X), index=X.index)
        projections.columns = ["PC{}".format(i+1) for i, _ in enumerate(projections.columns)]
        explained_variance_ratio = {}
        for i, col in enumerate(projections.columns):
            explained_variance_ratio[col] = pca.explained_variance_ratio_[i] * 100

        run_time_min, run_time_sec = divmod(time.time() - start_time, 60)
        print("PCA finished in {} minute(s) and "
              "{} second(s).".format(int(run_time_min),
                                     int(run_time_sec)))
        return projections, explained_variance_ratio

    @staticmethod
    def get_tsne_components(X, n=2):
        start_time = time.time()

        tsne = TSNE(n_components=n, verbose=1)
        embedding = pd.DataFrame(tsne.fit_transform(X), index=X.index)
        embedding.columns = ["TSNE{}".format(i+1) for i, _ in enumerate(embedding.columns)]

        run_time_min, run_time_sec = divmod(time.time() - start_time, 60)
        print("TSNE finished in {} minute(s) and "
              "{} second(s).".format(int(run_time_min),
                                     int(run_time_sec)))
        return embedding

    @staticmethod
    def get_umap_components(X, n=2):
        start_time = time.time()

        reducer = umap.UMAP(n_components=n, verbose=True)
        embedding = pd.DataFrame(reducer.fit_transform(X), index=X.index)
        embedding.columns = ["UMAP{}".format(i+1) for i, _ in enumerate(embedding.columns)]

        run_time_min, run_time_sec = divmod(time.time() - start_time, 60)
        print("UMAP finished in {} minute(s) and "
              "{} second(s).".format(int(run_time_min),
                                     int(run_time_sec)))
        return embedding

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
            fig.savefig(os.path.join(self.outdir, filename))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Meta table path: {}".format(self.meta_path))
        print("  > Meta table ID: {}".format(self.meta_id))
        print("  > Meta table color ID: {}".format(self.color_id))
        print("  > Algorithm: {}".format(self.algorithm))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("", flush=True)


if __name__ == '__main__':
    m = main()
    m.start()
