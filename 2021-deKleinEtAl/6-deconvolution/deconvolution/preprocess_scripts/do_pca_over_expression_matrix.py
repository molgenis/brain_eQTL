#!/usr/bin/env python3

"""
File:         do_pca_over_expression_matrix.py
Created:      20212/02/09
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
from functools import reduce
import argparse
import time
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.regression.linear_model import OLS
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Do PCA over Expression Matrix"
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
./do_pca_over_expression_matrix.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-EUR/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovarOLSPerDataset.CovariatesRemovedOLS.80PCAsOverSamplesRemoved.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex_EUR.txt.gz \
    -n 25 \
    -of 2022-02-09-MetaBrain_CortexEUR
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.n_components = getattr(arguments, 'n_components')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        if outdir is None:
            outdir = str(Path(__file__).parent.parent)
        self.plot_outdir = os.path.join(outdir, 'do_pca_over_expression_matrix', outfolder, 'plot')
        self.file_outdir = os.path.join(outdir, 'do_pca_over_expression_matrix', outfolder, 'data')
        for outdir in [self.plot_outdir, self.file_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        self.palette = {
            "AMPAD-MAYO-V2": "#9C9FA0",
            "CMC_HBCC_set2": "#0877B4",
            "GTEx": "#0FA67D",
            "AMPAD-ROSMAP-V2": "#6950A1",
            "BrainGVEX-V2": "#48B2E5",
            "TargetALS": "#D5C77A",
            "AMPAD-MSBB-V2": "#5CC5BF",
            "NABEC-H610": "#6D743A",
            "LIBD_1M": "#808080",
            "LIBD_5M": "#808080",
            "ENA": "#D46727",
            "LIBD_h650": "#808080",
            "GVEX": "#48B2E5",
            "NABEC-H550": "#6D743A",
            "CMC_HBCC_set3": "#0877B4",
            "UCLA_ASD": "#F36D2A",
            "CMC": "#EAE453",
            "CMC_HBCC_set1": "#0877B4",
            "Braineac": "#E49D26",
            "Bipseq_1M": "#000000",
            "Bipseq_h650": "#000000",
            "Brainseq": "#C778A6"
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
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the data matrix.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-n",
                            "--n_components",
                            type=int,
                            default=100,
                            help="The number of PCA components to calculate. "
                                 "Default: 100.")
        parser.add_argument("-od",
                            "--outdir",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of the output path.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # Construct the output filename.
        filename = os.path.basename(self.data_path).replace(".gz", "").replace(".txt", "")

        # Load sample-dataset file.
        print("Loading sample-to-dataset.")
        std_df = self.load_file(self.std_path, header=0, index_col=None)

        # Pre-process data.
        print("Pre-processing samples-to-dataset.")
        sample_to_dataset = dict(zip(std_df.iloc[:, 0], std_df.iloc[:, 1]))

        # Load data.
        print("Loading data.")
        df = self.load_file(self.data_path, header=0, index_col=0)

        print("PCA analysis.")
        self.pca(df=df,
                 filename=filename,
                 sample_to_dataset=sample_to_dataset,
                 file_appendix="",
                 plot_appendix="")

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def pca(self, df, filename, sample_to_dataset, file_appendix="", plot_appendix=""):
        # samples should be on the columns and genes on the rows.
        zscores = (df - df.mean(axis=0)) / df.std(axis=0)
        pca = PCA(n_components=self.n_components)
        pca.fit(zscores)
        components_df = pd.DataFrame(pca.components_)
        components_df.index = ["PC{}".format(i + 1) for i, _ in enumerate(components_df.index)]
        components_df.columns = df.columns

        print("\tSaving file.")
        self.save_file(df=components_df, outpath=os.path.join(self.file_outdir, "{}.PCAOverSamplesEigenvectors.txt.gz".format(filename, file_appendix)))

        print("\tPlotting PCA")
        plot_df = components_df.T
        plot_df["hue"] = plot_df.index.map(sample_to_dataset)
        self.plot(df=plot_df,
                  x="PC1",
                  y="PC2",
                  hue="hue",
                  palette=self.palette,
                  xlabel="PC1 [{:.2f}%]".format(pca.explained_variance_ratio_[0] * 100),
                  ylabel="PC2 [{:.2f}%]".format(pca.explained_variance_ratio_[1] * 100),
                  title="PCA - eigenvectors",
                  filename="eigenvectors_plot{}".format(plot_appendix))

    def plot(self, df, x="x", y="y", hue=None, palette=None, xlabel=None,
             ylabel=None, title="", filename="PCA_plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, gridspec_kw={"width_ratios": [0.9, 0.1]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        s=100,
                        linewidth=0,
                        legend=None,
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

        if palette is not None:
            handles = []
            for label, color in palette.items():
                if label in df[hue].values.tolist():
                    handles.append(mpatches.Patch(color=color, label=label))
            ax2.legend(handles=handles, loc="center")

        #fig.savefig(os.path.join(self.plot_outdir, "{}.pdf".format(filename)))
        fig.savefig(os.path.join(self.plot_outdir, "{}.png".format(filename)))
        plt.close()

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Data: {}".format(self.data_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > N-components: {}".format(self.n_components))
        print("  > Plot output directory: {}".format(self.plot_outdir))
        print("  > File output directory: {}".format(self.file_outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
