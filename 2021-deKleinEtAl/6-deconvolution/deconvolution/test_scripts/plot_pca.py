#!/usr/bin/env python3

"""
File:         plot_pca.py
Created:      2021/06/30
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
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Plot PCA"
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
./plot_pca.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2020-11-20-decon-QTL/cis/cortex/expression_table/2020-07-16-MetaBrainDeconQtlGenes.TMM.SampSelect.ZeroVarRemov.covRemoved.expAdded.txt -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR -p GTE-EUR- -o CortexEUR_noENA_noGVEX

./plot_pca.py -d ../../../../2020-11-10-DeconOptimizer/2020-11-10-decon-optimizer/preprocess_scripts/CortexEUR_noENA_noGVEX/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.CovariatesRemovedOLS.ProbesCentered.SamplesZTransformed.ExpAdded.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR -p GTE-EUR- -o CortexEUR_noENA_noGVEX
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.gte_prefix = getattr(arguments, 'gte_prefix')
        self.gte_path = getattr(arguments, 'gene_to_expression')
        self.outfile = getattr(arguments, 'outfile')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.file_cohort_dict = {
            "AMPAD-MAYO-V2": "MAYO",
            "CMC_HBCC_set2": "CMC HBCC",
            "GTEx": "GTEx",
            "AMPAD-ROSMAP-V2": "ROSMAP",
            "BrainGVEX-V2": "Brain GVEx",
            "TargetALS": "Target ALS",
            "AMPAD-MSBB-V2": "MSBB",
            "NABEC-H610": "NABEC",
            "LIBD_1M": "LIBD",
            "ENA": "ENA",
            "LIBD_h650": "LIBD",
            "GVEX": "GVEX",
            "NABEC-H550": "NABEC",
            "CMC_HBCC_set3": "CMC HBCC",
            "UCLA_ASD": "UCLA ASD",
            "CMC": "CMC",
            "CMC_HBCC_set1": "CMC HBCC"
        }

        self.palette = {
            "MAYO": "#9c9fa0",
            "CMC HBCC": "#0877b4",
            "GTEx": "#0fa67d",
            "ROSMAP": "#6950a1",
            "Brain GVEx": "#48b2e5",
            "Target ALS": "#d5c77a",
            "MSBB": "#5cc5bf",
            "NABEC": "#6d743a",
            "LIBD": "#e49d26",
            "ENA": "#d46727",
            "GVEX": "#000000",
            "UCLA ASD": "#f36d2a",
            "CMC": "#eae453",
            "NA": "#808080"
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
        parser.add_argument("-gte",
                            "--gene_to_expression",
                            type=str,
                            required=True,
                            help="The path to the gene-expression link files.")
        parser.add_argument("-p",
                            "--gte_prefix",
                            type=str,
                            required=True,
                            help="The gene-expression link file prefix.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=False,
                            default="",
                            help="The name of the output file.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # Loading samples.
        print("Loading samples.")
        gte_combined_df = None
        cohort_to_sample = {}
        for infile in glob.glob(os.path.join(self.gte_path, "{}*.txt".format(self.gte_prefix))):
            file = os.path.basename(infile).replace(".txt", "").replace(self.gte_prefix, "")
            gte_df = self.load_file(infile, header=None, index_col=None)
            gte_df["file"] = file
            if gte_combined_df is None:
                gte_combined_df = gte_df
            else:
                gte_combined_df = pd.concat([gte_combined_df, gte_df], axis=0, ignore_index=True)

            cohort_to_sample[file] = gte_df.iloc[:, 1]
        gte_combined_df["cohort"] = gte_combined_df.iloc[:, 2].map(self.file_cohort_dict)
        sample_to_cohort = dict(zip(gte_combined_df.iloc[:, 1], gte_combined_df.iloc[:, 3]))

        print("Loading file.")
        df = self.load_file(self.data_path, header=0, index_col=0)
        df = df.dropna(axis=1, how='all')
        print(df)

        print("Performing PCA.")
        centered = df - df.mean(axis=0)
        zscores = centered / centered.std(axis=0)
        pca = PCA(n_components=2)
        pca.fit(zscores)
        components_df = pd.DataFrame(pca.components_)
        components_df.index = ["Comp{}".format(i + 1) for i, _ in
                               enumerate(components_df.index)]
        components_df.columns = df.columns
        print(components_df)

        print("Plotting.")
        plot_df = components_df.T
        plot_df["cohort"] = plot_df.index.map(sample_to_cohort)
        plot_df["cohort"] = plot_df["cohort"].fillna('NA')
        self.plot(df=plot_df, x="Comp1", y="Comp2", hue="cohort",
                  palette=self.palette,
                  xlabel="PC1 [{:.2f}%]".format(pca.explained_variance_ratio_[0] * 100),
                  ylabel="PC2 [{:.2f}%]".format(pca.explained_variance_ratio_[1] * 100),
                  title="PCA - eigenvectors",
                  filename="PCA_plot_{}".format(self.outfile))

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
    def reverse_dict(dict):
        out_dict = {}
        seen_keys = set()
        for key, value in dict.items():
            if key in seen_keys:
                print("Key {} has muiltiple values.".format(key))
            seen_keys.add(key)

            if value in out_dict.keys():
                keys = out_dict[value]
                keys.append(key)
                out_dict[value] = keys
            else:
                out_dict[value] = [key]

        return out_dict

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

        fig.savefig(os.path.join(self.outdir, "{}.png".format(filename)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Data: {}".format(self.data_path))
        print("  > GtE path: {}".format(self.gte_path))
        print("  >   GtE prefix: {}".format(self.gte_prefix))
        print("  > Output directory {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
