#!/usr/bin/env python3

"""
File:         plot_decon_results_per_dataset.py
Created:      2022/01/17
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
import argparse
import itertools
import math
import os

# Third party imports.
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import upsetplot as up

# Local application imports.

# Metadata
__program__ = "Plot Decon-eQTL Results per Datasets"
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
./plot_decon_results_per_dataset.py -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/combine_gte_files/SampleToDataset.txt.gz -of 2021-12-07-CortexEUR-cis-InhibitorySummedWithOtherNeuron

./plot_decon_results_per_dataset.py -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz -of 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), "plot")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
            "Brainseq": "#C778A6",
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "OPC": "#009E73",
            "OPCs...COPs": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Pericytes": "#808080",
            "OtherNeuron": "#2690ce"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading std data.")
        std_df = self.load_file(self.std_path, header=0, index_col=None)

        dataset_sample_counts = list(zip(*np.unique(std_df.iloc[:, 1], return_counts=True)))
        dataset_sample_counts.sort(key=lambda x: -x[1])

        print("Loading decon-eQTL data")
        pvalue_data = {}
        fdr_data = {}
        beta_data = {}
        dataset_order = []
        colormap = {}
        for dataset, sample_size in dataset_sample_counts:
            dataset_result_path = os.path.join("decon_eqtl", self.outfolder + "-" + dataset, "deconvolutionResults.txt.gz")
            print(dataset_result_path)

            if os.path.exists(dataset_result_path):
                dataset_str = "{} [N={}]".format(dataset, sample_size)
                dataset_order.append(dataset_str)
                colormap[dataset_str] = self.palette[dataset]

                print("\tLoading results from dataset {}".format(dataset_str))
                decon_df = self.load_file(dataset_result_path)
                decon_fdr_df = decon_df.loc[:, [x for x in decon_df.columns if x.endswith("_FDR")]]

                decon_fdr_df = (decon_fdr_df <= 0.05).astype(int)
                decon_fdr_df.columns = [x.split("_")[0] for x in decon_fdr_df.columns]
                for ct in decon_fdr_df.columns:
                    fdr_df = decon_fdr_df.loc[:, [ct]]
                    fdr_df.columns = [dataset_str]
                    if ct in fdr_data.keys():
                        fdr_data[ct].append(fdr_df)
                    else:
                        fdr_data[ct] = [fdr_df]

                decon_pval_df = decon_df.loc[:, [x for x in decon_df.columns if x.endswith("_pvalue")]]
                decon_pval_df.columns = [x.split("_")[0] for x in decon_pval_df.columns]
                for ct in decon_pval_df.columns:
                    pvalue_df = decon_pval_df.loc[:, [ct]]
                    pvalue_df.columns = [dataset_str]
                    if ct in pvalue_data.keys():
                        pvalue_data[ct].append(pvalue_df)
                    else:
                        pvalue_data[ct] = [pvalue_df]

                counts = {}
                for column in decon_fdr_df.columns:
                    counts[column] = set(decon_fdr_df.loc[decon_fdr_df[column] == 1, :].index)
                self.plot_upsetplot(data=counts,
                                    filename="{}_ieQTLs_upsetplot_{}".format(self.outfolder, dataset))

                decon_beta_df = decon_df.loc[:, [x for x in decon_df.columns if x.endswith(":GT")]]
                decon_beta_df.columns = [x.split("_")[0].replace(":GT", "") for x in decon_beta_df.columns]
                for ct in decon_beta_df.columns:
                    beta_df = decon_beta_df.loc[:, [ct]]
                    beta_df.columns = [dataset_str]
                    if ct in beta_data.keys():
                        beta_data[ct].append(beta_df)
                    else:
                        beta_data[ct] = [beta_df]

        # for ct, df_list in beta_data.items():
        #     df = pd.concat(df_list, axis=1)
        #     print(df)
        #     exit()

        # for ct, df_list in pvalue_data.items():
        #     print("\t {}".format(ct))
        #
        #     df = pd.concat(df_list, axis=1)
        #     df.to_excel("{}_{}.xlsx".format(self.outfolder, ct))
        #     exit()

        print("Plotting upsetplots.")
        signif_df_list = []
        for ct, df_list in fdr_data.items():
            print("\t {}".format(ct))

            df = pd.concat(df_list, axis=1)
            # df.to_excel("{}_{}.xlsx".format(self.outfolder, ct))
            # continue

            counts = {}
            for column in df.columns:
                counts[column] = set(df.loc[df[column] == 1, :].index)

            self.plot_upsetplot(data=counts,
                                filename="{}_ieQTLs_upsetplot_{}".format(self.outfolder, ct))

            signif_df = df.sum(axis=0).to_frame()
            signif_df.columns = ["N"]
            signif_df["dataset"] = signif_df.index
            signif_df["cell type"] = ct
            signif_df_list.append(signif_df)
        # exit()

        signif_df = pd.concat(signif_df_list, axis=0)
        print(signif_df)

        dataset_sums_df = signif_df.groupby(['dataset']).sum()
        dataset_sums_df["group"] = ""
        dataset_sums_df["dataset"] = dataset_sums_df.index
        dataset_sums_df = dataset_sums_df.loc[dataset_order, :]
        print(dataset_sums_df)
        print(dataset_order)

        ct_sums_df = signif_df.groupby(['cell type']).sum()
        ct_sums = [[index, value[0]] for index, value in ct_sums_df.iterrows()]
        ct_sums_df["group"] = ""
        ct_sums_df["cell type"] = ct_sums_df.index
        ct_sums.sort(key=lambda x: -x[1])
        ct_order = [x[0] for x in ct_sums]
        ct_sums_df = ct_sums_df.loc[ct_order, :]
        print(ct_sums_df)
        print(ct_order)

        print("Plotting")
        self.plot_barplot(
            df=dataset_sums_df,
            group_column="group",
            groups=[""],
            x="N",
            y="dataset",
            xlabel="#ieQTLs (FDR<0.05)",
            ylabel="",
            palette=colormap,
            filename="{}_ieQTLs_barplot_per_ct_sum".format(self.outfolder)
        )
        self.plot_barplot(
            df=ct_sums_df,
            group_column="group",
            groups=[""],
            x="N",
            y="cell type",
            xlabel="#ieQTLs (FDR<0.05)",
            ylabel="",
            palette=self.palette,
            filename="{}_ieQTLs_barplot_per_dataset_sum".format(self.outfolder)
        )
        self.plot_barplot(
            df=signif_df,
            group_column="cell type",
            groups=ct_order,
            x="N",
            y="dataset",
            xlabel="#ieQTLs (FDR<0.05)",
            ylabel="",
            palette=colormap,
            filename="{}_ieQTLs_barplot_per_ct".format(self.outfolder)
        )
        self.plot_barplot(
            df=signif_df,
            group_column="dataset",
            groups=dataset_order,
            x="N",
            y="cell type",
            xlabel="#ieQTLs (FDR<0.05)",
            ylabel="",
            palette=self.palette,
            filename="{}_ieQTLs_barplot_per_dataset".format(self.outfolder)
        )

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def plot_upsetplot(self, data, filename):
        counts = self.count(data)
        counts = counts[counts > 0]
        up.plot(counts, sort_by='cardinality', show_counts=True)
        plt.savefig(os.path.join(self.outdir, "{}.png".format(filename)))
        plt.close()

    @staticmethod
    def count(input_data):
        combinations = []
        cols = list(input_data.keys())
        for i in range(1, len(cols) + 1):
            combinations.extend(list(itertools.combinations(cols, i)))

        indices = []
        data = []
        for combination in combinations:
            index = []
            for col in cols:
                if col in combination:
                    index.append(True)
                else:
                    index.append(False)

            background = set()
            for key in cols:
                if key not in combination:
                    work_set = input_data[key].copy()
                    background.update(work_set)

            overlap = None
            for key in combination:
                work_set = input_data[key].copy()
                if overlap is None:
                    overlap = work_set
                else:
                    overlap = overlap.intersection(work_set)

            duplicate_set = overlap.intersection(background)
            length = len(overlap) - len(duplicate_set)

            indices.append(index)
            data.append(length)

        s = pd.Series(data, index=pd.MultiIndex.from_tuples(indices, names=cols))
        s.name = "value"
        return s

    def plot_barplot(self, df, group_column, groups, x="x", y="y", xlabel="",
                         ylabel="", palette=None, filename=""):
        if df.shape[0] <= 2:
            return

        nplots = len(groups)
        ncols = math.ceil(np.sqrt(nplots))
        nrows = math.ceil(nplots / ncols)

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharey="all",
                                 figsize=(12 * ncols, 12 * nrows))
        sns.set(color_codes=True)

        row_index = 0
        col_index = 0
        for i in range(ncols * nrows):
            if nrows == 1 and ncols == 1:
                ax = axes
            elif nrows == 1 and ncols > 1:
                ax = axes[col_index]
            elif nrows > 1 and ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            if i < nplots:
                plot_df = df.loc[df[group_column] == groups[i], :].copy()
                plot_df.dropna(inplace=True)

                sns.despine(fig=fig, ax=ax)

                g = sns.barplot(x=x,
                                y=y,
                                hue=y,
                                palette=palette,
                                dodge=False,
                                data=plot_df,
                                orient="h",
                                ax=ax)
                g.legend_.remove()

                tmp_xlabel = ""
                if row_index == (nrows - 1):
                    tmp_xlabel = xlabel
                ax.set_xlabel(tmp_xlabel,
                              fontsize=20,
                              fontweight='bold')
                tmp_ylabel = ""
                if col_index == 0:
                    tmp_ylabel = ylabel
                ax.set_ylabel(tmp_ylabel,
                              fontsize=20,
                              fontweight='bold')

                ax.set_title(groups[i],
                             fontsize=25,
                             fontweight='bold')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        plt.tight_layout()
        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved: {}".format(outpath))

    def print_arguments(self):
        print("Arguments:")
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Output folder: {}".format(self.outfolder))
        print("  > Output directory {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
