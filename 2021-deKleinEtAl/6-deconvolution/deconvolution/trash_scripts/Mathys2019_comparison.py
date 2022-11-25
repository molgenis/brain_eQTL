#!/usr/bin/env python3

"""
File:         mathys2019_comparison.py
Created:      2020/09/16
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
import os

# Third party imports.
import pandas as pd
from scipy import stats
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Mathys 2019 Comparison"
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
        self.mathys_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/41586_2019_1195_MOESM8_ESM.csv"
        self.decon_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-16-decon-eQTL/cis/cortex/decon_out/deconvolutionResults.csv"
        self.gene_info_path = "/groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Loading Mathys et al 2019 supplementary tabel 6.")
        mathys_df = pd.read_csv(self.mathys_path, sep=",", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.mathys_path),
                                      mathys_df.shape))
        mathys_df = mathys_df.loc[
                    (mathys_df["subpopulation"].str.startswith("Ex")) | (
                        mathys_df["subpopulation"].str.startswith("In")), :]
        mathys_df.drop_duplicates(subset=['gene.name'], inplace=True)
        print(mathys_df)

        print("Loading gene info matrix.")
        gene_info_df = pd.read_csv(self.gene_info_path, sep="\t", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.gene_info_path),
                                      gene_info_df.shape))
        gene_dict = dict(zip(gene_info_df["ArrayAddress"], gene_info_df["Symbol"]))

        print("Loading decon-eQTL cell type interaction results.")
        decon_df = pd.read_csv(self.decon_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon_path),
                                      decon_df.shape))
        print(decon_df)
        decon_df.reset_index(drop=False, inplace=True)
        decon_df = decon_df.melt(id_vars="index", value_vars=[x for x in decon_df.columns if x.endswith("pvalue")])
        decon_df["FDR"] = multitest.multipletests(decon_df['value'], method='fdr_bh')[1]
        decon_df = decon_df.loc[decon_df['variable'] == "Neuron_pvalue", :]
        decon_df['Symbol'] = [gene_dict[x.split("_")[0]] for x in decon_df["index"]]
        decon_df = decon_df.groupby(decon_df["Symbol"], group_keys=False).apply(lambda x: x.loc[x["FDR"].idxmin()])
        decon_df.reset_index(drop=True, inplace=True)
        print(decon_df)
        print("Signif. interaction: {}, No signif. interaction: {}".format(decon_df.loc[decon_df["FDR"] < 0.05, :].shape[0], decon_df.loc[decon_df["FDR"] >= 0.05, :].shape[0]))

        df = mathys_df.merge(decon_df, left_on="gene.name", right_on="Symbol")
        print(df)

        size = df.loc[:, "subpopulation"].value_counts()
        interaction = df.loc[df["FDR"] < 0.05, "subpopulation"].value_counts()
        no_interaction = df.loc[df["FDR"] >= 0.05, "subpopulation"].value_counts()
        celltypes = list(df["subpopulation"].unique())

        total_inter = interaction.sum()
        total_no_inter = no_interaction.sum()

        data = []
        indices = []
        for ct in celltypes:
            n = 0
            if ct in size:
                n = size[ct]

            n_inter = 0
            if ct in interaction:
                n_inter = interaction[ct]

            n_no_inter = 0
            if ct in no_interaction:
                n_no_inter = no_interaction[ct]

            oog_interaction = (total_inter - n_inter)
            oog_no_interaction = (total_no_inter - n_no_inter)

            oddsratio, pvalue = stats.fisher_exact([[n_inter, oog_interaction], [n_no_inter, oog_no_interaction]])

            data.append([n, n_inter, n_no_inter, oog_interaction,
                         oog_no_interaction, oddsratio, pvalue])
            indices.append(ct)

        result_df = pd.DataFrame(data,
                                 columns=["N", "interaction",
                                          "no interaction",
                                          "out-of-group interaction",
                                          "out-of-group no interaction",
                                          "oddsratio", "pvalue"],
                                 index=indices)
        result_df["signif."] = result_df["oddsratio"] < 0.05
        result_df.sort_values(by="pvalue", ascending=True, inplace=True)
        print(result_df)

        result_df.to_csv(os.path.join(self.outdir, "Mathys2019Comparison.txt.gz"), sep="\t", header=True, index=True, compression="gzip")



# def start(self):
#     print("Loading Mathys et al 2019 supplementary tabel 6.")
#     mathys_df = pd.read_csv(self.mathys_path, sep=",", header=0)
#     print("\tLoaded dataframe: {} "
#           "with shape: {}".format(os.path.basename(self.mathys_path),
#                                   mathys_df.shape))
#     mathys_df = mathys_df.loc[(mathys_df["subpopulation"].str.startswith("Ex")) | (mathys_df["subpopulation"].str.startswith("In")), :]
#     mathys_df.drop_duplicates(subset=['gene.name'], inplace=True)
#     print(mathys_df)
#
#     print("Loading gene info matrix.")
#     gene_info_df = pd.read_csv(self.gene_info_path, sep="\t", header=0)
#     print("\tLoaded dataframe: {} "
#           "with shape: {}".format(os.path.basename(self.gene_info_path),
#                                   gene_info_df.shape))
#
#     data_df = mathys_df.merge(gene_info_df, left_on="gene.name", right_on="Symbol")
#     del mathys_df, gene_info_df
#     print(data_df)
#
#     print("Loading decon-eQTL cell type interaction results.")
#     decon_df = pd.read_csv(self.decon_path, sep="\t", header=0, index_col=0)
#     print("\tLoaded dataframe: {} "
#           "with shape: {}".format(os.path.basename(self.decon_path),
#                                   decon_df.shape))
#     decon_df['ArrayAddress'], decon_df['RsID'] = decon_df.index.str.split('_', 1).str
#     decon_df = decon_df.loc[:, ['Neuron_pvalue', 'Beta4_Neuron', 'Beta9_Neuron:GT', 'ArrayAddress', 'RsID']]
#     decon_df["Neuron_FDR"] = multitest.multipletests(decon_df['Neuron_pvalue'], method='fdr_bh')[1]
#     print(decon_df)
#
#     df = data_df.merge(decon_df, left_on="ArrayAddress", right_on="ArrayAddress")
#     del data_df, decon_df
#     print(df)
#
#     size = df.loc[:, "subpopulation"].value_counts()
#     unique = df.drop_duplicates(subset=['gene.name'])["subpopulation"].value_counts()
#     interaction = df.loc[df["Neuron_FDR"] < 0.05, "subpopulation"].value_counts()
#     no_interaction = df.loc[df["Neuron_FDR"] >= 0.05, "subpopulation"].value_counts()
#     celltypes = list(df["subpopulation"].unique())
#
#     total_inter = interaction.sum()
#     total_no_inter = no_interaction.sum()
#
#     data = []
#     indices = []
#     for ct in celltypes:
#         n = 0
#         if ct in size:
#             n = size[ct]
#
#         n_unique = 0
#         if ct in unique:
#             n_unique = unique[ct]
#
#         n_inter = 0
#         if ct in interaction:
#             n_inter = interaction[ct]
#
#         n_no_inter = 0
#         if ct in no_interaction:
#             n_no_inter = no_interaction[ct]
#
#         oog_interaction = (total_inter - n_inter)
#         oog_no_interaction = (total_no_inter - n_no_inter)
#
#         oddsratio, pvalue = stats.fisher_exact([[n_no_inter, oog_no_interaction], [n_inter, oog_interaction]])
#
#         data.append([n, n_unique, n_inter, n_no_inter, oog_interaction, oog_no_interaction, oddsratio, pvalue])
#         indices.append(ct)
#
#     result_df = pd.DataFrame(data,
#                              columns=["N", "N Unique", "interaction", "no interaction", "out-of-group interaction", "out-of-group no interaction", "oddsratio", "pvalue"],
#                              index=indices)
#     result_df["signif."] = result_df["oddsratio"] < 0.05
#     result_df.sort_index(inplace=True)
#     print(result_df)


if __name__ == '__main__':
    m = main()
    m.start()
