"""
File:         deconvolution_covariate_comparison.py
Created:      2020/04/15
Last Changed: 2020/06/19
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
import math
import os

# Third party imports.
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir, p_value_to_symbol


class DeconvolutionCovariateComparison:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'deconvolution_covariate_comparison')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.groups = dataset.get_groups()
        self.cellmap_methods = dataset.get_cellmap_methods()
        self.marker_genes = dataset.get_marker_genes()
        self.cov_df = dataset.get_cov_df()
        self.cmap = dataset.get_diverging_cmap()

    def start(self):
        print("Plotting deconvolution convariate comparison.")
        self.print_arguments()
        corr_df, pval_df = self.correlate(self.cov_df)

        prefixes = [x[0] for x in self.cellmap_methods]
        prefixes.append(self.marker_genes)
        prefixes = [x.replace("_", "") for x in prefixes]
        self.plot(corr_df, pval_df, self.groups, prefixes, self.cmap,
                  self.outdir, self.extension)

    @staticmethod
    def correlate(df):
        print("Creating the correlation matrix.")
        corr_df = pd.DataFrame(np.nan, index=df.index, columns=df.index)
        pval_df = pd.DataFrame("", index=df.index, columns=df.index)
        for i, row1 in enumerate(df.index):
            for j, row2 in enumerate(df.index):
                coef, p = stats.spearmanr(df.loc[row1, :], df.loc[row2, :])
                corr_df.loc[row1, row2] = coef
                # pval_df.loc[row1, row2] = p_value_to_symbol(p)
                pval_df.loc[row1, row2] = "{:.2e}".format(p)

        return corr_df, pval_df

    def plot(self, corr_df, pval_df, groups, prefixes, cmap, outdir, extension):
        print("Plotting")

        decon_groups = []
        for group in groups:
            for prefix in prefixes:
                if prefix in group[3]:
                    decon_groups.append(group)
                    break

        norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)

        heatmapkws = dict(square=False, cbar=False, cmap=cmap, fmt='',
                          linewidths=1.0, center=0, vmin=-1, vmax=1,
                          annot_kws={"size": 12, "color": "#808080"})

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=len(decon_groups), nrows=len(decon_groups),
                                 figsize=(29, 21))
        plt.subplots_adjust(left=0.2, right=0.87, bottom=0.2, top=0.95,
                            wspace=0.1, hspace=0.1)

        for i in range(len(decon_groups)):
            (ra, rb, ylabel, yremove) = decon_groups[i]
            for j in range(len(decon_groups)):
                (ca, cb, xlabel, xremove) = decon_groups[j]
                ax = axes[i, j]

                if i >= j:
                    print("\tPlotting axes[{}, {}]".format(i, j))
                    xticklabels = False
                    if i == len(decon_groups) - 1:
                        xticklabels = True
                    yticklabels = False
                    if j == 0:
                        yticklabels = True

                    data_subset = corr_df.iloc[ra:rb, ca:cb]
                    annot_subset = pval_df.iloc[ra:rb, ca:cb]

                    if xremove == yremove:
                        label = xremove.replace("_", "")
                        index_index = 1
                        if xremove.startswith("McKenzie"):
                            label = '_'.join(data_subset.index[0].split("_")[0:2])
                            index_index = 2
                        self.heatmap(data_subset.copy(), annot_subset.copy(), label, index_index, self.cmap, outdir, extension)

                    sns.heatmap(data_subset,
                                annot=annot_subset,
                                xticklabels=xticklabels,
                                yticklabels=yticklabels,
                                ax=ax,
                                **heatmapkws)

                    if xticklabels:
                        new_xticks = [x.replace(xremove, '').replace("_", " ")
                                      for x in data_subset.columns]
                        ax.set_xticklabels(new_xticks, fontsize=14, rotation=90)
                        ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')

                    if yticklabels:
                        new_yticks = [y.replace(yremove, '').replace("_", " ")
                                      for y in data_subset.index]
                        ax.set_yticklabels(new_yticks, fontsize=14, rotation=0)
                        ax.set_ylabel(ylabel, fontsize=16, fontweight='bold')
                else:
                    ax.set_axis_off()

        cax = fig.add_axes([0.9, 0.2, 0.01, 0.7])
        sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.align_ylabels(axes[:, 0])
        fig.align_xlabels(axes[len(decon_groups) - 1, :])
        fig.colorbar(sm, cax=cax)
        fig.suptitle('Deconvolution Covariate Correlations', fontsize=40,
                     fontweight='bold')
        fig.savefig(os.path.join(outdir,
                                 "deconvolution_covariate_comparison.{}".format(
                                     extension)))
        plt.close()

        mg = ["McKenzie_Astrocyte_AQP4", "McKenzie_EndothelialCell_KDR", "McKenzie_Microglia_TLR2", "McKenzie_Neuron_CCK", "McKenzie_Oligodendrocyte_CA2"]
        data_subset = corr_df.loc[mg, mg]
        annot_subset = pval_df.loc[mg, mg]
        label = "McKenzie_best"
        index_index = 1
        self.heatmap(data_subset, annot_subset, label, index_index, self.cmap, outdir, extension)
        exit()

        comparisons = [["McKenzie_Neuron_CCK", "CellMapPCA_Neuron_PC1", "CellMapNMF_Neuron_C1", "CellMapNNLS_Neuron"],
                       ["McKenzie_Oligodendrocyte_CA2", "CellMapPCA_Oligodendrocyte_PC1", "CellMapNMF_Oligodendrocyte_C1", "CellMapNNLS_Oligodendrocyte"],
                       ["McKenzie_EndothelialCell_KDR", "CellMapPCA_EndothelialCell_PC1", "CellMapNMF_EndothelialCell_C1", "CellMapNNLS_EndothelialCell"],
                       ["McKenzie_Microglia_TLR2", "CellMapPCA_Macrophage_PC1", "CellMapNMF_Macrophage_C1", "CellMapNNLS_Macrophage"],
                       ["McKenzie_Astrocyte_AQP4", "CellMapPCA_Astrocyte_PC1", "CellMapNMF_Astrocyte_C1", "CellMapNNLS_Astrocyte"]]
        for comparison, celltype in zip(comparisons, ["Neuron", "Oligodendrocyte", "EndothelialCell", "Microglia", "Astrocyte"]):
            data_subset = corr_df.loc[comparison, comparison]
            annot_subset = pval_df.loc[comparison, comparison]
            data_subset.index = ["Mcke", "PCA", "NMF", "NNLS"]
            data_subset.columns = ["Mcke", "PCA", "NMF", "NNLS"]
            annot_subset.index = ["Mcke", "PCA", "NMF", "NNLS"]
            annot_subset.columns = ["Mcke", "PCA", "NMF", "NNLS"]
            index_index = 0
            self.heatmap(data_subset, annot_subset, celltype, index_index,
                         self.cmap, outdir, extension)

    @staticmethod
    def heatmap(data_subset, annot_subset, label, index_index, cmap, outdir,
                extension):
        abbreviations = {"Neuron": "neuro", "Oligodendrocyte": "oligo",
                         "EndothelialCell": "endo", "Microglia": "micro",
                         "Macrophage": "macro", "Astrocyte": "astro",
                         "CellMapNNLS": "NNLS", "CellMapPCA_PC1": "PCA",
                         "CellMapNMF_C1": "NMF", "McKenzie": "McKe"}
        new_index = []
        for x in data_subset.index:
            x = x.split("_")[index_index]
            if x in abbreviations.keys():
                x = abbreviations[x]
            new_index.append(x)
        new_cols = []
        for x in data_subset.columns:
            x = x.split("_")[index_index]
            if x in abbreviations.keys():
                x = abbreviations[x]
            new_cols.append(x)

        data_subset.index = new_index
        data_subset.columns = new_cols

        for i in range(len(data_subset.index)):
            for j in range(len(data_subset.columns)):
                if i < j:
                    data_subset.iloc[i, j] = np.nan
                    annot_subset.iloc[i, j] = ""

        sns.set(color_codes=True)
        g = sns.clustermap(data_subset, cmap=cmap,
                           row_cluster=False, col_cluster=False,
                           yticklabels=True, xticklabels=True, square=True,
                           vmin=-1, vmax=1, annot=data_subset.round(2), fmt='',
                           annot_kws={"size": 16, "color": "#000000"},
                           figsize=(12, 12))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=25, rotation=0))
        plt.setp(
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(),
                                         fontsize=25, rotation=45))
        #g.cax.set_visible(False)
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "{}_clustermap.{}".format(label,
                                                                 extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Groups: {}".format(self.groups))
        print("  > CellMap Methods: {}".format(self.cellmap_methods))
        print("  > Marker Genes: {}".format(self.marker_genes))
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
