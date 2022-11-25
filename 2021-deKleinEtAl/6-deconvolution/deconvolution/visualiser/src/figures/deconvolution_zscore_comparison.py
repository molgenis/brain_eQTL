"""
File:         deconvolution_zscore_comparison.py
Created:      2020/04/07
Last Changed: 2020/06/03
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
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.
from general.utilities import prepare_output_dir, p_value_to_symbol


class DeconvolutionZscoreComparison:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir,
                                   'deconvolution_zscore_comparison')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.inter_df = dataset.get_inter_cov_zscore_df()
        self.celltypes = dataset.get_celltypes()
        self.cellmap_methods = dataset.get_cellmap_methods()
        self.marker_genes = dataset.get_marker_genes()
        self.color_map = dataset.get_colormap()

    def start(self):
        print("Plotting CellMap Methods vs Marker genes.")
        self.print_arguments()

        for i, method1 in enumerate(self.cellmap_methods):
            method1_df = self.inter_df.loc[self.inter_df.index.str.startswith(method1[0]), :]

            # Plot vs other methods.
            for j, method2 in enumerate(self.cellmap_methods):
                method2_df = self.inter_df.loc[self.inter_df.index.str.startswith(method2[0]), :]
                if i < j:
                    print("{} vs {}".format(method1[0], method2[0]))
                    self.compare_cellmap_method(method1, method1_df, method2,
                                                method2_df, self.celltypes,
                                                self.color_map, self.outdir,
                                                self.extension)

            # Plot vs marker genes.
            marker_df = self.inter_df.loc[self.inter_df.index.str.startswith(self.marker_genes), :]
            print("{} vs {}MarkerGenes".format(method1[0], self.marker_genes))
            self.compare_cellmap_with_marker_genes(method1, method1_df,
                                                   self.marker_genes, marker_df,
                                                   self.celltypes, self.color_map,
                                                   self.outdir, self.extension)

    @staticmethod
    def compare_cellmap_method(method1, method1_df, method2, method2_df,
                               celltypes, color_map, outdir, extension):
        nplots = len(method1_df.index)

        ncols = 1
        nrows = math.ceil(nplots / ncols)

        method1_name = method1[0].split("_")[0]
        method2_name = method2[0].split("_")[0]

        sns.set(rc={'figure.figsize': (12 * ncols, 9 * nrows)})
        sns.set_style("ticks")
        fig = plt.figure()
        grid = fig.add_gridspec(ncols=ncols, nrows=nrows)

        row_index = 0
        col_index = 0

        for celltype in celltypes:
            method_celltype = celltype
            if method_celltype == "Microglia":
                method_celltype = "Macrophage"

            # Get the data.
            method1_data = method1_df.loc[method1[0] + method_celltype + method1[1], :]
            method2_data = method2_df.loc[method2[0] + method_celltype + method2[1], :]

            df = pd.DataFrame({method1_name: method1_data, method2_name: method2_data})
            df.dropna(inplace=True)

            # Get the color.
            color = color_map[celltype]

            # Creat the subplot.
            ax = fig.add_subplot(grid[row_index, col_index])
            sns.despine(fig=fig, ax=ax)

            coef_str = "NA"
            p_str = "NA"
            if len(df.index) > 1:
                # Calculate the correlation.
                coef, p = stats.spearmanr(df[method1_name], df[method2_name])
                coef_str = "{:.2f}".format(coef)
                # p_str = p_value_to_symbol(p)
                p_str = "p = {:.2e}".format(p)

                # Plot.
                g = sns.regplot(x=method1_name,
                                y=method2_name,
                                data=df,
                                scatter_kws={'facecolors': '#000000',
#                                             'edgecolor': '#000000',
                                             'linewidth': 0,
                                             'alpha': 0.5},
                                line_kws={"color": color},
                                ax=ax
                                )

            # Add the text.
            ax.annotate(
                'r = {} [{}]'.format(coef_str, p_str),
                xy=(0.03, 0.94),
                xycoords=ax.transAxes,
                color=color,
                fontsize=18,
                fontweight='bold')

            ax.text(0.5, 1.05, method_celltype,
                    fontsize=26, weight='bold', ha='center', va='bottom',
                    color=color,
                    transform=ax.transAxes)

            ax.set_ylabel((method1[0] + method1[1]).replace("_", " "),
                          fontsize=18,
                          fontweight='bold')
            ax.set_xlabel((method2[0] + method2[1]).replace("_", " "),
                          fontsize=18,
                          fontweight='bold')

            ax.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
            ax.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)

            # Increment indices.
            col_index += 1
            if col_index == ncols:
                col_index = 0
                row_index += 1

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}_vs_{}.{}".format(method1_name.split("_")[0],
                                                      method2_name.split("_")[0],
                                                      extension)))
        plt.close()


    @staticmethod
    def compare_cellmap_with_marker_genes(method, method_df, marker_genes_prefix,
                                          marker_df, celltypes, color_map,
                                          outdir, extension):
        nplots = len(celltypes)

        method_name = method[0].split("_")[0]

        sns.set(rc={'figure.figsize': (12 * nplots, 9 * nplots)})
        sns.set_style("ticks")
        fig = plt.figure()
        grid = fig.add_gridspec(ncols=nplots, nrows=nplots)

        row_index = 0
        col_index = 0

        for celltype in celltypes:
            method_celltype = celltype
            if method_celltype == "Microglia":
                method_celltype = "Macrophage"

            # Get the method data.
            method_data = method_df.loc[method[0] + method_celltype + method[1], :]

            # Get the marker gene data.
            marker_data = marker_df.loc[marker_df.index.str.startswith(marker_genes_prefix + celltype), :]

            # Get the color.
            color = color_map[celltype]

            # Loop over each marker gene.
            for marker_gene, marker_data in marker_data.iterrows():
                df = pd.DataFrame({method_name: method_data, marker_gene: marker_data})
                df.dropna(inplace=True)

                # Creat the subplot.
                ax = fig.add_subplot(grid[row_index, col_index])
                sns.despine(fig=fig, ax=ax)

                coef_str = "NA"
                p_str = "NA"
                if len(df.index) > 1:
                    # Calculate the correlation.
                    coef, p = stats.spearmanr(df[method_name], df[marker_gene])
                    coef_str = "{:.2f}".format(coef)
                    # p_str = p_value_to_symbol(p)
                    p_str = "p = {:.2e}".format(p)

                    # Plot.
                    g = sns.regplot(x=method_name,
                                    y=marker_gene,
                                    data=df,
                                    scatter_kws={'facecolors': '#000000',
#                                                 'edgecolor': '#000000',
                                                 'linewidth': 0,
                                                 'alpha': 0.5},
                                    line_kws={"color": color},
                                    ax=ax
                                    )

                # Add the text.
                ax.annotate(
                    'r = {} [{}]'.format(coef_str, p_str),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color=color,
                    fontsize=18,
                    fontweight='bold')

                ax.text(0.5, 1.05,
                        '{} vs {} {}'.format(method_name.replace("_", " "),
                                             marker_gene.split("_")[0],
                                             marker_gene.split("_")[2]),
                        fontsize=26, weight='bold', ha='center', va='bottom',
                        color=color,
                        transform=ax.transAxes)

                ax.set_ylabel((method[0] + method_celltype + method[1]).replace("_", " "),
                              fontsize=18,
                              fontweight='bold')
                ax.set_xlabel(marker_gene,
                              fontsize=18,
                              fontweight='bold')

                ax.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
                ax.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)

                # Increment indices.
                col_index += 1
                if col_index == nplots:
                    col_index = 0
                    row_index += 1

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}_vs_{}"
                                 "MarkerGenes.{}".format(method_name.split("_")[0],
                                                         marker_genes_prefix,
                                                         extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > CellMap Methods: {}".format(self.cellmap_methods))
        print("  > Marker Genes: {}".format(self.marker_genes))
        print("  > Output directory: {}".format(self.outdir))
        print("")
