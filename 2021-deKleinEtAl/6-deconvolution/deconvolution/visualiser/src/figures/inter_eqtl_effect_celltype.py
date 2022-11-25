"""
File:         inter_eqtl_effect_celltype.py
Created:      2020/04/20
Last Changed: 2020/06/16
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
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib
import itertools
matplotlib.use('Agg')
import upsetplot as up
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.
from general.utilities import prepare_output_dir


class IntereQTLEffectCelltype:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_effect_celltype')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.eqtlinter_df = dataset.get_eqtl_and_interactions_df()
        self.celltypes = dataset.get_celltypes()
        self.cellmap_methods = dataset.get_cellmap_methods()
        self.marker_genes = dataset.get_marker_genes()
        self.z_score_cutoff = dataset.get_significance_cutoff()
        self.colormap = dataset.get_colormap()

    def start(self):
        print("Plotting cell type mediated interaction eQTLs")
        self.print_arguments()

        total = len(set(self.eqtlinter_df.index.values))

        print("Plotting deconvolution methods")
        overview_data = {}
        overview_overlap = {}
        celltype_mediated_eqtls = set()
        for (prefix, suffix) in self.cellmap_methods:
            name = prefix.replace("_", "") + suffix
            method_eqtls = set()
            data = {}

            for celltype in self.celltypes:
                method_celltype = celltype
                if method_celltype == "Microglia":
                    method_celltype = "Macrophage"

                df = self.eqtlinter_df.loc[:,
                     prefix + method_celltype + suffix].copy()
                eqtls = set(df.loc[df == 1].index.values)
                method_eqtls.update(eqtls)
                data[method_celltype] = eqtls
            celltype_mediated_eqtls.update(method_eqtls)

            # Plot.
            overlap_vals = self.get_overlap(data)
            self.upsetplot(data, prefix.replace("_", "") + suffix,
                           self.outdir, self.extension)
            self.plot_pie(total, len(method_eqtls), name, self.outdir,
                          self.extension)

            overview_data[name] = data
            overview_overlap[name] = overlap_vals

        print("Plotting {}".format(self.marker_genes))
        name = self.marker_genes.replace("_", "")
        data = {}
        for celltype in self.celltypes:
            df = self.eqtlinter_df.loc[:,
                 self.eqtlinter_df.columns.str.startswith(
                     self.marker_genes + celltype + "_")].copy()
            best_marker = df.sum(axis=0).idxmax()
            eqtls = set(df.loc[df[best_marker] == 1].index.values)
            data[celltype] = eqtls

            eqtls = set(df.loc[df.sum(axis=1) > 0].index.values)
            celltype_mediated_eqtls.update(eqtls)

        overlap_vals = self.get_overlap(data)
        self.upsetplot(data, name, self.outdir, self.extension)
        overview_data[name] = data
        overview_overlap[name] = overlap_vals

        print("Plotting overview data.")
        self.plot_overview(overview_data, overview_overlap, self.celltypes,
                           self.outdir, self.extension, self.colormap)

        print("Plotting all eQTLs.")
        part = len(celltype_mediated_eqtls)
        self.plot_pie(total, part, "all_cis-eQTLs", self.outdir, self.extension)

        print("Plotting {} marker genes separately".format(self.marker_genes))
        for celltype in self.celltypes:
            df = self.eqtlinter_df.loc[:,
                 self.eqtlinter_df.columns.str.startswith(
                     self.marker_genes + celltype + "_")].copy()
            data = {}
            for column in df:
                gene_name = column.split("_")[2]
                eqtls = set(df.loc[df[column] == 1, column].index.values)
                data[gene_name] = eqtls
            self.upsetplot(data, "{}{}".format(self.marker_genes, celltype),
                           self.outdir, self.extension)

        print("Plotting all methods combined")
        celltype_data = {}
        for celltype in self.celltypes:
            data = {}

            # Add the marker genes.
            marker_df = self.eqtlinter_df.loc[:,
                        self.eqtlinter_df.columns.str.startswith(
                            self.marker_genes + celltype + "_")].copy()
            best_marker = marker_df.sum(axis=0).idxmax()
            eqtls = set(marker_df.loc[marker_df[best_marker] == 1].index.values)
            data[self.marker_genes.replace("_", "")] = eqtls

            method_celltype = celltype
            title = celltype
            if method_celltype == "Microglia":
                method_celltype = "Macrophage"
                title = method_celltype + "_" + celltype

            # Add the deconvolution methods.
            for (prefix, suffix) in self.cellmap_methods:
                df = self.eqtlinter_df.loc[:,
                     prefix + method_celltype + suffix].copy()
                eqtls = set(df.loc[df == 1].index.values)
                data[prefix.replace("_", "") + suffix] = eqtls

            # Plot.
            self.upsetplot(data, title, self.outdir, self.extension)

            all_eQTLs_of_celltype = []
            for value in data.values():
                all_eQTLs_of_celltype.extend(value)
            celltype_data[celltype] = set(all_eQTLs_of_celltype)

        print("Plotting all celltypes combined")
        self.upsetplot(celltype_data, "Cell-Type_eQTLs", self.outdir, self.extension)

    @staticmethod
    def get_overlap(data):
        overlap = pd.DataFrame(np.nan, index=data.keys(), columns=data.keys())
        overlap_values = []
        for i, key1 in enumerate(overlap.index):
            for j, key2 in enumerate(overlap.columns):
                value1 = data[key1]
                value2 = data[key2]

                if len(value1) == 0 or len(value2) == 0:
                    continue

                n = len(set(value1).intersection(set(value2)))
                fraction = round(n / len(value1), 2)
                overlap.iloc[i, j] = fraction

                if i != j:
                    overlap_values.append([key1, key2, fraction])
        return overlap_values

    @staticmethod
    def plot_overview(data, overlap, celltypes, outdir, extension, colormap):
        boxplot_data = []
        for key, value in overlap.items():
            for row in value:
                boxplot_data.append([key] + row)
        boxplot_df = pd.DataFrame(boxplot_data,
                                  columns=["method", "from", "to", "overlap"])

        group_df = boxplot_df.copy()
        group_df = group_df.groupby("method").median()
        group_df.columns = ['median']
        group_df.index.name = None
        group_df.sort_values(['median'], inplace=True)

        boxplot_df["hue"] = "other"
        boxplot_df.loc[boxplot_df["method"] == group_df.index[0], "hue"] = "best"

        barplot_data = {}
        methods = data.keys()
        method_sum = {method: 0 for method in methods}
        for i, celltype in enumerate(celltypes):
            celltype_data = []
            for method in methods:
                method_celltype = celltype
                if celltype not in data[method].keys():
                    method_celltype = "Macrophage"

                old_sum = method_sum[method]
                new_sum = old_sum + len(data[method][method_celltype])
                method_sum[method] = new_sum
                celltype_data.append([method,
                                      new_sum,
                                      colormap[method_celltype]])
            barplot_data[i] = celltype_data

        sns.set(rc={'figure.figsize': (18, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(1, 2)
        plt.subplots_adjust(top=0.9, bottom=0.1, wspace=0.01, hspace=0.01)
        for ax in [ax1, ax2]:
            sns.despine(fig=fig, ax=ax)

        # Plot the bar graph.
        max_val = max(method_sum.values())
        major_ticks = 10 ** (math.floor(math.log10(max(max_val, 100))))
        minor_ticks = int(major_ticks / 2)
        for i in range(0, int(max_val) + (1 * major_ticks), minor_ticks):
            alpha = 0.025
            if i % major_ticks == 0:
                alpha = 0.15
            ax1.axhline(i, ls='-', color="#000000", alpha=alpha, zorder=-1)

        for i in range(len(celltypes) - 1, -1, -1):
            barplot_df = pd.DataFrame(barplot_data[i], columns=["method", "len", "color"])
            sns.barplot(x="method", y="len", palette=barplot_df["color"],
                        data=barplot_df, order=group_df.index, ax=ax1)

        ax1.text(0.5, 1.02,
                    'Total Interactions per Celltype',
                    fontsize=20, weight='bold', ha='center', va='bottom',
                    transform=ax1.transAxes)
        ax1.set_xlabel('',
                       fontsize=14,
                       fontweight='bold')
        ax1.set_ylabel('N cis-eQTL interactions',
                        fontsize=14,
                        fontweight='bold')

        handles = []
        for celltype in celltypes[::-1]:
            label = celltype
            if celltype == "Microglia":
                label = "Microglia/Macrophage"
            handles.append(mpatches.Patch(color=colormap[celltype], label=label))
        ax1.legend(handles=handles, loc="upper left")
        ax1.tick_params(labelsize=14)

        # Plot boxplot.
        for i in range(0, 100, 5):
            alpha = 0.025
            if i % 10 == 0:
                alpha = 0.15
            ax2.axhline(i / 100, ls='-', color="#000000", alpha=alpha,
                        zorder=-1)

        sns.boxplot(x="method", y="overlap", hue="hue", data=boxplot_df,
                    palette={"best": "#0072B2", "other": "#808080"},
                    order=group_df.index, showfliers=False, dodge=False,
                    ax=ax2)
        sns.swarmplot(x="method", y="overlap", data=boxplot_df, color=".25",
                      order=group_df.index, ax=ax2)
        ax2.get_legend().set_visible(False)

        ax2.text(0.5, 1.02,
                 'Average Overlap in Interactions\nBetween Cell Types',
                 fontsize=20, weight='bold', ha='center', va='bottom',
                 transform=ax2.transAxes)
        ax2.set_ylabel('overlap',
                       fontsize=14,
                       fontweight='bold')
        ax2.set_xlabel('',
                       fontsize=14,
                       fontweight='bold')
        ax2.tick_params(labelsize=14)
        ax2.set_ylim(0, 1)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "overview.{}".format(extension)))
        plt.close()

    def upsetplot(self, data, title, outdir, extension):
        counts = self.count(data)
        counts = counts[counts > 0]
        up.plot(counts, sort_by='cardinality', show_counts=True)
        plt.suptitle('{}'.format(title.replace("_", " ")), fontsize=18, fontweight='bold')
        plt.savefig(os.path.join(outdir, "{}_upsetplot.{}".format(title, extension)))
        plt.close()

    @staticmethod
    def count(input_data):
        combinations = []
        cols = list(input_data.keys())
        for i in range(1, len(cols) + 1):
            combinations.extend(list(itertools.combinations(cols, i)))

        abbreviations = {"Neuron": "neuro", "Oligodendrocyte": "oligo",
                         "EndothelialCell": "endo", "Microglia": "micro",
                         "Macrophage": "macro", "Astrocyte": "astro",
                         "CellMapNNLS": "NNLS", "CellMapPCA_PC1": "PCA",
                         "CellMapNMF_C1": "NMF", "McKenzie": "McKe"}
        abbr_cols = []
        for col in cols:
            if col in abbreviations.keys():
                abbr_cols.append(abbreviations[col])
            else:
                abbr_cols.append(col)

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

        s = pd.Series(data,
                      index=pd.MultiIndex.from_tuples(indices, names=abbr_cols))
        s.name = "value"
        return s

    @staticmethod
    def plot_pie(total, part, title, outdir, extension):
        labels = ['YES', 'NO']
        sizes = [part, total - part]
        explode = (0.1, 0)

        fig, ax = plt.subplots()
        _, _, autotexts = ax.pie(sizes,
                                 explode=explode,
                                 labels=labels,
                                 colors=["#0072B2", "#808080"],
                                 startangle=90,
                                 autopct=lambda pct: "{:.1f}%\n(n={:d})".format(
                                     pct, int(pct / 100. * sum(sizes))),
                                 textprops=dict(fontweight="bold",
                                                fontsize=20))
        for autotext in autotexts:
            autotext.set_color('white')
        ax.axis('equal')

        fig.suptitle("Cell-Type Mediated eQTLs", fontsize=20, fontweight='bold')
        fig.savefig(os.path.join(outdir, "{}_pie.{}".format(title, extension)))
        plt.show()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL and interaction matrix shape: {}".format(
            self.eqtlinter_df.shape))
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > CellMap Methods: {}".format(self.cellmap_methods))
        print("  > Marker Genes: {}".format(self.marker_genes))
        print("  > Output directory: {}".format(self.outdir))
        print("")
