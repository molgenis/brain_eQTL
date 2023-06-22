#!/usr/bin/env python3

"""
File:         marker_genes_selection.py
Created:      2020/04/08
Last Changed: 2020/04/26
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
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Marker Gene Selection"
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
        self.excel_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/41598_2018_27293_MOESM3_ESM.xlsx"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Load the excel file")
        df = pd.read_excel(self.excel_path, header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.excel_path),
                                      df.shape))

        # Remove unwanted celltype.
        df = df.loc[df["Celltype"] != "opc", :]

        # Remap.
        df["Celltype"] = df["Celltype"].map({"opc": "oligodendrocyte precursor cell",
                                             "oli": "oligodendrocyte",
                                             "neu": "neuron",
                                             "mic": "microglia",
                                             "end": "endothelial cell",
                                             "ast": "astrocyte"})

        colormap = {
            "neuron": "#b38d84",
            "oligodendrocyte": "#5d9166",
            "endothelial cell": "#f2a7a7",
            "microglia": "#e8c06f",
            "astrocyte": "#9b7bb8"}

        # Order.
        df.sort_values(["Celltype", "Cell_type_mentions"], inplace=True,
                       ascending=False)
        df['x'] = df.groupby('Celltype').cumcount()

        print("Plotting.")
        sns.set(rc={'figure.figsize': (10, 7.5)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        df1 = df.loc[df["x"] <= 5, :].copy()
        df2 = df.loc[df["x"] >= 5, :].copy()
        g = sns.lineplot(x="x", y="Cell_type_mentions", hue="Celltype",
                         palette=colormap, linewidth=2,
                         data=df1, ax=ax)
        sns.lineplot(x="x", y="Cell_type_mentions", hue="Celltype",
                     data=df2, legend=False, linewidth=2,
                     palette={i: "#808080" for i in df2["Celltype"].unique()},
                     ax=ax)
        g.set(yscale="log")

        ax.text(0.5, 1.08, "Marker Gene - Cell Type Mentions",
                fontsize=20, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02, "McKenzie et al. 2018",
                fontsize=18, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        g.set_ylabel("cell type mentions",
                     fontsize=16,
                     fontweight='bold')
        g.set_xlabel("index",
                     fontsize=16,
                     fontweight='bold')

        for i in range(0, 5, 1):
            ax.axhline(10**i, ls='-', color="#000000", alpha=0.15, zorder=-1)
        ax.axvline(5, ls='--', color="#000000", zorder=-1)

        handles = []
        for celltype, color in colormap.items():
            handles.append(mpatches.Patch(color=color, label=celltype))
        plt.legend(handles=handles)

        ax.tick_params(labelsize=12)
        plt.setp(ax.get_legend().get_texts(), fontsize='9')
        plt.setp(ax.get_legend().get_title(), fontsize='10', fontweight='bold')
        fig.savefig(os.path.join(self.outdir, "mckenzie_cell_type_mentions.png"))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
