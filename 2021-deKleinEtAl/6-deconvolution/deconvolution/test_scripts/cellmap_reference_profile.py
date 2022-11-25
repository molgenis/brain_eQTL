#!/usr/bin/env python3

"""
File:         cellmap_reference_profile.py
Created:      2020/04/06
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
import scipy.stats as stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "CellMap Profiles"
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
        self.profile_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Load the profile")
        df = pd.read_csv(self.profile_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.profile_path),
                                      df.shape))
        df.columns = [x.split("_")[1] for x in df]

        # Convert the profile expression CPM to z-scores.
        print("Normalize the data")
        normalized_df = self.normalize(df).T
        normalized_df.index.name = None
        print("\tNew shape: {}".format(normalized_df.shape))

        # Color map.
        colormap = {"Neuron": "#b38d84", "Oligodendrocyte": "#5d9166",
                    "EndothelialCell": "#f2a7a7", "Macrophage": "#e8c06f",
                    "Astrocyte": "#9b7bb8"}

        # Assign colors.
        col_colors = normalized_df.idxmax(axis=0).map(colormap)
        col_colors.name = "cell type"

        sns.set(color_codes=True)
        g = sns.clustermap(normalized_df, center=0, cmap="RdBu_r",
                           yticklabels=True, xticklabels=False,
                           col_colors=col_colors,
                           dendrogram_ratio=(.1, .1),
                           figsize=(12, 9))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=10))
        g.fig.subplots_adjust(bottom=0.05, top=0.7)
        plt.tight_layout()
        g.savefig(os.path.join(self.outdir, "cellmap_sc_reference_profile.png"))
        plt.close()

        sns.set(rc={'figure.figsize': (3, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        handles = []
        for celltype, color in colormap.items():
            handles.append(mpatches.Patch(color=color, label=celltype))
        ax.legend(handles=handles)
        ax.set_axis_off()
        fig.savefig(os.path.join(self.outdir, "legend.png"))
        plt.close()

    @staticmethod
    def normalize(df):
        df = df.loc[df.std(axis=1) > 0, :]
        df = df.subtract(df.mean(axis=1), axis=0).divide(df.std(axis=1), axis=0)
        return df


if __name__ == '__main__':
    m = main()
    m.start()
