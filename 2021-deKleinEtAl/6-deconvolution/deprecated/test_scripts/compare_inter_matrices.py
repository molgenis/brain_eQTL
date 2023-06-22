#!/usr/bin/env python3

"""
File:         compare_inter_matrices.py
Created:      2020/04/06
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
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Compare Interaction Matrices"
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
        # Interaction matrix files.
        self.eqtl_ia_filepath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/merge_groups/output/interaction_table.txt.gz"
        self.custom_ia_filepath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/output/interaction_table.txt.gz"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Loading translate dicts.")
        eia_df = pd.read_csv(self.eqtl_ia_filepath, sep="\t", header=0,
                             index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.eqtl_ia_filepath),
                                      eia_df.shape))
        cia_df = pd.read_csv(self.custom_ia_filepath, sep="\t", header=0,
                             index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.custom_ia_filepath),
                                      cia_df.shape))

        # Filter on overlapping rows and columns.
        print("Filtering matrices.")
        row_overlap = []
        for cov in eia_df.index:
            if cov in cia_df.index:
                row_overlap.append(cov)

        eia_col_overlap = []
        cia_col_overlap = []
        for sample in eia_df.columns:
            for sample2 in cia_df.columns:
                if sample.startswith(sample2):
                    eia_col_overlap.append(sample)
                    cia_col_overlap.append(sample2)
                    break

        eia_df = eia_df.loc[row_overlap, eia_col_overlap]
        cia_df = cia_df.loc[row_overlap, cia_col_overlap]

        print("EIA DF: {}\tCIA DF: {}".format(eia_df.shape, cia_df.shape))

        # Plot.
        print("Plotting")
        df = pd.DataFrame({"eQTLInteractionAnalyser": eia_df.melt()["value"],
                           "CustomInteractionAnalyser": cia_df.melt()["value"]})

        coef, p_value = stats.spearmanr(df["eQTLInteractionAnalyser"], df["CustomInteractionAnalyser"])

        fig, ax = plt.subplots()
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        g = sns.regplot(x="eQTLInteractionAnalyser",
                        y="CustomInteractionAnalyser",
                        data=df,
                        scatter_kws={'facecolors': '#000000',
                                     'edgecolor': '#000000',
                                     'alpha': 0.5},
                        line_kws={"color": "#D7191C"},
                        )
        g.set_title('Z-score regression (Spearman r = {:.2f}, '
                    'p = {:.2e})'.format(coef, p_value))
        g.set_ylabel('CustomInteractionAnalyser Z-score',
                     fontsize=8,
                     fontweight='bold')
        g.set_xlabel('eQTLInteractionAnalyser Z-score',
                     fontsize=8,
                     fontweight='bold')
        g.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        g.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        fig.savefig(os.path.join(self.outdir, "inter_matrix_comparison.png"))


if __name__ == '__main__':
    m = main()
    m.start()
