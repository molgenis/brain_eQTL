#!/usr/bin/env python3

"""
File:         cellmap_genes_per_celltype.py
Created:      2020/04/08
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
import scipy.stats as stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

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
        self.translate_path = "/groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz"
        self.expression_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/output/create_deconvolution_matrices/ct_profile_expr_table.txt.gz"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Load the profile")
        profile_df = pd.read_csv(self.profile_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.profile_path),
                                      profile_df.shape))

        # Load the translate file.
        print("Loading translate matrix.")
        trans_df = pd.read_csv(self.translate_path, sep="\t", header=0, index_col=None)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.translate_path),
                                      trans_df.shape))

        # # Filter the translate dict.
        # trans_df = trans_df.loc[:, ["ArrayAddress", "Symbol"]]
        # trans_df = trans_df.set_index("Symbol")
        # trans_df.index.name = "-"

        # Load the expression file.
        print("Loading expression matrix.")
        expr_df = pd.read_csv(self.expression_path, sep="\t", header=0,
                              index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.expression_path),
                                      expr_df.shape))

        # # Merge the two dataframes.
        # expr_df = expr_df.merge(trans_df, left_index=True, right_index=True)
        # expr_df = expr_df.set_index("Symbol")
        # expr_df.index.name = "-"

        # Subset the expression genes in the cellmap reference profile.
        expr_df = expr_df.loc[expr_df.index.isin(profile_df.index), :]
        self.clustermap(expr_df, "all", self.outdir)

        # Find the genes specific to each celltype.
        gene_celltypes = self.normalize(profile_df).idxmax(axis=1)

        for celltype in profile_df.columns:
            print("\tWorking on: {}".format(celltype))
            ct_genes = gene_celltypes[gene_celltypes == celltype].index
            ct_expr = expr_df.loc[expr_df.index.isin(ct_genes), :]

            self.clustermap(ct_expr, celltype, self.outdir)

    @staticmethod
    def clustermap(df, celltype, outdir):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           yticklabels=False, xticklabels=False,
                           figsize=(12, 9))
        plt.setp(g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_ymajorticklabels(),
            fontsize=10))
        g.fig.suptitle('CellMap Profile {}'.format(celltype))
        plt.tight_layout()
        g.savefig(os.path.join(outdir,
                               "{}_expresion_heatmap.png".format(celltype)))
        plt.close()

    @staticmethod
    def normalize(df):
        df = df.loc[df.std(axis=1) > 0, :]
        df = df.T
        df = (df-df.mean()) / df.std()
        df = df.T
        return df


if __name__ == '__main__':
    m = main()
    m.start()
