#!/usr/bin/env python3

"""
File:         compare_marker_genes_with_truth.py
Created:      2020/07/10
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
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Compare Marker Genes With Truth"
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
       self.mg_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_output/create_deconvolution_matrices/marker_genes.txt.gz"
       self.gt_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/AMP-AD/IHC_counts.txt.gz"
       #self.gt_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/AMP-AD/single_cell_counts.txt.gz"
       self.outdir = str(Path(__file__).parent.parent)
       self.palette = {
           "Neuron": "#0072B2",
           "Oligodendrocyte": "#009E73",
           "EndothelialCell": "#CC79A7",
           "Microglia": "#E69F00",
           "Macrophage": "#E69F00",
           "Astrocyte": "#D55E00",
           "Pericytes": "#808080"
       }

    def start(self):
        print("Load the marker-gene expression file.")
        mg_df = pd.read_csv(self.mg_path, sep="\t", header=0, index_col=0).T
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.mg_path),
                                      mg_df.shape))

        print("Load the IHC file.")
        gt_df = pd.read_csv(self.gt_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.gt_path),
                                      gt_df.shape))

        overlap = np.intersect1d(mg_df.index, gt_df.index)
        mg_df = mg_df.loc[overlap, :]
        gt_df = gt_df.loc[overlap, :]

        mg_df.index.name = "-"
        gt_df.index.name = "-"

        if not mg_df.index.equals(gt_df.index):
            print("Invalid order")
            exit()

        print(mg_df.shape)
        print(gt_df.shape)

        corr_data = []
        for mg, mg_epr in mg_df.T.iterrows():
            mg_name = " ".join(mg.split("_")[1:])
            for cell_type, ihc_counts in gt_df.T.iterrows():
                coef, p = stats.spearmanr(mg_epr, ihc_counts)
                corr_data.append([mg_name, cell_type, coef, p])
        corr_df = pd.DataFrame(corr_data, columns=["marker gene", "cell type", "coefficient", "p-value"])
        corr_df.sort_values(by="coefficient", ascending=False, inplace=True)
        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None):
            print(corr_df)


if __name__ == '__main__':
    m = main()
    m.start()
