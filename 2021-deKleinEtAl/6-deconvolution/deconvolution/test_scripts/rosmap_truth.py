#!/usr/bin/env python3

"""
File:         rosmap_truth.py
Created:      2020/06/30
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
import gzip
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "ROSMAP Turth"
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
        self.sc_umap_path = "/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/filtered_column_metadata.txt.gz"
        self.ids_path = "/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/meta/ROSMAP_IDkey.csv"
        self.gte_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_output/combine_gte_files/GTE_combined.txt.gz"
        self.ct_dict = {
            "Ex": "excitatory neurons",
            "In": "inhibitory neurons",
            "Ast": "astrocytes",
            "Oli": "oligodendrocytes",
            "Opc": "oligodendrocyte precursor cells",
            "Mic": "microglia",
            "End": "endothelial cells",
            "Per": "pericytes"
        }
        self.ct_ct_dict = {
            "Ex": "Neuron",
            "In": "Neuron",
            "Ast": "Astrocyte",
            "Oli": "Oligodendrocyte",
            "Opc": "Oligodendrocyte",
            "Mic": "Macrophage",
            "End": "EndothelialCell",
            "Per": "Pericytes"
        }

        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'AMP-AD')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        print("Load the sc-UMAP file")
        sc_umap_df = pd.read_csv(self.sc_umap_path, sep="\t", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.sc_umap_path),
                                      sc_umap_df.shape))
        print(sc_umap_df)

        print("Load the IDs file")
        id_df = pd.read_csv(self.ids_path, sep=",", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.ids_path),
                                      id_df.shape))
        print(id_df)
        wgs_dict = self.create_dict(id_df, 'projid', 'wgs_id')

        print("Load the GTE file")
        gte_df = pd.read_csv(self.gte_path, sep="\t")
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.gte_path),
                                      gte_df.shape))
        print(gte_df)
        gte_dict = self.create_dict(gte_df, 0, 1)

        # count the cells.
        composition = []
        cell_types = list(sc_umap_df['broad.cell.type'].unique())
        for sample in sc_umap_df['projid'].unique():
            subset = sc_umap_df.loc[sc_umap_df['projid'] == sample, 'broad.cell.type']
            counts = subset.value_counts()
            comp = [sample]
            for cell in cell_types:
                if cell in counts.index:
                    comp.append(counts[cell])
                else:
                    comp.append(0)

            composition.append(comp)
        comp_df = pd.DataFrame(composition, columns=['-'] + cell_types)
        comp_df.set_index('-', inplace=True)
        print(comp_df)
        print(comp_df.shape)

        # Create the output dataset.
        comp_df.columns = [self.ct_ct_dict[x] for x in comp_df.columns]
        comp_df = comp_df.groupby(level=0, axis=1).sum()
        comp_df = comp_df.divide(comp_df.sum(axis=1), axis=0)
        # comp_df['wgs_id'] = comp_df.index.map(wgs_dict)
        # comp_df['expr_id'] = comp_df['wgs_id'].map(gte_dict)
        # print(len(comp_df.index.unique()))
        # print(len(comp_df['wgs_id'].unique()))
        # print(len(comp_df['expr_id'].unique()))
        # print(comp_df)
        comp_df.index = comp_df.index.map(wgs_dict).map(gte_dict)
        comp_df = comp_df.loc[~comp_df.index.isnull(), :]
        print(comp_df)
        print(comp_df.mean(axis=0))

        # Save.
        comp_df.to_csv(os.path.join(self.outdir, 'ground_truth.txt.gz'),
                       compression="gzip",
                       sep="\t",
                       header=True,
                       index=True)

    @staticmethod
    def create_dict(data, key1, key2):
        str_key1 = key1
        if isinstance(key1, int):
            str_key1 = data.columns[key1]
        elif isinstance(key1, str):
            pass
        else:
            "unexpected input"

        str_key2 = key2
        if isinstance(key2, int):
            str_key2 = data.columns[key2]
        elif isinstance(key2, str):
            pass
        else:
            "unexpected input"

        df = data.loc[:, [str_key1, str_key2]].copy()
        df.dropna(inplace=True)
        return dict(zip(df.loc[:, str_key1], df.loc[:, str_key2]))


if __name__ == '__main__':
    m = main()
    m.start()
