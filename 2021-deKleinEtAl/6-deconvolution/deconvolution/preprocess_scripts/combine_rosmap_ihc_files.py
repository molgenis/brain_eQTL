#!/usr/bin/env python3

"""
File:         combine_rosmap_ihc_files.py
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
import glob
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Combine ROSMAP IHC Files"
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
        self.ihc_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/ROSMAP_IHC"
        self.ids_path = "/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/meta/ROSMAP_IDkey.csv"
        self.gte_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_output/combine_gte_files/GTE_combined.txt.gz"
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'AMP-AD')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.trans_dict = {
            "endo": "EndothelialCell",
            "microglia": "Macrophage",
            "astro": "Astrocyte",
            "neuro": "Neuron",
            "oligo": "Oligodendrocyte"
        }


    def start(self):
        combined = None
        for file in glob.glob(os.path.join(self.ihc_path, "IHC.*")):
            df = pd.read_csv(file, sep="\t", header=0)
            df.index = [os.path.basename(file).split(".")[1]]
            if combined is None:
                combined = df.T
            else:
                combined = combined.merge(df.T, left_index=True, right_index=True)
        combined.index.name = "-"
        combined.index = [int(x) for x in combined.index]
        combined.columns = [self.trans_dict[x] for x in combined.columns]
        print(combined)

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


        combined = combined.divide(combined.sum(axis=1), axis=0)
        combined.index = combined.index.map(wgs_dict).map(gte_dict)
        combined = combined.loc[~combined.index.isnull(), :]
        print(combined)
        print(combined.mean(axis=0))

        combined.to_csv(os.path.join(self.outdir, 'IHC_counts.txt.gz'),
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
