#!/usr/bin/env python3

"""
File:         combine_julienbryois2021_files.py
Created:      2022/01/12
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
__program__ = "Combine Julien Bryois et al. 2021 Files"
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
        self.input_directory = "/Users/mvochteloo/Documents/PhD/Data/JulienBryois2021"
        # self.folders = {"EndothelialCells": "Endothelial.cells"}
        self.folders = {"Astrocytes": "Astrocytes",
                        "EndothelialCells": "Endothelial.cells",
                        "ExcitatoryNeurons": "Excitatory.neurons",
                        "InhibitoryNeurons": "Inhibitory.neurons",
                        "Microglia": "Microglia",
                        "Oligodendrocytes": "Oligodendrocytes",
                        "OPCsCOPs": "OPCs...COPs",
                        "Pericytes": "Pericytes"}
        self.snp_pos_file = "snp_pos.txt"
        self.filter_path = "/Users/mvochteloo/Documents/PhD/Data/MetaBrain_eqtl_filter_list.txt"

        self.outdir = os.path.join(self.input_directory, "COMBINED")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        filter_df = pd.read_csv(self.filter_path, sep="\t", header=None, index_col=None)
        filter_list = set(filter_df.iloc[:, 0].values)

        combined_df_list = []
        for folder, filename in self.folders.items():
            df_list = []
            for i in range(1, 23):
                fpath = os.path.join(self.input_directory, folder, "{}.{}".format(filename, i))
                if not os.path.exists(fpath):
                    print("Filepath '{}' does not exist.".format(fpath))
                    exit()

                df = pd.read_csv(fpath, sep=" ", header=None, index_col=None, nrows=None)
                print("\tLoaded dataframe: {} "
                      "with shape: {}".format(os.path.basename(fpath),
                                              df.shape))
                df_list.append(df)

            df = pd.concat(df_list, axis=0)
            df.columns = ["Gene_id", "SNP_id", "Distance to TSS",
                          "Nominal p-value", "Beta"]

            df["gene"] = [x.split("_")[1] for x in df["Gene_id"]]
            df.index = df["gene"] + "_" + df["SNP_id"]
            df = df.loc[[x for x in df.index if x in filter_list], :]

            df = df.loc[:, ["Nominal p-value", "Beta"]]
            df.columns = ["{} p-value".format(folder), "{} beta".format(folder)]

            combined_df_list.append(df)

        df = pd.concat(combined_df_list, axis=1)

        snp_df = pd.read_csv(os.path.join(self.input_directory, self.snp_pos_file), sep="\t", header=0, index_col=None)
        snp_ae_trans_dict = dict(zip(snp_df["SNP"], snp_df["effect_allele"]))
        df["SNP"] = [x.split("_")[1] for x in df.index]
        df["effect_allele"] = df["SNP"].map(snp_ae_trans_dict)
        print(df)

        outpath = os.path.join(self.outdir, "JulienBryois2021SummaryStats.txt.gz")
        df.to_csv(outpath,
                  sep="\t",
                  index=True,
                  header=True,
                  compression="gzip")
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

if __name__ == '__main__':
    m = main()
    m.start()
