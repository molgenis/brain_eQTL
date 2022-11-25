#!/usr/bin/env python3

"""
File:         merge_mr_and_decon_results.py
Created:      2020/11/26
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

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Merge MR and decon results"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
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
        self.table_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/alltraits_MR_table_v3.xlsx"
        self.sheet_names = ['all_results', 'tophit_results']
        self.table_ea = "EA"
        self.decon_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults_alleles_FDR_betas.txt.gz"
        self.decon_ea = "DeconQTLAlleleAssessed"
        self.decon_drop = ["SNPType", "AlleleAssessed"]

        self.outfile = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/MRFindingsWithCTMediationScores.xlsx"

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading decon-QTL data.")
        decon_df = self.load_file(self.decon_path, index_col=None)
        decon_df["chrom"] = decon_df["SNPName"].str.split(":", expand=True)[0].astype(np.int64)
        decon_df["pos"] = decon_df["SNPName"].str.split(":", expand=True)[1].astype(np.int64)
        decon_df["SNP"] = decon_df["SNPName"].str.split(":", expand=True)[2]
        decon_df.drop(["SNPName"] + self.decon_drop, axis=1, inplace=True)
        decon_df.rename(columns={"ProbeName": "ENSG", "HGNCName": "gene"}, inplace=True)
        decon_df.columns = [x.replace("CellMapNNLS_", "") for x in decon_df.columns]
        beta_columns = [x for x in decon_df.columns if x.endswith("_Beta")]
        print(beta_columns)
        print(decon_df)

        print("")
        print("### Step2 ###")
        print("Loading MR results.")
        sheets_dict = self.load_file(self.table_path, skiprows=1, index_col=None)

        print("")
        print("### Step3 ###")
        print("Merging data.")
        with pd.ExcelWriter(self.outfile) as writer:
            for name, sheet in sheets_dict.items():
                if name in self.sheet_names:
                    print(sheet[["coloc.default.value", "coloc.Jimmy.value"]])

                    print("\tMerging sheet {} with decon-QTL data.".format(name))
                    df = sheet.merge(decon_df,
                                     left_on=["chrom", "pos", "SNP", "ENSG", "gene"],
                                     right_on=["chrom", "pos", "SNP", "ENSG", "gene"],
                                     how="left")
                    print(df[["coloc.default.value", "coloc.Jimmy.value"]])

                    # Flip.
                    df["deconQTLFlip"] = df[self.table_ea] != df[self.decon_ea]
                    print(df[[self.table_ea, self.decon_ea, "deconQTLFlip"] + beta_columns])
                    for beta_col in beta_columns:
                        df[beta_col] = df[beta_col] * df["deconQTLFlip"].map({True: -1, False: 1})
                    print(df[[self.table_ea, self.decon_ea, "deconQTLFlip"] + beta_columns])
                    df.drop([self.decon_ea], axis=1, inplace=True)

                    df.to_excel(writer, sheet_name=name, na_rep="NA")
                    print("\t\tSaving sheet {} with shape {}".format(name, df.shape))
        print("")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  skiprows=0, sheet_name=None):
        if path.endswith(".xlsx"):
            df = pd.read_excel(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows, skiprows=skiprows, sheet_name=sheet_name)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, skiprows=skiprows)
        return df

    def print_arguments(self):
        print("Arguments:")
        print("  > Table path: {}".format(self.table_path))
        print("  > Sheet names: {}".format(self.sheet_names))
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Output file: {}".format(self.outfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
