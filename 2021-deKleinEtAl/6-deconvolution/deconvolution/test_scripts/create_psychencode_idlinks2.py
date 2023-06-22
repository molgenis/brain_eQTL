#!/usr/bin/env python3

"""
File:         create_psychencode_idlinks2.py
Created:      2021/10/08
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
import numpy as np

# Local application imports.

# Metadata
__program__ = "Create PsychENCODE ID Links 2"
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


"""
Syntax:
./create_psychencode_idlinks2.py
"""


class main():
    def __init__(self):
        self.basedir = os.path.join(str(Path(__file__).parent.parent), 'PsychENCODE')
        self.link_file = "PsychENCODE_ID_links.txt.gz"
        self.samplelinks_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-25-samplelinks/"
        self.cortex_eur_std_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/OLD/ContainsDuplicateSamples/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz"

    def start(self):
        print("Step 1: loading data")
        link_df = self.load_file(path=os.path.join(self.basedir, self.link_file), header=0, index_col=None)
        print(link_df)

        gte_dfs = []
        for file_path in glob.glob(os.path.join(self.samplelinks_path, "*.txt")):
            tissue, _, etnicity, dataset, _ = [x.replace("-", "") for x in os.path.basename(file_path).split(".txt")]
            df = self.load_file(path=file_path, header=None, index_col=None)
            df.columns = ["genotype_id", "rnaseq_id"]

            df["tissue"] = tissue
            df["etnicity"] = etnicity
            df["dataset"] = dataset
            gte_dfs.append(df)

        gte_df = pd.concat(gte_dfs, axis=0, ignore_index=True)
        print(gte_df)
        self.save_file(df=gte_df, outpath=os.path.join(self.basedir, "gte_combined_df.txt.gz"), index=False)
        gte_dict = dict(zip(gte_df["genotype_id"].astype(str), gte_df["rnaseq_id"].astype(str)))
        etg_dict = dict(zip(gte_df["rnaseq_id"].astype(str), gte_df["genotype_id"].astype(str)))
        metabrain_gte_matches = set(gte_df["genotype_id"].astype(str) + "+" + gte_df["rnaseq_id"].astype(str))

        matches = []
        for _, row in link_df.iterrows():

            row_values = {}
            for name in ["PsychENCODE ID", "rnaseq_id", "genotype_id", "sample_full"]:
                value = row[name]
                if value in row_values:
                    names = row_values[value]
                    names.append(name)
                    row_values[value] = names
                else:
                    row_values[value] = [name]

            for value, names in row_values.items():
                value = str(value)
                if value in gte_dict.keys():
                    if value + "+" + gte_dict[value] not in metabrain_gte_matches:
                        continue
                    matches.append([row["PsychENCODE ID"], "/".join(names), "genotype_id", gte_dict[value], value])
                if value in etg_dict.keys():
                    if etg_dict[value] + "+" + value not in metabrain_gte_matches:
                        continue
                    matches.append([row["PsychENCODE ID"], "/".join(names), "rnaseq_id", value, etg_dict[value]])

        match_df = pd.DataFrame(matches, columns=["psychencode_id", "match_name", "match_type", "rnaseq_id", "genotype_id"])
        match_df["key"] = match_df["psychencode_id"].astype(str) + "+" + match_df["genotype_id"].astype(str) + "+" + match_df["rnaseq_id"].astype(str)
        match_df = match_df.groupby("key").first()
        match_df = match_df.merge(gte_df, on=["rnaseq_id", "genotype_id"], how="left")
        self.save_file(df=match_df, outpath=os.path.join(self.basedir, "PsychENCODE_ID_links2.txt.gz"), index=False)
        self.save_file(df=match_df[["psychencode_id", "rnaseq_id"]], outpath=os.path.join(self.basedir, "PsychENCODE-linkFile.txt.gz"), index=False)

        for column in ["match_name", "match_type", "tissue", "etnicity", "dataset"]:
            print(column)
            print(match_df[column].value_counts())

        found_samples = set(match_df["psychencode_id"])
        missing_samples = []
        for sample in link_df["PsychENCODE ID"].unique():
            if sample not in found_samples:
                missing_samples.append(sample)
        missing_df = link_df.loc[link_df["PsychENCODE ID"].isin(missing_samples), :]
        print(missing_df)
        self.save_file(df=link_df, outpath=os.path.join(self.basedir, "missing_PsychENCODE_IDs.txt.gz"), index=False)

        cortex_eur_df = match_df.loc[(match_df["tissue"] == "cortex") & (match_df["etnicity"] == "EUR"), :]
        cortex_eur_df_found_samples = set(cortex_eur_df["rnaseq_id"])
        print(cortex_eur_df)
        print(cortex_eur_df["dataset"].value_counts())

        cortex_eur_std_df = self.load_file(path=os.path.join(self.basedir, self.cortex_eur_std_path), header=0, index_col=None)
        cortex_eur_not_found = []
        for index, (sample, dataset) in cortex_eur_std_df.iterrows():
            if sample not in cortex_eur_df_found_samples:
                cortex_eur_not_found.append([sample, dataset])

        cortex_eur_not_found_df = pd.DataFrame(cortex_eur_not_found, columns=["sample", "dataset"])
        print(cortex_eur_not_found_df["dataset"].value_counts())
        self.save_file(df=cortex_eur_not_found_df, outpath=os.path.join(self.basedir, "cortex_eur_not_found.txt.gz"), index=False)
        exit()


    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  skiprows=0, sheet_name=None, low_memory=True):
        if path.endswith(".xlsx"):
            df = pd.read_excel(path, header=header, index_col=index_col,
                         nrows=nrows, skiprows=skiprows, sheet_name=sheet_name)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, skiprows=skiprows,
                             low_memory=low_memory)

        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))


if __name__ == '__main__':
    m = main()
    m.start()
