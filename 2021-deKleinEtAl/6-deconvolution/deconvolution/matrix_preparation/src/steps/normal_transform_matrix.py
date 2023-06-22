"""
File:         normal_transform_matrix.py
Created:      2020/10/27
Last Changed: 2020/10/28
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
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe


class NormalTransformMatrix:
    def __init__(self, settings, log, df, inpath, force, outdir):
        self.log = log
        self.df = df
        self.inpath = inpath
        self.force = force
        self.print_interval = 500

        # Prepare an output directories.
        outdir = os.path.join(outdir, 'normal_transform_matrix')
        prepare_output_dir(outdir)

        # Declare output path.
        self.outpath = os.path.join(outdir, os.path.basename(self.inpath))
        self.normalized_df = None

    def start(self):
        self.log.info("Starting normal transforming matrix.")
        self.print_arguments()

        # Check if output file exist.
        if not check_file_exists(self.outpath) or self.force:
            self.normalized_df = self.normal_transform()
            self.save()
        else:
            self.log.info("Skipping step.")

    def normal_transform(self):
        # loading deconvolution matrix.
        if self.df is None:
            self.log.info("Loading matrix.")
            self.df = load_dataframe(self.inpath,
                                     header=0,
                                     index_col=0,
                                     logger=self.log)

        new_data = []
        match = 0
        count = 0
        print("Processing data.")
        for i, (index, row) in enumerate(self.df.iterrows()):
            if (i == 0) or (i % self.print_interval == 0):
                print("\tprocessed {}/{} [{:.2f}%] lines".format(i, self.df.shape[0], (100 / self.df.shape[0])*i))
                print("\t\t{}/{} are identical".format(match, count))

            zscores1 = self.rank_method1(index, row)
            zscores2 = self.rank_method2(index, row)

            if np.array_equal(zscores1, zscores2):
                match += 1

            count += 1

        exit()
        return pd.DataFrame(new_data, index=self.df.index, columns=self.df.columns)

    @staticmethod
    def rank_method1(index, row):
        work_df = row.to_frame()

        # Get the rank of each value.
        tmp = work_df.copy()
        tmp.sort_values(by=index, ascending=True, inplace=True)
        tmp["rank"] = range(1, tmp.shape[0] + 1)
        tmp = tmp.drop([index], axis=1)

        # Combine the dataframes.
        work_df = work_df.merge(tmp, left_index=True, right_index=True)

        # work_df["pvalue"] = 1 - (work_df["rank"] / (work_df.shape[0] + 1))
        # work_df["zscore"] = stats.norm.isf(work_df["pvalue"])
        work_df["pvalue"] = (work_df["rank"] - 0.5) / work_df.shape[0]
        work_df["zscore"] = stats.norm.ppf(work_df["pvalue"])
        work_df.loc[work_df["pvalue"] > (1.0 - 1e-16), "zscore"] = -8.209536151601387
        work_df.loc[work_df["pvalue"] < 1e-323, "zscore"] = 38.44939448087599

        return work_df["zscore"].values

    @staticmethod
    def rank_method2(index, row):
        work_df = row.to_frame()
        work_df["rank"] = work_df.loc[:, index].rank(ascending=True)
        work_df["pvalue"] = (work_df["rank"] - 0.5) / work_df.shape[0]
        work_df["zscore"] = stats.norm.ppf(work_df["pvalue"])
        work_df.loc[work_df["pvalue"] > (1.0 - 1e-16), "zscore"] = -8.209536151601387
        work_df.loc[work_df["pvalue"] < 1e-323, "zscore"] = 38.44939448087599

        return work_df["zscore"].values

    def save(self):
        print("\tSaving matrix.")
        save_dataframe(df=self.normalized_df, outpath=self.outpath,
                       index=True, header=True, logger=self.log)

    def clear_variables(self):
        self.df = None
        self.inpath = None
        self.force = None

    def print_arguments(self):
        self.log.info("Arguments:")
        if self.df is not None:
            self.log.info("  > Input shape: {}".format(self.df.shape))
        else:
            self.log.info("  > Input file: {}".format(self.inpath))
        self.log.info("  > Output path: {}".format(self.outpath))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
