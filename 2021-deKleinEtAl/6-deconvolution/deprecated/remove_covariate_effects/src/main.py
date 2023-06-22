"""
File:         main.py
Created:      2020/10/14
Last Changed: 2020/10/20
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
import statsmodels.api as sm

# Local application imports.
from utilities import check_file_exists, prepare_output_dir, load_dataframe, construct_dict_from_df


class Main:
    def __init__(self, matrix, covariates, sample_dict):
        self.matrix_inpath = matrix
        self.covariates_inpath = covariates
        self.sample_dict_inpath = sample_dict
        self.print_interval = 500
        self.write_interval = 1000

        # Validate input.
        if not self.validate():
            exit()

        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Prepare an output directory.
        outdir = os.path.join(current_dir, 'output')
        prepare_output_dir(outdir)

        # Define the output file.
        self.outpath = os.path.join(outdir, os.path.basename(matrix))

    def validate(self):
        # Check if input files exist.
        for filepath in [self.matrix_inpath, self.covariates_inpath, self.sample_dict_inpath]:
            if filepath is not None and not check_file_exists(filepath):
                print("File {} does not exist".format(filepath))
                return False

        # Check if correct extension.
        if not self.matrix_inpath.endswith(".txt.gz"):
            print("Matrix input must be in .txt.gz format")
            return False

        return True

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        covariates_df, matrix_df, sample_dict = self.load_data()

        print("### Step2 ###")
        order = self.check_overlap(covariates_df, matrix_df, sample_dict)

        print("### Step3 ###")
        self.work(covariates_df.T, order)

        print("")
        print("Program completed.")

    def load_data(self):
        print("Loading covariates.")
        covariates_df = load_dataframe(self.covariates_inpath,
                                       index_col=0,
                                       header=0)

        print("Loading matrix header.")
        matrix_df = load_dataframe(self.matrix_inpath,
                                   index_col=0,
                                   header=0,
                                   nrows=0)

        print("Loading sample_dict.")
        sample_dict = None
        if self.sample_dict_inpath is not None:
            sample_dict_df = load_dataframe(self.sample_dict_inpath,
                                            index_col=None,
                                            header=0)

            sample_dict = construct_dict_from_df(sample_dict_df,
                                                 sample_dict_df.columns[0],
                                                 sample_dict_df.columns[1])

        return covariates_df, matrix_df, sample_dict

    @staticmethod
    def check_overlap(covariates_df, matrix_df, sample_dict):
        matrix_df_columns = [sample_dict[x] if x in sample_dict else x for x in matrix_df.columns]
        matrix_df.columns = matrix_df_columns

        order = []
        for sample in covariates_df.columns:
            if sample in matrix_df.columns:
                order.append(matrix_df_columns.index(sample) + 1)

        print("\t{}/{} [{:.2f}%] of columns are overlapping".format(len(order), len(matrix_df.columns), (100/len(matrix_df.columns)*len(order))))

        return order

    def work(self, covariates_df, order):
        print("Correcting data.")
        buffer = []
        with gzip.open(self.matrix_inpath, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % self.print_interval == 0):
                    print("\tprocessed {} lines.".format(i))

                if len(buffer) > self.write_interval:
                    self.write_buffer(self.outpath, buffer)
                    buffer = []

                splitted_line = np.array(line.decode().strip('\n').split('\t'))
                index = splitted_line[0]
                data = splitted_line[order]
                if i == 0:
                    # validate the order is correct.
                    if not np.array_equal(data, covariates_df.index.to_numpy()):
                        print("Error! Something went wrong in the sample order.")
                        exit()
                    buffer.append('\t'.join([index] + data.tolist()))
                else:
                    data = self.remove_covariates(np.asfarray(data, float), covariates_df)
                    buffer.append('\t'.join([index] + data.astype(str).tolist()))
        f.close()

        if len(buffer) > 0:
            self.write_buffer(self.outpath, buffer)

    @staticmethod
    def remove_covariates(y, X):
        ols = sm.OLS(y, X)
        try:
            ols_result = ols.fit()
            residuals = ols_result.resid.values
            return y.mean() + residuals

        except np.linalg.LinAlgError as e:
            print("\t\tError: {}".format(e))

        return [np.nan] * len(y)

    @staticmethod
    def write_buffer(filename, buffer):
        with gzip.open(filename, 'wb') as f:
            for line in buffer:
                f.write(line.encode())
        f.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Matrix file: {}".format(self.matrix_inpath))
        print("  > Covariates file: {}".format(self.covariates_inpath))
        print("  > Sample dict file: {}".format(self.sample_dict_inpath))
        print("  > Outpath {}".format(self.outpath))
        print("")
