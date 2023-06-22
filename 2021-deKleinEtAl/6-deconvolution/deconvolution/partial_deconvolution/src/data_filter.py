"""
File:         data_filter.py
Created:      2020/09/04
Last Changed: 2021/08/04
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
import pandas as pd
import numpy as np

# Local application imports.


class DataFilter:
    def __init__(self, settings, raw_expression):
        self.sample_to_dataset_path = settings.get_sample_to_dataset_path()
        self.sample_filter_path = settings.get_sample_filter_path()
        self.dataset_filter = settings.get_dataset_filter()
        self.dataset_correction = settings.get_dataset_correction()
        self.raw_expression = raw_expression

        self.sample_to_dataset_dict = None
        self.filtered_expression = None
        self.datasets = None
        self.reference_dataset = None

    def work(self):
        # Load the sample to dataset file and sample filter file.
        std_df = self.load_sample_to_dataset_file()
        sample_filter = self.load_sample_filter()

        self.set_sample_to_dataset_dict(std_df=std_df)

        # Construct list of requested samples.
        incl_samples = list(std_df.iloc[:, 0].values)
        if self.dataset_filter is not None:
            incl_samples = list(std_df.loc[std_df.iloc[:, 1].isin(self.dataset_filter), std_df.columns[0]].values)
        if sample_filter is not None:
            incl_samples = [x for x in incl_samples if x in sample_filter]

        # Filter the expression.
        self.filtered_expression = self.raw_expression.loc[:, incl_samples]

        # Create the cohorts matrix.
        if self.dataset_correction:
            self.datasets = self.create_datasets_df(std_df=std_df)
            self.reference_dataset = self.datasets.columns[0]

    def load_sample_to_dataset_file(self):
        print("Loading samples annotation")
        return pd.read_csv(self.sample_to_dataset_path,
                           header=0,
                           index_col=None,
                           sep="\t",
                           low_memory=False)

    def load_sample_filter(self):
        print("Loading sample filter")
        if self.sample_filter_path is None:
            return None
        else:
            df = pd.read_csv(self.sample_filter_path, sep="\t",
                             header=0, index_col=0,
                             nrows=1)
            return list(df.columns)

    @staticmethod
    def create_datasets_df(std_df):
        print("Creating datasets data frame.")

        dataset_sample_counts = list(zip(*np.unique(std_df.iloc[:, 1], return_counts=True)))
        dataset_sample_counts.sort(key=lambda x: -x[1])
        datasets = [csc[0] for csc in dataset_sample_counts]

        samples = list(std_df.iloc[:, 0].values)

        dataset_df = pd.DataFrame(0, index=samples, columns=datasets)
        for dataset in datasets:
            dataset_df.loc[std_df.loc[std_df.iloc[:, 1] == dataset, std_df.columns[0]], dataset] = 1
        dataset_df.index.name = "-"

        return dataset_df

    def set_sample_to_dataset_dict(self, std_df):
        self.sample_to_dataset_dict = dict(zip(std_df.iloc[:, 0], std_df.iloc[:, 1]))

    def get_sample_to_dataset_dict(self):
        return self.sample_to_dataset_dict

    def get_filtered_expression(self):
        return self.filtered_expression

    def get_datasets(self):
        return self.datasets

    def get_shape_diff(self):
        return self.raw_expression.shape[0] - self.filtered_expression.shape[0], self.raw_expression.shape[1] - self.filtered_expression.shape[1]

    def get_reference_dataset(self):
        return self.reference_dataset

    def print_info(self):
        print("Pre-filter shape:")
        print("\t{} rows and {} columns".format(self.raw_expression.shape[0], self.raw_expression.shape[1]))
        print("Post-filter shape:")
        print("\t{} rows and {} columns".format(self.filtered_expression.shape[0], self.filtered_expression.shape[1]))
        print("Shape difference:")
        print("\t{} rows and {} columns".format(self.get_shape_diff()[0], self.get_shape_diff()[1]))
        if self.dataset_correction:
            print("Reference dataset = {}".format(self.datasets.columns[0]))
        print("")
