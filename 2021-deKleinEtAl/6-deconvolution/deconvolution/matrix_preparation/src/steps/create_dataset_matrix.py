"""
File:         create_dataset_matrix.py
Created:      2020/10/08
Last Changed: 2021/07/28
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
import pandas as pd

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe, construct_dict_from_df


class CreateDatasetMatrix:
    def __init__(self, settings, log, dts_dict, sample_order, force, outdir):
        self.log = log
        self.dts_dict = dts_dict
        self.sample_order = sample_order
        self.force = force

        # Prepare an output directory.
        outdir = os.path.join(outdir, 'create_dataset_matrix')
        prepare_output_dir(outdir)
        self.outpath = os.path.join(outdir, "dataset_matrix.txt.gz")

        # Declare variables.
        self.sample_info_df = None
        self.dataset_df = None

    def start(self):
        self.log.info("Starting creating dataset matrix.")
        self.print_arguments()

        # Check if output file exist.
        if not check_file_exists(self.outpath) or self.force:
            # Create dataset dataframe.
            self.log.info("Constructing dataset matrix.")
            self.dataset_df = self.create_dataset_df(self.dts_dict,
                                                     self.sample_order)
            self.save()
        else:
            self.log.info("Skipping step.")

    def create_dataset_df(self, dts_dict, sample_order):
        dataset_sample_counts = [(key, len(value)) for key, value in dts_dict.items()]
        dataset_sample_counts.sort(key=lambda x: -x[1])
        datasets = [csc[0] for csc in dataset_sample_counts]

        dataset_df = pd.DataFrame(0, index=sample_order, columns=datasets)
        for dataset in datasets:
            dataset_df.loc[dts_dict[dataset], dataset] = 1
        dataset_df.index.name = "-"

        # Validate.
        if not dataset_df.sum(axis=0).all():
            self.log.error("\tSome samples do not have a cohort.")
            exit()

        return dataset_df

    def save(self):
        save_dataframe(df=self.dataset_df, outpath=self.outpath,
                       index=True, header=True, logger=self.log)

    def clear_variables(self):
        self.dts_dict = None
        self.sample_order = None
        self.force = None

    def get_dataset_file(self):
        return self.outpath

    def get_dataset_df(self):
        return self.dataset_df

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Output path: {}".format(self.outpath))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
