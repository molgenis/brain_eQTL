"""
File:         create_cohort_matrix.py
Created:      2020/10/08
Last Changed: 2020/11/20
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


class CreateCohortMatrix:
    def __init__(self, settings, log, sample_dict, sample_order, force, outdir):
        self.inpath = settings["info_datafile"]
        self.sample_id = settings["sample_id"]
        self.cohort_id = settings["cohort_id"]
        self.log = log
        self.sample_dict = sample_dict
        self.sample_order = sample_order
        self.force = force

        # Prepare an output directory.
        outdir = os.path.join(outdir, 'create_cohort_matrix')
        prepare_output_dir(outdir)
        self.outpath = os.path.join(outdir, "cohort_matrix.txt.gz")

        # Declare variables.
        self.sample_info_df = None
        self.cohort_df = None

    def start(self):
        self.log.info("Starting creating cohort matrix.")
        self.print_arguments()

        # Check if output file exist.
        if not check_file_exists(self.outpath) or self.force:
            # Load the sample info.
            self.log.info("Loading sample information matrix.")
            self.sample_info_df = load_dataframe(inpath=self.inpath,
                                                 header=0,
                                                 index_col=None,
                                                 low_memory=False,
                                                 logger=self.log)

            # Construct sample-cohort dict.
            self.log.info("Creating sample to cohort dict.")
            sample_cohort_dict = construct_dict_from_df(self.sample_info_df,
                                                        self.sample_id,
                                                        self.cohort_id)

            # Create cohort dataframe.
            self.log.info("Constructing cohort matrix.")
            self.cohort_df = self.create_cohort_df(self.sample_dict,
                                                   self.sample_order,
                                                   sample_cohort_dict)
            self.save()
        else:
            self.log.info("Skipping step.")

    def create_cohort_df(self, sample_dict, sample_order, sample_cohort_dict):
        cohort_df = pd.DataFrame(0,
                                 index=set(sample_cohort_dict.values()),
                                 columns=sample_order)

        for sample in sample_order:
            cohort = None
            if sample in sample_cohort_dict:
                cohort = sample_cohort_dict[sample]
            else:
                if sample_dict[sample] in sample_cohort_dict:
                    cohort = sample_cohort_dict[sample_dict[sample]]

            if cohort is not None:
                cohort_df.loc[cohort, sample] = 1
            else:
                self.log.warning("Sample {} has no cohort.".format(sample))

        # Validate.
        if not cohort_df.sum(axis=0).all():
            self.log.error("\tSome samples do not have a cohort.")
            exit()

        # Remove cohorts without a sample.
        cohort_df = cohort_df.loc[cohort_df.sum(axis=1) > 0, :]

        return cohort_df

    def save(self):
        save_dataframe(df=self.cohort_df, outpath=self.outpath,
                       index=True, header=True, logger=self.log)

    def clear_variables(self):
        self.inpath = None
        self.sample_dict = None
        self.sample_dict = None
        self.sample_order = None
        self.force = None

    def get_outpath(self):
        return self.outpath

    def get_sample_info_df(self):
        return self.sample_info_df

    def get_cohort_file(self):
        return self.outpath

    def get_cohort_df(self):
        return self.cohort_df

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Input file: {}".format(self.inpath))
        self.log.info("  > Sample ID: {}".format(self.sample_id))
        self.log.info("  > Cohort ID: {}".format(self.cohort_id))
        self.log.info("  > Output path: {}".format(self.outpath))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
