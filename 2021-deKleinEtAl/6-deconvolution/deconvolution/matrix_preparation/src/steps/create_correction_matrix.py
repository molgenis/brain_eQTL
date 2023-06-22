"""
File:         create_correction_matrix.py
Created:      2020/10/15
Last Changed: 2022/01/19
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
from functools import reduce

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, \
    save_dataframe


class CreateCorrectionMatrix:
    def __init__(self, settings, log, dataset_file, dataset_df, sample_dict,
                 sample_order, force, outdir):
        self.cov_file = settings["tech_covariates_datafile"]
        self.mds_file = settings["mds_datafile"]
        self.technical_covariates = settings["technical_covariates"]
        self.log = log
        self.dataset_file = dataset_file
        self.dataset_df = dataset_df
        self.sample_dict = sample_dict
        self.sample_order = sample_order
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_correction_matrix')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.outpath = os.path.join(self.outdir, "correction_table.txt.gz")

        # Create empty variable.
        self.correction_df = None
        self.n_tech_covs = 0
        self.n_mds = 0

    def start(self):
        self.log.info("Starting creating technical covariate file.")
        self.print_arguments()

        # Check if output file exist.
        if not check_file_exists(self.outpath) or self.force:
            self.correction_df = self.create_tech_covs_file()
            self.save()
        else:
            self.log.info("Skipping step.")

    def create_tech_covs_file(self):
        # Load the technical covariates.
        self.log.info("Loading technical covariates matrix.")
        tcov_df = load_dataframe(inpath=self.cov_file,
                                 header=0,
                                 index_col=0,
                                 logger=self.log)

        # Filter on samples and technical covariates.
        self.log.info("Filtering on samples and technical covariates.")
        tcov_df.index = [self.sample_dict[x] if x in self.sample_dict else x for x in tcov_df.index]
        tcov_df = tcov_df.loc[self.sample_order, :].copy()
        save_dataframe(df=tcov_df.T, outpath=os.path.join(self.outdir, "technical_covariates_table.txt.gz"),
                       index=True, header=True, logger=self.log)
        if self.technical_covariates:
            save_dataframe(df=tcov_df.loc[:, self.technical_covariates].T,
                           outpath=os.path.join(self.outdir, "technical_covariates_table_subset.txt.gz"),
                           index=True, header=True, logger=self.log)

        # Load the MDS components.
        self.log.info("Loading MDS matrix.")
        mds_df = load_dataframe(inpath=self.mds_file,
                                header=0,
                                index_col=0,
                                logger=self.log)

        # Filter on samples and technical covariates.
        self.log.info("Filtering on samples and technical covariates.")
        mds_df.index = [self.sample_dict[x] if x in self.sample_dict else x for x in mds_df.index]
        mds_df = mds_df.loc[self.sample_order, :].copy()

        save_dataframe(df=mds_df.T, outpath=os.path.join(self.outdir, "mds_covariates_table.txt.gz"),
                       index=True, header=True, logger=self.log)

        tmp_combined_df = tcov_df.merge(mds_df, left_index=True, right_index=True)
        save_dataframe(df=tmp_combined_df.T, outpath=os.path.join(self.outdir, "technical_and_mds_covariates_table.txt.gz"),
                       index=True, header=True, logger=self.log)

        # Loading cohort matrix.
        self.log.info("Loading dataset matrix.")
        if self.dataset_df is None:
            self.dataset_df = load_dataframe(self.dataset_file,
                                             header=0,
                                             index_col=0,
                                             logger=self.log)

        # merge.
        self.log.info("Merging matrices.")
        correction_df = reduce(lambda left, right: pd.merge(left,
                                                           right,
                                                           left_index=True,
                                                           right_index=True),
                             [tcov_df,
                              mds_df,
                              self.dataset_df])
        correction_df = correction_df.T
        correction_df.index.name = "-"
        self.log.info("\t Correction matrix shape: {}".format(correction_df.shape))

        # Validate sample order.
        if not correction_df.columns.equals(self.sample_order):
            correction_df = correction_df[self.sample_order]

        return correction_df

    def save(self):
        save_dataframe(df=self.correction_df, outpath=self.outpath,
                       index=True, header=True, logger=self.log)

    def clear_variables(self):
        self.dataset_file = None
        self.dataset_df = None
        self.cov_file = None
        self.force = None

    def get_correction_file(self):
        return self.outpath

    def get_correction_df(self):
        return self.correction_df

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Covariates input file: {}".format(self.cov_file))
        self.log.info("  > MDS input file: {}".format(self.mds_file))
        if self.dataset_df is not None:
            self.log.info("  > Cohort input shape: {}".format(self.dataset_df.shape))
        else:
            self.log.info("  > Cohort input file: {}".format(self.dataset_file))
        self.log.info("  > Technical covariates: {}".format(self.technical_covariates))
        self.log.info("  > Output path: {}".format(self.outpath))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
