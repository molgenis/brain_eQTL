"""
File:         correct_dataset_effects.py
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
import numpy as np
import pandas as pd
import statsmodels.api as sm

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe


class CorrectDatasetEffects:
    def __init__(self, settings, log, dataset_file, dataset_df, sign_expr_file,
                 sign_expr_df,  force, outdir):
        self.log = log
        self.dataset_file = dataset_file
        self.dataset_df = dataset_df
        self.sign_expr_file = sign_expr_file
        self.sign_expr_df = sign_expr_df
        self.force = force
        self.print_interval = 500

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'correct_dataset_effects')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.sign_expr_dc_outpath = os.path.join(self.outdir, "signature_expr_table_dc.txt.gz")

        # Create empty variable.
        self.expr_dc_df = None
        self.sign_expr_dc_df = None

    def start(self):
        self.log.info("Correcting expression data for dataset effects.")
        self.print_arguments()

        self.log.info("Correcting signature expression data.")
        if not check_file_exists(self.sign_expr_dc_outpath) or self.force:
            if self.dataset_df is None:
                self.dataset_df = load_dataframe(self.dataset_file,
                                                 header=0,
                                                 index_col=0,
                                                 logger=self.log)

            if self.sign_expr_df is None:
                self.sign_expr_df = load_dataframe(self.sign_expr_file,
                                                   header=0,
                                                   index_col=0,
                                                   logger=self.log)

            self.sign_expr_dc_df = self.dataset_correction(self.sign_expr_df, self.dataset_df)
            save_dataframe(df=self.sign_expr_dc_df, outpath=self.sign_expr_dc_outpath,
                           index=True, header=True, logger=self.log)
        else:
            self.log.info("\tSkipping step.")

    def dataset_correction(self, raw_expression, dataset_df):
        expression_df = raw_expression.dropna()

        new_expression_data = []
        for i, (index, expression) in enumerate(expression_df.iterrows()):
            if (i % self.print_interval == 0) or (i == (expression_df.shape[0] - 1)):
                self.log.info("\tProcessing {}/{} "
                              "[{:.2f}%]".format(i,
                                                 (expression_df.shape[0] - 1),
                                                 (100 / (expression_df.shape[0] - 1)) * i))

            if i == 0 and not expression.index.equals(dataset_df.index):
                self.log.error("Expression series does not match dataset matrix.")
                exit()

            ols = sm.OLS(expression, dataset_df)
            try:
                ols_result = ols.fit()
                residuals = ols_result.resid.values
                corrected_expression = expression.mean() + residuals
                new_expression_data.append(corrected_expression.tolist())

            except np.linalg.LinAlgError as e:
                self.log.error("\t\tError: {}".format(e))

        new_expression_df = pd.DataFrame(new_expression_data,
                                         index=expression_df.index,
                                         columns=expression_df.columns)

        return new_expression_df

    def clear_variables(self):
        self.dataset_df = None
        self.dataset_df = None
        self.sign_expr_file = None
        self.sign_expr_df = None
        self.force = None

    def get_expr_dc_df(self):
        return self.expr_dc_df

    def get_sign_expr_dc_file(self):
        return self.sign_expr_dc_outpath

    def get_sign_expr_dc_df(self):
        return self.sign_expr_dc_df

    def print_arguments(self):
        self.log.info("Arguments:")
        if self.dataset_df is not None:
            self.log.info("  > Dataset input shape: {}".format(self.dataset_df.shape))
        else:
            self.log.info("  > Dataset input file: {}".format(self.dataset_file))
        if self.sign_expr_df is not None:
            self.log.info("  > Signature expression input shape: {}".format(self.sign_expr_df.shape))
        else:
            self.log.info("  > Signature expression input file: {}".format(self.sign_expr_file))
        self.log.info("  > Output directory: {}".format(self.outdir))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
