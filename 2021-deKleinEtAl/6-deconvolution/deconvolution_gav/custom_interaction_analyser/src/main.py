"""
File:         main.py
Created:      2020/10/14
Last Changed: 2020/10/27
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
from datetime import datetime
import itertools
import pickle
import random
import time
import gzip
import os

# Third party imports.
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
from functools import reduce

# Local application imports.
from .storage_container import StorageContainer
from local_settings import LocalSettings
from utilities import check_file_exists, prepare_output_dir, load_dataframe


class Main:
    def __init__(self, input_folder, output_folder, settings_file, skip_rows, n_eqtls,
                 n_samples, verbose):
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir, output_folder)
        prepare_output_dir(self.outdir)

        # Safe settings.
        input_dir = settings.get_setting("input_dir")
        filenames = settings.get_setting("filenames")
        self.input = os.path.join(input_dir, input_folder)
        self.geno_inpath = os.path.join(input_dir, input_folder, filenames["genotype"])
        self.expr_inpath = os.path.join(input_dir, input_folder, filenames["expression"])
        self.tech_covs_inpath = os.path.join(input_dir, input_folder, filenames["technical_covariates"])
        self.covs_inpath = os.path.join(input_dir, input_folder, filenames["covariates"])

        self.correct_snp_tc_inter = settings.get_setting("correct_for_snp_tech_cov_interaction")
        self.pvalues_filename = settings.get_setting("real_pvalues_pickle_filename")
        self.coef_filename = settings.get_setting("coef_pickle_filename")
        self.std_err_filename = settings.get_setting("std_err_pickle_filename")
        self.perm_order_filename = settings.get_setting("permutations_order_pickle_filename")
        self.perm_pvalues_filename = settings.get_setting("permuted_pvalues_pickle_filename")
        self.n_perm = settings.get_setting("n_permutations")
        self.max_end_time = int(time.time()) + settings.get_setting("max_runtime_in_hours") * 60 * 60
        self.panic_time = self.max_end_time - (settings.get_setting("panic_time_in_min") * 60)
        self.skip_rows = skip_rows
        self.n_eqtls = n_eqtls
        self.n_samples = n_samples
        self.verbose = verbose

    def start(self):
        self.print_arguments()
        print("Starting Custom Interaction Analyser "
              "[{}]".format(datetime.now().strftime("%d-%m-%Y, %H:%M:%S")))

        # Start the timer.
        start_time = int(time.time())

        # Get the permutation orders.
        permutation_orders = None
        perm_orders_outfile = os.path.join(self.outdir,
                                           self.perm_order_filename + ".pkl")
        if check_file_exists(perm_orders_outfile):
            print("Loading permutation order")
            permutation_orders = self.load_pickle(perm_orders_outfile)

            # Validate the permutation orders for the given input.
            if len(permutation_orders) != (self.n_perm + 1):
                print("\tinvalid")
                permutation_orders = None

            if permutation_orders is not None:
                if permutation_orders[0] != None:
                    print("\tinvalid")
                    permutation_orders = None
                for order in permutation_orders[1:]:
                    if len(order) != self.n_samples:
                        print("\tinvalid")
                        permutation_orders = None
                        break

            print("\tvalid")

        if permutation_orders is None:
            print("Creating permutation order")
            permutation_orders = self.create_perm_orders()
            self.dump_pickle(permutation_orders, self.outdir,
                             self.perm_order_filename)

        # Start the work.
        print("Start the analysis", flush=True)
        storage = self.work(permutation_orders)

        print("Saving output files", flush=True)
        filename_suffix = "{}_{}".format(self.skip_rows, self.skip_rows + storage.get_n_rows())
        self.dump_pickle(storage.get_pvalues(),
                         self.outdir,
                         self.pvalues_filename,
                         filename_suffix=filename_suffix,
                         subdir=True, unique=True)
        self.dump_pickle(storage.get_perm_pvalues(),
                         self.outdir,
                         self.perm_pvalues_filename,
                         filename_suffix=filename_suffix,
                         subdir=True, unique=True)
        self.dump_pickle(storage.get_coefficients(),
                         self.outdir,
                         self.coef_filename,
                         filename_suffix=filename_suffix,
                         subdir=True, unique=True)
        self.dump_pickle(storage.get_std_errors(),
                         self.outdir,
                         self.std_err_filename,
                         filename_suffix=filename_suffix,
                         subdir=True, unique=True)

        # Print the process time.
        run_time = int(time.time()) - start_time
        run_time_min, run_time_sec = divmod(run_time, 60)
        run_time_hour, run_time_min = divmod(run_time_min, 60)
        print("Finished in  {} hour(s), {} minute(s) and "
              "{} second(s)".format(int(run_time_hour),
                                    int(run_time_min),
                                    int(run_time_sec)))
        print("Received {:.2f} analyses per minute".format((self.n_eqtls * (self.n_perm + 1)) /
                                                           (run_time / 60)))

        # Shutdown the manager.
        print("Shutting down manager [{}]".format(
            datetime.now().strftime("%d-%m-%Y, %H:%M:%S")), flush=True)

    def work(self, permutation_orders):
        # Load the data
        print("Loading data", flush=True)
        tech_covs_df = load_dataframe(self.tech_covs_inpath, header=0, index_col=0)
        covs_df = load_dataframe(self.covs_inpath, header=0, index_col=0)

        geno_df = load_dataframe(self.geno_inpath, header=0, index_col=0,
                                 skiprows=[i for i in
                                           range(1,  self.skip_rows + 1)],
                                 nrows=self.n_eqtls)
        expr_df = load_dataframe(self.expr_inpath, header=0, index_col=0,
                                 skiprows=[i for i in
                                           range(1,  self.skip_rows + 1)],
                                 nrows=self.n_eqtls)

        # Validate the dataframes match up.
        dfs = [tech_covs_df, covs_df, geno_df, expr_df]
        for (a, b) in list(itertools.combinations(dfs, 2)):
            if a is not None and b is not None and \
                    not a.columns.identical(b.columns):
                print("Order of samples are not identical.")
                exit()

        # Replace -1 with NaN in the genotype dataframe. This way we can
        # drop missing values.
        geno_df.replace(-1, np.nan, inplace=True)

        # Initialize the storage object.
        print("Creating storage object")
        storage = StorageContainer(colnames=covs_df.index.to_list())

        # Start working.
        print("Starting interaction analyser", flush=True)
        for row_index, eqtl_index in enumerate([i for i in
                                                range(self.skip_rows,
                                                      self.skip_rows +
                                                      self.n_eqtls)]):
            print("\tProcessing eQTL {}/{} "
                  "[{:.0f}%]".format(row_index + 1,
                                     self.n_eqtls,
                                     (100 / self.n_eqtls) * (row_index + 1)),
                  flush=True)
            start_time = time.time()

            # Get the complete genotype row for the permutation later.
            genotype_all = geno_df.iloc[row_index, :].copy()

            # Get the missing genotype indices.
            indices = np.arange(geno_df.shape[1])
            eqtl_indices = indices[~geno_df.iloc[row_index, :].isnull().values]

            # Subset the row and present samples for this eQTL.
            genotype = geno_df.iloc[row_index, eqtl_indices].copy()
            expression = expr_df.iloc[row_index, eqtl_indices].copy()
            technical_covs = tech_covs_df.iloc[:, eqtl_indices].copy()

            # Initialize variables.
            storage.add_row(eqtl_index, "{}_{}".format(genotype.name, expression.name))

            # Create the base model. Null model are all the technical
            # covariates multiplied with the genotype + the SNP.
            intercept = pd.DataFrame(1, index=genotype.index,
                                     columns=["intercept"])
            base_matrix = reduce(lambda left, right: pd.merge(left,
                                                              right,
                                                              left_index=True,
                                                              right_index=True),
                                 [intercept,
                                  genotype.to_frame(),
                                  technical_covs.T])
            if self.correct_snp_tc_inter:
                tech_inter_matrix = technical_covs.mul(genotype, axis=1)
                tech_inter_matrix.index = ["{}_X_SNP".format(x) for x in
                                           technical_covs.index]
                base_matrix = base_matrix.merge(tech_inter_matrix.T,
                                                left_index=True,
                                                right_index=True)

            # Regress out the base model from the expression values.
            expression_hat = self.remove_covariates(expression, base_matrix)

            # Loop over the covariates.
            for cov_index in range(len(covs_df.index)):
                if storage.has_error():
                    break

                # Get the covariate we are processing.
                covariate = covs_df.iloc[cov_index, eqtl_indices].copy()
                cov_name = covariate.name

                if self.verbose:
                    print("\t\tWorking on '{}'".format(cov_name), flush=True)

                # Create the null model.
                null_matrix = covariate.to_frame()
                n_null = null_matrix.shape[0]
                df_null, rss_null, _, _ = self.create_model(null_matrix,
                                                            expression_hat)

                if self.verbose:
                    print("\t\tn_null: {}\tdf_null: {}\trss_null: {}\t".format(n_null, df_null, rss_null))

                # Loop over each permutation sample order. The first order
                # is the normal order and the remainder are random shuffles.
                for order_id, sample_order in enumerate(permutation_orders):
                    if storage.has_error():
                        break

                    if self.verbose:
                        print("\t\t\tWorking on 'order_{}'".format(order_id),
                              flush=True)

                    # Reorder the covariate based on the sample order.
                    # Make sure the labels are in the same order, just
                    # shuffle the values.
                    covariate_all = covs_df.iloc[cov_index, :].copy()
                    if sample_order is not None:
                        covariate_all_index = covariate_all.index
                        covariate_all = covariate_all.reindex(covariate_all.index[sample_order])
                        covariate_all.index = covariate_all_index

                    # Calculate the interaction effect of the covariate of
                    # interest. Then drop the NA's from the interaction
                    # term.
                    inter_of_interest = covariate_all * genotype_all
                    inter_name = "{}_X_SNP".format(cov_name)
                    inter_of_interest.name = inter_name
                    inter_of_interest = inter_of_interest.iloc[eqtl_indices]

                    del covariate_all

                    # Check if the drop is identical (see above).
                    if not inter_of_interest.index.equals(null_matrix.index):
                        print("\t\t\tError in permutation reordering "
                              "(ID: {})".format(order_id), flush=True)
                        storage.set_error()
                        continue

                    # Create the alternative matrix and add the interaction
                    # term.
                    alt_matrix = null_matrix.copy()
                    alt_matrix = alt_matrix.merge(inter_of_interest.to_frame(),
                                                  left_index=True,
                                                  right_index=True)

                    del inter_of_interest

                    # Create the alternative model.
                    n_alt = alt_matrix.shape[0]
                    df_alt, rss_alt, coefficients_alt, std_errors_alt = self.create_model(alt_matrix,
                                                                                          expression_hat,
                                                                                          cols=[inter_name])

                    del alt_matrix

                    if self.verbose:
                        print("\t\t\tn_alt: {}\tdf_alt: {}\trss_alt: {}\talt_coefficients: {}\talt_std_erros: {}".format(n_alt, df_alt, rss_alt, coefficients_alt, std_errors_alt))

                    # Make sure the n's are identical.
                    if n_null != n_alt:
                        print("\t\t\tError due to unequal n_null and n_alt",
                              flush=True)
                        storage.set_error()
                        continue

                    # Safe the coefficient and std error.
                    storage.add_coefficient(order_id, coefficients_alt[inter_name])
                    storage.add_std_error(order_id, std_errors_alt[inter_name])

                    # Compare the null and alternative model.
                    fvalue = self.calc_f_value(rss_null, rss_alt,
                                               df_null, df_alt, n_null)
                    pvalue = self.get_p_value(fvalue, df_null, df_alt, n_null)

                    if self.verbose:
                        print("\t\t\tfvalue: {}\tpvalue: {}".format(fvalue, pvalue))

                    # Safe the p-values.
                    storage.add_pvalue(order_id, pvalue)

                    del fvalue, pvalue

                    # Check whether we are almost running out of time.
                    if time.time() > self.panic_time:
                        print("\tPanic!!!", flush=True)
                        return storage

            # Safe the results of the eQTL.
            storage.store_row()

            # Print the time.
            print("\t\tfinished in {:.4f} second(s).".format(time.time() - start_time, flush=True))

        return storage

    @staticmethod
    def load_pickle(fpath):
        with open(fpath, "rb") as f:
            content = pickle.load(f)
        f.close()

        return content

    @staticmethod
    def dump_pickle(content, directory, filename, filename_suffix="",
                    subdir=False, unique=False):
        if content is None or (isinstance(content, list) and len(content) == 0):
            return

        full_directory = directory
        if subdir:
            full_directory = os.path.join(directory, filename)

            if not os.path.exists(full_directory):
                os.makedirs(full_directory)

        full_filename = filename
        if filename_suffix != "":
            full_filename = "{}_{}".format(filename, filename_suffix)
            if unique:
                full_filename = "{}_{}_{}".format(filename,
                                                  filename_suffix,
                                                  int(time.time()))

        fpath = os.path.join(full_directory, full_filename + ".pkl")

        with open(fpath, "wb") as f:
            pickle.dump(content, f)
        f.close()

        print_str = os.path.join(os.path.basename(os.path.dirname(fpath)),
                                 os.path.basename(fpath))
        print("\tcreated {}".format(print_str))

    def create_perm_orders(self):
        sample_orders = [None]
        for _ in range(self.n_perm):
            sample_orders.append(self.create_random_shuffle())
        return sample_orders

    def create_random_shuffle(self):
        column_indices = [x for x in range(self.n_samples)]
        random.shuffle(column_indices)
        return column_indices

    @staticmethod
    def write_buffer(filename, buffer):
        with gzip.open(filename, 'wb') as f:
            for line in buffer:
                f.write(line.encode())
        f.close()

    @staticmethod
    def remove_covariates(y, X):
        ols = sm.OLS(y.values, X)
        try:
            ols_result = ols.fit()
            residuals = ols_result.resid.values
            return pd.Series(y.mean() + residuals, index=y.index, name=y.name)

        except np.linalg.LinAlgError as e:
            print("\t\tError: {}".format(e))

        return [np.nan] * len(y)

    @staticmethod
    def create_model(X, y, cols=None):
        # Perform the Ordinary least squares fit.
        ols = sm.OLS(y.values, X)
        try:
            ols_result = ols.fit()
        except np.linalg.LinAlgError as e:
            print("\t\tError: {}".format(e))
            return X.shape[1], np.nan, {col: np.nan for col in cols}, {col: np.nan for col in cols}

        # Extract the required info from the model.
        coefficients = {}
        std_errors = {}
        if cols is not None:
            for col in cols:
                if col in X.columns:
                    coef = ols_result.params[col]
                    std_err = ols_result.bse[col]
                else:
                    coef = np.nan
                    std_err = np.nan
                coefficients[col] = coef
                std_errors[col] = std_err

        return X.shape[1], ols_result.ssr, coefficients, std_errors

    @staticmethod
    def calc_f_value(rss1, rss2, df1, df2, n):
        if (rss1 == np.nan) or (rss2 == np.nan):
            return np.nan
        if df1 >= df2:
            return np.nan
        if df2 >= n:
            return np.nan
        if rss2 >= rss1:
            return 0

        return ((rss1 - rss2) / (df2 - df1)) / (rss2 / (n - df2))

    @staticmethod
    def get_p_value(f_value, df1, df2, n):
        # Lower and upper limit of stats.f.sf
        # stats.f.sf(1827.95, dfn=1, dfd=3661) = 5e-324
        # stats.f.sf(9.9e-12, dfn=1, dfd=3661) = 0.9999974567714613

        # Lower and upper limit of cdf
        # stats.f.cdf(69, dfn=1, dfd=3661) = 0.9999999999999999
        # stats.f.cdf(1e-320, dfn=1, dfd=3661) = 1.0730071046473278e-160

        if f_value == np.nan:
            return np.nan
        if df1 >= df2:
            return np.nan
        if df2 >= n:
            return np.nan

        return stats.f.sf(f_value, dfn=(df2 - df1), dfd=(n - df2))

    def print_arguments(self):
        panic_time_string = datetime.fromtimestamp(self.panic_time).strftime(
            "%d-%m-%Y, %H:%M:%S")
        end_time_string = datetime.fromtimestamp(self.max_end_time).strftime(
            "%d-%m-%Y, %H:%M:%S")
        print("Arguments:")
        print("  > Input directory: {}".format(self.input))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Genotype datafile: {}".format(self.geno_inpath))
        print("  > Expression datafile: {}".format(self.expr_inpath))
        print("  > Technical covariates datafile: {}".format(self.tech_covs_inpath))
        print("  > Covariates datafile: {}".format(self.covs_inpath))
        print("  > Panic datetime: {}".format(panic_time_string))
        print("  > Max end datetime: {}".format(end_time_string))
        print("  > Skip rows: {}".format(self.skip_rows))
        print("  > EQTLs: {}".format(self.n_eqtls))
        print("  > Samples: {}".format(self.n_samples))
        print("  > Permutations: {}".format(self.n_perm))
        print("  > Verbose: {}".format(self.verbose))
        print("", flush=True)
