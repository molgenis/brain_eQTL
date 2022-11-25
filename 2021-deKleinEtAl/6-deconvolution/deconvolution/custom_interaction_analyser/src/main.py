"""
File:         main.py
Created:      2020/04/23
Last Changed: 2020/06/09
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
from .storage import Storage
from general.local_settings import LocalSettings
from general.utilities import check_file_exists, prepare_output_dir
from general.df_utilities import load_dataframe


class Main:
    """
    Class for the manager.
    """

    def __init__(self, name, settings_file, skip_rows, n_eqtls, n_samples,
                 verbose):
        """
        Initializer of the class.

        :param name: string, the name of the base input/ouput directory.
        :param settings_file: string, the name of the settings file.
        :param skip_rows: int, the number of rows to skip.
        :param n_eqtls: int, the number of eqtls in the input files.
        :param n_samples: int, the number of samples in the input files.
        :param cores: int, the number of cores to use.
        :param verbose: boolean, whether or not to print all update info.
        :param include: boolean, whether or not to include the unfinished
                        wait_list.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir, name)
        prepare_output_dir(self.outdir)

        # Safe settings.
        input_dir = settings.get_setting("input_dir")
        filenames = settings.get_setting("filenames")
        self.geno_inpath = os.path.join(input_dir, name, filenames["genotype"])
        self.expr_inpath = os.path.join(input_dir, name, filenames["expression"])
        self.cov_inpath = os.path.join(input_dir, name, filenames["covariates"])
        self.drop_covs = settings.get_setting("drop_covariates")
        self.tech_covs = settings.get_setting("technical_covariates")
        self.cov_outdir = settings.get_setting("covariates_folder")
        self.tech_cov_outdir = settings.get_setting("technical_covariates_folder")
        self.perm_orders_filename = settings.get_setting("permutations_order_pickle_filename")
        self.pvalues_filename = settings.get_setting("actual_pvalues_pickle_filename")
        self.snp_tvalues_filename = settings.get_setting("snp_tvalues_pickle_filename")
        self.inter_tvalues_filename = settings.get_setting("inter_tvalues_pickle_filename")
        self.perm_pvalues_filename = settings.get_setting("permuted_pvalues_pickle_filename")
        self.n_permutations = settings.get_setting("n_permutations")
        self.max_end_time = int(time.time()) + settings.get_setting("max_runtime_in_hours") * 60 * 60
        self.panic_time = self.max_end_time - (settings.get_setting("panic_time_in_min") * 60)
        self.skip_rows = skip_rows
        self.n_eqtls = n_eqtls
        self.n_samples = n_samples
        self.verbose = verbose

    def start(self):
        """
        Method to start the manager.
        """
        self.print_arguments()
        print("Starting Custom Interaction Analyser "
              "[{}]".format(datetime.now().strftime("%d-%m-%Y, %H:%M:%S")))

        # Start the timer.
        start_time = int(time.time())

        # Get the permutation orders.
        permutation_orders = None
        perm_orders_outfile = os.path.join(self.outdir,
                                           self.perm_orders_filename + ".pkl")
        if check_file_exists(perm_orders_outfile):
            print("Loading permutation order")
            permutation_orders = self.load_pickle(perm_orders_outfile)

            # Validate the permutation orders for the given input.
            if len(permutation_orders) != (self.n_permutations + 1):
                print("\tinvalid")
                permutation_orders = None

            if permutation_orders is not None:
                for order in permutation_orders:
                    if len(order) != self.n_samples:
                        print("\tinvalid")
                        permutation_orders = None
                        break

            print("\tvalid")

        if permutation_orders is None:
            print("Creating permutation order")
            permutation_orders = self.create_perm_orders()
            self.dump_pickle(permutation_orders, self.outdir,
                             self.perm_orders_filename)

        # Start the work.
        print("Start the analyses", flush=True)
        storage = self.work(permutation_orders)
        tc_container = storage.get_tech_cov_container()
        c_container = storage.get_cov_container()

        print("Saving output files", flush=True)
        filename_suffix = "{}_{}".format(self.skip_rows, self.n_eqtls)
        for container, outdir in zip([tc_container, c_container], [self.tech_cov_outdir, self.cov_outdir]):
            full_outdir = os.path.join(self.outdir, outdir)
            prepare_output_dir(full_outdir)

            self.dump_pickle(container.get_pvalues(),
                             full_outdir,
                             self.pvalues_filename,
                             filename_suffix=filename_suffix,
                             subdir=True, unique=True)
            self.dump_pickle(container.get_snp_tvalues(),
                             full_outdir,
                             self.snp_tvalues_filename,
                             filename_suffix=filename_suffix,
                             subdir=True, unique=True)
            self.dump_pickle(container.get_inter_tvalues(),
                             full_outdir,
                             self.inter_tvalues_filename,
                             filename_suffix=filename_suffix,
                             subdir=True, unique=True)
            self.dump_pickle(container.get_perm_pvalues(),
                             full_outdir,
                             self.perm_pvalues_filename,
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
        print("Received {:.2f} analyses per minute".format((self.n_eqtls * (self.n_permutations + 1)) /
                                                           (run_time / 60)))

        # Shutdown the manager.
        print("Shutting down manager [{}]".format(
            datetime.now().strftime("%d-%m-%Y, %H:%M:%S")), flush=True)

    def work(self, permutation_orders):
        """
        Method that does the interaction analysis.

        :param storage: object, a storage object containing all results.
        """
        # Load the data
        print("Loading data", flush=True)
        cov_df = load_dataframe(self.cov_inpath, header=0, index_col=0)

        geno_df = load_dataframe(self.geno_inpath, header=0, index_col=0,
                                 skiprows=[i for i in
                                           range(1,  self.skip_rows + 1)],
                                 nrows=self.n_eqtls)
        expr_df = load_dataframe(self.expr_inpath, header=0, index_col=0,
                                 skiprows=[i for i in
                                           range(1,  self.skip_rows + 1)],
                                 nrows=self.n_eqtls)

        # Drop the covariates we don't want.
        if len(self.drop_covs) > 0:
            cov_df.drop(self.drop_covs, axis=0, inplace=True)

        # Split the covariate table into covariates of interest and technical
        # covariates.
        print("Extracting technical covariates data frame")
        tech_cov_df = cov_df.loc[self.tech_covs, :].copy()
        print("\tShape: {}".format(tech_cov_df.shape))

        # Replace -1 with NaN in the genotype dataframe. This way we can
        # drop missing values.
        geno_df.replace(-1, np.nan, inplace=True)

        # Initialize the storage object.
        print("Creating storage object")
        tech_cov_names = []
        cov_names = []
        for rowname in cov_df.index:
            if rowname in self.tech_covs:
                tech_cov_names.append(rowname)
            else:
                cov_names.append(rowname)
        storage = Storage(tech_covs=tech_cov_names, covs=cov_names)
        storage.print_info()

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

            # Get the complete genotype row for the permutation later.
            genotype_all = geno_df.iloc[row_index, :].copy()

            # Get the missing genotype indices.
            indices = np.arange(geno_df.shape[1])
            eqtl_indices = indices[~geno_df.iloc[row_index, :].isnull().values]

            # Subset the row and present samples for this eQTL.
            genotype = geno_df.iloc[row_index, eqtl_indices].copy()
            expression = expr_df.iloc[row_index, eqtl_indices].copy()
            technical_covs = tech_cov_df.iloc[:, eqtl_indices].copy()
            covariates = cov_df.iloc[:, eqtl_indices].copy()

            # Create the null model. Null model are all the technical
            # covariates multiplied with the genotype + the SNP.
            tech_inter_matrix = technical_covs.mul(genotype, axis=1)
            tech_inter_matrix.index = ["{}_X_SNP".format(x) for x in
                                       technical_covs.index]
            intercept = pd.DataFrame(1, index=genotype.index,
                                     columns=["intercept"])
            base_matrix = reduce(lambda left, right: pd.merge(left,
                                                              right,
                                                              left_index=True,
                                                              right_index=True),
                                 [intercept,
                                  genotype.to_frame(),
                                  technical_covs.T,
                                  tech_inter_matrix.T])

            # Initialize variables.
            storage.add_row(eqtl_index, genotype.name)

            # Loop over the covariates.
            for cov_index in range(len(cov_df.index)):
                if storage.has_error():
                    break

                # Get the covariate we are processing.
                covariate = covariates.iloc[cov_index, :]
                cov_name = covariate.name

                if self.verbose:
                    print("\t\tWorking on '{}'".format(cov_name), flush=True)

                # Add the covariate to the null matrix if it isn't already.
                null_matrix = base_matrix.copy()
                if cov_name not in null_matrix.columns:
                    covariate_df = covariate.copy()
                    null_matrix = null_matrix.merge(covariate_df.to_frame(),
                                                    left_index=True,
                                                    right_index=True)

                # Create the null model.
                n_null = null_matrix.shape[0]
                df_null, rss_null, _ = self.create_model(null_matrix,
                                                         expression)

                # if self.verbose:
                #     print("\t\tn_null: {}\tdf_null: {}\trss_null: {}\t".format(n_null, df_null, rss_null))

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
                    covariate_all = cov_df.iloc[cov_index, :].copy()
                    covariate_all_index = covariate_all.index
                    covariate_all = covariate_all.reindex(
                        covariate_all.index[sample_order])
                    covariate_all.index = covariate_all_index

                    # Calculate the interaction effect of the covariate of
                    # interest. Then drop the NA's from the interaction
                    # term.
                    inter_of_interest = covariate_all * genotype_all
                    inter_name = "{}_X_SNP".format(cov_name)
                    if inter_name in null_matrix.columns:
                        inter_name = inter_name + "_2"
                    inter_of_interest.name = inter_name
                    inter_of_interest = inter_of_interest.iloc[eqtl_indices]

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

                    # Create the alternative model.
                    n_alt = alt_matrix.shape[0]
                    df_alt, rss_alt, alt_tvalues = self.create_model(alt_matrix,
                                                                     expression,
                                                                     tvalue_cols=[genotype.name, inter_name])

                    # if self.verbose:
                    #     print("\t\t\tn_alt: {}\tdf_alt: {}\trss_alt: {}\talt_tvalues: {}".format(n_alt, df_alt, rss_alt, alt_tvalues))

                    # Safe the t-values.
                    storage.add_value(cov_name, order_id, "snp_tvalue", alt_tvalues[genotype.name])
                    storage.add_value(cov_name, order_id, "inter_tvalue", alt_tvalues[inter_name])

                    # Make sure the n's are identical.
                    if n_null != n_alt:
                        print("\t\t\tError due to unequal n_null and n_alt",
                              flush=True)
                        storage.set_error()
                        continue

                    # Compare the null and alternative model.
                    fvalue = self.calc_f_value(rss_null, rss_alt,
                                               df_null, df_alt, n_null)
                    pvalue = self.get_p_value(fvalue, df_null, df_alt, n_null)

                    # if self.verbose:
                    #     print("\t\t\tfvalue: {}\tpvalue: {}".format(fvalue, pvalue))

                    # Safe the p-values.
                    storage.add_value(cov_name, order_id, "pvalue", pvalue)

                    # Check whether we are almost running out of time.
                    if time.time() > self.panic_time:
                        print("\tPanic!!!", flush=True)
                        return storage

            # Safe the results of the eQTL.
            storage.store_row()

        return storage

    @staticmethod
    def load_pickle(fpath):
        """
        Method for loading data from a pickle file.

        :param fpath: string, the pickle input file.
        :return content: , the pickle content.
        """
        with open(fpath, "rb") as f:
            content = pickle.load(f)
        f.close()

        return content

    @staticmethod
    def dump_pickle(content, directory, filename, filename_suffix="",
                    subdir=False, unique=False):
        """
        Method for dumping data to a pickle file.

        :param content: , the pickle content.
        :param directory: string, the pickle output directory.
        :param filename: string, the pickle output file.
        :param filename_suffix: string, a suffix for the output filename.
        :param subdir: boolean, whether or not to put the file in equally
                       named subdirectory.
        :param unique: boolean, whether or not to make the filename unique.
        """
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
        """
        Method for creating x random shuffles of the column indices.

        return sample_orders: list, a nested list with sample orders.
        """
        sample_orders = [self.get_column_indices()]
        for _ in range(self.n_permutations):
            sample_orders.append(self.create_random_shuffle())
        return sample_orders

    def get_column_indices(self):
        """
        Method for getting a ordered list of sample indices.

        :return : list, ordered list with sample indices.
        """
        return [x for x in range(self.n_samples)]

    def create_random_shuffle(self):
        """
        Method for creating a random shuffle of the column indices.

        :return column_indices: list, a shuffles list with sample indices.
        """
        column_indices = self.get_column_indices()
        random.shuffle(column_indices)
        return column_indices

    @staticmethod
    def write_buffer(filename, buffer):
        """
        Method for writing a list of strings to a gzipped file,

        :param filename: string, the output file path.
        :param buffer: list, the lines of strings to write.
        """
        # Write output files.
        if os.path.isfile(filename):
            mode = 'ab'
        else:
            mode = 'wb'

        with gzip.open(filename, mode) as f:
            for line in buffer:
                f.write(line.encode())
        f.close()

    @staticmethod
    def create_model(X, y, tvalue_cols=None):
        """
        Method for creating a multilinear model.

        :param X: DataFrame, the matrix with rows as samples and columns as
                             dimensions.
        :param y: Series, the outcome values.
        :param tvalue_cols: list, the name(s) of the variable(s) to get
                            tvalue(s) for.
        :return df: int, the degrees of freedom of this model.
        :return ssr: float, the residual sum of squares of this fit.
        :return tvalue: dict, beta  / std error.
        """
        # # Filter on features with no standard deviation.
        # stds = X.std(axis=0) > 0
        # stds['intercept'] = True
        # X = X.loc[:, stds]

        df = X.shape[1]

        # Perform the Ordinary least squares fit.
        ols = sm.OLS(y.values, X)
        try:
            ols_result = ols.fit()
        except np.linalg.LinAlgError as e:
            print("\t\tError: {}".format(e))
            return df, np.nan, {col: np.nan for col in tvalue_cols}

        ssr = ols_result.ssr

        tvalues = {}
        if tvalue_cols is not None:
            for col in tvalue_cols:
                tvalue = 0
                if col in X.columns:
                    coef = ols_result.params[col]
                    std_err = ols_result.bse[col]
                    if std_err > 0:
                        tvalue = coef / std_err
                else:
                    tvalue = np.nan
                tvalues[col] = tvalue

        return df, ssr, tvalues

    @staticmethod
    def calc_f_value(rss1, rss2, df1, df2, n):
        """
        Method for comparing the risdual sum squared of two models using
        the F statistic.

        :param rss1: float, the residual sum of squares of the null model.
        :param rss2: float, the residual sum of squares of the alternative model.
        :param df1: int, the degrees of freedom of the null model.
        :param df2: int, the degrees of freedom of the alternative model.
        :param n: int, the number of samples in the model.
        :return : float, the p-value of the comparison.
        """
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
        """
        Method for getting the p-value corresponding to a F-distribution.

        :param f_value: float, the f-value.
        :param df1: int, the degrees of freedom of the null model.
        :param df2: int, the degrees of freedom of the alternative model.
        :param n: int, the number of samples in the model.
        :return : float, the p-value corresponding to the f-value.
        """
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
        """
        Method for printing the variables of the class.
        """
        panic_time_string = datetime.fromtimestamp(self.panic_time).strftime(
            "%d-%m-%Y, %H:%M:%S")
        end_time_string = datetime.fromtimestamp(self.max_end_time).strftime(
            "%d-%m-%Y, %H:%M:%S")
        print("Arguments:")
        print("  > Genotype datafile: {}".format(self.geno_inpath))
        print("  > Expression datafile: {}".format(self.expr_inpath))
        print("  > Covariates datafile: {}".format(self.cov_inpath))
        print("  > Drop covariates: {}".format(self.drop_covs))
        print("  > Technical covariates: {}".format(self.tech_covs))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Permutations: {}".format(self.n_permutations))
        print("  > Panic datetime: {}".format(panic_time_string))
        print("  > Max end datetime: {}".format(end_time_string))
        print("  > Skip rows: {}".format(self.skip_rows))
        print("  > EQTLs: {}".format(self.n_eqtls))
        print("  > Samples: {}".format(self.n_samples))
        print("  > Verbose: {}".format(self.verbose))
        print("", flush=True)
