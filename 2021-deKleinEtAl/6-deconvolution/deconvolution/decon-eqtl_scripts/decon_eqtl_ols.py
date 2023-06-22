#!/usr/bin/env python3

"""
File:         decon_eqtl_ols.py
Created:      2021/09/27
Last Changed: 2022/02/10
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
import argparse
import time
import os

# Third party imports.
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import betainc

# Local application imports.

"""
Syntax:
/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_ols.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/SampleToDataset.txt.gz -of CortexEUR-cis-NormalisedMAF5 -maf 0.05
"""

# Metadata
__program__ = "Decon-eQTL Reborn"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind. The methodology is based on Decon-eQTL: " \
                  "https://doi.org/10.1186/s12859-020-03576-5. Novel additions" \
                  "are the proper handling of missing values as well as " \
                  "calculating a permutation based FDR.".format(__program__,
                                                                __author__,
                                                                __license__)


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.genotype_na = getattr(arguments, 'genotype_na')
        self.nrows = getattr(arguments, 'rows')
        self.call_rate = getattr(arguments, 'call_rate')
        self.hw_pval = getattr(arguments, 'hardy_weinberg_pvalue')
        self.maf = getattr(arguments, 'minor_allele_frequency')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        if outdir is None:
            outdir = str(Path(__file__).parent)
        self.outdir = os.path.join(outdir, "decon_eqtl_ols", outfolder)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set other variables.
        self.alpha = 0.05
        self.print_interval = 30

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
        parser.add_argument("-cc",
                            "--cell_counts",
                            type=str,
                            required=True,
                            help="The path to the cell counts matrix.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-na",
                            "--genotype_na",
                            type=str,
                            required=False,
                            default=-1,
                            help="The genotype value that equals a missing "
                                 "value. Default: -1.")
        parser.add_argument("-r",
                            "--rows",
                            type=int,
                            default=None,
                            help="The number of rows to analyze.")
        parser.add_argument("-cr",
                            "--call_rate",
                            type=float,
                            required=False,
                            default=0.95,
                            help="The minimal call rate of a SNP (per dataset)."
                                 "Equals to (1 - missingness). "
                                 "Default: >= 0.95.")
        parser.add_argument("-hw",
                            "--hardy_weinberg_pvalue",
                            type=float,
                            required=False,
                            default=1e-4,
                            help="The Hardy-Weinberg p-value threshold."
                                 "Default: >= 1e-4.")
        parser.add_argument("-maf",
                            "--minor_allele_frequency",
                            type=float,
                            required=False,
                            default=0.01,
                            help="The MAF threshold. Default: >0.01.")
        parser.add_argument("-od",
                            "--outdir",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of the output path.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### STEP 1 ###")
        print("Loading genotype data and dataset info.")
        geno_df = self.load_file(self.geno_path, header=0, index_col=0, nrows=self.nrows)
        std_df = self.load_file(self.std_path, header=0, index_col=None)

        # Validate that the input data matches.
        self.validate_data(std_df=std_df,
                           geno_df=geno_df)

        print("Checking dataset sample sizes")
        dataset_sample_counts = list(zip(*np.unique(std_df.iloc[:, 1], return_counts=True)))
        dataset_sample_counts.sort(key=lambda x: -x[1])
        datasets = [x[0] for x in dataset_sample_counts]
        max_dataset_length = np.max([len(str(dataset[0])) for dataset in dataset_sample_counts])
        for dataset, sample_size in dataset_sample_counts:
            print("\t{:{}s}  {:,} samples".format(dataset, max_dataset_length, sample_size))
        if dataset_sample_counts[-1][1] <= 1:
            print("\t  One or more datasets have a smaller sample "
                  "size than recommended. Consider excluded these")
            exit()

        print("Calculating genotype call rate per dataset")
        geno_df, call_rate_df = self.calculate_call_rate(geno_df=geno_df,
                                                         std_df=std_df,
                                                         datasets=datasets)
        call_rate_n_skipped = (call_rate_df.min(axis=1) < self.call_rate).sum()
        if call_rate_n_skipped > 0:
            print("\t{:,} eQTLs have had dataset(s) filled with NaN "
                  "values due to call rate threshold ".format(call_rate_n_skipped))

        print("Calculating genotype stats for inclusing criteria")
        geno_stats_df = self.calculate_genotype_stats(df=geno_df)

        # Checking which eQTLs pass the requirements
        n_keep_mask = (geno_stats_df.loc[:, "N"] >= 6).to_numpy()
        hwpval_keep_mask = (geno_stats_df.loc[:, "HW pval"] >= self.hw_pval).to_numpy()
        maf_keep_mask = (geno_stats_df.loc[:, "MAF"] > self.maf).to_numpy()
        combined_keep_mask = n_keep_mask & hwpval_keep_mask & maf_keep_mask
        geno_n_skipped = np.size(combined_keep_mask) - np.sum(combined_keep_mask)
        if geno_n_skipped > 0:
            print("\t{:,} eQTL(s) failed the sample size threshold".format(np.size(n_keep_mask) - np.sum(n_keep_mask)))
            print("\t{:,} eQTL(s) failed the Hardy-Weinberg p-value threshold".format(np.size(hwpval_keep_mask) - np.sum(hwpval_keep_mask)))
            print("\t{:,} eQTL(s) failed the MAF threshold".format(np.size(maf_keep_mask) - np.sum(maf_keep_mask)))
            print("\t----------------------------------------")
            print("\t{:,} eQTL(s) are discarded in total".format(geno_n_skipped))

        print("  Saving results.")
        self.save_file(df=call_rate_df, outpath=os.path.join(self.outdir, "call_rate.txt.gz"))
        self.save_file(df=geno_stats_df, outpath=os.path.join(self.outdir, "geno_stats.txt.gz"))

        del call_rate_df, geno_stats_df

        if geno_n_skipped == self.nrows:
            print("Error, no valid eQTLs.")
            exit()

        print("")

        ########################################################################

        print("### STEP 2 ###")
        print("Loading other data.")
        expr_df = self.load_file(self.expr_path, header=0, index_col=0, nrows=self.nrows)
        cc_df = self.load_file(self.cc_path, header=0, index_col=0)

        # Transpose if need be. We want samples always as columns.
        if cc_df.shape[0] == geno_df.shape[1] or cc_df.shape[0] == expr_df.shape[1]:
            print("\t  Transposing covariate matrix.")
            cc_df = cc_df.T

        # Remove columns with all nan. This line is just because
        # the expression file I used somehow had one column with nan's.
        expr_df.dropna(axis='columns', how='all', inplace=True)

        # Select eQTL rows that meet requirements.
        geno_df = geno_df.loc[combined_keep_mask, :]
        expr_df = expr_df.loc[combined_keep_mask, :]

        print("\tValidating input.")
        self.validate_data(std_df=std_df,
                           geno_df=geno_df,
                           expr_df=expr_df,
                           cc_df=cc_df)
        print("", flush=True)

        del std_df

        ########################################################################

        print("### STEP 3 ###")
        print("Pre-processing data.")
        # Check if the requirements are met.
        if expr_df.min(axis=1).min() < 0:
            print("Error: expression matrix contains negative values.")
            exit()
        if cc_df.min(axis=1).min() < 0:
            print("Error: cell count matrix contains negative values.")
            exit()
        # TODO: check if the column sum is always == 1

        # Convert to numpy for speed.
        geno_m = geno_df.to_numpy(np.float64)
        expr_m = expr_df.to_numpy(np.float64)
        cc_m = cc_df.to_numpy(np.float64)

        # Replace missing values with nan
        geno_m[geno_m == self.genotype_na] = np.nan

        # Revert 2 ** value transform.
        expr_m = np.log2(expr_m)

        # Save properties.
        eqtl_indices = expr_df.index + "_" + geno_df.index
        cell_types_indices = cc_df.index.to_numpy(dtype=object)

        del geno_df, expr_df, cc_df

        # Print info.
        n_eqtls = geno_m.shape[0]
        n_samples = geno_m.shape[1]
        n_covariates = cc_m.shape[0]
        print("\tN-eQTLs: {:,}".format(n_eqtls))
        print("\tN-samples: {:,}".format(n_samples))
        print("\tN-covariates: {:,}".format(n_covariates))
        print("\tN-datasets: {:,}".format(len(datasets)))
        print("", flush=True)

        ########################################################################

        print("### STEP 4 ###")
        print("Analyzing eQTLs.")

        # Save the degrees of freedom the alternative model.
        n_columns = (n_covariates * 2)

        # Initializing output matrices / arrays.
        real_pvalues_m = np.empty((n_eqtls, n_covariates), dtype=np.float64)
        betas_alt_m = np.empty((n_eqtls, n_columns), dtype=np.float64)

        # Construct the X matrix.
        X = np.empty((n_samples, n_columns), np.float32)
        for cov_index in range(n_covariates):
            X[:, cov_index] = cc_m[cov_index, :]

        # Start loop.
        start_time = int(time.time())
        last_print_time = None
        for row_index in range(n_eqtls):
            # Print update for user.
            now_time = int(time.time())
            if n_eqtls > 1 and (last_print_time is None or (now_time - last_print_time) >= self.print_interval or row_index == (n_eqtls - 1)):
                print("\t[{}] {:,}/{:,} eQTLs analysed [{:.2f}%]".format(time.strftime('%H:%M:%S', time.gmtime(now_time - start_time)),
                                                                         row_index,
                                                                         (n_eqtls - 1),
                                                                         (100 / (n_eqtls - 1)) * row_index),
                      flush=True)
                last_print_time = now_time

            # Get the genotype.
            genotype = geno_m[row_index, :]

            # Construct the mask to remove missing values.
            sample_mask = ~np.isnan(genotype)
            n = np.sum(sample_mask)

            # Define the y vector.
            y = expr_m[row_index, sample_mask]

            # Fill in the X matrix.
            for cov_index in range(n_covariates):
                X[:, n_covariates + cov_index] = genotype * X[:, cov_index]

            # df = pd.DataFrame(X, columns=["Intercept", "Genotype"] + [x.split("_")[1] for x in cell_types_indices] + ["{}xGenotype".format(x.split("_")[1]) for x in cell_types_indices])
            # df["Expression"] = expr_m[row_index, :]
            # print(df)
            # df.to_csv(os.path.join(self.outdir, "ieqtl.csv"), sep=",", header=True, index=True)
            #
            # df = df.loc[sample_mask, :]
            # betas_alt, residuals_alt = self.fit_and_predict(X=df, y=y)
            # exit()

            # Model the alternative matrix (with the interaction term).
            # This is the matrix with expression ~ cc1 + cc2 + cc1 * geno +
            # cc2 * geno.
            betas_alt, residuals_alt = self.fit_and_predict(X=X[sample_mask, :], y=y)
            if betas_alt is None or residuals_alt is None:
                betas_alt_m[row_index, :] = np.nan
                real_pvalues_m[row_index, :] = np.nan
                continue
            rss_alt = self.calc_rss(residuals=residuals_alt)

            # Save the alternative model stats.
            betas_alt_m[row_index, :] = betas_alt

            # Remove one interaction column (cc * geno) one by one and
            # determine the significance of the change in residuals sum of
            # squares with a f-test.
            covariates_base_mask = np.ones(n_columns, dtype=bool)
            for cov_index in range(n_covariates):
                # Construct mask.
                covariates_mask = covariates_base_mask.copy()
                covariates_mask[n_covariates + cov_index] = False

                # Model the null matrix (without the interaction term). In
                # this model 1 (and only 1!) of the cc * geno terms is removed.
                _, residuals_null = self.fit_and_predict(X=X[sample_mask, :][:, covariates_mask], y=y)
                if residuals_null is None:
                    real_pvalues_m[row_index, cov_index] = np.nan
                    continue
                rss_null = self.calc_rss(residuals=residuals_null)

                # Calculate and save the p-value.
                p_value = self.calc_p_value(rss1=rss_null,
                                            rss2=rss_alt,
                                            df1=n_columns - 1,
                                            df2=n_columns,
                                            n=n)
                real_pvalues_m[row_index, cov_index] = p_value

        # Cap the p-values.
        real_pvalues_m[real_pvalues_m == 0] = 2.2250738585072014e-308

        # Print the number of significant hits.
        self.print_n_signif(m=real_pvalues_m, colnames=cell_types_indices, type="p-value")
        print("", flush=True)

        ########################################################################

        print("### STEP 5 ###")
        print("Saving results.")

        decon_df = pd.DataFrame(np.hstack((real_pvalues_m, betas_alt_m)),
                                index=eqtl_indices,
                                columns=["{}_pvalue".format(x) for x in cell_types_indices] +
                                        ["Beta{}_{}".format(i + 2, x) for i, x in enumerate(cell_types_indices)] +
                                        ["Beta{}_{}:GT".format(i + 2 + len(cell_types_indices), x) for i, x in enumerate(cell_types_indices)])
        print(decon_df)
        self.save_file(df=decon_df, outpath=os.path.join(self.outdir, "deconvolutionResults.txt.gz"))

        # Save the beta's.
        self.save_matrix(m=betas_alt_m, outpath=os.path.join(self.outdir, "betas_alternative_model.npy"))

        del decon_df, betas_alt_m

        print("", flush=True)

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def validate_data(std_df, geno_df=None, expr_df=None, cc_df=None):
        samples = std_df.iloc[:, 0].values.tolist()
        if geno_df is not None and geno_df.columns.tolist() != samples:
                print("\tThe genotype file header does not match "
                      "the sample-to-dataset link file")
                exit()

        if expr_df is not None and expr_df.columns.tolist() != samples:
                print("\tThe expression file header does not match "
                      "the sample-to-dataset link file")
                exit()

        if cc_df is not None and cc_df.columns.tolist() != samples:
                print("\tThe cell count file header does not match "
                      "the sample-to-dataset link file")
                exit()

    def calculate_call_rate(self, geno_df, std_df, datasets):
        # Calculate the fraction of NaNs per dataset.
        call_rate_df = pd.DataFrame(np.nan, index=geno_df.index, columns=["{} CR".format(dataset) for dataset in datasets])
        for dataset in datasets:
            sample_mask = (std_df.iloc[:, 1] == dataset).to_numpy()
            call_rate_s = (geno_df.loc[:, sample_mask.astype(bool)] != self.genotype_na).astype(int).sum(axis=1) / sample_mask.sum()
            call_rate_df.loc[:, "{} CR".format(dataset)] = call_rate_s

            # If the call rate is too high, replace all genotypes of that
            # dataset with missing.
            row_mask = call_rate_s < self.call_rate
            geno_df.loc[row_mask, sample_mask.astype(bool)] = self.genotype_na

        return geno_df, call_rate_df

    def calculate_genotype_stats(self, df):
        rounded_m = df.to_numpy(dtype=np.float64)
        rounded_m = np.rint(rounded_m)

        # Calculate the total samples that are not NaN.
        nan = np.sum(rounded_m == self.genotype_na, axis=1)
        n = rounded_m.shape[1] - nan

        # Count the genotypes.
        zero_a = np.sum(rounded_m == 0, axis=1)
        one_a = np.sum(rounded_m == 1, axis=1)
        two_a = np.sum(rounded_m == 2, axis=1)

        # Calculate the Hardy-Weinberg p-value.
        hwe_pvalues_a = self.calc_hwe_pvalue(obs_hets=one_a, obs_hom1=zero_a, obs_hom2=two_a)

        # Count the alleles.
        allele1_a = (zero_a * 2) + one_a
        allele2_a = (two_a * 2) + one_a

        # Calculate the MAF.
        maf = np.minimum(allele1_a, allele2_a) / (allele1_a + allele2_a)

        # Determine which allele is the minor allele.
        allele_m = np.column_stack((allele1_a, allele2_a))
        ma = np.argmin(allele_m, axis=1) * 2

        # Construct output data frame.
        output_df = pd.DataFrame({"N": n,
                                  "NaN": nan,
                                  "0": zero_a,
                                  "1": one_a,
                                  "2": two_a,
                                  "HW pval": hwe_pvalues_a,
                                  "allele 1": allele1_a,
                                  "allele 2": allele2_a,
                                  "MA": ma,
                                  "MAF": maf,
                                  })
        del rounded_m, allele_m

        return output_df

    @staticmethod
    def calc_hwe_pvalue(obs_hets, obs_hom1, obs_hom2):
        """
        exact SNP test of Hardy-Weinberg Equilibrium as described in Wigginton,
        JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
        Hardy-Weinberg Equilibrium. AJHG 76: 887-893
        """
        if not 'int' in str(obs_hets.dtype) or not 'int' in str(obs_hets.dtype) or not 'int' in str(obs_hets.dtype):
            obs_hets = np.rint(obs_hets)
            obs_hom1 = np.rint(obs_hom1)
            obs_hom2 = np.rint(obs_hom2)

        # Force homc to be the max and homr to be the min observed genotype.
        obs_homc = np.maximum(obs_hom1, obs_hom2)
        obs_homr = np.minimum(obs_hom1, obs_hom2)

        # Calculate some other stats we need.
        rare_copies = 2 * obs_homr + obs_hets
        l_genotypes = obs_hets + obs_homc + obs_homr
        n = np.size(obs_hets)

        # Get the distribution midpoint.
        mid = np.rint(rare_copies * (2 * l_genotypes - rare_copies) / (2 * l_genotypes)).astype(np.int)
        mid[mid % 2 != rare_copies % 2] += 1

        # Calculate the start points for the evaluation.
        curr_homr = (rare_copies - mid) / 2
        curr_homc = l_genotypes - mid - curr_homr

        # Calculate the left side.
        left_steps = np.floor(mid / 2).astype(int)
        max_left_steps = np.max(left_steps)
        left_het_probs = np.zeros((n, max_left_steps + 1), dtype=np.float64)
        left_het_probs[:, 0] = 1
        for i in np.arange(0, max_left_steps, 1, dtype=np.float64):
            prob = left_het_probs[:, int(i)] * (mid - (i * 2)) * ((mid - (i * 2)) - 1.0) / (4.0 * (curr_homr + i + 1.0) * (curr_homc + i + 1.0))
            prob[mid - (i * 2) <= 0] = 0
            left_het_probs[:, int(i) + 1] = prob

        # Calculate the right side.
        right_steps = np.floor((rare_copies - mid) / 2).astype(int)
        max_right_steps = np.max(right_steps)
        right_het_probs = np.zeros((n, max_right_steps + 1), dtype=np.float64)
        right_het_probs[:, 0] = 1
        for i in np.arange(0, max_right_steps, 1, dtype=np.float64):
            prob = right_het_probs[:, int(i)] * 4.0 * (curr_homr - i) * (curr_homc - i) / (((i * 2) + mid + 2.0) * ((i * 2) + mid + 1.0))
            prob[(i * 2) + mid >= rare_copies] = 0
            right_het_probs[:, int(i) + 1] = prob

        # Combine the sides.
        het_probs = np.hstack((np.flip(left_het_probs, axis=1), right_het_probs[:, 1:]))

        # Normalize.
        sum = np.sum(het_probs, axis=1)
        het_probs = het_probs / sum[:, np.newaxis]

        # Replace values higher then probability of obs_hets with 0.
        threshold_col_a = (max_left_steps - left_steps) + np.floor(obs_hets / 2).astype(int)
        threshold = np.array([het_probs[i, threshold_col] for i, threshold_col in enumerate(threshold_col_a)])
        het_probs[het_probs > threshold[:, np.newaxis]] = 0

        # Calculate the p-values.
        p_hwe = np.sum(het_probs, axis=1)
        p_hwe[p_hwe > 1] = 1

        return p_hwe

    @staticmethod
    def fit_and_predict(X, y):
        model = sm.OLS(y, X)
        try:
            results = model.fit()
            return results.params, results.resid
        except np.linalg.LinAlgError as e:
            print("\t\tError: {}".format(e))
            return None, None

    @staticmethod
    def calc_rss(residuals):
        res_squared = residuals * residuals
        return np.sum(res_squared)

    @staticmethod
    def calc_p_value(rss1, rss2, df1, df2, n):
        """
        last row == stats.f.sf(f_value, dfn=(df2 - df1), dfd=(n - df2)) but this
        is faster somehow.

        https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html

        1 - I(a,b,x) = I(b, a, 1-x)
        """
        if rss2 >= rss1:
            return 1
        dfn = df2 - df1
        dfd = n - df2
        f_value = ((rss1 - rss2) / dfn) / (rss2 / dfd)
        p_value = betainc(dfd / 2, dfn / 2, 1 - ((dfn * f_value) / ((dfn * f_value) + dfd)))
        return p_value

    def print_n_signif(self, m, colnames, type):
        print("\nN-interaction ({} <= {}):".format(type, self.alpha))
        n_hits_a = (m <= self.alpha).sum(axis=0)
        n_hits_total = np.sum(n_hits_a)
        cov_length = np.max([len(x) for x in colnames])
        hits_length = np.max([len(str(x)) for x in n_hits_a] + [len(str(n_hits_total))])
        for n_hits, cell_type in zip(n_hits_a, colnames):
            print("\t{:{}s}  {:{}d}".format(cell_type, cov_length, n_hits, hits_length))
        print("\t{}".format("".join(["-"] * cov_length)))
        print("\t{:{}s}  {:{}d}".format("total", cov_length, n_hits_total, hits_length))

        print("", flush=True)

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath), df.shape))

    @staticmethod
    def save_matrix(m, outpath):
        with open(outpath, 'wb') as f:
            np.save(f, m)
        f.close()

        print("\tSaved matrix: {} "
              "with shape: {}".format(os.path.basename(outpath), m.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > Genotype NaN: {}".format(self.genotype_na))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > N rows: {}".format(self.nrows))
        print("  > SNP call rate: >{}".format(self.call_rate))
        print("  > Hardy-Weinberg p-value: >{}".format(self.hw_pval))
        print("  > MAF: >{}".format(self.maf))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
