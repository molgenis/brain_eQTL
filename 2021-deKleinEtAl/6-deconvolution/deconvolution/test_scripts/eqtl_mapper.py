#!/usr/bin/env python3

"""
File:         eqtl_mapper.py
Created:      2022/04/01
Last Changed: 2022/04/06
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
from scipy.special import betainc
from statsmodels.stats import multitest

# Local application imports.

"""
Syntax:
./eqtl_mapper.py -h

./eqtl_mapper.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -al /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/expression_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected
    
./eqtl_mapper.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -al /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-DeconExpressionMatrix
    
### Harm-Jan files ###

./eqtl_mapper.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -al /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-HJFiles/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-DeconExpressionMatrix-HJFiles
    
"""

# Metadata
__program__ = "eQTL Mapper"
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
        self.alle_path = getattr(arguments, 'alleles')
        self.expr_path = getattr(arguments, 'expression')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.genotype_na = getattr(arguments, 'genotype_na')
        self.call_rate = getattr(arguments, 'call_rate')
        self.hw_pval = getattr(arguments, 'hardy_weinberg_pvalue')
        self.maf = getattr(arguments, 'minor_allele_frequency')
        self.nrows = getattr(arguments, 'rows')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        if outdir is None:
            outdir = str(Path(__file__).parent)
        self.outdir = os.path.join(outdir, "eqtl_mapper", outfolder)
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
        parser.add_argument("-al",
                            "--alleles",
                            type=str,
                            required=True,
                            help="The path to the alleles matrix.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            default=None,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-na",
                            "--genotype_na",
                            type=int,
                            required=False,
                            default=-1,
                            help="The genotype value that equals a missing "
                                 "value. Default: -1.")
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
        parser.add_argument("-r",
                            "--rows",
                            type=int,
                            default=None,
                            help="The number of rows to analyze.")
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

        if self.std_path is not None:
            std_df = self.load_file(self.std_path, header=0, index_col=None)

            # Validate that the input data matches.
            self.validate_data(std_df=std_df,
                               geno_df=geno_df)
        else:
            # Create sample-to-dataset file with all the samples having the
            # same dataset.
            std_df = pd.DataFrame({"sample": geno_df.columns, "dataset": "None"})

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
        cr_keep_mask = ~(geno_df == self.genotype_na).all(axis=1).to_numpy(dtype=bool)
        geno_stats_df = pd.DataFrame(np.nan, index=geno_df.index, columns=["N", "NaN", "0", "1", "2", "min GS", "HW pval", "allele1", "allele2", "MA", "MAF"])
        geno_stats_df["N"] = 0
        geno_stats_df["NaN"] = geno_df.shape[1]
        geno_stats_df.loc[cr_keep_mask, :] = self.calculate_genotype_stats(df=geno_df.loc[cr_keep_mask, :])

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

        # Add mask to genotype stats data frame.
        geno_stats_df["mask"] = 0
        geno_stats_df.loc[combined_keep_mask, "mask"] = 1

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
        alle_df = self.load_file(self.alle_path, header=0, index_col=0, nrows=self.nrows)

        # Select eQTL rows that meet requirements.
        geno_df = geno_df.loc[combined_keep_mask, :]
        alle_df = alle_df.loc[combined_keep_mask, :]
        expr_df = expr_df.loc[combined_keep_mask, :]

        print("\tValidating input.")
        self.validate_data(std_df=std_df,
                           geno_df=geno_df,
                           alle_df=alle_df,
                           expr_df=expr_df)
        print("", flush=True)
        del std_df

        ########################################################################

        print("### STEP 3 ###")
        print("Pre-processing data.")
        # Add the allele assed column.
        alle_df["AlleleAssessed"] = alle_df["Alleles"].str.split("/", n=None, expand=True)[1]
        alle_df.drop(["AltAllele"], axis=1, inplace=True)
        alle_df.reset_index(drop=True, inplace=True)

        # Convert to numpy for speed.
        geno_m = geno_df.to_numpy(np.float64)
        expr_m = expr_df.to_numpy(np.float64)

        # Replace missing values with nan
        geno_m[geno_m == self.genotype_na] = np.nan

        # Save properties.
        snps = list(geno_df.index)
        genes = list(expr_df.index)
        del geno_df, expr_df

        # Print info.
        n_eqtls = geno_m.shape[0]
        n_samples = geno_m.shape[1]
        print("Summary stats:")
        print("\tN-eQTLs: {:,}".format(n_eqtls))
        print("\tN-samples: {:,}".format(n_samples))
        print("\tN-datasets: {:,}".format(len(datasets)))
        print("", flush=True)

        ########################################################################

        print("### STEP 4 ###")
        print("Analyzing eQTLs.")

        # Initializing output matrices / arrays.
        eqtl_results_m = np.empty((n_eqtls, 8), dtype=np.float64)

        # Start loop.
        start_time = int(time.time())
        last_print_time = None
        for eqtl_index in range(n_eqtls):
            # Print update for user.
            now_time = int(time.time())
            if n_eqtls > 1 and (last_print_time is None or (now_time - last_print_time) >= self.print_interval or eqtl_index == (n_eqtls - 1)):
                print("\t[{}] {:,}/{:,} eQTLs analysed [{:.2f}%]".format(time.strftime('%H:%M:%S', time.gmtime(now_time - start_time)),
                                                                         eqtl_index,
                                                                         (n_eqtls - 1),
                                                                         (100 / (n_eqtls - 1)) * eqtl_index),
                      flush=True)
                last_print_time = now_time

            # Get the genotype.
            genotype = geno_m[eqtl_index, :]

            # Construct the mask to remove missing values.
            mask = ~np.isnan(genotype)
            n = np.sum(mask)

            # Create the matrix.
            X = np.empty((n, 2), np.float32)
            X[:, 0] = 1
            X[:, 1] = genotype[mask]

            # Get the expression.
            y = expr_m[eqtl_index, mask]

            # if each column is unique.
            if (np.min(np.std(X[:, 1:], axis=0)) == 0) or (
                    np.unique(X, axis=1).shape[1] != 4):
                # Save results.
                eqtl_results_m[eqtl_index, :] = np.array(
                    [n] + [np.nan] * 7)
                continue

            # Compute the rss for just the intercept.
            rss_model1 = self.calc_rss(y=y,
                                       y_hat=np.mean(y))

            # Compute the rss for the main eQTL effect.
            inv_m = self.inverse(X=X)
            betas = self.fit(X=X,
                             y=y,
                             inv_m=inv_m)
            rss_model2 = self.calc_rss(y=y,
                                       y_hat=self.predict(X=X,
                                                          betas=betas))
            std = self.calc_std(rss=rss_model2,
                                n=n,
                                df=2,
                                inv_m=inv_m)

            # Calculate eQTL p-value.
            p_value = self.calc_p_value(
                rss1=rss_model1,
                rss2=rss_model2,
                df1=1,
                df2=2,
                n=n
            )

            # Calculate the t-values.
            t_values = betas / std

            # Save results.
            eqtl_results_m[eqtl_index, :] = np.hstack((
                np.array([n]),
                betas,
                std,
                t_values,
                np.array([p_value]))
            )
        print("", flush=True)

        ########################################################################

        print("### STEP 5 ###")
        print("Saving results.")

        df = pd.DataFrame(eqtl_results_m,
                          columns=["N",
                                   "beta-intercept",
                                   "beta-genotype",
                                   "std-intercept",
                                   "std-genotype",
                                   "tvalue-intercept",
                                   "tvalue-genotype",
                                   "p-value"]
                          )
        df = pd.concat([alle_df, df], axis=1)
        df.insert(0, "ProbeName", genes)
        df.insert(0, "SNPName", snps)
        df["FDR"] = np.nan
        df.loc[~df["p-value"].isnull(), "FDR"] = \
        multitest.multipletests(df.loc[~df["p-value"].isnull(), "p-value"],
                                method='fdr_bh')[1]
        print("\t{:,} eQTLs (p-value <0.05)".format(df.loc[df["p-value"] < 0.05, :].shape[0]))
        print("\t{:,} eQTLs (BH-FDR <0.05)".format(df.loc[df["FDR"] < 0.05, :].shape[0]))

        # Save.
        self.save_file(df=df, outpath=os.path.join(self.outdir, "eQTLResults.txt.gz"))
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
    def validate_data(std_df, geno_df=None, alle_df=None, expr_df=None):
        samples = std_df.iloc[:, 0].values.tolist()
        if geno_df is not None and geno_df.columns.tolist() != samples:
            print("\tThe genotype file header does not match "
                  "the sample-to-dataset link file")
            exit()

        if geno_df is not None and alle_df is not None and geno_df.index.tolist() != alle_df.index.tolist():
            print("\tThe alleles file index does not match "
                  "the genotype file index")
            exit()

        if expr_df is not None and expr_df.columns.tolist() != samples:
            print("\tThe expression file header does not match "
                  "the sample-to-dataset link file")
            exit()

    def calculate_call_rate(self, geno_df, std_df, datasets):
        """
        Calculate the fraction of NaNs per dataset.
        """

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
                                  }, index=df.index)
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
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath), df.shape))

    @staticmethod
    def calc_rss(y, y_hat):
        res = y - y_hat
        res_squared = res * res
        return np.sum(res_squared)

    @staticmethod
    def inverse(X):
        X_square = X.T.dot(X)
        try:
            return np.linalg.inv(X_square)
        except np.linalg.LinAlgError:
            print("Warning: using pseudo-inverse")
            return np.linalg.pinv(X_square)

    def fit(self, X, y, inv_m=None):
        if inv_m is None:
            inv_m = self.inverse(X)
        return inv_m.dot(X.T).dot(y)

    @staticmethod
    def predict(X, betas):
        return np.dot(X, betas)

    @staticmethod
    def calc_std(rss, n, df, inv_m):
        return np.sqrt(rss / (n - df) * np.diag(inv_m))

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
        p_value = betainc(dfd / 2, dfn / 2,
                          1 - ((dfn * f_value) / ((dfn * f_value) + dfd)))
        if p_value == 0:
            p_value = 2.2250738585072014e-308
        return p_value

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Alleles path: {}".format(self.alle_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Genotype NaN: {}".format(self.genotype_na))
        print("  > SNP call rate: >{}".format(self.call_rate))
        print("  > Hardy-Weinberg p-value: >{}".format(self.hw_pval))
        print("  > MAF: >{}".format(self.maf))
        print("  > N rows: {}".format(self.nrows))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
