#!/usr/bin/env python3

"""
File:         correlate_align_metrics_with_genes.py
Created:      2021/12/07
Last Changed: 2022/03/31
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
from statsmodels.regression.linear_model import OLS

# Local application imports.

# Metadata
__program__ = "Correlate Alignment Metrics with Genes"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
Syntax: 
./correlate_align_metrics_with_genes.py -h

./correlate_align_metrics_with_genes.py -am /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis/create_correction_matrix/technical_covariates_table.txt.gz -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-EUR/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz -o 2021-12-07-CortexEUR-cis-round1

./correlate_align_metrics_with_genes.py -am /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/create_correction_matrix/technical_covariates_table.txt.gz -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-AFR/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz -o 2021-12-07-CortexAFR-cis-round3

./correlate_align_metrics_with_genes.py \
    -am /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-02-05-step0-correlate-covariates-with-expression/2020-02-05-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.txt.gz \
    -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-EURandAFR-noENA-noAMPAD/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz \
    -o 2022-03-31-CortexEURandAFR-noENA-noAMPAD-trans-round2
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.am_path = getattr(arguments, 'alignment_metrics')
        self.expr_path = getattr(arguments, 'expression')
        self.out_filename = getattr(arguments, 'outfile')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'correlate_align_metrics_with_genes')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
        parser.add_argument("-am",
                            "--alignment_metrics",
                            type=str,
                            required=True,
                            help="The path to RNAseq alignment metrics "
                                 "matrix.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            default="output",
                            help="The name of the outfile. Default: output.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # Load data.
        print("Loading data.")
        am_df = self.load_file(self.am_path, header=0, index_col=0)
        expr_df = self.load_file(self.expr_path, header=0, index_col=0, nrows=None)

        print("Pre-processing data.")
        # Make sure order is the same.
        samples = set(am_df.columns.tolist()).intersection(set(expr_df.columns.tolist()))
        am_df = am_df.loc[:, samples]
        expr_df = expr_df.loc[:, samples]

        # Filter on columns without std.
        am_df = am_df.loc[am_df.std(axis=1) > 0, :]
        expr_df = expr_df.loc[expr_df.std(axis=1) > 0, :]

        # Drop rows with NaN.
        am_df.dropna(inplace=True)

        # Remove columns that are a linear combination of others.
        am_df = self.remove_multicollinearity(am_df.T).T

        ### EUR cis ###
        # correction_df = am_df.loc[["PCT_MRNA_BASES",
        #                            "PCT_INTRONIC_BASES",
        #                            "PCT_USABLE_BASES",
        #                            "PCT_CODING_BASES",
        #                            "PCT_INTERGENIC_BASES",
        #                            "PCT_UTR_BASES",
        #                            "MEDIAN_3PRIME_BIAS",
        #                            "PCT_READS_ALIGNED_IN_PAIRS",
        #                            "X.GC_R1",
        #                            "PCT_RIBOSOMAL_BASES"
        #                            ], :].T
        # correction_df = am_df.loc[["PCT_MRNA_BASES",
        #                            "PCT_INTRONIC_BASES",
        #                            "PCT_USABLE_BASES",
        #                            "PCT_CODING_BASES",
        #                            "PCT_INTERGENIC_BASES",
        #                            "PCT_UTR_BASES",
        #                            "MEDIAN_3PRIME_BIAS",
        #                            "PCT_READS_ALIGNED_IN_PAIRS",
        #                            "X.GC_R1",
        #                            "PCT_RIBOSOMAL_BASES",
        #                            "deletion_length",
        #                            "MEAN_READ_LENGTH",
        #                            "avg_sequence_length_R1",
        #                            "MEDIAN_5PRIME_BIAS",
        #                            "PCT_PF_READS_IMPROPER_PAIRS",
        #                            "PF_READS_IMPROPER_PAIRS",
        #                            "insertion_length",
        #                            "PF_HQ_ALIGNED_Q20_BASES",
        #                            "PF_HQ_ERROR_RATE",
        #                            "PF_MISMATCH_RATE"
        #                            ], :].T

        ### AFR  cis ###
        # correction_df = am_df.loc[["PCT_INTRONIC_BASES",
        #                            "PCT_MRNA_BASES",
        #                            "PCT_INTERGENIC_BASES",
        #                            "PCT_USABLE_BASES",
        #                            "PCT_CODING_BASES",
        #                            "MEDIAN_3PRIME_BIAS",
        #                            "PCT_UTR_BASES",
        #                            "X.GC_R1",
        #                            "X.GC_R2",
        #                            "deletion_length"
        #                            ], :].T
        # correction_df = am_df.loc[["PCT_INTRONIC_BASES",
        #                            "PCT_MRNA_BASES",
        #                            "PCT_INTERGENIC_BASES",
        #                            "PCT_USABLE_BASES",
        #                            "PCT_CODING_BASES",
        #                            "MEDIAN_3PRIME_BIAS",
        #                            "PCT_UTR_BASES",
        #                            "X.GC_R1",
        #                            "X.GC_R2",
        #                            "deletion_length",
        #                            "MEDIAN_5PRIME_BIAS",
        #                            "PF_MISMATCH_RATE",
        #                            "PF_HQ_ERROR_RATE",
        #                            "PCT_CHIMERAS",
        #                            "PCT_R2_TRANSCRIPT_STRAND_READS",
        #                            "PCT_R1_TRANSCRIPT_STRAND_READS",
        #                            "avg_sequence_length_R1",
        #                            "PCT_PF_READS_IMPROPER_PAIRS",
        #                            "PCT_READS_ALIGNED_IN_PAIRS",
        #                            "MEAN_READ_LENGTH"
        #                            ], :].T

        # correction_df = am_df.loc[["PCT_RIBOSOMAL_BASES",
        #                            "PCT_MRNA_BASES",
        #                            "PCT_INTRONIC_BASES",
        #                            "PCT_USABLE_BASES",
        #                            "PCT_CODING_BASES",
        #                            "INTRONIC_BASES",
        #                            "PCT_UTR_BASES",
        #                            "RIBOSOMAL_BASES",
        #                            "MEDIAN_3PRIME_BIAS",
        #                            "unmapped_tooshort"
        #                            ], :].T
        #
        # print("Correcting gene expression data.")
        # expr_df = self.calculate_residuals(df=expr_df, correction_df=correction_df)

        # Safe the indices.
        metrics = am_df.index.tolist()
        genes = expr_df.index.tolist()

        # Convert to numpy.
        am_m = am_df.to_numpy()
        expr_m = expr_df.to_numpy()
        del am_df, expr_df

        print(am_m.shape)
        print(expr_m.shape)

        print("Correlating.")
        corr_m = np.corrcoef(am_m, expr_m)[:am_m.shape[0], am_m.shape[0]:]
        corr_df = pd.DataFrame(corr_m, index=metrics, columns=genes)
        print(corr_df)

        # print("Saving file.")
        # self.save_file(df=corr_df,
        #                outpath=os.path.join(self.outdir, "{}_correlations.txt.gz".format(self.out_filename)))

        print("Ranking top metrics.")
        summary_df = pd.DataFrame({"max abs correlation": corr_df.abs().max(axis=1),
                                   "mean abs correlation": corr_df.abs().mean(axis=1),
                                   "std abs correlation": corr_df.abs().std(axis=1)})
        summary_df.sort_values(by="max abs correlation", inplace=True, ascending=False)
        print(summary_df)

        print("Saving file.")
        self.save_file(df=summary_df,
                       outpath=os.path.join(self.outdir, "{}_summary.txt.gz".format(self.out_filename)))

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def remove_multicollinearity(self, df, threshold=1):
        indices = np.arange(df.shape[1])
        max_r2 = np.inf
        while len(indices) > 1 and max_r2 >= threshold:
            r2 = np.array([self.calc_ols_rsquared(df.iloc[:, indices], ix) for ix in range(len(indices))])
            max_r2 = max(r2)

            if max_r2 >= threshold:
                max_index = np.where(r2 == max_r2)[0][0]
                print("Removing {}".format(df.columns[max_index]))
                indices = np.delete(indices, max_index)

        return df.iloc[:, indices]

    @staticmethod
    def calc_ols_rsquared(df, idx):
        return OLS(df.iloc[:, idx], df.loc[:, np.arange(df.shape[1]) != idx]).fit().rsquared

    @staticmethod
    def calculate_residuals(df, correction_df):
        corrected_m = np.empty(df.shape, dtype=np.float64)
        last_print_time = None
        n_tests = df.shape[0]
        for i in range(n_tests):
            now_time = int(time.time())
            if last_print_time is None or (now_time - last_print_time) >= 10 or (i + 1) == n_tests:
                last_print_time = now_time
                print("\t{}/{} genes corrected [{:.2f}%]".format((i + 1), n_tests, (100 / n_tests) * (i + 1)))

            ols = OLS(df.iloc[i, :], correction_df)
            results = ols.fit()
            # print(results.summary())
            corrected_m[i, :] = results.resid

        return pd.DataFrame(corrected_m, index=df.index, columns=df.columns)

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Alignment metrics: {}".format(self.am_path))
        print("  > Expression: {}".format(self.expr_path))
        print("  > Output filename: {}".format(self.out_filename))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
