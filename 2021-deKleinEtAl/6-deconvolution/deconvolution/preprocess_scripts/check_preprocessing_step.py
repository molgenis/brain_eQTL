#!/usr/bin/env python3

"""
File:         check_preprocessing_step.py
Created:      2021/07/08
Last Changed:
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
import argparse
import os

# Third party imports.
import pandas as pd
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Check Preprocessing Step"
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
./check_preprocessing_step.py -d1 /groups/umcg-biogen/tmp01/output/2020-11-10-DeconOptimizer/preprocess_scripts/pre_process_expression_matrix/CortexEUR-cis-AllSoloCorrections/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.MDSAndCohortsRemovedOLS.txt.gz -d2 /groups/umcg-biogen/tmp01/output/2020-11-10-DeconOptimizer/preprocess_scripts/pre_process_expression_matrix/CortexEUR-cis-AllSoloCorrections/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CohortsRemovedOLS.txt.gz -n 1000

./check_preprocessing_step.py -d1 /groups/umcg-biogen/tmp01/output/2020-11-10-DeconOptimizer/preprocess_scripts/pre_process_expression_matrix/WrongMDSCorrection/CortexEUR-cis-NoENA-NoGVEX/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.MDSRemovedOLS.txt.gz -d2 /groups/umcg-biogen/tmp01/output/2020-11-10-DeconOptimizer/preprocess_scripts/pre_process_expression_matrix/CortexEUR-cis-NoENA-NoGVEX/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.MDSRemovedOLS.txt.gz -n 1000


# Log2 transformed:
./check_preprocessing_step.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-16-decon-eQTL/cis/cortex/expression_table/output_newCorrection/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.SampleSelection.ProbesWithZeroVarianceRemoved.txt.gz -n 100

# Covariate removed OLS:
./check_preprocessing_step.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.txt.gz -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-16-decon-eQTL/cis/cortex/expression_table/output_newCorrection/MetaBrain.allCohorts.2020-01-31.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.txt.gz -n 100

# Exp addded:
./check_preprocessing_step.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt.gz -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-16-decon-eQTL/cis/cortex/expression_table/MetaBrain.Cortex.newCovCorrection.allCohorts.2020-01-31.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.exponentAdded.txt -n 100

# Selected
./check_preprocessing_step.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/perform_deconvolution/deconvolution_table.txt -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt
# Average = 0.9842


./check_preprocessing_step.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2020-11-20-decon-QTL/cis/cortex/expression_table/2020-07-16-MetaBrainDeconQtlGenes.TMM.SampSelect.ZeroVarRemov.covRemoved.expAdded.txt -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt
# Average = 0.9978

./check_preprocessing_step.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt
# Average = 1.0000



# Log2 transformed: /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-16-decon-eQTL/cis/cortex/expression_table/output_newCorrection/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.SampleSelection.ProbesWithZeroVarianceRemoved.txt.gz
# Covariate removed OLS: /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-16-decon-eQTL/cis/cortex/expression_table/output_newCorrection/MetaBrain.allCohorts.2020-01-31.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.txt.gz
# Exp addded: /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-16-decon-eQTL/cis/cortex/expression_table/MetaBrain.Cortex.newCovCorrection.allCohorts.2020-01-31.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.exponentAdded.txt
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data1_path = getattr(arguments, 'data1')
        self.data2_path = getattr(arguments, 'data2')
        self.nrows = getattr(arguments, 'nrows')

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
        parser.add_argument("-d1",
                            "--data1",
                            type=str,
                            required=True,
                            help="The path to the first data matrix")
        parser.add_argument("-d2",
                            "--data2",
                            type=str,
                            required=True,
                            help="The path to the second data matrix")
        parser.add_argument("-n",
                            "--nrows",
                            type=int,
                            default=None,
                            help="The number of rows to load. Default: all.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")
        df1 = self.load_file(self.data1_path, header=0, index_col=0, nrows=self.nrows)
        df2 = self.load_file(self.data2_path, header=0, index_col=0, nrows=self.nrows)

        print("Checking overlap.")
        if not df1.index.equals(df2.index):
            row_overlap = set(df1.index).intersection(set(df2.index))
            df1 = df1.loc[row_overlap, :]
            df2 = df2.loc[row_overlap, :]

        if not df1.columns.equals(df2.columns):
            col_overlap = set(df1.columns).intersection(set(df2.columns))
            df1 = df1.loc[:, col_overlap]
            df2 = df2.loc[:, col_overlap]

        print(df1)
        print(df2)

        print("Correlating.")
        coef_data = []
        for i in range(df1.shape[0]):
            coef, _ = stats.pearsonr(df1.iloc[i, :], df2.iloc[i, :])
            coef_data.append(coef)

        s = pd.Series(coef_data, index=df1.index)
        print(s)
        print("")
        print("----------------")
        print("Average = {:.4f}".format(s.abs().mean()))

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def print_arguments(self):
        print("Arguments:")
        print("  > Data 1 path: {}".format(self.data1_path))
        print("  > Data 2 path: {}".format(self.data2_path))
        print("  > N rows: {}".format(self.nrows))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
