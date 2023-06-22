#!/usr/bin/env python3

"""
File:         select_and_reorder_matrix.py
Created:      2021/06/23
Last Changed: 2021/08/05
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
import gzip
import os

# Third party imports.
import pandas as pd
import numpy as np

# Local application imports.

# Metadata
__program__ = "Select and Reorder Matrix"
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
./select_and_reorder_matrix.py -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt.gz -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -n CortexEUR-cis

./select_and_reorder_matrix.py -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-HalfNormalised/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceHalfNormalised.ExpAdded.txt.gz -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -n CortexEUR-cis-HalfNormalised

./select_and_reorder_matrix.py -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt.gz -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -n CortexEUR-cis-Normalised

./select_and_reorder_matrix.py -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis-PrimaryeQTLs/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt.gz -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis-PrimaryeQTLs/create_matrices/genotype_table.txt.gz -n CortexEUR-cis-PrimaryeQTLs

### Cortex EUR ###

./select_and_reorder_matrix.py \
    -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/2021-12-07-CortexEUR-cis-Normalised/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt.gz \
    -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis/create_matrices/genotype_table.txt.gz \
    -n 2021-12-07-CortexEUR-cis-Normalised

### Cortex EUR trans ###

./select_and_reorder_matrix.py \
    -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-Normalised/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt.gz \
    -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -n 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected
    
./select_and_reorder_matrix.py \
    -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-Normalised/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.100PCsRemovedOLS.ForceNormalised.ExpAdded.txt.gz \
    -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -n 2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected
    
./select_and_reorder_matrix.py \
    -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-Normalised/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt.gz \
    -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -n 2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected
    
./select_and_reorder_matrix.py \
    -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-Normalised/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.80PCsRemovedOLS.ForceNormalised.ExpAdded.txt.gz \
    -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -n 2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected
    
### Cortex EUR trans Harm-Jan files ###

./select_and_reorder_matrix.py \
    -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-Normalised-HJFiles/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt.gz \
    -so /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -n 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected-HJFiles
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.expr_path = getattr(arguments, 'expression')
        self.sample_order_path = getattr(arguments, 'sample_order')
        self.name = getattr(arguments, 'name')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'select_and_reorder_matrix', self.name)
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
        parser.add_argument("-eq",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the replicating eQTLs")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
        parser.add_argument("-so",
                            "--sample_order",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to a matrix of which we take the"
                                 "sample order. Default: None.")
        parser.add_argument("-n",
                            "--name",
                            type=str,
                            required=True,
                            help="The name of the output directory")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading eQTL data.")
        eqtl_df = self.load_file(self.eqtl_path, header=0, index_col=None)
        genes_of_interest = eqtl_df["ProbeName"]

        print("Processing expression file.")
        expr_df = self.process_file(inpath=self.expr_path,
                                    filter_set=set(genes_of_interest))

        print("Reordering matrix.")
        missing = [gene for gene in genes_of_interest if gene not in expr_df.index]
        if len(missing) > 0:
            print("Missing {} genes: {}".format(len(missing), ", ".join(missing)))
            exit()
        expr_df = expr_df.loc[genes_of_interest, :]
        expr_df.columns.name = None
        expr_df.index.name = None

        if self.sample_order_path is not None:
            print("Loading sample order data.")
            sample_order_df = self.load_file(self.sample_order_path, header=0, index_col=0, nrows=1)
            sample_order = sample_order_df.columns.tolist()

            if len(set(sample_order).symmetric_difference(set(expr_df.columns.tolist()))) != 0:
                print("Sample order is not valid. Skipped step.")
            else:
                expr_df = expr_df.loc[:, sample_order]

        print("Saving output file.")
        self.save_file(df=expr_df,
                       outpath=os.path.join(self.outdir, os.path.basename(self.expr_path).replace(".gz", "")))

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
    def process_file(inpath, filter_set, flip_dict=None, sep="\t",
                     header=0, index_col=0):
        found_indices = set()
        info = []
        with gzip.open(inpath, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % 1000 == 0):
                    print("\t{} lines processed, {}/{} rows found [{:.2f}%]".format(i, len(info), len(filter_set), (100 / len(filter_set)) * len(info)))

                splitted_line = line.decode().strip('\n').split(sep)
                if i == 0:
                    info.append(splitted_line)
                else:
                    index = splitted_line[0]
                    if index not in filter_set or index in found_indices:
                        continue

                    if flip_dict is not None and flip_dict[index]:
                        values = [float(x) if x != '' else np.nan for x in splitted_line[1:]]
                        data = np.array(values)
                        data = 2 - data
                        info.append([index] + data.tolist())
                    else:
                        info.append(splitted_line)

                    found_indices.add(index)
        f.close()

        df = pd.DataFrame(info)
        df.index = df.iloc[:, index_col]
        df.index.name = "-"
        df = df.drop(df.columns[index_col], axis=1)
        df.columns = df.iloc[header, :]
        df = df.drop(df.index[header], axis=0)
        return df

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
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Sample order  path: {}".format(self.sample_order_path))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
