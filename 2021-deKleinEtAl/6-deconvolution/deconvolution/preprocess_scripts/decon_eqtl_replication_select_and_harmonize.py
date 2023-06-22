#!/usr/bin/env python3

"""
File:         decon_eqtl_replication_select_and_harmonize.py
Created:      2021/06/23
Last Changed: 2021/08/03
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
__program__ = "Decon-eQTL Select and Harmonize"
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
./decon_eqtl_replication_select_and_harmonize.py -re /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/replicateCortex-EUR/eQTLProbesFDR0.05-ProbeLevel.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_alleles.txt.gz -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/create_matrices/genotype_table.txt.gz -al /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/create_matrices/genotype_alleles.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexAFR-cis/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt.gz -of CortexAFR-cis-Replication-EUR

./decon_eqtl_replication_select_and_harmonize.py -re /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/replicateCortex-EUR/eQTLProbesFDR0.05-ProbeLevel.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_alleles.txt.gz -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/create_matrices/genotype_table.txt.gz -al /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/create_matrices/genotype_alleles.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexAFR-cis-Normalised/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt.gz -of CortexAFR-cis-Replication-EUR-Normalised

### Cortex AFR Replication ###

./decon_eqtl_replication_select_and_harmonize.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt.gz \
    -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_alleles.txt.gz \
    -of 2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved 
    
./decon_eqtl_replication_select_and_harmonize.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-02-09-CortexAFR-cis-replicationOfCortexEUR-NegativeToZero-DatasetAndRAMCorrected \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt.gz \
    -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_alleles.txt.gz \
    -of 2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved 
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_dir = getattr(arguments, 'data')
        self.expr_path = getattr(arguments, 'expression')
        self.discovery_alleles_path = getattr(arguments, 'discovery_alleles')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        if outdir is None:
            outdir = str(Path(__file__).parent.parent)
        self.outdir = os.path.join(outdir, 'decon_eqtl_replication_select_and_harmonize', outfolder)
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
        parser.add_argument("-da",
                            "--discovery_alleles",
                            type=str,
                            required=True,
                            help="The path to the discovery alleles matrix")
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to Matrix Preparation input dir.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
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

        print("Loading eQTL file.")
        eqtl_df = self.load_file(os.path.join(self.data_dir, "combine_eqtlprobes", "eQTLprobes_combined.txt.gz"), header=0, index_col=None)

        print("Loading genotype data files")
        geno_df = self.load_file(os.path.join(self.data_dir, "create_matrices", "genotype_table.txt.gz"), header=0, index_col=0)
        allele_df = self.load_file(os.path.join(self.data_dir, "create_matrices", "genotype_alleles.txt.gz"), header=0, index_col=0)

        print("\tChecking order")
        if eqtl_df["SNPName"].tolist() != list(geno_df.index):
            print("Genotype matrix does not match eQTL matrix.")
            exit()
        if eqtl_df["SNPName"].tolist() != list(allele_df.index):
            print("Alleles matrix does not match eQTL matrix.")
            exit()

        print("Preprocessing expression matrix")
        expr_df = self.load_file(inpath=self.expr_path, header=0, index_col=0)
        expr_df = expr_df.groupby(expr_df.index).first()
        for probe_name in eqtl_df["ProbeName"]:
            if probe_name not in expr_df.index:
                expr_df.loc[probe_name, :] = np.nan
        expr_df = expr_df.loc[eqtl_df["ProbeName"], :]
        expr_df.index.name = None

        print("\tChecking order")
        if eqtl_df["ProbeName"].tolist() != list(expr_df.index):
            print("Expression matrix does not match eQTL matrix.")
            exit()

        print("Filter eQTL file on present data.")
        geno_df[geno_df == -1] = np.nan
        geno_available_mask = np.logical_and(~geno_df.isnull().all(1).to_numpy(dtype=bool), ~allele_df.isnull().all(1).to_numpy(dtype=bool))
        expr_available_mask = ~expr_df.isnull().all(1).to_numpy(dtype=bool)
        mask = np.logical_and(geno_available_mask, expr_available_mask)

        eqtl_df = eqtl_df.loc[mask, :]
        geno_df = geno_df.loc[mask, :]
        allele_df = allele_df.loc[mask, :]
        expr_df = expr_df.loc[mask, :]

        print("Creating flip mask.")
        discovery_alleles_df = self.load_file(self.discovery_alleles_path, header=0, index_col=0)
        discovery_alleles_df = discovery_alleles_df.groupby(discovery_alleles_df.index).first()
        replication_alleles_df = allele_df.groupby(allele_df.index).first()
        ma_df = replication_alleles_df[["MinorAllele"]].merge(discovery_alleles_df[["MinorAllele"]], left_index=True, right_index=True)
        ma_df["flip"] = ma_df.iloc[:, 0] != ma_df.iloc[:, 1]
        ma_df = ma_df.loc[eqtl_df["SNPName"], :]
        flip_mask = ma_df["flip"].to_numpy(dtype=bool)
        del ma_df, replication_alleles_df

        print("Flipping genotype encoding.")
        allele_df = discovery_alleles_df.loc[allele_df.index, :]

        geno_m = geno_df.to_numpy(dtype=np.float64)
        missing_mask = geno_m == -1
        geno_m[flip_mask, :] = 2 - geno_m[flip_mask, :]
        geno_m[missing_mask] = -1
        geno_df = pd.DataFrame(geno_m, index=geno_df.index, columns=geno_df.columns)
        geno_df.index.name = None
        del geno_m, missing_mask

        self.save_file(df=eqtl_df, outpath=os.path.join(self.outdir, "eQTLprobes_combined.txt.gz"), index=False)
        self.save_file(df=geno_df, outpath=os.path.join(self.outdir, "genotype_table.txt.gz"))
        self.save_file(df=allele_df, outpath=os.path.join(self.outdir, "genotype_alleles.txt.gz"))
        self.save_file(df=expr_df, outpath=os.path.join(self.outdir, "expression_table.txt.gz"))

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
        print("  > Data directory: {}".format(self.data_dir))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Discovery alleles path: {}".format(self.discovery_alleles_path))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
