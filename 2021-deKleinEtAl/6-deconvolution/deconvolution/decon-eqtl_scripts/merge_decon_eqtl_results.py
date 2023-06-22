#!/usr/bin/env python3

"""
File:         merge_decon_eqtl_results.py
Created:      2022/02/10
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
from statsmodels.stats import multitest

# Third party imports.
import pandas as pd

# Local application imports.


"""
Syntax:
./merge_decon_eqtl_results.py \
    -w /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron \
    -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -a /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz
    
./merge_decon_eqtl_results.py \
    -w /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-02-09-CortexAFR-cis-replicationOfCortexEUR-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron \
    -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-02-09-CortexAFR-cis-replicationOfCortexEUR-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -a /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/genotype_alleles.txt.gz
    
### INHIBITORY SEPERATE ###

./merge_decon_eqtl_results.py \
    -w /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected \
    -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -a /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz
    
./merge_decon_eqtl_results.py \
    -w /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-14-CortexAFR-cis-replicationOfCortexEUR-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected \
    -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-02-09-CortexAFR-cis-replicationOfCortexEUR-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -a /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/genotype_alleles.txt.gz
    
### Trans ###

./merge_decon_eqtl_results.py \
    -w /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected \
    -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -a /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz
    
./merge_decon_eqtl_results.py \
    -w /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected \
    -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -a /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz
    
./merge_decon_eqtl_results.py \
    -w /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected \
    -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -a /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz
    
./merge_decon_eqtl_results.py \
    -w /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected \
    -e /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -a /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz
"""

# Metadata
__program__ = "Check Shuffle"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.work_dir = getattr(arguments, 'work_dir')
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.alleles_path = getattr(arguments, 'alleles')

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-w",
                            "--work_dir",
                            type=str,
                            required=True,
                            help="The path to the working directory")
        parser.add_argument("-e",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix.")
        parser.add_argument("-a",
                            "--alleles",
                            type=str,
                            required=True,
                            help="The path to the alleles matrix")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Load data.")
        decon_df = self.load_file(os.path.join(self.work_dir, "deconvolutionResults.txt.gz"), header=0, index_col=0)
        geno_stats_df = self.load_file(os.path.join(self.work_dir, "geno_stats.txt.gz"), header=0, index_col=0)
        eqtl_df = self.load_file(self.eqtl_path, header=0, index_col=None)
        alleles_df = self.load_file(self.alleles_path, header=0, index_col=0)

        print(decon_df)
        print(geno_stats_df)
        print(eqtl_df)
        print(alleles_df)

        print("Merge FDR with decon results.")
        decon_df = decon_df.loc[:, [col for col in decon_df.columns if not col.endswith("_FDR")]]

        if os.path.exists(os.path.join(self.work_dir, "BH_FDR.txt.gz")):
            bh_fdr_df = self.load_file(os.path.join(self.work_dir, "BH_FDR.txt.gz"), header=0, index_col=0)
            bh_fdr_df.columns = ["{} BH-FDR".format(col) for col in bh_fdr_df.columns]
            decon_df = decon_df.merge(bh_fdr_df, left_index=True, right_index=True)
        else:
            for column in decon_df.columns:
                if not column.endswith("_pvalue"):
                    continue

                decon_df[column.replace("_pvalue", " BH-FDR")] = multitest.multipletests(decon_df[column], method='fdr_bh')[1]

        if os.path.exists(os.path.join(self.work_dir, "EMP_FDR.txt.gz")):
            emp_fdr_df = self.load_file( os.path.join(self.work_dir, "EMP_FDR.txt.gz"), header=0, index_col=0)
            emp_fdr_df.columns = ["{} Perm-FDR".format(col) for col in emp_fdr_df.columns]
            decon_df = decon_df.merge(emp_fdr_df, left_index=True, right_index=True)
        print(decon_df)

        print("Pre-process")
        # Change the col names.
        new_columns = []
        for col in decon_df.columns:
            if "pvalue" in col:
                new_columns.append(col.replace("_", " "))
            elif "GT" in col:
                new_columns.append("{} interaction beta".format(col.split("_")[1].split(":")[0]))
            elif "Beta" in col:
                new_columns.append("{} beta".format(col.split("_")[1]))
            else:
                new_columns.append(col)
        decon_df.columns = new_columns

        # Pre-process the alleles data.
        snp_to_alles_dict = dict(zip(alleles_df.index, alleles_df["Alleles"]))

        # Flip the main z-score to the decon-eQTL allele assessed.
        zscore_column = "OverallZScore"
        if zscore_column not in eqtl_df.columns:
            if "MetaPZ" in eqtl_df.columns:
                zscore_column = "MetaPZ"
            else:
                print("OverallZScore / MetaPZ not in eQTL column.")
                exit()
        eqtl_df["Alleles"] = eqtl_df["SNPName"].map(snp_to_alles_dict)
        eqtl_df["DeconAlleleAssessed"] = eqtl_df["Alleles"].str.split("/", n=1, expand=True)[1]
        eqtl_df["eQTLZscoreFlip"] = eqtl_df["AlleleAssessed"] != eqtl_df["DeconAlleleAssessed"]
        eqtl_df["ZscoreFlipped"] = eqtl_df[zscore_column] * eqtl_df["eQTLZscoreFlip"].map({True: -1, False: 1})

        if "GeneSymbol" not in eqtl_df:
            gene_trans_df = pd.read_csv("/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.PsychENCODEGenesExtended.txt.gz",
                sep="\t", index_col=None)
            print(gene_trans_df)
            gene_trans_dict = dict(zip(gene_trans_df["ArrayAddress"], gene_trans_df["Symbol"]))
            eqtl_df["GeneSymbol"] = eqtl_df["ProbeName"].map(gene_trans_dict)

        # Select the columns we want.
        eqtl_df = eqtl_df.loc[:, ["ProbeName", "GeneSymbol", "SNPName", "Alleles", "DeconAlleleAssessed", "ZscoreFlipped"]].copy()
        eqtl_df.columns = ["Gene", "Gene symbol", "SNP", "Alleles", "Allele assessed", "Overall z-score"]
        eqtl_df.index = eqtl_df["Gene"] + "_" + eqtl_df["SNP"]

        # Merge eQTL with geno stats.
        geno_stats_df = geno_stats_df.loc[:, ["N", "HW pval", "MA", "MAF"]].copy()
        geno_stats_df.columns = ["N", "HW pval", "Minor allele", "MAF"]
        geno_stats_df = geno_stats_df.groupby(geno_stats_df.index).first()
        geno_stats_df.dropna(inplace=True)
        eqtl_df = eqtl_df.merge(geno_stats_df, left_on="SNP", right_index=True, how="left")
        eqtl_df = eqtl_df.loc[:, ["Gene", "Gene symbol", "SNP", "Alleles", "Allele assessed", "N", "HW pval", "Minor allele", "MAF", "Overall z-score"]]

        # Merge.
        df = eqtl_df.merge(decon_df, left_index=True, right_index=True, how="inner")
        print(df)

        print("Saving output")
        self.save_file(df=df,
                       outpath=os.path.join(self.work_dir, "merged_decon_results.txt.gz"),
                       index=False)
        self.save_file(df=df,
                       outpath=os.path.join(self.work_dir, "merged_decon_results.xlsx"),
                       index=False,
                       sheet_name="Decon-eQTL Cortex EUR")

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
    def save_file(df, outpath, header=True, index=True, sep="\t", na_rep="NA",
                  sheet_name="Sheet1"):
        if outpath.endswith('xlsx'):
            df.to_excel(outpath,
                        sheet_name=sheet_name,
                        na_rep=na_rep,
                        header=header,
                        index=index)
        else:
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
        print("  > Working directory: {}".format(self.work_dir))
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype alleles path: {}".format(self.alleles_path))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
