#!/usr/bin/env python3

"""
File:         pre_process_decon_expression_matrix.py
Created:      2021/07/06
Last Changed: 2022/04/04
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
from functools import reduce
import argparse
import time
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.regression.linear_model import OLS
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Pre-process Decon-eQTL Expression Matrix"
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
./pre_process_decon_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR -p GTE-EUR- -of CortexEUR

./pre_process_decon_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR -p GTE-EUR- -e ENA GVEX -of CortexEUR_noENA_noGVEX

./pre_process_decon_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR -p GTE-AFR- -e ENA -of CortexAFR_noENA

./pre_process_decon_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR -p GTE-EUR- -of CortexEUR-cis-Normalised

./pre_process_decon_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR -p GTE-AFR- -of CortexAFR-cis-Normalised

### NEW ###

./pre_process_decon_expression_matrix.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-EUR/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz \
    -ra /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis/create_correction_matrix/technical_covariates_table_top20.txt.gz \
    -m /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis/create_correction_matrix/mds_covariates_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz \
    -of 2021-12-07-CortexEUR-cis-Normalised

./pre_process_decon_expression_matrix.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-AFR/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz \
    -ra /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/create_correction_matrix/technical_covariates_table_top20.txt.gz \
    -m /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/create_correction_matrix/mds_covariates_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/combine_gte_files/SampleToDataset.txt.gz \
    -of 2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved
    
./pre_process_decon_expression_matrix.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-EURandAFR-noENA/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz \
    -ra /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/technical_covariates_table_subset.txt.gz \
    -m /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/mds_covariates_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-Normalised
    
./pre_process_decon_expression_matrix.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-EURandAFR-noENA/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz \
    -ra /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/technical_covariates_table_subset.txt.gz \
    -m /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/mds_covariates_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -rpc 100 \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-Normalised
    
./pre_process_decon_expression_matrix.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-EURandAFR-noENA/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz \
    -ra /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/technical_covariates_table_subset.txt.gz \
    -m /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/mds_covariates_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-Normalised
    
./pre_process_decon_expression_matrix.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-EURandAFR-noENA/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz \
    -ra /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/technical_covariates_table_subset.txt.gz \
    -m /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/mds_covariates_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -rpc 80 \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-Normalised
    
### Harm-Jan files ###

./pre_process_decon_expression_matrix.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/output-PCATitration-MDSCorrectedPerDsCovarOverall-cortex-EURandAFR-noENA/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz \
    -ra /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/technical_covariates_table_subset_hj.txt.gz \
    -m /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_correction_matrix/mds_covariates_table_hj.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-trans-Normalised-0PCs-HJFiles
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.rna_alignment_path = getattr(arguments, 'rna_alignment')
        self.sex_path = getattr(arguments, 'sex')
        self.mds_path = getattr(arguments, 'mds')
        self.expression_pcs_path = getattr(arguments, 'expression_pcs')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.remove_n_pcs = getattr(arguments, 'rm_n_pcs')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        if outdir is None:
            outdir = str(Path(__file__).parent.parent)
        self.plot_outdir = os.path.join(outdir, 'pre_process_decon_expression_matrix', outfolder, 'plot')
        self.file_outdir = os.path.join(outdir, 'pre_process_decon_expression_matrix', outfolder, 'data')
        for outdir in [self.plot_outdir, self.file_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        self.palette = {
            "AMPAD-MAYO-V2": "#9C9FA0",
            "CMC_HBCC_set2": "#0877B4",
            "GTEx": "#0FA67D",
            "AMPAD-ROSMAP-V2": "#6950A1",
            "BrainGVEX-V2": "#48B2E5",
            "TargetALS": "#D5C77A",
            "AMPAD-MSBB-V2": "#5CC5BF",
            "NABEC-H610": "#6D743A",
            "LIBD_1M": "#808080",
            "LIBD_5M": "#808080",
            "ENA": "#D46727",
            "LIBD_h650": "#808080",
            "GVEX": "#48B2E5",
            "NABEC-H550": "#6D743A",
            "CMC_HBCC_set3": "#0877B4",
            "UCLA_ASD": "#F36D2A",
            "CMC": "#EAE453",
            "CMC_HBCC_set1": "#0877B4",
            "Braineac": "#E49D26",
            "Bipseq_1M": "#000000",
            "Bipseq_h650": "#000000",
            "Brainseq": "#C778A6"
        }

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
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the data matrix.")
        parser.add_argument("-ra",
                            "--rna_alignment",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the RNAseq alignment metrics"
                                 " matrix.")
        parser.add_argument("-s",
                            "--sex",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the sex matrix.")
        parser.add_argument("-m",
                            "--mds",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the mds matrix.")
        parser.add_argument("-ep",
                            "--expression_pcs",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the expression pcs matrix.")
        parser.add_argument("-rpc",
                            "--rm_n_pcs",
                            type=int,
                            required=False,
                            default=0,
                            help="The number of expression PCs to remove. "
                                 "Default: 0.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
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

        # Construct the output filename.
        filename = os.path.basename(self.data_path).replace(".gz", "").replace(".txt", "")

        # Load sample-dataset file.
        print("Loading sample-to-dataset.")
        std_df = self.load_file(self.std_path, header=0, index_col=None)

        # Pre-process data.
        print("Pre-processing samples-to-dataset.")
        samples = std_df.iloc[:, 0].values.tolist()
        sample_to_dataset = dict(zip(std_df.iloc[:, 0], std_df.iloc[:, 1]))

        dataset_sample_counts = list(zip(*np.unique(std_df.iloc[:, 1], return_counts=True)))
        dataset_sample_counts.sort(key=lambda x: -x[1])
        datasets = [csc[0] for csc in dataset_sample_counts]
        print("\tDatasets: {} [N = {}]".format(", ".join(datasets), len(datasets)))

        dataset_s = std_df.copy()
        dataset_s.set_index(std_df.columns[0], inplace=True)
        dataset_df = pd.get_dummies(dataset_s, prefix="", prefix_sep="")
        dataset_df = dataset_df.loc[:, datasets]

        # Load data.
        print("Loading data.")
        df = self.load_file(self.data_path, header=0, index_col=0)

        print("Step 1: sample selection.")
        print("\tUsing {}/{} samples.".format(len(samples), df.shape[1]))
        df = df.loc[:, samples]

        print("Step 2: remove probes with zero variance.")
        mask = df.std(axis=1) != 0
        print("\tUsing {}/{} probes.".format(np.sum(mask), np.size(mask)))
        df = df.loc[mask, :]

        print("Step 3: log2 transform.")
        min_value = df.min(axis=1).min()
        if min_value <= 0:
            df = np.log2(df - min_value + 1)
        else:
            df = np.log2(df + 1)

        print("Step 4: PCA analysis.")
        _ = self.pca(df=df,
                     sample_to_dataset=sample_to_dataset,
                     plot_appendix="_1_Log2Transformed")

        print("Step 5: Construct correction matrix.")
        ram_df = None
        if self.rna_alignment_path is not None:
            ram_df = self.load_file(self.rna_alignment_path, header=0, index_col=0)
            ram_df = ram_df.loc[:, samples].T

        sex_df = None
        if self.sex_path is not None:
            sex_df = self.load_file(self.sex_path, header=0, index_col=0)
            sex_df = sex_df.loc[:, samples].T

        mds_df = None
        if self.mds_path is not None:
            mds_df = self.load_file(self.mds_path, header=0, index_col=0)
            mds_df = mds_df.loc[:, samples].T

        expr_pcs_df = None
        if self.expression_pcs_path is not None:
            expr_pcs_df = self.load_file(self.expression_pcs_path, header=0, index_col=0)
            expr_pcs_df = expr_pcs_df.loc[:, samples].T

        correction_df = self.prepare_correction_matrix(ram_df=ram_df,
                                                       sex_df=sex_df,
                                                       mds_df=mds_df,
                                                       expr_pcs_df=expr_pcs_df,
                                                       dataset_df=dataset_df)

        print("\tSaving file.")
        self.save_file(df=correction_df, outpath=os.path.join(self.file_outdir, "correction_matrix1.txt.gz"))

        print("Step 6: remove technical covariates OLS.")
        corrected_df = self.calculate_residuals(df=df, correction_df=correction_df)

        print("Step 7: PCA analysis.")
        pca_components_df = self.pca(df=corrected_df,
                                     sample_to_dataset=sample_to_dataset,
                                     n_components=self.remove_n_pcs,
                                     plot_appendix="_2_Log2Transformed_CovariatesRemovedOLS")

        prefn_corrected_df = None
        correction_df2 = None
        if self.remove_n_pcs > 0:
            correction_df2 = self.prepare_correction_matrix(expr_pcs_df=pca_components_df)

            print("\tSaving file.")
            self.save_file(df=correction_df2, outpath=os.path.join(self.file_outdir, "correction_matrix2.txt.gz"))

            prefn_corrected_df = corrected_df.copy()
            print(corrected_df)
            print(correction_df2)

        print("Step 8: force normalise.")
        normal_df = self.force_normalise(df=corrected_df)

        print("Step 9: PCA analysis.")
        self.pca(df=normal_df,
                 sample_to_dataset=sample_to_dataset,
                 plot_appendix="_3_Log2Transformed_CovariatesRemovedOLS_ForceNormalised")

        print("Step 10: exp added.")
        decon_df = np.power(2, normal_df)

        print("\tSaving file.")
        self.save_file(df=decon_df, outpath=os.path.join(self.file_outdir, "{}.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt.gz".format(filename)))

        print("Step 11: PCA analysis.")
        self.pca(df=decon_df,
                 sample_to_dataset=sample_to_dataset,
                 plot_appendix="_4_Log2Transformed_CovariatesRemovedOLS_ForceNormalised_ExpAdded")

        del corrected_df, normal_df, decon_df

        if self.remove_n_pcs > 0:
            print(prefn_corrected_df)
            print(correction_df2)

            print("Step 8: remove N expression PCs.")
            corrected_df = self.calculate_residuals(df=prefn_corrected_df,
                                                    correction_df=correction_df2)

            print("Step 9: force normalise.")
            normal_df = self.force_normalise(df=corrected_df)

            print("Step 10: exp added.")
            decon_df = np.power(2, normal_df)

            print("\tSaving file.")
            self.save_file(df=decon_df, outpath=os.path.join(self.file_outdir,
                                                             "{}.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.{}PCsRemovedOLS.ForceNormalised.ExpAdded.txt.gz".format(filename, self.remove_n_pcs)))

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

    @staticmethod
    def reverse_dict(dict):
        out_dict = {}
        seen_keys = set()
        for key, value in dict.items():
            if key in seen_keys:
                print("Key {} has muiltiple values.".format(key))
            seen_keys.add(key)

            if value in out_dict.keys():
                keys = out_dict[value]
                keys.append(key)
                out_dict[value] = keys
            else:
                out_dict[value] = [key]

        return out_dict

    def prepare_correction_matrix(self, ram_df=None, sex_df=None, mds_df=None,
                                  expr_pcs_df=None, dataset_df=None):
        ram_df_subset_df = None
        if ram_df is not None:
            # Remove columns without variance and filter the RNAseq alignment
            # metrics on VIF.
            ram_df_subset_df = ram_df.copy()
            ram_df_subset_df = self.remove_multicollinearity(ram_df_subset_df.loc[:, ram_df_subset_df.std(axis=0) != 0])

        # Merge the RNAseq alignment metrics with the sex and genotype
        # MDS components.
        correction_df = None
        if ram_df is not None:
            correction_df = ram_df

        if sex_df is not None:
            if correction_df is not None:
                correction_df = ram_df_subset_df.merge(sex_df, left_index=True, right_index=True)
            else:
                correction_df = sex_df

        if mds_df is not None:
            if correction_df is not None:
                correction_df = correction_df.merge(mds_df, left_index=True, right_index=True)
            else:
                correction_df = mds_df

        if expr_pcs_df is not None:
            if correction_df is not None:
                correction_df = correction_df.merge(expr_pcs_df, left_index=True, right_index=True)
            else:
                correction_df = expr_pcs_df

        # Add the dataset dummies but exclude the dataset with the highest
        # number of samples.
        if dataset_df is not None:
            if correction_df is not None:
                correction_df = correction_df.merge(dataset_df.iloc[:, 1:], left_index=True, right_index=True)
            else:
                correction_df = dataset_df.iloc[:, 1:]

        # Add intercept.
        correction_df.insert(0, "INTERCEPT", 1)
        correction_df.index.name = "-"

        return correction_df

    def remove_multicollinearity(self, df, threshold=0.9999):
        indices = np.arange(df.shape[1])
        max_vif = np.inf
        while len(indices) > 1 and max_vif > threshold:
            vif = np.array([self.calc_ols_rsquared(df.iloc[:, indices], ix) for ix in range(len(indices))])
            max_vif = max(vif)

            if max_vif > threshold:
                max_index = np.where(vif == max_vif)[0][0]
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
    def force_normalise(df):
        return pd.DataFrame(stats.norm.ppf((df.rank(axis=1, ascending=True) - 0.5) / df.shape[1]), index=df.index, columns=df.columns)

    def pca(self, df, sample_to_dataset, n_components=2, plot_appendix=""):
        if n_components < 2:
            n_components = 2

        # samples should be on the columns and genes on the rows.
        zscores = (df - df.mean(axis=0)) / df.std(axis=0)
        pca = PCA(n_components=n_components)
        pca.fit(zscores)
        expl_variance = {"PC{}".format(i+1): pca.explained_variance_ratio_[i] * 100 for i in range(2)}
        components_df = pd.DataFrame(pca.components_)
        components_df.index = ["Comp{}".format(i + 1) for i, _ in enumerate(components_df.index)]
        components_df.columns = df.columns
        components_df = components_df.T

        print("\tPlotting PCA.")
        plot_df = components_df.copy()
        plot_df["hue"] = plot_df.index.map(sample_to_dataset)
        self.plot(df=plot_df,
                  x="Comp1",
                  y="Comp2",
                  hue="hue",
                  palette=self.palette,
                  xlabel="PC1 [{:.2f}%]".format(expl_variance["PC1"]),
                  ylabel="PC2 [{:.2f}%]".format(expl_variance["PC2"]),
                  title="PCA - eigenvectors",
                  filename="eigenvectors_plot{}".format(plot_appendix))

        return components_df

    def plot(self, df, x="x", y="y", hue=None, palette=None, xlabel=None,
             ylabel=None, title="", filename="PCA_plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, gridspec_kw={"width_ratios": [0.9, 0.1]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        s=100,
                        linewidth=0,
                        legend=None,
                        palette=palette,
                        ax=ax1)

        ax1.set_title(title,
                      fontsize=20,
                      fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_xlabel(xlabel,
                       fontsize=14,
                       fontweight='bold')

        if palette is not None:
            handles = []
            for label, color in palette.items():
                if label in df[hue].values.tolist():
                    handles.append(mpatches.Patch(color=color, label=label))
            ax2.legend(handles=handles, loc="center")

        #fig.savefig(os.path.join(self.plot_outdir, "{}.pdf".format(filename)))
        fig.savefig(os.path.join(self.plot_outdir, "{}.png".format(filename)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Data: {}".format(self.data_path))
        print("  > RNAseq alignment metrics: {}".format(self.rna_alignment_path))
        print("  > Sex: {}".format(self.sex_path))
        print("  > MDS: {}".format(self.mds_path))
        print("  > Expression PCs: {}".format(self.expression_pcs_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Remove N-PCs: {}".format(self.remove_n_pcs))
        print("  > Plot output directory: {}".format(self.plot_outdir))
        print("  > File output directory: {}".format(self.file_outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
