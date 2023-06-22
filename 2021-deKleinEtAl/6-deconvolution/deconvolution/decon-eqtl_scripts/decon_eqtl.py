#!/usr/bin/env python3

"""
File:         decon_eqtl.py
Created:      2021/07/12
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
import itertools
import warnings
import argparse
import random
import time
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.special import betainc
from scipy.special import ndtri
from statsmodels.stats import multitest

# Local application imports.

"""
Syntax:
./decon_eqtl.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR-Normalised/genotype_table.txt -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexAFR-cis-Normalised/data/SampleToDataset.txt.gz -of CortexAFR-cis-Replication-EUR-NormalisedMAF5 -maf 0.05

/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/SampleToDataset.txt.gz -of CortexEUR-cis-NormalisedMAF5 -maf 0.05  

/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table_CNS7.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/SampleToDataset.txt.gz -of CortexEUR-cis-NormalisedMAF5-CNS7Profile -maf 0.05 

/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table_CNS7_NoPericytes.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/SampleToDataset.txt.gz -of CortexEUR-cis-NormalisedMAF5-LimitedConfigs-NewProfileNoPericytes -maf 0.05 -ac limited

/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/SampleToDataset.txt.gz -maf 0.05 -ac limited -of CortexEUR-cis-NormalisedMAF5-ALlConfigs

/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/SampleToDataset.txt.gz -maf 0.05 -ac limited -of CortexEUR-cis-NormalisedMAF5-LimitedConfigs 

/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table_PsychENCODE.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/SampleToDataset.txt.gz -maf 0.05 -ac limited -of CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile

./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/OLD/2021-09-22/CortexEUR-cis/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/OLD/2021-09-22/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/SampleToDataset.txt.gz \
    -of CortexEUR-cis-NormalisedMAF5-LimitedConfigs \
    -maf 5.05 \
    -ac limited \
    -of 2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-OldProfile

./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -ac limited \
    -of 2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDevCT

./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table_NoInhibitory.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -ac limited \
    -of 2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDevNoInhibitoryCT


### Cortex EUR ###

./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -ac limited \
    -of 2021-12-07-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron
    
./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -ac limited \
    -of 2022-01-25-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron
    
./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -ac limited \
    -of 2022-01-25-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-NegativeToZero-DatasetAndRAMCorrected
    
./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -of 2022-01-27-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron
    
./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -of 2022-05-04-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-test
    
### Cortex AFR Replication of EUR ###    

./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/expression_table.txt.gz \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -ac limited \
    -of 2021-12-22-CortexAFR-replicationOfCortexEUR20211207-ProbesWithZeroVarianceRemoved-InhibitorySummedWithOtherNeuron
    
./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/expression_table.txt.gz \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-02-09-CortexAFR-cis-replicationOfCortexEUR-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-02-09-CortexAFR-cis-replicationOfCortexEUR-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -of 2022-02-09-CortexAFR-cis-replicationOfCortexEUR-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron

./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/expression_table.txt.gz \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-02-09-CortexAFR-cis-replicationOfCortexEUR-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-02-09-CortexAFR-cis-replicationOfCortexEUR-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -of 2022-03-14-CortexAFR-cis-replicationOfCortexEUR-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected

### Trans ###

./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected
    
./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.80PCsRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected
    
./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected

./decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.100PCsRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -of 2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected

"""

# Metadata
__program__ = "Decon-eQTL"
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
        self.allele_configs = getattr(arguments, 'allele_configurations')
        self.call_rate = getattr(arguments, 'call_rate')
        self.hw_pval = getattr(arguments, 'hardy_weinberg_pvalue')
        self.maf = getattr(arguments, 'minor_allele_frequency')
        self.zscore_cutoff = getattr(arguments, 'zscore_cutoff')
        self.nrows = getattr(arguments, 'rows')
        self.n_permutations = getattr(arguments, 'permutations')
        self.permutation_index_offset = getattr(arguments, 'permutation_index_offset')
        self.leading_zeros = getattr(arguments, 'permutation_leading_zeros')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        if outdir is None:
            outdir = str(Path(__file__).parent)
        self.outdir = os.path.join(outdir, "decon_eqtl", outfolder)
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
                            default=None,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-na",
                            "--genotype_na",
                            type=int,
                            required=False,
                            default=-1,
                            help="The genotype value that equals a missing "
                                 "value. Default: -1.")
        parser.add_argument("-ac",
                            "--allele_configurations",
                            type=str,
                            choices=["complete", "limited"],
                            default="complete",
                            help="The number of allele encoding configurations"
                                 "to test. Either 'complete' for all or"
                                 "'limited' for restricting a max of 1 flipped"
                                 "allele. Default: complete.")
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
        parser.add_argument("-zc",
                            "--zscore_cutoff",
                            type=int,
                            required=False,
                            default=None,
                            help="Remove outliers with z-score >X on one or"
                                 "more cell fractions. Default: None.")
        parser.add_argument("-r",
                            "--rows",
                            type=int,
                            default=None,
                            help="The number of rows to analyze.")
        parser.add_argument("-p",
                            "--permutations",
                            type=int,
                            default=0,
                            help="The number of permutations to run. "
                                 "Default: 0.")
        parser.add_argument("-po",
                            "--permutation_index_offset",
                            type=int,
                            default=0,
                            help="The first index of the first permutation "
                                 "run. Default: 0.")
        parser.add_argument("-plz",
                            "--permutation_leading_zeros",
                            type=int,
                            default=0,
                            help="The number of leading zeros to print for the "
                                 "permutation output files. Default: 0.")
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
        print("Loading genotype data.")
        geno_df = self.load_file(self.geno_path, header=0, index_col=0, nrows=self.nrows)

        if self.std_path is not None:
            print("Loading dataset info.")
            std_df = self.load_file(self.std_path, header=0, index_col=None)

            # Validate that the input data matches.
            self.validate_data(std_df=std_df,
                               geno_df=geno_df)
        else:
            # Create sample-to-dataset file with all the samples having the
            # same dataset.
            std_df = pd.DataFrame({"sample": geno_df.columns, "dataset": "None"})

        cc_df = None
        zscore_mask = pd.Series(True, index=std_df.iloc[:, 0])
        if self.zscore_cutoff is not None:
            print("Loading cell fraction data.")
            cc_df = self.load_file(self.cc_path, header=0, index_col=0)

            # Transpose if need be. We want samples always as columns.
            if cc_df.shape[0] == geno_df.shape[1]:
                cc_df = cc_df.T

            # Validate that the input data matches.
            self.validate_data(std_df=std_df,
                               geno_df=geno_df,
                               cc_df=cc_df)

            # Transform to z-scores.
            cc_zscore_df = cc_df.subtract(cc_df.mean(axis=1), axis=0).divide(cc_df.std(axis=1), axis=0)
            zscore_mask = cc_zscore_df.abs().max(axis=0) < self.zscore_cutoff
            print("\t {:,} samples removed.".format(zscore_mask.size - zscore_mask.sum()))

            # Filter matrices.
            geno_df = geno_df.loc[:, zscore_mask]
            std_df = std_df.loc[zscore_mask.to_numpy(dtype=bool), :]
            cc_df = cc_df.loc[:, zscore_mask]

        print("Checking dataset sample sizes.")
        dataset_sample_counts = list(zip(*np.unique(std_df.iloc[:, 1], return_counts=True)))
        dataset_sample_counts.sort(key=lambda x: -x[1])
        datasets = [x[0] for x in dataset_sample_counts]
        max_dataset_length = np.max([len(str(dataset[0])) for dataset in dataset_sample_counts])
        for dataset, sample_size in dataset_sample_counts:
            print("\t{:{}s}  {:,} samples".format(dataset, max_dataset_length, sample_size))
        if dataset_sample_counts[-1][1] <= 1:
            print("\t  One or more datasets have a smaller sample "
                  "size than recommended. Consider excluded these.")
            exit()

        print("Calculating genotype call rate per dataset.")
        geno_df, call_rate_df = self.calculate_call_rate(geno_df=geno_df,
                                                         std_df=std_df,
                                                         datasets=datasets)
        call_rate_n_skipped = (call_rate_df.min(axis=1) < self.call_rate).sum()
        if call_rate_n_skipped > 0:
            print("\t{:,} eQTLs have had dataset(s) filled with NaN "
                  "values due to call rate threshold ".format(call_rate_n_skipped))

        print("Calculating genotype stats for inclusing criteria.")
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

        if self.permutation_index_offset == 0:
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
        if cc_df is None:
            print("Loading cell fraction data.")
            cc_df = self.load_file(self.cc_path, header=0, index_col=0)

            # Transpose if need be. We want samples always as columns.
            if cc_df.shape[0] == geno_df.shape[1]:
                cc_df = cc_df.T

        # Remove columns with all nan. This line is just because
        # the expression file I used somehow had one column with nan's.
        expr_df.dropna(axis='columns', how='all', inplace=True)

        # Select eQTL rows that meet requirements.
        geno_df = geno_df.loc[combined_keep_mask, :]
        expr_df = expr_df.loc[combined_keep_mask, :]

        # Remove samples if they did not pass z-score threshold.
        expr_df = expr_df.loc[:, zscore_mask]

        print("\tValidating input.")
        self.validate_data(std_df=std_df,
                           geno_df=geno_df,
                           expr_df=expr_df,
                           cc_df=cc_df)
        print("", flush=True)

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
        std_m = std_df.to_numpy(object)

        # Replace missing values with nan
        geno_m[geno_m == self.genotype_na] = np.nan

        # Save properties.
        eqtl_indices = expr_df.index + "_" + geno_df.index
        cell_types_indices = cc_df.index.to_numpy(dtype=object)

        del geno_df, expr_df, cc_df, std_df

        # Print info.
        n_eqtls = geno_m.shape[0]
        n_samples = geno_m.shape[1]
        n_covariates = cc_m.shape[0]
        print("Summary stats:")
        print("\tN-eQTLs: {:,}".format(n_eqtls))
        print("\tN-samples: {:,}".format(n_samples))
        print("\tN-covariates: {:,}".format(n_covariates))
        print("\tN-datasets: {:,}".format(len(datasets)))
        print("", flush=True)

        ########################################################################

        print("### STEP 4 ###")
        print("Analyzing eQTLs.")

        # Create a list of possible genotype encoding configuration. True means
        # we change the encoding by 2 - value. False means we do nothing.
        alt_model_configs = self.create_model_configs(n=n_covariates, type=self.allele_configs)
        null_model_configs = self.create_model_configs(n=n_covariates - 1, type=self.allele_configs)

        # Calculate and print some info about the analyses to be performed.
        n_configurations_alt = len(alt_model_configs)
        n_configurations_null = len(null_model_configs)
        n_models = n_eqtls * (n_configurations_alt + n_covariates * n_configurations_null)
        print("\tN-configurations (full model): {:,}".format(n_configurations_alt))
        print("\tN-configurations (cell type model): {:,}".format(n_configurations_null))
        print("\tN-models: {:,}".format(n_models))
        print("")

        # Initializing output matrices / arrays.
        real_pvalues_m = np.empty((n_eqtls, n_covariates), dtype=np.float64)
        betas_alt_m = np.empty((n_eqtls, n_covariates * 2), dtype=np.float64)
        rss_null_m = np.empty((n_eqtls, n_covariates), dtype=np.float64)

        # Save the degrees of freedom the alternative model.
        df = n_covariates * 2

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
            mask = ~np.isnan(genotype)
            n = np.sum(mask)

            # Model the alternative matrix (with the interaction term).
            # This is the matrix with expression ~ cc1 + cc2 + cc1 * geno +
            # cc2 * geno.
            betas_alt, rss_alt = \
                self.model(
                    genotype=genotype[mask],
                    expression=expr_m[row_index, mask],
                    cell_fractions=cc_m[:, mask],
                    configs=alt_model_configs,
                    n_samples=n,
                    n_covariates=n_covariates
                )

            # Save the alternative model stats.
            betas_alt_m[row_index, :] = betas_alt

            # Remove one interaction column (cc * geno) one by one and
            # determine the significance of the change in residuals sum of
            # squares with a f-test.
            for cov_index in range(n_covariates):
                # Model the null matrix (without the interaction term). In
                # this model 1 (and only 1!) of the cc * geno terms is removed.
                _, rss_null = \
                    self.model(
                        genotype=genotype[mask],
                        expression=expr_m[row_index, mask],
                        cell_fractions=cc_m[:, mask],
                        configs=null_model_configs,
                        n_samples=n,
                        n_covariates=n_covariates,
                        exclude=cov_index
                    )

                # Calculate and save the p-value.
                p_value = self.calc_p_value(rss1=rss_null,
                                            rss2=rss_alt,
                                            df1=df - 1,
                                            df2=df,
                                            n=n)
                real_pvalues_m[row_index, cov_index] = p_value

                # Save the RSS null for permutations later.
                rss_null_m[row_index, cov_index] = rss_null

        # Cap the p-values.
        real_pvalues_m[real_pvalues_m == 0] = 2.2250738585072014e-308

        # Print the number of significant hits.
        self.print_n_signif(m=real_pvalues_m, colnames=cell_types_indices, type="p-value")
        print("", flush=True)

        ########################################################################

        if self.permutation_index_offset == 0:
            print("### STEP 5 ###")
            print("Saving results.")

            decon_df = pd.DataFrame(np.hstack((real_pvalues_m, betas_alt_m)),
                                    index=eqtl_indices,
                                    columns=["{}_pvalue".format(x) for x in cell_types_indices] +
                                            ["Beta{}_{}".format(i + 1, x) for i, x in enumerate(cell_types_indices)] +
                                            ["Beta{}_{}:GT".format(len(cell_types_indices) + i + 1, x) for i, x in enumerate(cell_types_indices)])
            print(decon_df)

            if self.n_permutations == 0:
                for ct in cell_types_indices:
                    decon_df["{}_FDR".format(ct)] = multitest.multipletests(decon_df.loc[:, "{}_pvalue".format(ct)], method='fdr_bh')[1]
                fdr_m = decon_df.loc[:, ["{}_FDR".format(ct) for ct in cell_types_indices]].to_numpy()
                self.print_n_signif(m=fdr_m, colnames=cell_types_indices, type="FDR")
                del fdr_m

            self.save_file(df=decon_df, outpath=os.path.join(self.outdir, "deconvolutionResults.txt.gz"))

            # Save the beta's.
            self.save_matrix(m=betas_alt_m, outpath=os.path.join(self.outdir, "betas_alternative_model.npy"))

            del decon_df, betas_alt_m

            print("", flush=True)

        #######################################################################

        if self.n_permutations <= 0:
            exit()

        print("### STEP 6 ###")
        print("Performing permutations.")

        # Calculate and print some info about the analyses to be performed.
        n_permutation_values = n_eqtls * self.n_permutations
        n_perm_models = n_eqtls * (self.n_permutations * n_covariates * n_configurations_null)
        print("\tN-permutation p-values (per cell type): {:,}".format(n_permutation_values))
        print("\tN-models: {:,}".format(n_perm_models))
        print("")

        # Create the permutation orders.
        perm_order_m = self.create_perm_orders(
            n_permutations=self.n_permutations,
            n_samples=n_samples,
            datasets=datasets,
            std_m=std_m,
            seed_offset=self.permutation_index_offset)

        # Calculate the overlap between the genotype matrices
        perm_overlap_m = self.check_permutation_overlap(geno_m=geno_m,
                                                        perm_order_m=perm_order_m)

        geno_nan_mask = np.isnan(geno_m)
        nanfilled_geno_m = np.copy(geno_m)
        if np.sum(geno_nan_mask) > 0:
            # Calculate the dataset genotype means per eQTL and fill nan values
            # with the dataset mean.
            geno_dataset_mean_m = self.calculate_geno_mean_per_dataset(
                geno_m=geno_m,
                datasets=datasets,
                std_m=std_m)
            nanfilled_geno_m[geno_nan_mask] = geno_dataset_mean_m[geno_nan_mask]

        # Initializing output matrices / arrays.
        perm_pvalues_m = np.empty((n_eqtls, n_covariates, self.n_permutations), dtype=np.float64)
        perm_betas_alt_m = np.empty((n_eqtls, n_covariates, self.n_permutations, n_covariates * 2), dtype=np.float64)

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

            # Get the genotype arrays.
            genotype = geno_m[row_index, :]
            nanfilled_genotype = nanfilled_geno_m[row_index, :]

            # Construct the mask to remove missing values.
            mask = np.logical_and(~np.isnan(genotype), zscore_mask)
            n = np.sum(mask)

            # Loop over the covariates.
            for cov_index in range(n_covariates):

                # Perform n permutations.
                for perm_index in range(self.n_permutations):
                    # Calculate the interaction term for both allele encodings.
                    # Use the NaN-filled genotype array for this.
                    interaction_a = nanfilled_genotype * cc_m[cov_index, :]
                    interaction_flipped_a = (2 - nanfilled_genotype) * cc_m[cov_index, :]

                    # Shuffle the indices.
                    perm_order = perm_order_m[perm_index, :]
                    interaction_a = interaction_a[perm_order]
                    interaction_flipped_a = interaction_flipped_a[perm_order]

                    # Model the alternative matrix (with the interaction
                    # term) and shuffle the genotype of the interaction of
                    # interest.
                    perm_betas_alt, perm_rss_alt = \
                        self.model(
                            genotype=genotype[mask],
                            expression=expr_m[row_index, mask],
                            cell_fractions=cc_m[:, mask],
                            configs=alt_model_configs,
                            n_samples=n,
                            n_covariates=n_covariates,
                            shuffle_index=cov_index,
                            shuffle_inter=interaction_a[mask],
                            shuffle_inter_flipped=interaction_flipped_a[mask]
                        )

                    # Save the permuted alternative model stats.
                    perm_betas_alt_m[row_index, cov_index, perm_index, :] = perm_betas_alt

                    # Calculate and save the permutation p-value.
                    perm_pvalue = self.calc_p_value(
                        rss1=rss_null_m[row_index, cov_index],
                        rss2=perm_rss_alt,
                        df1=df - 1,
                        df2=df,
                        n=n)
                    perm_pvalues_m[row_index, cov_index, perm_index] = perm_pvalue

        print("", flush=True)

        # #######################################################################
        #
        # print("### STEP 7 ###")
        # print("Calculating permutation based FDR.")
        # perm_fdr_m = np.empty_like(real_pvalues_m)
        # perm_fdr_m[:] = np.nan
        # for cov_index in range(n_covariates):
        #     for row_index in range(n_eqtls):
        #         # Get the real p-value.
        #         real_pvalue = real_pvalues_m[row_index, cov_index]
        #
        #         # Get the rank of this p-value in both distributions.
        #         rank = np.sum(real_pvalues_m[:, cov_index] <= real_pvalue)
        #         perm_rank = np.sum(perm_pvalues_m[:, cov_index, :] <= real_pvalue)
        #
        #         # Calculate and safe the fdr.
        #         perm_fdr_m[row_index, cov_index] = (perm_rank / self.n_permutations) / rank
        #
        # # Cap the permutation FDR values.
        # perm_fdr_m[perm_fdr_m > 1] = 1
        # perm_fdr_m[perm_fdr_m == 0] = 2.2250738585072014e-308
        #
        # # Print the number of significant hits.
        # self.print_n_signif(m=perm_fdr_m, colnames=cell_types_indices, type="fdr-value")
        # print("", flush=True)
        #
        # ########################################################################

        print("### STEP 8 ###")
        print("Saving results.")

        file_suffix = "{}{}_until_{}{}".format("0" * self.leading_zeros,
                                               self.permutation_index_offset,
                                               "0" * self.leading_zeros,
                                               self.permutation_index_offset + self.n_permutations - 1)
        self.save_matrix(m=perm_order_m, outpath=os.path.join(self.outdir,
                                                              "perm_orders_{}.npy".format(file_suffix)))
        self.save_matrix(m=perm_overlap_m, outpath=os.path.join(self.outdir,
                                                                "perm_order_overlap_{}.npy".format(file_suffix)))
        self.save_matrix(m=perm_pvalues_m, outpath=os.path.join(self.outdir,
                                                                "permutation_pvalues_{}.npy".format(file_suffix)))
        self.save_matrix(m=perm_betas_alt_m, outpath=os.path.join(self.outdir,
                                                                  "permutation_betas_alternative_model_{}.npy".format(file_suffix)))

        lowest_pvalues_m = np.transpose(np.min(perm_pvalues_m, axis=0))
        lowest_pvalues_df = pd.DataFrame(lowest_pvalues_m,
                                         index=["permutation_{:06d}".format(self.permutation_index_offset + perm_index) for perm_index in range(self.n_permutations)],
                                         columns=["{}_pvalue".format(cell_type) for cell_type in cell_types_indices])
        self.save_file(df=lowest_pvalues_df, outpath=os.path.join(self.outdir, "lowest_permutation_pvalues_{}.txt.gz".format(file_suffix)))

        # perm_fdr_df = pd.DataFrame(perm_fdr_m, columns=["{}_FDR".format(cell_type) for cell_type in cell_types_indices])
        # self.save_file(df=perm_fdr_df, outpath=os.path.join(self.outdir, "permutation_FDR.txt.gz"))

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
                  "the sample-to-dataset link file.")
            exit()

        if expr_df is not None and expr_df.columns.tolist() != samples:
            print("\tThe expression file header does not match "
                  "the sample-to-dataset link file.")
            exit()

        if cc_df is not None and cc_df.columns.tolist() != samples:
            print("\tThe cell count file header does not match "
                  "the sample-to-dataset link file.")
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
    def create_model_configs(n, type):
        """
        Create the allele encoding configurations. In the Decon-eQTL article
        they restrict the configurations to max one opposite. This limits the
        configurations from k^2 to (2*k) + 2.

        Example for n = 3
        Output:
            [ False False False ]
            [ True False False ]
            [ False True False ]
            [ False False True ]
            [ True True True ]
            [ False True True ]
            [ True False True ]
            [ True True False ]

        """
        if type == "complete":
            return list(itertools.product([True, False], repeat=n))
        elif type == "limited":
            configurations = []

            false_array = np.zeros(n, dtype=bool)
            true_array = np.ones(n, dtype=bool)
            for value, array in zip([True, False], [false_array, true_array]):
                configurations.append(array)
                for i in range(n):
                    configuration = np.copy(array)
                    configuration[i] = value
                    configurations.append(configuration)

            return configurations
        else:
            print("Unexpected argument for allele encoding configurations.")
            exit()
            return None

    @staticmethod
    def create_perm_orders(n_permutations, n_samples, datasets, std_m,
                           seed_offset):
        """
        Shuffles an array of size n_samples into a random order (with seed)
        N times. However, this function only shuffles samples within the same
        dataset.
        """
        default_order = np.arange(n_samples, dtype=np.uint16)
        all_indices = set(default_order)
        perm_order_m = np.empty((n_permutations, n_samples), dtype=np.uint16)
        for i in range(perm_order_m.shape[0]):
            perm_order = np.copy(default_order)
            for dataset in datasets:
                dataset_mask = std_m[:, 1] == dataset
                if np.sum(dataset_mask) <= 1:
                    continue

                copy = perm_order[dataset_mask]
                random.Random(seed_offset + i).shuffle(copy)
                perm_order[dataset_mask] = copy

            if set(perm_order) != all_indices:
                print("Invalid permutation order.")
                exit()

            perm_order_m[i, :] = perm_order
        return perm_order_m

    @staticmethod
    def check_permutation_overlap(geno_m, perm_order_m):
        """
        Function to calculate how identical the permuted genotype matrix is to
        the original genotype matrix. Returns a matrix with n-eQTLs x n-permutations
        with the % overlap in the cells.
        """
        perm_overlap_m = np.empty((geno_m.shape[0], perm_order_m.shape[0]), dtype=np.float64)

        for perm_index in range(perm_order_m.shape[0]):
            perm_order = perm_order_m[perm_index, :]
            perm_geno_m = np.copy(geno_m)
            perm_geno_m = perm_geno_m[:, perm_order]
            perm_overlap_m[:, perm_index] = np.sum(np.logical_and(geno_m == perm_geno_m, ~np.isnan(geno_m)), axis=1) / geno_m.shape[1]

        return perm_overlap_m

    @staticmethod
    def calculate_geno_mean_per_dataset(geno_m, datasets, std_m):
        """
        Method for calculating the mean genotype per dataset per eQTL. Missing
        values are not included. Returns a matrix with n-eQTLs x n-samples
        with the mean genotype per dataset in the cells.
        """
        geno_dataset_mean_m = np.empty_like(geno_m)

        # Ignore RuntimeWarning: Mean of empty slice
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            # Fill in the missing values with the dataset mean.
            for dataset_index, dataset in enumerate(datasets):
                dataset_mask = std_m[:, 1] == dataset
                geno_dataset_m = geno_m[:, dataset_mask]
                dataset_mean_a = np.nanmean(geno_dataset_m, axis=1)
                geno_dataset_mean_m[:, dataset_mask] = np.tile(dataset_mean_a[:, np.newaxis], np.sum(dataset_mask))

        return geno_dataset_mean_m

    @staticmethod
    def model(genotype, expression, cell_fractions, configs, n_samples,
              n_covariates, exclude=None, shuffle_index=None,
              shuffle_inter=None, shuffle_inter_flipped=None):
        """
        Create the interaction model. Try different allele encodings and
        find the optimal configurations. Only the best configuration is stored
        and returned.
        """
        if shuffle_index is not None and (shuffle_inter is None or shuffle_inter_flipped is None):
            print("Both shuffle_inter and shuffle_inter_flipped arguments "
                  "are required if shuffle index is set.")
            exit()

        n_columns = n_covariates * 2
        if exclude is not None:
            n_columns -= 1

        # Create the model matrix. Leave the interaction columns blank for now.
        X = np.empty((n_samples, n_columns), dtype=np.float64)
        for cov_index in range(n_covariates):
            X[:, cov_index] = cell_fractions[cov_index, :]

        # Try different configurations for the genotype encoding.
        top_config = None
        top_betas = None
        top_rss = np.inf
        for config in configs:
            # Fill in the alternative matrix with the right configuration
            # of allele encoding.
            for cov_index, flip in enumerate(config):
                # If we exclude an interaction term we still have that cell
                # type fraction in the matrix. Therefore, when matching
                # interaction column position with cell fraction column
                # (cc_index) position we need to increment with 1 for all
                # cell fractions > exclude column.
                cc_index = cov_index
                if exclude is not None and cc_index >= exclude:
                    cc_index += 1

                # Use a shuffled genotype vector if we are doing a permutation
                # analysis.
                if cov_index == shuffle_index:
                    if flip:
                        X[:, n_covariates + cov_index] = shuffle_inter_flipped
                    else:
                        X[:, n_covariates + cov_index] = shuffle_inter
                else:
                    # Calculate genotype * cell fraction. If flip is true we
                    # change the allele encoding (0 = 2, 1 = 1, 2 = 0).
                    if flip:
                        X[:, n_covariates + cov_index] = (2 - genotype) * X[:, cc_index]
                    else:
                        X[:, n_covariates + cov_index] = genotype * X[:, cc_index]

            # Check if all values are positive.
            if np.min(X) < 0:
                print("Error: negative values in regression matrix.")
                exit()

            # Model the expression vector as non-negative linear combination of
            # the model matrix. Determine the beta's as well as the squared
            # Euclidean norm. Note: squared Eucledian norm is identical to
            # calculating the RSS.
            betas, rnorm = nnls(X, expression)
            rss_alt = rnorm * rnorm

            # Only safe the best configuration.
            if rss_alt < top_rss:
                top_config = config
                top_betas = betas
                top_rss = rss_alt

        # The beta's of the interaction terms are flipped if we
        # flipped the allele encoding. This makes it possible that some
        # betas are negative even though we use NNLS.
        flip_array = np.hstack((np.ones(n_covariates), np.vectorize({True: -1, False: 1}.get)(top_config)))
        top_betas = top_betas * flip_array

        # Insert NaN in betas if we excluded an interaction term.
        if exclude is not None:
            top_betas = np.insert(top_betas, n_covariates + exclude, np.nan)

        return top_betas, top_rss

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
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Genotype NaN: {}".format(self.genotype_na))
        print("  > Allele configs: {}".format(self.allele_configs))
        print("  > SNP call rate: >{}".format(self.call_rate))
        print("  > Hardy-Weinberg p-value: >{}".format(self.hw_pval))
        print("  > MAF: >{}".format(self.maf))
        print("  > Cell fraction z-score cut-off: >={}".format(self.zscore_cutoff))
        print("  > N rows: {}".format(self.nrows))
        print("  > N permutations: {}".format(self.n_permutations))
        print("  > Permutation index offset: {}".format(self.permutation_index_offset))
        print("  > Permutation leading zeros: {}".format(self.leading_zeros))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
