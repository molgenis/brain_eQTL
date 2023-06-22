#!/usr/bin/env python3

"""
File:         partial_deconvolution.py
Created:      2020/06/29
Last Changed: 2020/11/22
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

# Third party imports.

# Local application imports.
from partial_deconvolution.src.main import Main
from partial_deconvolution.src.settings import Settings
from partial_deconvolution.src.cmd_line_arguments import CommandLineArguments

# Metadata
__program__ = "Partial Deconvolution"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
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
./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/CellMap_brain_celltype_avgCPM.txt -t /groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o OLD_PROFILE_CORTEX_EUR_TMM_LOG2 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -os IHC_0CPM_LOG2_FILTERED_CC -dataset_correction -log2 -sum_to_one -visualise -e pdf

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/CellMap_brain_CNS7_avgCPM.txt -t /groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o NEW_PROFILE_CORTEX_EUR_TMM_LOG2 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -os IHC_0CPM_LOG2_FILTERED_CC -dataset_correction -log2 -sum_to_one -visualise -e png

### PSYCHENCODE ###

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/CellMap_brain_CNS7_NoPericytes_avgCPM.txt -t /groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o NEW_PROFILE_NOPERICYTES_CORTEX_EUR_TMM_LOG2 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -os IHC_0CPM_LOG2_FILTERED -log2 -sum_to_one -visualise -e png

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/CellMap_brain_CNS7_NoPericytes_PsychENCODEMarkerGenes_avgCPM.txt -t /groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o NEW_PROFILE_PSYCHENCODE_MARKER_GENES_NOPERICYTES_CORTEX_EUR_TMM_LOG2 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -os IHC_0CPM_LOG2_FILTERED -log2 -sum_to_one -visualise -e png

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE_ref_profile_unique_avg.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.PsychENCODEGenesExtended.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o PSYCHENCODE_PROFILE_CORTEX_EUR_TMM_LOG2 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -os IHC_0CPM_LOG2_FILTERED_CC -log2 -sum_to_one -visualise -e png


# NOVEL 22-11-2021
./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/DER-20_Single_cell_expression_processed_TPM_MGFIltered_AvgPerCT.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o PSYCHENCODE_PROFILE_CORTEX_EUR_TMM_LOG2 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -os IHC_0CPM_LOG2 -log2 -visualise -e png

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrain/data/geneCounts.TPM.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/DER-20_Single_cell_expression_processed_TPM_MGFIltered_AvgPerCT.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o PSYCHENCODE_PROFILE_TPM_LOG2 -os IHC_0CPM_LOG2 -log2 -visualise -e png

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrain/data/geneCounts.TPM.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/DER-20_Single_cell_expression_processed_TPM_MGFIltered_AvgPerCT.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/single_cell_counts_all_celltypes.txt.gz -o PSYCHENCODE_PROFILE_CORTEX_EUR_TPM_LOG2 -os SC_0CPM_LOG2 -log2 -visualise -e png

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrainAndPsychENCODE/data/geneCounts.TPM.MergedExonLength.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/DER-20_Single_cell_expression_processed_TPM_MGFIltered_AvgPerCT.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o PSYCHENCODE_PROFILE_METABRAIN_AND_PSYCHENCODE_EXON_TPM_LOG2 -os IHC_0CPM_LOG2 -log2 -visualise -e png

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrainAndPsychENCODE/data/geneCounts.TPM.MergedExonLength.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/DER-20_Single_cell_expression_processed_TPM_MGFIltered_AvgPerCT_NoDev.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o PSYCHENCODE_PROFILE_METABRAIN_AND_PSYCHENCODE_EXON_TPM_LOG2_NODEV -os IHC_0CPM_LOG2 -log2 -visualise -e png


### PSYCHENCODE ###

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/CellMap_brain_CNS7_avgCPM.txt -t /groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/single_cell_counts_all_celltypes.txt.gz -o NEW_PROFILE_CORTEX_EUR_TMM_LOG2 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -os SC_0CPM_LOG2_FILTERED_CC -dataset_correction -log2 -sum_to_one -visualise -e png

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/CellMap_brain_CNS7_NoPericytes_avgCPM.txt -t /groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/single_cell_counts_all_celltypes.txt.gz -o NEW_PROFILE_NOPERICYTES_CORTEX_EUR_TMM_LOG2 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -os SC_0CPM_LOG2_FILTERED_CC -dataset_correction -log2 -sum_to_one -visualise -e png

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/CellMap_brain_CNS7_NoPericytes_avgCPM.txt -t /groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/single_cell_counts_all_celltypes.txt.gz -o NEW_PROFILE_NOPERICYTES_CORTEX_AFR_TMM_LOG2 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/combine_gte_files/SampleToDataset.txt.gz -os SC_0CPM_LOG2_FILTERED_CC -dataset_correction -log2 -sum_to_one -visualise -e png

"""


if __name__ == '__main__':
    # Get the command line arguments.
    CLA = CommandLineArguments(program=__program__,
                               version=__version__,
                               description=__description__)
    SETTINGS = Settings(data_path=CLA.get_argument("data"),
                        signature_path=CLA.get_argument("signature"),
                        translate_path=CLA.get_argument("translate"),
                        ground_truth_path=CLA.get_argument("ground_truth"),
                        sample_to_dataset_path=CLA.get_argument("sample_to_dataset"),
                        sample_filter_path=CLA.get_argument("sample_filter"),
                        dataset_filter=CLA.get_argument("dataset_filter"),
                        min_expr=CLA.get_argument("min_expr"),
                        dataset_correction=CLA.get_argument("dataset_correction"),
                        normalize=CLA.get_argument("normalize"),
                        zscore=CLA.get_argument("zscore"),
                        log2=CLA.get_argument("log2"),
                        decon_method=CLA.get_argument("decon_method"),
                        extension=CLA.get_argument("extension")
                        )

    # Start the program.
    PROGRAM = Main(settings=SETTINGS,
                   outdir=CLA.get_argument("outdir"),
                   outsubdir=CLA.get_argument("outsubdir"),
                   visualise=CLA.get_argument("visualise"),
                   plot_ids=CLA.get_argument("plot_id"),)
    PROGRAM.start()
