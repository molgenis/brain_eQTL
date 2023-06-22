#!/usr/bin/env python3

"""
File:         pre_process_TPM_expression_matrix.py
Created:      2021/11/23
Last Changed: 2022/10/18
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
import glob
import math
import time
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Local application imports.

# Metadata
__program__ = "Pre-process TPM Expression Matrix"
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
./pre_process_TPM_expression_matrix2.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-01-31-raw-count-tables/ \
    -ra /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis/create_correction_matrix/technical_covariates_table_top20.txt.gz \
    -gi /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/gencode.v32.primary_assembly.annotation-genelengths.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex_EUR.txt.gz \
    -o 2022-01-19-MetaBrain-CortexEUR-NegativeToZero-DatasetAndRAMCorrected
    
./pre_process_TPM_expression_matrix2.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-01-31-raw-count-tables/ \
    -ra /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/create_correction_matrix/technical_covariates_table_top20.txt.gz \
    -gi /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/gencode.v32.primary_assembly.annotation-genelengths.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex_AFR.txt.gz \
    -o 2022-02-09-MetaBrain-CortexAFR-NegativeToZero-DatasetAndRAMCorrected
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.rna_alignment_path = getattr(arguments, 'rna_alignment')
        self.gene_info_path = getattr(arguments, 'gene_info')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        outdir = str(Path(__file__).parent.parent)
        self.plot_outdir = os.path.join(outdir, 'pre_process_TPM_expression_matrix', outfolder, 'plot')
        self.file_outdir = os.path.join(outdir, 'pre_process_TPM_expression_matrix', outfolder, 'data')
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

        self.marker_genes = ['ENSG00000081237', 'ENSG00000102468', 'ENSG00000178031', 'ENSG00000170075', 'ENSG00000183230', 'ENSG00000187800', 'ENSG00000162745', 'ENSG00000110975', 'ENSG00000152583', 'ENSG00000135919', 'ENSG00000113140', 'ENSG00000163421', 'ENSG00000187094', 'ENSG00000198879', 'ENSG00000133392', 'ENSG00000172164', 'ENSG00000124126', 'ENSG00000143119', 'ENSG00000170962', 'ENSG00000127249', 'ENSG00000116147', 'ENSG00000198822', 'ENSG00000146648', 'ENSG00000131095', 'ENSG00000169756', 'ENSG00000163630', 'ENSG00000183960', 'ENSG00000148482', 'ENSG00000204262', 'ENSG00000175164', 'ENSG00000198121', 'ENSG00000150556', 'ENSG00000198846', 'ENSG00000149564', 'ENSG00000156966', 'ENSG00000206190', 'ENSG00000008394', 'ENSG00000168481', 'ENSG00000237289', 'ENSG00000103089', 'ENSG00000088992', 'ENSG00000188906', 'ENSG00000187398', 'ENSG00000153395', 'ENSG00000165072', 'ENSG00000101311', 'ENSG00000184985', 'ENSG00000080293', 'ENSG00000125430', 'ENSG00000181234', 'ENSG00000165959', 'ENSG00000104313', 'ENSG00000171246', 'ENSG00000179796', 'ENSG00000125850', 'ENSG00000120549', 'ENSG00000150893', 'ENSG00000167755', 'ENSG00000158859', 'ENSG00000168081', 'ENSG00000136099', 'ENSG00000054598', 'ENSG00000153266', 'ENSG00000095713', 'ENSG00000115956', 'ENSG00000144040', 'ENSG00000164600', 'ENSG00000177469', 'ENSG00000180287', 'ENSG00000196660', 'ENSG00000168772', 'ENSG00000100362', 'ENSG00000087116', 'ENSG00000143858', 'ENSG00000157890', 'ENSG00000164124', 'ENSG00000092969', 'ENSG00000117009', 'ENSG00000105695', 'ENSG00000204287', 'ENSG00000136160', 'ENSG00000159788', 'ENSG00000242808', 'ENSG00000169562', 'ENSG00000127252', 'ENSG00000147246', 'ENSG00000183779', 'ENSG00000070808', 'ENSG00000137812', 'ENSG00000132938', 'ENSG00000164794', 'ENSG00000115593', 'ENSG00000125675', 'ENSG00000158089', 'ENSG00000166863', 'ENSG00000182836', 'ENSG00000021645', 'ENSG00000103740', 'ENSG00000137252', 'ENSG00000069122', 'ENSG00000072315', 'ENSG00000168461', 'ENSG00000104611', 'ENSG00000122966', 'ENSG00000187957', 'ENSG00000116761', 'ENSG00000196358', 'ENSG00000064270', 'ENSG00000006128', 'ENSG00000103489', 'ENSG00000168314', 'ENSG00000125730', 'ENSG00000196569', 'ENSG00000166448', 'ENSG00000139155', 'ENSG00000180801', 'ENSG00000185518', 'ENSG00000162631', 'ENSG00000155897', 'ENSG00000244242', 'ENSG00000152380', 'ENSG00000163032', 'ENSG00000068078', 'ENSG00000124493', 'ENSG00000188641', 'ENSG00000070371', 'ENSG00000116132', 'ENSG00000169184', 'ENSG00000198959', 'ENSG00000148704', 'ENSG00000108370', 'ENSG00000234965', 'ENSG00000171532', 'ENSG00000011465', 'ENSG00000182463', 'ENSG00000172137', 'ENSG00000153234', 'ENSG00000119938', 'ENSG00000144908', 'ENSG00000157570', 'ENSG00000182902', 'ENSG00000178235', 'ENSG00000172005', 'ENSG00000162687', 'ENSG00000160219', 'ENSG00000100033', 'ENSG00000196353', 'ENSG00000135643', 'ENSG00000103154', 'ENSG00000115756', 'ENSG00000152661', 'ENSG00000077943', 'ENSG00000111863', 'ENSG00000214548', 'ENSG00000144355', 'ENSG00000163285', 'ENSG00000132669', 'ENSG00000178878', 'ENSG00000101203', 'ENSG00000189056', 'ENSG00000128656', 'ENSG00000146530', 'ENSG00000152127', 'ENSG00000144619', 'ENSG00000164741', 'ENSG00000160200', 'ENSG00000165973', 'ENSG00000151012', 'ENSG00000162511', 'ENSG00000171444', 'ENSG00000166436', 'ENSG00000143013', 'ENSG00000123838', 'ENSG00000015133', 'ENSG00000145555', 'ENSG00000049323', 'ENSG00000151229', 'ENSG00000148180', 'ENSG00000106018', 'ENSG00000102755', 'ENSG00000164188', 'ENSG00000135577', 'ENSG00000166710', 'ENSG00000151789', 'ENSG00000134215', 'ENSG00000080493', 'ENSG00000101384', 'ENSG00000136960', 'ENSG00000163492', 'ENSG00000112309', 'ENSG00000113100', 'ENSG00000132470', 'ENSG00000198771', 'ENSG00000165246', 'ENSG00000110723', 'ENSG00000125869', 'ENSG00000066405', 'ENSG00000175352', 'ENSG00000168453', 'ENSG00000174607', 'ENSG00000038945', 'ENSG00000152784', 'ENSG00000140836', 'ENSG00000084636', 'ENSG00000150361', 'ENSG00000109255', 'ENSG00000168243', 'ENSG00000123243', 'ENSG00000022355', 'ENSG00000124205', 'ENSG00000100276', 'ENSG00000145703', 'ENSG00000165949', 'ENSG00000133789', 'ENSG00000166573', 'ENSG00000146122', 'ENSG00000091513', 'ENSG00000136541', 'ENSG00000141469', 'ENSG00000197430', 'ENSG00000162390', 'ENSG00000105855', 'ENSG00000172987', 'ENSG00000150722', 'ENSG00000187416', 'ENSG00000069018', 'ENSG00000170775', 'ENSG00000128573', 'ENSG00000118513', 'ENSG00000170381', 'ENSG00000012817', 'ENSG00000137573', 'ENSG00000114374', 'ENSG00000240891', 'ENSG00000171759', 'ENSG00000163132', 'ENSG00000127920', 'ENSG00000115844', 'ENSG00000107242', 'ENSG00000123901', 'ENSG00000162542', 'ENSG00000089199', 'ENSG00000108375', 'ENSG00000138759', 'ENSG00000102445', 'ENSG00000142156', 'ENSG00000113327', 'ENSG00000172508', 'ENSG00000119042', 'ENSG00000135424', 'ENSG00000145708', 'ENSG00000117707', 'ENSG00000141338', 'ENSG00000184368', 'ENSG00000075651', 'ENSG00000134207', 'ENSG00000136235', 'ENSG00000130529', 'ENSG00000077063', 'ENSG00000143248', 'ENSG00000146592', 'ENSG00000106976', 'ENSG00000172673', 'ENSG00000063438', 'ENSG00000081923', 'ENSG00000124302', 'ENSG00000182578', 'ENSG00000188848', 'ENSG00000140848', 'ENSG00000127152', 'ENSG00000096696', 'ENSG00000180537', 'ENSG00000182348', 'ENSG00000145864', 'ENSG00000213949', 'ENSG00000135744', 'ENSG00000160781', 'ENSG00000086205', 'ENSG00000140678', 'ENSG00000123560', 'ENSG00000103184', 'ENSG00000171860', 'ENSG00000129538', 'ENSG00000124749', 'ENSG00000074966', 'ENSG00000069431', 'ENSG00000168497', 'ENSG00000117266', 'ENSG00000132965', 'ENSG00000007372', 'ENSG00000019582', 'ENSG00000146090', 'ENSG00000138741', 'ENSG00000101489', 'ENSG00000132164', 'ENSG00000182841', 'ENSG00000152894', 'ENSG00000105880', 'ENSG00000124731', 'ENSG00000122584', 'ENSG00000157404', 'ENSG00000164161', 'ENSG00000165124', 'ENSG00000173391', 'ENSG00000163110', 'ENSG00000196083', 'ENSG00000115461', 'ENSG00000136167', 'ENSG00000197410', 'ENSG00000124145', 'ENSG00000151892', 'ENSG00000013297', 'ENSG00000197299', 'ENSG00000141404', 'ENSG00000152953', 'ENSG00000170458', 'ENSG00000175899', 'ENSG00000136531', 'ENSG00000113532', 'ENSG00000174948', 'ENSG00000180071', 'ENSG00000113504', 'ENSG00000164946', 'ENSG00000087494', 'ENSG00000158292', 'ENSG00000076513', 'ENSG00000164199', 'ENSG00000163531', 'ENSG00000196104', 'ENSG00000157005', 'ENSG00000072163', 'ENSG00000104888', 'ENSG00000150594', 'ENSG00000146469', 'ENSG00000174600', 'ENSG00000113494', 'ENSG00000073737', 'ENSG00000021300', 'ENSG00000079335', 'ENSG00000071967', 'ENSG00000134853', 'ENSG00000060656', 'ENSG00000128242', 'ENSG00000139354', 'ENSG00000171951', 'ENSG00000112902', 'ENSG00000164120', 'ENSG00000144406', 'ENSG00000133110', 'ENSG00000203867', 'ENSG00000091128', 'ENSG00000176641', 'ENSG00000069188', 'ENSG00000138356', 'ENSG00000171004', 'ENSG00000122012', 'ENSG00000134343', 'ENSG00000038295', 'ENSG00000198626', 'ENSG00000125844', 'ENSG00000135750', 'ENSG00000146477', 'ENSG00000160191', 'ENSG00000163873', 'ENSG00000118308', 'ENSG00000019991', 'ENSG00000221869', 'ENSG00000188039', 'ENSG00000145248', 'ENSG00000160293', 'ENSG00000060718', 'ENSG00000141449', 'ENSG00000138696', 'ENSG00000121742', 'ENSG00000181031', 'ENSG00000163518', 'ENSG00000106328', 'ENSG00000150938', 'ENSG00000153179', 'ENSG00000148908', 'ENSG00000198838', 'ENSG00000118432', 'ENSG00000115252', 'ENSG00000072657', 'ENSG00000168237', 'ENSG00000171502', 'ENSG00000106714', 'ENSG00000197256', 'ENSG00000116299', 'ENSG00000185634', 'ENSG00000106852', 'ENSG00000099139', 'ENSG00000118194', 'ENSG00000183166', 'ENSG00000150275', 'ENSG00000189108', 'ENSG00000173786', 'ENSG00000112149', 'ENSG00000182601', 'ENSG00000197971', 'ENSG00000103196', 'ENSG00000135046', 'ENSG00000198682', 'ENSG00000162494', 'ENSG00000106829', 'ENSG00000157851', 'ENSG00000186998', 'ENSG00000140563', 'ENSG00000150656', 'ENSG00000107249', 'ENSG00000112139', 'ENSG00000002745', 'ENSG00000131773', 'ENSG00000112964', 'ENSG00000172554', 'ENSG00000140379', 'ENSG00000137809', 'ENSG00000106078']

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
        parser.add_argument("-gi",
                            "--gene_info",
                            type=str,
                            required=True,
                            help="The path to the gene info matrix.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

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

        print("Loading data.")
        counts_df_list = []
        stats_df_list = []
        found_samples = set()
        special_samples_trans_dict = {"AN11864_ba41-42-22": "AN11864_ba41.42.22",
                                      "UMB1376_ba41-42-22": "UMB1376_ba41.42.22"}
        for filepath in glob.glob(os.path.join(self.data_path, "*.txt.gz")):
            df = self.load_file(filepath, header=0, index_col=0)
            df.dropna(how="all", inplace=True)

            filename = os.path.basename(filepath)
            if "TargetALS" in filename:
                df.columns = [re.sub("-", "_", colname) for colname in df.columns]
                df.columns = [re.sub("\\.", "_", colname) for colname in df.columns]
                df.columns = [re.search(".*(HRA_[0-9]+)", colname).group(1) for colname in df.columns]
            elif "Braineac" in filename:
                df.columns = [re.search(".*(A653.*)", colname).group(1) for colname in df.columns]
            elif "GTEx" in filename:
                df.columns = [re.search("(.*)_.*", colname).group(1) for colname in df.columns]
            elif "NABEC" in filename:
                df.columns = [re.search(".*_(.*)", colname).group(1) for colname in df.columns]
            elif "ENA" in filename:
                df.columns = [re.search(".*_(.*)", colname).group(1) for colname in df.columns]
            elif "BrainGVEx" in filename:
                df.columns = [re.sub("_", "-", colname) for colname in df.columns]
                df.columns = [re.sub("\\.", "-", colname) for colname in df.columns]
                df.columns = [re.sub("^X", "", colname) for colname in df.columns]
                df.columns = [re.search("^([0-9]+-[0-9]+)-[0-9]+-[0-9]+", colname).group(1) for colname in df.columns]
            elif "BipSeq" in filename:
                df.columns = [re.search("Br[0-9]+_(R[0-9]+)", colname).group(1) for colname in df.columns]
            elif "UCLA_ASD" in filename:
                df.columns = [re.search("[aA-zZ]+[0-9]+_(.+)", colname).group(1) for colname in df.columns]
                df.columns = [special_samples_trans_dict[colname] if colname in special_samples_trans_dict else colname for colname in df.columns]
            elif "CMC_HBCC" in filename:
                df.columns = [re.search("individualID.*_specimenID.(.*)", colname).group(1) for colname in df.columns]
            elif "CMC" in filename:
                df.columns = [re.search("^CMC_[aA-zZ]+_[0-9]+_(.*)", colname).group(1) for colname in df.columns]
            elif "MSBB" in filename:
                df.columns = [re.sub(".accepted_hits.sort.coordReadsPerGene.out.tab", "", colname) for colname in df.columns]
                df.columns = [re.search("AMPAD_MSSM_[0-9]+_(.*)", colname).group(1) for colname in df.columns]
            elif "ROSMAP" in filename:
                df.columns = [re.sub("^X", "", colname) for colname in df.columns]
                df.columns = [re.sub("ReadsPerGene.out.tab", "", colname) for colname in df.columns]
                df.columns = [re.search(".*_.*_(.*_.*)", colname).group(1) for colname in df.columns]
            elif "MayoCBE" in filename:
                df.columns = [re.sub("^X", "", colname) for colname in df.columns]
                df.columns = [re.search("[0-9]+_CER", colname).group(0) for colname in df.columns]
            elif "MayoTCX" in filename:
                df.columns = [re.sub("^X", "", colname) for colname in df.columns]
                df.columns = [re.search("[0-9]+_TCX", colname).group(0) for colname in df.columns]
            elif "Brainseq" in filename:
                # these did not get adjusted in the other talbes, so keep the same
                pass
            else:
                print("Unexpected input file.")
                exit()

            found_samples.update(set(df.columns))

            if "ENA" in filename:
                # Ends with
                # __no_feature
                # __ambiguous
                # __too_low_aQual
                # __not_aligned
                # __alignment_not_unique
                stats_df = df.iloc[(df.shape[0] - 5):, :]
                stats_df.index = ["N_noFeature", "N_ambiguous", "N_too_low_aQual", "N_unmapped", "N_multimapping"]

                counts_df_list.append(df.iloc[:(df.shape[0] - 5), :])
                stats_df_list.append(stats_df)
            else:
                # Starts with
                # N_unmapped
                # N_multimapping
                # N_noFeature
                # N_ambiguous
                stats_df_list.append(df.iloc[:4, :])
                counts_df_list.append(df.iloc[4:, :])

        missing_samples = [sample for sample in samples if sample not in found_samples]
        print("\t  Missing MetaBrain samples [N={}]: {}".format(len(missing_samples), ", ".join(missing_samples)))

        counts_df = pd.concat(counts_df_list, axis=1)
        counts_df.fillna(0, inplace=True)
        print(counts_df)

        print("Step 1: sample selection.")
        print("\tUsing {}/{} samples.".format(len(samples), counts_df.shape[1]))
        counts_df = counts_df.loc[:, samples]

        print("Step 2: remove probes with zero variance.")
        mask = counts_df.std(axis=1) != 0
        print("\tUsing {}/{} probes.".format(np.sum(mask), np.size(mask)))
        counts_df = counts_df.loc[mask, :]

        print("Step 3: remove samples with zero counts / variance.")
        mask = (counts_df.std(axis=0) != 0) & (counts_df.sum(axis=0) != 0)
        print("\tUsing {}/{} samples.".format(np.sum(mask), np.size(mask)))
        counts_df = counts_df.loc[:, mask]

        print("\tSaving data")
        self.save_file(df=counts_df, outpath=os.path.join(self.file_outdir,"geneCounts.txt.gz"))
        # counts_df = self.load_file(os.path.join(self.file_outdir, "geneCounts.txt.gz"), header=0, index_col=0)

        print("Step 4: PCA analysis.")
        self.pca(df=counts_df,
                 filename="GeneCounts",
                 sample_to_dataset=sample_to_dataset,
                 plot_appendix="_1_geneCounts")

        print("Loading gene length data")
        gene_info_df = self.load_file(self.gene_info_path, header=0, index_col=0)
        missing_genes = [gene for gene in counts_df.index if gene not in gene_info_df.index]
        if len(missing_genes) > 0:
            print("Warning: missing gene info for {} genes".format(len(missing_genes)))
            print(missing_genes)

        # Subset genes for which we have info.
        gene_overlap = set(counts_df.index).intersection(set(gene_info_df.index))
        counts_df = counts_df.loc[gene_overlap, :]

        print("Step 5: Calculating TPM values")
        # https://btep.ccr.cancer.gov/question/faq/what-is-the-difference-between-rpkm-fpkm-and-tpm/
        # Divide the read counts by the length of each gene in kilobases. This
        # gives you reads per kilobase (RPK).
        kilo_bases_s = gene_info_df.loc[gene_overlap, "MergedExonLength"] / 1e3
        rpk_df = counts_df.divide(kilo_bases_s, axis=0)

        # Count up all the RPK values in a sample and divide this number by
        # 1,000,000. This is your “per million” scaling factor.
        pm_scaling_factor = rpk_df.sum(axis=0) / 1e6

        # Divide the RPK values by the “per million” scaling factor.
        # This gives you TPM.
        tpm_df = rpk_df.divide(pm_scaling_factor, axis=1)
        del counts_df, kilo_bases_s, rpk_df, pm_scaling_factor

        print("\tSaving data")
        self.save_file(df=tpm_df, outpath=os.path.join(self.file_outdir,"geneCounts.TPM.MergedExonLength.txt.gz"))
        # tpm_df = self.load_file(os.path.join(self.file_outdir,"geneCounts.TPM.MergedExonLength.txt.gz"), header=0, index_col=0)

        print("Step 6: PCA analysis.")
        self.pca(df=tpm_df,
                 filename="geneCounts.TPM.MergedExonLength",
                 sample_to_dataset=sample_to_dataset,
                 plot_appendix="_2_TPM")

        print("Step 7: log2 transform.")
        min_value = tpm_df.min(axis=1).min()
        if min_value <= 0:
            tpm_df = np.log2(tpm_df - min_value + 1)
        else:
            tpm_df = np.log2(tpm_df + 1)

        print("\tSaving data")
        self.save_file(df=tpm_df, outpath=os.path.join(self.file_outdir,"geneCounts.TPM.MergedExonLength.Log2Transformed.txt.gz"))
        # tpm_df = self.load_file(os.path.join(self.file_outdir, "geneCounts.TPM.MergedExonLength.Log2Transformed.txt.gz"), header=0, index_col=0)

        print("Step 8: save mean and std per gene.")
        mean = tpm_df.mean(axis=1)
        std = tpm_df.std(axis=1)

        print("Step 9: PCA analysis.")
        self.pca(df=tpm_df,
                 filename="geneCounts.TPM.MergedExonLength.Log2Transformed",
                 sample_to_dataset=sample_to_dataset,
                 plot_appendix="_3_TPM_Log2Transformed")

        print("Step 10: Construct correction matrix.")
        ram_df = None
        if self.rna_alignment_path is not None:
            ram_df = self.load_file(self.rna_alignment_path,
                                    header=0,
                                    index_col=0)
            ram_df = ram_df.loc[:, samples].T

        correction_df = self.prepare_correction_matrix(dataset_df=dataset_df,
                                                       ram_df=ram_df)

        # correction_df = self.prepare_correction_matrix2(ram_df=ram_df)

        print("\tSaving file.")
        self.save_file(df=correction_df, outpath=os.path.join(self.file_outdir, "correction_matrix.txt.gz"))

        self.plot_scatterplot(df=tpm_df.T,
                              sa_df=std_df,
                              columns=tpm_df.index[:5],
                              filename="_1_TPM_Log2Transformed")

        print("Step 11: remove technical covariates OLS.")
        tpm_corrected_df = self.calculate_residuals(df=tpm_df,
                                                    correction_df=correction_df)

        print("\tSaving data")
        self.save_file(df=tpm_corrected_df, outpath=os.path.join(self.file_outdir,"geneCounts.TPM.MergedExonLength.Log2Transformed.CovariatesRemovedOLS.txt.gz"))

        print("Step 12: PCA analysis.")
        self.pca(df=tpm_corrected_df,
                 filename="geneCounts.TPM.MergedExonLength.Log2Transformed.CovariatesRemovedOLS",
                 sample_to_dataset=sample_to_dataset,
                 plot_appendix="_4_TPM_Log2Transformed_CovariatesRemovedOLS")

        self.plot_scatterplot(df=tpm_corrected_df.T,
                              sa_df=std_df,
                              columns=tpm_df.index[:5],
                              filename="_2_TPM_Log2Transformed_CovariatesRemovedOLS")

        print("Step 13: return distribution shape and location.")
        tpm_corrected_df = tpm_corrected_df.subtract(tpm_corrected_df.mean(axis=1), axis=0).mul(std / tpm_corrected_df.std(axis=1), axis=0).add(mean, axis=0)

        self.plot_scatterplot(df=tpm_corrected_df.T,
                              sa_df=std_df,
                              columns=tpm_df.index[:5],
                              filename="_3_TPM_Log2Transformed_CovariatesRemovedOLS_ScaleAndLocReturned")

        print("\tSaving data")
        self.save_file(df=tpm_corrected_df, outpath=os.path.join(self.file_outdir,"geneCounts.TPM.MergedExonLength.Log2Transformed.CovariatesRemovedOLS.ScaleAndLocReturned.txt.gz"))

        print("Step 14: PCA analysis.")
        self.pca(df=tpm_corrected_df,
                 filename="geneCounts.TPM.MergedExonLength.Log2Transformed.CovariatesRemovedOLS.ScaleAndLocReturned",
                 sample_to_dataset=sample_to_dataset,
                 plot_appendix="_5_TPM_Log2Transformed_CovariatesRemovedOLS_ScaleAndColReturned")

        print("Step 15: Replace negative with zero.")
        min_value = tpm_corrected_df.values.min()
        if min_value < 0:
            print("\tLowest value before: {}".format(min_value))
            tpm_corrected_df[tpm_corrected_df < 0] = 0
            print("\tLowest value after: {}".format(tpm_corrected_df.values.min()))

        self.plot_scatterplot(df=tpm_corrected_df.T,
                              sa_df=std_df,
                              columns=tpm_df.index[:5],
                              filename="_4_TPM_Log2Transformed_CovariatesRemovedOLS_ScaleAndLocReturned_NegativeToZero")

        self.save_file(df=tpm_corrected_df, outpath=os.path.join(self.file_outdir, "geneCounts.TPM.MergedExonLength.Log2Transformed.CovariatesRemovedOLS.ScaleAndLocReturned.NegativeToZero.txt.gz"))
        # tpm_corrected_df = self.load_file(os.path.join(self.file_outdir, "geneCounts.TPM.MergedExonLength.Log2Transformed.CovariatesRemovedOLS.ScaleAndLocReturned.NegativeToZero.txt.gz"), header=0, index_col=0)

        print("Step 14: PCA analysis.")
        self.pca(df=tpm_corrected_df,
                 filename="geneCounts.TPM.MergedExonLength.Log2Transformed.CovariatesRemovedOLS.ScaleAndLocReturned.NegativeToZero",
                 sample_to_dataset=sample_to_dataset,
                 plot_appendix="_6_TPM_Log2Transformed_CovariatesRemovedOLS_ScaleAndColReturned_NegativeToZero")

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

    def prepare_correction_matrix(self, dataset_df, ram_df=None):
        correction_df = dataset_df.iloc[:, 1:]
        if ram_df is not None:
            # Remove columns without variance and filter the RNAseq alignment
            # metrics on VIF.
            ram_df_subset_df = ram_df.copy()
            ram_df_subset_df = self.remove_multicollinearity(ram_df_subset_df.loc[:, ram_df_subset_df.std(axis=0) != 0])

            # Merge with correction data frame.
            correction_df = correction_df.merge(ram_df_subset_df, left_index=True, right_index=True)

        # Add intercept.
        correction_df.insert(0, "INTERCEPT", 1)
        correction_df.index.name = "-"

        return correction_df

    def prepare_correction_matrix2(self, ram_df):
        # Remove columns without variance and filter the RNAseq alignment
        # metrics on VIF.
        correction_df = ram_df.copy()
        correction_df = self.remove_multicollinearity(correction_df.loc[:, correction_df.std(axis=0) != 0])

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

    def pca(self, df, filename, sample_to_dataset, plot_appendix=""):
        # samples should be on the columns and genes on the rows.
        zscores = (df - df.mean(axis=0)) / df.std(axis=0)
        pca = PCA(n_components=25)
        pca.fit(zscores)
        components_df = pd.DataFrame(pca.components_)
        components_df.index = ["Comp{}".format(i + 1) for i, _ in enumerate(components_df.index)]
        components_df.columns = df.columns

        print("\tSaving file.")
        self.save_file(df=components_df, outpath=os.path.join(self.file_outdir, "{}.PCAOverSamplesEigenvectors.txt.gz".format(filename)))

        print("\tPlotting PCA")
        plot_df = components_df.T
        plot_df["cohort"] = plot_df.index.map(sample_to_dataset)
        self.plot(df=plot_df, x="Comp1", y="Comp2", hue="cohort", palette=self.palette,
                  xlabel="PC1 [{:.2f}%]".format(pca.explained_variance_ratio_[0] * 100),
                  ylabel="PC2 [{:.2f}%]".format(pca.explained_variance_ratio_[1] * 100),
                  title="PCA - eigenvectors",
                  filename="eigenvectors_plot{}".format(plot_appendix))

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

    def plot_scatterplot(self, df, sa_df, columns=None, filename="plot"):
        std_df = sa_df.copy()
        std_df.index = std_df["rnaseq_id"]
        std_df = std_df[["dataset"]]
        if columns is None:
            columns = df.columns.tolist()

        nplots = len(columns) + 1
        ncols = math.ceil(np.sqrt(nplots))
        nrows = math.ceil(nplots / ncols)

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 # sharex='col',
                                 # sharey='row',
                                 figsize=(12 * ncols, 12 * nrows))
        sns.set(color_codes=True)

        row_index = 0
        col_index = 0
        groups_present = set()
        for i in range(ncols * nrows):
            if nrows == 1:
                ax = axes[col_index]
            elif ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            if i < len(columns):
                sns.despine(fig=fig, ax=ax)

                pi = columns[i]

                # Merge.
                plot_df = df[[pi]].copy()
                plot_df.columns = ["y"]
                plot_df.dropna(inplace=True)

                # Add color
                hue = None
                palette = None
                if std_df is not None:
                    hue = std_df.columns[0]
                    palette = self.palette
                    plot_df = plot_df.merge(std_df, left_index=True, right_index=True)

                # set order.
                plot_df["x"] = df[pi].argsort()
                counter = 0
                for group in plot_df["dataset"].unique():
                    mask = plot_df["dataset"] == group
                    subset = plot_df.loc[mask, :]
                    plot_df.loc[mask, "x"] = subset["y"].argsort() + counter
                    counter += np.sum(mask)
                    groups_present.add(group)

                # Plot.
                self.plot_scatterplot_panel(ax=ax,
                                            df=plot_df,
                                            hue=hue,
                                            palette=palette,
                                            title=pi)
            else:
                ax.set_axis_off()

                if std_df is not None and self.palette is not None and i == (nplots - 1):
                    handles = []
                    for key, value in self.palette.items():
                        if key in groups_present:
                            handles.append(mpatches.Patch(color=value, label=key))
                    ax.legend(handles=handles, loc=4, fontsize=25)

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        fig.savefig(os.path.join(self.plot_outdir, "scatterplot{}.png".format(filename)))
        plt.close()

    @staticmethod
    def plot_scatterplot_panel(ax, df, x="x", y="y", hue=None, palette=None,
                     xlabel="", ylabel="", title=""):

        # Plot.
        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        palette=palette,
                        linewidth=0,
                        legend=False,
                        ax=ax)

        ax.set_title(title,
                     fontsize=40,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=20,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=20,
                      fontweight='bold')

    def print_arguments(self):
        print("Arguments:")
        print("  > Data: {}".format(self.data_path))
        print("  > RNAseq alignment metrics: {}".format(self.rna_alignment_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Plot output directory: {}".format(self.plot_outdir))
        print("  > File output directory: {}".format(self.file_outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
