#!/usr/bin/env python3

"""
File:         visualise_ct_mediated_eqtl_replication.py
Created:      2022/02/16
Last Changed: 2022/02/21
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
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Visualise CellType mediated eQTL Replication"
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
./visualise_ct_mediated_eqtl_replication.py \
    -eq ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined_withFDRCol.txt.gz \
    -ge ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -al ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -ex ../preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset_CF4SDFiltered.txt.gz \
    -d ../decon-eqtl_scripts/decon_eqtl/2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -rr sn_replication/cis/2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected/single_nucleus_replication.txt.gz \
    -br bryois_replication/bryois_replication.txt.gz \
    -i ENSG00000188732.11_7:23681366:rs4722244:C_T_Oligodendrocyte ENSG00000084628.10_1:31208167:rs7549197:T_C_Oligodendrocyte ENSG00000015592.16_8:27245507:rs17366947:A_G_Oligodendrocyte ENSG00000133805.15_11:10452342:rs11042811:C_T_Oligodendrocyte ENSG00000085117.12_11:44614997:rs7935166:G_A_Oligodendrocyte ENSG00000004468.13_4:15735725:rs4698412:G_A_Astrocyte  \
    -n 1350 \
    -e png pdf
    
./visualise_ct_mediated_eqtl_replication.py \
    -eq ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined_withFDRCol.txt.gz \
    -ge ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -al ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -ex ../preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std ../matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset_CF4SDFiltered.txt.gz \
    -d ../decon-eqtl_scripts/decon_eqtl/2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -rr sn_replication/cis/2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected/single_nucleus_replication.txt.gz \
    -br bryois_replication/bryois_replication.txt.gz \
    -i ENSG00000184293.7_12:9724600:rs7306304:G_A_Microglia ENSG00000203710.11_1:207577223:rs679515:T_C_Oligodendrocyte ENSG00000004468.13_4:15735725:rs4698412:G_A_Astrocyte ENSG00000078487.17_7:100419831:rs7783159:G_A_Oligodendrocyte \
    -n 1350 \
    -e png pdf
    
./visualise_ct_mediated_eqtl_replication.py \
    -eq ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ge ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -al ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -ex ../preprocess_scripts/select_and_reorder_matrix/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.100PCsRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -d ../decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -rr sn_replication/trans/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/single_nucleus_replication.txt.gz \
    -i ENSG00000110042.8_7:64546850:rs73130769:C_G_Oligodendrocyte ENSG00000110042.8_7:64508517:rs57196191:A_C_Oligodendrocyte ENSG00000110042.8_7:64480828:rs73127029:G_A_Oligodendrocyte ENSG00000110042.8_7:64563203:rs2862795:A_G_Oligodendrocyte \
    -e png pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.geno_path = getattr(arguments, 'genotype')
        self.alleles_path = getattr(arguments, 'alleles')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cellcount%')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.decon_path = getattr(arguments, 'decon')
        self.rosmap_replication_path = getattr(arguments, 'rosmap_replication')
        self.bryois_replication_path = getattr(arguments, 'bryois_replication')
        self.interest = getattr(arguments, 'interest')
        self.nrows = getattr(arguments, 'nrows')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'visualise_ct_mediated_eqtl_replication')

        self.palette = {
            2.: "#E69F00",
            1.: "#0072B2",
            0.: "#D55E00",
            "AST": "#D55E00",
            "END": "#CC79A7",
            "EX": "#56B4E9",
            "IN": "#0072B2",
            "MIC": "#E69F00",
            "OPC": "#009E73",
            "OLI": "#009E73",
            "PER": "#808080"
        }

        self.ct_abbrevations = {
            "Astrocyte": "AST",
            "EndothelialCell": "END",
            "Excitatory": "EX",
            "Inhibitory": "IN",
            "Microglia": "MIC",
            "OPC": "OPC",
            "OPCsCOPs": "OPC",
            "Oligodendrocyte": "OLI",
            "Pericytes": "PER",
            "Pericyte": "PER"
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
                            help="show program's version number and exit")
        parser.add_argument("-eq",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the eqtl matrix")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix")
        parser.add_argument("-al",
                            "--alleles",
                            type=str,
                            required=True,
                            help="The path to the alleles matrix")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix")
        parser.add_argument("-cc",
                            "--cellcount%",
                            type=str,
                            required=True,
                            help="The path to the cell count % matrix")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=False,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-d",
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix")
        parser.add_argument("-rr",
                            "--rosmap_replication",
                            type=str,
                            required=True,
                            help="The path to the ROSMAP replication "
                                 "matrix")
        parser.add_argument("-br",
                            "--bryois_replication",
                            type=str,
                            default=None,
                            help="The path to the Bryois replication "
                                 "matrix")
        parser.add_argument("-i",
                            "--interest",
                            nargs="+",
                            type=str,
                            required=True,
                            help="The IDs to plot.")
        parser.add_argument("-n",
                            "--nrows",
                            type=int,
                            required=False,
                            default=None,
                            help="Cap the number of runs to load. "
                                 "Default: None.")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        eqtl_df = self.load_file(self.eqtl_path, header=0, index_col=None, nrows=self.nrows)
        geno_df = self.load_file(self.geno_path, header=0, index_col=0, nrows=self.nrows)
        alleles_df = self.load_file(self.alleles_path, header=0, index_col=0, nrows=self.nrows)
        expr_df = self.load_file(self.expr_path, header=0, index_col=0, nrows=self.nrows)
        cc_df = self.load_file(self.cc_path, header=0, index_col=0)
        decon_df = self.load_file(self.decon_path, header=0, index_col=None, nrows=self.nrows)
        rosmap_df = self.load_file(self.rosmap_replication_path, header=0, index_col=None)

        bryois_df = None
        if self.bryois_replication_path is not None:
            bryois_df = self.load_file(self.bryois_replication_path, header=0, index_col=None)

        if self.std_path:
            std_df = self.load_file(self.std_path, header=0, index_col=None)

            print("Filter data")
            samples = list(std_df.iloc[:, 0])
            geno_df = geno_df.loc[:, samples]
            expr_df = expr_df.loc[:, samples]
            cc_df = cc_df.loc[samples, :]

        print("Validate data")
        probes = list(eqtl_df["ProbeName"])
        snps = list(eqtl_df["SNPName"])
        samples = list(expr_df.columns)
        if list(geno_df.index) != snps:
            print("Error, genotype does not match eQTL file.")
            exit()
        if list(geno_df.columns) != samples:
            print("Error, genotype does not match expression file.")
            exit()
        if list(alleles_df.index) != snps:
            print("Error, allele does not match eQTL file.")
            exit()
        if list(expr_df.index) != probes:
            print("Error, expression does not match eQTL file.")
            exit()
        if list(cc_df.index) != samples:
            print("Error, cell fractions does not match expression file.")
            exit()

        print("Pre-process data")
        decon_df.index = decon_df["Gene"] + "_" + decon_df["SNP"]
        rosmap_df.index = rosmap_df["Gene"] + "_" + rosmap_df["SNP"]
        print(decon_df)

        rosmap_ct = [col.replace("SN ", "").replace(" eQTL pvalue", "") for col in rosmap_df.columns if col.startswith("SN ") and col.endswith(" eQTL pvalue")]
        rosmap_ct.sort()

        bryois_ct = []
        if bryois_df is not None:
            bryois_df.index = bryois_df["Gene"] + "_" + bryois_df["SNP"]
            bryois_ct = [col.replace("Bryois ", "").replace(" pvalue", "") for
                         col in bryois_df.columns if
                         col.startswith("Bryois ") and col.endswith(" pvalue")]
            bryois_ct.sort()

        print("Visualizing data")
        nrows = len(self.interest)
        ncols = 4 if bryois_df is not None else 3

        sns.set(rc={'figure.figsize': (ncols * 12, nrows * 9)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='none',
                                 sharey='none')

        bulk_shared_ylim = [0, 0]
        ieqtl_shared_xlim = [0, 0]
        sn_shared_xlim = [0, 0]
        bryois_shared_xlim = [0, 0]

        row_count = 0
        for row_index, (_, row) in enumerate(eqtl_df.iterrows()):
            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            id = probe_name + "_" + snp_name

            interest_id = None
            for interest in self.interest:
                if interest.startswith(id):
                    interest_id = interest

            if interest_id is None:
                continue

            print("\tWorking on: {}\t{}\t{} [{}/{} "
                  "{:.2f}%]".format(snp_name, probe_name, hgnc_name,
                                    row_index + 1,
                                    eqtl_df.shape[0],
                                    (100 / eqtl_df.shape[0]) * (row_index + 1)))

            # Get the genotype / expression data.
            genotype = geno_df.iloc[[row_index], :].copy().T
            if genotype.columns != [snp_name]:
                print("\t\tGenotype file not in identical order as eQTL file.")
                exit()
            expression = expr_df.iloc[[row_index], :].copy().T
            if expression.columns != [probe_name]:
                print("\t\tExpression file not in identical order as eQTL file.")
                exit()
            data = genotype.merge(expression, left_index=True, right_index=True)
            data.columns = ["genotype", "expression"]
            data["group"] = data["genotype"].round(0)

            # Remove missing values.
            data = data.loc[(data['genotype'] >= 0.0) &
                            (data['genotype'] <= 2.0), :]

            # Get the allele data.
            alleles = alleles_df.iloc[row_index, :]
            if alleles.name != snp_name:
                print("\t\tAlleles file not in identical order as eQTL file.")
                exit()
            # A/T = 0.0/2.0
            # by default we assume T = 2.0 to be minor
            major_allele = alleles["Alleles"].split("/")[0]
            minor_allele = alleles["Alleles"].split("/")[1]

            # Check if we need to flip the genotypes.
            counts = data["group"].value_counts()
            for x in [0.0, 1.0, 2.0]:
                if x not in counts:
                    counts.loc[x] = 0
            zero_geno_count = (counts[0.0] * 2) + counts[1.0]
            two_geno_count = (counts[2.0] * 2) + counts[1.0]
            if two_geno_count > zero_geno_count:
                # Turns out that 0.0 was the minor.
                minor_allele = alleles["Alleles"].split("/")[0]
                major_allele = alleles["Alleles"].split("/")[1]
                data["genotype"] = 2.0 - data["genotype"]
                data["group"] = 2.0 - data["group"]

            allele_map = {0.0: "{}/{}".format(major_allele, major_allele),
                          1.0: "{}/{}".format(major_allele, minor_allele),
                          2.0: "{}/{}".format(minor_allele, minor_allele)}

            # Determine annotation stats.
            eqtl_pvalue = np.nan
            if "MetaP" in row:
                eqtl_pvalue = row["MetaP"]
            elif "PValue" in row:
                eqtl_pvalue = row["PValue"]
            eqtl_pvalue_str = "{:.2e}".format(eqtl_pvalue)
            if eqtl_pvalue == 0:
                eqtl_pvalue_str = "<{:.1e}".format(1e-308)
            eqtl_pearsonr, _ = stats.pearsonr(data["expression"], data["genotype"])
            minor_allele_frequency = min(zero_geno_count, two_geno_count) / (zero_geno_count + two_geno_count)

            # Fill the eQTL plot annotation.
            annot1 = ["N: {:,}".format(data.shape[0]),
                      "r: {:.2f}".format(eqtl_pearsonr),
                      "p-value: {}".format(eqtl_pvalue_str),
                      "MAF: {:.2f}".format(minor_allele_frequency)]

            # Plot the main eQTL effect.
            ylim = self.eqtl_plot(
                fig=fig,
                ax=axes[row_count, 0],
                df=data,
                x="group",
                y="expression",
                palette=self.palette,
                allele_map=allele_map,
                xlabel=snp_name,
                ylabel="{} expression".format(hgnc_name),
                annot=annot1,
                title="Bulk eQTL effect" if row_count == 0 else ""
            )

            if ylim[0] < bulk_shared_ylim[0]:
                bulk_shared_ylim[0] = ylim[0]
            if ylim[1] > bulk_shared_ylim[1]:
                bulk_shared_ylim[1] = ylim[1]

            # Add the cell type of interest
            cell_type = interest_id.replace(id + "_", "")
            cf = cc_df.loc[:, [cell_type]]
            cf.columns = ["context"]
            data = data.merge(cf, left_index=True, right_index=True)

            # Determine annotation stats.
            interaction_pvalue = decon_df.loc[id, "{} pvalue".format(cell_type)]
            interaction_pvalue_str = "{:.2e}".format(interaction_pvalue)
            if interaction_pvalue == 0:
                interaction_pvalue_str = "<{:.1e}".format(1e-308)

            # Fill the interaction plot annotation.
            annot2 = ["{} - {}".format(snp_name.split(":")[2], minor_allele),
                      "Interaction p-value: {}".format(interaction_pvalue_str),
                      "eQTL p-value: {}".format(eqtl_pvalue_str),
                      "MAF: {:.2f}".format(minor_allele_frequency)]

            # Plot the interaction eQTL.
            xlim, ylim = self.inter_plot(
                fig=fig,
                ax=axes[row_count, 1],
                df=data,
                x="context",
                y="expression",
                group="group",
                palette=self.palette,
                allele_map=allele_map,
                xlabel="{} proportion".format(cell_type),
                ylabel="",
                annot=annot2,
                title="Bulk ieQTL effect" if row_count == 0 else ""
            )

            if xlim[0] < ieqtl_shared_xlim[0]:
                ieqtl_shared_xlim[0] = xlim[0]
            if xlim[1] > ieqtl_shared_xlim[1]:
                ieqtl_shared_xlim[1] = xlim[1]
            if ylim[0] < bulk_shared_ylim[0]:
                bulk_shared_ylim[0] = ylim[0]
            if ylim[1] > bulk_shared_ylim[1]:
                bulk_shared_ylim[1] = ylim[1]

            # Plot the replication in ROSMAP.
            if id in rosmap_df.index:
                rosmap_forest_df = pd.DataFrame(np.nan, index=rosmap_ct, columns=["cell type", "affect allele", "x", "lower", "upper"])
                for ct in rosmap_ct:
                    beta = float(rosmap_df.loc[id, "SN {} eQTL beta".format(ct)])
                    se = float(rosmap_df.loc[id, "SN {} eQTL se".format(ct)])
                    if rosmap_df.loc[id, "Allele assessed"] != minor_allele:
                        beta = beta * -1
                    rosmap_forest_df.loc[ct, :] = [self.ct_abbrevations[ct], minor_allele, beta, beta - se, beta + se]

                xlim = self.stripplot(
                    fig=fig,
                    ax=axes[row_count, 2],
                    df=rosmap_forest_df,
                    x="x",
                    y="cell type",
                    palette=self.palette,
                    xlabel="eQTL beta" if row_count == (nrows - 1) else "",
                    title="ROSMAP SN eQTL effect" if row_count == 0 else ""
                )

                if xlim[0] < sn_shared_xlim[0]:
                    sn_shared_xlim[0] = xlim[0]
                if xlim[1] > sn_shared_xlim[1]:
                    sn_shared_xlim[1] = xlim[1]

                del rosmap_forest_df
            else:
                axes[row_count, 2].set_axis_off()

            # Plot the replication in ROSMAP.
            if bryois_df is not None:
                if id in bryois_df.index:
                    bryois_forest_df = pd.DataFrame(np.nan, index=bryois_ct, columns=["cell type", "affect allele", "x", "lower", "upper"])
                    for ct in bryois_ct:
                        beta = float(bryois_df.loc[id, "Bryois {} eQTL beta".format(ct)])
                        pvalue = float(bryois_df.loc[id, "Bryois {} pvalue".format(ct)])
                        maf = float(bryois_df.loc[id, "MetaBrain MAF"])
                        n = float(bryois_df.loc[id, "Bryois N"])

                        _, se = self.pvalue_to_beta_and_se(beta=beta,
                                                           pvalue=pvalue,
                                                           maf=maf,
                                                           n=n)
                        if bryois_df.loc[id, "Allele assessed"] != minor_allele:
                            beta = beta * -1
                        bryois_forest_df.loc[ct, :] = [self.ct_abbrevations[ct], minor_allele, beta, beta - se, beta + se]

                    xlim = self.stripplot(
                        fig=fig,
                        ax=axes[row_count, 3],
                        df=bryois_forest_df,
                        x="x",
                        y="cell type",
                        palette=self.palette,
                        xlabel="eQTL beta" if row_count == (nrows - 1) else "",
                        title="Bryois 2021 eQTL effect" if row_count == 0 else ""

                    )

                    if xlim[0] < bryois_shared_xlim[0]:
                        bryois_shared_xlim[0] = xlim[0]
                    if xlim[1] > bryois_shared_xlim[1]:
                        bryois_shared_xlim[1] = xlim[1]

                    del bryois_forest_df
                else:
                    axes[row_count, 3].set_axis_off()

            row_count += 1

        for row_count in range(nrows):
            axes[row_count, 0].set_ylim(bulk_shared_ylim[0], bulk_shared_ylim[1])
            axes[row_count, 1].set_xlim(ieqtl_shared_xlim[0], ieqtl_shared_xlim[1])
            axes[row_count, 1].set_ylim(bulk_shared_ylim[0], bulk_shared_ylim[1])
            axes[row_count, 2].set_xlim(sn_shared_xlim[0] - 0.05, sn_shared_xlim[1] + 0.05)
            if bryois_df is not None:
                axes[row_count, 3].set_xlim(bryois_shared_xlim[0] - 0.05, bryois_shared_xlim[1] + 0.05)

            axes[row_count, 1].set_yticks([])
            if row_count < (nrows - 1):
                axes[row_count, 1].set_xticks([])
                axes[row_count, 2].set_xticks([])
                if bryois_df is not None:
                    axes[row_count, 3].set_xticks([])

        for extension in self.extensions:
            outpath = os.path.join(self.outdir, "visualise_ct_mediated_eqtl_replication.{}".format(extension))
            print("Saving plot: {}".format(os.path.basename(outpath)))
            fig.savefig(outpath)
        plt.close()

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
    def eqtl_plot(fig, ax, df, x="x", y="y", palette=None, allele_map=None,
                  annot=None, xlabel="", ylabel="", title=""):
        sns.despine(fig=fig, ax=ax)

        sns.regplot(x=x,
                    y=y,
                    data=df,
                    scatter=False,
                    ci=None,
                    line_kws={"color": "#000000"},
                    ax=ax)
        sns.violinplot(x=x,
                       y=y,
                       data=df,
                       palette=palette,
                       cut=0,
                       zorder=-1,
                       ax=ax)
        plt.setp(ax.collections, alpha=.75)
        sns.boxplot(x=x,
                    y=y,
                    data=df,
                    whis=np.inf,
                    color="white",
                    zorder=-1,
                    ax=ax)

        if annot is not None:
            for i, annot_label in enumerate(annot):
                ax.annotate(annot_label,
                            xy=(0.03, 0.94 - (i * 0.04)),
                            xycoords=ax.transAxes,
                            color="#000000",
                            alpha=0.75,
                            fontsize=12,
                            fontweight='bold')

        if allele_map is not None:
            ax.set_xticks(range(3))
            ax.set_xticklabels([allele_map[0.0], allele_map[1.0], allele_map[2.0]])

        ax.set_title(title,
                     fontsize=16,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        return ax.get_ylim()

    @staticmethod
    def inter_plot(fig, ax, df, x="x", y="y", group="group", palette=None,
                   allele_map=None, annot=None, xlabel="", ylabel="",
                   title=""):
        if len(set(df[group].unique()).symmetric_difference({0, 1, 2})) > 0:
            return

        sns.despine(fig=fig, ax=ax)

        for i, group_id in enumerate([0, 1, 2]):
            subset = df.loc[df[group] == group_id, :].copy()
            allele = group_id
            if allele_map is not None:
                allele = allele_map[group_id]

            coef_str = "NA"
            r_annot_pos = (-1, -1)
            if len(subset.index) > 1:
                coef, p = stats.pearsonr(subset[y], subset[x])
                coef_str = "{:.2f}".format(coef)

                subset["intercept"] = 1
                betas = np.linalg.inv(subset[["intercept", x]].T.dot(subset[["intercept", x]])).dot(subset[["intercept", x]].T).dot(subset[y])
                subset["y_hat"] = np.dot(subset[["intercept", x]], betas)
                subset.sort_values(x, inplace=True)

                r_annot_pos = (subset.iloc[-1, :][x] + (subset[x].max() * 0.05), subset.iloc[-1, :]["y_hat"])

                sns.regplot(x=x, y=y, data=subset, ci=None,
                            scatter_kws={'facecolors': palette[group_id],
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": palette[group_id], "alpha": 0.75},
                            ax=ax
                            )

            ax.annotate(
                '{}\n{}'.format(allele, coef_str),
                xy=r_annot_pos,
                color=palette[group_id],
                alpha=0.75,
                fontsize=16,
                fontweight='bold')

        if annot is not None:
            for i, annot_label in enumerate(annot):
                ax.annotate(annot_label,
                            xy=(0.03, 0.94 - (i * 0.04)),
                            xycoords=ax.transAxes,
                            color="#000000",
                            alpha=0.75,
                            fontsize=12,
                            fontweight='bold')

        (xmin, xmax) = (df[x].min(), df[x].max())
        (ymin, ymax) = (df[y].min(), df[y].max())
        xmargin = (xmax - xmin) * 0.05
        ymargin = (ymax - ymin) * 0.05

        ax.set_xlim(xmin - xmargin, xmax + xmargin)
        ax.set_ylim(ymin - ymargin, ymax + ymargin)

        ax.set_title(title,
                     fontsize=16,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        return ax.get_xlim(), ax.get_ylim()

    @staticmethod
    def stripplot(fig, ax, df, x="x", y="y", lower="lower", upper="upper",
                  palette=None, xlabel="", ylabel="", title=""):
        sns.despine(fig=fig, ax=ax)

        ax.axvline(-1, ls='--', color="#000000", alpha=0.15, zorder=-1)
        ax.axvline(-0.5, ls='--', color="#000000", alpha=0.15, zorder=-1)
        ax.axvline(0, ls='--', color="#000000", alpha=0.3, zorder=-1, linewidth=3)
        ax.axvline(0.5, ls='--', color="#000000", alpha=0.15, zorder=-1)
        ax.axvline(1, ls='--', color="#000000", alpha=0.15, zorder=-1)

        df_m = df.melt(id_vars=[y], value_vars=[lower, upper])
        sns.pointplot(x="value",
                      y=y,
                      data=df_m,
                      linewidth=4,
                      join=False,
                      palette=palette,
                      ax=ax)

        sns.stripplot(x=x,
                      y=y,
                      data=df,
                      size=25,
                      dodge=False,
                      orient="h",
                      palette=palette,
                      linewidth=0,
                      edgecolor="w",
                      jitter=0,
                      ax=ax)

        ax.set_title(title,
                     fontsize=16,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        return ax.get_xlim()

    @staticmethod
    def pvalue_to_beta_and_se(beta, pvalue, maf, n):
        zscore = stats.norm.ppf(pvalue / 2)
        if beta > 0:
            zscore = zscore * -1

        if pvalue == 1:
            zscore = 4
        if pvalue == 0:
            zscore = -40.

        chi = zscore * zscore
        a = 2 * maf * (1 - maf) * (n + chi)
        beta = zscore / a ** (1/2)
        se = 1 / a ** (1/2)

        return beta, se

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Alleles path: {}".format(self.alleles_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell count % path: {}".format(self.cc_path))
        print("  > Sample-to-dataset file: {}".format(self.std_path))
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > ROSMAP replication path: {}".format(self.rosmap_replication_path))
        print("  > Bryois replication path: {}".format(self.bryois_replication_path))
        print("  > Interest: {}".format(self.interest))
        print("  > Nrows: {}".format(self.nrows))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
