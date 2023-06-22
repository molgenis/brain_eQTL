#!/usr/bin/env python3

"""
File:         visualise_interaction_eqtl.py
Created:      2022/04/07
Last Changed: 2022/04/13
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
import matplotlib.patches as mpatches
from scipy import stats
from statsmodels.regression.linear_model import OLS

# Local application imports.

# Metadata
__program__ = "Visualise Interaction eQTL"
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
./visualise_interaction_eqtl.py -h

./visualise_interaction_eqtl.py \
    -eq ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ge ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -al ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -ex ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -co ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/create_covs_matrix/ADstatus.txt.gz \
    -std ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -i Inhibitory_8:115627247:rs2737226:T_C_Alzheimer_disease \
    -n 332 \
    -e png
    
./visualise_interaction_eqtl.py \
    -eq ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ge ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -al ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -ex ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/expression_table.txt.gz \
    -co ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_covs_matrix/ADstatus.txt.gz \
    -std ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -i ENSG00000172137.19_7:12244161:rs1990622:A_G_Alzheimer_disease \
    -n 225 \
    -e png
    
./visualise_interaction_eqtl.py \
    -eq ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_eqtlprobes/eQTLprobes_combined.txt.gz \
    -ge ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -al ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -ex ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/expression_table.txt.gz \
    -co ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_covs_matrix/Age_ST10Removed.txt.gz \
    -std ../matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -i ENSG00000172137.19_7:12244161:rs1990622:A_G_Age \
    -n 225 \
    -e png
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.geno_path = getattr(arguments, 'genotype')
        self.alleles_path = getattr(arguments, 'alleles')
        self.expr_path = getattr(arguments, 'expression')
        self.cova_path = getattr(arguments, 'covariate')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.interest = getattr(arguments, 'interest')
        self.nrows = getattr(arguments, 'nrows')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'visualise_interaction_eqtl')

        self.group_palette = {
            2.: "#E69F00",
            1.: "#0072B2",
            0.: "#D55E00"
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
        parser.add_argument("-co",
                            "--covariate",
                            type=str,
                            required=True,
                            help="The path to the covariate matrix.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
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
        cova_df = self.load_file(self.cova_path, header=0, index_col=0)
        std_df = self.load_file(self.std_path, header=0, index_col=None)

        expr_df = expr_df.T

        print("Validate data")
        probes = list(eqtl_df["ProbeName"])
        snps = list(eqtl_df["SNPName"])
        samples = list(std_df["sample"])
        if list(geno_df.index) != snps:
            print("Error, genotype does not match eQTL file.")
            exit()
        if list(alleles_df.index) != snps:
            print("Error, allele does not match eQTL file.")
            exit()
        # if list(expr_df.index) != probes:
        #     print("Error, expression does not match eQTL file.")
        #     exit()
        if list(geno_df.columns) != samples:
            print("Error, genotype does not match samples.")
            exit()
        if list(expr_df.columns) != samples:
            print("Error, expression does not match samples.")
            exit()
        if list(cova_df.columns) != samples:
            print("Error, covariates does not match samples.")
            exit()

        # Create the dataset map.
        std_dict = dict(zip(std_df["sample"], std_df["dataset"]))

        print("Iterating over eQTLs.")
        for row_index, (_, row) in enumerate(eqtl_df.iterrows()):
            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            # id = probe_name + "_" + snp_name
            id = snp_name

            interest_id = None
            for interest in self.interest:
                if id in interest:
                # if interest.startswith(id):
                    interest_id = interest

            if interest_id is None:
                continue

            probe_name = interest_id.split("_")[0]
            hgnc_name = probe_name

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
            # expression = expr_df.iloc[[row_index], :].copy().T
            expression = expr_df.loc[[probe_name], :].copy().T
            if expression.columns != [probe_name]:
                print("\t\tExpression file not in identical order as eQTL file.")
                exit()
            data = genotype.merge(expression, left_index=True, right_index=True)
            data.columns = ["genotype", "expression"]
            data.insert(0, "intercept", 1)
            data["group"] = data["genotype"].round(0)

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

            # Add the covariate of interest.
            covariate = interest_id.replace(probe_name + "_", "").replace(snp_name + "_", "").replace("_", " ")
            covariate_df = cova_df.loc[[covariate], :].T
            covariate_df.columns = ["covariate"]
            data = data.merge(covariate_df, left_index=True, right_index=True)
            data["interaction"] = data["genotype"] * data["covariate"]

            # Add the datasets.
            data["dataset"] = data.index.map(std_dict)

            # Remove missing values.
            data = data.loc[(data['genotype'] != -1) & (data['covariate'] != -1), :]

            # Determine annotation stats.
            eqtl_pvalue = OLS(data["expression"], data[["intercept", "genotype"]]).fit().pvalues[1]
            eqtl_pearsonr, _ = stats.pearsonr(data["expression"], data["genotype"])
            minor_allele_frequency = min(zero_geno_count, two_geno_count) / (zero_geno_count + two_geno_count)

            # Fill the eQTL plot annotation.
            annot1 = ["N: {:,}".format(data.shape[0]),
                      "r: {:.2f}".format(eqtl_pearsonr),
                      "p-value: {:.2e}".format(eqtl_pvalue),
                      "MAF: {:.2f}".format(minor_allele_frequency)]

            # Plot the main eQTL effect.
            self.eqtl_plot(df=data,
                           x="group",
                           y="expression",
                           palette=self.group_palette,
                           allele_map=allele_map,
                           xlabel=snp_name,
                           ylabel="{} expression".format(hgnc_name),
                           annot=annot1,
                           title="eQTL",
                           filename="{}_{}_{}_{}".format(row_index, probe_name, hgnc_name, snp_name)
                           )

            # Determine annotation stats.
            ieqtl_pvalue = OLS(data["expression"], data[["intercept", "genotype", "covariate", "interaction"]]).fit().pvalues[3]

            # Fill the interaction plot annotation.
            annot2 = ["{} - {}".format(snp_name.split(":")[2], minor_allele),
                      "N: {:,}".format(data.shape[0]),
                      "interaction p-value: {:.2e}".format(ieqtl_pvalue),
                      "eQTL p-value: {:.2e}".format(eqtl_pvalue),
                      "MAF: {:.2f}".format(minor_allele_frequency)]

            if len(set(data["covariate"].unique())) <= 2:
                covs = list(data["covariate"].unique())
                covs.sort()
                for cov in covs:
                    cov_data = data.loc[data["covariate"] == cov, :].copy()
                    coef, _ = stats.pearsonr(cov_data["expression"], cov_data["genotype"])
                    annot2.append("{} r: {:.2f}".format(cov, coef))
                    del cov_data

            # Plot the interaction eQTL.
            self.inter_plot(df=data,
                            x="covariate",
                            y="expression",
                            group="group",
                            group_palette=self.group_palette,
                            allele_map=allele_map,
                            xlabel=covariate,
                            ylabel="{} cell fraction".format(hgnc_name),
                            annot=annot2,
                            title="ieQTL",
                            filename="{}_{}_{}_{}_{}".format(row_index, probe_name, hgnc_name, snp_name, covariate.replace(" ", ""))
                            )

            for dataset in data["dataset"].unique().tolist() + ["AMPAD", "NoAMPAD"]:
                subset = data.loc[data["dataset"] == dataset, :].copy()
                if dataset == "AMPAD":
                    subset = data.loc[data["dataset"].isin(["AMPAD-MAYO-V2", "AMPAD-ROSMAP-V2", "AMPAD-MSBB-V2"]), :].copy()
                if dataset == "NoAMPAD":
                    subset = data.loc[~data["dataset"].isin(["AMPAD-MAYO-V2", "AMPAD-ROSMAP-V2", "AMPAD-MSBB-V2"]), :].copy()

                # Determine annotation stats.
                eqtl_pvalue = OLS(subset["expression"], subset[["intercept", "genotype"]]).fit().pvalues[1]
                ieqtl_pvalue = np.nan
                if subset["covariate"].std() > 0:
                    ieqtl_pvalue = OLS(subset["expression"], subset[["intercept", "genotype", "covariate", "interaction"]]).fit().pvalues[3]

                counts = subset["group"].value_counts()
                for x in [0.0, 1.0, 2.0]:
                    if x not in counts:
                        counts.loc[x] = 0
                zero_geno_count = (counts[0.0] * 2) + counts[1.0]
                two_geno_count = (counts[2.0] * 2) + counts[1.0]
                minor_allele_frequency = min(zero_geno_count, two_geno_count) / ( zero_geno_count + two_geno_count)

                # Fill the interaction plot annotation.
                annot3 = [
                    "{} - {}".format(snp_name.split(":")[2], minor_allele),
                    "N: {:,}".format(subset.shape[0]),
                    "interaction p-value: {:.2e}".format(ieqtl_pvalue),
                    "eQTL p-value: {:.2e}".format(eqtl_pvalue),
                    "MAF: {:.2f}".format(minor_allele_frequency)]

                if len(set(subset["covariate"].unique())) <= 2:
                    covs = list(subset["covariate"].unique())
                    covs.sort()
                    for cov in covs:
                        cov_data = subset.loc[subset["covariate"] == cov, :].copy()
                        if cov_data.shape[0] < 2:
                            continue
                        coef, _ = stats.pearsonr(cov_data["expression"], cov_data["genotype"])
                        annot3.append("{} r: {:.2f}".format(cov, coef))
                        del cov_data

                self.inter_plot(df=subset,
                                x="covariate",
                                y="expression",
                                group="group",
                                group_palette=self.group_palette,
                                allele_map=allele_map,
                                xlabel=covariate,
                                ylabel="{} expression".format(hgnc_name),
                                annot=annot3,
                                title="ieQTL - {}".format(dataset),
                                filename="{}_{}_{}_{}_{}_{}".format(row_index,
                                                                    probe_name,
                                                                    hgnc_name,
                                                                    snp_name,
                                                                    covariate.replace(" ", ""),
                                                                    dataset)
                                )

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def eqtl_plot(self, df, x="x", y="y", palette=None, allele_map=None,
                  annot=None, xlabel="", ylabel="", title="",
                  filename="eqtl_plot"):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
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

        for extension in self.extensions:
            outpath = os.path.join(self.outdir, "{}.{}".format(filename, extension))
            print("\t\tSaving plot: {}".format(os.path.basename(outpath)))
            fig.savefig(outpath)
        plt.close()

    def inter_plot(self, df, x="x", y="y", group="group",
                   group_palette=None, allele_map=None,
                   annot=None, xlabel="", ylabel="", title="",
                   filename="ieqtl_plot"):
        if len(set(df[group].unique()).symmetric_difference({0, 1, 2})) > 0:
            return

        if len(set(df[x].unique())) <= 2:
            self.binary_inter_plot(df=df,
                                   x=x,
                                   y=y,
                                   group=group,
                                   group_palette=group_palette,
                                   allele_map=allele_map,
                                   annot=annot,
                                   xlabel=xlabel,
                                   ylabel=ylabel,
                                   title=title,
                                   filename=filename)
        else:
            self.numerical_inter_plot(df=df,
                                      x=x,
                                      y=y,
                                      group=group,
                                      group_palette=group_palette,
                                      allele_map=allele_map,
                                      annot=annot,
                                      xlabel=xlabel,
                                      ylabel=ylabel,
                                      title=title,
                                      filename=filename)

    def binary_inter_plot(self, df, x="x", y="y", group="group",
                   group_palette=None, allele_map=None,
                   annot=None, xlabel="", ylabel="", title="",
                   filename="ieqtl_plot"):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        sns.boxplot(x=x,
                    y=y,
                    hue=group,
                    data=df,
                    palette={key: "#FFFFFF" for key in self.group_palette.keys()},
                    ax=ax)
        sns.swarmplot(x=x,
                      y=y,
                      hue=group,
                      data=df,
                      palette=self.group_palette,
                      dodge=True,
                      ax=ax)
        ax.get_legend().remove()

        if group_palette is not None:
            handles = []
            for group, color in group_palette.items():
                label = group
                if allele_map is not None:
                    label = allele_map[group]
                handles.append(mpatches.Patch(color=color, label=label))
            ax.legend(handles=handles, loc=4)

        if annot is not None:
            for i, annot_label in enumerate(annot):
                ax.annotate(annot_label,
                            xy=(0.03, 0.94 - (i * 0.04)),
                            xycoords=ax.transAxes,
                            color="#000000",
                            alpha=0.75,
                            fontsize=12,
                            fontweight='bold')

        ax.set_title(title,
                     fontsize=16,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        for extension in self.extensions:
            outpath = os.path.join(self.outdir, "{}.{}".format(filename, extension))
            print("\t\tSaving plot: {}".format(os.path.basename(outpath)))
            fig.savefig(outpath)
        plt.close()

    def numerical_inter_plot(self, df, x="x", y="y", group="group",
                   group_palette=None, allele_map=None,
                   annot=None, xlabel="", ylabel="", title="",
                   filename="ieqtl_plot"):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        for i, group_id in enumerate([0, 1, 2]):
            subset = df.loc[df[group] == group_id, :].copy()
            allele = group_id
            if allele_map is not None:
                allele = allele_map[group_id]

            coef_str = "NA"
            if len(subset.index) > 1:
                coef, p = stats.pearsonr(subset[y], subset[x])
                coef_str = "{:.2f}".format(coef)

                sns.regplot(x=x, y=y, data=subset, ci=None,
                            scatter_kws={'facecolors': group_palette[group_id],
                                         'linewidth': 0,
                                         'alpha': 0.3},
                            line_kws={"color": group_palette[group_id],
                                      "alpha": 0.75},
                            ax=ax
                            )

            ax.annotate(
                '{}: r = {}'.format(allele, coef_str),
                xy=(0.03, 0.94 - (i * 0.04)),
                xycoords=ax.transAxes,
                color=group_palette[group_id],
                alpha=0.75,
                fontsize=12,
                fontweight='bold')

        if annot is not None:
            for i, annot_label in enumerate(annot):
                ax.annotate(annot_label,
                            xy=(0.03, 0.82 - (i * 0.04)),
                            xycoords=ax.transAxes,
                            color="#000000",
                            alpha=0.75,
                            fontsize=12,
                            fontweight='bold')

        ax.set_title(title,
                     fontsize=16,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        for extension in self.extensions:
            outpath = os.path.join(self.outdir, "{}.{}".format(filename, extension))
            print("\t\tSaving plot: {}".format(os.path.basename(outpath)))
            fig.savefig(outpath)
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Alleles path: {}".format(self.alleles_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Covariate path: {}".format(self.cova_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Interest: {}".format(self.interest))
        print("  > Nrows: {}".format(self.nrows))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
