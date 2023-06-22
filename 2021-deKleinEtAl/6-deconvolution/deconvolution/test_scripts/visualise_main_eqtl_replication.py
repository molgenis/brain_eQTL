#!/usr/bin/env python3

"""
File:         visualise_main_eqtl_replication.py
Created:      2021/08/03
Last Changed: 2021/09/02
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
from colour import Color
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
__program__ = "Visualise Main eQTL Replication"
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
./visualise_main_eqtl_replication.py -dq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/Iteration1/eQTLProbesFDR0.05-ProbeLevel.txt.gz -dgte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/GTE_combined.txt.gz -dg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/genotypedump/GenotypeData.txt.gz -de /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/output-cortex/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz -dcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -dn EUR -rq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/replicateCortex-EUR/eQTLProbesFDR0.05-ProbeLevel.txt.gz -rgte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/combine_gte_files/GTE_combined.txt.gz -rg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/genotypedump_EUR_SNPs/GenotypeData.txt.gz -re /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/output-cortex/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz -rcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/perform_deconvolution/deconvolution_table.txt.gz -rn AFR
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.d_eqtl_path = getattr(arguments, 'discovery_eqtl')
        self.d_gte_path = getattr(arguments, 'discovery_gte')
        self.d_geno_path = getattr(arguments, 'discovery_genotype')
        self.d_allele_path = getattr(arguments, 'discovery_allele')
        self.d_expr_path = getattr(arguments, 'discovery_expression')
        self.d_cc_path = getattr(arguments, 'discovery_cell_counts')
        self.d_name = getattr(arguments, 'discovery_name')
        self.r_geno_path = getattr(arguments, 'replication_genotype')
        self.r_eqtl_path = getattr(arguments, 'replication_eqtl')
        self.r_gte_path = getattr(arguments, 'replication_gte')
        self.r_allele_path = getattr(arguments, 'replication_allele')
        self.r_expr_path = getattr(arguments, 'replication_expression')
        self.r_cc_path = getattr(arguments, 'replication_cell_counts')
        self.r_name = getattr(arguments, 'replication_name')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'visualise_main_eqtl_replication')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.colormap = {
            "minor": "#E69F00",
            "center": "#0072B2",
            "major": "#D55E00"
        }

        # Create color map.
        self.palette = {0.0: "#D55E00",
                        1.0: "#0072B2",
                        2.0: "#E69F00"}

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
        parser.add_argument("-dq",
                            "--discovery_eqtl",
                            type=str,
                            required=True,
                            help="The path to the discovery eqtl matrix.")
        parser.add_argument("-dgte",
                            "--discovery_gte",
                            type=str,
                            required=True,
                            help="The path to the discovery gte matrix.")
        parser.add_argument("-dg",
                            "--discovery_genotype",
                            type=str,
                            required=True,
                            help="The path to the discovery genotype matrix.")
        parser.add_argument("-da",
                            "--discovery_allele",
                            type=str,
                            required=False,
                            help="The path to the discovery allele matrix.")
        parser.add_argument("-de",
                            "--discovery_expression",
                            type=str,
                            required=True,
                            help="The path to the discovery expression matrix.")
        parser.add_argument("-dcc",
                            "--discovery_cell_counts",
                            type=str,
                            required=True,
                            help="The path to the discovery cell counts matrix.")
        parser.add_argument("-dn",
                            "--discovery_name",
                            type=str,
                            required=True,
                            help="The name for the discovery analysis.")
        parser.add_argument("-rq",
                            "--replication_eqtl",
                            type=str,
                            required=True,
                            help="The path to the replication eqtl matrix.")
        parser.add_argument("-rgte",
                            "--replication_gte",
                            type=str,
                            required=True,
                            help="The path to the replication gte matrix.")
        parser.add_argument("-rg",
                            "--replication_genotype",
                            type=str,
                            required=True,
                            help="The path to the replication genotype matrix.")
        parser.add_argument("-ra",
                            "--replication_allele",
                            type=str,
                            required=False,
                            help="The path to the replication allele matrix.")
        parser.add_argument("-re",
                            "--replication_expression",
                            type=str,
                            required=True,
                            help="The path to the replication expression "
                                 "matrix.")
        parser.add_argument("-rcc",
                            "--replication_cell_counts",
                            type=str,
                            required=True,
                            help="The path to the replication cell counts "
                                 "matrix.")
        parser.add_argument("-rn",
                            "--replication_name",
                            type=str,
                            required=True,
                            help="The name for the replication analysis.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()
        #
        # print("Loading discovery data")
        d_eqtl_df = self.load_file(self.d_eqtl_path, header=0, index_col=None)
        # d_gte_df = self.load_file(path=self.d_gte_path, header=None, index_col=None)
        # d_geno_df = self.load_file(self.d_geno_path, header=0, index_col=0)
        # d_expr_df = self.load_file(self.d_expr_path, header=0, index_col=0)
        #
        # # Pre-process.
        d_eqtl_df.index = d_eqtl_df["ProbeName"] + "_" + d_eqtl_df["SNPName"]
        # d_geno_df = d_geno_df.groupby(d_geno_df.index).first()
        # d_expr_df = d_expr_df.groupby(d_expr_df.index).first()
        #
        # print("Loading replication data")
        r_eqtl_df = self.load_file(self.r_eqtl_path, header=0, index_col=None)
        # r_gte_df = self.load_file(path=self.r_gte_path, header=None, index_col=None)
        # r_geno_df = self.load_file(self.r_geno_path, header=0, index_col=0)
        # r_expr_df = self.load_file(self.r_expr_path, header=0, index_col=0)
        #
        # # Pre-process.
        r_eqtl_df.index = r_eqtl_df["ProbeName"] + "_" + r_eqtl_df["SNPName"]
        # r_geno_df = r_geno_df.groupby(r_geno_df.index).first()
        # r_expr_df = r_expr_df.groupby(r_expr_df.index).first()
        #
        # print("Reorder matrices")
        # d_alle_df = d_geno_df.loc[d_eqtl_df["SNPName"], ["Alleles", "MinorAllele"]]
        # d_geno_df = d_geno_df.loc[d_eqtl_df["SNPName"], d_gte_df.iloc[:, 0]]
        # d_expr_df = d_expr_df.loc[d_eqtl_df["ProbeName"], d_gte_df.iloc[:, 1]]
        #
        # r_alle_df = r_geno_df.loc[r_eqtl_df["SNPName"], ["Alleles", "MinorAllele"]]
        # r_geno_df = r_geno_df.loc[r_eqtl_df["SNPName"], r_gte_df.iloc[:, 0]]
        # r_expr_df = r_expr_df.loc[r_eqtl_df["ProbeName"], r_gte_df.iloc[:, 1]]
        #
        # print("Checking matrices")
        # if list(d_alle_df.index) != list(d_eqtl_df["SNPName"].values):
        #     print("Unequal input matrix.")
        #     exit()
        # if list(d_geno_df.index) != list(d_eqtl_df["SNPName"].values):
        #     print("Unequal input matrix.")
        #     exit()
        # if list(d_expr_df.index) != list(d_eqtl_df["ProbeName"].values):
        #     print("Unequal input matrix.")
        #     exit()
        # if list(r_alle_df.index) != list(r_eqtl_df["SNPName"].values):
        #     print("Unequal input matrix.")
        #     exit()
        # if list(r_geno_df.index) != list(r_eqtl_df["SNPName"].values):
        #     print("Unequal input matrix.")
        #     exit()
        # if list(r_expr_df.index) != list(r_eqtl_df["ProbeName"].values):
        #     print("Unequal input matrix.")
        #     exit()
        #
        # print("Saving matrices")
        # self.save_file(df=d_alle_df, outpath=os.path.join("eQTLReplicationData", "Original_EMP_EUR_alleles.txt.gz"))
        # self.save_file(df=d_geno_df, outpath=os.path.join("eQTLReplicationData", "Original_EMP_EUR_expression.txt.gz"))
        # self.save_file(df=d_expr_df, outpath=os.path.join("eQTLReplicationData", "Original_EMP_EUR_genotype.txt.gz"))
        # self.save_file(df=r_alle_df, outpath=os.path.join("eQTLReplicationData", "Original_EMP_AFR_alleles.txt.gz"))
        # self.save_file(df=r_geno_df, outpath=os.path.join("eQTLReplicationData", "Original_EMP_AFR_expression.txt.gz"))
        # self.save_file(df=r_expr_df, outpath=os.path.join("eQTLReplicationData", "Original_EMP_AFR_genotype.txt.gz"))

        print("Loading matrices")
        d_alle_df = self.load_file(os.path.join("eQTLReplicationData", "Original_EMP_EUR_alleles.txt.gz"), header=0, index_col=0)
        d_geno_df = self.load_file(os.path.join("eQTLReplicationData", "Original_EMP_EUR_expression.txt.gz"), header=0, index_col=0)
        d_expr_df = self.load_file(os.path.join("eQTLReplicationData", "Original_EMP_EUR_genotype.txt.gz"), header=0, index_col=0)
        r_alle_df = self.load_file(os.path.join("eQTLReplicationData", "Original_EMP_AFR_alleles.txt.gz"), header=0, index_col=0)
        r_geno_df = self.load_file(os.path.join("eQTLReplicationData", "Original_EMP_AFR_expression.txt.gz"), header=0, index_col=0)
        r_expr_df = self.load_file(os.path.join("eQTLReplicationData", "Original_EMP_AFR_genotype.txt.gz"), header=0, index_col=0)

        print(d_alle_df)
        print(d_geno_df)
        print(d_expr_df)

        print(r_alle_df)
        print(r_geno_df)
        print(r_expr_df)

        ########################################################################

        print("Modelling discovery expression ~ genotype")
        d_coefs_a = []
        d_indices_a = []
        d_ma_a = []
        for i, index in enumerate(d_eqtl_df.index):
            genotype = d_geno_df.iloc[i, :].to_numpy()
            expression = d_expr_df.iloc[i, :].to_numpy()
            mask = genotype != -1
            # mask = np.where((genotype == 0) | (genotype == 1) | (genotype == 2))
            coef = np.nan
            if np.std(genotype[mask]) != 0 and np.std(expression[mask]) != 0:
                coef, _ = stats.pearsonr(genotype[mask], expression[mask])
            d_coefs_a.append(coef)
            d_indices_a.append(index)
            d_ma_a.append(d_alle_df.iloc[i, :]["MinorAllele"])
        d_coefs_df = pd.DataFrame({"x": d_coefs_a, "discovery MA": d_ma_a}, index=d_indices_a)
        del d_coefs_a, d_indices_a
        print(d_coefs_df)

        print("Getting Meta-Beta values")
        d_meta_beta = []
        for index, row in d_eqtl_df.iterrows():
            meta_beta_se_str = row['Meta-Beta (SE)']
            d_meta_beta.append([float(str(meta_beta_se_str).split(" (")[0]), row['AlleleAssessed']])
        d_meta_beta_df = pd.DataFrame(d_meta_beta, index=d_eqtl_df.index, columns=["discovery Meta-Beta", "AlleleAssessed"])
        print(d_meta_beta_df)

        print("Comparing correlation to Meta-Beta")
        discovery_df = d_coefs_df.merge(d_meta_beta_df, left_index=True, right_index=True)
        discovery_df["match"] = (discovery_df["discovery MA"] == discovery_df["AlleleAssessed"]).map({True: 1, False: -1})
        discovery_df["x"] = discovery_df["x"] * discovery_df["match"]
        discovery_df["z-score"] = (discovery_df["x"] - discovery_df["x"].mean()) / discovery_df["x"].std()
        self.plot_replication(df=discovery_df,
                              x="discovery Meta-Beta",
                              y="z-score",
                              xlabel="Meta-Beta",
                              ylabel="CorrCoef. Z-score",
                              name="corrcoef_vs_metabeta_{}".format(self.d_name),
                              title=self.d_name)

        ########################################################################

        print("Modelling replication expression ~ genotype")
        r_coefs_a = []
        r_indices_a = []
        r_ma_a = []
        for i, index in enumerate(r_eqtl_df.index):
            genotype = r_geno_df.iloc[i, :].to_numpy()
            expression = r_expr_df.iloc[i, :].to_numpy()
            mask = genotype != -1
            # mask = np.where((genotype == 0) | (genotype == 1) | (genotype == 2))
            coef = np.nan
            if np.std(genotype[mask]) != 0 and np.std(expression[mask]) != 0:
                coef, _ = stats.pearsonr(genotype[mask], expression[mask])
            r_coefs_a.append(coef)
            r_indices_a.append(index)
            r_ma_a.append(r_alle_df.iloc[i, :]["MinorAllele"])
        r_coefs_df = pd.DataFrame({"y": r_coefs_a, "replication MA": r_ma_a}, index=r_indices_a)
        del r_coefs_a, r_indices_a
        print(r_coefs_df)

        print("Getting Meta-Beta values")
        r_meta_beta = []
        for index, row in r_eqtl_df.iterrows():
            meta_beta_se_str = row['Meta-Beta (SE)']
            r_meta_beta.append([float(str(meta_beta_se_str).split(" (")[0]), row['AlleleAssessed']])
        r_meta_beta_df = pd.DataFrame(r_meta_beta, index=r_eqtl_df.index, columns=["replication Meta-Beta", "AlleleAssessed"])
        print(r_meta_beta_df)

        print("Comparing correlation to Meta-Beta")
        replication_df = r_coefs_df.merge(r_meta_beta_df, left_index=True, right_index=True)
        replication_df["match"] = (replication_df["replication MA"] == replication_df["AlleleAssessed"]).map({True: 1, False: -1})
        replication_df["y"] = replication_df["y"] * replication_df["match"]
        replication_df["z-score"] = (replication_df["y"] - replication_df["y"].mean()) / replication_df["y"].std()
        self.plot_replication(df=replication_df,
                              x="replication Meta-Beta",
                              y="z-score",
                              xlabel="Meta-Beta",
                              ylabel="CorrCoef. Z-score",
                              name="corrcoef_vs_metabeta_{}".format(self.r_name),
                              title=self.r_name)

        ########################################################################

        print("Combining data")
        plot_df = d_coefs_df.merge(r_coefs_df, left_index=True, right_index=True)
        plot_df.dropna(inplace=True)

        print("\tFlipping effects")
        plot_df["match"] = (plot_df["discovery MA"] == plot_df["replication MA"]).map({True: 1, False: -1})
        plot_df["y"] = plot_df["y"] * plot_df["match"]

        print("Plotting comparison")
        self.plot_replication(df=plot_df,
                              xlabel=self.d_name,
                              ylabel=self.r_name,
                              name="main_eqtl_replication_plot_{}_vs_{}".format(self.d_name, self.r_name),
                              title="expression ~ genotype")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=None, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
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
    def plot_eqtl(df, palette, ax, title="", xlabel="", ylabel="", annotate=None):
        # Calculate the correlation.
        coef, _ = stats.pearsonr(df["genotype"], df["expression"])

        # Plot the scatter / box plot.
        sns.regplot(x="genotype", y="expression", data=df,
                    scatter=False,
                    line_kws={"color": "#000000"},
                    ax=ax
                    )
        sns.boxplot(x="round_geno", y="expression", data=df,
                    palette=palette,
                    showfliers=False,
                    zorder=1,
                    ax=ax)

        ax.annotate(
            'N = {:,}'.format(df.shape[0]),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=12,
            fontweight='bold')
        ax.annotate(
            'r = {:.2f}'.format(coef),
            xy=(0.03, 0.90),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=12,
            fontweight='bold')

        if annotate is not None:
            for i, (label, value, rounding) in enumerate(annotate):
                ax.annotate(
                    '{} = {:{}}'.format(label, value, rounding),
                    xy=(0.03, 0.86 - (0.04 * i)),
                    xycoords=ax.transAxes,
                    color="#000000",
                    alpha=0.75,
                    fontsize=12,
                    fontweight='bold')

        ax.set_title(title,
                     fontsize=22,
                     weight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

    @staticmethod
    def plot_inter_eqtl(df, palette, ax, title="", xlabel="", ylabel="",
                        annotate=None):

        for i, genotype in enumerate([0.0, 1.0, 2.0]):
            subset = df.loc[df["round_geno"] == genotype, :].copy()
            color = palette[genotype]
            coef = np.nan
            if len(subset.index) > 1:
                # Calculate the correlation.
                coef, _ = stats.pearsonr(subset["cell count"], subset["expression"])

                # Plot the scatter / box plot.
                sns.regplot(x="cell count", y="expression", data=subset,
                            scatter_kws={'facecolors': color,
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": color, "alpha": 0.75},
                            ax=ax
                            )

            ax.annotate(
                '{}: r = {:.2f} [N = {:,}]'.format(genotype, coef, subset.shape[0]),
                xy=(0.03, 0.94 - (0.04 * i)),
                xycoords=ax.transAxes,
                color=color,
                alpha=0.75,
                fontsize=12,
                fontweight='bold')

        if annotate is not None:
            for i, (label, value, rounding) in enumerate(annotate):
                ax.annotate(
                    '{} = {:{}}'.format(label, value, rounding),
                    xy=(0.03, 0.82 - (0.04 * i)),
                    xycoords=ax.transAxes,
                    color="#000000",
                    alpha=0.75,
                    fontsize=12,
                    fontweight='bold')

        ax.set_title(title,
                     fontsize=22,
                     weight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

    def plot_replication(self, df, x="x", y="y", hue=None, xlabel="",
                         ylabel="", name="", title=""):
        if df.shape[0] <= 2:
            return

        facecolors = "#000000"
        if hue is not None:
            facecolors = df[hue]

        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(12, 12))
        sns.set(color_codes=True)

        sns.despine(fig=fig, ax=ax)

        lower_quadrant = df.loc[(df[x] < 0) & (df[y] < 0), :]
        upper_quadrant = df.loc[(df[x] > 0) & (df[y] > 0), :]
        concordance = (100 / df.shape[0]) * (lower_quadrant.shape[0] + upper_quadrant.shape[0])

        coef, _ = stats.pearsonr(df[y], df[x])

        sns.regplot(x=x, y=y, data=df, ci=None,
                    scatter_kws={'facecolors': facecolors,
                                 'linewidth': 0,
                                 'alpha': 0.75},
                    line_kws={"color": "#0072B2",
                              'linewidth': 5},
                    ax=ax)

        ax.annotate(
            'N = {}'.format(df.shape[0]),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=18,
            fontweight='bold')
        ax.annotate(
            'r = {:.2f}'.format(coef),
            xy=(0.03, 0.90),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=18,
            fontweight='bold')
        ax.annotate(
            'concordance = {:.0f}%'.format(concordance),
            xy=(0.03, 0.86),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=18,
            fontweight='bold')

        ax.axhline(0, ls='--', color="#000000", zorder=-1)
        ax.axvline(0, ls='--', color="#000000", zorder=-1)

        ax.set_xlabel(xlabel,
                      fontsize=20,
                      fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=20,
                      fontweight='bold')
        ax.set_title(title,
                     fontsize=25,
                     fontweight='bold')

        outpath = os.path.join(self.outdir, "{}.png".format(name))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved: {}".format(outpath))

    def print_arguments(self):
        print("Arguments:")
        print("  > {} [discovery]:".format(self.d_name))
        print("    > eQTL: {}".format(self.d_eqtl_path))
        print("    > GtE: {}".format(self.d_gte_path))
        print("    > Genotype: {}".format(self.d_geno_path))
        print("    > Allele: {}".format(self.d_allele_path))
        print("    > Expression: {}".format(self.d_expr_path))
        print("    > Cell counts: {}".format(self.d_cc_path))
        print("  > {} [replication]:".format(self.r_name))
        print("    > eQTL: {}".format(self.r_eqtl_path))
        print("    > GtE: {}".format(self.r_gte_path))
        print("    > Genotype: {}".format(self.r_geno_path))
        print("    > Allele: {}".format(self.r_allele_path))
        print("    > Expression: {}".format(self.r_expr_path))
        print("    > Cell counts: {}".format(self.r_cc_path))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
