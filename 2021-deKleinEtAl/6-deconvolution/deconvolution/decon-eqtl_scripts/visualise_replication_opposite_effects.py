#!/usr/bin/env python3

"""
File:         visualise_replication_opposite_effects.py
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
__program__ = "Visualise Replication Opposite Effects"
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
./visualise_replication_opposite_effects.py -dq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/Iteration1/eQTLProbesFDR0.05-ProbeLevel.txt.gz -dg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -de /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt -dcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -dn EUR -rq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/replicateCortex-EUR/eQTLProbesFDR0.05-ProbeLevel.txt.gz -rg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR/genotype_table.txt -re /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt -rcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/perform_deconvolution/deconvolution_table.txt.gz -rn AFR -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_replication_plot/data/replication_opposite_efects_table.txt.gz

./visualise_replication_opposite_effects.py -dq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/Iteration1/eQTLProbesFDR0.05-ProbeLevel.txt.gz -dg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -de /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -dcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -dn EUR -rq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/replicateCortex-EUR/eQTLProbesFDR0.05-ProbeLevel.txt.gz -rg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR-Normalised/genotype_table.txt -re /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -rcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/perform_deconvolution/deconvolution_table.txt.gz -rn AFR -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_replication_plot/data/replication_opposite_efects_table.txt.gz

./visualise_replication_opposite_effects.py -dq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/Iteration1/eQTLProbesFDR0.05-ProbeLevel.txt.gz -dg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -de /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/expression_table.txt.gz -dcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -dn EUR -rq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/replicateCortex-EUR/eQTLProbesFDR0.05-ProbeLevel.txt.gz -rg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR/genotype_table.txt -re /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/create_matrices/expression_table.txt.gz -rcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/perform_deconvolution/deconvolution_table.txt.gz -rn AFR -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_replication_plot/data/replication_opposite_efects_table.txt.gz

./visualise_replication_opposite_effects.py -dq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/Iteration1/eQTLProbesFDR0.05-ProbeLevel.txt.gz -dg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -de /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/expression_table.txt.gz -dcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -dn EUR -rq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/replicateCortex-EUR/eQTLProbesFDR0.05-ProbeLevel.txt.gz -rg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/create_matrices/genotype_table.txt.gz -re /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/create_matrices/expression_table.txt.gz -rcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/perform_deconvolution/deconvolution_table.txt.gz -rn AFR -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_replication_plot/data/replication_opposite_efects_table.txt.gz

./visualise_replication_opposite_effects.py -dq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/Iteration1/eQTLProbesFDR0.05-ProbeLevel.txt.gz -dg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/genotypedump/GenotypeData.txt.gz -de /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/output-cortex/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz -dcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -dn EUR -rq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/replicateCortex-EUR/eQTLProbesFDR0.05-ProbeLevel.txt.gz -rg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/genotypedump_EUR_SNPs/GenotypeData.txt.gz -re /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/output-cortex/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz -rcc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/perform_deconvolution/deconvolution_table.txt.gz -rn AFR -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_replication_plot/data/replication_opposite_efects_table.txt.gz
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.d_eqtl_path = getattr(arguments, 'discovery_eqtl')
        self.d_geno_path = getattr(arguments, 'discovery_genotype')
        self.d_allele_path = getattr(arguments, 'discovery_allele')
        self.d_expr_path = getattr(arguments, 'discovery_expression')
        self.d_cc_path = getattr(arguments, 'discovery_cell_counts')
        self.d_name = getattr(arguments, 'discovery_name')
        self.r_geno_path = getattr(arguments, 'replication_genotype')
        self.r_eqtl_path = getattr(arguments, 'replication_eqtl')
        self.r_allele_path = getattr(arguments, 'replication_allele')
        self.r_expr_path = getattr(arguments, 'replication_expression')
        self.r_cc_path = getattr(arguments, 'replication_cell_counts')
        self.r_name = getattr(arguments, 'replication_name')
        self.interest_path = getattr(arguments, 'interest')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'visualise_opposite_effects')
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
        parser.add_argument("-i",
                            "--interest",
                            type=str,
                            required=True,
                            help="The path to the interest matrix.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading interest")
        interest_df = self.load_file(self.interest_path, header=0, index_col=0)
        # nrows = max(max(interest_df["discovery index"]), max(interest_df["replication index"])) + 1
        nrows = None

        print("Loading discovery data")
        d_eqtl_df = self.load_file(self.d_eqtl_path, header=0, index_col=None, nrows=nrows)
        d_eqtl_df.index = d_eqtl_df["ProbeName"] + "_" + d_eqtl_df["SNPName"]
        d_geno_df = self.load_file(self.d_geno_path, header=0, index_col=0, nrows=nrows)
        d_geno_df = d_geno_df.groupby(d_geno_df.index).first()
        # d_allele_df = self.load_file(self.d_allele_path, header=0, index_col=0, nrows=nrows)
        d_expr_df = self.load_file(self.d_expr_path, header=0, index_col=0, nrows=nrows)
        d_expr_df = d_expr_df.groupby(d_expr_df.index).first()
        d_cc_df = self.load_file(self.d_cc_path, header=0, index_col=0)

        print("Loading replication data")
        r_eqtl_df = self.load_file(self.r_eqtl_path, header=0, index_col=None, nrows=nrows)
        r_eqtl_df.index = r_eqtl_df["ProbeName"] + "_" + r_eqtl_df["SNPName"]
        r_geno_df = self.load_file(self.r_geno_path, header=0, index_col=0, nrows=nrows)
        r_geno_df = r_geno_df.groupby(r_geno_df.index).first()
        # r_allele_df = self.load_file(self.r_allele_path, header=0, index_col=0, nrows=nrows)
        r_expr_df = self.load_file(self.r_expr_path, header=0, index_col=0, nrows=nrows)
        r_expr_df = r_expr_df.groupby(r_expr_df.index).first()
        r_cc_df = self.load_file(self.r_cc_path, header=0, index_col=0)

        print("Looping over interest")
        for _, (index, cell_type, d_pvalue, d_fdr, d_beta, r_pvalue, r_beta, r_fdr) in interest_df.iterrows():
            gene = index.split("_")[0]
            snp = "_".join(index.split("_")[1:])
            print("\tPlotting {} - {} - {}".format(gene, snp, cell_type))

            # Initialize plot.
            sns.set(rc={'figure.figsize': (24, 18)})
            sns.set_style("ticks")
            fig, axes = plt.subplots(ncols=2, nrows=2)

            for eqtl_ax, inter_ax, eqtl_df, geno_df, expr_df, cc_df, pvalue, fdr, \
                beta, name in ([*axes[:, 0], d_eqtl_df, d_geno_df, d_expr_df, d_cc_df, d_pvalue, d_fdr, d_beta, self.d_name],
                               [*axes[:, 1], r_eqtl_df, r_geno_df, r_expr_df, r_cc_df, r_pvalue, r_fdr, r_beta, self.r_name]):
                sns.despine(fig=fig, ax=eqtl_ax)
                sns.despine(fig=fig, ax=inter_ax)

                eqtl_pvalue = eqtl_df.loc[gene + "_" + snp, "PValue"]
                eqtl_fdr = eqtl_df.loc[gene + "_" + snp, "FDR"]

                genotype = geno_df.loc[snp, :]
                expression = expr_df.loc[gene, :]
                cell_count = cc_df.loc[:, cell_type]

                plot_df = pd.DataFrame({"genotype": genotype, "expression": expression, "cell count": cell_count})
                plot_df["round_geno"] = np.rint(genotype)
                plot_df = plot_df.loc[plot_df["genotype"] != -1, :]

                self.plot_eqtl(df=plot_df,
                               palette=self.palette,
                               ax=eqtl_ax,
                               title=name,
                               xlabel=snp,
                               ylabel=gene,
                               annotate=[("eQTL p-value", eqtl_pvalue, ".2e"), ("eQTL FDR", eqtl_fdr, ".2e")])
                self.plot_inter_eqtl(df=plot_df,
                                     palette=self.palette,
                                     ax=inter_ax,
                                     title="",
                                     xlabel=cell_type,
                                     ylabel=gene,
                                     annotate=[("Decon-eQTL p-value", pvalue, ".2e"),
                                               ("Decon-eQTL FDR", fdr, ".2e"),
                                               ("Decon-eQTL beta", beta, ".2f")])

            outpath = os.path.join(self.outdir, "opposite_effects_plot_{}_vs_{}_{}_{}_{}.png".format(self.d_name, self.r_name, gene, snp, cell_type))
            fig.savefig(outpath)
            plt.close()
            print("\tSaved: {}".format(outpath))

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=None, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

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

    def print_arguments(self):
        print("Arguments:")
        print("  > {} [discovery]:".format(self.d_name))
        print("    > eQTL: {}".format(self.d_eqtl_path))
        print("    > Genotype: {}".format(self.d_geno_path))
        print("    > Allele: {}".format(self.d_allele_path))
        print("    > Expression: {}".format(self.d_expr_path))
        print("    > Cell counts: {}".format(self.d_cc_path))
        print("  > {} [replication]:".format(self.r_name))
        print("    > eQTL: {}".format(self.r_eqtl_path))
        print("    > Genotype: {}".format(self.r_geno_path))
        print("    > Allele: {}".format(self.r_allele_path))
        print("    > Expression: {}".format(self.r_expr_path))
        print("    > Cell counts: {}".format(self.r_cc_path))
        print("  > Interest: {}".format(self.interest_path))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
