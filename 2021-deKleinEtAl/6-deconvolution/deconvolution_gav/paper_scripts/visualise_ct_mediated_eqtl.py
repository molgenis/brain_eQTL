#!/usr/bin/env python3

"""
File:         visualise_ct_mediated_eqtl.py
Created:      2020/09/23
Last Changed: 2021/02/22
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
import pandas as pd
import seaborn as sns
import matplotlib
from statsmodels.stats import multitest
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Visualise CellType mediated eQTL"
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
./visualise_ct_mediated_eqtl.py -eq ../matrix_preparation/cortex_eur_cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ge ../matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -al ../matrix_preparation/cortex_eur_cis/create_matrices/genotype_alleles.txt.gz -ex ../matrix_preparation/cortex_eur_cis/create_matrices/expression_table.txt.gz -cc ../matrix_preparation/cortex_eur_cis/perform_deconvolution/deconvolution_table.txt.gz -d ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -i CLECL1 CYP24A1 -n 600 -e pdf
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
        self.decon_path = getattr(arguments, 'decon')
        self.alpha = getattr(arguments, 'alpha')
        self.interest = getattr(arguments, 'interest')
        self.nrows = getattr(arguments, 'nrows')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'visualise_ct_mediated_eqtl')

        self.colormap = {
            "minor": "#E69F00",
            "center": "#0072B2",
            "major": "#D55E00",
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00"
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

        # Create color map.
        self.group_color_map, self.value_color_map = self.create_color_map(self.colormap)

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
        parser.add_argument("-d",
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix")
        parser.add_argument("-a",
                            "--alpha",
                            type=float,
                            required=False,
                            default=0.05,
                            help="The significance cut-off. Default: 0.05.")
        parser.add_argument("-i",
                            "--interest",
                            nargs="+",
                            type=str,
                            required=False,
                            default=None,
                            help="The HGNCSymbols to plot. Default: none.")
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

    @staticmethod
    def create_color_map(colormap):
        major_small = list(Color(colormap["major"]).range_to(Color("#FFFFFF"), 12))[2]
        center_small = list(Color(colormap["center"]).range_to(Color("#FFFFFF"), 12))[2]

        palette = list(Color(colormap["major"]).range_to(Color(major_small), 50)) + \
                  list(Color(major_small).range_to(Color(colormap["center"]), 50)) + \
                  list(Color(colormap["center"]).range_to(Color(center_small), 50)) + \
                  list(Color(center_small).range_to(Color(colormap["minor"]), 51))
        colors = [str(x).upper() for x in palette]
        values = [x / 100 for x in list(range(201))]
        group_color_map = {0.0: colormap["major"], 1.0: colormap["center"], 2.0: colormap["minor"]}
        value_color_map = {}
        for val, col in zip(values, colors):
            value_color_map[val] = col
        return group_color_map, value_color_map

    def start(self):
        self.print_arguments()
        eqtl_df, geno_df, alleles_df, expr_df, cc_df, decon_df = self.load()
        _, decon_fdr_df = self.bh_correct(decon_df)
        self.plot_eqtls(eqtl_df, geno_df, alleles_df, expr_df, cc_df, decon_fdr_df)

    def load(self):
        print("Loading input files.")
        eqtl_df = self.load_file(self.eqtl_path, index_col=None, nrows=self.nrows)
        geno_df = self.load_file(self.geno_path, nrows=self.nrows)
        alleles_df = self.load_file(self.alleles_path, nrows=self.nrows)
        expr_df = self.load_file(self.expr_path, nrows=self.nrows)
        cc_df = self.load_file(self.cc_path)
        decon_df = self.load_file(self.decon_path)

        return eqtl_df, geno_df, alleles_df, expr_df, cc_df, decon_df

    @staticmethod
    def bh_correct(pvalue_df):
        df = pvalue_df.copy()
        pval_data = []
        fdr_data = []
        indices = []
        for col in df.columns:
            if col.endswith("_pvalue"):
                pval_data.append(df.loc[:, col])
                fdr_data.append(multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
                indices.append(col.replace("_pvalue", ""))
        pval_df = pd.DataFrame(pval_data, index=indices, columns=df.index)
        fdr_df = pd.DataFrame(fdr_data, index=indices, columns=df.index)

        return pval_df.T, fdr_df.T

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def plot_eqtls(self, eqtl_df, geno_df, alleles_df, expr_df, cc_df, decon_df):
        print("Plotting interaction eQTL plots.")

        print("Iterating over eQTLs.")
        for i, (index, row) in enumerate(eqtl_df.iterrows()):
            # Extract the usefull information from the row.
            fdr_value = row["FDR"]
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            eqtl_type = row["CisTrans"]

            if self.interest is not None and hgnc_name not in self.interest:
                continue

            print("\tWorking on: {}\t{}\t{} [{}/{} "
                  "{:.2f}%]".format(snp_name, probe_name, hgnc_name,
                                    i + 1,
                                    eqtl_df.shape[0],
                                    (100 / eqtl_df.shape[0]) * (i + 1)))

            # Get the genotype / expression data.
            genotype = geno_df.iloc[i, :].T.to_frame()
            if genotype.columns != [snp_name]:
                print("\t\tGenotype file not in identical order as eQTL file.")
                exit()
            expression = expr_df.iloc[i, :].T.to_frame()
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
            alleles = alleles_df.iloc[i, :]
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
            data["alleles"] = data["group"].map(allele_map)

            # Determine the minor allele frequency.
            minor_allele_frequency = min(zero_geno_count, two_geno_count) / (
                        zero_geno_count + two_geno_count)

            # Add the color.
            data["round_geno"] = data["genotype"].round(2)
            data["value_hue"] = data["round_geno"].map(self.value_color_map)
            data["group_hue"] = data["group"].map(self.group_color_map)

            # Check if the SNP has an interaction effect.
            interaction_effect = decon_df.loc["{}_{}".format(probe_name, snp_name), :]
            if not interaction_effect.name.startswith(probe_name) or not interaction_effect.name.endswith(snp_name):
                print("\t\tCannot find probe-snp combination in decon results.")
                exit()
            interaction_effect = interaction_effect.to_frame()
            interaction_effect.columns = ["FDR"]
            interaction_effect = interaction_effect.loc[
                                 interaction_effect["FDR"] < self.alpha, :]
            interaction_effect = interaction_effect.reindex(
                interaction_effect["FDR"].abs().sort_values(
                    ascending=True).index)

            # Prepare output directory.
            if len(interaction_effect.index) > 0:
                eqtl_interaction_outdir = os.path.join(self.outdir,
                                                       "{}_{}_{}_{}".format(
                                                           index, snp_name,
                                                           probe_name,
                                                           hgnc_name))
                if not os.path.exists(eqtl_interaction_outdir):
                    os.makedirs(eqtl_interaction_outdir)

                self.plot_simple_eqtl(i, fdr_value, snp_name, probe_name,
                                      hgnc_name, eqtl_type, data,
                                      minor_allele, minor_allele_frequency,
                                      allele_map, self.group_color_map,
                                      eqtl_interaction_outdir,
                                      self.extensions)

                for index2, (fdr,) in interaction_effect.iterrows():
                    print("\t\t{}: FDR = {}".format(index2, fdr))
                    eqtl_data = data.copy()
                    cov_data = cc_df.loc[:, index2].to_frame()
                    eqtl_data = eqtl_data.merge(cov_data, left_index=True,
                                                right_index=True)

                    self.plot_inter_eqtl(i, snp_name, probe_name, hgnc_name,
                                         eqtl_type, eqtl_data, index2, fdr,
                                         allele_map, self.group_color_map,
                                         eqtl_interaction_outdir,
                                         self.extensions)
            else:
                print("\t\tNo significant interactions.")

    @staticmethod
    def plot_simple_eqtl(count, fdr_value, snp_name, probe_name, hgnc_name,
                         eqtl_type, df, minor_allele, minor_allele_frequency,
                         allele_map, group_color_map, outdir, extensions):
        # Calculate the correlation.
        coef, p = stats.spearmanr(df["genotype"],
                                  df["expression"])

        # Prepare the figure.
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        # Plot the scatter / box plot.
        sns.regplot(x="genotype", y="expression", data=df,
                    scatter=False,
                    # scatter_kws={'facecolors': df['value_hue'],
                    #              'edgecolors': "#808080"},
                    line_kws={"color": "#000000"},
                    ax=ax
                    )
        sns.boxplot(x="group", y="expression", data=df,
                    palette=group_color_map,
                    showfliers=False,
                    zorder=1,
                    # boxprops=dict(alpha=.3),
                    ax=ax)
        # plt.setp(ax.artists, edgecolor='k', facecolor='w')
        # plt.setp(ax.lines, color='k')

        # Set the other aesthetics.
        ax.set_xticks(range(3))
        ax.set_xticklabels([allele_map[0.0], allele_map[1.0], allele_map[2.0]])
        ax.text(0.5, 1.06,
                # '{} {}-eQTL [{}]'.format(hgnc_name, eqtl_type,
                #                          p_value_to_symbol(p_value)),
                '{} {}-eQTL [FDR = {:.2e}]'.format(hgnc_name, eqtl_type, fdr_value),
                fontsize=22, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                # 'r = {:.2f} [{}]    minor allele frequency '
                # '{} = {:.2f}'.format(coef,
                #                      p_value_to_symbol(p),
                #                      minor_allele,
                #                      minor_allele_frequency),
                'r = {:.2f} [p = {:.2e}]    minor allele frequency '
                '{} = {:.2f}'.format(coef,
                                     p,
                                     minor_allele,
                                     minor_allele_frequency),
                fontsize=14, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)
        ax.set_ylabel('{} ({}) expression'.format(probe_name, hgnc_name),
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel('SNP {}'.format(snp_name),
                      fontsize=14,
                      fontweight='bold')

        # Safe the plot.
        for extension in extensions:
            filename = "{}_simple_eqtl_{}_{}_{}.{}".format(count,
                                                           snp_name,
                                                           probe_name,
                                                           hgnc_name,
                                                           extension)
            print("\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(outdir, filename), dpi=300)
        plt.close()

    @staticmethod
    def plot_inter_eqtl(count, snp_name, probe_name, hgnc_name, eqtl_type, df,
                        cov_name, fdr, allele_map, group_color_map, outdir,
                        extensions):
        # calculate axis limits.
        ymin_value = df["expression"].min()
        ymin = ymin_value - abs(ymin_value * 0.2)
        ymax_value = df["expression"].max()
        ymax = ymax_value + abs(ymax_value * 0.6)

        xmin_value = df[cov_name].min()
        xmin = xmin_value - max(abs(xmin_value * 0.1), 0.05)
        xmax_value = df[cov_name].max()
        xmax = xmax_value + max(abs(xmax_value * 0.1), 0.05)

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        label_pos = {0.0: 0.94, 1.0: 0.90, 2.0: 0.86}
        for i, genotype in enumerate([1.0, 0.0, 2.0]):
            subset = df.loc[df["round_geno"] == genotype, :].copy()
            color = group_color_map[genotype]
            allele = allele_map[genotype]

            coef_str = "NA"
            p_str = "NA"
            if len(subset.index) > 1:
                # Regression.
                coef, p = stats.spearmanr(subset["expression"],
                                          subset[cov_name])
                coef_str = "{:.2f}".format(coef)
                # p_str = p_value_to_symbol(p)
                p_str = "p = {:.2e}".format(p)

                # Plot.
                sns.regplot(x=cov_name, y="expression", data=subset,
                            scatter_kws={'facecolors': subset['value_hue'],
#                                         'edgecolors': subset['group_hue'],
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": color, "alpha": 0.75},
                            ax=ax
                            )

            # Add the text.
            ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
            ax.annotate(
                '{}: r = {} [{}]'.format(allele, coef_str, p_str),
                # xy=(0.03, 0.94 - ((i / 100) * 4)),
                xy=(0.03, label_pos[genotype]),
                xycoords=ax.transAxes,
                color=color,
                alpha=0.75,
                fontsize=12,
                fontweight='bold')

        ax.text(0.5, 1.06,
                '{} {}-eQTL Interaction with {} '.format(hgnc_name,
                                                         eqtl_type,
                                                         cov_name),
                fontsize=18, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                'SNPName: {}  ProbeName: {}  '
                'FDR: {:.2e}'.format(snp_name, probe_name, fdr),
                fontsize=12, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel('{} ({}) expression'.format(probe_name, hgnc_name),
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(cov_name,
                      fontsize=14,
                      fontweight='bold')

        # ax.axvline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)
        # ax.axhline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)

        # Safe the plot.
        plt.tight_layout()
        for extension in extensions:
            filename = "{}_inter_eqtl_{}_{}_{}_{}.{}".format(
                count,
                snp_name,
                probe_name,
                hgnc_name,
                cov_name,
                extension)
            print("\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(outdir, filename), dpi=300)
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Alleles path: {}".format(self.alleles_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell count % path: {}".format(self.cc_path))
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Alpha: {}".format(self.alpha))
        print("  > Interest: {}".format(self.interest))
        print("  > Nrows: {}".format(self.nrows))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
