#!/usr/bin/env python3

"""
File:         visualise_sn_eqtl.py
Created:      2021/01/28
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
import math
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
__program__ = "Visualise Single-Nucleus eQTL"
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
./visualise_sn_eqtl.py -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/cis_100Perm/MASK/eQTLsFDR-ProbeLevel.txt.gz -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/genotypedump/GenotypeData.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/MASK_expression.txt -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/ROSMAP-scRNAseq-genometoexpressioncoupling.txt -i STMN4 NKAIN1 FAM221A SLC25A27 -e pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cell_types = getattr(arguments, 'cell_types')
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.gte_path = getattr(arguments, 'geno_to_expr_coupling')
        self.alpha = getattr(arguments, 'alpha')
        self.interest = getattr(arguments, 'interest')
        self.nrows = getattr(arguments, 'nrows')
        self.extensions = getattr(arguments, 'extension')

        # Replace default cell types.
        if self.cell_types == "all":
            self.cell_types = ["EX", "OLI", "IN", "AST", "OPC", "MIC", "PER", "END"]

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'visualise_sn_eqtl')

        self.colormap = {
            "minor": "#E69F00",
            "center": "#0072B2",
            "major": "#D55E00",
            "EX": "#0072B2",
            "IN": "#0072B2",
            "OLI": "#009E73",
            "OPC": "#009E73",
            "END": "#CC79A7",
            "MIC": "#E69F00",
            "AST": "#D55E00",
            "PER": "#F0E442"
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

        # Create color map.
        self.group_color_map, self.value_color_map = self.create_color_map(
            self.colormap)

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
        parser.add_argument("-ct",
                            "--cell_types",
                            nargs="*",
                            type=str,
                            required=False,
                            default="all",
                            choices=["EX", "IN", "OLI", "OPC", "END", "MIC",
                                     "AST"],
                            help="The cell type to visualise. Default: all.")
        parser.add_argument("-eq",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the eqtl matrix. Note: use "
                                 "'MASK' to declare an file specific to each"
                                 "cell type.")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix. Note: use "
                                 "'MASK' to declare an file specific to each"
                                 "cell type.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix. Note: use "
                                 "'MASK' to declare an file specific to each"
                                 "cell type.")
        parser.add_argument("-gte",
                            "--geno_to_expr_coupling",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to a genotype (left column) to "
                                 "expression (right column) coupling table. "
                                 "Note: use 'MASK' to declare an file specific "
                                 "to each cell type.")
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
                            required=True,
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
        major_small = \
        list(Color(colormap["major"]).range_to(Color("#FFFFFF"), 12))[2]
        center_small = \
        list(Color(colormap["center"]).range_to(Color("#FFFFFF"), 12))[2]

        palette = list(
            Color(colormap["major"]).range_to(Color(major_small), 50)) + \
                  list(Color(major_small).range_to(Color(colormap["center"]),
                                                   50)) + \
                  list(Color(colormap["center"]).range_to(Color(center_small),
                                                          50)) + \
                  list(Color(center_small).range_to(Color(colormap["minor"]),
                                                    51))
        colors = [str(x).upper() for x in palette]
        values = [x / 100 for x in list(range(201))]
        group_color_map = {0.0: colormap["major"], 1.0: colormap["center"],
                           2.0: colormap["minor"]}
        value_color_map = {}
        for val, col in zip(values, colors):
            value_color_map[val] = col
        return group_color_map, value_color_map

    def start(self):
        self.print_arguments()

        print("### Loading global data ###")
        geno_df = None
        alleles_df = None
        if "MASK" not in self.geno_path:
            print("### Loading genotype file ###")
            geno_df = self.load_file(self.geno_path, index_col=0,
                                     nrows=self.nrows)

            # Split the dump file.
            alleles_df = geno_df[["Alleles", "MinorAllele"]].copy()
            geno_df.drop(["Alleles", "MinorAllele"], axis=1, inplace=True)

        gte_dict = None
        if self.gte_path is not None and "MASK" not in self.gte_path:
            gte_df = self.load_file(self.gte_path)
            gte_dict = dict(zip(gte_df.iloc[:, 1], gte_df.iloc[:, 0]))

        print("### Iterating through cell types ###")
        interest_corr = {x: [] for x in self.interest}
        for ct in self.cell_types:
            ct_outdir = os.path.join(self.outdir, ct)

            print("  ### Loading {} input files ###".format(ct))
            if "MASK" in self.geno_path:
                geno_df = self.load_file(self.geno_path.replace("MASK", ct),
                                         index_col=0, nrows=self.nrows)

                # Split the dump file.
                alleles_df = geno_df[["Alleles", "MinorAllele"]].copy()
                geno_df.drop(["Alleles", "MinorAllele"], axis=1, inplace=True)

            if self.gte_path is not None and "MASK" in self.gte_path:
                gte_df = self.load_file(self.gte_path.replace("MASK", ct))
                gte_dict = dict(zip(gte_df.iloc[:, 1], gte_df.iloc[:, 0]))

            eqtl_df = self.load_file(self.eqtl_path.replace("MASK", ct),
                                     nrows=self.nrows)
            expr_df = self.load_file(self.expr_path.replace("MASK", ct),
                                     index_col=0, nrows=self.nrows)

            print("  ### PLotting {} eQTLs ###".format(ct))
            # Plot.
            for i, (index, row) in enumerate(eqtl_df.iterrows()):
                # Extract the usefull information from the row.
                fdr_value = row["FDR"]
                snp_name = row["SNPName"]
                probe_name = row["ProbeName"]
                hgnc_name = row["HGNCName"]
                eqtl_type = row["CisTrans"]

                if hgnc_name not in self.interest:
                    continue

                print("\tWorking on: {}\t{}\t{} [{}/{} "
                      "{:.2f}%]".format(snp_name, probe_name, hgnc_name,
                                        i + 1,
                                        eqtl_df.shape[0],
                                        (100 / eqtl_df.shape[0]) * (i + 1)))

                # Prep output directory.
                if not os.path.exists(ct_outdir):
                    os.makedirs(ct_outdir)

                # Get the genotype / expression data.
                genotype = geno_df.loc[[snp_name], :].copy().T
                if genotype.columns != [snp_name]:
                    print("\t\tCould not extract valid genotype data for '{}'".format(snp_name))
                    exit()
                expression = expr_df.loc[[probe_name], :].copy().T
                if expression.columns != [probe_name]:
                    print("\t\tCould not extract valid expression data for '{}'".format(probe_name))
                    exit()

                # Apply genotype to expression coupling.
                if gte_dict is not None:
                    mask = []
                    new_index = []
                    for index in expression.index:
                        if int(index) in gte_dict.keys():
                            mask.append(True)
                            new_index.append(gte_dict[int(index)])
                        else:
                            mask.append(False)
                    expression = expression.loc[mask, :]
                    expression.index = new_index

                # Merge the data.
                data = genotype.merge(expression, left_index=True,
                                      right_index=True)
                data.columns = ["genotype", "expression"]
                data["group"] = data["genotype"].round(0)

                # Remove missing values.
                data = data.loc[(data['genotype'] >= 0.0) &
                                (data['genotype'] <= 2.0), :]

                # Get the allele data.
                alleles = alleles_df.loc[snp_name, :]
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

                allele_map = {
                    0.0: "{}/{}".format(major_allele, major_allele),
                    1.0: "{}/{}".format(major_allele, minor_allele),
                    2.0: "{}/{}".format(minor_allele, minor_allele)}
                data["alleles"] = data["group"].map(allele_map)

                # Determine the minor allele frequency.
                minor_allele_frequency = min(zero_geno_count,
                                             two_geno_count) / (
                                                 zero_geno_count + two_geno_count)

                # Add the color.
                data["round_geno"] = data["genotype"].round(2)
                data["value_hue"] = data["round_geno"].map(self.value_color_map)
                data["group_hue"] = data["group"].map(self.group_color_map)

                if fdr_value < self.alpha:
                    print("\t\tPlotting: {}\t{}\t{}".format(snp_name, probe_name, hgnc_name))
                    self.plot_simple_eqtl(i, ct, fdr_value, snp_name, probe_name,
                                          hgnc_name, eqtl_type, data,
                                          minor_allele, minor_allele_frequency,
                                          allele_map, ct_outdir)

                # Calculate correlation
                coef, p = stats.spearmanr(data["genotype"], data["expression"])
                lower, upper = self.calc_pearson_ci(coef, data.shape[0])

                # Saving data.
                data_list = interest_corr[hgnc_name]
                data_list.append([minor_allele, ct, coef, lower, upper])
                interest_corr[hgnc_name] = data_list

        # # Combine data.
        # print("Plotting combined eQTL")
        # for name, info in interest_data.items():
        #     print("\t\tPlotting: {}".format(name))
        #     self.plot_combined_eqtl(name, info)

        print("Plotting forest plot")
        for name, info in interest_corr.items():
            print("\t\tPlotting: {}".format(name))
            df = pd.DataFrame(info, columns=["affect allele", "cell type", "mean", "lower", "upper"])
            counts = df["affect allele"].value_counts()
            counts.sort_values(inplace=True)
            majority_allele = counts.index[0]
            if len(counts.index) > 1:
                tmp = df.copy()
                for index, row in tmp.iterrows():
                    if row["affect allele"] != majority_allele:
                        df.at[index, "affect allele"] = majority_allele
                        df.at[index, "mean"] = df.at[index, "mean"] * -1
                        df.at[index, "lower"] = df.at[index, "lower"] * -1
                        df.at[index, "upper"] = df.at[index, "upper"] * -1
                del tmp
            df.set_index("cell type", inplace=True)
            for ct in self.cell_types:
                if ct not in df.index:
                    df.loc[ct, :] = [np.nan, np.nan, np.nan, np.nan]
            df = df.loc[self.cell_types, :]
            df.reset_index(inplace=True, drop=False)
            self.plot_stripplot(df, name)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=None, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def plot_simple_eqtl(self, count, cell_type, fdr_value, snp_name, probe_name, hgnc_name,
                         eqtl_type, df, minor_allele, minor_allele_frequency,
                         allele_map, outdir=None):
        if outdir is None:
            outdir = self.outdir

        # Prepare the figure.
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        title = '{}-SN {} {}-eQTL [FDR = {:.2e}]'.format(cell_type, hgnc_name, eqtl_type, fdr_value)
        subtitle_appendix = 'minor allele frequency {} = {:.2f}'.format(minor_allele, minor_allele_frequency)
        ylabel = '{} ({}) expression'.format(probe_name, hgnc_name)
        xlabel = 'SNP {}'.format(snp_name)

        self.plot_eqtl(df, allele_map, self.group_color_map, ax,
                       title=title, subtitle_appendix=subtitle_appendix,
                       xlabel=xlabel, ylabel=ylabel)

        # Safe the plot.
        for extension in self.extensions:
            filename = "{}_{}_simple_eqtl_{}_{}_{}.{}".format(cell_type,
                                                              count,
                                                              snp_name,
                                                              probe_name,
                                                              hgnc_name,
                                                              extension)
            print("\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(outdir, filename), dpi=300)
        plt.close()

    def plot_combined_eqtl(self, name, info, outdir=None):
        if outdir is None:
            outdir = self.outdir

        # Define the dimensions of the grid.
        ncols = 2
        nrows = math.ceil(len(info) / ncols)

        # Initialize plot.
        sns.set(rc={'figure.figsize': (12*ncols, 9*nrows)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(ncols=ncols,
                                 nrows=nrows)

        row_index = 0
        col_index = 0
        for ct, value in info.items():
            # Prep axes.
            ax = axes[row_index, col_index]
            sns.despine(fig=fig, ax=ax)

            title_prefix = '{}-SN'.format(ct)
            xlabel = 'SNP {}'.format(value["snp_name"])
            ylabel = '{} expression'.format(value["hgnc_name"])
            if value is not None:
                # Plot.
                self.plot_eqtl(value["data"],
                               value["allele_map"],
                               self.group_color_map,
                               ax,
                               title=title_prefix + '    FDR = {:.2e} ({})'.format(value["fdr_value"], self.p_value_to_symbol(value["fdr_value"])),
                               xlabel=xlabel,
                               ylabel=ylabel)
            else:
                ax.set_title('{}    FDR = NA    r = NA'.format(title_prefix),
                             fontsize=22, weight='bold', ha='center',
                             va='bottom',
                             transform=ax.transAxes)
                ax.set_xlabel(xlabel,
                              fontsize=14,
                              fontweight='bold')
                ax.set_ylabel(ylabel,
                              fontsize=14,
                              fontweight='bold')

            # Increment indices.
            col_index += 1
            if col_index == ncols:
                col_index = 0
                row_index += 1

        # Safe the plot.
        for extension in self.extensions:
            filename = "{}_combined_eqtl.{}".format(name, extension)
            print("\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(outdir, filename), dpi=300)
        plt.close()

    @staticmethod
    def plot_eqtl(df, allele_map, group_color_map, ax, title="",
                  subtitle_appendix="", xlabel="", ylabel=""):

        # Calculate the correlation.
        coef, p = stats.spearmanr(df["genotype"],
                                  df["expression"])

        # Plot the scatter / box plot.
        sns.regplot(x="genotype", y="expression", data=df,
                    scatter=False,
                    line_kws={"color": "#000000"},
                    ax=ax
                    )
        sns.boxplot(x="group", y="expression", data=df,
                    palette=group_color_map,
                    showfliers=False,
                    zorder=1,
                    ax=ax)

        # Set the other aesthetics.
        ax.set_xticks(range(3))
        ax.set_xticklabels([allele_map[0.0], allele_map[1.0], allele_map[2.0]])

        if subtitle_appendix == "":
            ax.set_title('{}    r = {:.2f}'.format(title, coef),
                         fontsize=22, weight='bold', ha='center', va='bottom',
                         transform=ax.transAxes)
        else:
            ax.text(0.5, 1.06,
                    title,
                    fontsize=22, weight='bold', ha='center', va='bottom',
                    transform=ax.transAxes)
            ax.text(0.5, 1.02,
                    'r = {:.2f}    {}'.format(coef, subtitle_appendix),
                    fontsize=14, alpha=0.75, ha='center', va='bottom',
                    transform=ax.transAxes)
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

    @staticmethod
    def p_value_to_symbol(p_value):
        output = ""
        try:
            if p_value > 0.05:
                output = "NS"
            elif 0.05 >= p_value > 0.01:
                output = "*"
            elif 0.01 >= p_value > 0.001:
                output = "**"
            elif 0.001 >= p_value > 0.0001:
                output = "***"
            elif p_value <= 0.0001:
                output = "****"
        except TypeError:
            pass

        return output

    @staticmethod
    def calc_pearson_ci(coef, n):
        stderr = 1.0 / math.sqrt(n - 3)
        delta = 1.96 * stderr
        lower = math.tanh(math.atanh(coef) - delta)
        upper = math.tanh(math.atanh(coef) + delta)
        return lower, upper

    def plot_stripplot(self, df, name):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        df_m = df.melt(id_vars=["cell type"], value_vars=["lower", "upper"])
        sns.pointplot(x="value",
                      y="cell type",
                      data=df_m,
                      join=False,
                      palette=self.colormap,
                      ax=ax)

        sns.stripplot(x="mean",
                      y="cell type",
                      data=df,
                      size=25,
                      dodge=False,
                      orient="h",
                      palette=self.colormap,
                      linewidth=1,
                      edgecolor="w",
                      jitter=0,
                      ax=ax)

        ax.set_ylabel('',
                      fontsize=12,
                      fontweight='bold')
        ax.set_xlabel('Spearman r',
                      fontsize=12,
                      fontweight='bold')
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=10)

        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        ax.set_xlim([-1, 1])

        for extension in self.extensions:
            filename = "{}_stripplot.{}".format(name, extension)
            print("\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(self.outdir, filename), dpi=300)
        plt.close()


    def print_arguments(self):
        print("Arguments:")
        print("  > Cell types: {}".format(self.cell_types))
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Genotype to expression coupling path: {}".format(self.gte_path))
        print("  > Alpha: {}".format(self.alpha))
        print("  > Interest: {}".format(self.interest))
        print("  > Nrows: {}".format(self.nrows))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
