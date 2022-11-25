#!/usr/bin/env python3

"""
File:         visualise_inter_eqtl.py
Created:      2020/10/27
Last Changed:
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
import gzip
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
__program__ = "Visualise Inter eQTL"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.geno_path = getattr(arguments, 'genotype')
        self.alleles_path = getattr(arguments, 'alleles')
        self.expr_path = getattr(arguments, 'expression')
        self.cov_path = getattr(arguments, 'covariates')
        self.decon_path = getattr(arguments, 'decon')
        self.alpha = getattr(arguments, 'alpha')
        self.interest = getattr(arguments, 'interest')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'inter_eQTLs')

        self.colormap = {
            "minor": "#E69F00",
            "center": "#0072B2",
            "major": "#D55E00",
            "low": 246,
            "hig": 24,
            "male": "#56B4E9",
            "female": "#CC79A7",
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
        self.group_color_map, self.value_color_map = self.create_color_map(
            self.colormap)
        self.sex_color_map = {"Male": self.colormap["male"],
                              "Female": self.colormap["female"]}

    def create_argument_parser(self):
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
                            help="The path to the expression matrix")
        parser.add_argument("-c",
                            "--covariates",
                            type=str,
                            required=True,
                            help="The path to the covariates matrix")
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
                            required=True,
                            help="The covariates to plot")
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

        print("### Step 1 ###")
        print("Loading eQTL data.")
        eqtl_df = self.load_file(self.eqtl_path, index_col=None)

        print("### Step 2 ###")
        print("Loading covariate / deconvolution / eQTL data.")
        cov_indices, cov_df = self.search_file(self.cov_path,
                                               interest=self.interest)
        decon_indices, decon_df = self.search_file(self.decon_path,
                                                   interest=self.interest)

        print("### Step 3 ###")
        print("Preprocessing data.")
        decon_df = decon_df.T
        snp_names = []
        probe_names = []
        for index in decon_df.index:
            snp_names.append("_".join(index.split("_")[:-1]))
            probe_names.append(index.split("_")[-1])
        decon_df["SNPName"] = snp_names
        decon_df["ProbeName"] = probe_names
        decon_df.reset_index(inplace=True, drop=True)

        if cov_indices != decon_indices:
            print("Covariate and deconvolution results do not match.")
            exit()

        print("### Step 4 ###")
        print("Combine deconvolution and eQTL data frames.")
        full_eqtl_df = eqtl_df.merge(decon_df, left_on=["SNPName", "ProbeName"],
                                     right_on=["SNPName", "ProbeName"])
        print("\tSHape: {}".format(full_eqtl_df.shape))

        print("### Step 5 ###")
        print("Loading other data.")
        geno_df = self.load_file(self.geno_path, nrows=full_eqtl_df.shape[0])
        alleles_df = self.load_file(self.alleles_path, nrows=full_eqtl_df.shape[0])
        expr_df = self.load_file(self.expr_path, nrows=full_eqtl_df.shape[0])

        print("### Step 2 ###")
        self.plot(full_eqtl_df, geno_df, alleles_df, expr_df, cov_df)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def search_file(path, sep="\t", header_index=0, index_col=0, nrows=None,
                    interest=None, print_interval=500):
        columns = []
        indices_int = []
        indices_str = []
        data = []

        if interest is None:
            return None
        else:
            search_list = set(interest)

        print("Searching file.")
        with gzip.open(path, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % print_interval == 0):
                    print("\tprocessed {} lines".format(i))
                if len(search_list) == 0:
                    break

                splitted_line = line.decode().strip('\n').split(sep)
                index = None
                if index_col is not None:
                    index = splitted_line[index_col]

                data_start_index = 0
                if index_col is not None:
                    data_start_index = index_col + 1

                content = splitted_line[data_start_index:]

                if header_index is not None and i == header_index:
                    columns = content
                else:
                    if interest is None or index in interest:
                        indices_int.append(i)
                        indices_str.append(index)
                        data.append(np.array(content, dtype=float))

                        search_list.remove(index)

                if nrows is not None and i > nrows:
                    break

        f.close()

        df = pd.DataFrame(data, columns=columns, index=indices_str)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return indices_int, df

    def plot(self, eqtl_df, geno_df, alleles_df, expr_df, cc_df):
        print("Plotting interaction eQTL plots.")

        print("Iterating over interest covariates.")
        for i, interest in enumerate(self.interest):
            print("\tWorking on: {}\t{}/{} [{:.2f}%]".format(interest, i, len(interest), (100/len(interest)*i)))

            # Get the order.
            interest_pvalues = eqtl_df[[interest]].copy()
            interest_pvalues.sort_values(by=interest, ascending=True, inplace=True)

            # Prepare output directory.
            interest_outdir = os.path.join(self.outdir, interest)
            if not os.path.exists(interest_outdir):
                os.makedirs(interest_outdir)

            print("\tIterating over eQTLs.")
            count = 0
            for index, (inter_pvalue, ) in interest_pvalues.iterrows():
                if inter_pvalue >= self.alpha:
                    continue
                if count > 5:
                    break

                # Extract the usefull information from the row.
                row = eqtl_df.iloc[index, :]
                p_value = row["PValue"]
                snp_name = row["SNPName"]
                probe_name = row["ProbeName"]
                hgnc_name = row["HGNCName"]
                eqtl_type = row["CisTrans"]

                print("\t\tWorking on: {}\t{}\t{}\t{}".format(index, snp_name, probe_name, hgnc_name))

                # Get the genotype / expression data.
                genotype = geno_df.iloc[index, :].T.to_frame()
                if genotype.columns != [snp_name]:
                    print("\t\t\tGenotype file not in identical order as eQTL file.")
                    exit()
                expression = expr_df.iloc[index, :].T.to_frame()
                if expression.columns != [probe_name]:
                    print("\t\t\tExpression file not in identical order as eQTL file.")
                    exit()
                data = genotype.merge(expression, left_index=True, right_index=True)
                data.columns = ["genotype", "expression"]
                data["group"] = data["genotype"].round(0)

                # Remove missing values.
                data = data.loc[(data['genotype'] >= 0.0) &
                                (data['genotype'] <= 2.0), :]

                # Get the allele data.
                alleles = alleles_df.iloc[index, :]
                if alleles.name != snp_name:
                    print("\t\t\tAlleles file not in identical order as eQTL file.")
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

                interest_eqtl_outdir = os.path.join(interest_outdir,
                                                    "{}_{}_{}_{}_{}".format(
                                                           count, index, snp_name,
                                                           probe_name, hgnc_name))

                if not os.path.exists(interest_eqtl_outdir):
                    os.makedirs(interest_eqtl_outdir)

                self.plot_simple_eqtl(index, p_value, snp_name, probe_name,
                                      hgnc_name, eqtl_type, data,
                                      minor_allele, minor_allele_frequency,
                                      allele_map, self.group_color_map,
                                      interest_eqtl_outdir,
                                      self.extensions)


                eqtl_data = data.copy()
                cov_data = cc_df.loc[interest, :].to_frame()
                eqtl_data = eqtl_data.merge(cov_data,
                                            left_index=True,
                                            right_index=True)

                self.plot_inter_eqtl(index, snp_name, probe_name, hgnc_name,
                                     eqtl_type, eqtl_data, interest, inter_pvalue,
                                     allele_map, self.group_color_map,
                                     interest_eqtl_outdir,
                                     self.extensions)

                count += 1

    @staticmethod
    def plot_simple_eqtl(count, p_value, snp_name, probe_name, hgnc_name,
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
        ax.text(0.5, 1.06,
                '{} {}-eQTL [p = {:.2e}]'.format(hgnc_name, eqtl_type, p_value),
                fontsize=22, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
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
            print("\t\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(outdir, filename), dpi=300)
        plt.close()


    @staticmethod
    def plot_inter_eqtl(count, snp_name, probe_name, hgnc_name, eqtl_type, df,
                        cov_name, pvalue, allele_map, group_color_map, outdir,
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
                p_str = "p = {:.2e}".format(p)

                # Plot.
                sns.regplot(x=cov_name, y="expression", data=subset,
                            scatter_kws={'facecolors': subset['value_hue'],
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": color, "alpha": 0.75},
                            ax=ax
                            )

            # Add the text.
            ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
            ax.annotate(
                '{}: r = {} [{}]'.format(allele, coef_str, p_str),
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
                'P-value: {:.2e}'.format(snp_name, probe_name, pvalue),
                fontsize=12, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel('{} ({}) expression'.format(probe_name, hgnc_name),
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(cov_name,
                      fontsize=14,
                      fontweight='bold')

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
            print("\t\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(outdir, filename), dpi=300)
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Alleles path: {}".format(self.alleles_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Covariates path: {}".format(self.cov_path))
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Alpha: {}".format(self.alpha))
        print("  > Interest: {}".format(self.interest))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
