#!/usr/bin/env python3

"""
File:         sn_simple_eqtl_plot.py
Created:      2020/11/27
Last Changed: 2020/12/02
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
__program__ = "SN Simple eQTL Plot"
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

./sn_simple_eqtl_plot.py
"""


class main():
    def __init__(self):

        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.expr_file_suffix = getattr(arguments, 'expression_suffix')
        self.eqtl_type = getattr(arguments, 'eqtl_type')
        self.extensions = getattr(arguments, 'extension')
        self.interest = getattr(arguments, 'interest')
        self.plot_cc_corr = getattr(arguments, 'plot_cc_corr')
        self.filter = getattr(arguments, 'filter')

        # Set variables.
        self.palette = {
            "In": "#56B4E9",
            "Ex": "#0072B2",
            "Oli": "#009E73",
            "End": "#CC79A7",
            "Mic": "#E69F00",
            "Ast": "#D55E00",
            "Per": "#808080",
            "Opc": "#F0E442"
        }

        # Set file paths.
        self.eqtl_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/{}_100Perm/ROSMAP-scRNAseq-snpProbe-{}.txt".format(self.eqtl_type, self.eqtl_type)
        self.signif_eqtl_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/{}_100Perm/COMBINED/eQTLProbesFDR0.05-ProbeLevel.txt.gz".format(self.eqtl_type)
        self.geno_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/genotypedump/GenotypeData.txt.gz"
        self.sn_expr_folder = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/"
        self.cc_path = "/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/2020-10-22-MVochteloo-Copy/cell_counts.txt"
        self.gte_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/ROSMAP-scRNAseq-genometoexpressioncoupling.txt"
        self.gene_info_path = "/groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz"

        # Construct output directories.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "sn_simple_eqtl_plot", self.eqtl_type)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.cc_corr_outdir = os.path.join(self.outdir, "cf_corr")
        if self.plot_cc_corr:
            if not os.path.exists(self.cc_corr_outdir):
                os.makedirs(self.cc_corr_outdir)

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
        parser.add_argument("-et",
                            "--eqtl_type",
                            type=str,
                            choices=["cis", "trans"],
                            help="The eQTL type to plot.")
        parser.add_argument("-es",
                            "--expression_suffix",
                            type=str,
                            default="_expression.txt",
                            help="The suffix of the expression file. Default: "
                                 "<cell type>_expression.txt")
        parser.add_argument("-i",
                            "--interest",
                            nargs="+",
                            type=str,
                            default=None,
                            help="The HGNC names to print the info of, "
                                 "default: None.")
        parser.add_argument("-e",
                            "--extension",
                            type=str,
                            choices=["png", "pdf"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")
        parser.add_argument("-plot_cc_corr",
                            action='store_true',
                            help="Plot correlation with cell counts."
                                 " Default: False.")
        parser.add_argument("-filter",
                            action='store_true',
                            help="Filter on only significant eQTLs."
                                 " Default: False.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data")
        eqtl_df = self.load_file(self.eqtl_path, index_col=None, header=None)
        geno_df = self.load_file(self.geno_path)
        cc_df = self.load_file(self.cc_path)
        gte_df = self.load_file(self.gte_path, index_col=None, header=None)
        gene_info_df = self.load_file(self.gene_info_path)

        sn_expr_data = {}
        for sn_expr_path in glob.glob(os.path.join(self.sn_expr_folder, "*{}".format(self.expr_file_suffix))):
            cell_type = os.path.basename(sn_expr_path).replace(self.expr_file_suffix, "")
            sn_expr_data[cell_type] = self.load_file(sn_expr_path)

        signif_eqtl_df = None
        if self.filter:
            signif_eqtl_df = self.load_file(self.signif_eqtl_path, index_col=None)

        print("### Step3 ###")
        print("Preprocess")
        cc_df = cc_df / cc_df.sum(axis=0)
        cc_df.index = [x.upper() for x in cc_df.index]

        gte_dict = dict(zip(gte_df.iloc[:, 0], gte_df.iloc[:, 1]))
        del gte_df

        gene_dict = dict(zip(gene_info_df["ArrayAddress"], gene_info_df["Symbol"]))
        del gene_info_df

        print("### Step3 ###")
        print("Filter data")
        alleles_df = geno_df[["Alleles", "MinorAllele"]].copy()

        mask = []
        colnames = []
        for col in geno_df.columns:
            if col in gte_dict.keys():
                mask.append(True)
                colnames.append(str(gte_dict[col]))
            else:
                mask.append(False)

        geno_df = geno_df.loc[:, mask].copy()
        geno_df.columns = colnames
        geno_df.index.name = "-"

        print("### Step4 ###")
        print("Visualise snp vs cc")
        plotted_snps = set()
        corr_data = []
        for i, (_, (snp_name, probe_name)) in enumerate(eqtl_df.iterrows()):
            hgnc_name = gene_dict[probe_name]

            if (self.interest is not None and hgnc_name not in self.interest) or (signif_eqtl_df is not None and signif_eqtl_df.loc[(signif_eqtl_df["SNPName"] == snp_name) & (signif_eqtl_df["ProbeName"] == probe_name) & (signif_eqtl_df["HGNCName"] == hgnc_name), :].shape[0] == 0):
                continue

            eqtl_outdir = os.path.join(self.outdir, "{}_{}_{}_{}".format(i, snp_name, probe_name, hgnc_name))
            if not os.path.exists(eqtl_outdir):
                os.makedirs(eqtl_outdir)

            print("\tWorking on: {}\t{}\t{}".format(snp_name, probe_name, hgnc_name))

            # Get the genotype data.
            genotype = geno_df.loc[[snp_name], :].copy()
            genotype = genotype.T
            genotype = genotype.loc[(genotype[snp_name] >= 0.0) &
                                    (genotype[snp_name] <= 2.0), :]

            # Get the allele data.
            alleles = alleles_df.loc[snp_name, :].copy()
            major_allele = alleles["Alleles"].split("/")[0]
            minor_allele = alleles["Alleles"].split("/")[1]
            allele_map = {0.0: "{}/{}".format(major_allele, major_allele),
                          1.0: "{}/{}".format(major_allele, minor_allele),
                          2.0: "{}/{}".format(minor_allele, minor_allele)}

            for ct, expr_df in sn_expr_data.items():
                if signif_eqtl_df is not None and signif_eqtl_df.loc[(signif_eqtl_df["SNPName"] == snp_name) & (signif_eqtl_df["ProbeName"] == probe_name) & (signif_eqtl_df["HGNCName"] == hgnc_name) & (signif_eqtl_df["CellType"] == ct), :].shape[0] == 0:
                    continue

                # Get the expression data.
                if probe_name not in expr_df.index:
                    print("\t\tCell type {} does not have expression for probe {}".format(ct, probe_name))
                    continue
                expression = expr_df.loc[[probe_name], :].copy()
                data = expression.T.merge(genotype.copy(), left_index=True, right_index=True)

                self.plot(df=data, x=snp_name, y=probe_name, alleles=allele_map,
                          xlabel="SNP {}".format(snp_name),
                          ylabel="{} ({}) expression".format(probe_name,
                                                             hgnc_name),
                          title="{} {}-eQTL [Data: {}]".format(hgnc_name, self.eqtl_type, ct),
                          file_suffix="_{}data".format(ct),
                          outdir=eqtl_outdir)

                if self.plot_cc_corr and "{}_{}".format(snp_name, ct) not in plotted_snps:
                    data = data.merge(cc_df.loc[[ct], :].T, left_index=True, right_index=True)

                    self.plot(df=data, x=snp_name, y=ct, alleles=allele_map,
                              xlabel="SNP {}".format(snp_name),
                              ylabel="{} cell proportion".format(ct),
                              title="{} vs {} correlation".format(snp_name, ct),
                              file_suffix="_cellCounts",
                              outdir=self.cc_corr_outdir)

                    plotted_snps.add("{}_{}".format(snp_name, ct))

                    coef, p = stats.spearmanr(data[snp_name], data[ct])
                    corr_data.append([snp_name, ct, coef, p])
        corr_df = pd.DataFrame(corr_data, columns = ["SNPName", "CellType", "Coef", "p"])
        corr_df["AbsCoef"] = corr_df["Coef"].abs()
        corr_df.sort_values(by="AbsCoef", ascending=False, inplace=True)
        print(corr_df)
        print(corr_df["AbsCoef"].mean())
        corr_df.to_csv(os.path.join(self.cc_corr_outdir, "correlations.txt.gz"), compression="gzip", header=True, index=False, sep="\t")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def plot(self, df, x, y, alleles, xlabel="", ylabel="", title="", file_suffix="", outdir=None):
        if outdir is None:
            outdir = self.outdir

        # Calculate the correlation.
        coef, p = stats.spearmanr(df[x], df[y])

        df["round_{}".format(x)] = df[x].round(2)
        counts = df["round_{}".format(x)].value_counts()
        for genotype in [0.0, 1.0, 2.0]:
            if genotype not in counts:
                counts[genotype] = 0

        # Prepare the figure.
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        # Plot the scatter / box plot.
        sns.regplot(x=x, y=y, data=df,
                    # scatter=False,
                    scatter_kws={"color": "#808080"},
                    line_kws={"color": "#000000"},
                    ax=ax
                    )
        sns.boxplot(x="round_{}".format(x), y=y, data=df,
                    palette={0.0: "#D55E00", 1.0: "#0072B2", 2.0: "#E69F00"},
                    showfliers=False,
                    zorder=1,
                    # boxprops=dict(alpha=.3),
                    ax=ax)

        ax.set_xticks(range(3))
        ax.set_xticklabels([alleles[0.0], alleles[1.0], alleles[2.0]])

        # Set the other aesthetics.
        ax.text(0.5, 1.06,
                title,
                fontsize=22, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                'r = {:.2f} [p = {:.2e}]    N = {}:{}, {}:{}, {}:{}'.format(coef, p, alleles[0.0], counts[0.0], alleles[1.0], counts[1.0], alleles[2.0], counts[2.0]),
                fontsize=14, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        # Safe the plot.
        for extension in self.extensions:
            filename = "sn_simple_eqtl_{}_{}{}.{}".format(x, y, file_suffix, extension)
            print("\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(outdir, filename))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL type: {}".format(self.eqtl_type))
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression folder: {}".format(self.sn_expr_folder))
        print("  > Expression file suffix: {}".format(self.expr_file_suffix))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > GTE coupling path: {}".format(self.gte_path))
        print("  > Gene info path: {}".format(self.gene_info_path))
        print("  > Interest: {}".format(self.interest))
        print("  > Extension: {}".format(self.extensions))
        print("  > Plot cell count: {}".format(self.plot_cc_corr))
        print("  > Filter FDR <0.05: {}".format(self.filter))
        print("  > Output directory: {}".format(self.outdir))
        if self.plot_cc_corr:
            print("  > Cell count correlations directory: {}".format(self.cc_corr_outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
