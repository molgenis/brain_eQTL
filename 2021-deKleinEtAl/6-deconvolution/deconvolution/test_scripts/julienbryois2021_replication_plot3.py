#!/usr/bin/env python3

"""
File:         julienbryois2021_replication_plot3.py
Created:      2022/01/24
Last Changed: 2022/02/11
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
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from statsmodels.stats import multitest
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

"""
Syntax:
./julienbryois2021_replication_plot3.py -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-19-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron-ExprMatrixCovariatesRemovedOLS-NegativeToZero/deconvolutionResults.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_alleles.txt.gz -dn MetaBrain_Decon-eQTL -o 2022-01-19-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron-ExprMatrixCovariatesRemovedOLS-NegativeToZero

./julienbryois2021_replication_plot3.py \
    -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron/deconvolutionResults_BH_FDR.txt.gz \
    -dg /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron/geno_stats.txt.gz \
    -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_alleles.txt.gz \
    -dn MetaBrain_Decon-eQTL_BH \
    -o 2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron-BH_FDR

"""

# Metadata
__program__ = "Julien Bryois et al. 2021 Replication Plot 3"
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
        self.discovery_data_path = getattr(arguments, 'discovery_data_path')
        self.discovery_genotype_stats_path = getattr(arguments, 'discovery_genotype_stats_path')
        self.discovery_alleles_path = getattr(arguments, 'discovery_alleles_path')
        self.discovery_name = getattr(arguments, 'discovery_name')
        self.out_filename = getattr(arguments, 'outfile')
        self.extensions = getattr(arguments, 'extension')

        self.replication_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/julienbryois2021"
        self.replication_name = "JulienBryois2021"

        # Set variables.
        outdir = os.path.join(str(Path(__file__).parent.parent), 'julienbryois2021_replication_plot')
        self.outdir = os.path.join(outdir, 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.colormap = {
            "Excitatory": "#56B4E9",
            "ExcitatoryNeurons": "#56B4E9",
            "Inhibitory": "#0072B2",
            "InhibitoryNeurons": "#0072B2",
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "Oligodendrocytes": "#009E73",
            "OPC": "#009E73",
            "OPCs...COPs": "#009E73",
            "OPCsCOPs": "#009E73",
            "EndothelialCell": "#CC79A7",
            "EndothelialCells": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Astrocytes": "#D55E00",
            "Pericytes": "#808080",
            "OtherNeuron": "#2690ce"
        }

        self.bryois_ct_trans = {
            'Astrocytes': 'Astrocytes',
            'Endothelial cells': 'EndothelialCells',
            'Excitatory neurons': 'ExcitatoryNeurons',
            'Inhibitory neurons': 'InhibitoryNeurons',
            'Microglia': 'Microglia',
            'OPCs / COPs': 'Oligodendrocytes',
            'Oligodendrocytes': 'OPCsCOPs',
            'Pericytes': 'Pericytes'
        }

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
        parser.add_argument("-dd",
                            "--discovery_data_path",
                            type=str,
                            required=True,
                            help="The path to the discovery deconvolution "
                                 "results matrix")
        parser.add_argument("-dg",
                            "--discovery_genotype_stats_path",
                            type=str,
                            required=True,
                            help="The path to the discovery genotypes matrix")
        parser.add_argument("-da",
                            "--discovery_alleles_path",
                            type=str,
                            required=True,
                            help="The path to the discovery alleles matrix")
        parser.add_argument("-dn",
                            "--discovery_name",
                            type=str,
                            default="one",
                            help="The name for the discovery data")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the outfile.")
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
        discovery_df = self.load_file(self.discovery_data_path, header=0, index_col=None)
        discovery_df["SNPName"] = ["_".join(x.split("_")[1:]) for x in discovery_df.iloc[:, 0]]
        discovery_cell_types = [x.replace("_pvalue", "") for x in discovery_df if x.endswith("_pvalue")]

        discovery_genotype_stats_df = self.load_file(self.discovery_genotype_stats_path, header=0, index_col=0)
        maf_dict = dict(zip(discovery_genotype_stats_df.index, discovery_genotype_stats_df["MAF"]))
        maf_dict2 = dict(zip([x.split(":")[2] for x in discovery_genotype_stats_df.index], discovery_genotype_stats_df["MAF"]))
        n_dict = dict(zip(discovery_genotype_stats_df.index, discovery_genotype_stats_df["N"]))
        discovery_df["MAF"] = discovery_df["SNPName"].map(maf_dict)
        discovery_df["N"] = discovery_df["SNPName"].map(n_dict)

        discovery_alleles_df = self.load_file(self.discovery_alleles_path, header=0, index_col=None)
        discovery_alleles_df.columns = ["SNPName", "Alleles", "MinorAllele"]
        discovery_alleles_df["DeconAllele"] = [x.split("/")[1] for x in discovery_alleles_df["Alleles"]]
        snp_to_da_dict = dict(zip(discovery_alleles_df["SNPName"], discovery_alleles_df["DeconAllele"]))

        discovery_df["DeconAllele"] = discovery_df["SNPName"].map(snp_to_da_dict)
        discovery_df.index = ["{}_{}".format(x.split(".")[0], x.split(":")[2]) for x in discovery_df.iloc[:, 0]]
        discovery_df = discovery_df.iloc[:, 1:]
        discovery_df.reset_index(drop=False, inplace=True)

        discovery_fdr_df = discovery_df.loc[:, [x for x in discovery_df.columns if x.endswith("_FDR")]]
        discovery_fdr_df = (discovery_fdr_df <= 0.05).astype(int)
        discovery_unique_ieqtls_df = discovery_fdr_df.loc[discovery_fdr_df.sum(axis=1) == 1, :]

        replication_df = self.load_file(os.path.join(self.replication_path, "JulienBryois2021SummaryStats.txt.gz"), header=0, index_col=0, nrows=None)
        replication_cell_types = [x.replace(" p-value", "") for x in replication_df if x.endswith(" p-value")]
        replication_df.reset_index(drop=False, inplace=True)
        replication_df["N"] = 196
        replication_df["MAF"] = replication_df["SNP"].map(maf_dict2)

        # replication_signif_df = self.load_file(os.path.join(self.replication_path, "media-3-filtered.txt.gz"), header=0, index_col=0, nrows=None)
        # replication_signif_df["cell_type"] = replication_signif_df["cell_type"].map(self.bryois_ct_trans)
        # replication_signif_df.reset_index(drop=False, inplace=True)
        #
        # count_bryois_eqtls = dict(zip(*np.unique(replication_signif_df["cell_type"], return_counts=True)))

        print("Calculating replication stats")
        indices = []
        # columns = ["{}\n[N={:,}]".format(replication_ct, count_bryois_eqtls[replication_ct]) for replication_ct in replication_cell_types]
        columns = replication_cell_types
        color_data = []
        annot_data = []
        for discovery_ct in discovery_cell_types:
            discovery_pvalue_column = "{}_pvalue".format(discovery_ct)
            discovery_fdr_column = "{}_FDR".format(discovery_ct)
            discovery_beta_column = None
            for col in discovery_df:
                if discovery_ct in col and col.endswith(":GT"):
                    discovery_beta_column = col
                    break

            # Select the discovery ieQTLs.
            ieqtls = list(discovery_unique_ieqtls_df.loc[discovery_unique_ieqtls_df[discovery_fdr_column] == 1, :].index)
            # ieqtls = list(discovery_df.loc[discovery_df[discovery_fdr_column] <= 0.05, :].index)
            discovery_subset_df = discovery_df.loc[ieqtls, ["index", discovery_pvalue_column, discovery_beta_column, "DeconAllele", "N", "MAF"]].copy()
            discovery_subset_df.columns = ["index", "discovery p-value", "discovery beta", "discovery allele", "discovery n", "discovery MAF"]

            self.pvalue_to_zscore(df=discovery_subset_df,
                                  beta_col="discovery beta",
                                  p_col="discovery p-value",
                                  prefix="discovery ")
            self.zscore_to_beta(df=discovery_subset_df,
                                z_col="discovery z-score",
                                maf_col="discovery MAF",
                                n_col="discovery n",
                                prefix="discovery zscore-to-beta ")

            # Save index.
            indices.append("{}\n[N={:,}]".format(discovery_ct, len(ieqtls)))

            ct_color_data = [np.nan] * len(replication_cell_types)
            ct_annot_data = [""] * len(replication_cell_types)

            if len(ieqtls) > 0:
                ct_color_data = []
                ct_annot_data = []
                for replication_ct in replication_cell_types:
                    rb = np.nan
                    annot_str = ""

                    # Select the replication eQTLs.
                    replication_subset_df = replication_df.loc[:, ["index", "{} p-value".format(replication_ct), "{} beta".format(replication_ct), "effect_allele", "N", "MAF"]].copy()
                    replication_subset_df.columns = ["index", "replication p-value", "replication beta", "replication allele", "replication n", "replication MAF"]
                    replication_subset_df.dropna(inplace=True)

                    if replication_subset_df.shape[0] > 0:
                        # Merge together.
                        subset_df = discovery_subset_df.merge(replication_subset_df, on=["index"], how="inner")
                        n_overlap = subset_df.shape[0]

                        # Flip direction.
                        subset_df["flip"] = (subset_df["discovery allele"] == subset_df["replication allele"]).replace({False: -1, True: 1})
                        subset_df["replication beta"] = subset_df["replication beta"] * subset_df["flip"]

                        # Calculate se.
                        self.pvalue_to_zscore(df=subset_df,
                                              beta_col="replication beta",
                                              p_col="replication p-value",
                                              prefix="replication ")
                        self.zscore_to_beta(df=subset_df,
                                            z_col="replication z-score",
                                            maf_col="replication MAF",
                                            n_col="replication n",
                                            prefix="replication zscore-to-beta ")

                        # Calculate the p1.
                        subset_df.sort_values(by="discovery p-value", inplace=True)
                        pi1 = self.calculate_p1(p=subset_df["replication p-value"])

                        # Calculate the rb.
                        rb_est = self.calculate_rb(b1=subset_df["discovery zscore-to-beta beta"],
                                                   se1=subset_df["discovery zscore-to-beta se"],
                                                   b2=subset_df["replication zscore-to-beta beta"],
                                                   se2=subset_df["replication zscore-to-beta se"],
                                                   )
                        rb = rb_est[0]

                        # Filter on significant eQTLs.
                        # subset_df = subset_df.loc[subset_df["replication adjusted p-value"] <= 0.05, :]
                        subset_df["replication FDR"] = multitest.multipletests(subset_df["replication p-value"], method='fdr_bh')[1]
                        subset_df = subset_df.loc[subset_df["replication FDR"] <= 0.05, :]
                        n_replicating = subset_df.shape[0]

                        concordance = np.nan
                        replication_rate = 0
                        if n_replicating > 0:
                            # Calculate concordance.
                            lower_quadrant = subset_df.loc[(subset_df["discovery beta"] < 0) & (subset_df["replication beta"] < 0), :]
                            upper_quadrant = subset_df.loc[(subset_df["discovery beta"] > 0) & (subset_df["replication beta"] > 0), :]
                            concordance = (100 / subset_df.shape[0]) * (lower_quadrant.shape[0] + upper_quadrant.shape[0])

                            # Calculate replication rate.
                            replication_rate = (100 / n_overlap) * n_replicating

                        # Construc the annotation string.
                        annot_str = "Overlap {:,}\n\u03C01 {:.2f}\nRb: {:.2f}\nRepl. {:,} [{:.0f}%]\nConc. {:.0f}%".format(n_overlap, pi1, rb, n_replicating, replication_rate, concordance)

                    # Save.
                    ct_color_data.append(rb)
                    ct_annot_data.append(annot_str)

            # Save.
            color_data.append(ct_color_data)
            annot_data.append(ct_annot_data)

        color_df = pd.DataFrame(color_data, index=indices, columns=columns)
        annot_df = pd.DataFrame(annot_data, index=indices, columns=columns)
        print(color_df)
        print(annot_df)

        # Plot.
        self.plot_heatmap(df=color_df,
                          annot_df=annot_df)

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
    def pvalue_to_zscore(df, beta_col, p_col, prefix=""):
        p_values = df[p_col].to_numpy()
        zscores = stats.norm.ppf(p_values / 2)
        mask = np.ones_like(p_values)
        mask[df[beta_col] > 0] = -1
        df["{}z-score".format(prefix)] = zscores * mask
        df.loc[df[p_col] == 1, "{}z-score".format(prefix)] = 0
        df.loc[df[p_col] == 0, "{}z-score".format(prefix)] = -40.

    @staticmethod
    def zscore_to_beta(df, z_col, maf_col, n_col, prefix=""):
        chi = df[z_col] * df[z_col]
        a = 2 * df[maf_col] * (1 - df[maf_col]) * (df[n_col] + chi)
        df["{}beta".format(prefix)] = df[z_col] / a ** (1/2)
        df["{}se".format(prefix)] = 1 / a ** (1/2)

    @staticmethod
    def zscore_to_maf(df, z_col, beta_col, n_col, prefix=""):
        df["{}MAF".format(prefix)] = (2 * df[beta_col]**2 - np.sqrt(4 * df[beta_col]**4 - (8 * df[beta_col]**2 * df[z_col]**2)/(df[n_col] + df[z_col]**2))) / (4 * df[beta_col]**2)

    @staticmethod
    def calculate_p1(p):
        importr("qvalue")
        pvals = robjects.FloatVector(p)
        lambda_seq = robjects.FloatVector([x for x in np.arange(0.05, 1, 0.05) if p.max() > x])
        pi0est = robjects.r['pi0est'](pvals, lambda_seq)
        return 1 - np.array(pi0est.rx2('pi0'))[0]

    @staticmethod
    def calculate_rb(b1, se1, b2, se2, theta=0):
        robjects.r("source('Rb.R')")
        b1 = robjects.FloatVector(b1)
        se1 = robjects.FloatVector(se1)
        b2 = robjects.FloatVector(b2)
        se2 = robjects.FloatVector(se2)
        calcu_cor_true = robjects.globalenv['calcu_cor_true']
        rb = calcu_cor_true(b1, se1, b2, se2, theta)
        return np.array(rb)[0]

    def plot_heatmap(self, df, annot_df, xlabel="", ylabel="", appendix=""):
        fig, axes = plt.subplots(nrows=2,
                                 ncols=2,
                                 figsize=(1 * df.shape[1] + 10, 1 * df.shape[0] + 10),
                                 gridspec_kw={"width_ratios": [0.2, 0.8],
                                              "height_ratios": [0.8, 0.2]})
        sns.set(color_codes=True)

        annot_df.fillna("", inplace=True)

        row_index = 0
        col_index = 0
        for _ in range(4):
            ax = axes[row_index, col_index]
            if row_index == 0 and col_index == 1:
                sns.heatmap(df, cmap="Blues",
                            square=True, annot=annot_df, fmt='',
                            cbar=True, annot_kws={"size": 10},
                            ax=ax)

                plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
                plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=45))

                ax.set_xlabel(xlabel, fontsize=14)
                ax.xaxis.set_label_position('top')

                ax.set_ylabel(ylabel, fontsize=14)
                ax.yaxis.set_label_position('right')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > 1:
                col_index = 0
                row_index += 1

        plt.tight_layout()
        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}_corr_heatmap{}.{}".format(self.out_filename, appendix, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Discovery:")
        print("    > Data path: {}".format(self.discovery_data_path))
        print("    > Genotype stats path: {}".format(self.discovery_genotype_stats_path))
        print("    > Alleles path: {}".format(self.discovery_alleles_path))
        print("    > Name: {}".format(self.discovery_name))
        print("  > Replication:")
        print("    > File path: {}".format(self.replication_path))
        print("    > Name: {}".format(self.replication_name))
        print("  > Output filename: {}".format(self.out_filename))
        print("  > Outpath {}".format(self.outdir))
        print("  > Extensions: {}".format(self.extensions))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
