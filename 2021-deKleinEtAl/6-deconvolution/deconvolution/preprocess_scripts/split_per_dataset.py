#!/usr/bin/env python3

"""
File:         split_per_dataset.py
Created:      2022/01/17
Last Changed: 2022/02/10
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
import argparse
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

"""
Syntax:
./split_per_dataset.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/combine_gte_files/SampleToDataset.txt.gz \
    -of 2021-12-07-CortexEUR-cis-InhibitorySummedWithOtherNeuron
    
./split_per_dataset.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -of 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron
"""

# Metadata
__program__ = "Split per Dataset"
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
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), "split_per_dataset", self.outfolder)
        self.plot_outdir = os.path.join(self.outdir, "plot")
        for outdir in [self.outdir, self.plot_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        self.decon_basedir = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/"
        self.job_outdir = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/jobs"
        self.job_output_outdir = os.path.join(self.job_outdir, "output")
        self.time = "05:55:00"
        self.python_executable = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl.py"
        self.genotype_na = -1
        self.allele_configs = "limited"
        self.call_rate = 0.95
        self.hw_pval = 1e-4
        self.maf = 0.05

        for outdir in [self.job_outdir, self.job_output_outdir]:
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
            "Brainseq": "#C778A6",
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "OPC": "#009E73",
            "OPCs...COPs": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Pericytes": "#808080",
            "OtherNeuron": "#2690ce"
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
                            help="show program's version number and exit.")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
        parser.add_argument("-cc",
                            "--cell_counts",
                            type=str,
                            required=True,
                            help="The path to the cell counts matrix.")
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

        print("Loading data.")
        geno_df = self.load_file(self.geno_path, header=0, index_col=0)
        expr_df = self.load_file(self.expr_path, header=0, index_col=0)
        cc_df = self.load_file(self.cc_path, header=0, index_col=0)
        std_df = self.load_file(self.std_path, header=0, index_col=None)

        print("Validating input.")
        self.validate_data(std_df=std_df,
                           geno_df=geno_df,
                           expr_df=expr_df,
                           cc_df=cc_df)

        print("Subsetting dataset")
        cc_dfm_list = []
        for dataset in std_df.iloc[:, 1].unique():
            mask = (std_df.iloc[:, 1] == dataset).to_numpy(dtype=bool)
            print("\tSubsetting {} with {} samples.".format(dataset, np.sum(mask)))

            dataset_outdir = os.path.join(self.outdir, dataset)
            if not os.path.exists(dataset_outdir):
                os.makedirs(dataset_outdir)

            geno_path = os.path.join(dataset_outdir, "genotype_table.txt.gz")
            expr_path = os.path.join(dataset_outdir, "expression_table.txt.gz")
            cc_path = os.path.join(dataset_outdir, "cell_count_table.txt.gz")

            self.save_file(df=geno_df.loc[:, mask], outpath=geno_path)
            self.save_file(df=expr_df.loc[:, mask], outpath=expr_path)
            self.save_file(df=cc_df.loc[mask, :], outpath=cc_path)

            corr_df = self.correlate(df=cc_df.loc[mask, :])
            self.plot_heatmap(df=corr_df,
                              annot_df=corr_df.round(2),
                              xlabel=dataset,
                              ylabel=dataset,
                              name=dataset)

            if dataset in ["AMPAD-ROSMAP-V2", "CMC"]:
                cc_dfm = cc_df.loc[mask, :].melt()
                cc_dfm["dataset"] = dataset
                cc_dfm_list.append(cc_dfm)

            self.create_job_file(job_name=self.outfolder,
                                 dataset=dataset,
                                 geno_path=geno_path,
                                 expr_path=expr_path,
                                 cc_path=cc_path)

            print("")

        print("Print cell fraction violin plot")
        cc_dfm = pd.concat(cc_dfm_list, axis=0)
        print(cc_dfm)

        self.plot_boxplot(df=cc_dfm,
                          x="variable",
                          y="value",
                          hue="dataset",
                          palette=self.palette,
                          name="cell_count")

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
    def validate_data(std_df, geno_df=None, expr_df=None, cc_df=None):
        samples = std_df.iloc[:, 0].values.tolist()
        if geno_df is not None and geno_df.columns.tolist() != samples:
                print("\tThe genotype file header does not match "
                      "the sample-to-dataset link file")
                exit()

        if expr_df is not None and expr_df.columns.tolist() != samples:
                print("\tThe expression file header does not match "
                      "the sample-to-dataset link file")
                exit()

        if cc_df is not None and cc_df.index.tolist() != samples:
                print("\tThe cell count file index does not match "
                      "the sample-to-dataset link file")
                exit()

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath), df.shape))

    @staticmethod
    def correlate(df):
        corr_df = pd.DataFrame(np.nan, index=df.columns, columns=df.columns)

        for i, x_col in enumerate(df.columns):
            for j, y_col in enumerate(df.columns):
                if i < j:
                    continue
                corr_data = df.loc[:, [x_col, y_col]].copy()
                corr_data.dropna(inplace=True)

                coef = np.nan
                if np.min(corr_data.std(axis=0)) > 0:
                    coef, _ = stats.pearsonr(corr_data.iloc[:, 1], corr_data.iloc[:, 0])

                corr_df.loc[x_col, y_col] = coef

        return corr_df

    @staticmethod
    def mask_non_significant(df, pvalue_df, a=0.05):
        signif_df = df.copy()
        for i in range(signif_df.shape[0]):
            for j in range(signif_df.shape[1]):
                if np.isnan(pvalue_df.iloc[i, j]) or pvalue_df.iloc[i, j] > a:
                    signif_df.iloc[i, j] = 0

        return signif_df

    def plot_heatmap(self, df, annot_df, xlabel="", ylabel="", vmin=-1, vmax=1,
                     name=""):
        cmap = sns.diverging_palette(246, 24, as_cmap=True)

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
                sns.heatmap(df, cmap=cmap, vmin=vmin, vmax=vmax, center=0,
                            square=True, annot=annot_df, fmt='',
                            cbar=False, annot_kws={"size": 14, "color": "#000000"},
                            ax=ax)

                plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
                plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=90))

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
        fig.savefig(os.path.join(self.plot_outdir, "{}_corr_heatmap_Pearson.png".format(name)))
        plt.close()

    def create_job_file(self, job_name, dataset, geno_path, expr_path, cc_path,
                        leading_zeros=0, n_permutations=0,
                        permutation_index_offset=0, cpus=1, mem=4, nodes=1,
                        qos="regular"):
        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}-{}".format(job_name, dataset),
                 "#SBATCH --output={}".format(os.path.join(self.job_output_outdir, job_name + "_" + dataset + ".out")),
                 "#SBATCH --error={}".format(os.path.join(self.job_output_outdir, job_name + "_" + dataset + ".out")),
                 "#SBATCH --time={}".format(self.time),
                 "#SBATCH --cpus-per-task {}".format(cpus),
                 "#SBATCH --mem {}gb".format(mem),
                 "#SBATCH --nodes {}".format(nodes),
                 "#SBATCH --qos={}".format(qos),
                 "",
                 "module load Python/3.7.4-GCCcore-7.3.0-bare",
                 "source $HOME/env/bin/activate",
                 "",
                 "{} \\".format(self.python_executable),
                 "     -ge {} \\".format(geno_path),
                 "     -ex {} \\".format(expr_path),
                 "     -cc {} \\".format(cc_path),
                 "     -na {} \\".format(self.genotype_na),
                 "     -ac {} \\".format(self.allele_configs),
                 "     -cr {} \\".format(self.call_rate),
                 "     -hw {} \\".format(self.hw_pval),
                 "     -maf {} \\".format(self.maf),
                 "     -p {} \\".format(n_permutations),
                 "     -po {} \\".format(permutation_index_offset),
                 "     -plz {} \\".format(leading_zeros),
                 "     -od {} \\".format(self.decon_basedir),
                 "     -of {}-{}".format(self.outfolder, dataset),
                 "",
                 "deactivate",
                 ""]

        jobfile_path = os.path.join(self.job_outdir, job_name + "_" + dataset + ".sh")
        with open(jobfile_path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(jobfile_path)))

    def plot_boxplot(self, df, x="variable", y="value", hue=None, xlabel="",
                       ylabel="", name="", order=None, palette=None):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        sns.violinplot(x=x,
                       y=y,
                       hue=hue,
                       data=df,
                       order=order,
                       palette=palette,
                       cut=0,
                       dodge=True,
                       ax=ax)

        plt.setp(ax.collections, alpha=.75)

        sns.boxplot(x=x,
                    y=y,
                    hue=hue,
                    data=df,
                    order=order,
                    whis=np.inf,
                    color="white",
                    dodge=True,
                    ax=ax)

        if ax.get_legend() is not None:
            ax.get_legend().remove()

        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')

        plt.tight_layout()
        fig.savefig(os.path.join(self.plot_outdir, "{}_boxplot.png".format(name)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
