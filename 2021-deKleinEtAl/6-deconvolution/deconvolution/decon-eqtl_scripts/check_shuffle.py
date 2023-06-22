#!/usr/bin/env python3

"""
File:         check_shuffle.py
Created:      2021/07/28
Last Changed: 2021/09/07
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
import sys
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.stats import multitest

# Local application imports.


"""
Syntax:
./check_shuffle.py -id /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts -if CortexEUR-cis-WithPermutations -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis/data/SampleToDataset.txt.gz

./check_shuffle.py -id /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts -if CortexEUR-cis-HalfNormalizedMAF5 -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis/data/SampleToDataset.txt.gz

./check_shuffle.py -id /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts -if CortexEUR-cis-NormalisedMAF5 -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis/data/SampleToDataset.txt.gz
"""

# Metadata
__program__ = "Check Shuffle"
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
        indir = getattr(arguments, 'indir')
        infolder = getattr(arguments, 'infolder')
        self.geno_path = getattr(arguments, 'genotype')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.n_files = getattr(arguments, 'n_files')

        # Set variables.
        if indir is None:
            indir = str(Path(__file__).parent.parent)
        self.indir = os.path.join(indir, "decon_eqtl_with_permutation_fdr", infolder)
        self.outdir = os.path.join(self.indir, "plot")
        for dir in [self.indir, self.outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        self.dataset_cohort_dict = {
            "AMPAD-MAYO-V2": "MAYO",
            "CMC_HBCC_set2": "CMC HBCC",
            "GTEx": "GTEx",
            "AMPAD-ROSMAP-V2": "ROSMAP",
            "BrainGVEX-V2": "Brain GVEx",
            "TargetALS": "Target ALS",
            "AMPAD-MSBB-V2": "MSBB",
            "NABEC-H610": "NABEC",
            "LIBD_1M": "LIBD",
            "ENA": "ENA",
            "LIBD_h650": "LIBD",
            "GVEX": "GVEX",
            "NABEC-H550": "NABEC",
            "CMC_HBCC_set3": "CMC HBCC",
            "UCLA_ASD": "UCLA ASD",
            "CMC": "CMC",
            "CMC_HBCC_set1": "CMC HBCC"
        }

        self.cohort_palette = {
            "MAYO": "#9c9fa0",
            "CMC HBCC": "#0877b4",
            "GTEx": "#0fa67d",
            "ROSMAP": "#6950a1",
            "Brain GVEx": "#48b2e5",
            "Target ALS": "#d5c77a",
            "MSBB": "#5cc5bf",
            "NABEC": "#6d743a",
            "LIBD": "#e49d26",
            "ENA": "#d46727",
            "GVEX": "#000000",
            "UCLA ASD": "#f36d2a",
            "CMC": "#eae453",
            "NA": "#808080"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-id",
                            "--indir",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of the input path.")
        parser.add_argument("-if",
                            "--infolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the input folder.")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=False,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-n",
                            "--n_files",
                            type=int,
                            required=False,
                            default=None,
                            help="The number of files to load. Default: None.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading color info.")
        col_colors = None
        if self.std_path is not None:
            std_df = self.load_file(self.std_path, header=0, index_col=None)
            col_colors = [self.cohort_palette[self.dataset_cohort_dict[dataset]] for dataset in std_df.iloc[:, 1]]

        ########################################################################

        print("Plotting permutation order.")
        perm_order_m_list = []
        perm_order_inpaths = glob.glob(os.path.join(self.indir, "perm_orders_*"))
        perm_order_inpaths.sort(key=self.natural_keys)
        for perm_order_inpath in perm_order_inpaths:
            perm_order_m_list.append(self.load_matrix(perm_order_inpath))
        perm_order_m = np.vstack(perm_order_m_list)
        del perm_order_m_list
        print("\tShape: {}".format(perm_order_m.shape))

        self.plot_clustermap(m=perm_order_m, col_colors=col_colors,
                             filename="perm_order_m")
        del perm_order_m
        print("")

        ########################################################################

        print("Loading permutation order overlap")
        perm_order_overlap_m_list = []
        perm_order_overlap_inpaths = glob.glob(os.path.join(self.indir, "perm_order_overlap_*"))
        perm_order_overlap_inpaths.sort(key=self.natural_keys)
        for perm_order_overlap_inpath in perm_order_overlap_inpaths:
            perm_order_overlap_m_list.append(self.load_matrix(perm_order_overlap_inpath))
        perm_order_overlap_m = np.hstack(perm_order_overlap_m_list)
        del perm_order_overlap_m_list
        print("\tShape: {}".format(perm_order_overlap_m.shape))

        print("\tplotting distribution")
        self.distplot(x=perm_order_overlap_m.flatten(), xlabel="%overlap", filename="overlap")
        print("")

        ########################################################################

        print("Loading genotype matrix and calculating MAF")
        maf_a = self.load_file(os.path.join(self.indir, "MAF.txt.gz"), header=0, index_col=0).to_numpy().flatten()
        print("\tShape: {}".format(maf_a.shape))

        print("\tplotting distribution")
        self.distplot(x=maf_a, xlabel="MAF", filename="MAF")
        print("")

        ########################################################################

        print("Loading permutation p-value data")
        perm_pvalues_m_list = []
        perm_pvalues_inpaths = glob.glob(os.path.join(self.indir, "permutation_pvalues_*"))
        perm_pvalues_inpaths.sort(key=self.natural_keys)
        for perm_pvalues_inpath in perm_pvalues_inpaths:
            perm_pvalues_m_list.append(self.load_matrix(perm_pvalues_inpath))
        perm_pvalues_m = np.dstack(perm_pvalues_m_list)
        del perm_pvalues_m_list
        print("\tShape: {}".format(perm_pvalues_m.shape))

        print("\tplotting distributions")
        for cov_index in range(perm_pvalues_m.shape[1]):
            self.distplot(x=perm_pvalues_m[:, cov_index, :].flatten(), xlabel="permutation p-value", filename="perm_pvalues_cov{}".format(cov_index))
        print("")

        ########################################################################

        print("Loading real p-value data")
        decon_df = self.load_file(os.path.join(self.indir, "deconvolutionResults.txt.gz"), header=0, index_col=0)
        pvalue_columns = [x for x in decon_df.columns if x.endswith("_pvalue")]
        real_pvalues_m = decon_df.loc[:, pvalue_columns].to_numpy()

        print("\tplotting distribution")
        self.distplot(x=real_pvalues_m.flatten(), xlabel="real p-value", filename="real_pvalues")
        print("")

        ########################################################################

        print("Creating overview plots")

        # Align the shape of the matrices.
        perm_order_overlap_m = np.repeat(perm_order_overlap_m[:, np.newaxis, :], perm_pvalues_m.shape[1], axis=1)
        maf_m = np.tile(maf_a[:, np.newaxis, np.newaxis], (perm_order_overlap_m.shape[1], perm_order_overlap_m.shape[2]))
        real_pvalues_m = np.repeat(real_pvalues_m[:, :, np.newaxis], perm_pvalues_m.shape[2], axis=2)

        print("\tConstruct plot DF")
        print(real_pvalues_m.shape, perm_pvalues_m.shape, perm_order_overlap_m.shape, maf_m.shape)
        plot_df = pd.DataFrame({"real p-value": real_pvalues_m.flatten(),
                                "perm. p-value": perm_pvalues_m.flatten(),
                                "%overlap": perm_order_overlap_m.flatten(),
                                "MAF": maf_m.flatten()})
        plot_df["log10 real p-value"] = np.log10(plot_df["real p-value"] + 2.2250738585072014e-308)
        plot_df["log10 perm. p-value"] = np.log10(plot_df["perm. p-value"] + 2.2250738585072014e-308)

        print("\tPlotting single comparisons")
        self.single_regplot(df=plot_df,
                            x="log10 real p-value",
                            y="log10 perm. p-value",
                            xlabel="log10 nominal p-value",
                            ylabel="log10 perm p-value",
                            filename="nominal_vs_permuted_pvalues_log10")
        self.single_regplot(df=plot_df,
                            x="real p-value",
                            y="perm. p-value",
                            xlabel="nominal p-value",
                            ylabel="perm p-value",
                            filename="nominal_vs_permuted_pvalues")
        self.single_regplot(df=plot_df,
                            x="log10 perm. p-value",
                            y="%overlap",
                            xlabel="perm. p-value",
                            ylabel="%overlap",
                            filename="perm_pval_vs_overlap")
        self.single_regplot(df=plot_df,
                            x="log10 perm. p-value",
                            y="MAF",
                            xlabel="perm. p-value",
                            ylabel="MAF",
                            filename="perm_pval_vs_MAF")
        self.single_regplot(df=plot_df,
                            x="MAF",
                            y="%overlap",
                            xlabel="MAF",
                            ylabel="%overlap",
                            filename="MAF_vs_overlap")

        print("\tPlotting grouped comparisons")
        self.regplot(df=plot_df,
                     x="perm. p-value",
                     y="%overlap",
                     group="real p-value",
                     xlabel="perm p-value",
                     ylabel="%identical genotype",
                     filename="perm_pval_vs_overlap_panels")
        self.regplot(df=plot_df,
                     x="perm. p-value",
                     y="MAF",
                     group="real p-value",
                     xlabel="perm p-value",
                     ylabel="MAF",
                     filename="perm_pval_vs_MAF_panels")
        print("")

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
    def load_matrix(inpath):
        m = np.load(inpath)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      m.shape))
        return m

    @staticmethod
    def natural_keys(text):
        return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', text)]

    def distplot(self, x, xlabel="", title="", filename="distribution"):
        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(12, 12))

        sns.despine(fig=fig, ax=ax)

        sns.kdeplot(x, shade=True, color="#808080", ax=ax, cut=0, zorder=-1)
        ax.axvline(x.mean(), ls='--', color="#808080", zorder=-1)

        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel("density",
                      fontsize=14,
                      fontweight='bold')
        ax.set_title(title,
                     fontsize=18,
                     fontweight='bold')

        ax.annotate(
            'N = {:,}'.format(np.size(x)),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=15,
            fontweight='bold')

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def single_regplot(self, df, x, y, xlabel="", ylabel="", title="", filename="plot"):
        sns.set(rc={'figure.figsize': (12, 12)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1,
                                       ncols=2,
                                       gridspec_kw={"width_ratios": [0.85, 0.15]})
        sns.despine(fig=fig, ax=ax1)

        sns.regplot(x=x, y=y, data=df, ci=None,
                    scatter_kws={'facecolors': "#000000",
                                 'linewidth': 0,
                                 'alpha': 0.25},
                    line_kws={"color": "#b22222"},
                    ax=ax1)

        ax1.set_xlabel(xlabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_title(title,
                      fontsize=18,
                      fontweight='bold')

        # Change margins.
        xlim = ax1.get_xlim()
        ylim = ax1.get_ylim()

        xmargin = (xlim[1] - xlim[0]) * 0.05
        ymargin = (ylim[1] - ylim[0]) * 0.05

        ax1.set_xlim(xlim[0] - xmargin, xlim[1] + xmargin)
        ax1.set_ylim(ylim[0] - ymargin, ylim[1] + ymargin)

        # Set annotation.
        ax2.set_axis_off()
        coef, _ = stats.pearsonr(df[y], df[x])
        ax2.annotate(
            'r = {:.2f}'.format(coef),
            xy=(0, 0.94),
            xycoords=ax2.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold')
        ax2.annotate(
            'total N = {:,}'.format(df.shape[0]),
            xy=(0, 0.91),
            xycoords=ax2.transAxes,
            color="#404040",
            fontsize=14,
            fontweight='bold')

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def regplot(self, df, x, y, group, xlabel="", ylabel="", filename="regplot"):
        tmp_df = df.copy()
        offset = 2.2250738585072014e-308

        tmp_df.loc[:, [x, group]] = np.log10(tmp_df.loc[:, [x, group]] + offset)
        cutoff = np.log10(0.05 + offset)
        tmp_df["hue"] = "#808080"
        tmp_df.loc[tmp_df[x] < cutoff, "hue"] = "#0072B2"
        groups = np.array([-100, -75, -50, -25, -10, -5, -1, 0, None])

        ngroups = np.size(groups)
        ncols = int(np.ceil(np.sqrt(ngroups)))
        nrows = int(np.ceil(ngroups / ncols))

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharey="all",
                                 figsize=(12 * ncols, 12 * nrows))
        sns.set(color_codes=True)

        row_index = 0
        col_index = 0
        for i in range(ncols * nrows):
            if nrows == 1 and ncols == 1:
                ax = axes
            elif nrows == 1 and ncols > 1:
                ax = axes[col_index]
            elif nrows > 1 and ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            if i < ngroups:
                tmp_plot_df = tmp_df.copy()
                title = "All"
                if groups[i] is not None:
                    mask = tmp_df[group] <= groups[i]
                    title = "{} <= {}".format(groups[i], group)
                    if i > 0:
                        mask = np.logical_and(tmp_df[group] > groups[i - 1], tmp_df[group] <= groups[i])
                        title = "{} < {} <= {}".format(groups[i - 1], group, groups[i])
                    tmp_plot_df = tmp_plot_df.loc[mask, :]

                sns.despine(fig=fig, ax=ax)

                sns.regplot(x=x, y=y, data=tmp_plot_df, ci=None,
                            scatter_kws={'facecolors': tmp_plot_df["hue"],
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": "#000000"},
                            ax=ax)

                ax.axvline(cutoff, ls='--', color="#000000", zorder=-1)

                coef, _ = stats.pearsonr(tmp_plot_df[y], tmp_plot_df[x])
                ax.annotate(
                    'r = {:.2f}'.format(coef),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=15,
                    fontweight='bold')
                ax.annotate(
                    'total N = {:,}'.format(tmp_plot_df.shape[0]),
                    xy=(0.03, 0.9),
                    xycoords=ax.transAxes,
                    color="#404040",
                    fontsize=15,
                    fontweight='bold')
                ax.annotate(
                    'N = {:,}'.format(tmp_plot_df[tmp_plot_df["hue"] == "#0072B2"].shape[0]),
                    xy=(0.03, 0.86),
                    xycoords=ax.transAxes,
                    color="#0072B2",
                    fontsize=15,
                    fontweight='bold')
                ax.annotate(
                    'N = {:,}'.format(tmp_plot_df[tmp_plot_df["hue"] == "#808080"].shape[0]),
                    xy=(0.03, 0.82),
                    xycoords=ax.transAxes,
                    color="#808080",
                    fontsize=15,
                    fontweight='bold')

                ax.set_xlabel("log10 {}".format(xlabel),
                              fontsize=14,
                              fontweight='bold')
                ax.set_ylabel("{}".format(ylabel),
                              fontsize=14,
                              fontweight='bold')
                ax.set_title(title,
                             fontsize=18,
                             fontweight='bold')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def plot_clustermap(self, m, col_colors=None, filename="clustermap"):
        sys.setrecursionlimit(100000)

        sns.set(color_codes=True)
        g = sns.clustermap(m, cmap="Blues",
                           col_colors=col_colors,
                           row_cluster=True, col_cluster=False,
                           yticklabels=False, xticklabels=False,
                           figsize=(12, 9))
        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        plt.tight_layout()
        g.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > N-files: {}".format(self.n_files))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
