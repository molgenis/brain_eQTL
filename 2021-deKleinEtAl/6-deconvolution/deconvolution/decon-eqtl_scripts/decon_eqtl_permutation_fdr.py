#!/usr/bin/env python3

"""
File:         decon_eqtl_permutation_fdr.py
Created:      2021/06/07
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
from pathlib import Path
import argparse
import glob
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy import special
from scipy.optimize import minimize
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from statsmodels.stats import multitest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from venn import venn

# Local application imports.

"""
Syntax:
./decon_eqtl_permutation_fdr.py \
    -id /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts \
    -if 2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron
    
./decon_eqtl_permutation_fdr.py \
    -id /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts \
    -if 2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected 
"""

# Metadata
__program__ = "Decon-eQTL Permutation FDR"
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
        self.alpha = 0.05

        # Set variables.
        if indir is None:
            indir = str(Path(__file__).parent.parent)
        self.indir = os.path.join(indir, "decon_eqtl", infolder)
        self.outdir = os.path.join(self.indir, "plot")
        for dir in [self.indir, self.outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

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

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        ########################################################################

        # Load the nominal data.
        print("Loading nominal p-value data")
        nominal_pvalues_df = self.load_file(os.path.join(self.indir, "deconvolutionResults.txt.gz"), header=0, index_col=0)
        pvalue_columns = [x for x in nominal_pvalues_df.columns if x.endswith("_pvalue")]
        nominal_pvalues_m = nominal_pvalues_df.loc[:, pvalue_columns].to_numpy()
        rownames = nominal_pvalues_df.index.tolist()
        colnames = [x.replace("_pvalue", "") for x in pvalue_columns]
        del nominal_pvalues_df
        print("\tShape: {}".format(nominal_pvalues_m.shape))

        # Plotting.
        for i in range(nominal_pvalues_m.shape[1]):
            self.distplot(a=nominal_pvalues_m[:, i],
                          xlabel="nominal p-value",
                          title=colnames[i],
                          filename="{}_nominal_pval_distribution".format(colnames[i].lower()))
        print("")

        ########################################################################

        print("Loading permutation p-value data")
        perm_pvalues_m_list = []
        perm_pvalues_inpaths = glob.glob(os.path.join(self.indir, "permutation_pvalues_*"))
        perm_pvalues_inpaths.sort(key=self.natural_keys)
        for perm_pvalues_inpath in perm_pvalues_inpaths:
            perm_pvalues_m_list.append(self.load_matrix(perm_pvalues_inpath))
        perm_pvalues_m = np.dstack(perm_pvalues_m_list)
        n_permutations = perm_pvalues_m.shape[2]
        del perm_pvalues_m_list
        print("\tShape: {}".format(perm_pvalues_m.shape))

        # Plotting.
        for i in range(perm_pvalues_m.shape[1]):
            self.distplot(a=perm_pvalues_m[:, i, :],
                          xlabel="permuted p-value",
                          title=colnames[i],
                          filename="{}_permuted_pval_distribution".format(colnames[i].lower()))
        print("")

        ########################################################################

        print("Calculating lowest p-value per cell type per permutations")
        lowest_perm_pvalues_m = np.empty((perm_pvalues_m.shape[2], perm_pvalues_m.shape[1]), dtype=np.float64)
        for cov_index in range(perm_pvalues_m.shape[1]):
            for perm_index in range(perm_pvalues_m.shape[2]):
                lowest_perm_pvalues_m[perm_index, cov_index] = np.min(perm_pvalues_m[:, cov_index, perm_index])
        print("\tShape: {}".format(lowest_perm_pvalues_m.shape))

        # Plotting.
        for i in range(lowest_perm_pvalues_m.shape[1]):
            self.distplot(a=lowest_perm_pvalues_m[:, i],
                          xlabel="permuted p-value",
                          title=colnames[i],
                          filename="{}_lowest_permuted_pval_distribution".format(colnames[i].lower()))
        print("")

        ########################################################################

        print("Analyzing per cell type")
        bh_fdr_m = np.empty_like(nominal_pvalues_m, dtype=np.float64)
        emp_fdr_m = np.empty_like(nominal_pvalues_m, dtype=np.float64)
        qvalues_m = np.empty_like(nominal_pvalues_m, dtype=np.float64)
        per_eqtl_qvalues_m = np.empty_like(nominal_pvalues_m, dtype=np.float64)
        for cov_index in range(nominal_pvalues_m.shape[1]):
            stats_m = np.empty((nominal_pvalues_m.shape[0], 9), dtype=np.float64)

            column = colnames[cov_index]
            print("  Analyzing: {}".format(column))

            # Extract the data.
            nominal_pvalues = nominal_pvalues_m[:, cov_index]
            perm_pvalues = perm_pvalues_m[:, cov_index].flatten()
            lowest_perm_pvalues = lowest_perm_pvalues_m[:, cov_index]
            print("\tnominal p-values: {}".format(nominal_pvalues.shape))
            print("\tPermutation p-values: {}".format(perm_pvalues.shape))
            print("\tLowest permutation p-values: {}".format(lowest_perm_pvalues.shape))
            print("")

            print("\tMethod 1: Benjamini-Hochberg FDR")
            bh_fdr = multitest.multipletests(nominal_pvalues, method='fdr_bh')[1]

            print("\tMethod 2: EMP style FDR")
            ranks = np.array([np.sum(nominal_pvalues <= nominal_pvalue) for nominal_pvalue in nominal_pvalues])
            perm_ranks = np.array([np.sum(perm_pvalues <= nominal_pvalue) for nominal_pvalue in nominal_pvalues])
            emp_fdr = (perm_ranks / n_permutations) / ranks
            emp_fdr[emp_fdr > 1] = 1

            print("\tMethod 3: Fast-QTL style FDR")
            print("\t  Fitting beta function.")
            (a, b), nfev, nit = self.fit_beta_distribution(p=lowest_perm_pvalues)
            print(a, b, nfev, nit)

            # Plotting distributions.
            self.histplot(nominal=nominal_pvalues,
                          permuted=perm_pvalues,
                          lowest_permuted=lowest_perm_pvalues,
                          beta_parameters=(a, b),
                          nit=nit,
                          nfev=nfev,
                          filename="{}_distribution".format(column))

            # Calculating adjusted p-values.
            print("\tCalculating adjusted p-values.")
            adj_pvalues = stats.beta.cdf(nominal_pvalues, a, b)

            print("\tCalculating q-values.")
            qvalues = self.qvalues(p=adj_pvalues)

            print("\tMethod 4: per eQTL FDR")
            print("\t  Calculating adjusted p-values.")
            per_eqtl_ranks = np.empty(nominal_pvalues_m.shape[0], dtype=np.float64)
            per_eqtl_adj_pvalues = np.empty(nominal_pvalues_m.shape[0], dtype=np.float64)
            for eqtl_index in range(nominal_pvalues_m.shape[0]):
                nominal_pvalue = nominal_pvalues_m[eqtl_index, cov_index]
                perm_pvalues = perm_pvalues_m[eqtl_index, cov_index].flatten()
                rank = np.sum(perm_pvalues <= nominal_pvalue)
                adj_pvalue = (rank + 0.5) / (n_permutations + 1)
                # if adj_pvalue < nominal_pvalue:
                #     adj_pvalue = nominal_pvalue

                per_eqtl_ranks[eqtl_index] = rank
                per_eqtl_adj_pvalues[eqtl_index] = adj_pvalue

            print("\t  Calculating q-values.")
            per_eqtl_qvalues = self.qvalues(p=per_eqtl_adj_pvalues)

            if column.startswith("Excitatory"):
                print("\tPlotting.")
                for log10 in [(False, False), (True, True)]:
                    outdir = self.outdir
                    if log10[0] or log10[1]:
                        outdir = os.path.join(self.outdir, "log10")
                    if not os.path.exists(outdir):
                        os.makedirs(outdir)

                    self.regplot(x=nominal_pvalues, y=bh_fdr,
                                 xlabel="p-values",
                                 ylabel="BH fdr-values",
                                 log10=log10,
                                 filename="{}_nominal_pval_vs_bh_fdr".format(column),
                                 outdir=outdir)

                    self.regplot(x=nominal_pvalues, y=ranks,
                                 xlabel="nominal p-values",
                                 ylabel="ranks",
                                 log10=(log10[0], False),
                                 cutoff=(True, False),
                                 filename="{}_nominal_pval_vs_ranks".format(column),
                                 outdir=outdir)
                    self.regplot(x=nominal_pvalues, y=perm_ranks,
                                 xlabel="nominal p-values",
                                 ylabel="perm_ranks",
                                 log10=(log10[0], False),
                                 cutoff=(True, False),
                                 filename="{}_nominal_pval_vs_perm_ranks".format(column),
                                 outdir=outdir)
                    self.regplot(x=nominal_pvalues, y=emp_fdr,
                                 xlabel="p-values",
                                 ylabel="EMP fdr-values",
                                 log10=log10,
                                 filename="{}_nominal_pval_vs_emp_fdr".format(column),
                                 outdir=outdir)

                    self.regplot(x=nominal_pvalues, y=adj_pvalues,
                                 xlabel="nominal p-values",
                                 ylabel="adj. p-values",
                                 log10=log10,
                                 filename="{}_nominal_pval_vs_adj_pval".format(column),
                                 outdir=outdir)
                    self.regplot(x=adj_pvalues, y=qvalues,
                                 xlabel="adj. p-values",
                                 ylabel="q-values",
                                 log10=log10,
                                 filename="{}_adj_pval_vs_qval".format(column),
                                 outdir=outdir)
                    self.regplot(x=adj_pvalues, y=perm_ranks,
                                 xlabel="adj. p-values",
                                 ylabel="perm_ranks",
                                 log10=log10,
                                 filename="{}_adj_pval_vs_perm_ranks".format(column),
                                 outdir=outdir)
                    self.regplot(x=nominal_pvalues, y=qvalues,
                                 xlabel="p-values",
                                 ylabel="q-values",
                                 log10=log10,
                                 filename="{}_nominal_pval_vs_qval".format(column),
                                 outdir=outdir)

                    self.regplot(x=nominal_pvalues, y=per_eqtl_adj_pvalues,
                                 xlabel="nominal p-values",
                                 ylabel="adj. p-values",
                                 log10=log10,
                                 filename="{}_per_eqtl_nominal_pval_vs_adj_pval".format(column),
                                 outdir=outdir)
                    self.regplot(x=per_eqtl_adj_pvalues, y=per_eqtl_qvalues,
                                 xlabel="adj. p-values",
                                 ylabel="q-values",
                                 log10=log10,
                                 filename="{}_per_eqtl_adj_pval_vs_qval".format(column),
                                 outdir=outdir)
                    self.regplot(x=nominal_pvalues, y=per_eqtl_qvalues,
                                 xlabel="p-values",
                                 ylabel="q-values",
                                 log10=log10,
                                 filename="{}_per_eqtl_nominal_pval_vs_qval".format(column),
                                 outdir=outdir)

                    self.regplot(x=qvalues, y=bh_fdr,
                                 xlabel="q-values",
                                 ylabel="BH fdr-values",
                                 log10=log10,
                                 filename="{}_qval_vs_bh_fdr".format(column),
                                 outdir=outdir)
                    self.regplot(x=qvalues, y=emp_fdr,
                                 xlabel="q-values",
                                 ylabel="EMP fdr-values",
                                 log10=log10,
                                 filename="{}_qval_vs_emp_fdr".format(column),
                                 outdir=outdir)
                    self.regplot(x=qvalues, y=per_eqtl_qvalues,
                                 xlabel="q-values",
                                 ylabel="per eQTL qvalues",
                                 log10=log10,
                                 filename="{}_qval_vs_per_eqtl_qval".format(column),
                                 outdir=outdir)
                    self.regplot(x=bh_fdr, y=emp_fdr,
                                 xlabel="BH fdr-values",
                                 ylabel="EMP fdr-values",
                                 log10=log10,
                                 filename="{}_bh_fdr_vs_emp_fdr".format(column),
                                 outdir=outdir)
                    self.regplot(x=bh_fdr, y=per_eqtl_qvalues,
                                 xlabel="BH fdr-values",
                                 ylabel="per eQTL qvalues",
                                 log10=log10,
                                 filename="{}_bh_fdr_vs_per_eqtl_qval".format(column),
                                 outdir=outdir)
                    self.regplot(x=emp_fdr, y=per_eqtl_qvalues,
                                 xlabel="EMP fdr-values",
                                 ylabel="per eQTL qvalues",
                                 log10=log10,
                                 filename="{}_emp_fdr_vs_per_eqtl_qval".format(column),
                                 outdir=outdir)

            # Saving results.
            bh_fdr_m[:, cov_index] = bh_fdr
            emp_fdr_m[:, cov_index] = emp_fdr
            qvalues_m[:, cov_index] = qvalues
            per_eqtl_qvalues_m[:, cov_index] = per_eqtl_qvalues

            # Saving stats.
            stats_m[:, 0] = nominal_pvalues
            stats_m[:, 1] = ranks
            stats_m[:, 2] = perm_ranks
            stats_m[:, 3] = emp_fdr
            stats_m[:, 4] = adj_pvalues
            stats_m[:, 5] = qvalues
            stats_m[:, 6] = per_eqtl_ranks
            stats_m[:, 7] = per_eqtl_adj_pvalues
            stats_m[:, 8] = per_eqtl_qvalues

            # Save thats data frame.
            stats_df = pd.DataFrame(stats_m, index=rownames, columns=["nominal p-value", "rank", "permutation rank", "EMP-FDR", "adj. p-value", "q-value", "permutation rank (per eQTL)", "adj. p-value (per eQTL)", "q-value (per eQTL)"])
            self.save_file(df=stats_df, outpath=os.path.join(self.indir, "{}_stats_df.txt.gz".format(column)))
            print("")
            del stats_m, stats_df
        print("")

        ########################################################################

        print("Creating VENN diagram")
        venn_dict = {}
        for m, name in zip([bh_fdr_m, emp_fdr_m, qvalues_m, per_eqtl_qvalues_m], ["BH_FDR", "EMP_FDR", "qvalues", "per_eQTL"]):
            signif_hits = set()
            for j, colname in enumerate(colnames):
                signif_hits.update(set(["{}_{}".format(rownames[i], colname) for i in range(len(rownames)) if m[i, j] <= 0.05]))
            venn_dict["{} [N = {:,}]".format(name, len(signif_hits))] = signif_hits

        self.vennplot(data=venn_dict,
                      title="interaction overlap",
                      filename="signif_interaction_overlap")

        print("Saving data frames")
        for m, name in zip([bh_fdr_m, emp_fdr_m, qvalues_m, per_eqtl_qvalues_m], ["BH_FDR", "EMP_FDR", "qvalues", "per_eQTL"]):
            self.print_n_signif(m=m, colnames=colnames, type=name)
            df = pd.DataFrame(m, columns=colnames, index=rownames)
            self.save_file(df=df, outpath=os.path.join(self.indir, "{}.txt.gz".format(name)))

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

    def distplot(self, a, xlabel="", title="", filename="distribution"):
        data = a.flatten()
        data = data[~np.isnan(data)]


        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(12, 12))

        sns.despine(fig=fig, ax=ax)

        sns.kdeplot(data, shade=True, color="#808080", ax=ax, cut=0, zorder=-1)
        ax.axvline(np.mean(data), ls='--', color="#808080", zorder=-1)

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
            'N = {:,}'.format(np.size(data)),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=15,
            fontweight='bold')

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\t  Saved figure: {} ".format(os.path.basename(outpath)))

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

    def fit_beta_distribution(self, p, a_bnd=(0.1, 10), b_bnd=(1, 1000000)):
        a, b = self.beta_distribution_initial_guess(x=p)
        x0 = np.array([min(max(a, a_bnd[0]), a_bnd[1]), min(max(b, b_bnd[0]), b_bnd[1])])
        res = minimize(self.beta_distribution_mle_function,
                       x0=x0,
                       args=(p, ),
                       method='nelder-mead',
                       bounds=(a_bnd, b_bnd),
                       options={"maxiter": 1000, "disp": True})
        return res.x, res.nfev, res.nit

    @staticmethod
    def beta_distribution_initial_guess(x):
        """
        https://stats.stackexchange.com/questions/13245/which-is-a-good-tool-to-compute-parameters-for-a-beta-distribution
        """
        mean = np.mean(x)
        var = np.var(x)
        a = mean * ((mean * (1 - mean) / var) - 1)
        b = (1 - mean) * ((mean * (1 - mean) / var) - 1)
        return a, b

    @staticmethod
    def beta_distribution_mle_function(x, p):
        k, n = x
        ll = (k - 1) * np.sum(np.log(p)) + (n - 1) * np.sum(np.log(1 - p)) - np.size(p) * special.betaln(k, n)
        return -1 * ll

    @staticmethod
    def qvalues(p):
        qvalue = importr("qvalue")
        pvals = robjects.FloatVector(p)
        qobj = robjects.r['qvalue'](pvals)
        return np.array(qobj.rx2('qvalues'))

    def print_n_signif(self, m, colnames, type):
        print("\nN-interaction ({} <= {}):".format(type, self.alpha))
        n_hits_a = (m <= self.alpha).sum(axis=0)
        n_hits_total = np.sum(n_hits_a)
        cov_length = np.max([len(x) for x in colnames])
        hits_length = np.max([len(str(x)) for x in n_hits_a] + [len(str(n_hits_total))])
        for n_hits, cell_type in zip(n_hits_a, colnames):
            print("\t{:{}s}  {:{}d}".format(cell_type, cov_length, n_hits, hits_length))
        print("\t{}".format("".join(["-"] * cov_length)))
        print("\t{:{}s}  {:{}d}".format("total", cov_length, n_hits_total, hits_length))

        print("", flush=True)

    def histplot(self, nominal, permuted, lowest_permuted, beta_parameters, nit, nfev, filename="plot"):
        sns.set(rc={'figure.figsize': (27, 12)})
        sns.set_style("ticks")
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3)

        for ax, data, log10_transform, add_beta_func, title in zip((ax1, ax2, ax3),
                                                                   (nominal, permuted, lowest_permuted),
                                                                   (False, False, True),
                                                                   (False, False, True),
                                                                   ("nominal", "permuted", "lowest + beta fit")):
            sns.despine(fig=fig, ax=ax1)

            plot_data = np.copy(data)
            xlabel = "p-value"
            if log10_transform:
                plot_data = -1 * np.log10(plot_data + 2.2250738585072014e-308)
                xlabel = "-log10 p-value"

            # plot in sections.
            sns.histplot(data=plot_data, ax=ax, color="black")

            ax.annotate(
                'N = {:,}'.format(np.size(plot_data)),
                xy=(0.7, 0.94),
                xycoords=ax.transAxes,
                color="#000000",
                fontsize=12)

            if add_beta_func:
                x_lim = ax.get_xlim()
                y_lim = ax.get_ylim()

                # Fit.
                a, b = beta_parameters
                x = np.linspace(np.min(data), np.max(data), 1000)
                y = stats.beta.pdf(x, a, b)

                # Rescale.
                (min_val, max_val) = y_lim
                min_y = np.min(y)
                max_y = np.max(y)
                y = (max_val - min_val) * ((y - min_y) / (max_y - min_y)) + min_val

                # Log transform.
                if log10_transform:
                    x = -1 * np.log10(x + 2.2250738585072014e-308)
                    y = -1 * np.log10(y + 2.2250738585072014e-308)

                ax.plot(x, y, 'r--', label='beta pdf')
                ax.set_xlim(x_lim)
                ax.set_ylim(y_lim)

                ax.annotate(
                    'a = {:.2e}'.format(a),
                    xy=(0.7, 0.90),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=12)
                ax.annotate(
                    'b = {:.2e}'.format(b),
                    xy=(0.7, 0.86),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=12)
                ax.annotate(
                    'Iterations: {}'.format(nit),
                    xy=(0.7, 0.82),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=12)
                ax.annotate(
                    'Function evaluations: {}'.format(nfev),
                    xy=(0.7, 0.78),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=12)

            ax.set_title(title,
                         fontsize=18,
                         fontweight='bold')
            ax.set_ylabel("count",
                          fontsize=14,
                          fontweight='bold')
            ax.set_xlabel(xlabel,
                          fontsize=14,
                          fontweight='bold')

        fig.savefig(os.path.join(self.outdir, "{}.png".format(filename)))
        plt.close()

    def regplot(self, x, y, xlabel="", ylabel="", title="", filename="plot", log10=(True, True), cutoff=(True, True), outdir=None):
        if outdir is None:
            outdir = self.outdir

        offset = 2.2250738585072014e-308

        sns.set(rc={'figure.figsize': (12, 12)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1,
                                       ncols=2,
                                       gridspec_kw={"width_ratios": [0.85, 0.15]})
        sns.despine(fig=fig, ax=ax1)

        x_cutoff = -np.Inf
        y_cutoff = -np.Inf
        if cutoff[0]:
            x_cutoff = 0.05
        if cutoff[1]:
            y_cutoff = 0.05

        if log10[0]:
            x = np.log10(x + offset)
            x_cutoff = np.log10(x_cutoff + offset)
        if log10[1]:
            y = np.log10(y + offset)
            y_cutoff = np.log10(y_cutoff + offset)

        plot_df = pd.DataFrame({"x": x, "y": y})
        plot_df["hue"] = "#808080"
        plot_df.loc[(plot_df["x"] < x_cutoff) & (
                    plot_df["y"] >= y_cutoff), "hue"] = "#0072B2"
        plot_df.loc[(plot_df["x"] >= x_cutoff) & (
                    plot_df["y"] < y_cutoff), "hue"] = "#D55E00"
        plot_df.loc[(plot_df["x"] < x_cutoff) & (
                    plot_df["y"] < y_cutoff), "hue"] = "#009E73"

        sns.regplot(x="x", y="y", data=plot_df, ci=None,
                    scatter_kws={'facecolors': plot_df["hue"],
                                 'linewidth': 0,
                                 'alpha': 0.75},
                    line_kws={"color": "#000000"},
                    ax=ax1)

        if y_cutoff != -np.Inf:
            ax1.axhline(y_cutoff, ls='--', color="#000000", zorder=-1)
        if x_cutoff != -np.Inf:
            ax1.axvline(x_cutoff, ls='--', color="#000000", zorder=-1)

        ax1.set_xlabel("{}{}".format("log10 " if log10[0] else "", xlabel),
                       fontsize=14,
                       fontweight='bold')
        ax1.set_ylabel("{}{}".format("log10 " if log10[1] else "", ylabel),
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

        new_xlim = (xlim[0] - xmargin, xlim[1] + xmargin)
        new_ylim = (ylim[0] - ymargin, ylim[1] + ymargin)

        ax1.set_xlim(new_xlim[0], new_xlim[1])
        ax1.set_ylim(new_ylim[0], new_ylim[1])

        # Set diagonal.
        min_pos = min(new_xlim[0], new_ylim[0])
        max_pos = max(new_xlim[1], new_ylim[1])
        ax1.plot([min_pos, max_pos],
                 [min_pos, max_pos],
                 ls='--', color="#000000", zorder=-1)

        # Set annotation.
        ax2.set_axis_off()
        coef, _ = stats.pearsonr(plot_df["y"], plot_df["x"])
        ax2.annotate(
            'r = {:.2f}'.format(coef),
            xy=(0, 0.94),
            xycoords=ax2.transAxes,
            color="#000000",
            fontsize=14,
            fontweight='bold')
        ax2.annotate(
            'total N = {:,}'.format(plot_df.shape[0]),
            xy=(0, 0.91),
            xycoords=ax2.transAxes,
            color="#404040",
            fontsize=14,
            fontweight='bold')
        ax2.annotate(
            'N = {:,}'.format(plot_df[plot_df["hue"] == "#009E73"].shape[0]),
            xy=(0, 0.88),
            xycoords=ax2.transAxes,
            color="#009E73",
            fontsize=14,
            fontweight='bold')
        ax2.annotate(
            'N = {:,}'.format(plot_df[plot_df["hue"] == "#0072B2"].shape[0]),
            xy=(0, 0.85),
            xycoords=ax2.transAxes,
            color="#0072B2",
            fontsize=14,
            fontweight='bold')
        ax2.annotate(
            'N = {:,}'.format(plot_df[plot_df["hue"] == "#D55E00"].shape[0]),
            xy=(0, 0.82),
            xycoords=ax2.transAxes,
            color="#D55E00",
            fontsize=14,
            fontweight='bold')
        ax2.annotate(
            'N = {:,}'.format(plot_df[plot_df["hue"] == "#808080"].shape[0]),
            xy=(0, 0.79),
            xycoords=ax2.transAxes,
            color="#808080",
            fontsize=14,
            fontweight='bold')

        file_appendix = ""
        if log10[0] or log10[1]:
            file_appendix += "_log10"
            if log10[0]:
                file_appendix += "_xTrue"
            if log10[1]:
                file_appendix += "_yTrue"

        outpath = os.path.join(outdir, "{}{}.png".format(filename, file_appendix))
        fig.savefig(outpath)
        plt.close()
        print("\tsaved plot: {}".format(outpath))

    def vennplot(self, data, title, filename):
        sns.set(rc={'figure.figsize': (12, 12)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        venn(data, ax=ax)

        ax.set_title(title,
                     fontsize=35,
                     fontweight='bold')

        fig.savefig(os.path.join(self.outdir, "{}.png".format(filename)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
