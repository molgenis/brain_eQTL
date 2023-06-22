#!/usr/bin/env python3

"""
File:         plot_cell_fraction_for_snp.py
Created:      2022/05/02
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
__program__ = "Plot Cell Fraction for SNP"
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
./plot_cell_fraction_for_snp.py \
    -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/2022-03-28-Cortex-cis-NegativeToZero-DatasetAndRAMCorrected-DeconvolutionTable-Transposed.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/data/MetaBrain_STD_cortex.txt.gz \
    -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -s 7:12244161:rs1990622:A_G \
    -e png pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.geno_path = getattr(arguments, 'genotype')
        self.snps = getattr(arguments, 'snps')
        self.extensions = getattr(arguments, 'extensions')
        self.genotype_na = -1
        self.call_rate = 0.95
        self.hw_pval = 1e-4
        self.maf = 0.01

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "cell_fraction_per_snp")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.cohort_info = [
            # ("GTE-AFR-AMPAD-MSBB-V2", "AMPAD-MSBB (AFR; n=TMP)"),
            ("GTE-EUR-AMPAD-MAYO-V2", "AMPAD-MAYO (EUR; n=TMP)"),
            ("GTE-EUR-AMPAD-MSBB-V2", "AMPAD-MSBB (EUR; n=TMP)"),
            ("GTE-EUR-AMPAD-ROSMAP-V2", "AMPAD-ROSMAP (EUR; n=TMP)"),
            ("GTE-EUR-BrainGVEX-V2", "BrainGVEX (EUR; n=TMP)"),
            ("GTE-AFR-CMC", "CMC (AFR; n=TMP)"),
            ("GTE-EUR-CMC", "CMC (EUR; n=TMP)"),
            # ("GTE-AFR-CMC_HBCC_set1", "CMC_HBCC_set1 (AFR; n=TMP)"),
            # ("GTE-AFR-CMC_HBCC_set2", "CMC_HBCC_set2 (AFR; n=TMP)"),
            ("GTE-EUR-CMC_HBCC_set2", "CMC_HBCC_set2 (EUR; n=TMP)"),
            ("GTE-EUR-CMC_HBCC_set3", "CMC_HBCC_set3 (EUR; n=TMP)"),
            ("GTE-EUR-GTEx", "GTEx (EUR; n=TMP)"),
            ("GTE-EUR-GVEX", "GVEX (EUR; n=TMP)"),
            ("GTE-AFR-LIBD_1M", "LIBD_1M (AFR; n=TMP)"),
            ("GTE-EUR-LIBD_1M", "LIBD_1M (EUR; n=TMP)"),
            ("GTE-AFR-LIBD_h650", "LIBD_h650 (AFR; n=TMP)"),
            ("GTE-EUR-LIBD_h650", "LIBD_h650 (EUR; n=TMP)"),
            ("GTE-EUR-NABEC-H550", "NABEC-H550 (EUR; n=TMP)"),
            ("GTE-EUR-NABEC-H610", "NABEC-H610 (EUR; n=TMP)"),
            ("GTE-EUR-UCLA_ASD", "UCLA_ASD (EUR; n=TMP)"),
            ("GTE-EUR-TargetALS", "TargetALS (EUR; n=TMP)")
        ]
        self.cohort_dict = {a: b for a, b in self.cohort_info}

        # Color map.
        self.colormap = {
            "Astrocyte": "#D55E00",
            "EndothelialCell": "#CC79A7",
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Microglia": "#E69F00",
            "Oligodendrocyte": "#009E73",
            "OtherNeuron": "#2690ce"
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fracttions matrix")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=False,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix.")
        parser.add_argument("-s",
                            "--snps",
                            nargs="*",
                            type=str,
                            required=True,
                            help="The snps to plot.")
        parser.add_argument("-e",
                            "--extensions",
                            type=str,
                            nargs="+",
                            default=["png"],
                            choices=["eps", "pdf", "pgf", "png", "ps", "raw", "rgba", "svg", "svgz"],
                            help="The output file format(s), default: ['png']")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data.")
        cf_df = self.load_file(self.cf_path, header=0, index_col=0).T
        std_df = self.load_file(self.std_path, header=0, index_col=None)
        geno_df = self.load_file(self.geno_path, header=0, index_col=0)
        geno_df = geno_df.groupby(geno_df.index).first()

        print("### Step2 ###")
        print("Pre-processing data frame")
        std_df.index = std_df["rnaseq_id"]
        df = cf_df.merge(std_df, left_index=True, right_index=True)
        df["cohort_id"] = "GTE-" + df["ethnicity"] + "-" + df["dataset"]
        print(df["cohort_id"].value_counts())
        print(df)

        print("### Step3 ###")
        print("Checking genotype statistics")
        geno_df = self.check_genotype_matrix(geno_df=geno_df,
                                             std_df=std_df)

        print("### Step4 ###")
        print("Analyzing data.")
        for snp in self.snps:
            print("\tSNP: {}".format(snp))
            outdir = os.path.join(self.outdir, snp)
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            samples = geno_df.columns[geno_df.loc[snp, :] != self.genotype_na]
            snp_df = df.loc[samples, :].copy()
            counts = snp_df["cohort_id"].value_counts()
            cohort_label_trans_dict = {}
            for _, row in snp_df.iterrows():
                cohort_label_trans_dict[row["cohort_id"]] = self.cohort_dict[row["cohort_id"]].replace("TMP", str(counts[row["cohort_id"]]))
            snp_df["cohort_label"] = snp_df["cohort_id"].map(cohort_label_trans_dict)

            order = []
            for cohort_id, _ in self.cohort_info:
                if cohort_id in cohort_label_trans_dict:
                    order.append(cohort_label_trans_dict[cohort_id])

            snp_df_meta = snp_df.copy()
            snp_df_meta["cohort_id"] = "Meta-Analysis"
            snp_df_meta["cohort_label"] = "Meta-Analysis (n={:,})".format(snp_df.shape[0])
            order.append("Meta-Analysis (n={:,})".format(snp_df.shape[0]))
            snp_df = pd.concat([snp_df, snp_df_meta], axis=0)
            print(snp_df)

            print("\tPlotting")
            self.create_violin_plot(
                df=snp_df,
                y="cohort_label",
                order=order,
                columns=cf_df.columns,
                palette=self.colormap,
                title=snp,
                outdir=outdir
            )

            print("\tCalculating means")
            means = {}
            for celltype in cf_df.columns:
                data = []
                for cohort in order:
                    cohort_data = snp_df.loc[snp_df["cohort_label"] == cohort, celltype]
                    data.append(cohort_data.mean())
                s = pd.Series(data, index=order)
                means[celltype] = s

            print("\tCalculating t-tests")
            ttest_results = {}
            foldc_results = {}
            n_tests = 0
            for celltype in cf_df.columns:
                t_test_data = []
                fold_change_data = []
                for i, cohort1 in enumerate(order):
                    cohort1_data = snp_df.loc[snp_df["cohort_label"] == cohort1, celltype].copy()
                    t_test_row = []
                    fold_change_row = []
                    for j, cohort2 in enumerate(order):
                        if i <= j:
                            t_test_row.append(np.nan)
                            fold_change_row.append(np.nan)
                            continue

                        cohort2_data = snp_df.loc[snp_df["cohort_label"] == cohort2, celltype].copy()
                        t2, p2 = stats.ttest_ind(cohort1_data, cohort2_data)
                        n_tests += 1

                        t_test_row.append(p2)
                        fold_change_row.append(means[celltype][cohort1] / means[celltype][cohort2])
                    t_test_data.append(t_test_row)
                    fold_change_data.append(fold_change_row)

                ct_ttest_df = pd.DataFrame(t_test_data,
                                           index=order,
                                           columns=order)
                ct_foldc_df = pd.DataFrame(fold_change_data,
                                           index=order,
                                           columns=order)

                appendix_df = pd.DataFrame([None, n_tests, 0.05 / n_tests],
                                           index=["", "Nr Tests", "Bonferroni threshold"],
                                           columns=[ct_ttest_df.columns[0]])
                ct_ttest_df = pd.concat([ct_ttest_df, appendix_df], axis=0)

                ttest_results[celltype] = ct_ttest_df
                foldc_results[celltype] = ct_foldc_df

            print("\tSaving results")
            with pd.ExcelWriter(os.path.join(outdir, "ttest_matrix.xlsx")) as writer:
                for ct, ttest_df in ttest_results.items():
                    ttest_df.to_excel(writer,
                                      sheet_name=ct,
                                      na_rep="",
                                      index=True)

            with pd.ExcelWriter(os.path.join(outdir, "fold_change_matrix.xlsx")) as writer:
                for ct, foldc_df in foldc_results.items():
                    foldc_df.to_excel(writer,
                                      sheet_name=ct,
                                      na_rep="",
                                      index=True)

    @staticmethod
    def load_file(path, header, index_col, sep="\t", nrows=None, low_memory=True):
        if path.endswith(".pkl"):
            df = pd.read_pickle(path)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def check_genotype_matrix(self, geno_df, std_df):
        for dataset in std_df.iloc[:, 1].unique():
            print("  {}".format(dataset))

            # Get the samples of this dataset.
            samples = [sample for sample in std_df.loc[std_df.iloc[:, 1] == dataset, :].iloc[:, 0].tolist() if sample in geno_df.columns]
            if len(samples) == 0:
                continue
            geno_dataset_df = geno_df.loc[:, samples].copy()

            # Check the present samples.
            n_a = np.sum(geno_dataset_df != self.genotype_na, axis=1)
            n_mask = np.isnan(n_a)

            # Check the call rate.
            call_rate_s = (geno_dataset_df != self.genotype_na).astype(int).sum(axis=1) / len(samples)
            cr_mask = (call_rate_s < self.call_rate).to_numpy(dtype=bool)

            # Check the MAF.
            rounded_m = geno_dataset_df.copy().to_numpy(dtype=np.float64)
            rounded_m = np.rint(rounded_m)
            zero_a = np.sum(rounded_m == 0, axis=1)
            one_a = np.sum(rounded_m == 1, axis=1)
            two_a = np.sum(rounded_m == 2, axis=1)
            allele1_a = (zero_a * 2) + one_a
            allele2_a = (two_a * 2) + one_a

            maf_a = np.empty_like(n_a, dtype=np.float64)
            maf_a[:] = np.nan
            maf_a[n_a > 0] = np.minimum(allele1_a[n_a > 0], allele2_a[n_a > 0]) / (allele1_a[n_a > 0] + allele2_a[n_a > 0])
            maf_mask = maf_a <= self.maf

            # Check the HWE p-value.
            hwe_pvalues_a = np.empty_like(n_a, dtype=np.float64)
            hwe_pvalues_a[:] = np.nan
            hwe_pvalues_a[n_a > 0] = self.calc_hwe_pvalue(obs_hets=one_a[n_a > 0],
                                                          obs_hom1=zero_a[n_a > 0],
                                                          obs_hom2=two_a[n_a > 0])
            hwe_mask = hwe_pvalues_a < self.hw_pval

            output_df = pd.DataFrame({"N": n_a,
                                      "0": zero_a,
                                      "1": one_a,
                                      "2": two_a,
                                      "HW pval": hwe_pvalues_a,
                                      "allele1": allele1_a,
                                      "allele2": allele2_a,
                                      "MAF": maf_a
                                      }, index=geno_df.index)
            print(output_df.loc[self.snps, :])

            # Set samples that did not meet requirements to NaN.
            mask = n_mask | cr_mask | maf_mask | hwe_mask
            if np.sum(mask) > 0:
                print("\t  {:,} eQTL(s) failed sample size threshold".format(np.sum(n_mask)))
                print("\t  {:,} eQTL(s) failed the call rate threshold".format(np.sum(cr_mask)))
                print("\t  {:,} eQTL(s) failed the MAF threshold".format(np.sum(maf_mask)))
                print("\t  {:,} eQTL(s) failed the Hardy-Weinberg p-value threshold".format(np.sum(hwe_mask)))
                print("\t  ----------------------------------------")
                print("\t  {:,} eQTL(s) are discarded in total".format(np.sum(mask)))
            geno_df.loc[mask, samples] = self.genotype_na

        return geno_df

    @staticmethod
    def calc_hwe_pvalue(obs_hets, obs_hom1, obs_hom2):
        """
        exact SNP test of Hardy-Weinberg Equilibrium as described in Wigginton,
        JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
        Hardy-Weinberg Equilibrium. AJHG 76: 887-893

        Adapted by M.Vochteloo to work on matrices.
        """
        if not 'int' in str(obs_hets.dtype) or not 'int' in str(obs_hets.dtype) or not 'int' in str(obs_hets.dtype):
            obs_hets = np.rint(obs_hets)
            obs_hom1 = np.rint(obs_hom1)
            obs_hom2 = np.rint(obs_hom2)

        # Force homc to be the max and homr to be the min observed genotype.
        obs_homc = np.maximum(obs_hom1, obs_hom2)
        obs_homr = np.minimum(obs_hom1, obs_hom2)

        # Calculate some other stats we need.
        rare_copies = 2 * obs_homr + obs_hets
        l_genotypes = obs_hets + obs_homc + obs_homr
        n = np.size(obs_hets)

        # Get the distribution midpoint.
        mid = np.rint(rare_copies * (2 * l_genotypes - rare_copies) / (2 * l_genotypes)).astype(np.int)
        mid[mid % 2 != rare_copies % 2] += 1

        # Calculate the start points for the evaluation.
        curr_homr = (rare_copies - mid) / 2
        curr_homc = l_genotypes - mid - curr_homr

        # Calculate the left side.
        left_steps = np.floor(mid / 2).astype(int)
        max_left_steps = np.max(left_steps)
        left_het_probs = np.zeros((n, max_left_steps + 1), dtype=np.float64)
        left_het_probs[:, 0] = 1
        for i in np.arange(0, max_left_steps, 1, dtype=np.float64):
            prob = left_het_probs[:, int(i)] * (mid - (i * 2)) * ((mid - (i * 2)) - 1.0) / (4.0 * (curr_homr + i + 1.0) * (curr_homc + i + 1.0))
            prob[mid - (i * 2) <= 0] = 0
            left_het_probs[:, int(i) + 1] = prob

        # Calculate the right side.
        right_steps = np.floor((rare_copies - mid) / 2).astype(int)
        max_right_steps = np.max(right_steps)
        right_het_probs = np.zeros((n, max_right_steps + 1), dtype=np.float64)
        right_het_probs[:, 0] = 1
        for i in np.arange(0, max_right_steps, 1, dtype=np.float64):
            prob = right_het_probs[:, int(i)] * 4.0 * (curr_homr - i) * (curr_homc - i) / (((i * 2) + mid + 2.0) * ((i * 2) + mid + 1.0))
            prob[(i * 2) + mid >= rare_copies] = 0
            right_het_probs[:, int(i) + 1] = prob

        # Combine the sides.
        het_probs = np.hstack((np.flip(left_het_probs, axis=1), right_het_probs[:, 1:]))

        # Normalize.
        sum = np.sum(het_probs, axis=1)
        het_probs = het_probs / sum[:, np.newaxis]

        # Replace values higher then probability of obs_hets with 0.
        threshold_col_a = (max_left_steps - left_steps) + np.floor(obs_hets / 2).astype(int)
        threshold = np.array([het_probs[i, threshold_col] for i, threshold_col in enumerate(threshold_col_a)])
        het_probs[het_probs > threshold[:, np.newaxis]] = 0

        # Calculate the p-values.
        p_hwe = np.sum(het_probs, axis=1)
        p_hwe[p_hwe > 1] = 1

        return p_hwe

    def create_violin_plot(self, df, y="y", columns=None, order=None, title="",
                           palette=None, outdir=None):
        ncols = 1
        if columns is not None:
            ncols = len(columns)

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=ncols,
                                 nrows=1,
                                 sharey="all",
                                 figsize=(12 * ncols, 19))

        for ax, x in zip(axes, columns):
            color = "#000000"
            if palette is not None and x in palette:
                color = palette[x]

            sns.despine(fig=fig, ax=ax)
            sns.violinplot(x=x,
                           y=y,
                           data=df,
                           order=order,
                           color=color,
                           cut=0,
                           dodge=False,
                           ax=ax)

            ax.set_ylabel("",
                          fontsize=20,
                          fontweight='bold')

        # Add the main title.
        fig.suptitle(title,
                     fontsize=30,
                     color="#000000",
                     weight='bold')

        for extension in self.extensions:
            filename = "cell_fraction_violinplot.{}".format(extension)
            if outdir is not None:
                outpath = os.path.join(outdir, filename)
            else:
                outpath = filename
            fig.savefig(outpath)
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Cell fraction file: {}".format(self.cf_path))
        print("  > Sample-to-dataset file: {}".format(self.std_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > SNPs: {}".format(", ".join(self.snps)))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
