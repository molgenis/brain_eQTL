"""
File:         combine_and_plot.py
Created:      2020/03/30
Last Changed: 2020/06/09
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
from bisect import bisect_left
from colour import Color
from itertools import groupby, count
import pickle
import math
import glob
import time
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.local_settings import LocalSettings
from general.df_utilities import save_dataframe, get_basename


class CombineAndPlot:
    """
    Main: this class combines the output of different jobs of the manager into
        one result.
    """
    def __init__(self, name, settings_file):
        """
        Initializer of the class.

        :param name: string, the name of the base input/ouput directory.
        :param settings_file: string, the name of the settings file.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Safe arguments.
        self.n_permutations = settings.get_setting("n_permutations")

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir, name)

        # Get the needed settings.
        self.cov_outdir = settings.get_setting("covariates_folder")
        self.tech_cov_outdir = settings.get_setting("technical_covariates_folder")
        self.pvalues_outfile = settings.get_setting("actual_pvalues_pickle_filename")
        self.snp_tvalues_outfile = settings.get_setting("snp_tvalues_pickle_filename")
        self.inter_tvalues_outfile = settings.get_setting("inter_tvalues_pickle_filename")
        self.perm_pvalues_outfile = settings.get_setting("permuted_pvalues_pickle_filename")

    def start(self):
        """
        Method to start the combiner.
        """
        print("Starting interaction analyser - combine and plot.")
        self.print_arguments()

        # Start the timer.
        start_time = time.time()

        for outdir in [self.cov_outdir, self.tech_cov_outdir]:
            full_outdir = os.path.join(self.outdir, outdir)
            if not os.path.exists(full_outdir):
                print("Error, output directory does not exist.")

            self.work(full_outdir)

        # Print the time.
        run_time_min, run_time_sec = divmod(time.time() - start_time, 60)
        run_time_hour, run_time_min = divmod(run_time_min, 60)
        print("finished in  {} hour(s), {} minute(s) and "
              "{} second(s).".format(int(run_time_hour),
                                     int(run_time_min),
                                     int(run_time_sec)), flush=True)

    def work(self, workdir):

        # pvalue_df = pd.read_csv(os.path.join(workdir, "pvalue_table.txt.gz"), sep="\t", header=0, index_col=0)
        # perm_fdr_df = pd.read_csv(os.path.join(workdir, "perm_fdr_table.txt.gz"), sep="\t", header=0, index_col=0)
        # bh_fdr_df = pd.read_csv(os.path.join(workdir, "bh_fdr_table.txt.gz"), sep="\t", header=0, index_col=0)

        # Combine the pickle files.
        print("Loading pvalue data.", flush=True)
        pcolumns, pvalues_data = self.combine_pickles(workdir,
                                                      self.pvalues_outfile,
                                                      columns=True)

        # Create a pandas dataframe from the nested list.
        print("Creating p-values dataframe.", flush=True)
        pvalue_df = self.create_df(pvalues_data, pcolumns)
        save_dataframe(df=pvalue_df,
                       outpath=os.path.join(workdir,
                                            "pvalue_table.txt.gz"),
                       header=True, index=True)

        # pvalue_df = pd.read_csv(os.path.join(workdir, "pvalue_table.txt.gz"),
        #                         sep="\t", header=0, index_col=0)
        # with open(os.path.join(workdir, "perm_pvalues.pkl"), "rb") as f:
        #     perm_pvalues = pickle.load(f)
        # f.close()

        # Get the pvalues from the dataframe.
        pvalues = pvalue_df.melt()["value"].values

        print("Loading permutation pvalue data.", flush=True)
        _, perm_pvalues = self.combine_pickles(workdir,
                                               self.perm_pvalues_outfile)
        # with open(os.path.join(workdir, "perm_pvalues.pkl"), "wb") as f:
        #      pickle.dump(perm_pvalues, f)
        # f.close()

        # Visualise distributions.
        print("Visualizing distributions.", flush=True)
        self.plot_distributions(perm_pvalues, pvalues, workdir)

        # return

        print("Loading SNP tvalue data.", flush=True)
        snp_tcolumns, snp_tvalues_data = self.combine_pickles(workdir,
                                                              self.snp_tvalues_outfile,
                                                              columns=True)

        # Create a pandas dataframe from the nested list.
        print("Creating SNP t-values dataframe.", flush=True)
        snp_tvalue_df = self.create_df(snp_tvalues_data, snp_tcolumns)
        save_dataframe(df=snp_tvalue_df,
                       outpath=os.path.join(workdir, "snp_tvalue_table.txt.gz"),
                       header=True, index=True)

        print("Loading inter tvalue data.", flush=True)
        inter_tcolumns, inter_tvalues_data = self.combine_pickles(workdir,
                                                                  self.inter_tvalues_outfile,
                                                                  columns=True)

        # Create a pandas dataframe from the nested list.
        print("Creating inter t-values dataframe.", flush=True)
        inter_tvalue_df = self.create_df(inter_tvalues_data, inter_tcolumns)
        save_dataframe(df=inter_tvalue_df,
                       outpath=os.path.join(workdir, "inter_tvalue_table.txt.gz"),
                       header=True, index=True)

        # Create a dataframe with z-scores.
        print("Creating Z-score dataframe.", flush=True)
        zscore_df = self.create_zscore_df(pvalue_df)
        save_dataframe(df=zscore_df,
                       outpath=os.path.join(workdir,
                                            "interaction_table.txt.gz"),
                       header=True, index=True)

        # Sort the lists.
        print("Sorting p-values.", flush=True)
        perm_pvalues = sorted(perm_pvalues)
        pvalues = sorted(pvalues)

        # Create the FDR dataframes.
        print("Creating permutation FDR dataframe.", flush=True)
        perm_fdr_df, perm_cutoff = self.create_perm_fdr_df(pvalue_df,
                                                           pvalues,
                                                           perm_pvalues,
                                                           self.n_permutations)
        perm_n_signif = self.count_n_significant(pvalues, perm_cutoff)
        print("\tPermutation FDR: {} p-values < signif. cutoff "
              "{:.2e} [{:.2f}%]".format(perm_n_signif, perm_cutoff,
                                        (100 / len(pvalues)) * perm_n_signif))
        # Write the output file.
        save_dataframe(df=perm_fdr_df,
                       outpath=os.path.join(workdir,
                                            "perm_fdr_table.txt.gz"),
                       header=True, index=True)

        print("Creating Benjamini-Hochberg FDR dataframe.", flush=True)
        bh_fdr_df, bh_cutoff = self.create_bh_fdr_df(pvalue_df, pvalues)
        bh_n_signif = self.count_n_significant(pvalues, bh_cutoff)
        print("\tBH FDR: {} p-values < signif. cutoff "
              "{:.2e} [{:.2f}%]".format(bh_n_signif, bh_cutoff,
                                        (100 / len(pvalues)) * bh_n_signif))
        save_dataframe(df=bh_fdr_df,
                       outpath=os.path.join(workdir,
                                            "bh_fdr_table.txt.gz"),
                       header=True, index=True)

        #return

        # Compare the two pvalue scores.
        print("Creating score visualisation [1/2].", flush=True)
        self.compare_pvalue_scores(pvalue_df, perm_fdr_df, bh_fdr_df,
                                   workdir)

    @staticmethod
    def combine_pickles(indir, filename, columns=False):
        """
        Method for combining the pickled lists into one nested list.

        :param indir: string, the input directory containing the pickle files.
        :param filename: string, the prefix name of the input file.
        :param columns: boolean, whether or not each pickle file has a column.
        :return col_list: list, the columns of the content.
        :return data: list, a nested list with the combined content.
        """
        # Declare variables.
        col_list = None
        data = []

        # Combine the found files.
        for i, fpath in enumerate(glob.glob(os.path.join(indir, filename,
                                                         filename + "*.pkl"))):
            with open(fpath, "rb") as f:
                try:
                    content = pickle.load(f)
                    if columns and i == 0:
                        col_list = content[0]
                    data.extend(content[1:])
                    # print("\tLoaded list: {} with length: {}".format(
                    #     get_basename(fpath), len(content)))
                except EOFError:
                    print("\tEOFError in: {} ".format(get_basename(fpath)))
            f.close()

        return col_list, data

    def create_df(self, data, columns):
        """
        Method for creating a pandas dataframe from a nested list.

        :param data: list, a nested list with the data.
        :param columns: list, the column names.
        :return df: DataFrame, the created pandas dataframe.
        """
        # Crate the data frame.
        df = pd.DataFrame(data, columns=columns)
        order_col = df.columns[0]
        df.sort_values(by=order_col, inplace=True)
        print("\tInput shape: {}".format(df.shape))

        # Check what eQTLs are present, missing, duplicated.
        reference = set(np.arange(df[order_col].min(), df[order_col].max() + 1))
        present = set(df[order_col])
        print("\tPresent indices: {}".format(self.group_consecutive_numbers(present)))
        missing = list(set(df[order_col]).symmetric_difference(reference))
        print("\tMissing indices: {}".format(self.group_consecutive_numbers(missing)))
        duplicated = list(df.loc[df.duplicated([order_col]), order_col].values)
        print("\tDuplicate indices: {}".format(self.group_consecutive_numbers(duplicated)))

        # Remove duplicate eQTLs.
        df.drop_duplicates(subset=order_col, keep="first", inplace=True)

        # Insert missing eQTLs.
        if len(missing) > 0:
            print("\tInserting missing indices")
            missing_df = pd.Series(missing).to_frame()
            missing_df.columns = [order_col]
            df = pd.concat([df, missing_df], axis=0).reset_index(drop=True)
            df.sort_values(by=order_col, inplace=True)

        # Remove the eQTL order column and set the SNPName as index.
        df.drop([order_col], axis=1, inplace=True)
        df.set_index("-", inplace=True)
        df = df.T
        print("\tOutput shape: {}".format(df.shape))

        return df

    @staticmethod
    def group_consecutive_numbers(numbers):
        groups = groupby(numbers, key=lambda item, c=count(): item - next(c))
        tmp = [list(g) for k, g in groups]
        return [str(x[0]) if len(x) == 1 else "{}-{}".format(x[0], x[-1]) for x in tmp]

    def create_zscore_df(self, df):
        """
        Method for converting a dataframe of p-values to z-scores.

        :param df: DataFrame, a dataframe containing p-values.
        :return zscore_df: DataFrame, a dataframe containing z-scores
        """
        count = 1
        total = df.shape[0] * df.shape[1]

        zscore_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                if (count == 1) or (count % int(total / 10) == 0) or (count == total):
                    print("\tworking on {}/{} [{:.2f}%]".format(count,
                                                                total,
                                                                (100 / total) * count))
                pvalue = df.iloc[row_index, col_index]
                zscore = self.get_z_score(pvalue)
                zscore_df.iloc[row_index, col_index] = zscore
                count += 1

        return zscore_df

    @staticmethod
    def get_z_score(p_value):
        """
        Method for converting a p-value to a z-score.

        :param p-value: float, the p-value to convert.
        :param z-score: float, the corresponding z-score.
        """
        if p_value > (1.0 - 1e-16):
            p_value = (1.0 - 1e-16)
        if p_value < 1e-323:
            p_value = 1e-323

        # The lower and upper limit of stats.norm.sf
        # stats.norm.isf((1 - 1e-16)) = -8.209536151601387
        # stats.norm.isf(1e-323) = 38.44939448087599
        return stats.norm.isf(p_value)

    @staticmethod
    def plot_distributions(perm_pvalues, pvalues, outdir):
        """
        Method for visualizing the distribution of the null and alternative
        p-values.

        :param perm_pvalues: list, the sorted null model p-values.
        :param pvalues: list, the sorted alternative model p-values.
        :param outdir: string, the output directory for the image.
        """
        # Create bins.
        bins = np.linspace(0, 1, 50)
        perm_pvalues_bins = pd.cut(perm_pvalues, bins=bins).value_counts().to_frame()
        pvalues_bins = pd.cut(pvalues, bins=bins).value_counts().to_frame()
        df = perm_pvalues_bins.merge(pvalues_bins, left_index=True, right_index=True)
        df.columns = ["perm_pvalues", "pvalues"]
        df = df.divide(df.sum())
        lowest_row = df[df["pvalues"] == df["pvalues"].min()]
        scaler = lowest_row.iloc[0, 1] / lowest_row.iloc[0, 0]
        df["baseline"] = (df["perm_pvalues"] * scaler)
        df["index"] = (bins[:-1] + bins[1:]) / 2

        # Set the step sizes.
        multiplication = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1]
        step_sizes = sorted([(x * y) for x in [1, 2, 2.5, 5]
                             for y in multiplication])

        # Plot.
        sns.set(rc={'figure.figsize': (18, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

        sns.barplot(x="index", y="perm_pvalues", color='cornflowerblue', data=df, ax=ax1)
        sns.barplot(x="index", y="pvalues", color='firebrick', data=df, ax=ax2)
        sns.barplot(x="index", y="pvalues", color='firebrick', data=df, ax=ax3)
        sns.barplot(x="index", y="baseline", color='cornflowerblue', data=df, ax=ax3)

        for ax, title, ymax in zip([ax1, ax2, ax3],
                                   ["Null", "Alternative", "Combined"],
                                   [df["perm_pvalues"].max(), df["pvalues"].max(), df["pvalues"].max()]):
            sns.despine(fig=fig, ax=ax)

            # Determine the step size.
            distances = pd.Series([abs(x - (ymax / 10)) for x in step_sizes])
            index = distances.idxmin()
            if isinstance(index, list):
                index = index[0]
            step_size = int(step_sizes[index] * 1000)

            for i in range(0, 1001, step_size):
                alpha = 0.025
                if i % (step_size * 2) == 0:
                    alpha = 0.15
                if i / 1000 < ymax:
                    ax.axhline(i / 1000, ls='-', color="#000000",
                               alpha=alpha, zorder=-1)

            ax.text(0.5, 1.02,
                    '{} Distribution'.format(title),
                    fontsize=16, weight='bold', ha='center', va='bottom',
                    transform=ax.transAxes)
            ax.set_xlabel('p-value',
                          fontsize=12,
                          fontweight='bold')
            ax.set_ylabel('proportion',
                          fontsize=12,
                          fontweight='bold')

            labels = ax.get_xticklabels()
            new_labels = []
            for i, label in enumerate(labels):
                if (i == 0) or (i % 3 == 0) or (i == (len(labels) - 1)):
                    new_labels.append("{:.2f}".format(float(label.get_text())))
                else:
                    new_labels.append("")
            ax.set_xticklabels(new_labels, rotation=45)

        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "pvalues_distributions.png"))
        plt.close()

    @staticmethod
    def create_perm_fdr_df(df, pvalues, perm_pvalues, n_perm):
        """
        Method for creating the permutation False Discovery Rate dataframe.

        FDR = (permutation rank / number of permutations) / actual rank

        :param df: DataFrame, the alternative p-value dataframe.
        :param perm_pvalues: list, the sorted null model p-values.
        :param pvalues: list, the sorted alternative model p-values.
        :param n_perm: int, the number of permutations performed.
        :return fdr_df: DataFrame, the permutation FDR dataframe.
        """
        count = 1
        total = df.shape[0] * df.shape[1]

        max_signif_pvalue = -np.inf

        fdr_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                if (count == 1) or (count % int(total / 10) == 0) or (count == total):
                    print("\tworking on {}/{} [{:.2f}%]".format(count,
                                                                total,
                                                                (100 / total) * count))
                pvalue = df.iloc[row_index, col_index]
                rank = bisect_left(pvalues, pvalue)
                perm_rank = bisect_left(perm_pvalues, pvalue)
                if (rank > 0) and (perm_rank > 0):
                    fdr_value = (perm_rank / n_perm) / rank
                    if fdr_value > 1:
                        fdr_value = 1
                else:
                    fdr_value = 0
                fdr_df.iloc[row_index, col_index] = fdr_value

                if fdr_value < 0.05 and pvalue > max_signif_pvalue:
                    max_signif_pvalue = pvalue

                count += 1

        return fdr_df, max_signif_pvalue

    @staticmethod
    def create_bh_fdr_df(df, pvalues):
        """
        Method for creating the Benjamini-Hochberg False Discovery Rate
        dataframe.

        FDR = p-value * (# p-values / rank)

        :param df: DataFrame, the alternative p-value dataframe.
        :param pvalues: list, the sorted alternative model p-values.
        :return fdr_df: DataFrame, the Benjamini-Hochberg FDR dataframe.
        """
        count = 1
        total = df.shape[0] * df.shape[1]

        max_signif_pvalue = -np.inf

        m = np.count_nonzero(~np.isnan(pvalues))
        fdr_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        fdr_values = []
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                if (count == 1) or (count % int(total / 10) == 0) or (count == total):
                    print("\tworking on {}/{} [{:.2f}%]".format(count,
                                                                total,
                                                                (100 / total) * count))
                pvalue = df.iloc[row_index, col_index]
                rank = bisect_left(pvalues, pvalue) + 1
                fdr_value = pvalue * (m / rank)
                if fdr_value > 1:
                    fdr_value = 1
                fdr_values.append((rank, row_index, col_index, fdr_value))
                fdr_df.iloc[row_index, col_index] = fdr_value

                if fdr_value < 0.05 and pvalue > max_signif_pvalue:
                    max_signif_pvalue = pvalue

                count += 1

        # Make sure the BH FDR is a monotome function. This goes through
        # the FDR values backwords and make sure that the next FDR is always
        # lower or equal to the previous.
        fdr_values = sorted(fdr_values, key=lambda x: x[0], reverse=True)
        prev_fdr_value = None
        for (rank, row_index, col_index, fdr_value) in fdr_values:
            if prev_fdr_value is not None and fdr_value > prev_fdr_value:
                fdr_df.iloc[row_index, col_index] = prev_fdr_value
                prev_fdr_value = prev_fdr_value
            else:
                prev_fdr_value = fdr_value

        return fdr_df, max_signif_pvalue

    @staticmethod
    def count_n_significant(sorted_values, threshold):
        """
        Method to count the number of values in a sorted list lower than
        a certain threshhold.

        :param sorted_values: list, sorted list of numbers.
        :param threshold: float, the cutoff to count the n values below.
        :return i: int, the number of values < threshold
        """
        count = 0
        for value in sorted_values:
            if value != np.nan:
                count += 1
                if value > threshold:
                    break

        return count

    def compare_pvalue_scores(self, pvalue_df, perm_fdr, bh_fdr, outdir,
                              a=0.05):
        """
        Method for comparing the different p-value dataframes.

        :param pvalue_df: DataFrame, the alternative p-value dataframe.
        :param perm_fdr: DataFrame, the permutation FDR dataframe.
        :param bh_fdr: DataFrame, the Benjamini-Hochberg FDR dataframe.
        :param outdir: string, the output directory for the image.
        :param a: float, the singificance cut-off.
        """
        # Combine the dataframes into one.
        df = pd.DataFrame({"p-value": pvalue_df.melt()["value"],
                           "permutation FDR": perm_fdr.melt()["value"],
                           "BH FDR": bh_fdr.melt()["value"]})
        df[(df > 1.0) | (df < 0.0)] = np.nan
        df.dropna(how="any", inplace=True)
        df = df.T

        bins = np.linspace(0, 1.0, 15)

        # Plot
        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=3, nrows=4, figsize=(12, 12),
                                 gridspec_kw={"height_ratios": [0.01, 0.33, 0.33, 0.33]})

        for i in range(3):
            ax = axes[0, i]
            ax.set_axis_off()

        for index1, (label1, row1) in enumerate(df.iterrows()):
            for index2, (label2, row2) in enumerate(df.iterrows()):
                ax = axes[index2+1, index1]
                sns.despine(fig=fig, ax=ax)
                print("\tPlotting axes[{}, {}]".format(index2, index1))

                plot_df = pd.DataFrame({label2: row2.copy(), label1: row1.copy()})

                if index2 != index1:
                    if index2 < index1:
                        plot_df = plot_df.loc[plot_df.max(axis=1) <= a, :]

                    # Calculate the lower and upper bound of the axis.
                    xy_min = plot_df.values.min()
                    xy_max = plot_df.values.max()

                    # Calculate the difference between x and y and create a
                    # colormap accordingly.
                    precision = 3
                    diff = round((plot_df[label2] - plot_df[label1]) ** 2, precision)
                    max_col_value = math.ceil(diff.max() * (10 ** precision))
                    colormap = self.create_color_map(length=max_col_value, precision=precision)
                    colors = diff.map(colormap, na_action=Color(rgb=(0.698, 0.133, 0.133)))

                    # Plot.
                    # print("[{}, {}]\tx: {}\ty: {}".format(index2+1, index1, label1, label2))
                    ax.scatter(x=label1, y=label2, data=plot_df, c=colors.values)
                    ax.plot([xy_min, xy_max], [xy_min, xy_max], ls="--", c=".3")

                    if index2 >= index1:
                        for i in range(int(xy_min * 100), int(xy_max * 100) + 10, 10):
                            alpha = 0.025
                            if i % 20 == 0:
                                alpha = 0.15
                            ax.axvline(i / 100, ls='-', color="#000000", alpha=alpha, zorder=-1)
                            ax.axhline(i / 100, ls='-', color="#000000", alpha=alpha, zorder=-1)
                    else:
                        for i in range(0, 65, 5):
                            alpha = 0.025
                            if i % 10 == 0:
                                alpha = 0.15
                            ax.axvline(i / 1000, ls='-', color="#000000", alpha=alpha, zorder=-1)
                            ax.axhline(i / 1000, ls='-', color="#000000", alpha=alpha, zorder=-1)

                    # ax.set_ylim(xy_min, xy_max)
                    # ax.set_xlim(xy_min, xy_max)

                    ax.set_xlabel(label1, fontsize=14, fontweight='bold')
                    ax.set_ylabel(label2, fontsize=14, fontweight='bold')
                    ax.tick_params(labelsize=10)

                else:
                    # print("[{}, {}]\tdata: {}".format(index2 + 1, index1, label1))
                    value_bins = pd.cut(plot_df[label1], bins=bins).value_counts().to_frame()
                    value_bins.columns = ["count"]
                    value_bins = value_bins.divide(value_bins.sum())
                    value_bins["index"] = pd.IntervalIndex(value_bins.index).mid
                    value_bins = value_bins.round(2)

                    sns.barplot(x="index", y="count", color='#696969',
                                data=value_bins, ax=ax)
                    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

                    for i in range(0, 101, 5):
                        alpha = 0.025
                        if i % 10 == 0:
                            alpha = 0.15
                        ax.axhline(i / 100, ls='-', color="#000000",
                                   alpha=alpha, zorder=-1)

                    ax.set_xlabel(label1, fontsize=14, fontweight='bold')
                    ax.set_ylabel("ratio", fontsize=14, fontweight='bold')
                    ax.tick_params(labelsize=10)

        fig.align_ylabels(axes[:, 0])
        fig.align_xlabels(axes[2, :])
        fig.suptitle('Multiple Testing Comparison', fontsize=20, fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "pvalue_vs_fdr_comparison.png"))
        plt.close()

    @staticmethod
    def create_color_map(length, precision):
        """
        Method for creating a color gradient for the p-values.
        """
        palette = list(Color("#D3D3D3").range_to(Color("#000000"), length + 2))
        colors = [x.rgb for x in palette]
        values = [x / (10 ** precision) for x in list(range(0, length + 2))]
        value_color_map = {x: y for x, y in zip(values, colors)}
        return value_color_map

    def print_arguments(self):
        """
        Method for printing the variables of the class.
        """
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("  > N permutations: {}".format(self.n_permutations))
        print(
            "  > Actual P-values output file: {}*.pkl".format(
                self.pvalues_outfile))
        print("  > Permutated P-values output file: {}*.pkl".format(
            self.perm_pvalues_outfile))
        print("")
