"""
File:         inter_pvalue_boxplot.py
Created:      2020/05/12
Last Changed: 2020/06/19
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
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class InterPvalueBoxplot:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'inter_pvalue_boxplot')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.inter_cov_pvalue_df = dataset.get_inter_cov_pvalue_df()
        self.inter_tech_cov_pvalue_df = dataset.get_inter_tech_cov_pvalue_df()
        self.z_score_cutoff = dataset.get_significance_cutoff()

    @staticmethod
    def percentile(n):
        def percentile_(x):
            return np.percentile(x, n)

        percentile_.__name__ = 'percentile_%s' % n
        return percentile_

    def start(self):
        print("Plotting p-values as boxplot.")
        self.print_arguments()

        cov_df = self.inter_cov_pvalue_df.T.melt()
        cov_df["type"] = "Covariate"
        tech_cov_df = self.inter_tech_cov_pvalue_df.T.melt()
        tech_cov_df["type"] = "Technical Covariate"

        df = pd.concat([cov_df, tech_cov_df], axis=0)
        df.dropna(inplace=True)

        group_df = df.copy()
        group_df = group_df.groupby("variable").agg({"value": ["median", self.percentile(25)]})
        group_df.columns = ['median', 'percentile(25)']
        group_df.sort_values(['median', 'percentile(25)'], ascending=[True, True], inplace=True)
        print(group_df)

        interest = group_df.loc[(group_df["median"] == 1.0) & (group_df["percentile(25)"] > 0.99), :].index

        for variable in interest:
            tmp = df.loc[df["variable"] == variable, :].copy()
            all = tmp.shape[0]
            tmp2 = tmp.loc[tmp["value"] == 1.0, :].copy()
            one = tmp2.shape[0]
            print("Variable: {}\tall: {}\tnone: {}\t[{:.2f}%]".format(variable, all, one, ((100 / all) * one)))

        self.plot(df, group_df.index, self.outdir, self.extension)

    @staticmethod
    def plot(data, order, outdir, extension):
        sns.set(rc={'figure.figsize': (
        12, max(0.2 * len(data["variable"].unique()), 9))})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        for i in range(0, 100, 5):
            alpha = 0.025
            if i % 10 == 0:
                alpha = 0.15
            ax.axvline(i / 100, ls='-', color="#000000", alpha=alpha, zorder=-1)

        sns.boxplot(x="value", y="variable", hue="type", data=data,
                    palette={"Covariate": "#0072B2",
                             "Technical Covariate": "#D55E00"},
                    order=order, showfliers=False)

        ax.text(0.5, 1.02,
                'Covariate Interactions P-Values',
                fontsize=22, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.set_ylabel('',
                      fontsize=18,
                      fontweight='bold')
        ax.set_xlabel('p-value',
                      fontsize=12,
                      fontweight='bold')
        ax.tick_params(labelsize=10)
        ax.set_xlim(0, 1)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "pvalues_boxplot.{}".format(extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Cov interaction p-value matrix shape: {}".format(
            self.inter_cov_pvalue_df.shape))
        print("  > Tech.cov interaction p-value matrix shape: {}".format(
            self.inter_tech_cov_pvalue_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
