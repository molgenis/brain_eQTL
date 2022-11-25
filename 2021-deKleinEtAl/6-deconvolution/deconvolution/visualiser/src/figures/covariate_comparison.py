"""
File:         covariate_comparison.py
Created:      2020/04/15
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
import scipy.stats as stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir, p_value_to_symbol


class CovariateComparison:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'covariate_comparison')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.groups = dataset.get_groups()
        self.cov_df = dataset.get_cov_df()
        self.cmap = dataset.get_diverging_cmap()

    def start(self):
        print("Plotting convariate comparison.")
        self.print_arguments()
        corr_df, pval_df = self.correlate(self.cov_df)
        self.plot(corr_df, pval_df, self.groups, self.cmap, self.outdir, self.extension)

    @staticmethod
    def correlate(df):
        print("Creating the correlation matrix.")
        corr_df = pd.DataFrame(np.nan, index=df.index, columns=df.index)
        pval_df = pd.DataFrame("", index=df.index, columns=df.index)
        for i, row1 in enumerate(df.index):
            for j, row2 in enumerate(df.index):
                if i >= j:
                    coef, p = stats.spearmanr(df.loc[row1, :], df.loc[row2, :])
                    corr_df.loc[row1, row2] = coef
                    #pval_df.loc[row1, row2] = p_value_to_symbol(p)
                    pval_df.loc[row1, row2] = "{:.2e}".format(p)
                    if (coef == 1.0) and (row1 != row2):
                        print("{} = {}".format(row1, row2))

        return corr_df, pval_df

    @staticmethod
    def plot(corr_df, pval_df, groups, cmap, outdir, extension):
        print("Plotting")

        gridspec_kw = {"height_ratios": [x[1] - x[0] for x in groups],
                       "width_ratios": [x[1] - x[0] for x in groups]}

        norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)

        heatmapkws = dict(square=False, cbar=False, cmap=cmap, fmt='',
                          linewidths=1.0, center=0, vmin=-1, vmax=1,
                          annot_kws={"size": 12, "color": "#808080"})

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=len(groups), nrows=len(groups),
                                 figsize=(101, 73), gridspec_kw=gridspec_kw)
        plt.subplots_adjust(left=0.2, right=0.87, bottom=0.2, top=0.95,
                            wspace=0.1, hspace=0.1)

        for i in range(len(groups)):
            (ra, rb, ylabel, yremove) = groups[i]
            for j in range(len(groups)):
                (ca, cb, xlabel, xremove) = groups[j]
                ax = axes[i, j]

                if i >= j:
                    print("\tPlotting axes[{}, {}]".format(i, j))
                    xticklabels = False
                    if i == (len(groups) - 1):
                        xticklabels = True
                    yticklabels = False
                    if j == 0:
                        yticklabels = True

                    data_subset = corr_df.iloc[ra:rb, ca:cb]
                    anmnot_subset = pval_df.iloc[ra:rb, ca:cb]

                    sns.heatmap(data_subset,
                                annot=anmnot_subset,
                                xticklabels=xticklabels,
                                yticklabels=yticklabels,
                                ax=ax,
                                **heatmapkws)

                    if xticklabels:
                        new_xticks = [x.replace(xremove, '').replace("_", " ") for x in data_subset.columns]
                        ax.set_xticklabels(new_xticks, fontsize=14, rotation=90)
                        ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')

                    if yticklabels:
                        new_yticks = [y.replace(yremove, '').replace("_", " ") for y in data_subset.index]
                        ax.set_yticklabels(new_yticks, fontsize=14, rotation=0)
                        ax.set_ylabel(ylabel, fontsize=16, fontweight='bold')

                else:
                    ax.set_axis_off()

        cax = fig.add_axes([0.9, 0.2, 0.01, 0.7])
        sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, cax=cax)

        fig.align_ylabels(axes[:, 0])
        fig.align_xlabels(axes[len(groups) - 1, :])
        fig.suptitle('Covariate Correlations', fontsize=40, fontweight='bold')
        fig.savefig(os.path.join(outdir, "covariate_comparison.{}".format(extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Groups: {}".format(self.groups))
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
