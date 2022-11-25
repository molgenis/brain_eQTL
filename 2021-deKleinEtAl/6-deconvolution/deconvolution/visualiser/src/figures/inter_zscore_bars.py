"""
File:         inter_zscores_bars.py
Created:      2020/03/16
Last Changed: 2020/06/03
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
import math
import os

# Third party imports.
import numpy as np
import scipy.stats as stats
from colour import Color
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class InterZscoreBars:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'inter_zscore_bars')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.inter_cov_zscore_df = dataset.get_inter_cov_zscore_df()
        self.inter_tech_cov_zscore_df = dataset.get_inter_tech_cov_zscore_df()
        self.z_score_cutoff = dataset.get_significance_cutoff()

    def start(self):
        print("Plotting interaction matrix z-scores as barplot.")
        self.print_arguments()

        print("Plotting covariates of interest")
        self.plot_bars(self.inter_cov_zscore_df, fontsize=5, outdir="covariates")

        print("Plotting technical covariates")
        self.plot_bars(self.inter_tech_cov_zscore_df, fontsize=10, outdir="technical_covariates")

    def plot_bars(self, df, outdir, fontsize=10):
        full_outdir = os.path.join(self.outdir, outdir)
        if not os.path.exists(full_outdir):
            os.makedirs(full_outdir)

        print("\tPlotting negative z-scores.")
        data = self.get_counts(df, lower_cutoff=0, upper_cutoff=np.inf)
        self.create_plots(data, full_outdir, fontsize,
                          "Negative ", "[<0]")

        print("\tPlotting positive z-scores.")
        data = self.get_counts(df, lower_cutoff=-np.inf, upper_cutoff=0)
        self.create_plots(data, full_outdir, fontsize,
                          "Positive ", "[>0]")

        pos_z_score_cutoff = abs(self.z_score_cutoff)
        neg_z_score_cutoff = -1 * pos_z_score_cutoff

        print("\tPlotting significant z-scores.")
        data = self.get_counts(df,
                               lower_cutoff=-np.inf,
                               upper_cutoff=pos_z_score_cutoff)
        self.create_plots(data, full_outdir, fontsize, "Significant ",
                          "[>{:.2f}]".format(pos_z_score_cutoff))

    def get_counts(self, data, lower_cutoff=0, upper_cutoff=0):
        df = data.copy()
        df[(df > lower_cutoff) & (df < upper_cutoff)] = 0
        sums = df.pow(2).sum(axis=1).to_frame().reset_index()
        sums.columns = ["index", "counts"]
        sums["color"] = self.create_color_map(sums["counts"])
        sums.sort_values(by=['counts'], ascending=False, inplace=True)

        return sums

    def create_plots(self, data, outdir, fontsize=10,
                     title_prefix="", subtitle_suffix=""):
        max_val = max(data["counts"])
        self.plot(data.copy(), max_val, outdir, self.extension, fontsize=fontsize,
                  title_prefix=title_prefix, subtitle_suffix=subtitle_suffix)
        self.plot(data.copy(), max_val, outdir, self.extension, top=10,
                  fontsize=14, title_prefix=title_prefix,
                  subtitle_suffix=subtitle_suffix)
        self.plot(data.copy(), max_val, outdir, self.extension, bottom=10,
                  fontsize=14, title_prefix=title_prefix,
                  subtitle_suffix=subtitle_suffix)

    @staticmethod
    def create_color_map(values, size=101):
        """
        """
        palette = list(Color("#8ABBDB").range_to(Color("#344A5A"), size))
        colors = [str(x).upper() for x in palette]

        start = 0
        stop = math.ceil(max(values))
        if start == stop:
            return ["#8ABBDB"] * len(values)
        step = (stop - start) / size

        color_list = []
        for value in values:
            index = int(value / step) - 1
            color_list.append(colors[index])
        return color_list

    @staticmethod
    def plot(df, max_val, outdir, extension, top=None, bottom=None,
             fontsize=10, title_prefix="", subtitle_suffix=""):
        # Prepare the data.
        subtitle = ''
        file_appendix = ''
        if top is not None:
            df = df.loc[df["counts"] > 0, :]
            if len(df.index) == 0:
                return
            df = df.head(top)
            subtitle = 'top {} covariates'.format(top)
            file_appendix = '_top{}'.format(top)
        if bottom is not None:
            df = df.loc[df["counts"] > 0, :]
            if len(df.index) == 0:
                return
            df = df.tail(bottom)
            subtitle = 'bottom {} covariates'.format(bottom)
            file_appendix = '_bottom{}'.format(bottom)

        # Plot.
        sns.set(rc={'figure.figsize': (12, max(0.1 * df.shape[0], 9))})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        major_ticks = 10 ** (math.floor(math.log10(max(max_val, 100))))
        minor_ticks = int(major_ticks / 2)
        for i in range(0, int(max(df["counts"])) + (1 * major_ticks), minor_ticks):
            alpha = 0.025
            if i % major_ticks == 0:
                alpha = 0.15
            ax.axvline(i, ls='-', color="#000000", alpha=alpha, zorder=-1)

        g = sns.barplot(x="counts", y="index", data=df,
                        palette=df["color"], orient="h")

        g.text(0.5, 1.05,
               '{}Covariate Interactions'.format(title_prefix),
               fontsize=22, weight='bold', ha='center', va='bottom',
               transform=ax.transAxes)
        g.text(0.5, 1.02, "{} {}".format(subtitle, subtitle_suffix),
               fontsize=12, alpha=0.75, ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('covariate',
                     fontsize=18,
                     fontweight='bold')
        g.set_xlabel('sum(z-score^2)',
                     fontsize=18,
                     fontweight='bold')
        ax.tick_params(labelsize=fontsize)
        ax.set_yticks(range(len(df.index)))
        ax.set_yticklabels(df["index"])
        ax.set_xlim(0, max_val)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}_cov_zscores_barplot{}.{}".format(title_prefix.replace(" ", "").lower(),
                                                                       file_appendix, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Cov interaction z-score matrix shape: {}".format(self.inter_cov_zscore_df.shape))
        print("  > Tech.cov interaction z-score matrix shape: {}".format(self.inter_tech_cov_zscore_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
