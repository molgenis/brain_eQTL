"""
File:         inter_clustermap.py
Created:      2020/03/16
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
import sys
import os

# Third party imports.
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class InterClusterMap:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'inter_clustermap')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.inter_cov_zscore_df = dataset.get_inter_cov_zscore_df()
        #self.inter_tech_cov_zscore_df = dataset.get_inter_tech_cov_zscore_df()
        self.inter_cov_inter_tvalue_df = dataset.get_inter_cov_inter_tvalue_df()
        #self.inter_tech_cov_inter_tvalue_df = dataset.get_inter_tech_cov_inter_tvalue_df()
        self.cmap = dataset.get_diverging_cmap()

    def start(self):
        sys.setrecursionlimit(10000)

        print("Plotting interaction clustermap")
        self.print_arguments()

        self.visualize_matrix(self.inter_cov_zscore_df, self.outdir, outfile_prefix="cov_zscore")
        #self.visualize_matrix(self.inter_tech_cov_zscore_df, self.outdir, outfile_prefix="tech_cov_zscore")
        self.visualize_matrix(self.inter_cov_inter_tvalue_df, self.outdir, outfile_prefix="cov_tvalue")
        #self.visualize_matrix(self.inter_tech_cov_inter_tvalue_df, self.outdir, outfile_prefix="tech_cov_tvalue")

    def visualize_matrix(self, df, outdir, outfile_prefix="", vmin=None, vmax=None):
        clean_df = df.dropna(axis=1)

        self.plot(df=clean_df, cmap=self.cmap, outdir=outdir,
                  extension=self.extension,
                  outfile_prefix=outfile_prefix, vmin=vmin, vmax=vmax)

        if vmin is None:
            vmin = clean_df.values.min()
        if vmax is None:
            vmax = clean_df.values.max()

        self.plot_colorbar(vmin=vmin, vmax=vmax, cmap=self.cmap, outdir=outdir,
                           name_prefix=outfile_prefix, extension=self.extension)

    @staticmethod
    def plot(df, cmap, outdir, extension, outfile_prefix="", vmin=None, vmax=None):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap=cmap,
                           vmin=vmin, vmax=vmax,
                           yticklabels=True, xticklabels=False,
                           dendrogram_ratio=(.1, .1),
                           figsize=(12, (.2 * (len(df.index)))))
        g.cax.set_visible(False)
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=10))
        g.savefig(os.path.join(outdir, "{}_clustermap.{}".format(outfile_prefix, extension)))
        plt.close()

    @staticmethod
    def plot_colorbar(vmin, vmax, cmap, outdir, name_prefix, extension):
        a = np.array([[vmin, 0, vmax]])
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(3, 9))
        img = sns.heatmap(a, center=0, vmin=vmin, vmax=vmax, cmap=cmap,
                          ax=ax1, cbar_ax=ax2)
        ax1.remove()
        plt.savefig(os.path.join(outdir, "{}_colorbar.{}".format(name_prefix, extension)), bbox_inches='tight')

    def print_arguments(self):
        print("Arguments:")
        print("  > Cov interaction z-score matrix shape: {}".format(self.inter_cov_zscore_df.shape))
        #print("  > Tech.cov interaction z-score matrix shape: {}".format(self.inter_tech_cov_zscore_df.shape))
        print("  > Cov interaction t-value matrix shape: {}".format(self.inter_cov_inter_tvalue_df.shape))
        #print("  > Tech.cov interaction t-value matrix shape: {}".format(self.inter_tech_cov_inter_tvalue_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
