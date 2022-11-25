"""
File:         covariate_clustermap.py
Created:      2020/06/02
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
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class CovariateClustermap:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'covariate_clustermap')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.cov_df = dataset.get_cov_df()
        self.cmap = dataset.get_diverging_cmap()

    def start(self):
        print("Plotting convariate comparison.")
        self.print_arguments()
        norm_df = self.normalize(self.cov_df)
        self.plot(norm_df, self.cmap, self.outdir, self.extension)

    @staticmethod
    def normalize(df):
        out_df = df.copy()
        out_df = out_df.loc[out_df.std(axis=1) > 0, :]
        out_df = out_df.subtract(out_df.mean(axis=1), axis=0).divide(out_df.std(axis=1), axis=0)
        return out_df

    @staticmethod
    def plot(df, cmap, outdir, extension):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0.5, cmap=cmap,
                           yticklabels=True, xticklabels=False,
                           row_cluster=False, col_cluster=True,
                           figsize=(12, (.2 * (len(df.index)))))
        plt.setp(g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=10))
        g.savefig(os.path.join(outdir, "covariate_clustermap.{}".format(extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
