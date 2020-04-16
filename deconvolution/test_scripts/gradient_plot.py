#!/usr/bin/env python3

"""
File:         gradient_plot.py
Created:      2020/04/15
Last Changed: 2020/04/16
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
import os

# Third party imports.
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Gradient Plot"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
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
       self.cov_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/output/create_cov_matrix/covariates_table.txt.gz"
       self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Load the covariate file.")
        cov_df = pd.read_csv(self.cov_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.cov_path),
                                      cov_df.shape))

        # Define the axis.
        xaxis = "CellMapNMF_Astrocyte_C1"
        yaxis = "CellMapNNLS_Astrocyte"
        zaxis = ""

        # # Define the axis.
        # xaxis = "Comp5"
        # yaxis = "Comp6"
        # zaxis = "SEX"

        # Get the data.
        xdata = cov_df.loc[xaxis, :]
        ydata = cov_df.loc[yaxis, :]
        # zdata = cov_df.loc[zaxis, :]
        # zdata = zdata.map({0: "Male", 1: "Female"})
        #
        # # Define colors.
        # colormap = {"Male": "#03165E", "Female": "#DC106C"}

        # De-mean.
        xdata = xdata - xdata.mean()
        ydata = ydata - ydata.mean()

        # Plot.
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")

        ax.axvline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)
        ax.axhline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)

        g = sns.scatterplot(x=xdata,
                            y=ydata,
                            # hue=zdata,
                            # palette=colormap,
                            alpha=0.5)

        g.set_title("")
        g.set_ylabel("{} - average".format(yaxis),
                     fontsize=10,
                     fontweight='bold')
        g.set_xlabel("{} - average".format(xaxis),
                     fontsize=10,
                     fontweight='bold')
        # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        # plt.setp(ax.get_legend().get_texts(), fontsize='10')
        # plt.setp(ax.get_legend().get_title(), fontsize='12')
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_{}_{}.png".format(xaxis,
                                                                    yaxis,
                                                                    zaxis)))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
