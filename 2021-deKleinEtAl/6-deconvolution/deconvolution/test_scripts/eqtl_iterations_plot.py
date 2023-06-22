#!/usr/bin/env python3

"""
File:         eqtl_iterations_plot.py
Created:      2020/04/27
Last Changed: 2020/04/28
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
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "eQTL Iterations Plot"
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
        self.input_dir = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/"
        self.subdir_regex = "Iteration*"
        self.filename = "eQTLProbesFDR0.05-ProbeLevel.txt.gz"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Loading files.")
        counts = []
        for subdir in glob.glob(
                os.path.join(self.input_dir, self.subdir_regex)):
            iteration = int(
                os.path.basename(os.path.normpath(subdir)).replace("Iteration",
                                                                   ""))
            filepath = os.path.join(subdir, self.filename)

            if os.path.exists(filepath):
                tmp = pd.read_csv(filepath, sep="\t", header=0, index_col=False)
                counts.append([iteration, tmp.shape[0]])

        print("Creating dataframe(s).")
        df = pd.DataFrame(counts, columns=["iteration", "count"])
        df = df.loc[df["count"] > 0, :]
        df.sort_values(by=['iteration'], inplace=True)

        print(df)

        cutoff = 4
        colors = (['#6495ED'] * cutoff) + (['#808080'] * (df.shape[0] - cutoff))

        print("Plotting.")
        sns.set(rc={'figure.figsize': (10, 7.5)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        g = sns.barplot(x="iteration", y="count", palette=colors, data=df)
        g.set(yscale="log")

        g.text(0.5, 1.025,
               'Number of cis-eQTLs per n$^{th}$ Degree effect',
               fontsize=20, weight='bold', ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('number of cis-eQTLs',
                     fontsize=16,
                     fontweight='bold')
        g.set_xlabel('iteration',
                     fontsize=16,
                     fontweight='bold')

        # for i in range(0, 5, 1):
        #     ax.axhline(10 ** i, ls='-', color="#000000", alpha=0.15, zorder=-1)
        for i in range(0, 13000, 1000):
            alpha = 0.025
            if i % 2000 == 0:
                alpha = 0.15
            ax.axhline(i, ls='-', color="#000000", alpha=alpha, zorder=-1)

        plt.legend(handles=[mpatches.Patch(color='#6495ED', label='included'),
                            mpatches.Patch(color='#808080', label='excluded')])

        ax.tick_params(labelsize=12)
        ax.set_xticks(range(len(df.index)))
        ax.set_xticklabels(df["iteration"])
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "ciseqtls_per_iteration_new.png"))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
