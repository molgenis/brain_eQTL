#!/usr/bin/env python3

"""
File:         replication_plot.py
Created:      2020/06/10
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
import os

# Third party imports.
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from colour import Color
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Replication Plot"
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
       self.eur_afr_result = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/identify_ct_mediated_eqtls/cis_output/all.txt"
       self.eur_result = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/identify_ct_mediated_eqtls/cis_new_output/all.txt"
       self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Load the files.")
        eur_afr_df = pd.read_csv(self.eur_afr_result, sep="\t", header=0, nrows=5256)
        print("\tLoaded dataframe: EUR_AFR "
              "with shape: {}".format(eur_afr_df.shape))

        eur_df = pd.read_csv(self.eur_result, sep="\t", header=0, nrows=893)
        print("\tLoaded dataframe: EUR "
              "with shape: {}".format(eur_df.shape))

        # Create ID column.
        eur_afr_df.index = eur_afr_df["SNPName"] + eur_afr_df["HGNCName"] + eur_afr_df["Covariate"]
        eur_df.index = eur_df["SNPName"] + eur_df["HGNCName"] + eur_df["Covariate"]

        # Find overlap.
        eur_afr_hits = set(eur_afr_df.index)
        eur_hits = set(eur_df.index)
        overlap = eur_afr_hits.intersection(eur_hits)
        print("EUR_AFR: {}\tEUR: {}\tOverlap: {}\t".format(len(eur_afr_hits),
                                                           len(eur_hits),
                                                           len(overlap)))

        # Subset.
        eur_afr_subset = eur_afr_df.loc[overlap, :]
        eur_subset = eur_df.loc[overlap, :]

        # Plot.
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")

        ax.axvline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)
        ax.axhline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)

        g = sns.scatterplot(x=eur_afr_subset["Inter"],
                            y=eur_subset["Inter"],
                            legend=False,
                            alpha=0.5)

        g.set_title("Interaction Replication")
        g.set_ylabel("EUR",
                     fontsize=10,
                     fontweight='bold')
        g.set_xlabel("EUR and AFR",
                     fontsize=10,
                     fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "interaction_replication_plot.png"))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
