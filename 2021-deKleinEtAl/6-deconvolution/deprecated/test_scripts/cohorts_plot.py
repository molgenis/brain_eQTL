#!/usr/bin/env python3

"""
File:         cohorts_plot.py
Created:      2020/04/29
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
__program__ = "Cohorts Plot"
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
        self.cohorts = ["AMPAD-MSBB-V2-AFR", "CMC-AFR", "LIBD_1M-AFR",
                        "LIBD_h650-AFR", "AMPAD-MAYO-V2-EUR",
                        "AMPAD-MSBB-V2-EUR", "AMPAD-ROSMAP-V2-EUR",
                        "BrainGVEX-V2-EUR", "CMC-EUR", "GTEx-EUR", "GVEx",
                        "LIBD_1M-EUR", "LIBD_h650-EUR", "NABEC-H550-EUR",
                        "NABEC-H610-EUR", "TargetALS-EUR", "UCLA_ASD-EUR",
                        "ENA-EU"]
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Loading file.")
        df = pd.read_csv(self.cov_path, sep="\t", header=0, index_col=0)
        df = df.loc[df.index.isin(self.cohorts), :]

        print("Count samples.")
        counts = df.sum(axis=1).to_frame()

        eu_data = []
        afr_data = []
        other_data = []
        for index, row in counts.iterrows():
            if "EU" in index or "EUR" in index:
                new_index = index.replace("-EUR", "").replace("-EU", "")
                eu_data.append([new_index] + list(row))
            elif "AFR" in index:
                new_index = index.replace("-AFR", "")
                afr_data.append([new_index] + list(row))
            else:
                other_data.append([index] + list(row))

        eu_df = pd.DataFrame(eu_data, columns=["cohort", "EU"])
        eu_df.set_index("cohort", inplace=True)
        afr_df = pd.DataFrame(afr_data, columns=["cohort", "AFR"])
        afr_df.set_index("cohort", inplace=True)
        oth_df = pd.DataFrame(other_data, columns=["cohort", "Other"])
        oth_df.set_index("cohort", inplace=True)

        tmp = oth_df.merge(eu_df, how="outer", left_index=True, right_index=True)
        df = tmp.merge(afr_df, how="outer", left_index=True, right_index=True)
        del tmp
        df.fillna(0, inplace=True)
        col_sums = df.sum()
        df["sum"] = df.sum(axis=1)
        df.index = [x.replace("-", " ").replace("_", " ").replace("V2", "(V2)") for x in df.index]
        df.sort_index(inplace=True)
        print(df)

        group_fraction = {}
        for index, value in col_sums.items():
            group_fraction[index] = (value / sum(col_sums)) * 100

        pallete = {"Other": "#808080",
                   "EU": "cornflowerblue",
                   "AFR": "firebrick"}

        print("Plotting.")
        sns.set(rc={'figure.figsize': (10, 7.5)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        sns.barplot(x="sum", y=df.index, color="firebrick",
                    data=df, orient="h", ax=ax)
        sns.barplot(x="EU", y=df.index, color="cornflowerblue",
                    data=df, orient="h", ax=ax)
        g = sns.barplot(x="Other", y=df.index, color="#808080",
                        data=df, orient="h", ax=ax)

        g.text(0.5, 1.025,
               'Cohort sizes',
               fontsize=20, weight='bold', ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('',
                     fontsize=16,
                     fontweight='bold')
        g.set_xlabel('sample size',
                     fontsize=16,
                     fontweight='bold')

        for i in range(0, 900, 50):
            alpha = 0.025
            if i % 100 == 0:
                alpha = 0.15
            ax.axvline(i, ls='-', color="#000000", alpha=alpha, zorder=-1)

        handles = []
        for name, color in pallete.items():
            handles.append(mpatches.Patch(color=color, label="{} [{:.0f}%]".format(name, group_fraction[name])))
        ax.legend(handles=handles)

        ax.tick_params(labelsize=12)
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "cohort_sample_sizes.png"))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
