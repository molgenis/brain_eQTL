#!/usr/bin/env python3

"""
File:         compare_inter_matrices.py
Created:      2020/04/06
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
import random
import os
import math

# Third party imports.
from matplotlib import pyplot as plt
from matplotlib import collections as mc
import seaborn as sns
import numpy as np
import pandas as pd
import scipy.stats as st

# Local application imports.

# Metadata
__program__ = "XKCD Plot"
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
        self.outdir = str(Path(__file__).parent.parent)
        self.n_points = 100
        self.eqtl_effect_size = 1
        self.eqtl_direction = "positive"
        self.interaction_effect_size = 1
        self.interaction_direction = "positive"
        self.color_map = {0.0: "#ABDCA2", 1.0: "#89A3D1", 2.0: "#F08C72"}

    def start(self):
        plt.xkcd()
        fig, (ax1, ax2) = plt.subplots(1, 2,
                                       gridspec_kw={
                                           "width_ratios": [0.45, 0.55]},
                                       figsize=(12, 7.5))
        plt.subplots_adjust(top=0.9, bottom=0.1, wspace=0.2, hspace=0.2)

        boxplot_data = self.generate_eqtl_data(self.eqtl_effect_size,
                                               self.eqtl_direction,
                                               self.n_points)
        boxplot_data_melt = boxplot_data.melt()
        boxplot_data_melt["hue"] = boxplot_data_melt["variable"].map(
            self.color_map)

        sns.regplot(x="variable", y="value", data=boxplot_data_melt,
                    scatter=False,
                    line_kws={'color': '#000000', 'alpha': 0.5},
                    ax=ax1
                    )
        sns.boxplot(x="variable", y="value", data=boxplot_data_melt,
                    palette=self.color_map,
                    zorder=-1,
                    boxprops=dict(alpha=.3),
                    ax=ax1)
        ax1.set_xticks(range(3))
        ax1.set_xticklabels(["C/C", "C/T", "T/T"])
        ax1.set_xlabel("")
        ax1.set_ylabel("expression")
        ax1.set_yticks([])
        ax1.set_ylim(0.5, 3.5)

        interaction_data = self.add_inter_data(boxplot_data,
                                               self.interaction_effect_size,
                                               self.interaction_direction)

        ends = []
        for index in boxplot_data.columns:
            X = np.linspace(0, 1, self.n_points)
            y = interaction_data.loc[:, index].copy()
            data = pd.DataFrame({"x": X,
                                 "y": y})
            data["hue"] = self.color_map[index]

            slope, intercept, _, _, _ = st.linregress(X, y)
            ends.append(intercept + slope)

            sns.regplot(x="x", y="y", data=data,
                        scatter_kws={'facecolors': data['hue'],
                                     'edgecolors': data['hue'],
                                     'alpha': 0.75},
                        line_kws={"color": self.color_map[index],
                                  "alpha": 0.75},
                        ax=ax2

                        )

        y_min = interaction_data.values.min() - 0.25
        y_max = interaction_data.values.max() + 0.25

        ax2.set_xticks([])
        ax2.set_xlabel("context")
        ax2.set_ylabel("")
        ax2.set_yticks([])
        ax2.set_ylim(y_min, y_max)

        for ax in [ax1, ax2]:
            self.remove_axis(ax)

        y_pos = boxplot_data.values.max()
        print(y_pos)

        ax1.annotate('EXPRESSION INCREASES\nBASED ON GENOTYPE', xy=(1.5, y_pos),
                     arrowprops=dict(arrowstyle='->'),
                     xytext=(0, y_pos+0.5))

        x_pos = 1.02
        a = 0.02
        b = 0.05
        lines = [[(x_pos, ends[0] + b), (x_pos, ends[1] - b)],
                 [(x_pos - a, ends[0] + b), (x_pos + a, ends[0] + b)],
                 [(x_pos - a, ends[1] - b), (x_pos + a, ends[1] - b)],
                 [(x_pos, ends[1] + b), (x_pos, ends[2] - b)],
                 [(x_pos - a, ends[1] + b), (x_pos + a, ends[1] + b)],
                 [(x_pos - a, ends[2] - b), (x_pos + a, ends[2] - b)],
                 ]
        for (x1, y1), (x2, y2) in lines:
            ax2.plot((x1, x2), (y1, y2), c='black')

        ax2.annotate('SMALL EQTL\nEFFECT', xy=(0.1, 0.1),
                     arrowprops=dict(arrowstyle='->'),
                     xytext=(0.1, 0.6))
        ax2.annotate('LARGE EQTL\nEFFECT', xy=(0.95, ends[2]-0.05),
                     arrowprops=dict(arrowstyle='->'),
                     xytext=(0.5, ends[2]))
        ax2.annotate('EFFECT SIZE IS\nINFLUENCED BY',
                     xy=(x_pos + 0.01, (ends[1] + ends[2]) / 2),
                     arrowprops=dict(arrowstyle='->'),
                     xytext=(x_pos + 0.1, ends[1] + b))
        ax2.annotate('CONTEXT', xy=(x_pos + 0.01, (ends[0] + ends[1]) / 2),
                     arrowprops=dict(arrowstyle='->'),
                     xytext=(x_pos + 0.1, ends[1] - b))

        # a = 0.1
        # ax2.annotate('EQTL effect\ndepends on',
        #              xy=(1-(2*a), ends[2]+(1*a)),
        #              arrowprops=dict(arrowstyle='->'),
        #              xytext=(0.1, y_max*0.78))
        # ax2.annotate('context',
        #              xy=(0+a, 0+(2*a)),
        #              arrowprops=dict(arrowstyle='->'),
        #              xytext=(0.1, y_max*0.72))

        plt.show()

    @ staticmethod
    def generate_eqtl_data(effect_size, direction, n_points, mu=1, sigma=0.3):
        order = None
        if direction == "positive":
            order = [1, 1.5, 2]
        elif direction == "negative":
            order = [2, 1.5, 1]
        else:
            print("unexpected direction")
            exit()
        order = [x + random.uniform(-0.25, 0.25) for x in order]

        df = pd.DataFrame({0.0: np.random.normal(mu * effect_size * order[0],
                                                 random.uniform(0, sigma),
                                                 n_points),
                           1.0: np.random.normal(mu * effect_size * order[1],
                                                 random.uniform(0, sigma),
                                                 n_points),
                           2.0: np.random.normal(mu * effect_size * order[2],
                                                 random.uniform(0, sigma),
                                                 n_points)})
        print(df)

        return df

    @staticmethod
    def add_inter_data(df, effect_size, direction):
        n_points = len(df.index)
        start = None
        stop = None
        if direction == "positive":
            start = 0
            stop = effect_size
        elif direction == "negative":
            start = effect_size
            stop = 0
        else:
            print("unexpected direction")
            exit()

        means = df.mean(axis=0)
        order = [means.idxmin(), 1.0, means.idxmax()]

        data_new = []
        for col, mul in zip(order, [-1, 0, 1]):
            y_data = df.loc[:, col]
            y_data_demean = (y_data - y_data.mean()) + 1
            y_data_demean.sort_values(inplace=True)

            inter_baseline = np.linspace(start, stop, n_points) * (
                        mul + random.uniform(-0.2, 0.2))
            inter_jitter = np.random.normal(0, 0.1, n_points)
            inter_values = inter_baseline + inter_jitter

            data_new.append(
                [x * y for x, y in zip(y_data_demean, inter_values)])

        df = pd.DataFrame(data_new, index=order).T
        print(df)

        return df

    @staticmethod
    def remove_axis(ax, sides=None):
        if sides is None:
            sides = ['right', 'top']
        if sides is "all":
            sides = ['right', 'left', 'bottom', 'top']
        for side in sides:
            ax.spines[side].set_color('none')


if __name__ == "__main__":
    m = main()
    m.start()
