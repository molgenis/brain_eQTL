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
        self.n_points = 200
        self.eqtl_effect_size = 1
        self.eqtl_direction = "positive"
        self.interaction_effect_size = 0.5
        self.interaction_direction = "positive"
        self.color_map = {0.0: "#ABDCA2", 1.0: "#89A3D1", 2.0: "#F08C72"}

    def start(self):
        plt.xkcd()

        fig, ax = plt.subplots(figsize=(12, 7.5))
        plt.subplots_adjust(top=0.9, bottom=0.1, wspace=0.2, hspace=0.2)
        self.remove_axis(ax)

        boxplot_data = self.generate_eqtl_data(self.eqtl_effect_size,
                                               self.eqtl_direction,
                                               self.n_points)

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
                        ax=ax

                        )

        y_min = interaction_data.values.min() - 0.25
        y_max = interaction_data.values.max() + 1

        ax.set_xticks([])
        ax.set_xlabel("context")
        ax.set_ylabel("expression")
        ax.set_yticks([])
        ax.set_xlim(-0.1, 1.75)
        ax.set_ylim(y_min, y_max)

        # this is an inset axes over the main axes
        ax1 = plt.axes([0.1, 0.65, .25, .25])
        self.remove_axis(ax1)
        self.plot_boxplot(ax1, boxplot_data, self.color_map, demean=True)
        ax2 = plt.axes([0.7, 0.65, .25, .25])
        self.remove_axis(ax2)
        self.plot_boxplot(ax2, boxplot_data, self.color_map)

        print(y_max)
        ax.annotate('',
                    xy=(0.1, 0.25),
                    arrowprops=dict(arrowstyle='->'),
                    xytext=(0.25, ends[2] * 0.65))
        ax.annotate('',
                    xy=(1.05, 0.25),
                    arrowprops=dict(arrowstyle='->'),
                    xytext=(1.4, ends[2] * 0.65))

        plt.show()

    @ staticmethod
    def generate_eqtl_data(effect_size, direction, n_points, mu=1, sigma=0.5):
        order = None
        if direction == "positive":
            order = [1, 1.5, 2]
        elif direction == "negative":
            order = [2, 1.5, 1]
        else:
            print("unexpected direction")
            exit()
        order = [x + random.uniform(-0.5, 0.5) for x in order]

        df = pd.DataFrame({0.0: np.random.normal(mu * effect_size * order[0],
                                                 random.uniform(0.1, sigma),
                                                 n_points),
                           1.0: np.random.normal(mu * effect_size * order[1],
                                                 random.uniform(0.1, sigma),
                                                 n_points),
                           2.0: np.random.normal(mu * effect_size * order[2],
                                                 random.uniform(0.1, sigma),
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

    @staticmethod
    def plot_boxplot(ax, data, color_map, demean=False):
        df = data.copy()
        if demean:
            df = df - df.mean()

        df_melt = df.melt()
        df_melt["hue"] = df_melt["variable"].map(color_map)

        sns.regplot(x="variable", y="value", data=df_melt,
                    scatter=False,
                    line_kws={'color': '#000000', 'alpha': 0.5},
                    ax=ax
                    )
        sns.boxplot(x="variable", y="value", data=df_melt,
                    palette=color_map,
                    zorder=-1,
                    boxprops=dict(alpha=.3),
                    ax=ax)
        ax.set_xticks(range(3))
        ax.set_xticklabels(["C/C", "C/T", "T/T"])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_yticks([])

if __name__ == "__main__":
    m = main()
    m.start()
