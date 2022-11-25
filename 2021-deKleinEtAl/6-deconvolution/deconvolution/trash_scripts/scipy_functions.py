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
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Scipy Function Visualized"
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

    def start(self):
        # Get the f-values
        rss = list(np.arange(0, 1, 1e-3)) + \
              list(np.arange(1, 5, 1e-2)) + \
              list(np.arange(5, 15, 1e-1)) + \
              list(np.arange(15, 1812, 1)) + \
              list(np.arange(1812, 1822, 1e-1)) + \
              list(np.arange(1822, 1827, 1e-2)) + \
              list(np.arange(1827, 1827.95, 1e-3))
        rss_tmp = list(np.arange(0, 1827.95, 10))
        df1 = 40
        df2 = 41
        n = 3703
        print("Total number of values: {}.".format(len(rss)))

        delta_rss = []
        f_values = []
        unique_f_values = []
        for i, val1 in enumerate(rss_tmp):
            if i % 100 == 0:
                print("\t{}/{}".format(i, len(rss_tmp)))
            for val2 in rss_tmp:
                f_value = self.get_f_value(val1, val2, df1, df2, n)
                if isinstance(f_value, float):
                    delta_rss.append(val1 - val2)
                    f_values.append(f_value)
                    if f_value not in unique_f_values:
                        unique_f_values.append(f_value)
        print("Total number of unique F-values: {}.".format(len(unique_f_values)))

        # Visualize.
        print("Visualise F-value Distribution")
        self.plot_distribution(outdir=self.outdir, x=delta_rss, y=f_values,
                               xlab="RSS1 - RSS2", ylab="F-value",
                               title="RSS to F-value Distribution",
                               filename="rss_to_fval_dist")

        exit()

        # Get the p-values.
        p_values = []
        f_values = rss
        for i, fval in enumerate(f_values):
            if i % 100 == 0:
                print("\t{}/{}".format(i, len(f_values)))
            p_values.append(self.get_p_value(fval, df1, df2, n))

        # Visualize.
        print("Visualise F-value to P-value Distribution")
        self.plot_distribution(outdir=self.outdir, x=f_values, y=p_values,
                               xlab="F-value", ylab="P-value", xlim=(-1, 5),
                               title="F-value to P-value Distribution",
                               filename="fval_to_pval_dist")

        # Get the z-scores.
        z_scores = []
        p_values = [1e-323, 1e-322, 1e-321, 1e-320, 1e-319, 1e-318, 1e-317, 1e-316, 1e-315, 1e-314, 1e-313, 1e-312, 1e-311, 1e-310, 1e-309, 1e-308, 1e-307, 1e-306, 1e-305, 1e-304, 1e-303, 1e-302, 1e-301, 1e-300, 1e-299, 1e-298, 1e-297, 1e-296, 1e-295, 1e-294, 1e-293, 1e-292, 1e-291, 1e-290, 1e-289, 1e-288, 1e-287, 1e-286, 1e-285, 1e-284, 1e-283, 1e-282, 1e-281, 1e-280, 1e-279, 1e-278, 1e-277, 1e-276, 1e-275, 1e-274, 1e-273, 1e-272, 1e-271, 1e-270, 1e-269, 1e-268, 1e-267, 1e-266, 1e-265, 1e-264, 1e-263, 1e-262, 1e-261, 1e-260, 1e-259, 1e-258, 1e-257, 1e-256, 1e-255, 1e-254, 1e-253, 1e-252, 1e-251, 1e-250, 1e-249, 1e-248, 1e-247, 1e-246, 1e-245, 1e-244, 1e-243, 1e-242, 1e-241, 1e-240, 1e-239, 1e-238, 1e-237, 1e-236, 1e-235, 1e-234, 1e-233, 1e-232, 1e-231, 1e-230, 1e-229, 1e-228, 1e-227, 1e-226, 1e-225, 1e-224, 1e-223, 1e-222, 1e-221, 1e-220, 1e-219, 1e-218, 1e-217, 1e-216, 1e-215, 1e-214, 1e-213, 1e-212, 1e-211, 1e-210, 1e-209, 1e-208, 1e-207, 1e-206, 1e-205, 1e-204, 1e-203, 1e-202, 1e-201, 1e-200, 1e-199, 1e-198, 1e-197, 1e-196, 1e-195, 1e-194, 1e-193, 1e-192, 1e-191, 1e-190, 1e-189, 1e-188, 1e-187, 1e-186, 1e-185, 1e-184, 1e-183, 1e-182, 1e-181, 1e-180, 1e-179, 1e-178, 1e-177, 1e-176, 1e-175, 1e-174, 1e-173, 1e-172, 1e-171, 1e-170, 1e-169, 1e-168, 1e-167, 1e-166, 1e-165, 1e-164, 1e-163, 1e-162, 1e-161, 1e-160, 1e-159, 1e-158, 1e-157, 1e-156, 1e-155, 1e-154, 1e-153, 1e-152, 1e-151, 1e-150, 1e-149, 1e-148, 1e-147, 1e-146, 1e-145, 1e-144, 1e-143, 1e-142, 1e-141, 1e-140, 1e-139, 1e-138, 1e-137, 1e-136, 1e-135, 1e-134, 1e-133, 1e-132, 1e-131, 1e-130, 1e-129, 1e-128, 1e-127, 1e-126, 1e-125, 1e-124, 1e-123, 1e-122, 1e-121, 1e-120, 1e-119, 1e-118, 1e-117, 1e-116, 1e-115, 1e-114, 1e-113, 1e-112, 1e-111, 1e-110, 1e-109, 1e-108, 1e-107, 1e-106, 1e-105, 1e-104, 1e-103, 1e-102, 1e-101, 1e-100, 1e-99, 1e-98, 1e-97, 1e-96, 1e-95, 1e-94, 1e-93, 1e-92, 1e-91, 1e-90, 1e-89, 1e-88, 1e-87, 1e-86, 1e-85, 1e-84, 1e-83, 1e-82, 1e-81, 1e-80, 1e-79, 1e-78, 1e-77, 1e-76, 1e-75, 1e-74, 1e-73, 1e-72, 1e-71, 1e-70, 1e-69, 1e-68, 1e-67, 1e-66, 1e-65, 1e-64, 1e-63, 1e-62, 1e-61, 1e-60, 1e-59, 1e-58, 1e-57, 1e-56, 1e-55, 1e-54, 1e-53, 1e-52, 1e-51, 1e-50, 1e-49, 1e-48, 1e-47, 1e-46, 1e-45, 1e-44, 1e-43, 1e-42, 1e-41, 1e-40, 1e-39, 1e-38, 1e-37, 1e-36, 1e-35, 1e-34, 1e-33, 1e-32, 1e-31, 1e-30, 1e-29, 1e-28, 1e-27, 1e-26, 1e-25, 1e-24, 1e-23, 1e-22, 1e-21, 1e-20, 1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3] + list(np.arange(1e-3, (1-1e-2), 1e-2)) + [(1-1e-2), (1-1e-3), (1-1e-4), (1-1e-5), (1-1e-6), (1-1e-7), (1-1e-8), (1-1e-9), (1-1e-10), (1-1e-11), (1-1e-12), (1-1e-13), (1-1e-14), (1-1e-15), (1-1e-16)]
        for i, pval in enumerate(p_values):
            if i % 100 == 0:
                print("\t{}/{}".format(i, len(p_values)))
            new_z_score = self.get_z_score(pval)
            if isinstance(new_z_score, float):
                z_scores.append(new_z_score)
            else:
                break

        # Visualize.
        print("Visualise P-value to Z-score Distribution")
        self.plot_distribution(outdir=self.outdir, x=p_values, y=z_scores,
                               xlab="P-value", ylab="Z-score",
                               title="P-value to Z-score Distribution",
                               filename="pval_to_zscore_dist")

    @staticmethod
    def get_f_value(rss1, rss2, df1, df2, n):
        if df1 >= df2:
            return np.nan
        if df2 >= n:
            return np.nan
        if (rss2 / (n - df2)) <= 0:
            return np.nan

        return ((rss1 - rss2) / (df2 - df1)) / (rss2 / (n - df2))

    @staticmethod
    def get_p_value(f_value, df1, df2, n):
        # stats.f.sf(1827.95, dfn=1, dfd=3661) = 5e-324
        # stats.f.sf(9.9e-12, dfn=1, dfd=3661) = 0.9999974567714613

        # stats.f.cdf(69, dfn=1, dfd=3661) = 0.9999999999999999
        # stats.f.cdf(1e-320, dfn=1, dfd=3661) = 1.0730071046473278e-160

        return stats.f.sf(f_value, dfn=(df2 - df1), dfd=(n - df2))

    @staticmethod
    def get_z_score(p_value):
        if p_value > (1.0 - 1e-16):
            p_value = 1e-16
        if p_value < 1e-323:
            p_value = 1e-323

        # stats.norm.isf((1 - 1e-16)) = -8.209536151601387
        # stats.norm.isf(1e-323) = 38.44939448087599
        return stats.norm.isf(p_value)

    @staticmethod
    def plot_distribution(outdir, x, y, xlim=None, ylim=None,
                          xlab="", ylab="", title="", filename=""):
        df = pd.DataFrame({"x": x,
                           "y": y})
        df.dropna(inplace=True)

        fig, ax = plt.subplots()
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        g = sns.scatterplot(x=x,
                            y=y,
                            data=df,
                            facecolors='#000000',
                            edgecolor='#000000',
                            alpha=0.5)
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)

        g.set_title(title)
        g.set_ylabel(ylab,
                     fontsize=8,
                     fontweight='bold')
        g.set_xlabel(xlab,
                     fontsize=8,
                     fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, filename + ".png"))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
