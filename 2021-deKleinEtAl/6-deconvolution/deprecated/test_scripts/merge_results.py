#!/usr/bin/env python3

"""
File:         merge_results.py
Created:      2020/07/09
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
import json
import glob
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Merge Results"
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
        self.input_directory = "/home/mvochteloo/Documents/MetaBrain/Deconvolution_tests/"

    def start(self):
        data = []

        for fpath in glob.glob(os.path.join(self.input_directory, "*")):
            split_path = os.path.basename(fpath).split("_")
            cohort = split_path[0]
            if "SINGLE_CELL" in fpath or "CELLMAP" in fpath:
                input_file = "_".join(split_path[1:-3])
                ground_truth = "_".join(split_path[-3:-1])
                method = split_path[-1]
            elif "IHC" in fpath:
                input_file = "_".join(split_path[1:-2])
                ground_truth = split_path[-2]
                method = split_path[-1]
            else:
                continue

            settings = json.loads(open(os.path.join(fpath, "settings.json")).read())
            data.append([input_file, ground_truth, method, settings["avg_residuals"], settings["comparison_rss"]])

        df = pd.DataFrame(data, columns=["Input", "GroundTruth", "Method", "Avg.Residuals", "RSS"])
        df.sort_values(by="RSS", ascending=True, inplace=True)
        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None):
            print(df)
        df.round(2).to_csv(os.path.join(self.input_directory, "partial_deconvolution_comparison_results.csv"), sep="\t", header=True, index=False)

        for method in df["Method"].unique():
            subset = df.loc[df["Method"] == method, :].copy()
            print(method)
            print(subset.groupby(["Input", "GroundTruth"]).mean())

        for input in df["Input"].unique():
            subset = df.loc[df["Input"] == input, :].copy()
            print(input)
            print(subset.groupby(["GroundTruth", "Method"]).mean())

        for gt in df["GroundTruth"].unique():
            subset = df.loc[df["GroundTruth"] == gt, :].copy()
            print(gt)
            print(subset.groupby(["Input", "Method"]).mean())



if __name__ == "__main__":
    m = main()
    m.start()
