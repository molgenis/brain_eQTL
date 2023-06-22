#!/usr/bin/env python3

"""
File:         combine_network_output.py
Created:      2020/05/19
Last Changed: 2022/02/10
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
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Combine Network Output"
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
        self.work_dir = "/home/mvochteloo/Documents/Master/Graduation/Results/cis_output/TestScripts/network_output"
        self.outdir = os.path.join(self.work_dir, 'combined')
        self.nrows = 20

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):

        covariate_data = {}
        for fpath in glob.glob(os.path.join(self.work_dir, "*.txt")):
            filename = os.path.basename(fpath)
            (network, _, covariate, group, info) = filename.split(".")[0].split("-")

            data = pd.read_csv(fpath, sep="\t", skiprows=12, header=0)
            data = data.loc[data["p-value"] <= 0.05, :]
            gsn_values = list(data["gene_set_name"].values)

            if covariate not in covariate_data.keys():
                covariate_data[covariate] = {}

            if group in covariate_data[covariate].keys():
                (cov_data, cov_cols) = covariate_data[covariate][group]
                cov_data.append(gsn_values)
                cov_cols.append("{}{}".format(network, info.upper()))
                covariate_data[covariate][group] = (cov_data, cov_cols)
            else:
                covariate_data[covariate][group] = ([gsn_values],
                                                    ["{}{}".format(network, info.upper())])

        for key, value in covariate_data.items():
            for group, (data, cols) in value.items():
                print("COVARIATE: {}\tGROUP: {}".format(key, group))
                df = pd.DataFrame(data, index=cols).T
                df = df.reindex(sorted(df.columns), axis=1)
                subset = df.iloc[:self.nrows, :]
                with pd.option_context('display.max_rows', None,
                                       'display.max_columns', None):
                    print(subset)
                subset.to_csv(os.path.join(self.outdir, "{}_{}.txt".format(key, group)),
                              sep="\t", header=True, index=False)


if __name__ == '__main__':
    m = main()
    m.start()
