#!/usr/bin/env python3

"""
File:         diverging_cmap.py
Created:      2020/06/19
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

# Third party imports.
import seaborn as sns
import numpy as np

# Local application imports.

# Metadata
__program__ = "Divering CMAP"
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
        self.correct_low = (0, 114/255, 178/255)
        self.correct_high = (185/255, 39/255, 50/255)

    def start(self):
        low_end_diff = np.inf
        high_end_diff = np.inf
        best_low = None
        best_high = None
        for i in range(359):
            pallete = sns.diverging_palette(i, 24, n=9)
            low_end = pallete[0]
            high_end = pallete[-1]
            low_end_difference = sum([abs(self.correct_low[i] - low_end[i]) for i in range(3)])
            high_end_difference = sum([abs(self.correct_high[i] - high_end[i]) for i in range(3)])
            if low_end_difference < low_end_diff:
                low_end_diff = low_end_difference
                best_low = i
            if high_end_difference < high_end_diff:
                high_end_diff = high_end_difference
                best_high = i
        print(best_low, best_high)


if __name__ == "__main__":
    m = main()
    m.start()
