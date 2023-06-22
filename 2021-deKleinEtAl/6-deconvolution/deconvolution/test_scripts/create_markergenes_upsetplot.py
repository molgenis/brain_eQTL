#!/usr/bin/env python3

"""
File:         create_markergenes_upsetplot.py
Created:      2021/01/13
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
import itertools
import os

# Third party imports.
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import upsetplot as up
import matplotlib.pyplot as plt

# Local application imports.


# Metadata
__program__ = "Create Marker Genes UpsetPlot"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
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
        self.data_path = "/Users/mvochteloo/Documents/PhD/Data/PyshENCODE/Derived/SC_Decomp/DER-19_Single_cell_markergenes_TPM.xlsx"

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'create_markergenes_upsetplot')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.trans_dict = {
            'Adult-Ex1': 'Adult-Ex',
            'Adult-Ex2': 'Adult-Ex',
            'Adult-Ex3': 'Adult-Ex',
            'Adult-Ex4': 'Adult-Ex',
            'Adult-Ex5': 'Adult-Ex',
            'Adult-Ex6': 'Adult-Ex',
            'Adult-Ex7': 'Adult-Ex',
            'Adult-Ex8': 'Adult-Ex',
            'Adult-In1': 'Adult-In',
            'Adult-In2': 'Adult-In',
            'Adult-In3': 'Adult-In',
            'Adult-In4': 'Adult-In',
            'Adult-In5': 'Adult-In',
            'Adult-In6': 'Adult-In',
            'Adult-In7': 'Adult-In',
            'Adult-In8': 'Adult-In',
            'Adult-Micro': 'Adult-Micro',
            'Adult-OPC': 'Adult-OPC',
            'Adult-Endo': 'Adult-Endo',
            'Adult-Astro': 'Adult-Astro',
            'Adult-Oligo': 'Adult-Oligo',
            'Adult-OtherNeuron': 'Adult-OtherNeuron',
            'Dev-replicating': 'Dev-replicating',
            'Dev-quiescent': 'Dev-quiescent'
        }

    def start(self):
        print("Loading data.")
        df = pd.read_excel(self.data_path, header=0, index_col=None, sheet_name="Sheet2")
        print(df["CellType"].unique())
        df["group"] = df["CellType"].map(self.trans_dict)
        print(df)

        print("Preprocessing data.")
        data = {}
        for ct in df["group"].unique():
            data[ct] = set(df.loc[df["group"] == ct, "GeneName"].tolist())
        print(data)
        counts = self.count(data)
        counts = counts[counts > 0]
        print(counts)

        print("Creating plot.")
        up.plot(counts, sort_by='cardinality', show_counts=True)
        plt.savefig(os.path.join(self.outdir, "markergenes_upsetplot.png"))
        plt.close()

    @staticmethod
    def count(input_data):
        combinations = []
        cols = list(input_data.keys())
        for i in range(1, len(cols) + 1):
            combinations.extend(list(itertools.combinations(cols, i)))

        abbreviations = {"CellMapNNLS_Neuron": "neuro",
                         "CellMapNNLS_Oligodendrocyte": "oligo",
                         "CellMapNNLS_EndothelialCell": "endo",
                         "CellMapNNLS_Macrophage": "macro",
                         "CellMapNNLS_Astrocyte": "astro",
                         "CellMapNNLS_Inhibitory": "inhib",
                         "CellMapNNLS_Microglia": "micro",
                         "CellMapNNLS_Excitatory": "excit",
                         "CellMapNNLS_Pericytes": "peri"}
        abbr_cols = []
        for col in cols:
            if col in abbreviations.keys():
                abbr_cols.append(abbreviations[col])
            else:
                abbr_cols.append(col)

        indices = []
        data = []
        for combination in combinations:
            index = []
            for col in cols:
                if col in combination:
                    index.append(True)
                else:
                    index.append(False)

            background = set()
            for key in cols:
                if key not in combination:
                    work_set = input_data[key].copy()
                    background.update(work_set)

            overlap = None
            for key in combination:
                work_set = input_data[key].copy()
                if overlap is None:
                    overlap = work_set
                else:
                    overlap = overlap.intersection(work_set)

            duplicate_set = overlap.intersection(background)
            length = len(overlap) - len(duplicate_set)

            indices.append(index)
            data.append(length)

        s = pd.Series(data,
                      index=pd.MultiIndex.from_tuples(indices, names=abbr_cols))
        s.name = "value"
        return s


if __name__ == '__main__':
    m = main()
    m.start()
