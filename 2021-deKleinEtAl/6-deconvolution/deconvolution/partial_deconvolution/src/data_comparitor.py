"""
File:         data_comparitor.py
Created:      2020/06/30
Last Changed: 2020/11/22
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

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.


class DataComparitor:
    def __init__(self, settings, deconvolution, ground_truth):
        self.deconvolution = deconvolution
        self.ground_truth = ground_truth

        self.comparison = None
        self.n_samples = None
        self.rss = None

    def work(self):
        if self.ground_truth is None:
            return

        decon_df = self.deconvolution.copy()
        truth_df = self.ground_truth.copy()

        overlap = np.intersect1d(decon_df.index, truth_df.index)
        decon_df = decon_df.loc[overlap, :]
        truth_df = truth_df.loc[overlap, :]

        decon_df.index.name = "-"
        truth_df.index.name = "-"

        if not decon_df.index.equals(truth_df.index):
            print("Invalid order")
            exit()

        print(decon_df)
        print(truth_df)

        # Sum the Ex and Ih cell types.
        trans_dict = {
            "Adult-Astrocytes": "Astrocyte",
            "Adult-Endothelial": "EndothelialCell",
            "Adult-Ex": "Neuron",
            "Adult-In": "Neuron",
            "Adult-Microglia": "Macrophage",
            "Adult-OPC": "Oligodendrocyte",
            "Adult-Oligo": "Oligodendrocyte",
            "Adult-OtherNeuron": "Neuron",
            "Dev-Quiescent": "Quiescent",
            "Dev-Replicating": "Quiescent",
            # "Adult-Astrocytes": "Astrocyte",
            # "Adult-Endothelial": "EndothelialCell",
            # "Adult-Ex": "Excitatory",
            # "Adult-In": "Inhibitory",
            # "Adult-Microglia": "Microglia",
            # "Adult-OPC": "OPC",
            # "Adult-Oligo": "Oligodendrocyte",
            # "Adult-OtherNeuron": "Neuron",
            # "Dev-Quiescent": "Quiescent",
            # "Dev-Replicating": "Quiescent",
        }
        decon_df = decon_df.T
        decon_df["cell type"] = [trans_dict[y] for y in [''.join([i for i in x if not i.isdigit()]) for x in decon_df.index]]
        decon_df = decon_df.groupby("cell type").sum()
        decon_df = decon_df.T
        decon_df.columns.name = None

        #
        # # cell_types = np.intersect1d(decon_df.columns, truth_df.columns)
        # cell_type_matches = (("Microglia", "Macrophage"),
        #                      ("EndothelialCell", "EndothelialCell"),
        #                      ("Oligodendrocyte", "Oligodendrocyte"),
        #                      ("Excitatory", "Neuron"),
        #                      ("Inhibitory", "Neuron"),
        #                      (("Excitatory", "Inhibitory"), "Neuron"),
        #                      ("Astrocyte", "Astrocyte"))
        # data = []
        # for i, sample in enumerate(overlap):
        #     decon = decon_df.iloc[i, :]
        #     truth = truth_df.iloc[i, :]
        #     #for cell in cell_types:
        #     #    data.append([sample, decon[cell], truth[cell], cell])
        #     for ct1, ct2 in cell_type_matches:
        #         if isinstance(ct1, str):
        #             value1 = decon[ct1]
        #         else:
        #             value1 = 0
        #             for sub_ct in ct1:
        #                 value1 += decon[sub_ct]
        #
        #         if isinstance(ct2, str):
        #             value2 = truth[ct2]
        #         else:
        #             value2 = 0
        #             for sub_ct in ct2:
        #                 value2 += truth[sub_ct]
        #
        #         name1 = ct1
        #         name2 = ct2
        #         if not isinstance(ct1, str):
        #             name1 = "+".join([sub_ct for sub_ct in ct1])
        #         if not isinstance(ct2, str):
        #             name2 = "+".join([sub_ct for sub_ct in ct1])
        #
        #         cell = "{}/{}".format(name1, name2)
        #         if ct1 == ct2:
        #             cell = ct1
        #
        #         data.append([sample, value1, value2, cell])
        # df = pd.DataFrame(data, columns=['-', 'x', 'y', 'hue'])
        # df.set_index('-', inplace=True)

        cell_types = np.intersect1d(decon_df.columns, truth_df.columns)
        data = []
        for i, sample in enumerate(overlap):
            decon = decon_df.iloc[i, :]
            truth = truth_df.iloc[i, :]
            for cell in cell_types:
                data.append([sample, decon[cell], truth[cell], cell])
        df = pd.DataFrame(data, columns=['-', 'x', 'y', 'hue'])
        df.set_index('-', inplace=True)

        self.comparison = df
        self.n_samples = len(overlap)
        self.rss = sum((df['x'] - df['y'])**2)

    def get_comparison(self):
        return self.comparison

    def get_n_samples(self):
        return self.n_samples

    def get_rss(self):
        return self.rss

    def print_info(self):
        if self.comparison is not None:
            print("Comparison matrix:\t{} rows and {} columns".format(self.comparison.shape[0], self.comparison.shape[1]))
        if self.n_samples is not None:
            print("N samples:\t{}".format(self.n_samples))
        if self.rss is not None:
            print("Residuals sum of squares:\t{:.2f}".format(self.rss))
        print("")
