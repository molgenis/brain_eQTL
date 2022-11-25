"""
File:         eqtl.py
Created:      2020/03/19
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

# Third party imports.
import numpy as np

# Local application imports.


class Eqtl:
    def __init__(self, name, snp_index, genotype, expression):
        self.name = name
        self.snp_index = snp_index

        df = self.combine_data(genotype, expression)
        self.samples, self.sample_indices = self.filter_data(df)
        self.n_missing = genotype.size - len(self.samples)
        del df

        self.group = None

    @staticmethod
    def combine_data(genotype, expression):
        genotype_df = genotype.T
        expression_df = expression.T
        data = genotype_df.merge(expression_df,
                                 left_index=True,
                                 right_index=True)
        data.columns = ["genotype", "expression"]
        del genotype_df, expression_df
        return data

    @staticmethod
    def filter_data(df):
        tmp_df = df.copy()
        tmp_df.replace(-1, np.nan, inplace=True)
        indices = np.arange(tmp_df.shape[0])
        sample_indices = indices[~tmp_df.isnull().any(axis=1).values]
        tmp_df.dropna(inplace=True)
        samples = tmp_df.index.to_list()
        del tmp_df
        return samples, sample_indices

    def get_name(self):
        return self.name

    def get_snp_index(self):
        return self.snp_index

    def get_samples(self):
        return self.samples

    def get_sample_indices(self):
        return self.sample_indices

    def get_n_samples(self):
        return len(self.samples)

    def get_n_missing(self):
        return self.n_missing
