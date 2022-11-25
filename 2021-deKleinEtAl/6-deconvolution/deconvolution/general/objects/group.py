"""
File:         group.py
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


class Group:
    def __init__(self, id, samples):
        self.id = "group_{}".format(id)
        self.samples = samples

        self.snp_indices = np.array([], dtype=np.int16)
        self.sample_indices = None
        self.eqtls = []

    def add_eqtl(self, eqtl):
        self.eqtls.append(eqtl)
        self.add_snp_index(eqtl.get_snp_index())
        self.set_sample_indices(eqtl.get_sample_indices())

    def add_snp_index(self, index):
        self.snp_indices = np.append(self.snp_indices, index)

    def set_sample_indices(self, indices):
        if self.sample_indices is None:
            self.sample_indices = indices

    def get_id(self):
        return self.id

    def get_samples(self):
        return self.samples

    def get_snp_indices(self):
        return self.snp_indices

    def get_sample_indices(self):
        return self.sample_indices

    def get_eqtls(self):
        return self.eqtls

    def get_n_eqtls(self):
        return len(self.eqtls)

    def get_n_samples(self):
        return len(self.sample_indices)

    def matches(self, sample_indices):
        return np.array_equal(self.sample_indices, sample_indices)
