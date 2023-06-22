"""
File:         storage_container.py
Created:      2020/10/14
Last Changed: 2020/10/21
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


# Local application imports.


class StorageContainer:
    def __init__(self, colnames):
        # Initialize the result list.
        self.pvalues_buffer = [[-1, "-"] + colnames]
        self.coefficients_buffer = [[-1, "-"] + colnames]
        self.std_errors_buffer = [[-1, "-"] + colnames]
        self.perm_pvalues_buffer = []

        # Initialize a new row.
        self.pvalues = None
        self.perm_pvalues = None
        self.coefficients = None
        self.std_errors = None

        # Initialize variables.
        self.n_rows = 0
        self.error = False

    def add_row(self, eqtl_index, genotype_name):
        self.pvalues = [eqtl_index, genotype_name]
        self.coefficients = [eqtl_index, genotype_name]
        self.std_errors = [eqtl_index, genotype_name]

        self.perm_pvalues = []

    def store_row(self):
        self.pvalues_buffer.append(self.pvalues)
        self.perm_pvalues_buffer.extend(self.perm_pvalues)
        self.coefficients_buffer.append(self.coefficients)
        self.std_errors_buffer.append(self.std_errors)

        self.pvalues = None
        self.perm_pvalues = None
        self.coefficients = None
        self.std_errors = None

        self.n_rows = self.n_rows + 1

    def add_pvalue(self, order_id, value):
        if order_id == 0:
            self.pvalues.append(value)
        else:
            self.perm_pvalues.append(value)

    def add_coefficient(self, order_id, value):
        if order_id == 0:
            self.coefficients.append(value)
        else:
            pass

    def add_std_error(self, order_id, value):
        if order_id == 0:
            self.std_errors.append(value)
        else:
            pass

    def get_pvalues(self):
        return self.pvalues_buffer

    def get_perm_pvalues(self):
        return self.perm_pvalues_buffer

    def get_coefficients(self):
        return self.coefficients_buffer

    def get_std_errors(self):
        return self.std_errors_buffer

    def has_error(self):
        return self.error

    def set_error(self):
        self.error = True

    def get_n_rows(self):
        return self.n_rows
