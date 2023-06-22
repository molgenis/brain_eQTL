#!/usr/bin/env python3

"""
File:         artificial_bulk_deconvolution.py
Created:      2020/11/09
Last Changed: 2020/11/12
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
from src.main import Main
from src.cmd_line_arguments import CommandLineArguments

# Metadata
__program__ = "Artificial Bulk Deconvolution"
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
if __name__ == '__main__':
    # Get the command line arguments.
    CLA = CommandLineArguments(program=__program__,
                               version=__version__,
                               description=__description__)
    INPUT_DIR_PATH = CLA.get_argument("input_dir_path")
    INPUT_SUFFIX = CLA.get_argument("input_suffix")
    CELL_COUNTS_PATH = CLA.get_argument("cell_counts_path")
    REFERENCE_PROFILE = CLA.get_argument("reference_profile_path")
    GENE_INFO_PATH = CLA.get_argument("gene_info_path")
    COMBINE = CLA.get_argument("combine")
    N_SAMPLES = CLA.get_argument("n_samples")

    # Start the program.
    PROGRAM = Main(input_path=INPUT_DIR_PATH,
                   input_suffix=INPUT_SUFFIX,
                   cell_counts_path=CELL_COUNTS_PATH,
                   ref_profile_path=REFERENCE_PROFILE,
                   gene_info_path=GENE_INFO_PATH,
                   combine=COMBINE,
                   n_samples=N_SAMPLES)
    PROGRAM.start()