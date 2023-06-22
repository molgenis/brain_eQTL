#!/usr/bin/env python3

"""
File:         remove_covariate_effects.py
Created:      2020/10/19
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

# Local application imports.
from remove_covariate_effects.src.main import Main
from remove_covariate_effects.src.cmd_line_arguments import CommandLineArguments

# Metadata
__program__ = "Remove Covariate Effects"
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
    MATRIX = CLA.get_argument("matrix")
    COVARIATES = CLA.get_argument("covariates")
    SAMPLE_DICT = CLA.get_argument("sample_dict")

    # Start the program.
    PROGRAM = Main(matrix=MATRIX,
                   covariates=COVARIATES,
                   sample_dict=SAMPLE_DICT)
    PROGRAM.start()