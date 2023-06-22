#!/usr/bin/env python3

"""
File:         custom_interaction_analyser.py
Created:      2020/10/15
Last Changed: 2020/10/27
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
from custom_interaction_analyser.src.main import Main
from custom_interaction_analyser.src.combiner import Combine
from custom_interaction_analyser.src.cmd_line_arguments import \
    CommandLineArguments

# Metadata
__program__ = "Custom Interaction Analyser"
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
    INPUT = CLA.get_argument("input")
    OUTPUT = CLA.get_argument("output")
    SETTINGS_FILE = CLA.get_argument("settings")
    COMBINE = CLA.get_argument("combine")

    if OUTPUT is None:
        OUTPUT = INPUT

    if COMBINE:
        ALPHA = CLA.get_argument("alpha")
        FORCE = CLA.get_argument("force")

        # Start the program.
        PROGRAM = Combine(input_folder=OUTPUT,
                          settings_file=SETTINGS_FILE,
                          alpha=ALPHA,
                          force=FORCE)
        PROGRAM.start()
    else:
        SKIP_ROWS = CLA.get_argument("skip_rows")
        N_EQTLS = CLA.get_argument("n_eqtls")
        N_SAMPLES = CLA.get_argument("n_samples")
        VERBOSE = CLA.get_argument("verbose")

        # Start the program.
        PROGRAM = Main(input_folder=INPUT,
                       output_folder=OUTPUT,
                       settings_file=SETTINGS_FILE,
                       skip_rows=SKIP_ROWS,
                       n_eqtls=N_EQTLS,
                       n_samples=N_SAMPLES,
                       verbose=VERBOSE)
        PROGRAM.start()
