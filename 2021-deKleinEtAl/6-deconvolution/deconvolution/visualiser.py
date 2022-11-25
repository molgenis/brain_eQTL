#!/usr/bin/env python3

"""
File:         visualiser.py
Created:      2020/03/13
Last Changed: 2020/06/05
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
from visualiser.src.main import Main
from visualiser.src.cmd_line_arguments import CommandLineArguments

# Metadata
__program__ = "Visualiser"
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
if __name__ == '__main__':
    # Get the command line arguments.
    CLA = CommandLineArguments(program=__program__,
                               version=__version__,
                               description=__description__)
    NAME = CLA.get_argument("name")
    SETTINGS_FILE = CLA.get_argument("settings")
    ALPHA = CLA.get_argument("alpha")
    PLOTS = CLA.get_argument("plots")
    TOP = CLA.get_argument("top")
    INTEREST = CLA.get_argument("interest")
    EXTENSION = CLA.get_argument("extension")
    VALIDATE = CLA.get_argument("validate")

    # Start the program.
    PROGRAM = Main(name=NAME,
                   settings_file=SETTINGS_FILE,
                   alpha=ALPHA,
                   plots=PLOTS,
                   top=TOP,
                   interest=INTEREST,
                   extension=EXTENSION,
                   validate=VALIDATE)
    PROGRAM.start()
