"""
File:         cmd_line_arguments.py
Created:      2020/10/08
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
import argparse

# Third party imports.

# Local application imports.


class CommandLineArguments:
    def __init__(self, program, version, description):
        # Safe variables.
        self.program = program
        self.version = version
        self.description = description

        # Get the arguments.
        parser = self.create_argument_parser()
        self.arguments = parser.parse_args()
        self.print_arguments()

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=self.program,
                                         description=self.description)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(self.program,
                                                   self.version),
                            help="show program's version number and exit")
        parser.add_argument("-n",
                            "--name",
                            type=str,
                            required=True,
                            help="The name of the input/output directory.")
        parser.add_argument("-s",
                            "--settings",
                            type=str,
                            default="default_settings",
                            help="The settings input file (without '.json'), "
                                 "default: 'default_settings'.")
        parser.add_argument("-f",
                            "--force_steps",
                            nargs="+",
                            type=str,
                            default=None,
                            choices=["all",
                                     "combine_gte_files",
                                     "combine_eqtlprobes",
                                     "create_cohort_matrix",
                                     "create_matrices",
                                     "correct_cohort_effects",
                                     "perform_deconvolution",
                                     "cerate_tech_covs_matrix",
                                     "create_covs_matrix",
                                     "create_extra_covs_matrix",
                                     "normal_transform_matrix"],
                            help="The steps to force the program to redo, "
                                 "default: None. Note, all dependend steps"
                                 "are forced too.")
        parser.add_argument("-ecm",
                            "--extra_cov_matrix",
                            type=str,
                            default=None,
                            help="The path to matrices to prepare as"
                                 "additional covariate matrix. Note: assumes"
                                 "tab separated with header and index."
                                 "Default: None.")
        parser.add_argument("-clear_log",
                            action='store_true',
                            help="Clear already existing log files. default: "
                                 "'False'.")

        return parser

    def print_arguments(self):
        for arg in vars(self.arguments):
            print("Input argument '{}' "
                  "has value '{}'.".format(arg, getattr(self.arguments, arg)))

    def get_argument(self, arg_key):
        if self.arguments is not None and arg_key in self.arguments:
            value = getattr(self.arguments, arg_key)
        else:
            value = None

        return value

    def get_all_arguments(self):
        return self.arguments
