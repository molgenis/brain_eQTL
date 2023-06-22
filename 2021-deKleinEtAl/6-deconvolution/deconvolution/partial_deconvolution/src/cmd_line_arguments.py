"""
File:         cmd_line_arguments.py
Created:      2020/06/29
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
        self.print()

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
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the bulk expression input file.")
        parser.add_argument("-s",
                            "--signature",
                            type=str,
                            required=True,
                            help="The path to the signature matrix file.")
        parser.add_argument("-t",
                            "--translate",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the gene to ensembl ID translate "
                                 "file.")
        parser.add_argument("-g",
                            "--ground_truth",
                            type=str,
                            default=None,
                            help="A matrix with the ground truth values. Must "
                                 "be in <type>_counts.<extension> format."
                                 "Default: None'.")
        parser.add_argument("-o",
                            "--outdir",
                            type=str,
                            default="output",
                            help="The name of the output directory. " \
                                 "Default: 'output'.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=False,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-sf",
                            "--sample_filter",
                            type=str,
                            default=None,
                            help="A path to a file containing which first line "
                                 "(tab seperated) contains the sample "
                                 "identifiers to include.")
        parser.add_argument("-df",
                            "--dataset_filter",
                            type=str,
                            nargs="+",
                            required=False,
                            default=None,
                            help="The name of the dataset(s) to include, "
                                 "default: include all. (space = _)")
        parser.add_argument("-os",
                            "--outsubdir",
                            type=str,
                            default=None,
                            help="The name of the output subdirectory. " \
                                 "Default: None.")
        parser.add_argument("-m",
                            "--min_expr",
                            type=int,
                            default=0,
                            help="The minimal expression value per gene."
                                 " Default: 0.")
        parser.add_argument("-dataset_correction",
                            action='store_true',
                            help="Correct for dataset means in the expression."
                                 " Default: False.")
        parser.add_argument("-n",
                            "--normalize",
                            type=str,
                            choices=["rows", "columns", "both"],
                            default=False,
                            help="Divide each value by row row / column sum."
                                 " Default: False.")
        parser.add_argument("-zscore",
                            action='store_true',
                            help="Z-score transform the profile values."
                                 " Default: False.")
        parser.add_argument("-log2",
                            action='store_true',
                            help="Log2 transform the profile values."
                                 " Default: False.")
        parser.add_argument("-dm",
                            "--decon_method",
                            type=str,
                            choices=["NNLS"],
                            default="NNLS",
                            help="The deconvolution method to use. "
                                 "Default: 'NNLS'.")
        parser.add_argument("-visualise",
                            action='store_true',
                            help="Whether or not to visualise the data."
                                 " Default: False.")
        parser.add_argument("-pid",
                            "--plot_id",
                            type=str,
                            nargs="+",
                            required=False,
                            help="The column names to plot in the -sa / "
                                 "--sample file.")
        parser.add_argument("-e",
                            "--extension",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default="png",
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser

    def print(self):
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
