"""
File:         cmd_line_arguments.py
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
import argparse

# Third party imports.

# Local application imports.


class CommandLineArguments:
    """
    CommandLineArguments: class that processing command line arguments.
    """

    def __init__(self, program, version, description):
        """
        Initializer method.

        :param program: string, name of the program.
        :param version: float, the version of the program.
        :param description: string, description of the program.
        """
        # Safe variables.
        self.program = program
        self.version = version
        self.description = description

        # Get the arguments.
        parser = self.create_argument_parser()
        self.arguments = parser.parse_args()
        self.print()

    def create_argument_parser(self):
        """
        Function to create the command line flags.

        :return: parser object,
        """
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
        parser.add_argument("-a",
                            "--alpha",
                            type=float,
                            default=0.05,
                            help="The significance cut-off,"
                                 "default: 0.05.")
        parser.add_argument("-p",
                            "--plots",
                            nargs="+",
                            type=str,
                            default=["all"],
                            choices=["covariate_clustermap",
                                     "covariate_comparison",
                                     "covariates_explained_by_others",
                                     "deconvolution_covariate_comparison",
                                     "deconvolution_zscore_comparison",
                                     "simple_eqtl_effect",
                                     "inter_clustermap",
                                     "inter_zscore_bars",
                                     "inter_pvalue_boxplot",
                                     "inter_zscore_dist",
                                     "inter_zscore_marker_genes",
                                     "inter_eqtl_zscore_bars",
                                     "inter_eqtl_effect",
                                     "inter_eqtl_effect_deconvolution",
                                     "inter_eqtl_effect_marker_vs_comp",
                                     "inter_eqtl_effect_celltype",
                                     "inter_eqtl_celltype_details"
                                     ],
                            help="The name of the figures to be created, "
                                 "default: 'all'.")
        parser.add_argument("-t",
                            "--top",
                            type=int,
                            default=1,
                            help="The number of top eQTLs to visualise, "
                                 "default: 1.")
        parser.add_argument("-i",
                            "--interest",
                            nargs="+",
                            type=int,
                            default=None,
                            help="The indices of the eQTLS to visualise,"
                                 "default: None. If set, -t / --top is "
                                 "discarded.")
        parser.add_argument("-e",
                            "--extension",
                            type=str,
                            default="png",
                            choices=["eps", "pdf", "pgf", "png", "ps", "raw", "rgba", "svg", "svgz"],
                            help="The output file format, default: 'png'")
        parser.add_argument("-validate",
                            action='store_true',
                            help="Validate that the input matrices match "
                                 "with each other and then quit, default: "
                                 "'False'.")

        return parser

    def print(self):
        """
        Method that prints the input settings.

        :return:
        """
        for arg in vars(self.arguments):
            print("Input argument '{}' "
                  "has value '{}'.".format(arg, getattr(self.arguments, arg)))

    def get_argument(self, arg_key):
        """
        Method to return the value of a command line key.

        :param arg_key: str, key in parse command line
        :return: int/str/dict, value of arg_key. If the key did not
                 exist, return None.
        """
        if self.arguments is not None and arg_key in self.arguments:
            value = getattr(self.arguments, arg_key)
        else:
            value = None

        return value

    def get_all_arguments(self):
        """
        Method to return all command line arguments.

        :return self.arguments: dict, dictionary containing all command
                line arguments.
        """
        return self.arguments
