#!/usr/bin/env python3

"""
File:         create_CIA_jobs.py
Created:      2020/10/16
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
from __future__ import print_function
from pathlib import Path
import argparse
import os

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Create CIA jobs"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.job = getattr(arguments, 'job').upper()
        self.input_folder = getattr(arguments, 'input').lower()
        self.output_folder = getattr(arguments, 'output').lower()
        self.settings = getattr(arguments, 'settings').lower()
        self.exclude = getattr(arguments, 'exclude')
        self.start_index = getattr(arguments, 'first')
        self.stop_index = getattr(arguments, 'last')
        self.batch_size = getattr(arguments, 'batch')
        self.n_samples = getattr(arguments, 'n_samples')
        self.cores = 1
        self.mem = getattr(arguments, 'mem')

        if self.output_folder is None:
            self.output_folder = self.input_folder

        # Set the variables.
        self.outdir = os.path.join(Path(__file__).parent.absolute(), self.job)
        self.log_file_outdir = os.path.join(self.outdir, 'output')
        time = getattr(arguments, 'time').lower()
        self.time = None
        if time == "short":
            self.time = "05:59:00"
        elif time == "medium":
            self.time = "23:59:00"
        elif time == "long":
            self.time = "6-23:49:00"

        for outdir in [self.outdir, self.log_file_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-j",
                            "--job",
                            type=str,
                            required=True,
                            help="The name of the job.")
        parser.add_argument("-i",
                            "--input",
                            type=str,
                            default=None,
                            help="The name of the input directory. "
                                 "Default: 'output'.")
        parser.add_argument("-o",
                            "--output",
                            type=str,
                            default=None,
                            help="The name of the output directory. "
                                 " Default: same as -i / --input.")
        parser.add_argument("-s",
                            "--settings",
                            type=str,
                            required=True,
                            help="The settings input file (without '.json'), "
                                 "default: 'default_settings'.")
        parser.add_argument("-f",
                            "--first",
                            type=int,
                            default=0,
                            help="The eQTL index of the first one to analyse, "
                                 "default: 0.")
        parser.add_argument("-l",
                            "--last",
                            type=int,
                            required=True,
                            help="The eQTL index of the last one to analyse")
        parser.add_argument("-b",
                            "--batch",
                            type=int,
                            default=50,
                            help="The number of eQTLs per job. Default: 50.")
        parser.add_argument("-ns",
                            "--n_samples",
                            type=int,
                            required=True,
                            help="The number of samples.")
        parser.add_argument("-t",
                            "--time",
                            type=str,
                            default="short",
                            choices=["short", "medium", "long"],
                            help="The time required.")
        parser.add_argument("-c",
                            "--cores",
                            type=int,
                            default=1,
                            help="The number of cores required.")
        parser.add_argument("-m",
                            "--mem",
                            type=int,
                            default=2,
                            help="The memory required.")
        parser.add_argument("-e",
                            "--exclude",
                            type=str,
                            default=None,
                            help="The name of node to exclude,"
                                 "default: None.")


        return parser.parse_args()

    def start(self):
        start_indices = [i for i in range(self.start_index, self.stop_index, self.batch_size)]
        for job_id, start_index in enumerate(start_indices):
            batch_size = self.batch_size
            if job_id == (len(start_indices) - 1):
                batch_size = self.stop_index - start_index
            self.write_job_file(job_id, start_index, batch_size)

    def write_job_file(self, job_id, start_index, batch_size):
        skip_rows = ""
        if start_index > 0:
            skip_rows = " -sr {}".format(start_index)

        job_name = "{}{}".format(self.job.upper(), job_id)
        out_filepath = os.path.join(self.log_file_outdir, job_name + ".out")
        bash_filepath = os.path.join(self.outdir, job_name + ".sh")

        exclude = ""
        if self.exclude is not None:
            exclude = "#SBATCH --exclude={}".format(self.exclude)

        lines = ["#!/bin/bash\n",
                 "#SBATCH --job-name={}\n".format(job_name),
                 "#SBATCH --output={}\n".format(out_filepath),
                 "#SBATCH --error={}\n".format(out_filepath),
                 "#SBATCH --time={}\n".format(self.time),
                 "#SBATCH --cpus-per-task={}\n".format(self.cores),
                 "#SBATCH --mem={}gb\n".format(self.mem),
                 "#SBATCH --nodes=1\n",
                 "#SBATCH --open-mode=append\n",
                 "#SBATCH --export=NONE\n",
                 "#SBATCH --get-user-env=L\n",
                 "{}\n".format(exclude),
                 "\n",
                 "module load Python/3.6.3-foss-2015b\n",
                 "source $HOME/venv/bin/activate\n",
                 "\n",
                 "python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/custom_interaction_analyser.py -i {} -o {} -s {}{} -ne {} -ns {}\n".format(self.input_folder, self.output_folder, self.settings, skip_rows, batch_size, self.n_samples),
                 "\n",
                 "deactivate\n"]

        with open(bash_filepath, "w") as f:
            for line in lines:
                f.write(line)
        f.close()


if __name__ == '__main__':
    m = main()
    m.start()
