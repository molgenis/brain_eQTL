#!/usr/bin/env python3

"""
File:         create_jobs.py
Created:      2020/06/05
Last Changed: 2020/06/08
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
__program__ = "Create jobs"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.job = getattr(arguments, 'job').upper()
        self.name = getattr(arguments, 'name').lower()
        self.settings = getattr(arguments, 'settings').lower()
        self.disease = ' '.join(getattr(arguments, 'disease'))
        self.alpha = getattr(arguments, 'alpha')
        self.exclude = getattr(arguments, 'exclude')
        self.types = getattr(arguments, 'types')

        # Set the variables.
        self.outdir = Path(__file__).parent.absolute()
        self.log_file_outdir = os.path.join(self.outdir, 'output')

        if not os.path.exists(self.log_file_outdir):
            os.makedirs(self.log_file_outdir)

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
                            help="The name of the job")
        parser.add_argument("-n",
                            "--name",
                            type=str,
                            required=True,
                            help="The name of the input/output directory.")
        parser.add_argument("-s",
                            "--settings",
                            type=str,
                            required=True,
                            help="The settings input file")
        parser.add_argument("-d",
                            "--disease",
                            nargs="+",
                            type=str,
                            default="",
                            help="The name of the disease to analyse,"
                                 "default: '' (i.e. SNPs).")
        parser.add_argument("-a",
                            "--alpha",
                            type=float,
                            default=0.05,
                            help="The significance cut-off,"
                                 "default: 0.05.")
        parser.add_argument("-e",
                            "--exclude",
                            type=str,
                            default=None,
                            help="The name of name to exclude,"
                                 "default: None.")
        parser.add_argument("-t",
                            "--types",
                            nargs="+",
                            type=str,
                            default=["all"],
                            choices=["prepare",
                                     "combine",
                                     "identify",
                                     "visualise"],
                            help="The job files to create.")

        return parser.parse_args()

    def start(self):
        if 'prepare' in self.types or 'all' in self.types:
            self.create_prepare_job()

        if 'combine' in self.types or 'all' in self.types:
            self.create_combine_job()

        if 'identify' in self.types or 'all' in self.types:
            self.create_identify_job()

        if 'visualise' in self.types or 'all' in self.types:
            self.create_visualise_job()

    def create_prepare_job(self):
        job_name = "{}_{}".format(self.job, "prepare")

        disease_str = ""
        settings_str = self.settings
        if self.disease != "":
            disease_str = " -d {} ".format(self.disease)
            settings_str = "disease_{}".format(self.settings)

        header = self.create_header(job_name, cpus="2", mem="64")
        content = ['python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation.py -n {} -s cis_{}{}\n'.format(self.name, settings_str, disease_str)]
        footer = self.create_footer()

        fpath = os.path.join(self.outdir, job_name + ".sh")
        lines = header + content + footer
        self.create_file(fpath, lines)

    def create_combine_job(self):
        job_name = "{}_{}".format(self.job, "combine")

        header = self.create_header(job_name, cpus="4", mem="16")
        content = ['python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser.py -n {} -s {} -combine \n'.format(self.name, self.settings)]
        footer = self.create_footer()

        fpath = os.path.join(self.outdir, job_name + ".sh")
        lines = header + content + footer
        self.create_file(fpath, lines)

    def create_identify_job(self):
        job_name = "{}_{}".format(self.job, "identify")

        header = self.create_header(job_name, cpus="4", mem="16")
        content = ['python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/identify_ct_mediated_eqtls.py -n {} -s {} -a {} -e png pdf\n'.format(self.name, self.settings, self.alpha)]
        footer = self.create_footer()

        fpath = os.path.join(self.outdir, job_name + ".sh")
        lines = header + content + footer
        self.create_file(fpath, lines)

    def create_visualise_job(self):
        job_name = "{}_{}".format(self.job, "visualise")

        header = self.create_header(job_name, cpus="2", mem="4")
        content = ['EXTENSIONS=("png" "pdf")\n',
                   '\n',
                   'for EXT in ${EXTENSIONS[@]}; do\n',
                   '\tpython3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/visualiser.py -n {} -s {} -a {} -p covariate_clustermap covariate_comparison covariates_explained_by_others deconvolution_covariate_comparison deconvolution_zscore_comparison inter_clustermap inter_eqtl_effect_celltype inter_pvalue_boxplot inter_zscore_bars -e $EXT\n'.format(self.name, self.settings, self.alpha),
                   '\tpython3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/visualiser.py -n {} -s {} -a {} -p simple_eqtl_effect inter_eqtl_zscore_bars inter_eqtl_effect inter_eqtl_celltype_details -i -e $EXT\n'.format(self.name, self.settings, self.alpha),
                   'done\n']
        footer = self.create_footer()

        fpath = os.path.join(self.outdir, job_name + ".sh")
        lines = header + content + footer
        self.create_file(fpath, lines)

    def create_header(self, job_name, time="05:59:00", cpus="1", mem="4"):
        out_filepath = os.path.join(self.log_file_outdir, job_name + ".out")

        exclude_str = ""
        if self.exclude is not None:
            exclude_str = "#SBATCH --exclude={}".format(self.exclude.lower())

        lines = ["#!/bin/bash\n",
                 "#SBATCH --job-name={}\n".format(job_name),
                 "#SBATCH --output={}\n".format(out_filepath),
                 "#SBATCH --error={}\n".format(out_filepath),
                 "#SBATCH --time={}\n".format(time),
                 "#SBATCH --cpus-per-task={}\n".format(cpus),
                 "#SBATCH --mem={}gb\n".format(mem),
                 "#SBATCH --nodes=1\n",
                 "#SBATCH --open-mode=append\n",
                 "#SBATCH --export=NONE\n",
                 "#SBATCH --get-user-env=L\n",
                 "{}\n".format(exclude_str),
                 "\n",
                 "module load Python/3.6.3-foss-2015b\n",
                 "source $HOME/venv/bin/activate\n",
                 "\n"]

        return lines

    @staticmethod
    def create_footer():
        return ["\n", "deactivate\n"]

    @staticmethod
    def create_file(fpath, lines):
        with open(fpath, "w") as f:
            for line in lines:
                f.write(line)
        f.close()


if __name__ == '__main__':
    m = main()
    m.start()
