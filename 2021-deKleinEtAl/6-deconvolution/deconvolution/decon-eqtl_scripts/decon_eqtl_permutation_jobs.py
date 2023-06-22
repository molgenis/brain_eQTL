#!/usr/bin/env python3

"""
File:         decon_eqtl_permutation_jobs.py
Created:      2021/07/23
Last Changed: 2022/01/25
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
import argparse
import os

# Third party imports.
import numpy as np

# Local application imports.

"""
Syntax:
./decon_eqtl_permutation_jobs.py \
    -py /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -od /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts \
    -of 2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron \
    -j 125 \
    -t 23:59:00
    
./decon_eqtl_permutation_jobs.py \
    -py /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -maf 0.05 \
    -zc 4 \
    -ac complete \
    -od /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts \
    -of 2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected \
    -j 50 \
    -t 6-23:59:00

"""

# Metadata
__program__ = "Create Decon-eQTL With Permutation Jobs"
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
        self.python_executable = getattr(arguments, 'python_executable')
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.genotype_na = getattr(arguments, 'genotype_na')
        self.allele_configs = getattr(arguments, 'allele_configurations')
        self.call_rate = getattr(arguments, 'call_rate')
        self.hw_pval = getattr(arguments, 'hardy_weinberg_pvalue')
        self.maf = getattr(arguments, 'minor_allele_frequency')
        self.zscore_cutoff = getattr(arguments, 'zscore_cutoff')
        self.outdir = getattr(arguments, 'outdir')
        self.outfolder = getattr(arguments, 'outfolder')
        self.n_permutations = getattr(arguments, 'permutations')
        self.n_jobs = getattr(arguments, 'jobs')
        self.time = getattr(arguments, 'time')

        self.base_outdir = os.path.join(self.outdir, "decon_eqtl", self.outfolder)
        self.job_outdir = os.path.join(self.base_outdir, "jobs")
        self.job_output_outdir = os.path.join(self.job_outdir, "output")
        for outdir in [self.job_outdir, self.job_output_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        if self.n_permutations % self.n_jobs != 0:
            print("-p / --permutations must be a multiplication of -j / --jobs.")
            exit()

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-py",
                            "--python_executable",
                            type=str,
                            required=True,
                            help="The path to the executable.")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
        parser.add_argument("-cc",
                            "--cell_counts",
                            type=str,
                            required=True,
                            help="The path to the cell counts matrix.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-na",
                            "--genotype_na",
                            type=str,
                            required=False,
                            default=-1,
                            help="The genotype value that equals a missing "
                                 "value. Default: -1.")
        parser.add_argument("-ac",
                            "--allele_configurations",
                            type=str,
                            choices=["complete", "limited"],
                            default="complete",
                            help="The number of allele encoding configurations"
                                 "to test. Either 'complete' for all or"
                                 "'limited' for restricting a max of 1 flipped"
                                 "allele. Default: complete.")
        parser.add_argument("-cr",
                            "--call_rate",
                            type=float,
                            required=False,
                            default=0.95,
                            help="The minimal call rate of a SNP (per dataset)."
                                 "Equals to (1 - missingness). "
                                 "Default: >= 0.95.")
        parser.add_argument("-hw",
                            "--hardy_weinberg_pvalue",
                            type=float,
                            required=False,
                            default=1e-4,
                            help="The Hardy-Weinberg p-value threshold."
                                 "Default: >= 1e-4.")
        parser.add_argument("-maf",
                            "--minor_allele_frequency",
                            type=float,
                            required=False,
                            default=0.01,
                            help="The MAF threshold. Default: >0.01.")
        parser.add_argument("-zc",
                            "--zscore_cutoff",
                            type=int,
                            required=False,
                            default=None,
                            help="Remove outliers with z-score >X on one or"
                                 "more cell fractions. Default: None.")
        parser.add_argument("-od",
                            "--outdir",
                            type=str,
                            required=True,
                            help="The name of the output path.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=True,
                            help="The name of the output folder.")
        parser.add_argument("-p",
                            "--permutations",
                            type=int,
                            default=1000,
                            help="The number of permutations to run. "
                                 "Default: 0.")
        parser.add_argument("-j",
                            "--jobs",
                            type=int,
                            default=100,
                            help="The number of jobs to make. "
                                 "Default: 100.")
        parser.add_argument("-t",
                            "--time",
                            type=str,
                            default="05:59:00",
                            choices=["05:59:00", "23:59:00", "6-23:59:00"],
                            help="The time for each job.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Creating standard job files")
        self.create_job_file(job_name=self.outfolder)

        print("Creating permutation job files:")
        n_permutations = int(self.n_permutations / self.n_jobs)
        value_length = len(str(self.n_permutations - n_permutations))
        for permutation_index_offset in np.arange(0, self.n_permutations, n_permutations):
            job_name = "{}_{:0{}d}_{:0{}d}".format(self.outfolder, permutation_index_offset, value_length, permutation_index_offset + n_permutations - 1, value_length)
            self.create_job_file(job_name=job_name,
                                 leading_zeros=value_length - len(str(permutation_index_offset)),
                                 n_permutations=n_permutations,
                                 permutation_index_offset=permutation_index_offset,
                                 outdir=self.job_outdir,
                                 )

    def create_job_file(self, job_name, leading_zeros=0, n_permutations=0,
                        permutation_index_offset=0, cpus=1, mem=4, nodes=1,
                        qos="regular", outdir=None):
        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(job_name),
                 "#SBATCH --output={}".format(os.path.join(self.job_output_outdir, job_name + ".out")),
                 "#SBATCH --error={}".format(os.path.join(self.job_output_outdir, job_name + ".out")),
                 "#SBATCH --time={}".format(self.time),
                 "#SBATCH --cpus-per-task {}".format(cpus),
                 "#SBATCH --mem {}gb".format(mem),
                 "#SBATCH --nodes {}".format(nodes),
                 "#SBATCH --qos={}".format(qos),
                 "",
                 "module load Python/3.7.4-GCCcore-7.3.0-bare",
                 "source $HOME/env/bin/activate",
                 "",
                 "{} \\".format(self.python_executable),
                 "     -ge {} \\".format(self.geno_path),
                 "     -ex {} \\".format(self.expr_path),
                 "     -cc {} \\".format(self.cc_path),
                 "     -std {} \\".format(self.std_path),
                 "     -na {} \\".format(self.genotype_na),
                 "     -ac {} \\".format(self.allele_configs),
                 "     -cr {} \\".format(self.call_rate),
                 "     -hw {} \\".format(self.hw_pval),
                 "     -maf {} \\".format(self.maf),
                 "     -zc {} \\".format(self.zscore_cutoff),
                 "     -p {} \\".format(n_permutations),
                 "     -po {} \\".format(permutation_index_offset),
                 "     -plz {} \\".format(leading_zeros),
                 "     -od {} \\".format(self.outdir),
                 "     -of {}".format(self.outfolder),
                 "",
                 "deactivate",
                 ""]

        if outdir is None:
            outdir = self.base_outdir

        jobfile_path = os.path.join(outdir, job_name + ".sh")
        with open(jobfile_path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(jobfile_path)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Python executable path: {}".format(self.python_executable))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > N permutations: {}".format(self.n_permutations))
        print("  > Genotype NaN: {}".format(self.genotype_na))
        print("  > Allele configs: {}".format(self.allele_configs))
        print("  > MAF: {}%".format(self.maf))
        print("  > Z-score cut-off: {}".format(self.zscore_cutoff))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Output folder: {}".format(self.outfolder))
        print("  > Jobs output folder: {}".format(self.job_outdir))
        print("  > N jobs: {}".format(self.n_jobs))
        print("  > Time: {}".format(self.time))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
