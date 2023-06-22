#!/usr/bin/env python3

"""
File:         decon_eqtl_zscore_jobs.py
Created:      2022/01/24
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
from __future__ import print_function
import argparse
import os

# Third party imports.

# Local application imports.

"""
Syntax:
./decon_eqtl_zscore_jobs.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -of 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron
    
./decon_eqtl_zscore_jobs.py \
    -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/create_matrices/genotype_table.txt.gz \
    -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/2021-12-07-CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt \
    -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/perform_deconvolution/deconvolution_table.txt.gz \
    -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected/combine_gte_files/SampleToDataset.txt.gz \
    -of 2022-02-27-CortexEUR-cis-ForceNormalised-MAF5-NegativeToZero-DatasetAndRAMCorrected
"""

# Metadata
__program__ = "Create Decon-eQTL Z-score Jobs"
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
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.outfolder = getattr(arguments, 'outfolder')

        self.decon_basedir = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/"
        self.job_outdir = os.path.join(self.decon_basedir, "jobs", self.outfolder)
        self.job_output_outdir = os.path.join(self.job_outdir, "output")
        self.time = "05:55:00"
        self.python_executable = os.path.join(self.decon_basedir, "decon_eqtl.py")
        self.genotype_na = -1
        self.allele_configs = "limited"
        self.call_rate = 0.95
        self.hw_pval = 1e-4
        self.maf = 0.05

        for outdir in [self.job_outdir, self.job_output_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

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
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        for zscore in [0, 1, 2, 3, 4, 5]:
            self.create_job_file(job_name=self.outfolder,
                                 zscore=zscore)

    def create_job_file(self, job_name, zscore, leading_zeros=0,
                        n_permutations=0, permutation_index_offset=0, cpus=1,
                        mem=4, nodes=1, qos="regular"):
        name = "{}-GT{}SD".format(job_name, zscore)

        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(name),
                 "#SBATCH --output={}".format(os.path.join(self.job_output_outdir, name + ".out")),
                 "#SBATCH --error={}".format(os.path.join(self.job_output_outdir, name + ".out")),
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
                 "     -zc {} \\".format(zscore),
                 "     -p {} \\".format(n_permutations),
                 "     -po {} \\".format(permutation_index_offset),
                 "     -plz {} \\".format(leading_zeros),
                 "     -od {} \\".format(self.decon_basedir),
                 "     -of {}".format(name),
                 "",
                 "deactivate",
                 ""]

        jobfile_path = os.path.join(self.job_outdir, name + ".sh")
        with open(jobfile_path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(jobfile_path)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Output folder: {}".format(self.outfolder))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
