#!/opt/conda/envs/eQTLGenPopAssign/bin/python
"""
File:         print_WGS_VCF_filter_overview.py
Created:      2022/10/20
Last Changed: 2022/11/11
Author:       M.Vochteloo

Copyright (C) 2022 M.Vochteloo
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
import gzip
import os
import re

# Third party imports.
import pandas as pd
import numpy as np

# Local application imports.

# Metadata
__program__ = "Print WGS VCF Filter Overview"
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

"""
Syntax: 

### AMP-AD ###
./print_WGS_VCF_filter_overview.py \
    --workdir /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-11-03-FilteredGenotypes \
    --vcf_file_format NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_CHR.recalibrated_variants.vcf.gz \
    --exclude Y others
"""


CHROMOSOMES = [str(chr) for chr in range(1, 23)] + ["X", "Y", "others"]


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.workdir = getattr(arguments, 'workdir')
        self.vcf_file_format = getattr(arguments, 'vcf_file_format')
        self.exclude = getattr(arguments, 'exclude')
        self.chromosome = getattr(arguments, 'chr')

        self.chromosomes = [CHR for CHR in CHROMOSOMES if (self.exclude is None) or (CHR not in self.exclude)]

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add other arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("--workdir",
                            type=str,
                            required=True,
                            help="The path to the working directory")
        parser.add_argument("--vcf_file_format",
                            type=str,
                            required=True,
                            help="The file format of the VCF files. CHR"
                                 "will be replaced with the chromosome number.")
        parser.add_argument("--chr",
                            type=str,
                            required=True,
                            help="the chromosome number.")
        parser.add_argument("--exclude",
                            nargs="*",
                            type=str,
                            choices=CHROMOSOMES,
                            default=None,
                            help="Exclude certain chromosomes.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Parsing chromosomes:")
        data = []
        prev_filter_thresh = None
        prev_failed_var_thresh = None

        row = [self.chromosome, False, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, np.nan, False, ""]

        job_logfile_path = os.path.join("custom_vcf_filter.log")
        if os.path.exists(job_logfile_path):
            row[1] = True
            row[15:] = self.read_job_logfile(filepath=job_logfile_path)

        # Parse thresholds
        job_file_path_prefix = os.path.join(self.workdir, self.vcf_file_format.replace(".vcf.gz", "-filtered"))
        job_vcffile_path = job_file_path_prefix + ".vcf.gz"
        job_logfile_path = job_file_path_prefix + ".log.gz"
        print("JOB LOG FILE")
        print(job_logfile_path)
        if os.path.exists(job_vcffile_path) and os.path.exists(job_logfile_path):
            filter_thresh = self.parse_thresholds(
                filepath=job_vcffile_path
            )
            # first check the settings that matter.
            failed_var_thresh = dict((key, filter_thresh[key]) for key in ('tresh_MAF', 'tresh_CR', 'tresh_HWE'))
            if prev_failed_var_thresh is not None:
                if failed_var_thresh != prev_failed_var_thresh:
                    print("Error, not all VCF files are filtered with "
                          "the same MAF, CR, or HWE settings.")
                    exit()

            # give warning if other settings differ.
            if prev_filter_thresh is not None:
                if filter_thresh != prev_filter_thresh:
                    print("Warning, not all VCF files are filtered with "
                          "the same settings.")

            row[2:15] = self.read_filter_logfile(
                filepath=job_logfile_path,
                thresh_maf=failed_var_thresh["tresh_MAF"],
                thresh_cr=failed_var_thresh["tresh_CR"],
                thresh_hwe=failed_var_thresh["tresh_HWE"]
            )

            filter_thresh["Chromosome"] = self.chromosome
            filter_thresh_df = pd.DataFrame.from_dict({key: [value] for key, value in filter_thresh.items()})
            self.save_file(df=filter_thresh_df, outpath=os.path.join(self.workdir, "VCFFilterSettings.txt"))

        data.append(row)

        df = pd.DataFrame(data, columns=["Chromosome", "Started", "MultiAllelic", "IndelBelowVQSR", "IndelNonPass", "SNVBelowVQSR", "SNVNonPass", "IncorrectInbreedingCoeff", "BelowInbreedingCoeff", "NoGTCol", "FailedPrefilterVarStats", "FailedCR", "FailedMAF", "FailedHWE", "PASSQC", "Parsed", "Written", "PctKept", "Finished", "Elapsed"])
        df.iloc[:, 2:17] = df.iloc[:, 2:17].astype(int)

        self.save_file(
            df=df,
            outpath=os.path.join(self.workdir, "VCFFilterSummaryStats.txt")
        )

        for column in df.columns:
            if column in ["Chromosome", "Started", "PctKept", "Finished", "Elapsed"]:
                continue
            df[column] = df[column].map('{:,.0f}'.format)

        print("\nVCF filter overview:")
        print(df.loc[:, ["Chromosome", "Started", "Parsed", "Written", "PctKept", "Finished", "Elapsed"]])
        print("\nFailed variant stats:")
        print(df.loc[:, ["Chromosome", "MultiAllelic", "IndelBelowVQSR", "IndelNonPass", "SNVBelowVQSR", "SNVNonPass", "IncorrectInbreedingCoeff", "BelowInbreedingCoeff", "NoGTCol", "FailedPrefilterVarStats", "FailedCR", "FailedMAF", "FailedHWE", "PASSQC"]])

    @staticmethod
    def parse_thresholds(filepath):
        filter_thresh = {}
        prev_line = None
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    settings_line = "".join(prev_line.split(" ")[1:]).replace("\n", "")
                    settings = [(argument, value) for argument, value in [x.split("=") for x in settings_line.split(";")]]
                    for argument, value in settings:
                        try:
                            value = float(value)
                        except ValueError:
                            try:
                                value = bool(value)
                            except ValueError:
                                pass
                        filter_thresh[argument] = value
                    break

                if not line.startswith("##"):
                    break

                prev_line = line
        f.close()

        if "tresh_MAF" not in filter_thresh or "tresh_CR" not in filter_thresh or "tresh_HWE" not in filter_thresh:
            print("Error, missing required thresholds MAF, CR, or HWE in filter settings.")
            exit()

        return filter_thresh

    @staticmethod
    def read_job_logfile(filepath):
        with open(filepath, 'r') as f:
            parsed = 0
            written = 0
            pct_kept = np.nan
            finished = False
            elapsed = ""
            for line in f:
                line = line.replace(",", "")
                if re.match("([0-9]+) lines parsed ([0-9]+) written", line):
                    match = re.search("([0-9]+) lines parsed ([0-9]+) written", line)
                    parsed = int(match.group(1))
                    written = int(match.group(2))
                    pct_kept = np.round((written / parsed) * 100, 1)
                elif re.match("Done. How about that!", line):
                    finished = True
                elif re.match("[0-9]{7}.batch. ([0-9]{2}:[0-9]{2}:[0-9]{2})", line):
                    match = re.search("[0-9]{7}.batch. ([0-9]{2}:[0-9]{2}:[0-9]{2})", line)
                    elapsed = match.group(1)
                else:
                    pass
        f.close()

        return [parsed, written, pct_kept, finished, elapsed]

    @staticmethod
    def read_filter_logfile(filepath, thresh_maf=0.05, thresh_cr=0.99,
                            thresh_hwe=1E-6):
        multi_allelic = 0
        indel_below_vqsr = 0
        indel_non_pass = 0
        snv_below_vqsr = 0
        snv_non_pass = 0
        incorrect_inbreeding_coeff = 0
        below_inbreeding_coeff = 0
        no_gt_col = 0
        failed_pre_filter_var_stats = 0
        n_failed_cr = 0
        n_failed_maf = 0
        n_failed_hwe = 0
        pass_qc = 0

        try:
            with gzip.open(filepath, 'rt') as f:
                for line in f:
                    (_, reason, stats, sec) = line.split("\t")
                    if reason == "MultiAllelic":
                        multi_allelic += 1
                    elif reason == "IndelBelowVQSR":
                        indel_below_vqsr += 1
                    elif reason == "IndelNonPass":
                        indel_non_pass += 1
                    elif reason == "SNVBelowVQSR":
                        snv_below_vqsr += 1
                    elif reason == "SNVNonPass":
                        snv_non_pass += 1
                    elif reason.startswith("IncorrectInbreedingCoeff"):
                        incorrect_inbreeding_coeff += 1
                    elif reason.startswith("BelowInbreedingCoeff"):
                        below_inbreeding_coeff += 1
                    elif reason == "NoGTCol":
                        no_gt_col += 1
                    elif reason == "FailedPrefilterVarStats":
                        failed_pre_filter_var_stats += 1

                        splitted_stats = stats.split(";")
                        cr = float(splitted_stats[0].replace("CR=", ""))
                        maf = float(splitted_stats[1].replace("MAF=", ""))
                        hwe = float(splitted_stats[2].replace("HWE=", ""))

                        if cr <= thresh_cr:
                            n_failed_cr += 1
                        if maf <= thresh_maf:
                            n_failed_maf += 1
                        if hwe != -1 and hwe <= thresh_hwe:
                            n_failed_hwe += 1
                    elif reason == "PASSQC":
                        pass_qc += 1
                    else:
                        pass
            f.close()
        except EOFError:
            pass

        return [multi_allelic, indel_below_vqsr, indel_non_pass,
                snv_below_vqsr, snv_non_pass, incorrect_inbreeding_coeff,
                below_inbreeding_coeff, no_gt_col, failed_pre_filter_var_stats,
                n_failed_cr, n_failed_maf, n_failed_hwe, pass_qc]

    @staticmethod
    def save_file(df, outpath, header=True, index=False, sep=","):
        print(df)
        df.to_csv(outpath, sep=sep, index=index, header=header)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Working directory: {}".format(self.workdir))
        print("  > vcf_file_format: {}".format(self.vcf_file_format))
        print("  > Chromosome: {}".format(self.chromosome))
        print("  > Exclude: {}".format(self.exclude))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
