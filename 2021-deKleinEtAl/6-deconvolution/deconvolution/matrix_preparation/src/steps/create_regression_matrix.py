"""
File:         create_reg_matrix.py
Created:      2020/03/16
Last Changed: 2020/03/19
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
import gzip
import os

# Third party imports.
import numpy as np
from scipy import stats

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists


class CreateRegressionMatrix:
    def __init__(self, settings, eqtl_df, geno_df, alleles_df, expr_df, force,
                 outdir):
        """
        The initializer for the class.

        :param settings: string, the settings.
        :param eqtl_df: DataFrame, the eQTL probes data.
        :param geno_df: DataFrame, the genotype data.
        :param alleles_df: DataFrame, the alleles data.
        :param expr_df: DataFrame, the expression data.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.eqtl_df = eqtl_df
        self.geno_df = geno_df
        self.alleles_df = alleles_df
        self.expr_df = expr_df
        self.force = force

        # Prepare an output directories.
        outdir = os.path.join(outdir, 'create_regression_matrix')
        prepare_output_dir(outdir)
        self.outpath = os.path.join(outdir, "regression_table.txt.gz")

    def start(self):
        print("Starting creating regression file.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.outpath) and not self.force:
            print("Skipping step.")
            return

        # Remove the output files.
        if os.path.isfile(self.outpath):
            print("Removing file: {}.".format(self.outpath))
            os.remove(self.outpath)

        # Prepare string buffer.
        regr_str_buffer = [
            "snp\tprobe\talleles\tminor_allele\tallele_assessed\tflipped\tslope"
            "\tintercept\tcorr_coeff\tp_value\tstd_err\toveral_z_score"
            "\tz_score_estimate\n"]

        # Correlating.
        print("Correlating:")
        nrows = self.eqtl_df.shape[0]
        for i, row in self.eqtl_df.iterrows():
            if (i % 250 == 0) or (i == (nrows - 1)):
                print("\t Processing {}/{} [{:.2f}%]".format(i, nrows,
                                                             (100 / nrows) * i
                                                             )
                      )

                # Write output files.
                self.write_buffer(self.outpath, regr_str_buffer)
                regr_str_buffer = []

            # Extract the usefull information.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            overal_z_score = row["OverallZScore"]
            allele_assessed = row["AlleleAssessed"]

            # Get the data.
            genotype = self.geno_df.iloc[i, :].T.to_frame()
            if snp_name != genotype.columns[0]:
                print("SNPName does not match in genotype subset.")
                break
            expression = self.expr_df.iloc[i, :].T.to_frame()
            if snp_name != expression.columns[0]:
                print("SNPName does not match in expression subset.")
                break
            data = genotype.merge(expression, left_index=True, right_index=True)
            data.columns = ["genotype", "expression"]

            # Remove missing values.
            data.replace(-1, np.nan, inplace=True)
            data.dropna(inplace=True)

            # Determine the alleles.
            (alleles, minor_allele) = self.alleles_df.iloc[i, :]

            # Determine whether to flip or not.
            flipped = False
            if allele_assessed != alleles.split("/")[1]:
                data['genotype'] = 2.0 - data['genotype']
                flipped = True

            # # Naive flip method.
            # if allele_assessed != minor_allele:
            #     data['genotype'] = 2.0 - data['genotype']

            # Calculate the correlation.
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                data["genotype"],
                data["expression"])

            # Calculate the z-score estimate.
            z_score_estimate = slope / std_err

            # # Naive flip method 2.0.
            # if allele_assessed != alleles.split("/")[1]:
            #     z_score_estimate = z_score_estimate * -1

            # Add to the buffer.
            regr_str = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n" \
                .format(snp_name,
                        probe_name,
                        alleles,
                        minor_allele,
                        allele_assessed,
                        flipped,
                        slope,
                        intercept,
                        r_value,
                        p_value,
                        std_err,
                        overal_z_score,
                        z_score_estimate
                        )
            regr_str_buffer.append(regr_str)

        # Write output files.
        if regr_str_buffer:
            self.write_buffer(self.outpath, regr_str_buffer)

    @staticmethod
    def write_buffer(filename, buffer):
        """
        Method for writing a list of strings to a gzipped file,

        :param filename: string, the output file path.
        :param buffer: list, the lines of strings to write.
        """
        # Write output files.
        if os.path.isfile(filename):
            mode = 'ab'
        else:
            mode = 'wb'

        with gzip.open(filename, mode) as f:
            for line in buffer:
                f.write(line.encode())
        f.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype matrix shape: {}".format(self.geno_df.shape))
        print("  > Alleles matrix shape: {}".format(self.alleles_df.shape))
        print("  > Expression matrix shape: {}".format(self.expr_df.shape))
        print("  > Output path: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
