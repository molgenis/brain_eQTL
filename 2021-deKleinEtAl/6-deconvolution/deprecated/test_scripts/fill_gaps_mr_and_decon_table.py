#!/usr/bin/env python3

"""
File:         fill_gaps_mr_and_decon_table.py
Created:      2021/02/03
Last Changed: 2021/02/22
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
import random
import re

# Third party imports.
import numpy as np
import pandas as pd
import requests

# Local application imports.

# Metadata
__program__ = "Fill Gaps MR and Decon Table"
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
        self.token = getattr(arguments, 'token')

        # Set the other variables.
        self.table_path = "../data/SupplementaryTable8-MR_findings_passing_suggestive_threshold.xlsx"
        self.sheet_names = ['TopMR']
        self.table_ea = "EA"

        self.decon_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults_alleles_FDR_betas.txt.gz"
        self.decon_ea = "DeconQTLAlleleAssessed"
        self.decon_drop = ["SNPType", "AlleleAssessed"]

        self.minimal_ld = 0.8

        self.outfile = "../data/SupplementaryTable8-MR_findings_passing_suggestive_threshold_gaps_filled.xlsx"

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-t",
                            "--token",
                            type=str,
                            required=True,
                            help="The LDPair API token.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading decon-QTL data.")
        decon_df = self.load_file(self.decon_path, index_col=None)
        decon_df["chrom"] = decon_df["SNPName"].str.split(":", expand=True)[0].astype(np.int64)
        decon_df["pos"] = decon_df["SNPName"].str.split(":", expand=True)[1].astype(np.int64)
        decon_df["SNP"] = decon_df["SNPName"].str.split(":", expand=True)[2]
        decon_df.drop(["SNPName"] + self.decon_drop, axis=1, inplace=True)
        decon_df.rename(columns={"ProbeName": "ENSG", "HGNCName": "gene"}, inplace=True)
        decon_df.columns = [x.replace("CellMapNNLS_", "") for x in decon_df.columns]
        print(decon_df)

        print("")
        print("### Step2 ###")
        print("Loading MR results.")
        sheets_dict = self.load_file(self.table_path, skiprows=1, index_col=None)

        print("")
        print("### Step3 ###")
        print("Filling gaps.")
        n_gaps = 0
        n_filled = 0
        with pd.ExcelWriter(self.outfile) as writer:
            eqtl_snp_data = []
            ld_r2_data = []
            flip_data = []
            for name, sheet in sheets_dict.items():
                if name in self.sheet_names:
                    tmp_sheet = sheet.copy()
                    for i, (index, row) in enumerate(tmp_sheet.iterrows()):
                        decon_cells = row[26:36]
                        if decon_cells.isnull().values.any():
                            print("\teQTL {} - {} on row {} has no decon-QTL values.".format(row["SNP"], row["ENSG"], i))
                            n_gaps += 1

                            # Get the instrument SNP
                            instrument_snp = row["SNP"]

                            # Look the gene up in the decon results.
                            ensg_decon_df = decon_df.loc[decon_df["ENSG"] == row["ENSG"], :].copy()

                            # Test which SNP is in highest LD.
                            best_match = None
                            r2_options = []
                            max_r2 = 0
                            matching = None
                            for _, decon_row in ensg_decon_df.iterrows():
                                info = self.get_ldpair_info(instrument_snp, decon_row["SNP"])
                                if info["WARNING"] is not np.nan:
                                    print(info)
                                r2_options.append(decon_row["SNP"] + " = " + str(info["R2"]))
                                if info["WARNING"] is np.nan and info["R2"] != np.nan and info["R2"] > max_r2:
                                    best_match = decon_row.copy()
                                    max_r2 = info["R2"]
                                    matching = info["matching"]

                            print("\t\tLD R^2 options: {}".format(", ".join(r2_options)))
                            if max_r2 < self.minimal_ld:
                                print("\t\tNo SNP match")
                                eqtl_snp_data.append(np.nan)
                                ld_r2_data.append(np.nan)
                                flip_data.append(np.nan)
                            else:
                                print("\t\tBest SNP match: {} - {} with LD R^2: {}\tMatching: {}-{}".format(best_match["SNP"], best_match["ENSG"], max_r2, row[self.table_ea], best_match[self.decon_ea]))
                                n_filled += 1

                                mr_ea_decon_aa = [row[self.table_ea], best_match[self.decon_ea]]

                                # Flip.
                                flip = False
                                if mr_ea_decon_aa not in matching:
                                    tmp_best_match = best_match.copy()
                                    for name, value in tmp_best_match.iteritems():
                                        if name.endswith("_Beta"):
                                            best_match[name] = value * -1
                                    flip = True

                                # Safe best match.
                                eqtl_snp_data.append(best_match["SNP"])
                                ld_r2_data.append(max_r2)

                                # Insert decon values into excel sheet.
                                sheet.iloc[i, 26:36] = best_match[3:13]
                                flip_data.append(flip)
                        else:
                            eqtl_snp_data.append(row["SNP"])
                            ld_r2_data.append(1)
                            flip_data.append(row["deconQTLFlip"])

                    sheet["deconQTLFlip"] = flip_data
                    sheet.insert(26, "eQTL SNP", eqtl_snp_data)
                    sheet.insert(27, "LD R^2", ld_r2_data)

                    print("\t\tSaving sheet {} with shape {}".format(name, sheet.shape))
                    sheet.to_excel(writer, sheet_name=name, na_rep="NA")

        print("Filled {}/{} gaps.".format(n_filled, n_gaps))

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  skiprows=0, sheet_name=None):
        if path.endswith(".xlsx"):
            df = pd.read_excel(path, header=header, index_col=index_col,
                         nrows=nrows, skiprows=skiprows, sheet_name=sheet_name)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, skiprows=skiprows)
        return df

    def get_ldpair_info(self, rs1, rs2):
        url = 'https://ldlink.nci.nih.gov/LDlinkRest/ldpair?var1={}&var2={}&pop=EUR&token={}'.format(rs1, rs2, self.token)
        r = requests.get(url)

        interest = ["D'", "R2", "Chi-sq", "p-value", "WARNING"]
        info = {x: np.nan for x in interest}
        snp_alleles = []
        haplo_count = 0
        matching = []
        for line in r.text.split("\n"):
            line = line.strip()
            if re.match("[ATCG]_[ATCG]:", line):
                alleles = line.split(":")[0].split("_")

                if haplo_count == 0 or haplo_count == 1:
                    snp_alleles.append(alleles[1])
                elif haplo_count == 2 or haplo_count == 3:
                    snp_alleles.append(alleles[0])

                haplo_count += 1

            for variable in interest:
                if line.startswith("{}: ".format(variable)):
                    value_str = line.replace("{}: ".format(variable), "").replace("<", "").replace(">", "")

                    if variable != "WARNING":
                        value = np.nan
                        try:
                            value = float(value_str)
                        except AttributeError:
                            pass
                    else:
                        value = value_str
                    info[variable] = value

            if re.match("{}\([ATCG]\) allele is correlated with {}\([ATCG]\) allele".format(rs1, rs2), line):
                part1 = re.search("{}(\([ATCG]\))".format(rs1), line).group(1)[1]
                part2 = re.search("{}(\([ATCG]\))".format(rs2), line).group(1)[1]
                matching.append([part1, part2])
            elif re.match("{} and {} are in linkage equilibrium".format(rs1, rs2), line):
                matching = None

        info[rs1] = snp_alleles[2:4]
        info[rs2] = snp_alleles[0:2]
        info["matching"] = matching

        return info

    def print_arguments(self):
        print("Arguments:")
        print("  > Table path: {}".format(self.table_path))
        print("  > Sheet names: {}".format(self.sheet_names))
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Output file: {}".format(self.outfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
