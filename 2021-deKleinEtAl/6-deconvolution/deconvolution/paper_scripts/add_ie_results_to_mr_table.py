#!/usr/bin/env python3

"""
File:         add_ie_results_to_mr_table.py
Created:      2021/02/026
Last Changed: 2022/03/14
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
import time
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
import requests

# Local application imports.

"""
Syntax:
./add_ie_results_to_mr_table.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -t
"""

# Metadata
__program__ = "Add Interaction eQTL results to MR Table"
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
        self.decon_path = getattr(arguments, 'decon_path')
        self.token = getattr(arguments, 'token')

        # Set variables.
        self.table_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/2022-03-14-Supplementary_Table_12_MR_findings_passing_suggestive_threshold.xlsx"
        self.sheet_name = 'TopMR'
        self.table_ea_column = "Effect Allele"
        self.table_decon_snp_column = "eQTL SNP"
        self.table_ld_r2_column = "LD R-squared"
        self.table_n_ieqtls_column = "Number of cell types with significant interaction (FDR<0.05)"
        self.minimal_ld = 0.8

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'add_ie_results_to_mr_table')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
        parser.add_argument("-d",
                            "--decon_path",
                            type=str,
                            required=True,
                            help="The path to the deconvolution "
                                 "results matrix.")
        parser.add_argument("-t",
                            "--token",
                            type=str,
                            required=True,
                            help="The LDPair API token.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")
        table_df = self.load_file(self.table_path, header=0, index_col=None, skiprows=1, sheet_name=self.sheet_name)
        decon_df = self.load_file(self.decon_path, header=0, index_col=None)

        print("Pre-process data.")
        table_df.dropna(axis=0, how="all", inplace=True)
        table_df.index = table_df["Ensembl Gene ID"].astype(str) + "_" + table_df["Chromosome"].astype(int).astype(str) + ":" + table_df["Position"].astype(int).astype(str) + ":" + table_df["SNP"].astype(str)
        print(table_df)
        print(list(table_df.columns))

        snp_info_df = decon_df["SNP"].str.split(":", n=None, expand=True)
        snp_info_df.columns = ["Chromosome", "Position", "SNP", "Alleles"]
        snp_info_df = pd.concat([decon_df[["Gene"]], snp_info_df], axis=1)
        index = snp_info_df["Gene"] + "_" + snp_info_df["Chromosome"] + ":" + snp_info_df["Position"] + ":" + snp_info_df["SNP"]
        decon_df.index = index
        decon_df.columns = [col.replace("beta", "Beta").replace("BH-FDR", "FDR") for col in decon_df.columns]
        print(decon_df)
        print(list(decon_df.columns))

        print("Checking if all columns are present")
        for col in [self.table_ea_column, self.table_decon_snp_column, self.table_ld_r2_column, self.table_n_ieqtls_column]:
            if col not in table_df.columns:
                print("Error, column '{}' is missing".format(col))
                exit()
        cell_types = [x.split(" ")[0] for x in decon_df.columns if "pvalue" in x]
        for ct in cell_types:
            if "{} Beta".format(ct) not in decon_df.columns:
                print("Error, column '{} Beta' is missing".format(ct))
                exit()
            if "{} FDR".format(ct) not in decon_df.columns:
                print("Error, column '{} FDR' is missing".format(ct))
                exit()
        print("\tValid.")

        print("Checking overlap")
        table_indices = set(table_df.index)
        decon_indices = set(decon_df.index)
        overlap = len(table_indices.intersection(decon_indices))
        print("\t{} / {} [{:.2f}%] of the rows have an exact match".format(overlap, len(table_indices), (100 / len(table_indices)) * overlap))
        #
        # print("Constructing the MR SNP link table.")
        # n_rows = table_df.shape[0]
        # n_exact_match = 0
        # n_proxy_found = 0
        # n_missing = 0
        # n_excluded = 0
        # links = []
        # last_print_time = None
        # for i, (index, row) in enumerate(table_df.iterrows()):
        #     # Update user on progress.
        #     now_time = int(time.time())
        #     if (last_print_time is None or (now_time - last_print_time) >= 30 or i == (n_rows - 1)):
        #         print("\t{:,}/{:,} rows processed [{:.2f}%]".format(i, (n_rows - 1), (100 / (n_rows - 1)) * i))
        #         last_print_time = now_time
        #
        #     id = row["ID"]
        #     if index in decon_indices:
        #         links.append([id, index, index, row["SNP"], 1, 1])
        #         n_exact_match += 1
        #     else:
        #         # Check if the instrument snp is in high ld with the ieQTL snp.
        #         gene = row["Ensembl Gene ID"]
        #         instrument_snp = row["SNP"]
        #         effect_allele = row[self.table_ea_column]
        #         if str(row["Proxy SNP"]).startswith("rs"):
        #             instrument_snp = row["Proxy SNP"]
        #
        #         # Get all ieQTLs tested for this gene.
        #         gene_decon_df_subset = decon_df.loc[decon_df["Gene"] == gene, ["Gene", "SNP", "Allele assessed"]].copy()
        #         if len(gene_decon_df_subset.index) > 0:
        #             gene_decon_df_subset["SNP"] = gene_decon_df_subset["SNP"].str.split(":", n=None, expand=True)[2]
        #
        #             match_index = None
        #             match_snp = None
        #             match_aa = None
        #             match_r2 = 0
        #             match_allele_combinations = None
        #             for ld_option_index, ld_option_row in gene_decon_df_subset.iterrows():
        #                 info = self.get_ldpair_info(instrument_snp, ld_option_row["SNP"])
        #
        #                 if info["WARNING"] is np.nan and info["R2"] != np.nan and info["R2"] > match_r2:
        #                     match_index = ld_option_index
        #                     match_snp = ld_option_row["SNP"]
        #                     match_aa = ld_option_row["Allele assessed"]
        #                     match_r2 = info["R2"]
        #                     match_allele_combinations = info["matching"]
        #
        #             if match_r2 >= self.minimal_ld:
        #                 flip = 1
        #                 if [effect_allele, match_aa] not in match_allele_combinations:
        #                     flip = -1
        #                 links.append([id, index, match_index, match_snp, match_r2, flip])
        #                 n_proxy_found += 1
        #             else:
        #                 links.append([id, index, np.nan, np.nan, np.nan, np.nan])
        #                 n_missing += 1
        #         else:
        #             links.append([id, index, np.nan, np.nan, np.nan, np.nan])
        #             n_excluded += 1
        #
        # print("\tN-rows with exact match: {:,} [{:.2f}%]".format(n_exact_match, (100 / n_rows) * n_exact_match))
        # print("\tN-rows with proxy SNP found: {:,} [{:.2f}%]".format(n_proxy_found, (100 / n_rows) * n_proxy_found))
        # print("\tN-rows with NaN: {:,} [{:.2f}%]".format(n_missing, (100 / n_rows) * n_missing))
        # print("\tN-rows whose genes were excluded: {:,} [{:.2f}%]".format(n_excluded, (100 / n_rows) * n_excluded))
        # print("")
        #
        # print("Saving link file")
        # links_df = pd.DataFrame(links, columns=["ID", "table index", "decon index", self.table_decon_snp_column, self.table_ld_r2_column, "flip"])
        # print(links_df)
        # self.save_file(df=links_df,
        #                outpath=os.path.join(self.outdir, "links.txt.gz"),
        #                header=True,
        #                index=False)

        links_df = self.load_file(os.path.join(self.outdir, "links.txt.gz"), header=0, index_col=None)
        print(links_df)

        print("Post-processing")
        # Check indices are correct.
        if (list(table_df["ID"]) != list(links_df["ID"])) or (list(table_df.index) != list(links_df["table index"])):
            print("Error, table indices do not match link dataframe.")
            exit()

        # Reorder the decon frame.
        links_df.index = links_df["table index"]
        decon_subset_df = links_df.merge(decon_df, left_on="decon index", right_index=True, how="left")

        # Flip the beta's.
        for ct in cell_types:
            decon_subset_df["{} Beta".format(ct)] = decon_subset_df["{} Beta".format(ct)] * decon_subset_df["flip"]

        # Calculate the number of significant interactions.
        decon_subset_df[self.table_n_ieqtls_column] = (decon_subset_df.loc[:, ["{} FDR".format(ct) for ct in cell_types]] <= 0.05).sum(axis=1)

        # Change the index to the table index.
        table_df.index = table_df["ID"]
        decon_subset_df.index = decon_subset_df["ID"]

        # Subset the right columns.
        columns_of_interest = [self.table_decon_snp_column, self.table_ld_r2_column] + ["{} Beta".format(ct) for ct in cell_types] + ["{} FDR".format(ct) for ct in cell_types] + [self.table_n_ieqtls_column]
        decon_subset_df = decon_subset_df.loc[:, columns_of_interest]

        # Set the existing columns to nan.
        for col in columns_of_interest:
            if col in table_df.columns:
                table_df[col] = np.nan
            else:
                print("Error, column '{}' not found in table".format(col))

        print(table_df[columns_of_interest])
        print(decon_subset_df)

        print("Merging tables")
        table_df.update(decon_subset_df)
        print(table_df)

        print("")
        print("### Step7 ###")
        print("Save data.")
        # self.save_file(table_df,
        #                outpath=os.path.join(self.outdir, os.path.basename(self.table_path).replace(".xlsx", ".txt.gz")),
        #                index=False)
        self.save_file(table_df,
                       outpath=os.path.join(self.outdir, os.path.basename(self.table_path)),
                       sheet_name=self.sheet_name,
                       na_rep="NA",
                       index=False)

    @staticmethod
    def load_file(path, header, index_col, sep="\t", nrows=None,
                  skiprows=0, sheet_name=None):
        if path.endswith(".xlsx"):
            df = pd.read_excel(path, header=header, index_col=index_col,
                         nrows=nrows, skiprows=skiprows, sheet_name=sheet_name)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t", na_rep="NA",
                  sheet_name="Sheet1"):
        if outpath.endswith('xlsx'):
            df.to_excel(outpath,
                        sheet_name=sheet_name,
                        na_rep=na_rep,
                        header=header,
                        index=index)
        else:
            compression = 'infer'
            if outpath.endswith('.gz'):
                compression = 'gzip'

            df.to_csv(outpath, sep=sep, index=index, header=header,
                      compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

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
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Table path: {}".format(self.table_path))
        print("  > Sheet name: {}".format(self.sheet_name))
        print("  > Table affect allele column: {}".format(self.table_ea_column))
        print("  > Minimal LD: {}".format(self.minimal_ld))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
