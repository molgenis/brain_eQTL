#!/usr/bin/env python3

"""
File:         add_ie_results_to_mr_table.py
Created:      2021/02/026
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
from pathlib import Path
import argparse
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
from statsmodels.stats import multitest
import requests

# Local application imports.

# Metadata
__program__ = "Merge MR and decon results"
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

        # Set variables.
        self.input_file = "Supplementary Table 11 - MR findings passing suggestive threshold.xlsx"
        self.indir = "../data/"
        self.sheet_name = 'TopMR'
        self.table_ea = "Effect Allele"
        self.decon_path = "../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv"
        self.alleles_path = "../matrix_preparation/cortex_eur_cis/create_matrices/genotype_alleles.txt.gz"
        self.drop_columns = ['eQTL SNP',
                             'LD R-squared',
                             'Astrocyte Beta',
                             'EndothelialCell Beta',
                             'Macrophage Beta',
                             'Neuron Beta',
                             'Oligodendrocyte Beta',
                             'Astrocyte FDR',
                             'EndothelialCell FDR',
                             'Macrophage FDR',
                             'Neuron FDR',
                             'Oligodendrocyte FDR',
                             'Number of cell types with significant interaction (FDR<0.05)',]
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
        parser.add_argument("-t",
                            "--token",
                            type=str,
                            required=True,
                            help="The LDPair API token.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("")
        print("### Step1 ###")
        print("Loading MR table.")
        df = self.load_file(os.path.join(self.indir, self.input_file), skiprows=1, index_col=None, sheet_name=self.sheet_name)
        order = df.columns
        df.drop(self.drop_columns, axis=1, inplace=True)

        print("")
        print("### Step2 ###")
        print("Loading deconvolution table.")
        decon_df, indices = self.preprocess_decon_df(self.load_file(self.decon_path))

        # Add alleles.
        alleles_df = self.load_file(self.alleles_path, index_col=None)
        alleles_df.columns = ["SNP", "Alleles", "MinorAllele"]
        snp_to_alles_dict = dict(zip(alleles_df["SNP"], alleles_df["Alleles"]))
        decon_df["DeconAlleles"] = decon_df["SNPName"].map(snp_to_alles_dict)
        decon_df["DeconQTLAlleleAssessed"] = decon_df["DeconAlleles"].str.split("/", n=1, expand=True)[1]
        del alleles_df

        print("")
        print("### Step3 ###")
        print("Set identical index.")
        df.index = df["Ensembl Gene ID"].astype(str) + df["SNP"].astype(str)
        decon_df.index = decon_df["ProbeName"].astype(str) + decon_df["SNP"].astype(str)

        overlap = len(set(df.index).intersection(set(decon_df.index)))
        print("  {} / {} [{:.2f}%] of the rows have a match".format(overlap, df.shape[0], (100/df.shape[0])*overlap))

        print("")
        print("### Step5 ###")
        print("Construct merge data frame.")
        merge_data = []
        filled = 0
        for i, (index, row) in enumerate(df.iterrows()):
            # if row["Gene"] not in ["CYP24A1", "CLECL1"]:
            #     continue

            print("\tWorking on {}/{} {:.2f}%]".format(i, df.shape[0], (100 / df.shape[0])*i))

            instrument_snp = row["SNP"]

            # check if index is in decon table.
            if index in decon_df.index:
                print("\t  Found match in deconvolution table.")

                decon_data = decon_df.loc[index, :]

                beta_values = decon_data.loc[indices["beta"]].copy()
                if decon_data["DeconQTLAlleleAssessed"] != row[self.table_ea]:
                    beta_values = beta_values * -1

                fdr_values = decon_data.loc[indices["FDR"]]
                n_signif = fdr_values[fdr_values < 0.05].count()

                merge_data.append([decon_data["SNP"], 1] + list(beta_values.values) + list(fdr_values.values) + [n_signif])
                filled += 1
            else:
                # Check if there is a decon result with a different SNP.
                gene_decon_df_subset = decon_df.loc[decon_df["ProbeName"] == row["Ensembl Gene ID"],:].copy()

                if len(gene_decon_df_subset.index) > 0:
                    # Test if there is a SNP in LD.
                    best_match = None
                    max_r2 = 0
                    for _, ld_option in gene_decon_df_subset.iterrows():
                        info = self.get_ldpair_info(instrument_snp,
                                                    ld_option["SNP"])
                        if info["WARNING"] is np.nan and info["R2"] != np.nan and info["R2"] > max_r2:
                            info["data"] = ld_option
                            best_match = info
                            max_r2 = info["R2"]

                    if max_r2 >= self.minimal_ld:
                        # LD is high enough.
                        print("\t  Found SNP in deconvolution table that is in high LD with instrument SNP.")

                        decon_data = best_match["data"]

                        mr_ea_decon_aa = [row[self.table_ea], decon_data["DeconQTLAlleleAssessed"]]

                        beta_values = decon_data.loc[indices["beta"]].copy()
                        if mr_ea_decon_aa not in best_match["matching"]:
                            beta_values = beta_values * -1

                        fdr_values = decon_data.loc[indices["FDR"]]
                        n_signif = fdr_values[fdr_values < 0.05].count()

                        merge_data.append([decon_data["SNP"], best_match["R2"]] + list(beta_values.values) + list(fdr_values.values) + [
                                              n_signif])
                        filled += 1
                        continue
                    else:
                        # No LD match found.
                        pass

                print("\t  No match.")
                merge_data.append([np.nan] * 12 + [0])
        exit()
        print("  {} / {} [{:.2f}%] of the rows are matched.".format(filled, df.shape[0], (100/df.shape[0])*filled))
        merge_df = pd.DataFrame(merge_data, index=df.index, columns=self.drop_columns)
        merge_df.to_pickle(os.path.join(self.outdir, "merge_df.pkl"))

        print("")
        print("### Step6 ###")
        print("Combine data.")
        combined_df = pd.concat([df, merge_df], axis=1)
        combined_df = combined_df.loc[:, order]
        combined_df.reset_index(drop=True, inplace=True)
        print(combined_df)

        print("")
        print("### Step7 ###")
        print("Save data.")
        self.save_file(combined_df, os.path.join(self.outdir, self.input_file), sheet_name=self.sheet_name, na_rep="NA", index=False)

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

    @staticmethod
    def save_file(df, outpath, sheet_name=None, na_rep=None, header=True, index=True):
        df.to_excel(outpath, sheet_name=sheet_name, na_rep=na_rep, header=header, index=index)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def preprocess_decon_df(self, df):
        beta_indices = []
        fdr_indices = []

        # Extract the p_values and convert to fdr.
        tmp = df.copy()
        data = []
        indices = []
        for col in tmp.columns:
            if col.endswith("_pvalue"):
                data.append(multitest.multipletests(tmp.loc[:, col], method='fdr_bh')[1])
                indice_name = col.replace("CellMapNNLS_", "").replace("_pvalue", "") + " FDR"
                fdr_indices.append(indice_name)
                indices.append(indice_name)
            elif col.endswith(":GT"):
                data.append(tmp.loc[:, col])
                indice_name = col.split("_")[2].split(":")[0] + " Beta"
                beta_indices.append(indice_name)
                indices.append(indice_name)
        fdr_beta_df = pd.DataFrame(data, index=indices, columns=tmp.index).T

        del tmp

        # Split the index.
        index_data = []
        for index in fdr_beta_df.index:
            probe_name = index.split("_")[0]
            snp_name = "_".join(index.split("_")[1:])
            chromosome, position, rs, _ = snp_name.split(":")
            index_data.append([probe_name, snp_name, chromosome, position, rs])
        index_df = pd.DataFrame(index_data, index=fdr_beta_df.index, columns=["ProbeName", "SNPName", "Chromosome", "Position", "SNP"])
        fdr_beta_df = index_df.merge(fdr_beta_df, left_index=True, right_index=True)

        return fdr_beta_df, {"beta": beta_indices, "FDR": fdr_indices, "all": indices}

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
        print("  > Table path: {}".format(os.path.join(self.indir, self.input_file)))
        print("  > Sheet name: {}".format(self.sheet_name))
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Alleles path: {}".format(self.alleles_path))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
