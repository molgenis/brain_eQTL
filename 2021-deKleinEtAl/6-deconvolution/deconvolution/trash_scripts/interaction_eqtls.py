#!/usr/bin/env python3

"""
File:         interaction_eqtls.py
Created:      2020/05/01
Last Changed: 2020/05/22
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
import os

# Third party imports.
import scipy.stats as st
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Interaction eQTLs"
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
        self.input = {'cis_new_output': 3.60e-04,
                      'alzheimer_cis_new_output': 0.05,
                      'schizophrenia_cis_new_output': 0.05,
                      'parkinson_cis_new_output': 0.05,
                      'depression_cis_new_output': 0.05,
                      'ms_cis_new_output': 0.05}
        self.cols_of_interest = ["SEX", "CellMapNNLS_Astrocyte",
                                 "CellMapNNLS_EndothelialCell",
                                 "CellMapNNLS_Macrophage", "CellMapNNLS_Neuron",
                                 "CellMapNNLS_Oligodendrocyte"]
        self.maf_cutoff = 0.05
        self.max_url_len = 8190
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'interaction_eqtls')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        for dir, pval_cutoff in self.input.items():
            if dir != "cis_new_output":
                continue
            print("Directory: {}".format(dir))

            z_score_cutoff = abs(st.norm.isf(pval_cutoff))

            # Define filenames.
            inter_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/{}/covariates/interaction_table.txt.gz".format(
                dir)
            tvalue_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/{}/covariates/inter_tvalue_table.txt.gz".format(
                dir)
            eqtl_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/{}/combine_eqtlprobes/eQTLprobes_combined.txt.gz".format(
                dir)
            geno_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/{}/create_matrices/genotype_table.txt.gz".format(
                dir)
            alleles_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/{}/create_matrices/genotype_alleles.txt.gz".format(
                dir)
            outdir = os.path.join(self.outdir, dir)

            if not os.path.exists(outdir):
                os.makedirs(outdir)

            print("Loading dataframes.")
            inter_df = pd.read_csv(inter_path, sep="\t", header=0, index_col=0)
            print("\tInteraction matrix: {}".format(inter_df.shape))

            tvalue_df = pd.read_csv(tvalue_path, sep="\t", header=0,
                                    index_col=0)
            print("\tT-value matrix: {}".format(tvalue_df.shape))

            eqtl_df = pd.read_csv(eqtl_path, sep="\t", header=0)
            print("\teQTL matrix: {}".format(eqtl_df.shape))

            geno_df = pd.read_csv(geno_path, sep="\t", header=0, index_col=0)
            print("\tgenotype matrix: {}".format(geno_df.shape))

            alleles_df = pd.read_csv(alleles_path, sep="\t", header=0,
                                     index_col=0)
            print("\talleles matrix: {}".format(alleles_df.shape))

            # Check if the files math up.
            for i in range(inter_df.shape[1]):
                snp_name = eqtl_df["SNPName"][i]
                if not inter_df.columns[i].startswith(snp_name) or not \
                        tvalue_df.columns[i].startswith(snp_name):
                    print("Input files do not match (1).")
                    exit()
                if not geno_df.index[i] == snp_name:
                    print("Input files do not match (2).")
                    exit()
                if not alleles_df.index[i] == snp_name:
                    print("Input files do not match (3).")
                    exit()

            # Create a flip column.
            print("Getting genotype info.")
            eqtl_geno_info = []
            for i, (_, row) in enumerate(eqtl_df.iterrows()):
                geno_flip = 1
                eqtl_flip = 1

                # Determine the minor allele genotype and letter.
                genotype = geno_df.iloc[i, :].round(0)
                (geno_alleles, _) = alleles_df.iloc[i, :]
                counts = genotype.value_counts()
                for x in [0.0, 1.0, 2.0]:
                    if x not in counts:
                        counts.loc[x] = 0

                zero_geno_count = (counts[0.0] * 2) + counts[1.0]
                two_geno_count = (counts[2.0] * 2) + counts[1.0]
                minor_allele_letter = geno_alleles[-1]
                if two_geno_count > zero_geno_count:
                    geno_flip = -1
                    minor_allele_letter = geno_alleles[0]

                # Check if the eQTL direction is measured towards the minor
                # allele.
                if row["AlleleAssessed"] != minor_allele_letter:
                    eqtl_flip = -1

                sum = zero_geno_count + two_geno_count
                minor_allele_frequency = min(zero_geno_count,
                                             two_geno_count) / sum

                eqtl_geno_info.append(
                    [geno_flip, eqtl_flip, sum / 2, minor_allele_letter, minor_allele_frequency])
            eqtl_geno_info_df = pd.DataFrame(eqtl_geno_info,
                                             columns=["geno_flip", "eqtl_flip",
                                                      "N", "MA", "MAF"])

            print("Geno info:")
            for label, col in zip(["genotype", "eQTL"],
                                  ["geno_flip", "eqtl_flip"]):
                flip_counts = eqtl_geno_info_df[col].value_counts()
                for x in [-1, 1]:
                    if x not in flip_counts:
                        flip_counts.loc[x] = 0
                print("\t{}:\tYES: {}\tNO: {}".format(label, flip_counts[-1],
                                                      flip_counts[1]))

            flip_df = eqtl_geno_info_df[["geno_flip", "eqtl_flip"]].copy()
            flip_df.drop_duplicates(inplace=True)
            n_total = eqtl_geno_info_df.shape[0]
            n_not_agreed = flip_df.shape[0]
            print("\tFlip concordance:\t{:.4f}% [{} not]".format(
                (((n_total - n_not_agreed) / n_total) * 100), n_not_agreed))

            n_maf_droppped = len(eqtl_geno_info_df.loc[
                                 eqtl_geno_info_df["MAF"] < self.maf_cutoff,
                                 :].index)
            print("\tMAF below cutoff:\t{}".format(n_maf_droppped))

            # Combine the data.
            df = eqtl_df.merge(eqtl_geno_info_df, left_index=True,
                               right_index=True)

            print("Working on column:")
            include = 3
            interesting_indices = set()
            combined_df = None
            for col in self.cols_of_interest:
                print("\t{}".format(col))

                working_df = df.copy()
                working_df["inter_zscore"] = inter_df.loc[col, :].values
                working_df["unflipped_inter_tvalue"] = tvalue_df.loc[col,
                                                       :].values
                working_df["inter_tvalue"] = working_df[
                                                 "unflipped_inter_tvalue"] * \
                                             working_df["geno_flip"]
                working_df["main_zscore"] = working_df["OverallZScore"] * \
                                            working_df["eqtl_flip"]
                working_df["cell"] = col

                working_df = working_df.loc[
                             (working_df.loc[:, "MAF"] > self.maf_cutoff) &
                             (working_df.loc[:,
                              "inter_zscore"] > z_score_cutoff), :]
                if col != "SEX":
                    working_df = working_df.loc[
                                 working_df.loc[:, "inter_tvalue"] > 0, :]

                if len(working_df.index) <= 0:
                    continue

                # print(working_df.loc[working_df["HGNCName"] == "HLA-DQA2", ["SNPName", "ProbeName", "HGNCName", "N", "MA", "MAF", "main_zscore", "inter_tvalue"]])
                #continue

                if combined_df is None:
                    working_df.reset_index(inplace=True, drop=False)
                    combined_df = working_df
                else:
                    working_df.reset_index(inplace=True, drop=False)
                    combined_df = pd.concat([combined_df, working_df], axis=0, ignore_index=True)

                # print(working_df.loc[working_df["HGNCName"] == "PIWIL2", ["eqtl_flip", "OverallZScore", "main_zscore", "geno_flip", "unflipped_inter_tvalue", "inter_tvalue", "inter_tvalue"]])

                positive_interaction = working_df.loc[
                                       working_df.loc[:, "inter_tvalue"] > 0,
                                       :].copy()
                negative_interaction = working_df.loc[
                                       working_df.loc[:, "inter_tvalue"] < 0,
                                       :].copy()

                for direction_df, filename, sort in zip(
                        [positive_interaction, negative_interaction],
                        ["pos_inter", "neg_inter"], [False, True]):
                    if len(df.index) <= 0:
                        continue

                    up_regulated = direction_df.loc[
                                   direction_df.loc[:, "main_zscore"] > 0,
                                   :].copy()
                    up_regulated.sort_values(by="inter_tvalue", ascending=sort,
                                             inplace=True)
                    if len(up_regulated.index) > 1:
                        up_filename = os.path.join(outdir,
                                                   "{}_{}_up_regulated.txt".format(
                                                       col, filename))
                        self.save(up_regulated, up_filename, self.max_url_len,
                                  z_score_cutoff)
                        interesting_indices.update(
                            up_regulated.head(include).index)

                    down_regulated = direction_df.loc[
                                     direction_df.loc[:, "main_zscore"] < 0,
                                     :].copy()
                    down_regulated.sort_values(by="inter_tvalue",
                                               ascending=sort, inplace=True)
                    if len(down_regulated.index) > 1:
                        down_filename = os.path.join(outdir,
                                                     "{}_{}_down_regulated.txt".format(
                                                         col, filename))
                        self.save(down_regulated, down_filename,
                                  self.max_url_len,
                                  z_score_cutoff)
                        interesting_indices.update(
                            down_regulated.head(include).index)

            interesting_indices = list(interesting_indices)
            interesting_indices.sort()
            interesting_indices_str = ' '.join(
                [str(x) for x in interesting_indices])
            print("Interesting indices:")
            print(interesting_indices_str)

            map_dict = {"CellMapNNLS_Astrocyte": "astro",
                        "CellMapNNLS_EndothelialCell": "endo",
                        "CellMapNNLS_Macrophage": "macro",
                        "CellMapNNLS_Neuron": "neuron",
                        "CellMapNNLS_Oligodendrocyte": "oligo"}

            # print(combined_df.loc[combined_df["HGNCName"] == "HLA-DQA2", ["index", "SNPName",  "ProbeName", "HGNCName", "N", "MA", "MAF", "main_zscore", "inter_tvalue"]])
            if combined_df is not None and len(combined_df.index) > 0:
                combined_df = combined_df.loc[combined_df["cell"] != "SEX", :]
                if len(combined_df.index) > 0:
                    combined_df["cell"] = combined_df["cell"].map(map_dict)
                    combined_df.sort_values(by="inter_tvalue", ascending=False, inplace=True)
                    combined_df.set_index("index", inplace=True)
                    self.save(combined_df, os.path.join(outdir, "all.txt"),
                              self.max_url_len, z_score_cutoff, include_cell=True)

    @staticmethod
    def save(df, outfile, max_url_len, z_score_cutoff, include_cell=False):
        print(
            "\tWriting output file: {}\tlen: {}".format(outfile, len(df.index)))
        with open(outfile, 'w') as f:
            f.write("Index\tSNPName\tProbeName\tHGNCName\tN\tMAF\teQTL\tInter")
            if include_cell:
                f.write("\tCell")
            f.write("\n")

            url_string = ""
            url_genes = []
            for i, row in df.iterrows():
                f.write("{}\t{}\t{}\t{}\t{}\t{:.2f}\t"
                        "{:.2f}\t{:.2f}".format(i,
                                                row["SNPName"],
                                                row["ProbeName"],
                                                row["HGNCName"],
                                                row["N"],
                                                row["MAF"],
                                                row["main_zscore"],
                                                row["inter_tvalue"]
                                                ))
                if include_cell:
                    f.write("\t{}".format(row["cell"]))
                f.write("\n")
                if (len(url_string) + len(row["ProbeName"])) < max_url_len:
                    if row["HGNCName"] not in url_genes:
                        url_string += row["ProbeName"]
                        url_genes.append(row["HGNCName"])

            f.write("\n")

            f.write("Z-score cutoff: {}\n".format(z_score_cutoff))
            f.write("N genes: {}\n".format(len(df.index)))
            f.write("\n")

            unique_snps = set(df["SNPName"])
            f.write("Number of unique SNPs: {}\n".format(len(unique_snps)))
            f.write("Unique SNPs:\n{}\n".format(', '.join(unique_snps)))
            f.write("\n")

            unique_genes = set(df["HGNCName"])
            f.write("Number of unique genes: {}\n".format(len(unique_genes)))
            f.write("Unique genes:\n{}\n".format(', '.join(unique_genes)))
            f.write("\n")

            f.write("Number of URL genes: {}\n".format(len(url_genes)))
            f.write("URL genes:\n{}\n".format(', '.join(url_genes)))
            f.write("\n")

        f.close()


if __name__ == '__main__':
    m = main()
    m.start()
