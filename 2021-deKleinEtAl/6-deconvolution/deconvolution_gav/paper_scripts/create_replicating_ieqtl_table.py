#!/usr/bin/env python3

"""
File:         create_replicating_ieqtl_table.py
Created:      2021/01/27
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

# Third party imports.
import pandas as pd
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Create Replicating ieQTL Table"
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
        self.eqtl_type = getattr(arguments, 'type')

        # Declare input files.
        if self.eqtl_type == "cis":
            self.bulk_eqtl_infile = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/Iteration1/eQTLProbesFDR0.05-ProbeLevel.txt.gz"
        elif self.eqtl_type == "trans":
            self.bulk_eqtl_infile = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/trans/2020-05-26-Cortex-EUR-AFR-noENA-noPCA/Iteration1/eQTLs-crossMappingEQTLsRemoved-FDR0.05.txt.gz"
        self.bulk_alleles_infile = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_{}/create_matrices/genotype_alleles.txt.gz".format(self.eqtl_type)
        self.bulk_decon_infile = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/2020-11-20-decon-QTL/{}/cortex/decon_out/deconvolutionResults.csv".format(self.eqtl_type)
        self.sn_eqtl_infolder = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/{}_100Perm/".format(self.eqtl_type)
        self.sn_eqtl_filename = "eQTLsFDR-ProbeLevel.txt.gz"
        self.eqtl_ea_col = "AlleleAssessed"

        self.cell_types = [("AST", "Astrocyte"),
                           ("END", "EndothelialCell"),
                           ("MIC", "Microglia"),
                           ("EX", "ExNeuron"),
                           ("IN", "InNeuron"),
                           ("OLI", "Oligodendrocyte")]

        # Set variables.
        outdir = os.path.join(str(Path(__file__).parent.parent), 'create_replicating_ieqtl_table')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # Define the output file.
        self.outpath = os.path.join(outdir, "replicating_{}_ieqtl_table.xlsx".format(self.eqtl_type))

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
        parser.add_argument("-type",
                            type=str,
                            required=False,
                            choices=["cis", "trans"],
                            default="cis",
                            help="The type of eQTLs to plot.")

        return parser.parse_args()

    def start(self):
        print("### Loading bulk data ###")
        bulk_eqtl_df = pd.read_csv(self.bulk_eqtl_infile,
                                   sep="\t",
                                   header=0,
                                   index_col=None)
        print("\tBulk eQTL data frame: {}".format(bulk_eqtl_df.shape))
        print(bulk_eqtl_df)
        tmp = bulk_eqtl_df[["ProbeName", "SNPName"]].copy()
        tmp["key"] = tmp["ProbeName"] + "_" + tmp["SNPName"]
        if len(tmp["key"].unique()) != tmp.shape[0]:
            print("Error: eqtl keys are not unique!")
            exit()

        bulk_decon_df = pd.read_csv(self.bulk_decon_infile,
                                    sep="\t",
                                    header=0,
                                    index_col=0)
        print("\tBulk decon data frame: {}".format(bulk_eqtl_df.shape))
        print(bulk_decon_df)

        print("### Pre-processing bulk decon data ###")
        buk_decon_fdr_beta_df = self.preprocess_decon_df(bulk_decon_df)
        print(buk_decon_fdr_beta_df)

        print("### Merging bulk data ###")
        bulk_df = bulk_eqtl_df.merge(buk_decon_fdr_beta_df, on=["ProbeName", "SNPName"])
        print(bulk_df)

        # Add alleles.
        alleles_data = pd.read_csv(self.bulk_alleles_infile, sep="\t", header=0)
        alleles_data.columns = ["SNP", "Alleles", "MinorAllele"]
        snp_to_alles_dict = dict(zip(alleles_data["SNP"], alleles_data["Alleles"]))
        bulk_df["DeconAlleles"] = bulk_df["SNPName"].map(snp_to_alles_dict)
        bulk_df["DeconQTLAlleleAssessed"] = bulk_df["DeconAlleles"].str.split("/", n=1, expand=True)[1]
        bulk_df.drop(["DeconAlleles"], axis=1, inplace=True)

        print("### Flipping beta's ###")
        bulk_df["deconQTLFlip"] = bulk_df[self.eqtl_ea_col] != bulk_df["DeconQTLAlleleAssessed"]
        for col in bulk_df:
            if col.endswith("_beta"):
                bulk_df.loc[:, col] = bulk_df.loc[:, col] * bulk_df["deconQTLFlip"].map({True: -1, False: 1})
        bulk_df.drop(["deconQTLFlip", "DeconQTLAlleleAssessed"], axis=1, inplace=True)

        print("### Drop column with one distinct value ###")
        for col in bulk_df.columns:
            values = bulk_df[col].unique()
            if len(values) == 1:
                bulk_df.drop(col, inplace=True, axis=1)
                print("\tColumn '{}' has one distinct value: '{}'".format(col, values[0]))
        print(bulk_df)

        print("### Adding single-nucleus data ###")
        for col_index, (abbreviation, full_name) in enumerate(self.cell_types):
            # Load the data.
            sn_df = pd.read_csv(os.path.join(self.sn_eqtl_infolder, abbreviation, self.sn_eqtl_filename),
                                sep="\t",
                                header=0,
                                index_col=None)
            print("\tSingle-nucleus {} eQTL data frame: {}".format(full_name, sn_df.shape))

            # Select the columns of interest
            sn_subset = sn_df[["ProbeName", "SNPName", self.eqtl_ea_col, "FDR", "OverallZScore"]]
            sn_subset.columns = ["ProbeName", "SNPName", "sn_{}_ea_col".format(full_name), "sn_{}_eQTL_FDR".format(full_name), "sn_{}_eqtl_OverallZscore".format(full_name)]

            # Merge.
            bulk_df = bulk_df.merge(sn_subset, on=["ProbeName", "SNPName"], how="left")

        print(bulk_df)
        print("### Flipping z-score's ###")
        for _, full_name in self.cell_types:
            bulk_df["sn_{}_flip".format(full_name)] = bulk_df[self.eqtl_ea_col] != bulk_df["sn_{}_ea_col".format(full_name)]

            bulk_df.loc[:, "sn_{}_eqtl_OverallZscore".format(full_name)] = bulk_df.loc[:, "sn_{}_eqtl_OverallZscore".format(full_name)] * bulk_df["sn_{}_flip".format(full_name)].map({True: -1, False: 1})

            bulk_df.drop(["sn_{}_ea_col".format(full_name), "sn_{}_flip".format(full_name)], axis=1, inplace=True)
        print(bulk_df.columns)

        if bulk_df["PValue"].min() < 2.2250738585072014e-308:
            bulk_df["PValue"] = bulk_df["PValue"].astype(str)
        print(bulk_df.dtypes)

        bulk_df.columns = [x.replace("_", " ") for x in bulk_df.columns]

        print("### Saving data frame ###")
        bulk_df.to_excel(self.outpath, sheet_name="Bulk ieQTLs replication in SN", na_rep="NA", index=False)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(self.outpath),
                                      bulk_df.shape))

    def preprocess_decon_df(self, df):
        # Extract the p_values and convert to fdr.
        tmp = df.copy()
        data = []
        indices = []
        for col in tmp.columns:
            if col.endswith("_pvalue"):
                data.append(multitest.multipletests(tmp.loc[:, col], method='fdr_bh')[1])
                indices.append("bulk_" + col.replace("CellMapNNLS_", "").replace("_pvalue", "") + "_interaction_FDR")
            elif col.endswith(":GT"):
                data.append(tmp.loc[:, col])
                indices.append("bulk_" + col.split("_")[2].split(":")[0] + "_interaction_beta")
        fdr_beta_df = pd.DataFrame(data, index=indices, columns=tmp.index).T
        fdr_beta_df = fdr_beta_df.reindex(sorted(fdr_beta_df.columns), axis=1)

        del tmp

        # Split the index.
        probe_names = []
        snp_names = []
        for index in fdr_beta_df.index:
            probe_names.append(index.split("_")[0])
            snp_names.append("_".join(index.split("_")[1:]))
        fdr_beta_df.insert(0, "ProbeName", probe_names)
        fdr_beta_df.insert(1, "SNPName", snp_names)

        return fdr_beta_df


if __name__ == '__main__':
    m = main()
    m.start()
