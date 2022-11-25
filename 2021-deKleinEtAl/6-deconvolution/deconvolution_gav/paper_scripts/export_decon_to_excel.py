#!/usr/bin/env python3

"""
File:         export_decon_to_excel.py
Created:      2021/02/16
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
import re
import os

# Third party imports.
import pandas as pd

# Local application imports.

"""
Syntax:
./export_decon_to_excel.py
"""


# Metadata
__program__ = "Export Decon to Excel"
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
        self.data = {"cortex": {"decon": "../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv",
                                "eqtl": "../matrix_preparation/cortex_eur_cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz",
                                "alleles": "../matrix_preparation/cortex_eur_cis/create_matrices/genotype_alleles.txt.gz"},
                     "cerebellum": {"decon": "../2020-11-20-decon-QTL/cis/cerebellum/decon_out/deconvolutionResults.csv",
                                    "eqtl": "../matrix_preparation/cerebellum_eur_cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz",
                                    "alleles": "../matrix_preparation/cerebellum_eur_cis/create_matrices/genotype_alleles.txt.gz"}}
        self.gene_info_path = "/groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz"

        # Set variables.
        outdir = os.path.join(str(Path(__file__).parent.parent), 'export_decon_to_excel')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        self.outfile = os.path.join(outdir, "Supplementary Table 7 - Deconvoluted eQTLs.xlsx")

    def start(self):
        self.print_arguments()

        print("Load gene info.")
        gene_info_df = self.load_file(self.gene_info_path)
        gene_dict = dict(zip(gene_info_df["ArrayAddress"], gene_info_df["Symbol"]))
        del gene_info_df

        with pd.ExcelWriter(self.outfile) as writer:
            for tissue, data in self.data.items():
                decon_path = data["decon"]
                eqtl_path = data["eqtl"]
                alleles_path = data["alleles"]
                print("Tissue: {}".format(tissue))
                print("  > Decon file: {}".format(decon_path))
                print("  > eQTL file: {}".format(eqtl_path))
                print("  > Alleles file: {}".format(alleles_path))

                # Load decon data.
                decon_df = self.load_file(decon_path)

                # Pre-process decon.
                probe_names = []
                snp_names = []
                for index in decon_df.index:
                    probe_names.append(index.split("_")[0])
                    snp_names.append("_".join(index.split("_")[1:]))
                decon_df["Gene"] = probe_names
                decon_df["SNP"] = snp_names

                decon_columns = []
                new_columns = []
                for col in decon_df.columns:
                    parts = re.split('_|:', col)
                    if "pvalue" in col:
                        new_columns.append("{} {}".format(parts[1], parts[2]))
                    elif "GT" in col:
                        new_columns.append("{} Beta:{}".format(parts[2], parts[3]))
                    elif "Beta" in col:
                        new_columns.append("{} Beta".format(parts[2]))
                    else:
                        new_columns.append(col)
                decon_df.columns = new_columns

                # Load alleles data.
                alleles_data = self.load_file(alleles_path, index_col=None)
                alleles_data.columns = ["SNP", "Alleles", "MinorAllele"]
                snp_to_alles_dict = dict(zip(alleles_data["SNP"], alleles_data["Alleles"]))
                del alleles_data

                # Load eQTL data.
                eqtl_df = self.load_file(eqtl_path, index_col=None)
                eqtl_df = eqtl_df.loc[:, ["ProbeName", "SNPName"]]
                eqtl_df["HGNCName"] = eqtl_df["ProbeName"].map(gene_dict)
                eqtl_df = eqtl_df.loc[:, ["ProbeName", "HGNCName", "SNPName"]]
                eqtl_df.columns = ["Gene", "Gene symbol", "SNP"]
                eqtl_df["Alleles"] = eqtl_df["SNP"].map(snp_to_alles_dict)
                eqtl_df["Allele Assessed"] = eqtl_df["Alleles"].str.split("/", n=1, expand=True)[1]

                # Merge.
                df = eqtl_df.merge(decon_df, on=["Gene", "SNP"])
                print(df)

                # Save
                df.to_excel(writer, sheet_name=tissue, na_rep="NA", index=False)
                print("Saving sheet {} with shape {}".format(tissue, df.shape))
                print("")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def print_arguments(self):
        print("Arguments:")
        print("  > Output file: {}".format(self.outfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
