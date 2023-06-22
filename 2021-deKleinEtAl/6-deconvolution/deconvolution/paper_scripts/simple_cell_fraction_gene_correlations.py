#!/usr/bin/env python3

"""
File:         simple_cell_fraction_gene_correlations.py
Created:      2021/02/02
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
import glob
import os

# Third party imports.
import pandas as pd
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Simple Cell Fractions Gene Correlations"
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

"""
Syntax:
./simple_cell_fraction_gene_correlations.py -cf ../matrix_preparation/cortex_eur_trans/perform_deconvolution/deconvolution_table.txt -ge ../matrix_preparation/cortex_eur_trans/create_matrices/expression_table.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/trans/2020-05-26-Cortex-EUR-AFR-noENA-noPCA/gwasupdate/Iteration1-dsZscores -op cortex_eur_trans_
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.ge_path = getattr(arguments, 'gene_expression')
        self.gte_folder = getattr(arguments, 'gte_folder')
        self.outfile_prefix = getattr(arguments, 'outfile_prefix')
        self.cohort_info = ["GTE-AFR-AMPAD-MSBB-V2",
                            "GTE-EUR-AMPAD-MAYO-V2",
                            "GTE-EUR-AMPAD-MSBB-V2",
                            "GTE-EUR-AMPAD-ROSMAP-V2",
                            "GTE-EUR-BrainGVEX-V2",
                            "GTE-AFR-CMC",
                            "GTE-EUR-CMC",
                            "GTE-AFR-CMC_HBCC_set1",
                            "GTE-AFR-CMC_HBCC_set2",
                            "GTE-EUR-CMC_HBCC_set2",
                            "GTE-EUR-CMC_HBCC_set3",
                            "GTE-EUR-GTEx",
                            "GTE-EUR-GVEX",
                            "GTE-AFR-LIBD_1M",
                            "GTE-EUR-LIBD_1M",
                            "GTE-AFR-LIBD_h650",
                            "GTE-EUR-LIBD_h650",
                            "GTE-EUR-NABEC-H550",
                            "GTE-EUR-NABEC-H610",
                            "GTE-EUR-UCLA_ASD",
                            "GTE-EUR-TargetALS"
                            ]

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   "simple_gene_cellfraction_corr")

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
                            help="show program's version number and exit")
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fracttions matrix")
        parser.add_argument("-ge",
                            "--gene_expression",
                            type=str,
                            required=True,
                            help="The path to the gene expression matrix")
        parser.add_argument("-gte",
                            "--gte_folder",
                            type=str,
                            required=True,
                            help="The path to the folder containg 'GTE-*' files.")
        parser.add_argument("-op",
                            "--outfile_prefix",
                            type=str,
                            required=False,
                            default="",
                            help="A prefix for the output files.")

        return parser.parse_args()

    def start(self):
        print("Loading data.")
        cf_df = self.load_file(self.cf_path, low_memory=False)
        ge_df = self.load_file(self.ge_path)

        print("Pre-process")
        gte_combined = None
        for fpath in glob.glob(os.path.join(self.gte_folder, 'GTE-*')):
            gte_filename = os.path.basename(fpath).split(".")[0]
            if gte_filename not in self.cohort_info:
                print("GTE file '{}' will be skipped.".format(gte_filename))
                continue
            gte_df = self.load_file(fpath, index_col=None, header=None)
            gte_df.columns = ["gene_id", "expr_id"]
            gte_df["cohort"] = gte_filename
            if gte_combined is None:
                gte_combined = gte_df
            else:
                gte_combined = pd.concat([gte_combined, gte_df], axis=0)
        gte_combined.set_index("expr_id", inplace=True)

        # Subset.
        cf_df = cf_df.loc[gte_combined.index, :]
        ge_df = ge_df.loc[:, gte_combined.index]

        if not cf_df.index.equals(ge_df.columns):
            print("The expressiom file columns do not match the cell "
                  "type fractions file.")
            exit()
        cf_df = cf_df.T
        ge_df = ge_df.drop_duplicates()

        print(cf_df)
        print(ge_df)

        print("Correlate data.")
        corr_coef_df, _ = self.correlate_dataframes(ge_df, cf_df)
        self.save_dataframe(corr_coef_df, "all_datasets_coefficients")

        print("Correlate data per dataset.")
        for cohort in gte_combined["cohort"].unique():
            print("Processing dataset '{}'".format(cohort))
            cohort_gte = gte_combined.loc[gte_combined["cohort"] == cohort, :]
            cohort_cf = cf_df.loc[:, cohort_gte.index].copy()
            cohort_ge = ge_df.loc[:, cohort_gte.index].copy()

            cohort_corr_coef_df, _ = self.correlate_dataframes(cohort_ge,
                                                               cohort_cf)
            self.save_dataframe(cohort_corr_coef_df,
                                "{}_dataset_coefficients".format(cohort))

    @staticmethod
    def correlate_dataframes(df1, df2):
        all_coefficients = []
        all_pvalues = []
        for _, row1 in df1.iterrows():
            coefficients = []
            pvalues = []
            for i, (_, row2) in enumerate(df2.iterrows()):
                coef, p = stats.spearmanr(row1, row2)
                coefficients.append(coef)
                pvalues.append(p)
            all_coefficients.append(coefficients)
            all_pvalues.append(pvalues)

        corr_coef_df = pd.DataFrame(all_coefficients,
                                    index=df1.index,
                                    columns=df2.index)

        pvalue_df = pd.DataFrame(all_pvalues,
                                 index=df1.index,
                                 columns=df2.index)

        return corr_coef_df, pvalue_df

    def save_dataframe(self, df, name):
        df.to_csv(os.path.join(self.outdir,
                               '{}{}.txt.gz'.format(self.outfile_prefix, name)),
                  sep="\t", index=True, header=True, compression="gzip")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  low_memory=True):
        if path.endswith(".pkl"):
            df = pd.read_pickle(path)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df


if __name__ == '__main__':
    m = main()
    m.start()
