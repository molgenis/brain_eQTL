"""
File:         mask_matrices.py
Created:      2020/03/12
Last Changed: 2020/03/13
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
import os

# Third party imports.
import pandas as pd

# Local application imports.
from src.utilities import prepare_output_dir, check_file_exists
from src.df_utilities import save_dataframe


class MaskMatrices:
    def __init__(self, settings, geno_df, alleles_df, expr_df, cov_df,
                 force, outdir):
        """
        The initializer for the class.

        :param settings: string, the settings.
        :param geno_df: DataFrame, the genotype data.
        :param alleles_df: DataFrame, the alleles data.
        :param expr_df: DataFrame, the expression data.
        :param cov_df: DataFrame, the covariate data.
        :param marker_file: string, path to the marker file.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.geno_df = geno_df
        self.alleles_df = alleles_df
        self.expr_df = expr_df
        self.cov_df = cov_df
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'mask_matrices')
        prepare_output_dir(self.outdir)

        # Define the output names.
        self.eqtl_translate_outpath = os.path.join(self.outdir,
                                                   "eqtl_translate_table.txt.gz")
        self.sample_translate_outpath = os.path.join(self.outdir,
                                                     "sample_translate_table.txt.gz")
        self.cov_translate_outpath = os.path.join(self.outdir,
                                                  "cov_translate_table.txt.gz")
        self.geno_outpath = os.path.join(self.outdir, "genotype_table.txt.gz")
        self.alleles_outpath = os.path.join(self.outdir,
                                            "genotype_alleles.txt.gz")
        self.expr_outpath = os.path.join(self.outdir, "expression_table.txt.gz")
        self.cov_outpath = os.path.join(self.outdir, "covariates-cortex.txt.gz")

    def start(self):
        print("Starting creating masked files.")
        self.print_arguments()

        # Get the sizes.
        (n_eqtls, n_samples) = self.geno_df.shape
        n_covs = self.cov_df.shape[0]

        # Create masks.
        eqtl_mask = ["eqtl_" + str(x) for x in range(n_eqtls)]
        sample_mask = ["sample_" + str(x) for x in range(n_samples)]
        cov_mask = ["cov_" + str(x) for x in range(n_covs)]

        # Create translate dicts.
        print("Creating translation files.")
        if not check_file_exists(self.eqtl_translate_outpath) or self.force:
            eqtl_translate = pd.DataFrame({'unmasked': list(self.geno_df.index),
                                           'masked': eqtl_mask})
            save_dataframe(outpath=self.eqtl_translate_outpath,
                           df=eqtl_translate,
                           index=False, header=True)
            del eqtl_translate

        if not check_file_exists(self.sample_translate_outpath) or self.force:
            sample_translate = pd.DataFrame(
                {'unmasked': list(self.geno_df.columns),
                 'masked': sample_mask})
            save_dataframe(outpath=self.sample_translate_outpath,
                           df=sample_translate,
                           index=False, header=True)
            del sample_translate

        if not check_file_exists(self.cov_translate_outpath) or self.force:
            cov_translate = pd.DataFrame({'unmasked': list(self.cov_df.index),
                                          'masked': cov_mask})
            save_dataframe(outpath=self.cov_translate_outpath, df=cov_translate,
                           index=False, header=True)
            del cov_translate

        # Start masking the dataframes.
        print("Start masking files.")
        if not check_file_exists(self.geno_outpath) or self.force:
            self.geno_df.index = eqtl_mask
            self.geno_df.columns = sample_mask
            save_dataframe(outpath=self.geno_outpath, df=self.geno_df,
                           index=True, header=True)

        if not check_file_exists(self.alleles_outpath) or self.force:
            self.alleles_df.index = eqtl_mask
            save_dataframe(outpath=self.alleles_outpath, df=self.alleles_df,
                           index=True, header=True)

        if not check_file_exists(self.expr_outpath) or self.force:
            self.expr_df.index = eqtl_mask
            self.expr_df.columns = sample_mask
            save_dataframe(outpath=self.expr_outpath, df=self.expr_df,
                           index=True, header=True)

        if not check_file_exists(self.cov_outpath) or self.force:
            self.cov_df.index = cov_mask
            self.cov_df.columns = sample_mask
            save_dataframe(outpath=self.cov_outpath, df=self.cov_df,
                           index=True, header=True)

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype matrix shape: {}".format(self.geno_df.shape))
        print("  > Alleles matrix shape: {}".format(self.alleles_df.shape))
        print("  > Expression matrix shape: {}".format(self.expr_df.shape))
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Force: {}".format(self.force))
        print("")
