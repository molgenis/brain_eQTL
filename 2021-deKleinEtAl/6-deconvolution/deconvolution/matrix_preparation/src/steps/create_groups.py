"""
File:         create_groups.py
Created:      2020/03/12
Last Changed: 2020/03/20
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
import pickle
import os

# Third party imports.

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists, get_basename
from general.df_utilities import save_dataframe


class CreateGroups:
    def __init__(self, settings, eqtl_df, geno_df, alleles_df, expr_df, cov_df,
                 groups_file, force, outdir):
        """
        The initializer for the class.

        :param settings: string, the settings.
        :param eqtl_df: DataFrame, the eQTL probes data.
        :param geno_df: DataFrame, the genotype data.
        :param alleles_df: DataFrame, the alleles data.
        :param expr_df: DataFrame, the expression data.
        :param cov_df: DataFrame, the covariate data.
        :param groups_file: string, path to the groups file.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.eqtl_df = eqtl_df
        self.geno_df = geno_df
        self.alleles_df = alleles_df
        self.expr_df = expr_df
        self.cov_df = cov_df
        self.force = force

        # Load the groups.
        with open(groups_file, "rb") as f:
            groups_data = pickle.load(f)

        # Remove uninteresting groups.
        self.groups = self.filter_groups(groups_data,
                                         settings["min_eqtl_in_group"],
                                         settings["min_samples_in_group"])
        del groups_data

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_groups')
        prepare_output_dir(self.outdir)

    @staticmethod
    def filter_groups(groups_data, min_eqtl, min_samples):
        groups = {}
        incl_eqtls = 0
        total_eqtls = 0
        for group in groups_data:
            group_n_eqtls = group.get_n_eqtls()
            total_eqtls += group_n_eqtls
            if (group_n_eqtls >= min_eqtl) and \
                    (group.get_n_samples() >= min_samples):
                incl_eqtls += group_n_eqtls
                groups[group.get_id()] = group
        print("Filter completed, {}/{} [{:.2f}%] of the groups "
              "passed the threshold.".format(len(groups),
                                             len(groups_data),
                                             (100 / len(groups_data)) *
                                             len(groups)))
        print("Total eQTL's when combining groups: "
              "{}/{} [{:.2f}%].".format(incl_eqtls,
                                        total_eqtls,
                                        (100 / total_eqtls) * incl_eqtls))

        return groups

    def start(self):
        print("Creating groups.")
        for i, (group_id, group_obj) in enumerate(self.groups.items()):
            print("  Working on: {:10s} [{}/{} "
                  "{:.2f}%]".format(group_id, i + 1, len(self.groups),
                                    (100 / len(self.groups)) * (i + 1)))

            # Create the group dir.
            group_dir = os.path.join(self.outdir, group_id)
            prepare_output_dir(group_dir)

            # Define the output names.
            group_object = os.path.join(group_dir,
                                        "group.pkl")
            eqtl_outpath = os.path.join(group_dir,
                                        "eqtl_table.txt.gz")
            geno_outpath = os.path.join(group_dir,
                                        "genotype_table.txt.gz")
            alleles_outpath = os.path.join(group_dir,
                                           "genotype_alleles.txt.gz")
            expr_outpath = os.path.join(group_dir,
                                        "expression_table.txt.gz")
            cov_outpath = os.path.join(group_dir,
                                       "covariates_table.txt.gz")

            # Check if output file exist, if not, create it.
            if not check_file_exists(group_object) or self.force:
                with open(group_object, "wb") as f:
                    pickle.dump(group_obj, f)
                print("\tSaved group object: "
                      "{}".format(get_basename(group_object)))

            # Get the group indices.
            snp_mask = group_obj.get_snp_indices()
            sample_mask = group_obj.get_sample_indices()

            # Check if output file exist, if not, create it.
            if not check_file_exists(eqtl_outpath) or self.force:
                group_eqtl = self.eqtl_df.iloc[snp_mask, :].copy()
                save_dataframe(outpath=eqtl_outpath, df=group_eqtl,
                               index=False, header=True)
                del group_eqtl

            if not check_file_exists(geno_outpath) or self.force:
                group_geno = self.geno_df.iloc[snp_mask, sample_mask].copy()
                save_dataframe(outpath=geno_outpath, df=group_geno,
                               index=True, header=True)
                del group_geno

            if not check_file_exists(alleles_outpath) or self.force:
                group_alleles = self.alleles_df.iloc[snp_mask, :].copy()
                save_dataframe(outpath=alleles_outpath, df=group_alleles,
                               index=True, header=True)
                del group_alleles

            if not check_file_exists(expr_outpath) or self.force:
                group_expr = self.expr_df.iloc[snp_mask, sample_mask].copy()
                save_dataframe(outpath=expr_outpath, df=group_expr,
                               index=True, header=True)
                del group_expr

            if not check_file_exists(cov_outpath) or self.force:
                group_cov = self.cov_df.iloc[:, sample_mask].copy()
                save_dataframe(outpath=cov_outpath, df=group_cov,
                               index=True, header=True)
                del group_cov

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype matrix shape: {}".format(self.geno_df.shape))
        print("  > Alleles matrix shape: {}".format(self.alleles_df.shape))
        print("  > Expression matrix shape: {}".format(self.expr_df.shape))
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Force: {}".format(self.force))
        print("")
