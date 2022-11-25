"""
File:         dataset.py
Created:      2020/03/16
Last Changed: 2020/06/19
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
import itertools
import os

# Third party imports.
import pandas as pd
import seaborn as sns
import scipy.stats as stats

# Local application imports.
from general.df_utilities import load_dataframe


class Dataset:
    def __init__(self, name, settings, alpha=0.05, nrows=None, interest=None):
        self.input_dir = os.path.join(settings.get_setting("input_dir"), name)
        filenames = settings.get_setting("filenames")
        self.eqtl_filename = filenames["eqtl"]
        self.geno_filename = filenames["genotype"]
        self.alleles_filename = filenames["alleles"]
        self.expr_filename = filenames["expression"]
        self.cov_filename = filenames["covariates"]
        self.markers_filename = filenames["markers"]

        self.inter_input_dir = os.path.join(settings.get_setting("interaction_input_dir"), name)
        inter_subdirs = settings.get_setting("interaction_input_subfolders")
        self.inter_cov_subdir = inter_subdirs["covariates_of_interest"]
        self.inter_tech_cov_subdir = inter_subdirs["technical_covariates"]
        inter_filenames = settings.get_setting("interaction_filenames")
        self.pvalue_filename = inter_filenames["pvalues"]
        self.zscore_filename = inter_filenames["zscores"]
        self.snp_tvalue_filename = inter_filenames["snp_tvalues"]
        self.inter_tvalue_filename = inter_filenames["inter_tvalues"]

        self.groups = settings.get_setting("groups")
        self.celltypes = settings.get_setting("celltypes")
        self.colormap = settings.get_setting("colormap")
        self.cellmap_methods = settings.get_setting("cellmap_method_prefix_and_suffix")
        self.marker_genes = settings.get_setting("marker_genes_prefix")
        self.signif_cutoff = stats.norm.isf(alpha)
        self.interest = interest
        if nrows == -1:
            nrows = None
        if nrows is not None and nrows <= 0:
            print("Unexpected argument for -n / --n_eqtls: '{}'".format(nrows))
            exit()
        if self.interest is not None:
            nrows = max(self.interest) + 1
        self.nrows = nrows

        # Declare empty variables.
        self.eqtl_df = None
        self.geno_df = None
        self.alleles_df = None
        self.expr_df = None
        self.cov_df = None
        self.inter_cov_pvalue_df = None
        self.inter_tech_cov_pvalue_df = None
        self.inter_cov_zscore_df = None
        self.inter_tech_cov_zscore_df = None
        self.inter_cov_snp_tvalue_df = None
        self.inter_tech_cov_snp_tvalue_df = None
        self.inter_cov_inter_tvalue_df = None
        self.inter_tech_cov_inter_tvalue_df = None
        self.eqtl_and_interactions_df = None
        self.marker_df = None

    def load_all(self):
        print("Loading all dataframes for validation.")
        self.nrows = None
        self.get_celltypes()
        self.get_eqtl_df()
        self.get_geno_df()
        self.get_alleles_df()
        self.get_expr_df()
        self.get_cov_df()
        self.get_inter_cov_pvalue_df()
        self.get_inter_tech_cov_pvalue_df()
        self.get_inter_cov_zscore_df()
        self.get_inter_tech_cov_zscore_df()
        self.get_inter_cov_snp_tvalue_df()
        self.get_inter_tech_cov_snp_tvalue_df()
        self.get_inter_cov_inter_tvalue_df()
        self.get_inter_tech_cov_inter_tvalue_df()
        self.get_eqtl_and_interactions_df()
        self.get_marker_df()
        print("Validation finished.")
        exit()

    def get_groups(self):
        return self.groups

    def get_celltypes(self):
        return self.celltypes

    def get_colormap(self):
        return self.colormap

    def get_diverging_cmap(self):
        return sns.diverging_palette(self.colormap["low"], self.colormap["high"], as_cmap=True)

    def get_cellmap_methods(self):
        return self.cellmap_methods

    def get_marker_genes(self):
        return self.marker_genes

    def get_significance_cutoff(self):
        return self.signif_cutoff

    def get_eqtl_df(self):
        if self.eqtl_df is None:
            eqtl_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                         self.eqtl_filename),
                                     header=0,
                                     index_col=False,
                                     nrows=self.nrows)
            if self.interest is not None:
                eqtl_df = eqtl_df.iloc[self.interest, :]
            self.eqtl_df = eqtl_df

            self.validate()
        return self.eqtl_df

    def get_geno_df(self):
        if self.geno_df is None:
            geno_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                         self.geno_filename),
                                     header=0,
                                     index_col=0,
                                     nrows=self.nrows)
            if self.interest is not None:
                geno_df = geno_df.iloc[self.interest, :]
            self.geno_df = geno_df

            self.validate()
        return self.geno_df

    def get_alleles_df(self):
        if self.alleles_df is None:
            alleles_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                            self.alleles_filename),
                                        header=0,
                                        index_col=0,
                                        nrows=self.nrows)
            if self.interest is not None:
                alleles_df = alleles_df.iloc[self.interest, :]
            self.alleles_df = alleles_df

            self.validate()
        return self.alleles_df

    def get_expr_df(self):
        if self.expr_df is None:
            expr_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                         self.expr_filename),
                                     header=0,
                                     index_col=0,
                                     nrows=self.nrows)

            if self.interest is not None:
                expr_df = expr_df.iloc[self.interest, :]
            self.expr_df = expr_df

            self.validate()
        return self.expr_df

    def get_cov_df(self):
        if self.cov_df is None:
            self.cov_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                             self.cov_filename),
                                         header=0,
                                         index_col=0)
            self.validate()
        return self.cov_df

    def get_inter_cov_pvalue_df(self):
        if self.inter_cov_pvalue_df is None:
            inter_cov_pvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_cov_subdir,
                                    self.pvalue_filename),
                header=0,
                index_col=0)
            if self.interest is not None:
                inter_cov_pvalue_df = inter_cov_pvalue_df.iloc[:, self.interest]
            self.inter_cov_pvalue_df = inter_cov_pvalue_df

            self.validate()
        return self.inter_cov_pvalue_df

    def get_inter_tech_cov_pvalue_df(self):
        if self.inter_tech_cov_pvalue_df is None:
            inter_tech_cov_pvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_tech_cov_subdir,
                                    self.pvalue_filename),
                header=0,
                index_col=0)
            if self.interest is not None:
                inter_tech_cov_pvalue_df = inter_tech_cov_pvalue_df.iloc[:, self.interest]
            self.inter_tech_cov_pvalue_df = inter_tech_cov_pvalue_df

            self.validate()
        return self.inter_tech_cov_pvalue_df

    def get_inter_cov_zscore_df(self):
        if self.inter_cov_zscore_df is None:
            inter_cov_zscore_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_cov_subdir,
                                    self.zscore_filename),
                header=0,
                index_col=0)
            if self.interest is not None:
                inter_cov_zscore_df = inter_cov_zscore_df.iloc[:,
                                           self.interest]
            self.inter_cov_zscore_df = inter_cov_zscore_df

            self.validate()
        return self.inter_cov_zscore_df

    def get_inter_tech_cov_zscore_df(self):
        if self.inter_tech_cov_zscore_df is None:
            inter_tech_cov_zscore_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_tech_cov_subdir,
                                    self.zscore_filename),
                header=0,
                index_col=0)
            if self.interest is not None:
                inter_tech_cov_zscore_df = inter_tech_cov_zscore_df.iloc[:,
                                           self.interest]
            self.inter_tech_cov_zscore_df = inter_tech_cov_zscore_df

            self.validate()
        return self.inter_tech_cov_zscore_df

    def get_inter_cov_snp_tvalue_df(self):
        if self.inter_cov_snp_tvalue_df is None:
            inter_cov_snp_tvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_cov_subdir,
                                    self.snp_tvalue_filename),
                header=0,
                index_col=0)
            if self.interest is not None:
                inter_cov_snp_tvalue_df = inter_cov_snp_tvalue_df.iloc[:, self.interest]
            self.inter_cov_snp_tvalue_df = inter_cov_snp_tvalue_df

            self.validate()
        return self.inter_cov_snp_tvalue_df

    def get_inter_tech_cov_snp_tvalue_df(self):
        if self.inter_tech_cov_snp_tvalue_df is None:
            inter_tech_cov_snp_tvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_tech_cov_subdir,
                                    self.snp_tvalue_filename),
                header=0,
                index_col=0)
            if self.interest is not None:
                inter_tech_cov_snp_tvalue_df = inter_tech_cov_snp_tvalue_df.iloc[:, self.interest]
            self.inter_tech_cov_snp_tvalue_df = inter_tech_cov_snp_tvalue_df

            self.validate()
        return self.inter_tech_cov_snp_tvalue_df

    def get_inter_cov_inter_tvalue_df(self):
        if self.inter_cov_inter_tvalue_df is None:
            inter_cov_inter_tvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_cov_subdir,
                                    self.inter_tvalue_filename),
                header=0,
                index_col=0)
            if self.interest is not None:
                inter_cov_inter_tvalue_df = inter_cov_inter_tvalue_df.iloc[:, self.interest]
            self.inter_cov_inter_tvalue_df = inter_cov_inter_tvalue_df

            self.validate()
        return self.inter_cov_inter_tvalue_df

    def get_inter_tech_cov_inter_tvalue_df(self):
        if self.inter_tech_cov_inter_tvalue_df is None:
            inter_tech_cov_inter_tvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_tech_cov_subdir,
                                    self.inter_tvalue_filename),
                header=0,
                index_col=0)
            if self.interest is not None:
                inter_tech_cov_inter_tvalue_df = inter_tech_cov_inter_tvalue_df.iloc[:, self.interest]
            self.inter_tech_cov_inter_tvalue_df = inter_tech_cov_inter_tvalue_df

            self.validate()
        return self.inter_tech_cov_inter_tvalue_df

    def get_eqtl_and_interactions_df(self):
        # Get the complete input dataframes.
        df1 = load_dataframe(inpath=os.path.join(self.input_dir,
                                                 self.eqtl_filename),
                             header=0, index_col=False)
        df2 = load_dataframe(inpath=os.path.join(self.inter_input_dir,
                                                 self.inter_cov_subdir,
                                                 self.zscore_filename),
                             header=0, index_col=0).T

        # Check if the files math up.
        if df1.shape[0] != df2.shape[0]:
            print("Input files do not match (1).")
            exit()
        for i in range(df1.shape[0]):
            if not df2.index[i].startswith(df1["SNPName"][i]):
                print("Input files do not match (2).")
                exit()

        # Reset the indices.
        df1.reset_index(drop=True, inplace=True)
        df2.reset_index(drop=True, inplace=True)

        # Replace the z-scores with 1's and 0's (significant vs not-siginifcant)
        df2[df2 <= self.signif_cutoff] = 0
        df2[df2 > self.signif_cutoff] = 1
        df2 = df2.fillna(0).astype('int8')

        # Combine.
        self.eqtl_and_interactions_df = pd.concat([df1, df2], axis=1)
        self.eqtl_and_interactions_df.index = self.eqtl_and_interactions_df.index.astype(str) + "_" + self.eqtl_and_interactions_df["SNPName"] + "_" + self.eqtl_and_interactions_df["ProbeName"]
        del df1, df2

        return self.eqtl_and_interactions_df

    def get_marker_df(self):
        if self.marker_df is None:
            self.marker_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                                self.markers_filename),
                                            header=0,
                                            index_col=False)
            self.validate()
        return self.marker_df

    def validate(self):
        if self.eqtl_df is not None:
            if self.geno_df is not None:
                reference = self.eqtl_df["SNPName"].copy()
                reference.rename("-", inplace=True)
                alternative = pd.Series(self.geno_df.index, index=reference.index, name="-")
                if not alternative.equals(reference):
                    print("Order of SNPs in eqtl_df and geno_df "
                          "are not identical.")
                    exit()
            if self.expr_df is not None:
                reference = self.eqtl_df["ProbeName"].copy()
                reference.rename("-", inplace=True)
                alternative = pd.Series(self.expr_df.index, index=reference.index, name="-")
                if not alternative.equals(reference):
                    print("Order of Probes in eqtl_df and expr_df "
                          "are not identical.")

        if self.geno_df is not None and self.alleles_df is not None:
            if not self.geno_df.index.identical(self.alleles_df.index):
                print("Order of SNPs in geno_df and alleles_df are not "
                      "identical.")

        dfs = [self.geno_df, self.expr_df, self.cov_df]
        for (a, b) in list(itertools.combinations(dfs, 2)):
            if a is not None and b is not None and \
                    not a.columns.identical(b.columns):
                print("Order of samples are not identical.")
                exit()

        print("\tValid.")
