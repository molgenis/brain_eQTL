"""
File:         create_cov_matrices.py
Created:      2020/03/12
Last Changed: 2020/06/03
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
from functools import reduce
import numpy as np
import pandas as pd

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists
from general.df_utilities import load_dataframe, save_dataframe


class CreateCovMatrix:
    def __init__(self, settings, marker_file, celltype_pcs, celltype_cs,
                 deconvolution, sample_order, force, outdir):
        """
        The initializer for the class.

        :param settings: string, the settings.
        :param marker_file: string, path to the marker file.
        :param celltype_pcs: DataFrame, the first component from PCA of each
                            celltype expression.
        :param celltype_cs: DataFrame, the first component from NMF of each
                            celltype expression.
        :param deconvolution: DataFrame, the estimated cell count proportions
                              of each celltype per sample.
        :param sample_order: list, order of samples.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.cov_file = settings["covariate_datafile"]
        self.tech_covs = settings["technical_covariates"]
        self.cohorts = settings["cohorts"]
        self.ref_cohort = settings["reference_cohort"]
        self.pheno_file = settings["phenotype_datafile"]
        self.eig_file = settings["eigenvectors_datafile"]
        self.n_eigen = settings["num_eigenvectors"]
        self.eig_bef_cov_corr_file = settings["eigenvectors_before_cov_corr_datafile"]
        self.marker_file = marker_file
        self.sample_order = sample_order
        self.celltype_pcs = celltype_pcs
        self.celltype_cs = celltype_cs
        self.deconvolution = deconvolution
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_cov_matrix')
        prepare_output_dir(self.outdir)
        self.outpath = os.path.join(self.outdir, "covariates_table.txt.gz")

        # Variables.
        self.covariates = None
        self.sex_dict = {"M": 0, "F": 1, np.nan: -1}

    def start(self):
        print("Starting creating covariate file.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.outpath) and not self.force:
            print("Skipping step, loading result.")
            self.covariates = load_dataframe(inpath=self.outpath, header=0,
                                             index_col=0)
        else:
            self.covariates = self.combine_files()
            self.save()

    def combine_files(self):
        # read the covariates file.
        print("Loading covariate matrix.")
        cov_df = load_dataframe(inpath=self.cov_file, header=0, index_col=0)
        tech_cov_df = cov_df[self.tech_covs].copy()
        cohorts_df = cov_df[self.cohorts].copy()
        del cov_df

        # validate the cohorts.
        print("Validating cohorts.")
        colsums = cohorts_df.sum(axis=1)
        cohorts_df[self.ref_cohort] = 0
        cohorts_df.loc[colsums == 0, self.ref_cohort] = 1
        if not cohorts_df.sum(axis=1).all():
            print("\tSome samples do not have a cohort.")
            exit()
        else:
            print("\tValid.")

        # read the phenotype file.
        print("Loading phenotype matrix.")
        pheno_df = load_dataframe(inpath=self.pheno_file, header=0,
                                  index_col=4, low_memory=False)

        # Combine the two gender columns, keep 'sex.by.expression' as main
        # gender ans use 'Gender' when no information is available.
        pheno_df = pheno_df.loc[:, ["Gender", "sex.by.expression"]]
        pheno_df.replace("no expression available", np.nan, inplace=True)
        pheno_df["SEX"] = pheno_df['sex.by.expression'].combine_first(
            pheno_df['Gender'])
        gender_df = pheno_df["SEX"].to_frame()
        del pheno_df
        gender_df = gender_df.replace({"SEX": self.sex_dict})

        # read the eigenvectors file.
        print("Loading eigenvectors matrix.")
        eigen_df = load_dataframe(self.eig_file, header=0, index_col=0)
        eigen_df = eigen_df.loc[:, ["Comp{}".format(x) for x in range(1, self.n_eigen + 1)]]

        # read the eigenvectors before covariate correction file.
        print("Loading eigenvectors before cov. correction matrix.")
        cov_cor_df = load_dataframe(self.eig_bef_cov_corr_file, header=0, index_col=0)
        cov_cor_df.columns = ["PC1-before-cov-correction", "PC2-before-cov-correction"]

        # read the marker genes expression file.
        print("Loading marker genes matrix.")
        marker_df = load_dataframe(self.marker_file, header=0, index_col=0)
        marker_df.sort_index(inplace=True)
        marker_df.drop_duplicates(inplace=True)
        marker_df = marker_df.T

        # merge.
        print("Merging matrices.")
        comb_cov = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True),
                          [tech_cov_df, cohorts_df, gender_df, eigen_df, cov_cor_df, marker_df, self.celltype_pcs.T, self.celltype_cs.T, self.deconvolution])
        comb_cov = comb_cov.T
        comb_cov = comb_cov[self.sample_order]
        comb_cov.index.name = "-"
        print("\tShape: {}".format(comb_cov.shape))

        # Remove old dataframes.
        del tech_cov_df, cohorts_df, gender_df, eigen_df, cov_cor_df, marker_df

        return comb_cov

    def save(self):
        save_dataframe(df=self.covariates, outpath=self.outpath,
                       index=True, header=True)

    def clear_variables(self):
        self.cov_file = None
        self.cov_exclude = None
        self.cohorts = None
        self.pheno_file = None
        self.eig_file = None
        self.n_eigen = None
        self.eig_bef_cov_corr_file = None
        self.marker_file = None
        self.sample_order = None
        self.force = None

    def get_covariates(self):
        return self.covariates

    def get_sex_dict(self):
        return self.sex_dict

    def get_outpath(self):
        return self.outpath

    def print_arguments(self):
        print("Arguments:")
        print("  > Covariates input file: {}".format(self.cov_file))
        print("  > Technical Covarates: {}".format(self.tech_covs))
        print("  > Cohorts: {}".format(self.cohorts))
        print("  > Reference cohort: {}".format(self.ref_cohort))
        print("  > Phenotype input file: {}".format(self.pheno_file))
        print("  > Eigenvectors input file: {}".format(self.eig_file))
        print("  > N. Eigenvectors: {}".format(self.n_eigen))
        print("  > Eigenvec before cov. corr. input file: {}".format(self.eig_bef_cov_corr_file))
        print("  > Markers input file: {}".format(self.marker_file))
        print("  > Celltype PCs: {}".format(self.celltype_pcs.shape))
        print("  > Deconvolution: {}".format(self.deconvolution.shape))
        print("  > Output path: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
