"""
File:         perform_celltype_factorization.py
Created:      2020/04/07
Last Changed: 2020/07/15
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
import numpy as np
from sklearn.decomposition import PCA, NMF

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists
from general.df_utilities import load_dataframe, save_dataframe


class PerformCelltypeFactorization:
    def __init__(self, settings, profile_file, profile_df, ct_expr_file,
                 force, outdir):
        """
        The initializer for the class.


        :param settings: string, the settings.
        :param profile_file: string, the datafile contaioning the celltype
                             profile.
        :param profile_df: DataFrame, the celltype profile.
        :param ct_expr_file: string, the datafile containing expression
                             of the celltype profiles.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.profile_file = profile_file
        self.profile_df = profile_df
        self.ct_expr_file = ct_expr_file
        self.force = force

        # Prepare an output directory.
        self.outdir = os.path.join(outdir, 'perform_celltype_factorization')
        prepare_output_dir(self.outdir)
        self.pca_outpath = os.path.join(self.outdir, "celltype_pca.txt.gz")
        self.nmf_outpath = os.path.join(self.outdir, "celltype_nmf.txt.gz")

        # Create empty variables.
        self.celltype_expression = None
        self.celltype_pcs = None
        self.celltype_cs = None

    def start(self):
        print("Starting factorization of celltype profile expression.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.pca_outpath) and check_file_exists(self.nmf_outpath) and not self.force:
            print("Skipping step, loading result.")
            self.celltype_pcs = load_dataframe(inpath=self.pca_outpath,
                                               header=0, index_col=0)
            self.celltype_cs = load_dataframe(inpath=self.nmf_outpath,
                                              header=0, index_col=0)
        else:
            self.celltype_expression, self.celltype_pcs, self.celltype_cs = self.perform_matrix_factorization()
            self.save()

    def perform_matrix_factorization(self):
        # Load the expression data.
        print("Loading celltype expression data.")
        ct_expr_df = load_dataframe(inpath=self.ct_expr_file,
                                    header=0, index_col=0)

        if self.profile_df is None:
            # Load the celltype profile file.
            print("Loading cell type profile matrix.")
            self.profile_df = load_dataframe(self.profile_file,
                                             header=0, index_col=0)

        # Find the genes specific to each celltype.
        gene_celltypes = self.normalize(self.profile_df).idxmax(axis=1)

        # Construct a dataframe of the first component of each celltype
        # subset expression profile.
        pca_data = []
        print("Performing PCA")
        for celltype in self.profile_df.columns:
            print("\tWorking on: {}".format(celltype))
            ct_genes = gene_celltypes[gene_celltypes == celltype].index
            ct_expr = ct_expr_df.loc[ct_expr_df.index.isin(ct_genes), :]
            print("\t  N = {}".format(len(ct_expr.index)))

            # perform PCA over the expression of these genes.
            print("\t  PCA")
            pca_component = self.get_first_pca_component(ct_expr)
            pca_component_values = [x[0] for x in list(pca_component)]
            pca_data.append(pca_component_values)

        # Create the data frame.
        celltype_pcs = pd.DataFrame(pca_data,
                                    index=["{}PCA_{}_PC1".format(*x.split("_")) for x in self.profile_df.columns],
                                    columns=ct_expr_df.columns)

        # Shift the expression to be all positive.
        shifted_ct_expr = ct_expr_df.copy()
        if ct_expr_df.values.min() < 0:
            shifted_ct_expr = self.perform_shift(ct_expr_df)

        # Construct a dataframe of the first component of each celltype
        # subset expression profile.
        nmf_data = []
        print("Performing NMF")
        for celltype in self.profile_df.columns:
            print("\tWorking on: {}".format(celltype))
            ct_genes = gene_celltypes[gene_celltypes == celltype].index
            ct_expr = shifted_ct_expr.loc[shifted_ct_expr.index.isin(ct_genes), :]
            print("\t  N = {}".format(len(ct_expr.index)))

            # perform NMF over the expression of these genes.
            print("\t  NMF")
            nmf_component = self.get_first_nmf_component(ct_expr)
            nmf_component_values = [x[0] for x in list(nmf_component)]
            nmf_data.append(nmf_component_values)

        # Create the data frame.
        celltype_cs = pd.DataFrame(nmf_data,
                                   index=["{}NMF_{}_C1".format(*x.split("_")) for x in self.profile_df.columns],
                                   columns=shifted_ct_expr.columns)

        return ct_expr_df, celltype_pcs, celltype_cs

    def save(self):
        save_dataframe(df=self.celltype_pcs, outpath=self.pca_outpath,
                       index=True, header=True)
        save_dataframe(df=self.celltype_cs, outpath=self.nmf_outpath,
                       index=True, header=True)

    @staticmethod
    def normalize(df):
        out_df = df.copy()
        out_df = out_df.loc[out_df.std(axis=1) > 0, :]
        out_df = out_df.subtract(out_df.mean(axis=1), axis=0).divide(out_df.std(axis=1), axis=0)
        return out_df

    @staticmethod
    def get_first_pca_component(X):
        corr_matrix = np.dot(X.T, X) / (X.shape[0] - 1)

        pca = PCA(n_components=1)
        pca.fit(corr_matrix)
        print("\t\tExplained variance ratio: {:.2f}".format(pca.explained_variance_ratio_[0]))
        print("\t\tSingular values: {:.2f}".format(pca.singular_values_[0]))
        return pca.transform(corr_matrix)

    @staticmethod
    def perform_shift(df):
        return df + abs(df.values.min())

    @staticmethod
    def get_first_nmf_component(X):
        corr_matrix = np.dot(X.T, X) / (X.shape[0] - 1)

        nmf = NMF(n_components=1)
        nmf.fit(corr_matrix)
        print("\t\tReconstruction error: {:.2f}".format(nmf.reconstruction_err_))
        print("\t\tNumber of iterations: {}".format(nmf.n_iter_))
        return nmf.transform(corr_matrix)

    def get_celltype_expression(self):
        return self.celltype_expression

    def get_celltype_pcs(self):
        return self.celltype_pcs

    def get_celltype_cs(self):
        return self.celltype_cs

    def clear_variables(self):
        self.profile_file = None
        self.profile = None
        self.ct_expr_file = None
        self.force = None

    def print_arguments(self):
        print("Arguments:")
        if self.profile_df is not None:
            print("  > Celltype profile: {}".format(self.profile_df.shape))
        else:
            print("  > Celltype profile input file: {}".format(self.profile_file))
        print("  > Celltype expression input path: {}".format(self.ct_expr_file))
        print("  > Celltype PCA output file: {}".format(self.pca_outpath))
        print("  > Celltype NMF output file: {}".format(self.nmf_outpath))
        print("  > Force: {}".format(self.force))
        print("")
