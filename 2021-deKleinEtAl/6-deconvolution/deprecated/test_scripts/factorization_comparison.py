#!/usr/bin/env python3

"""
File:         factorization_comparison.py
Created:      2020/06/16
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
import os

# Third party imports.
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, NMF
from sklearn.metrics import explained_variance_score
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Factorization Comparison"
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
        self.expression_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_output/create_deconvolution_matrices/ct_profile_expr_table.txt.gz"
        self.profile_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt"
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'factorization_comparison')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        #Load the expression file.
        print("Load the expression file")
        expr_df = pd.read_csv(self.expression_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.expression_path),
                                      expr_df.shape))

        print("Load the profile")
        profile_df = pd.read_csv(self.profile_path, sep="\t", header=0,
                                 index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.profile_path),
                                      profile_df.shape))

        # Filter on overlap.
        overlap = np.intersect1d(profile_df.index, expr_df.index)
        profile_df = profile_df.loc[overlap, :]
        expr_df = expr_df.loc[overlap, :]

        # Set index name identical.
        profile_df.index.name = "-"
        expr_df.index.name = "-"

        # Check if identical.
        if not profile_df.index.equals(expr_df.index):
            print("Invalid order.")
            exit()

        # Find the genes specific to each celltype.
        gene_celltypes = self.normalize(profile_df).idxmax(axis=1)

        # Construct a dataframe of the first component of each celltype
        # subset expression profile.
        pca_data = []
        pca_var_explained = {}
        print("Performing PCA")
        for celltype in profile_df.columns:
            print("\tWorking on: {}".format(celltype))
            ct_genes = gene_celltypes[gene_celltypes == celltype].index
            ct_expr = expr_df.loc[expr_df.index.isin(ct_genes), :]
            print("\t  N = {}".format(len(ct_expr.index)))

            # # Create PCA plot.
            # self.plot_pca(ct_expr, celltype, self.outdir)

            # perform PCA over the expression of these genes.
            print("\t  PCA")
            pca_component, explained_var = self.get_first_pca_component(ct_expr)
            pca_component_values = [x[0] for x in list(pca_component)]
            pca_data.append(pca_component_values)
            pca_var_explained[celltype] = explained_var

        # Create the data frame.
        celltype_pcs = pd.DataFrame(pca_data,
                                    index=profile_df.columns,
                                    columns=expr_df.columns)
        self.plot_distributions(celltype_pcs, "PCA", self.outdir)

        # Shift the expression to be all positive.
        shifted_ct_expr = expr_df.copy()
        shifted_ct_expr = shifted_ct_expr + abs(shifted_ct_expr.values.min())

        # Construct a dataframe of the first component of each celltype
        # subset expression profile.
        nmf_data = []
        nmf_recon_error = {}
        print("Performing NMF")
        for celltype in profile_df.columns:
            print("\tWorking on: {}".format(celltype))
            ct_genes = gene_celltypes[gene_celltypes == celltype].index
            ct_expr = shifted_ct_expr.loc[shifted_ct_expr.index.isin(ct_genes), :]
            print("\t  N = {}".format(len(ct_expr.index)))

            # # Create NMF plot.
            # self.plot_nmf(ct_expr, celltype, self.outdir)

            # perform NMF over the expression of these genes.
            print("\t  NMF")
            nmf_component, recon_error = self.get_first_nmf_component(ct_expr)
            nmf_component_values = [x[0] for x in list(nmf_component)]
            nmf_data.append(nmf_component_values)
            nmf_recon_error[celltype] = recon_error

        # Create the data frame.
        celltype_cs = pd.DataFrame(nmf_data,
                                   index=profile_df.columns,
                                   columns=shifted_ct_expr.columns)
        self.plot_distributions(celltype_cs, "NMF", self.outdir)

        # # plot.
        # for celltype in profile_df.columns:
        #     df = pd.DataFrame({"PCA": celltype_pcs.loc[celltype, :], "NMF": celltype_cs.loc[celltype, :]})
        #     print(df)
        #     info = {"X": ("PCA", "variance explained: {:.2f}".format(pca_var_explained[celltype])),
        #             "Y": ("NMF", "reconstruction error: {:.2f}".format(nmf_recon_error[celltype]))}
        #     self.plot(df, info, celltype, self.outdir)

    @staticmethod
    def normalize(df):
        out_df = df.copy()
        out_df = out_df.loc[out_df.std(axis=1) > 0, :]
        out_df = out_df.subtract(out_df.mean(axis=1), axis=0).divide(
            out_df.std(axis=1), axis=0)
        return out_df

    def plot_pca(self, X, cell, outdir):
        corr_matrix = np.dot(X.T, X) / (X.shape[0] - 1)

        pca = PCA(n_components=2)
        pca.fit(corr_matrix)
        projections = pd.DataFrame(pca.transform(corr_matrix))
        projections.columns = ["PC1", "PC2"]
        info = {"X": ("PC1", " [{:.0f}%]".format(pca.explained_variance_ratio_[0]*100)),
                "Y": ("PC2", " [{:.0f}%]".format(pca.explained_variance_ratio_[1]*100))}
        self.plot(projections, info, cell, outdir)

    @staticmethod
    def get_first_pca_component(X):
        corr_matrix = np.dot(X.T, X) / (X.shape[0] - 1)

        pca = PCA(n_components=1)
        pca.fit(corr_matrix)
        print("\t\tExplained variance ratio: {:.2f}".format(pca.explained_variance_ratio_[0]))
        print("\t\tSingular values: {:.2f}".format(pca.singular_values_[0]))
        return pca.transform(corr_matrix), pca.explained_variance_ratio_[0]

    @staticmethod
    def get_first_nmf_component(X):
        corr_matrix = np.dot(X.T, X) / (X.shape[0] - 1)

        nmf = NMF(n_components=1)
        nmf.fit(corr_matrix)
        print("\t\tReconstruction error: {:.2f}".format(nmf.reconstruction_err_))
        print("\t\tNumber of iterations: {}".format(nmf.n_iter_))
        return nmf.transform(corr_matrix), nmf.reconstruction_err_

    def plot_nmf(self, X, cell, outdir):
        corr_matrix = np.dot(X.T, X) / (X.shape[0] - 1)

        nmf = NMF(n_components=2)
        nmf.fit(corr_matrix)
        projections = pd.DataFrame(nmf.transform(corr_matrix))
        projections.columns = ["LF1", "LF2"]
        info = {"X": ("LF1", ""),
                "Y": ("LF2", "")}
        self.plot(projections, info, cell, outdir)

    def plot_distributions(self, data, method, outdir):
        df = data.copy()
        df.index = [x.split("_")[-1] for x in df.index]
        dfm = df.T.melt()
        print(dfm)

        sns.set(style="ticks", color_codes=True)
        g = sns.FacetGrid(dfm, col='variable', sharex=True, sharey=True)
        g.map(sns.distplot, 'value')
        g.map(self.vertical_mean_line, 'value')
        g.set_titles('{col_name}')
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "{}_distributions.png".format(method)))
        plt.close()

    @staticmethod
    def vertical_mean_line(x, **kwargs):
        plt.axvline(x.mean(), ls="--", c="black")

    @staticmethod
    def plot(df, info, cell, outdir):
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")

        g = sns.scatterplot(x=info["X"][0],
                            y=info["Y"][0],
                            data=df,
                            legend=False,
                            alpha=0.5)

        g.set_title("")
        g.set_ylabel(' '.join(info["Y"]),
                     fontsize=10,
                     fontweight='bold')
        g.set_xlabel(' '.join(info["X"]),
                     fontsize=10,
                     fontweight='bold')

        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "{}_{}_vs_{}.png".format(cell, info["X"][0], info["Y"][0])))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
