#!/usr/bin/env python3

"""
File:         direction_test.py
Created:      2020/05/22
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
import scipy.stats as stats
import statsmodels.api as sm
from functools import reduce
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Direction Test"
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
        self.geno_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_output/create_matrices/genotype_table.txt.gz"
        self.expr_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_nocovcorr/create_matrices/expression_table.txt.gz"
        self.cov_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_output/create_cov_matrix/covariates_table.txt.gz"
        self.index = 0
        self.drop_covs = ["AMPAD-ROSMAP-V2-EUR"]
        self.tech_covs = [
            "PCT_CODING_BASES",
            "PCT_MRNA_BASES",
            "PCT_INTRONIC_BASES",
            "MEDIAN_3PRIME_BIAS",
            "PCT_USABLE_BASES",
            "PCT_INTERGENIC_BASES",
            "PCT_UTR_BASES",
            "PCT_READS_ALIGNED_IN_PAIRS",
            "PCT_CHIMERAS",
            "PF_READS_IMPROPER_PAIRS",
            "PF_HQ_ALIGNED_Q20_BASES",
            "PF_HQ_ALIGNED_BASES",
            "PCT_PF_READS_IMPROPER_PAIRS",
            "PF_READS_ALIGNED",
            "avg_mapped_read_length",
            "avg_input_read_length",
            "uniquely_mapped",
            "total_reads",
            "Total.Sequences_R1",
            "MDS1",
            "MDS2",
            "MDS3",
            "MDS4",
            "AMPAD-MAYO-V2-EUR",
            "AMPAD-MSBB-V2-EUR",
            "BrainGVEX-V2-EUR",
            "CMC_HBCC_set2-EUR",
            "CMC_HBCC_set3-EUR",
            "CMC-EUR",
            "ENA-EUR",
            "GTEx-EUR",
            "GVEX-EUR",
            "LIBD_1M-EUR",
            "LIBD_h650-EUR",
            "NABEC-H550-EUR",
            "NABEC-H610-EUR",
            "TargetALS-EUR",
            "UCLA_ASD-EUR"
        ]
        self.covs = ["SEX", "CellMapNNLS_Astrocyte",
                     "CellMapNNLS_EndothelialCell",
                     "CellMapNNLS_Macrophage",
                     "CellMapNNLS_Neuron",
                     "CellMapNNLS_Oligodendrocyte"]
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Loading dataframes")
        geno_df = pd.read_csv(self.geno_path, sep="\t", header=0, index_col=0,
                              nrows=self.index + 1)
        print("\tGenotype matrix: {}".format(geno_df.shape))
        expr_df = pd.read_csv(self.expr_path, sep="\t", header=0, index_col=0,
                              nrows=self.index + 1)
        print("\tExpression matrix: {}".format(expr_df.shape))
        cov_df = pd.read_csv(self.cov_path, sep="\t", header=0, index_col=0)
        print("\tCovariate matrix: {}".format(cov_df.shape))

        # Drop the covariates we don't want.
        if len(self.drop_covs) > 0:
            cov_df.drop(self.drop_covs, axis=0, inplace=True)

        tech_cov_df = cov_df.loc[self.tech_covs, :].copy()
        cov_df = cov_df.loc[self.covs, :]
        geno_df.replace(-1, np.nan, inplace=True)

        # Get the missing genotype indices.
        indices = np.arange(geno_df.shape[1])
        eqtl_indices = indices[~geno_df.iloc[self.index, :].isnull().values]

        # Subset the row and present samples for this eQTL.
        genotype = geno_df.iloc[self.index, eqtl_indices].copy()
        expression = expr_df.iloc[self.index, eqtl_indices].copy()
        technical_covs = tech_cov_df.iloc[:, eqtl_indices].copy()
        covariates = cov_df.iloc[:, eqtl_indices].copy()

        self.plot_eqtl(genotype, expression, geno_df.index[self.index], self.outdir)

        # FLIP.
        # genotype = 2.0 - genotype

        # Create the null model. Null model are all the technical
        # covariates multiplied with the genotype + the SNP.
        tech_inter_matrix = technical_covs.mul(genotype, axis=1)
        tech_inter_matrix.index = ["{}_X_SNP".format(x) for x in
                                   technical_covs.index]
        intercept = pd.DataFrame(1, index=genotype.index,
                                 columns=["intercept"])
        base_matrix = reduce(lambda left, right: pd.merge(left,
                                                          right,
                                                          left_index=True,
                                                          right_index=True),
                             [intercept,
                              genotype.to_frame(),
                              technical_covs.T,
                              tech_inter_matrix.T])

        # Check the eQTL effect.
        print(base_matrix)
        print(genotype.name)
        print(base_matrix[["intercept", genotype.name]])
        df_eqtl, rss_eqtl, info = self.create_model(base_matrix[["intercept", genotype.name]],
                                                    expression)
        print("\ttdf_eqtl: {}\trss_eqtl: {}".format(df_eqtl, rss_eqtl))
        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None):
            print(info)

        # Loop over the covariates.
        for cov_index in range(len(cov_df.index)):
            # Get the covariate we are processing.
            covariate = covariates.iloc[cov_index, :]
            cov_name = covariate.name

            if cov_name != "CellMapNNLS_Oligodendrocyte":
                continue
            print("Cov: {}".format(cov_name))

            # Add the covariate to the null matrix if it isn't already.
            null_matrix = base_matrix.copy()
            if cov_name not in null_matrix.columns:
                covariate_df = covariate.copy()
                null_matrix = null_matrix.merge(covariate_df.to_frame(),
                                                left_index=True,
                                                right_index=True)

            print(null_matrix)

            # Create the null model.
            n_null = null_matrix.shape[0]
            df_null, rss_null, info = self.create_model(null_matrix,
                                                        expression)
            print("\tn_null: {}\tdf_null: {}\trss_null: {}".format(n_null,
                                                                   df_null,
                                                                   rss_null))
            with pd.option_context('display.max_rows', None,
                                   'display.max_columns', None):
                print(info)

            # Calculate the interaction effect of the covariate of
            # interest. Then drop the NA's from the interaction
            # term.
            inter_of_interest = covariate * genotype
            inter_name = "{}_X_SNP".format(cov_name)
            if inter_name in null_matrix.columns:
                inter_name = inter_name + "_2"
            inter_of_interest.name = inter_name

            # Create the alternative matrix and add the interaction
            # term.
            alt_matrix = null_matrix.copy()
            alt_matrix = alt_matrix.merge(inter_of_interest.to_frame(),
                                          left_index=True,
                                          right_index=True)

            # Create the alternative model.
            n_alt = alt_matrix.shape[0]
            df_alt, rss_alt, info = self.create_model(alt_matrix,
                                                      expression)
            print("\tn_alt: {}\tdf_alt: {}\trss_alt: {}".format(n_alt,
                                                                df_alt,
                                                                rss_alt))
            with pd.option_context('display.max_rows', None,
                                   'display.max_columns', None):
                print(info)

            # Make sure the n's are identical.
            if n_null != n_alt:
                print("\t\t\tError due to unequal n_null and n_alt",
                      flush=True)
                continue

            # Compare the null and alternative model.
            fvalue = self.calc_f_value(rss_null, rss_alt,
                                       df_null, df_alt, n_null)
            pvalue = self.get_p_value(fvalue, df_null, df_alt, n_null)
            zscore = self.get_z_score(pvalue)
            print("\tfvalue: {}\tpvalue: {}\tzscore: {}".format(fvalue, pvalue, zscore))

    @staticmethod
    def plot_eqtl(genotype, expression, snp_name, outdir):
        group = genotype.round(0)
        counts = group.value_counts()

        # Prepare the figure.
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        # Plot the scatter / box plot.
        sns.regplot(x=genotype, y=expression,
                    ax=ax
                    )
        sns.boxplot(x=group, y=expression,
                    zorder=-1,
                    ax=ax)
        plt.setp(ax.artists, edgecolor='k', facecolor='w')
        plt.setp(ax.lines, color='k')
        ax.set_xticks(range(3))
        ax.set_xticklabels(["{} [{}]".format(0.0, counts[0.0]),
                            "{} [{}]".format(1.0, counts[1.0]),
                            "{} [{}]".format(2.0, counts[2.0])])

        # Safe the plot.
        fig.savefig(os.path.join(outdir, "{}.png".format(snp_name.replace(":", "_"))))
        plt.close()

    @staticmethod
    def create_model(X, y):
        df = X.shape[1]

        # Perform the Ordinary least squares fit.
        ols = sm.OLS(y.values, X)
        try:
            ols_result = ols.fit()
        except np.linalg.LinAlgError as e:
            print("\t\tError: {}".format(e))
            return df, np.nan, None

        ssr = ols_result.ssr

        info = []
        for col in X.columns:
            coef = ols_result.params[col]
            std_err = ols_result.bse[col]
            if std_err > 0:
                tvalue = coef / std_err
            else:
                tvalue = np.nan
            info.append([col, coef, std_err, tvalue])
        info_df = pd.DataFrame(info)
        info_df.columns = ["term", "b", "se", "t-value"]

        return df, ssr, info_df

    @staticmethod
    def calc_f_value(rss1, rss2, df1, df2, n):
        if df1 >= df2:
            return np.nan
        if df2 >= n:
            return np.nan
        if rss2 >= rss1:
            return 0

        return ((rss1 - rss2) / (df2 - df1)) / (rss2 / (n - df2))

    @staticmethod
    def get_p_value(f_value, df1, df2, n):
        if f_value == np.nan:
            return np.nan
        if df1 >= df2:
            return np.nan
        if df2 >= n:
            return np.nan

        return stats.f.sf(f_value, dfn=(df2 - df1), dfd=(n - df2))

    @staticmethod
    def get_z_score(p_value):
        if p_value > (1.0 - 1e-16):
            p_value = (1.0 - 1e-16)
        if p_value < 1e-323:
            p_value = 1e-323
        return stats.norm.isf(p_value)


if __name__ == "__main__":
    m = main()
    m.start()
