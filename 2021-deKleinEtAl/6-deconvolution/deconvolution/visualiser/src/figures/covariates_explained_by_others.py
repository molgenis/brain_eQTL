"""
File:         covariates_explained_by_others.py
Created:      2020/04/15
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
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from colour import Color
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class CovariatesExplainedByOthers:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type format.
        """
        self.outdir = os.path.join(outdir, 'covariates_explained_by_others')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.groups = dataset.get_groups()
        self.cov_df = dataset.get_cov_df()
        self.colormap = self.create_color_map()
        self.tech_covs = ["PCT_CODING_BASES",
                          "PCT_MRNA_BASES",
                          "PCT_INTRONIC_BASES",
                          "MEDIAN_3PRIME_BIAS",
                          "PCT_USABLE_BASES",
                          "PCT_INTERGENIC_BASES",
                          "PCT_UTR_BASES",
                          #"PF_HQ_ALIGNED_READS",
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
                          "UCLA_ASD-EUR",
#                          "AMPAD-ROSMAP-V2-EUR"
                          ]

    @staticmethod
    def create_color_map():
        """
        """
        palette = list(Color("#8ABBDB").range_to(Color("#344A5A"), 101))
        colors = [str(x).upper() for x in palette]
        values = [x / 100 for x in list(range(101))]
        color_map = {}
        for val, col in zip(values, colors):
            color_map[val] = col
        return color_map

    def start(self):
        print("Plotting if covariates are explained by other covariates.")
        self.print_arguments()

        null_matrix = self.cov_df.loc[self.tech_covs, :].copy()

        r2_df, coef_df = self.model(null_matrix.T, self.cov_df)
        self.plot_bars(r2_df, self.groups, self.outdir, self.extension)

        subset = coef_df.loc[:, r2_df.loc[r2_df["value"] == 1.0, "index"].values].copy()
        self.plot_clustermap(subset, self.outdir, self.extension)

        for index, row in subset.T.iterrows():
            formula = ["{} = ".format(index)]
            row.dropna(inplace=True)
            for name, value in row.iteritems():
                if round(value) != 0:
                    formula.append("{} x {}".format(value, name))
            print(''.join(formula))

    def model(self, null_matrix, cov_df):
        print("Modelling each covariate with a linear model of the null matrix.")
        r2_data = []
        coef_df = None
        for index, y in cov_df.iterrows():
            X = null_matrix.copy()
            if index in X.columns:
                X.drop([index], axis=1, inplace=True)
            score, coefficients = self.create_model(X, y)
            color = self.colormap[round(score, 2)]
            r2_data.append([index, score, color])
            coefficients_df = pd.DataFrame([coefficients], index=[index], columns=X.columns).T
            for col in null_matrix.columns:
                if col not in coefficients_df.index:
                    coefficients_df.loc[col, index] = np.nan
            if coef_df is None:
                coef_df = coefficients_df
            else:
                coef_df = coef_df.merge(coefficients_df, left_index=True, right_index=True)

        r2_df = pd.DataFrame(r2_data, columns=["index", "value", "color"])

        return r2_df, coef_df

    @staticmethod
    def create_model(X, y):
        """
        Method for creating a multilinear model.

        :param X: DataFrame, the matrix with rows as samples and columns as
                             dimensions.
        :param y: Series, the outcome values.
        :return degrees_freedom: int, the degrees of freedom of this model.
        :return residual_squared_sum: float, the residual sum of squares of this
                                      fit.
        """
        # Create the model.
        regressor = LinearRegression()
        regressor.fit(X, y)
        y_hat = regressor.predict(X)

        # Calculate the statistics of the model.
        score = r2_score(y, y_hat)

        return score, regressor.coef_

    @staticmethod
    def plot_bars(df, groups, outdir, extension):
        print("Plotting")

        gridspec_kw = {"height_ratios": [x[1] - x[0] for x in groups],
                       "width_ratios": [0.3, 0.7]}

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=2, nrows=len(groups),
                                 figsize=(9, 28), gridspec_kw=gridspec_kw)
        plt.subplots_adjust(top=0.95, bottom=0.05, wspace=0.1, hspace=0.2)

        for i in range(len(groups)):
            print("\tPlotting axes[{}, 1]".format(i))
            axes[i, 0].set_axis_off()
            ax = axes[i, 1]
            sns.despine(fig=fig, ax=ax)
            (a, b, ylabel, remove) = groups[i]
            xlabel = ""
            if i == (len(groups) - 1):
                xlabel = "MLR R2"

            subset = df.iloc[a:b, :]

            for i in range(0, 105, 5):
                alpha = 0.025
                if i % 10 == 0:
                    alpha = 0.15
                ax.axvline(i / 100, ls='-', color="#000000",
                           alpha=alpha, zorder=-1)

            sns.barplot(x="value", y="index", data=df.iloc[a:b, :],
                        palette=subset["color"], orient="h", ax=ax)

            new_ylabels = [x.replace(remove, '').replace("_", " ") for x in
                           subset["index"]]
            ax.set_yticklabels(new_ylabels, fontsize=10)
            ax.set_ylabel(ylabel, fontsize=16, fontweight='bold')
            ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')
            ax.set(xlim=(0, 1))

        fig.align_ylabels(axes[:, 1])
        fig.suptitle('Variance Explained by other Covariates', fontsize=25,
                     fontweight='bold')
        fig.savefig(os.path.join(outdir, "covariates_explained_by_others.{}".format(extension)))
        plt.close()

    @staticmethod
    def plot_clustermap(df, outdir, extension):

        for col in df.columns:
            data = df.loc[:, col].copy()
            data = data.to_frame()
            data.dropna(inplace=True)

            sns.set(rc={'figure.figsize': (12, 9)})
            sns.set_style("ticks")
            fig, ax = plt.subplots()
            sns.despine(fig=fig, ax=ax)

            g = sns.barplot(x=col, y=data.index, data=data,
                            orient="h")
            g.set_title(col)
            g.set_ylabel('covariate',
                         fontsize=18,
                         fontweight='bold')
            g.set_xlabel('coefficient',
                         fontsize=18,
                         fontweight='bold')
            ax.tick_params(labelsize=10)
            ax.set_yticks(range(len(data.index)))
            ax.set_yticklabels(data.index)
            plt.tight_layout()
            fig.savefig(os.path.join(outdir,
                                     "{}_coef_bars.{}".format(col,
                                                              extension)))
            plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Groups: {}".format(self.groups))
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
