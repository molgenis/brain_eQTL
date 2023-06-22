#!/usr/bin/env python3

"""
File:         gene_expresion_cohort_correction.py
Created:      2020/08/31
Last Changed: 2020/09/02
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
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.api as sm

# Local application imports.

# Metadata
__program__ = "Gene Expression Cohort Correction"
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
        self.expression_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz"
        self.info_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-03-09.brain.phenotypes.txt"
        self.palette = {
            "AMP-AD": "#A9AAA9",
            "Braineac": "#E7A023",
            "Brainseq": "#CC79A9",
            "CMC": "#EEE643",
            "GTEx": "#1A9F74",
            "NABEC": "#767833",
            "TargetALS": "#DDCC78",
            "ENA": "#D56228",
            "PsychEncode": "#5AB4E5",
            # "BipSeq": "#050708",
            # "Braineac": "#E7A023",
            # "BrainGVEx": "##5AB4E5",
            # "Brainseq": "#CC79A9",
            # "CMC": "#EEE643",
            # "CMC HBCC": "#0873B4",
            # "MSBB": "#050708",
            # "NABEC": "#767833",
            # "ROSMAP": "#42AB9A",
            # "TargetALS": "#DDCC78",
            # "UCLA ASD": "#F26322",
            # "ENA": "#D56228",
            # "GTEx": "#1A9F74",
            # "Mayo": "#A9AAA9",
        }
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'gene_expression')
        self.nrows = 1

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        print("Load the expression")
        expr_df = pd.read_csv(self.expression_path, sep="\t", header=0,
                              index_col=0, nrows=self.nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.expression_path),
                                      expr_df.shape))

        print("Load the information")
        info_df = pd.read_csv(self.info_path, sep="\t", header=0,
                              index_col=4, low_memory=False)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.info_path),
                                      info_df.shape))
        sample_cohort_dict = dict(zip(info_df.index, info_df["cohort"]))
        sample_metacohort_dict = dict(zip(info_df.index, info_df["MetaCohort"]))
        cohort_metacohort_dict = dict(zip(info_df["cohort"], info_df["MetaCohort"]))
        # cohorts = {}
        # for key, value in metacohort_dict.items():
        #     if value in cohorts.keys():
        #         bla = cohorts[value]
        #         bla.append(key)
        #         cohorts[value] = bla
        #     else:
        #         cohorts[value] = [key]
        # for key, value in cohorts.items():
        #     print(key, value)
        # exit()

        print("Preprocessing the data.")

        # Get the overlap.
        overlap = []
        cohorts = {}
        for rnaseq_id, cohort in sample_cohort_dict.items():
            if cohort in cohorts.keys():
                samples = cohorts[cohort]
                samples.append(rnaseq_id)
                cohorts[cohort] = samples
            else:
                cohorts[cohort] = [rnaseq_id]

            if rnaseq_id in expr_df.columns:
                overlap.append(rnaseq_id)

        # Reorder.
        expr_df = expr_df.loc[:, overlap].T

        print("\tWorking with {} sample(s), {} gene(s), and {} cohort(s).".format(
            len(overlap), expr_df.shape[1], len(cohorts.keys())
        ))

        print("Creating cohort dataframe.")

        # Create a cohort dataframe.
        cohort_df = pd.DataFrame(0, index=overlap, columns=cohorts.keys())
        for cohort in cohort_df.columns:
            cohort_df.loc[cohort_df.index.isin(cohorts[cohort]), cohort] = 1

        print("Plotting cohort dataframe.")
        self.plot_cohorts(cohort_df.T, sample_metacohort_dict,
                          cohort_metacohort_dict, self.palette, self.outdir)

        print("Plotting genes.")
        for gene, expression in expr_df.T.iterrows():
            print("\tGene: {}".format(gene))
            df = expression.to_frame()
            df.reset_index(drop=False, inplace=True)
            df.columns = ["Sample", "Expression"]
            # df["Cohort"] = df["Sample"].map(sample_cohort_dict)
            df["MetaCohort"] = df["Sample"].map(sample_metacohort_dict)
            df["Color"] = df["MetaCohort"].map(self.palette)
            print(df)
            df["Residuals"] = self.remove_batch_effects(expression, cohort_df)
            df["CorrectedExpression"] = df["Expression"].mean() + df["Residuals"]

            self.plot_expression(gene, df, self.palette, self.outdir)

    @staticmethod
    def remove_batch_effects(expression, cohorts):
        ols = sm.OLS(expression, cohorts)
        try:
            ols_result = ols.fit()
            return ols_result.resid.values
        except np.linalg.LinAlgError as e:
            print("\t\tError: {}".format(e))
            return None

    @staticmethod
    def plot_cohorts(df, smc_dict, cmc_dict, pallete, outdir):
        row_colors = [pallete[cmc_dict[x]] for x in df.index]
        col_colors = [pallete[smc_dict[x]] for x in df.columns]

        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0.5, cmap="Greys",
                           yticklabels=False, xticklabels=False,
                           row_cluster=False, col_cluster=False,
                           row_colors=row_colors, col_colors=col_colors,
                           figsize=(12, 6))
        plt.setp(g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=10))
        g.savefig(os.path.join(outdir, "cohort_heatmap.png"))
        plt.close()

    @staticmethod
    def plot_expression(gene, df, palette, outdir):
        print(df)

        # Plot.
        fig, (ax1, ax2, ax3) = plt.subplots(ncols=1,
                                            nrows=3,
                                            gridspec_kw={"height_ratios": [0.48, 0.48, 0.04]},
                                            figsize=(15, 8))
        sns.despine(fig=fig, ax=ax1)
        sns.despine(fig=fig, ax=ax2)
        ax3.set_axis_off()

        sns.set()
        sns.set_style("ticks")

        ax1.set_title("Gene: {}".format(gene))

        for column, ax, legend, xlabel in [("Expression", ax1, False, ""),
                                           ("CorrectedExpression", ax2, "brief", "sample")]:
            sns.scatterplot(x=df.index,
                            y=df[column],
                            hue=df["MetaCohort"],
                            palette=palette,
                            linewidth=0,
                            legend=legend,
                            ax=ax)

            ax.axhline(df[column].mean(), ls='--', color="#b22222", zorder=-1)

            ax.set_ylabel("expression",
                          fontsize=10,
                          fontweight='bold')
            ax.set_xlabel(xlabel,
                          fontsize=10,
                          fontweight='bold')

            if legend == "brief":
                plt.setp(ax.get_legend().get_texts(), fontsize='2')
                plt.setp(ax.get_legend().get_title(), fontsize='4')
                ax.legend(bbox_to_anchor=(0.5, -0.45), loc='lower center',
                           ncol=5)

        fig.savefig(os.path.join(outdir, "{}_raw.png".format(gene)))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
