#!/usr/bin/env python3

"""
File:         gtex_nnls_deconvolution.py
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
import math
import ast
import re
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy.optimize import nnls
import scipy.cluster.hierarchy as sch
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "GTEx NNLS Deconvolution"
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
        self.expression_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/gtex_tissue/SamplesZTransformed.txt.gz"
        # self.expression_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_output/create_deconvolution_matrices/ct_profile_expr_table.txt.gz"
        self.profile_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt"
        self.palette = {
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Average": "#E8E8E8"
        }
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'gtex_tissue')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):

        print("Load the expression file.")
        expr_df = pd.read_csv(self.expression_path, sep="\t", header=[0, 1, 2], index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.expression_path),
                                      expr_df.shape))
        column_list = []
        for column in expr_df.columns:
            if isinstance(column, str):
                column_list.append(
                    ast.literal_eval(re.findall(r'(\(.*?\)).', column)[0]))
            else:
                column_list.append(column)
        expr_df.columns = pd.MultiIndex.from_tuples(column_list)

        print("Load the profile file.")
        profile_df = pd.read_csv(self.profile_path, sep="\t", header=0,
                                 index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.profile_path),
                                      profile_df.shape))
        profile_df.columns = [x.split("_")[1] for x in profile_df.columns]
        profile_df = self.remove_non_informative_genes(profile_df)
        profile_index = profile_df.index

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

        print("Profile shape: {}".format(profile_df.shape))
        print("Expression shape: {}".format(expr_df.shape))

        # Determine the order.
        print("Performing hiarachical clustering.")
        linkage = sch.linkage(self.zscore_transform(profile_df))
        deno = sch.dendrogram(linkage, orientation='right')
        order = deno['leaves']

        # Filter order.
        new_profile_index = []
        for x in profile_index:
            if x in overlap:
                new_profile_index.append(x)

        # Heatmap of the profile.
        # self.plot_profile(profile_df, "default", self.outdir)
        # self.plot_profile(self.normalize(profile_df), "normalized", self.outdir)
        # self.plot_profile(self.zscore_transform(profile_df), "zscore_transformed", self.outdir)
        # self.plot_profile(self.log2_transform(profile_df), "log2_transformed", self.outdir)
        # self.plot_profile(self.normalize(self.log2_transform(profile_df)), "log2_transformed_and_normalized", self.outdir)
        #exit()

        # # Normalize per cell type.
        # profile_df = self.normalize(profile_df)

        # # Normalize the cellmap profile.
        # profile_df = self.zscore_transform(profile_df)

        # Log2 transform the cellmap profile.
        profile_df = self.log2_transform(profile_df)

        print("Shifting data to be positive.")
        # Shift the data to be positive.
        if profile_df.values.min() < 0:
            print("\tshifting profile")
            profile_df = profile_df + abs(profile_df.values.min())

        if expr_df.values.min() < 0:
            print("\tshifting expression")
            expr_df = expr_df + abs(expr_df.values.min())

        # print("Rescale.")
        # # Rescale the cellmap profile.
        # scaled_expr_df = []
        # for ((_, a), (_, b)) in zip(expr_df.iterrows(), profile_df.iterrows()):
        #     scaled_expr_df.append(self.rescale(a, b))
        # expr_df = pd.DataFrame(scaled_expr_df, columns=expr_df.columns)

        # Plotting profile.
        # print("Plotting distributions.")
        # tmp = profile_df.copy()
        # tmp_melt = tmp.melt()
        # self.plot_distributions(tmp_melt, "profile", self.outdir)
        # self.plot_boxplot(tmp_melt, "profile", self.palette, self.outdir)
        # self.plot_profile_stripplot(profile_df, self.palette, self.outdir, order=new_profile_index)
        # self.plot_profile_boxplot(profile_df, self.palette, self.outdir)

        # Plot the profile.
        self.plot_profile(profile_df, "decon", self.outdir)

        # Perform deconvolution per sample.
        print("Performing deconvolution.")
        decon_data = []
        residuals_data = []
        count = 0
        for index, sample in expr_df.T.iterrows():
            proportions, rnorm = self.nnls(profile_df, sample)
            decon_data.append(proportions)
            residuals_data.append(rnorm)

            # if count < 5:
            #     self.plot_sample_bars(profile_df.copy(), sample.copy(),
            #                      proportions.copy(), order, self.outdir,
            #                      sample_id + "_" + tissue)
            #     count += 1

        decon_df = pd.DataFrame(decon_data,
                                index=pd.MultiIndex.from_tuples(expr_df.columns, names=['sample', 'broad_region', 'specific_region']),
                                columns=profile_df.columns)
        print(decon_df)
        print(decon_df.mean(axis=0))

        # # Make the weights sum up to 1.
        # decon_df = self.sum_to_one(decon_df)
        # min_fraction = 0
        # max_fraction = 1

        # Plot the fractions.
        count = 0
        for sample in decon_df.index.get_level_values('sample').unique():
            df = decon_df.xs(sample, level='sample', drop_level=True)
            if len(df.index) > 4 and count < 5:
                self.plot_sample_stripplot(df, self.palette, self.outdir, sample)
                count += 1

        # Melt the data.
        decon_df.reset_index(drop=False, inplace=True)
        decon_melt = decon_df.melt(id_vars=["broad_region"], value_vars=["Astrocyte", "EndothelialCell", "Macrophage", "Neuron", "Oligodendrocyte"])

        # Create boxplot per broad region.
        for broad_region in decon_melt["broad_region"].unique():
            df = decon_melt.loc[decon_melt["broad_region"] == broad_region, :].copy()
            print(broad_region)
            print(df.groupby(by="variable").mean())
            print("")
            if len(df.index) > 10:
                self.plot_boxplot(df,
                                  "{}_celltype_fraction".format(broad_region),
                                  self.palette,
                                  self.outdir)

        print("Plotting distributions.")
        self.plot_distributions(decon_melt, "celltype_fraction", self.outdir)
        self.plot_boxplot(decon_melt, "celltype_fraction", self.palette, self.outdir)

        # Construct a Series for the residuals.
        residuals_df = pd.DataFrame({"variable": np.nan,
                                     "value": residuals_data},
                                    index=expr_df.columns)
        print(residuals_df)
        print("Average residual: {:.2f}".format(residuals_df["value"].mean()))

    @staticmethod
    def remove_non_informative_genes(df, min_max_cpm=5):
        print("Filtering genes.")
        out_df = df.copy()
        pre_shape = out_df.shape
        out_df = out_df.loc[(out_df.std(axis=1) != 0) & (out_df.max(axis=1) > min_max_cpm), :]
        post_shape = out_df.shape
        print("\tdropped {} rows and {} columns".format(pre_shape[0] - post_shape[0], pre_shape[1] - post_shape[1]))
        return out_df

    @staticmethod
    def normalize(df):
        print("Normalize profile.")
        out_df = df.copy()
        out_df = out_df.divide(out_df.sum(axis=1), axis=0)
        out_df = out_df.divide(out_df.sum(axis=0), axis=1)
        return out_df

    @staticmethod
    def zscore_transform(df):
        print("Z-score transform profile.")
        out_df = df.copy()
        out_df = out_df.loc[out_df.std(axis=1) > 0, :]
        out_df = out_df.subtract(out_df.mean(axis=1), axis=0).divide(
            out_df.std(axis=1), axis=0)
        return out_df

    @staticmethod
    def log2_transform(df):
        print("Log2 transform profile.")
        out_df = df.copy()
        out_df = out_df + 1
        out_df = np.log2(out_df)
        return out_df

    @staticmethod
    def rescale(a, b):
        c = ((max(b) - min(a)) / (max(a) - min(a)) * (a - min(a))) + min(a)
        return c

    @staticmethod
    def nnls(A, b):
        return nnls(A, b)

    @staticmethod
    def sum_to_one(X):
        print("Sum-to-one weights.")
        return X.divide(X.sum(axis=1), axis=0)

    @staticmethod
    def plot_profile(df, title, outdir):
        sns.set(color_codes=True)
        g = sns.clustermap(df.T, center=0, cmap="RdBu_r",
                           yticklabels=True, xticklabels=False,
                           dendrogram_ratio=(.1, .1),
                           figsize=(12, 9))
        plt.setp(g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=10))
        g.fig.subplots_adjust(bottom=0.05, top=0.7)
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "{}_profile.png".format(title)))
        plt.close()

    @staticmethod
    def plot_profile_stripplot(data, palette, outdir, order=None):
        df = data.copy()
        if order is not None:
            df = df.loc[order, :]
        df = df.iloc[:25, :]
        value_vars = df.columns
        df.reset_index(inplace=True, drop=False)
        dfm = df.melt(id_vars=[data.index.name], value_vars=value_vars)

        sns.set(rc={'figure.figsize': (8, 12)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        g = sns.stripplot(x="value",
                          y=data.index.name,
                          hue="variable",
                          data=dfm,
                          size=15,
                          dodge=False,
                          orient="h",
                          palette=palette,
                          linewidth=1,
                          edgecolor="w",
                          jitter=0,
                          ax=ax)

        g.set_ylabel('',
                     fontsize=12,
                     fontweight='bold')
        g.set_xlabel('value',
                     fontsize=12,
                     fontweight='bold')
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=10)

        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "profile_stripplot.png"))
        plt.close()

    @staticmethod
    def plot_profile_boxplot(data, palette, outdir):
        df = data.copy()
        dfm = df.melt()

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.boxplot(x="variable", y="value", data=dfm, palette=palette)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "profile_boxplot.png".format()))
        plt.close()

    @staticmethod
    def plot_sample_bars(X, y, coefficients, order, outdir, name):
        X = X.iloc[order, :]
        y = y.iloc[order].to_frame()
        y.columns = ["Sample"]

        vmax = X.values.max()

        ncols = len(X.columns) * 2 + 6
        width_ratios = []
        for i in range(ncols):
            width = 1 / ncols
            if (i == 0) or (i % 2 == 0):
                width = 0.01
            width_ratios.append(width)

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=ncols,
                                 nrows=1,
                                 gridspec_kw={"width_ratios": width_ratios},
                                 figsize=(18, 5))

        i = 0
        combined = pd.Series(0, index=y.index)
        for i, (_, row) in enumerate(X.T.iterrows()):
            coefficient_prefix = '+ '
            coefficient_suffix = ' *'
            cbar = False
            if i == 0:
                coefficient_prefix = ''
            if i == len(X.columns) - 1:
                cbar = True

            weight_ax = axes[i*2]
            vector_ax = axes[(i*2)+1]

            weight_ax.text(0.5, 0.5, '{}{:.2f}{}'.format(coefficient_prefix, coefficients[i], coefficient_suffix),
                           {'size': 16, 'weight': 'bold'},
                           horizontalalignment='center',
                           verticalalignment='center',
                           transform=weight_ax.transAxes)
            weight_ax.set_axis_off()

            g = sns.heatmap(row.to_frame(), cmap="Blues",
                            vmin=0, vmax=vmax, cbar=cbar,
                            yticklabels=False, xticklabels=True, ax=vector_ax)
            vector_ax.set_ylabel('')

            combined = combined + (coefficients[i] * row)

        combined = combined.to_frame()
        combined.columns = ["Sum"]

        vmin = max(combined.values.min(), y.values.min()) - 2
        vmax = max(combined.values.max(), y.values.max())

        equals_ax = axes[(i+1)*2]
        equals_ax.text(0.5, 0.5, '=',
                       {'size': 16, 'weight': 'bold'},
                       horizontalalignment='center',
                       verticalalignment='center',
                       transform=equals_ax.transAxes)
        equals_ax.set_axis_off()

        combined_ax = axes[((i+1)*2)+1]
        g = sns.heatmap(combined, cmap="Greens",
                        vmin=vmin, vmax=vmax, cbar=False,
                        yticklabels=False, xticklabels=True, ax=combined_ax)
        combined_ax.set_ylabel('')

        mius_ax = axes[(i+2)*2]
        mius_ax.text(0.5, 0.5, '-',
                     {'size': 16, 'weight': 'bold'},
                     horizontalalignment='center',
                     verticalalignment='center',
                     transform=mius_ax.transAxes)
        mius_ax.set_axis_off()

        sample_ax = axes[((i+2)*2)+1]
        g = sns.heatmap(y, cmap="Greens",
                        vmin=vmin, vmax=vmax, cbar=True,
                        yticklabels=False, xticklabels=True, ax=sample_ax)
        sample_ax.set_ylabel('')

        diff = (combined.iloc[:, 0] - y.iloc[:, 0]).to_frame()
        diff.columns = ["\u0394"]

        equals_ax = axes[(i+3)*2]
        equals_ax.text(0.5, 0.5, '=',
                       {'size': 16, 'weight': 'bold'},
                       horizontalalignment='center',
                       verticalalignment='center',
                       transform=equals_ax.transAxes)
        equals_ax.set_axis_off()

        diff_ax = axes[((i+3)*2)+1]
        g = sns.heatmap(diff, center=0, cmap="RdBu_r",
                        yticklabels=False, xticklabels=True, ax=diff_ax)
        diff_ax.set_ylabel('')

        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "sample_{}_barplot.png".format(name.lower())))
        plt.close()

    @staticmethod
    def plot_sample_stripplot(df, palette, outdir, sample_id):
        print(df)

        average = df.mean(axis=0).to_frame()
        average.columns = ["value"]
        average["hue"] = "Average"
        average.reset_index(inplace=True, drop=False)

        xmin = df.values.min() - 0.1
        xmax = df.values.max() + 0.1

        ncols = min(len(df.index), 3)
        nrows = math.ceil(len(df.index) / 3)

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(ncols*8, nrows*5))
        plt.subplots_adjust(top=0.95, bottom=0.05, wspace=0.1, hspace=0.2)

        row_id = 0
        col_id = 0
        count = 0
        for i, ((broad_region, specific_region), row) in enumerate(df.iterrows()):
            count += 1
            weights = row.copy()
            weights = weights.to_frame()
            weights.columns = ["value"]
            weights["hue"] = weights.index
            weights.reset_index(inplace=True, drop=False)
            data = pd.concat([average, weights], ignore_index=True)

            ax = axes[row_id, col_id]
            sns.despine(fig=fig, ax=ax)

            g = sns.stripplot(x="value",
                              y="index",
                              hue="hue",
                              data=data,
                              size=25,
                              dodge=False,
                              orient="h",
                              palette=palette,
                              linewidth=1,
                              edgecolor="w",
                              jitter=0,
                              ax=ax)

            ax.legend_.remove()
            ax.set_xlim(xmin, xmax)

            ax.text(0.5, 1.05, broad_region.replace("_", " "),
                    fontsize=20, weight='bold', ha='center', va='bottom',
                    transform=ax.transAxes)
            ax.text(0.5, 0.98, specific_region.replace("_", " "),
                    fontsize=14, weight='light', style="italic", ha='center', va='bottom',
                    transform=ax.transAxes)

            g.set_ylabel('',
                         fontsize=12,
                         fontweight='bold')
            g.set_xlabel('',
                         fontsize=12,
                         fontweight='bold')
            ax.tick_params(axis='x', labelsize=10)
            ax.tick_params(axis='y', labelsize=14)

            ax.xaxis.grid(False)
            ax.yaxis.grid(True)

            if col_id != 0:
                ax.set_yticklabels(["" for _ in row.index])

            col_id += 1
            if col_id == 3:
                col_id = 0
                row_id += 1

        for i in range(count, (ncols * nrows)):
            axes[row_id, col_id].set_axis_off()
            col_id += 1
            if col_id == 3:
                col_id = 0
                row_id += 1

        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "sample_{}_stripplot.png".format(sample_id)))
        plt.close()

    def plot_distributions(self, df, name, outdir):
        sns.set(style="ticks", color_codes=True)
        g = sns.FacetGrid(df, row='broad_region', col='variable', sharex=True, sharey=True)
        g.map(sns.distplot, 'value')
        g.map(self.vertical_mean_line, 'value')
        for i, ax in enumerate(g.axes.flat):
            title = ax.get_title().split("|")
            broad_region = title[0].split(" = ")[-1]
            variable = title[1].split(" = ")[-1]
            ax.set_title(variable, fontweight='bold')
            if i % 5 == 0:
                ax.set_ylabel(broad_region, fontweight='bold')
            else:
                ax.set_ylabel("", fontweight='bold')
        #g.set_titles('{col_name}')
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "{}_distributions.png".format(name)))
        plt.close()

    @staticmethod
    def vertical_mean_line(x, **kwargs):
        plt.axvline(x.mean(), ls="--", c="black")

    @staticmethod
    def plot_boxplot(df, name, palette, outdir):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        for i in range(0, 200, 25):
            if (i / 100) > df['value'].max():
                break
            alpha = 0.025
            if i % 50 == 0:
                alpha = 0.15
            ax.axhline(i / 100, ls='-', color="#000000", alpha=alpha, zorder=-1)

        sns.boxplot(x="variable", y="value", data=df, palette=palette)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "{}_boxplot.png".format(name)))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
