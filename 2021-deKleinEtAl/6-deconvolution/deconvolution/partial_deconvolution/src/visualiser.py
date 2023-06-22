"""
File:         visualiser.py
Created:      2020/06/29
Last Changed: 2021/10/19
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
import os

# Third party imports.
import scipy.cluster.hierarchy as sch
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.


class Visualiser:
    def __init__(self, settings, signature, expression, deconvolution,
                 ground_truth, comparison):
        self.outdir = os.path.join(settings.get_outsubdir_path(), "plot")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.extension = settings.get_extension()
        self.title = settings.get_title()
        self.subtitle = settings.get_subtitle()
        self.ground_truth_type = settings.get_ground_truth_type()
        self.signature = signature
        self.expression = expression
        self.deconvolution = deconvolution
        self.ground_truth = ground_truth
        self.comparison = comparison
        order = self.set_order()
        if order is None:
            order = [x for x in range(self.signature.shape[0])]
        self.order = order
        self.palette = {
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Neuron": "#0072B2",
            "OtherNeuron": "#2690ce",
            "Oligodendrocyte": "#009E73",
            "OPC": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Endothelial": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Astrocytes": "#D55E00",
            "Pericytes": "#808080",
            "Microglia/Macrophage": "#E69F00",
            "Excitatory/Neuron": "#56B4E9",
            "Inhibitory/Neuron": "#0072B2",
            "Excitatory+Inhibitory/Neuron": "#BEBEBE",
            "-": "#000000",
            'Adult-Ex1': '#56B4E9',
            'Adult-Ex2': '#56B4E9',
            'Adult-Ex3': '#56B4E9',
            'Adult-Ex4': '#56B4E9',
            'Adult-Ex5': '#56B4E9',
            'Adult-Ex6': '#56B4E9',
            'Adult-Ex7': '#56B4E9',
            'Adult-Ex8': '#56B4E9',
            'Adult-In1': '#0072B2',
            'Adult-In2': '#0072B2',
            'Adult-In3': '#0072B2',
            'Adult-In4': '#0072B2',
            'Adult-In5': '#0072B2',
            'Adult-In6': '#0072B2',
            'Adult-In7': '#0072B2',
            'Adult-In8': '#0072B2',
            'Adult-Microglia': '#E69F00',
            'Adult-OPC': '#1b8569',
            'Adult-Endothelial': '#CC79A7',
            'Adult-Astrocytes': '#D55E00',
            'Adult-Oligo': '#009E73',
            'Adult-OtherNeuron': '#2690ce',
            'Dev-Replicating': '#000000',
            'Dev-Quiescent': '#808080'
        }

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

    def set_order(self):
        order = None
        try:
            linkage = sch.linkage(self.signature)
            deno = sch.dendrogram(linkage, orientation='right')
            order = deno['leaves']
        except RecursionError:
            print("RecursionError: maximum recursion depth exceeded while getting the str of an object")
            pass

        return order

    def plot_profile_clustermap(self):
        df = self.signature.copy()
        self.create_clustermap(df.T, name='signature')

    def plot_profile_correlations(self):
        df = self.signature.copy()
        correlation_m = np.corrcoef(df.T)
        correlation_m[np.triu_indices(correlation_m.shape[0])] = np.nan
        self.create_clustermap(pd.DataFrame(correlation_m, index=df.columns, columns=df.columns),
                               vmin=-1,
                               vmax=1,
                               xticklabels=True,
                               yticklabels=True,
                               row_cluster=False,
                               col_cluster=False,
                               add_annot=True,
                               name='signature_correlations')

    def plot_profile_boxplot(self):
        df = self.signature.copy()
        dfm = df.melt()
        self.create_boxplot(dfm, name='signature', order=df.columns.tolist())

    def plot_deconvolution_clustermap(self):
        df = self.deconvolution.copy()
        self.create_clustermap(df.T, name='deconvolution')

    def plot_deconvolution_correlations(self):
        df = self.deconvolution.copy()
        correlation_m = np.corrcoef(df.T)
        correlation_m[np.triu_indices(correlation_m.shape[0])] = np.nan
        self.create_clustermap(pd.DataFrame(correlation_m, index=df.columns, columns=df.columns),
                               vmin=-1,
                               vmax=1,
                               xticklabels=True,
                               yticklabels=True,
                               row_cluster=False,
                               col_cluster=False,
                               add_annot=True,
                               name='deconvolution_correlations')

    def plot_deconvolution_per_sample(self, n=1):
        sign_df = self.signature.copy()
        expr_df = self.expression.copy()
        decon_df = self.deconvolution.copy()

        for i in range(n):
            self.visualise_sample_deconvolution(X=sign_df,
                                                y=expr_df.iloc[:, i],
                                                coefficients=decon_df.iloc[i, :],
                                                name=expr_df.columns[i])

    def plot_deconvolution_distribution(self):
        df = self.deconvolution.copy()
        dfm = df.melt()
        self.create_distribution(dfm, name='deconvolution')

    def plot_deconvolution_boxplot(self):
        df = self.deconvolution.copy()
        dfm = df.melt()
        self.create_boxplot(dfm, name='deconvolution', order=df.columns.tolist())

    def plot_ground_truth_distribution(self):
        if self.ground_truth is None:
            return

        df = self.ground_truth.copy()
        dfm = df.melt()
        self.create_distribution(dfm, name='ground_truth')

    def plot_ground_truth_boxplot(self):
        if self.ground_truth is None:
            return

        df = self.ground_truth.copy()
        dfm = df.melt()
        self.create_boxplot(dfm, name='ground_truth', order=df.columns.tolist())

    def plot_prediction_comparison(self):
        if self.comparison is None:
            return

        df = self.comparison.copy()
        self.create_regression(df,
                               "NNLS predictions",
                               "{} counts".format(self.ground_truth_type),
                               '{}_comparison'.format(self.ground_truth_type))

    def plot_recon_accuracy_boxplot(self, s):
        df = s.to_frame()
        df.columns = ["-"]
        dfm = df.melt()
        self.create_boxplot(dfm, name='reconstruction_accuracy')

    def plot_violin_comparison(self, label, trans_dict):
        df = self.deconvolution.copy()
        df = df.reset_index().melt(id_vars=["index"])
        df[label] = df["index"].map(trans_dict).astype(str)
        self.create_catplot(df, label)

    def plot_profile_stripplot(self, n=25):
        df = self.signature.copy()
        df = df.iloc[self.order, :]
        id_vars = df.index.name
        value_vars = df.columns
        df = df.iloc[:n, :]
        df.reset_index(inplace=True, drop=False)
        dfm = df.melt(id_vars=[id_vars], value_vars=value_vars)

        sns.set(rc={'figure.figsize': (8, 12)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        g = sns.stripplot(x="value",
                          y=id_vars,
                          hue="variable",
                          data=dfm,
                          size=15,
                          dodge=False,
                          orient="h",
                          palette=self.palette,
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
        fig.savefig(os.path.join(self.outdir, "profile_stripplot.{}".format(self.extension)))
        plt.close()

    def visualise_sample_deconvolution(self, X, y, coefficients, name):
        X = X.iloc[self.order, :]
        y = y.iloc[self.order].to_frame()
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

            weight_ax = axes[i * 2]
            vector_ax = axes[(i * 2) + 1]

            weight_ax.text(0.5, 0.5, '{}{:.2f}{}'.format(coefficient_prefix,
                                                         coefficients[i],
                                                         coefficient_suffix),
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

        equals_ax = axes[(i + 1) * 2]
        equals_ax.text(0.5, 0.5, '=',
                       {'size': 16, 'weight': 'bold'},
                       horizontalalignment='center',
                       verticalalignment='center',
                       transform=equals_ax.transAxes)
        equals_ax.set_axis_off()

        combined_ax = axes[((i + 1) * 2) + 1]
        g = sns.heatmap(combined, cmap="Greens",
                        vmin=vmin, vmax=vmax, cbar=False,
                        yticklabels=False, xticklabels=True, ax=combined_ax)
        combined_ax.set_ylabel('')

        mius_ax = axes[(i + 2) * 2]
        mius_ax.text(0.5, 0.5, '-',
                     {'size': 16, 'weight': 'bold'},
                     horizontalalignment='center',
                     verticalalignment='center',
                     transform=mius_ax.transAxes)
        mius_ax.set_axis_off()

        sample_ax = axes[((i + 2) * 2) + 1]
        g = sns.heatmap(y, cmap="Greens",
                        vmin=vmin, vmax=vmax, cbar=True,
                        yticklabels=False, xticklabels=True, ax=sample_ax)
        sample_ax.set_ylabel('')

        diff = (combined.iloc[:, 0] - y.iloc[:, 0]).to_frame()
        diff.columns = ["\u0394"]

        equals_ax = axes[(i + 3) * 2]
        equals_ax.text(0.5, 0.5, '=',
                       {'size': 16, 'weight': 'bold'},
                       horizontalalignment='center',
                       verticalalignment='center',
                       transform=equals_ax.transAxes)
        equals_ax.set_axis_off()

        diff_ax = axes[((i + 3) * 2) + 1]
        g = sns.heatmap(diff, center=0, cmap="RdBu_r",
                        yticklabels=False, xticklabels=True, ax=diff_ax)
        diff_ax.set_ylabel('')

        plt.tight_layout()
        fig.savefig(
            os.path.join(self.outdir, "sample_{}_deconvolution.{}".format(name.lower(), self.extension)))
        plt.close()

    def create_clustermap(self, df, name, xticklabels=False, yticklabels=True,
                          row_cluster=True, col_cluster=True, vmin=None,
                          vmax=None, add_annot=False):
        annot = None
        if add_annot:
            annot = df.round(2)

        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           row_cluster=row_cluster, col_cluster=col_cluster,
                           yticklabels=yticklabels, xticklabels=xticklabels,
                           vmin=vmin, vmax=vmax, annot=annot,
                           annot_kws={"size": 10},
                           dendrogram_ratio=(.1, .1),
                           figsize=(12, 12))
        if yticklabels:
            plt.setp(g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=10))

        if xticklabels:
            plt.setp(g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=10))

        g.fig.subplots_adjust(bottom=0.05, top=0.7)
        plt.tight_layout()
        g.savefig(os.path.join(self.outdir, "{}_clustermap.{}".format(name, self.extension)))
        plt.close()

    def create_catplot(self, df, name):
        sns.set(style="ticks")
        order = list(df[name].unique())
        order.sort()
        g = sns.catplot(x=name, y="value", col="variable", col_wrap=3,
                        data=df, kind="violin", order=order)
        for axes in g.axes.flat:
            axes.set_xticklabels(axes.get_xticklabels(),
                                 rotation=65,
                                 horizontalalignment='right')
        g.savefig(os.path.join(self.outdir, "{}_catplot.{}".format(name, self.extension)))
        plt.close()

    def create_boxplot(self, df, xlabel="", ylabel="", name="", order=None):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, gridspec_kw={"height_ratios": [0.9, 0.1]})
        sns.despine(fig=fig, ax=ax1)
        ax2.set_axis_off()

        sns.boxplot(x="variable", y="value", data=df, order=order,
                    palette=self.palette, ax=ax1)
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
        ax1.set_xlabel(xlabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=14,
                       fontweight='bold')
        fig.savefig(os.path.join(self.outdir, "{}_boxplot.{}".format(name, self.extension)))
        plt.close()

    def create_distribution(self, df, name=""):
        sns.set(style="ticks", color_codes=True)
        g = sns.FacetGrid(df, col='variable', sharex=True, sharey=True)
        g.map(sns.distplot, 'value')
        g.map(self.vertical_mean_line, 'value')
        g.set_titles('{col_name}')
        plt.tight_layout()
        g.savefig(os.path.join(self.outdir, "{}_distributions.{}".format(name, self.extension)))
        plt.close()

    @staticmethod
    def vertical_mean_line(x, **kwargs):
        plt.axvline(x.mean(), ls="--", c="black")

    def create_regression(self, df, xlabel="", ylabel="", name=""):
        df['color'] = df['hue'].map(self.palette)

        total_coef, total_p = stats.pearsonr(df["x"], df["y"])

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        groups = df['hue'].unique()

        coefs = {}
        for group in groups:
            subset = df.loc[df['hue'] == group, :].copy()
            color = self.palette[group]

            #coef, p = stats.spearmanr(subset["x"], subset["y"])
            coef, p = stats.pearsonr(subset["x"], subset["y"])
            coefs[group] = "r = {:.2f}".format(coef)

            sns.regplot(x="x", y="y", data=subset,
                        scatter_kws={'facecolors': subset['color'],
                                     'edgecolors': "#808080"},
                        line_kws={"color": color},
                        ax=ax
                        )

        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
        ax.plot([-0.1, 1.1], [-0.1, 1.1], ls="--", c=".3")

        ax.text(0.5, 1.1, self.title,
                fontsize=18, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02, self.subtitle + ", r = {:.2f}".format(total_coef),
                fontsize=14, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        handles = []
        for key, color in self.palette.items():
            if key in groups and key in coefs.keys():
                handles.append(mpatches.Patch(color=color, label="{} [{}]".format(key, coefs[key])))
        ax.legend(handles=handles)

        fig.savefig(os.path.join(self.outdir, "{}_regression.{}".format(name, self.extension)))
        plt.close()