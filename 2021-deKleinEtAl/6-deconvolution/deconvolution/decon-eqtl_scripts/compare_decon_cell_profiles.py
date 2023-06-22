#!/usr/bin/env python3

"""
File:         compare_decon_cell_profiles.py
Created:      2021/09/17
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
from colour import Color
import argparse
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Compare Decon Cell Profiles"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
Syntax:
./compare_decon_cell_profiles.py -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-Normalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceNormalised.ExpAdded.txt -cc1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-NormalisedMAF5/deconvolutionResults.txt.gz -n1 original -cc2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table_CNS7.txt.gz -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-NormalisedMAF5-CNS7Profile/deconvolutionResults.txt.gz -n2 CNS7 -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_cohort_matrix/sample_to_dataset.txt.gz
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path1 = getattr(arguments, 'cell_counts1')
        self.decon_path1 = getattr(arguments, 'decon1')
        self.name1 = getattr(arguments, 'name1')
        self.cc_path2 = getattr(arguments, 'cell_counts2')
        self.decon_path2 = getattr(arguments, 'decon2')
        self.name2 = getattr(arguments, 'name2')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.n = getattr(arguments, 'number_of_plots')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'compare_decon_cell_profiles')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.ct_links = {
            "Astrocyte": "Astrocyte",
            "EndothelialCell": "EndothelialCell",
            "Macrophage": "Microglia",
            "Neuron": "Neuron",
            "Oligodendrocyte": "Oligodendrocyte",
            "Excitatory": "Neuron",
            "Inhibitory": "Neuron",
            "Pericytes": "Pericytes"
        }

        self.colormap = {
            "minor": "#E69F00",
            "center": "#0072B2",
            "major": "#D55E00"
        }
        self.palette = {0.0: "#D55E00",
                        1.0: "#0072B2",
                        2.0: "#E69F00"}
        self.ct_colormap = {
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Neuron": "#2690ce",
            "Oligodendrocyte": "#009E73",
            "OPC": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Pericytes": "#808080"
        }
        self.dataset_cohort_dict = {
            "AMPAD-MAYO-V2": "MAYO",
            "CMC_HBCC_set2": "CMC HBCC",
            "GTEx": "GTEx",
            "AMPAD-ROSMAP-V2": "ROSMAP",
            "BrainGVEX-V2": "Brain GVEx",
            "TargetALS": "Target ALS",
            "AMPAD-MSBB-V2": "MSBB",
            "NABEC-H610": "NABEC",
            "LIBD_1M": "LIBD",
            "ENA": "ENA",
            "LIBD_h650": "LIBD",
            "GVEX": "GVEX",
            "NABEC-H550": "NABEC",
            "CMC_HBCC_set3": "CMC HBCC",
            "UCLA_ASD": "UCLA ASD",
            "CMC": "CMC",
            "CMC_HBCC_set1": "CMC HBCC"
        }
        self.cohort_palette = {
            "MAYO": "#9c9fa0",
            "CMC HBCC": "#0877b4",
            "GTEx": "#0fa67d",
            "ROSMAP": "#6950a1",
            "Brain GVEx": "#48b2e5",
            "Target ALS": "#d5c77a",
            "MSBB": "#5cc5bf",
            "NABEC": "#6d743a",
            "LIBD": "#e49d26",
            "ENA": "#d46727",
            "GVEX": "#000000",
            "UCLA ASD": "#f36d2a",
            "CMC": "#eae453",
            "NA": "#808080"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-eq",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the eqtl matrix.")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
        parser.add_argument("-cc1",
                            "--cell_counts1",
                            type=str,
                            required=True,
                            help="The path to the first cell counts matrix.")
        parser.add_argument("-d1",
                            "--decon1",
                            type=str,
                            required=True,
                            help="The path to the first decon matrix.")
        parser.add_argument("-n1",
                            "--name1",
                            type=str,
                            required=True,
                            help="The name for the first analysis.")
        parser.add_argument("-cc2",
                            "--cell_counts2",
                            type=str,
                            required=True,
                            help="The path to the second cell counts matrix.")
        parser.add_argument("-d2",
                            "--decon2",
                            type=str,
                            required=True,
                            help="The path to the second decon matrix.")
        parser.add_argument("-n2",
                            "--name2",
                            type=str,
                            required=True,
                            help="The name for the second analysis.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=False,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-n",
                            "--number_of_plots",
                            type=int,
                            default=5,
                            help="The number of plots. Default: 5.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading color data")
        std_df = self.load_file(self.std_path, header=0, index_col=None)
        sample_to_dataset_dict = dict(zip(std_df.iloc[:, 0], std_df.iloc[:, 1]))
        sample_color_dict = {sample: self.cohort_palette[self.dataset_cohort_dict[sample_to_dataset_dict[sample]]] for sample in std_df.iloc[:, 0]}
        del std_df

        print("Loading decon data")
        decon_df1 = self.load_file(self.decon_path1, header=0, index_col=0)
        decon_df2 = self.load_file(self.decon_path2, header=0, index_col=0)

        # Find rows with top opposite effects.
        decon_df1["id"] = decon_df1.index
        decon_df1_m = decon_df1.melt(id_vars=["id"], value_vars=[x for x in decon_df1.columns if x.endswith("_pvalue")])
        decon_df1_m.columns = ["ID", "cell type", "value1"]
        decon_df1_m["cell type"] = [x.replace("_pvalue", "") for x in decon_df1_m["cell type"]]
        # decon_df1_m.columns = ["ID", "index", "cell type1", "value1"]
        # decon_df1_m["cell type1"] = [x.split("_")[1] for x in decon_df1_m["cell type1"]]
        # decon_df1_m["CT link"] = decon_df1_m["cell type1"].map(self.ct_links)
        # print(decon_df1_m)

        decon_df2["id"] = decon_df2.index
        decon_df2_m = decon_df2.melt(id_vars=["id"], value_vars=[x for x in decon_df2.columns if x.endswith("_pvalue")])
        decon_df2_m.columns = ["ID", "cell type", "value2"]
        decon_df2_m["cell type"] = [x.replace("_pvalue", "") for x in decon_df2_m["cell type"]]
        # decon_df2_m.columns = ["ID", "index", "cell type2", "value2"]
        # decon_df2_m["cell type2"] = [x.split("_")[1] for x in decon_df2_m["cell type2"]]
        # decon_df2_m["CT link"] = decon_df2_m["cell type2"].map(self.ct_links)
        # print(decon_df2_m)

        #decon_df_m = decon_df1_m.merge(decon_df2_m, on=["ID", "index", "CT link"])
        decon_df_m = decon_df1_m.merge(decon_df2_m, on=["ID", "cell type"])
        decon_df_m["abs diff"] = (decon_df_m["value1"] - decon_df_m["value2"]).abs()
        decon_df_m.sort_values(by="abs diff", ascending=False, inplace=True)
        interest_df = decon_df_m.head(self.n)
        print(interest_df)
        del decon_df1_m, decon_df2_m

        print("Loading eQTL data")
        eqtl_df = self.load_file(self.eqtl_path, header=0, index_col=None)
        eqtl_df.index = eqtl_df["ProbeName"] + "_" + eqtl_df["SNPName"]

        print("Loading other data")
        nrows = None
        geno_df = self.load_file(self.geno_path, header=0, index_col=0, nrows=nrows)
        geno_df = geno_df.groupby(geno_df.index).first()
        expr_df = self.load_file(self.expr_path, header=0, index_col=0, nrows=nrows)
        expr_df = expr_df.groupby(expr_df.index).first()
        cc_df1 = self.load_file(self.cc_path1, header=0, index_col=0)
        cc_df2 = self.load_file(self.cc_path2, header=0, index_col=0)

        print("Plotting")
        for _, (id, cell_type, value1, value2, _) in interest_df.iterrows():
            print("\tPlotting '{}'".format(id))

            # Split the ID.
            gene = id.split("_")[0]
            snp = id.replace(gene + "_", "")

            # Get the eqtl data.
            eqtl_data = eqtl_df.loc[id, :]
            genotype = geno_df.loc[snp, :]
            expression = expr_df.loc[gene, :]

            # Get the interaction data.
            decon1 = decon_df1.loc[id, :]
            decon2 = decon_df2.loc[id, :]
            cell_count1 = cc_df1.loc[:, cell_type]
            cell_count2 = cc_df2.loc[:, cell_type]

            # Construct the plot data frame.
            plot_df = pd.DataFrame({"genotype": genotype,
                                    "expression": expression,
                                    "cell count1": cell_count1,
                                    "cell count2": cell_count2})
            plot_df["round_geno"] = np.rint(genotype)
            plot_df["hue"] = plot_df.index.map(sample_color_dict)
            plot_df = plot_df.loc[plot_df["genotype"] != -1, :]

            # Prepare bargraph data.
            decon1_betas_df, decon1_pvalues_df = self.prepare_bargraph_data(decon1)
            decon2_betas_df, decon2_pvalues_df = self.prepare_bargraph_data(decon2)

            # Plot.
            sns.set(rc={'figure.figsize': (24, 27)})
            sns.set_style("ticks")
            fig, axes = plt.subplots(ncols=2, nrows=3, gridspec_kw={"height_ratios": [0.4, 0.2, 0.4]},)
            for ax in axes:
                sns.despine(fig=fig, ax=ax)

            # Plot eQTL.
            self.plot_eqtl(df=plot_df,
                           palette=self.palette,
                           ax=axes[0, 0],
                           title="Main eQTL effect",
                           xlabel=snp,
                           ylabel=gene,
                           annotate=[("eQTL p-value", eqtl_data["PValue"], ".2e"),
                                     ("Iteration", eqtl_data["Iteration"], ".0f")])

            # Plot cell fraction comparison.
            self.plot_regplot(df=plot_df,
                              x="cell count1",
                              y="cell count2",
                              color="hue",
                              ax=axes[0, 1],
                              xlabel=self.name1,
                              ylabel=self.name2,
                              title=cell_type + " cell fractions"
                              )

            # Plot Decon-QTL beta's 1.
            self.plot_barplot(df=decon1_betas_df,
                              x="abs beta",
                              y="index",
                              hue="cell type",
                              yticklabels="label",
                              palette=self.ct_colormap,
                              xlabel="|beta|",
                              title="{} Decon-QTL output".format(self.name1),
                              annotate=[(index, value, ".2e") for _, (_, value, index) in decon1_pvalues_df.iterrows()],
                              ax=axes[1, 0])

            # Plot Decon-QTL beta's 2.
            self.plot_barplot(df=decon2_betas_df,
                              x="abs beta",
                              y="index",
                              hue="cell type",
                              yticklabels="label",
                              palette=self.ct_colormap,
                              xlabel="|beta|",
                              title="{} Decon-QTL output".format(self.name2),
                              annotate=[(index, value, ".2e") for _, (_, value, index) in decon2_pvalues_df.iterrows()],
                              ax=axes[1, 1])

            # Plot interaction 1.
            self.plot_inter_eqtl(df=plot_df,
                                 x="cell count1",
                                 y="expression",
                                 palette=self.palette,
                                 ax=axes[2, 0],
                                 title="{} interaction".format(self.name1),
                                 xlabel=cell_type,
                                 ylabel=gene,
                                 annotate=[
                                     ("Decon-eQTL p-value", decon1[cell_type + "_pvalue"], ".2e"),
                                     ("Decon-eQTL beta", decon1[[col for col in decon_df1.columns if col.endswith(cell_type + ":GT")][0]], ".2f")])

            # Plot interaction 2.
            self.plot_inter_eqtl(df=plot_df,
                                 x="cell count2",
                                 y="expression",
                                 palette=self.palette,
                                 ax=axes[2, 1],
                                 title="{} interaction".format(self.name2),
                                 xlabel=cell_type,
                                 ylabel=gene,
                                 annotate=[
                                     ("Decon-eQTL p-value", decon2[cell_type + "_pvalue"], ".2e"),
                                     ("Decon-eQTL beta", decon2[[col for col in decon_df2.columns if col.endswith(cell_type + ":GT")][0]], ".2f")])

            outpath = os.path.join(self.outdir, "compare_decon_cell_profile_{}_{}_{}.png".format(gene, snp, cell_type))
            fig.savefig(outpath)
            plt.close()
            print("\tSaved: {}".format(outpath))

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=None, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def prepare_bargraph_data(decon_s):
        decon_betas_df = decon_s.loc[[x for x in decon_s.index if "Beta" in x]].to_frame()
        decon_betas_df.reset_index(drop=False, inplace=True)
        decon_betas_df.columns = ["variable", "beta"]
        decon_betas_df["index"] = np.arange(1, decon_betas_df.shape[0] + 1)
        decon_betas_df["abs beta"] = decon_betas_df["beta"].abs()
        decon_betas_df["label"] = [x.split("_")[2] for x in decon_betas_df["variable"]]
        decon_betas_df["cell type"] = [x.replace(":GT", "") for x in decon_betas_df["label"]]

        pvalue_cols = [x for x in decon_s.index if x.endswith("_pvalue")]
        decon_pvalues_df = decon_s.loc[pvalue_cols].to_frame()
        decon_pvalues_df.reset_index(drop=False, inplace=True)
        decon_pvalues_df.columns = ["variable", "p-value"]
        decon_pvalues_df["index"] = np.arange(len(pvalue_cols), len(pvalue_cols) * 2)

        return decon_betas_df, decon_pvalues_df

    @staticmethod
    def plot_eqtl(df, palette, ax, title="", xlabel="", ylabel="", annotate=None):
        # Calculate the correlation.
        coef, _ = stats.pearsonr(df["genotype"], df["expression"])

        # Plot the scatter / box plot.
        sns.regplot(x="genotype", y="expression", data=df,
                    scatter=False,
                    line_kws={"color": "#000000"},
                    ax=ax
                    )
        sns.boxplot(x="round_geno", y="expression", data=df,
                    palette=palette,
                    showfliers=False,
                    zorder=1,
                    ax=ax)

        ax.annotate(
            'N = {:,}'.format(df.shape[0]),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=12,
            fontweight='bold')
        ax.annotate(
            'r = {:.2f}'.format(coef),
            xy=(0.03, 0.90),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=14,
            fontweight='bold')

        if annotate is not None:
            for i, (label, value, rounding) in enumerate(annotate):
                ax.annotate(
                    '{} = {:{}}'.format(label, value, rounding),
                    xy=(0.03, 0.86 - (0.04 * i)),
                    xycoords=ax.transAxes,
                    color="#000000",
                    alpha=0.75,
                    fontsize=14,
                    fontweight='bold')

        ax.set_title(title,
                     fontsize=22,
                     weight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

    @staticmethod
    def plot_barplot(ax, df, x="x", y="y", hue="hue", yticklabels=None, palette=None,
                     xlabel="", ylabel="", title="", annotate=None):

        # Plot.
        g = sns.barplot(x=x, y=y, hue=hue, palette=palette, dodge=False,
                        data=df, orient="h", ax=ax)
        g.legend_.remove()

        if annotate is not None:
            xlim = ax.get_xlim()
            for (ypos, value, rounding) in annotate:
                g.text(xlim[1] * 0.9,
                       ypos,
                       '{:{}}'.format(value, rounding),
                       color="#000000",
                       alpha=1,
                       fontsize=12,
                       ha="center",
                       va="center"
                       )

        if yticklabels is not None:
            ax.set_yticklabels(df[yticklabels])

        ax.set_title(title,
                     fontsize=22,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

    @staticmethod
    def plot_regplot(ax, df, x="x", y="y", color=None, xlabel="", ylabel="",
                     title="", annotate=None):
        # Set color
        point_color = "#000000"
        reg_color = "#b22222"
        if color is not None:
            point_color = df[color]
            reg_color = "#000000"

        # Plot.
        sns.regplot(x=x, y=y, data=df, ci=95,
                    scatter_kws={'facecolors': point_color,
                                 'linewidth': 0,
                                 'alpha': 0.75},
                    line_kws={"color": reg_color},
                    ax=ax
                    )

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.set_xlim(-0.01, xlim[1])
        ax.set_ylim(-0.01, ylim[1])
        max_pos = max(xlim[1], ylim[1]) + 0.05
        ax.plot([-0.01, max_pos], [-0.01, max_pos], ls='--', color="#000000", zorder=-1)

        # Regression.
        coef, _ = stats.spearmanr(df[y], df[x])

        # Add the text.
        ax.annotate(
            'r = {:.2f}'.format(coef),
            xy=(0.03, 0.9),
            xycoords=ax.transAxes,
            color=reg_color,
            alpha=0.75,
            fontsize=12,
            fontweight='bold')

        if annotate is not None:
            for i, (label, value, rounding) in enumerate(annotate):
                ax.annotate(
                    '{} = {:{}}'.format(label, value, rounding),
                    xy=(0.03, 0.86 - (0.04 * i)),
                    xycoords=ax.transAxes,
                    color="#000000",
                    alpha=0.75,
                    fontsize=14,
                    fontweight='bold')

        ax.set_title(title,
                     fontsize=22,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

    @staticmethod
    def plot_inter_eqtl(df, palette, ax, x="x", y="y",
                        title="", xlabel="", ylabel="",
                        annotate=None):

        for i, genotype in enumerate([0.0, 1.0, 2.0]):
            subset = df.loc[df["round_geno"] == genotype, :].copy()
            color = palette[genotype]
            coef = np.nan
            if len(subset.index) > 1:
                # Calculate the correlation.
                coef, _ = stats.pearsonr(subset[x], subset[y])

                # Plot the scatter / box plot.
                sns.regplot(x=x, y=y, data=subset,
                            scatter_kws={'facecolors': color,
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": color, "alpha": 0.75},
                            ax=ax
                            )

            ax.annotate(
                '{}: r = {:.2f} [N = {:,}]'.format(genotype, coef, subset.shape[0]),
                xy=(0.03, 0.94 - (0.04 * i)),
                xycoords=ax.transAxes,
                color=color,
                alpha=0.75,
                fontsize=14,
                fontweight='bold')

        if annotate is not None:
            for i, (label, value, rounding) in enumerate(annotate):
                ax.annotate(
                    '{} = {:{}}'.format(label, value, rounding),
                    xy=(0.03, 0.82 - (0.04 * i)),
                    xycoords=ax.transAxes,
                    color="#000000",
                    alpha=0.75,
                    fontsize=14,
                    fontweight='bold')

        ax.set_title(title,
                     fontsize=22,
                     weight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL: {}".format(self.eqtl_path))
        print("  > Genotype: {}".format(self.geno_path))
        print("  > Expression: {}".format(self.expr_path))
        print("  > {} [cell profile 1]:".format(self.name1))
        print("    > Cell counts: {}".format(self.cc_path1))
        print("    > Decon: {}".format(self.decon_path1))
        print("  > {} [cell profile 2]:".format(self.name2))
        print("    > Decon: {}".format(self.decon_path2))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > N: {}".format(self.n))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
