"""
File:         inter_eqtl_effect.py
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
from colour import Color
import os

# Third party imports.
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.
from general.utilities import prepare_output_dir, p_value_to_symbol


class IntereQTLEffect:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_effect')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.eqtl_df = dataset.get_eqtl_df()
        self.geno_df = dataset.get_geno_df()
        self.expr_df = dataset.get_expr_df()
        self.alleles_df = dataset.get_alleles_df()
        self.cov_df = dataset.get_cov_df()
        self.inter_df = dataset.get_inter_cov_zscore_df()
        self.z_score_cutoff = dataset.get_significance_cutoff()
        colormap = dataset.get_colormap()

        # Create color map.
        self.group_color_map, self.value_color_map = self.create_color_map(colormap)
        self.sex_color_map = {"Male": colormap["male"], "Female": colormap["female"]}

    def start(self):
        print("Plotting interaction eQTL plots.")
        self.print_arguments()

        print("Iterating over eQTLs.")
        for i, (index, row) in enumerate(self.eqtl_df.iterrows()):
            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            eqtl_type = row["CisTrans"]

            # if hgnc_name != "CLECL1":
            #     continue

            print("\tWorking on: {}\t{}\t{} [{}/{} "
                  "{:.2f}%]".format(snp_name, probe_name, hgnc_name,
                                    i + 1,
                                    self.eqtl_df.shape[0],
                                    (100 / self.eqtl_df.shape[0]) * (i + 1)))

            # Get the genotype / expression data.
            genotype = self.geno_df.iloc[i, :].T.to_frame()
            expression = self.expr_df.iloc[i, :].T.to_frame()
            data = genotype.merge(expression, left_index=True, right_index=True)
            data.columns = ["genotype", "expression"]
            data["group"] = data["genotype"].round(0)

            # Remove missing values.
            data = data.loc[(data['genotype'] >= 0.0) &
                            (data['genotype'] <= 2.0), :]

            # Get the allele data.
            (alleles, _) = self.alleles_df.iloc[i, :]
            # A/T = 0.0/2.0
            # by default we assume T = 2.0 to be minor
            minor_allele = alleles[-1]
            major_allele = alleles[0]

            # Check if we need to flip the genotypes.
            counts = data["group"].value_counts()
            for x in [0.0, 1.0, 2.0]:
                if x not in counts:
                    counts.loc[x] = 0
            zero_geno_count = (counts[0.0] * 2) + counts[1.0]
            two_geno_count = (counts[2.0] * 2) + counts[1.0]
            if two_geno_count > zero_geno_count:
                # Turns out that 0.0 was the minor.
                minor_allele = alleles[0]
                major_allele = alleles[-1]
                data["genotype"] = 2.0 - data["genotype"]
                data["group"] = 2.0 - data["group"]

            allele_map = {0.0: "{}/{}".format(major_allele, major_allele),
                          1.0: "{}/{}".format(major_allele, minor_allele),
                          2.0: "{}/{}".format(minor_allele, minor_allele)}
            data["alleles"] = data["group"].map(allele_map)

            # Add the color.
            data["round_geno"] = data["genotype"].round(2)
            data["value_hue"] = data["round_geno"].map(self.value_color_map)
            data["group_hue"] = data["group"].map(self.group_color_map)

            # Check if the SNP has an interaction effect.
            interaction_effect = self.inter_df.iloc[:, i].to_frame()
            interaction_effect.columns = ["zscore"]
            interaction_effect = interaction_effect.loc[
                                 interaction_effect["zscore"] > abs(
                                     self.z_score_cutoff), :]
            interaction_effect = interaction_effect.reindex(
                interaction_effect["zscore"].abs().sort_values(
                    ascending=False).index)

            # Prepare output directory.
            if len(interaction_effect.index) > 0:
                eqtl_interaction_outdir = os.path.join(self.outdir,
                                                       "{}_{}_{}_{}".format(
                                                           index, snp_name,
                                                           probe_name,
                                                           hgnc_name))
                if not os.path.exists(eqtl_interaction_outdir):
                    os.makedirs(eqtl_interaction_outdir)

                count = 0
                for index2, (row,) in interaction_effect.iterrows():
                    if index2 not in ["SEX", "CellMapNNLS_Astrocyte",
                                      "CellMapNNLS_EndothelialCell",
                                      "CellMapNNLS_Macrophage",
                                      "CellMapNNLS_Neuron",
                                      "CellMapNNLS_Oligodendrocyte"]:
                        continue

                    eqtl_data = data.copy()
                    cov_data = self.cov_df.loc[index2].to_frame()
                    eqtl_data = eqtl_data.merge(cov_data, left_index=True,
                                                right_index=True)

                    if len(eqtl_data[index2].value_counts().index) == 2:
                        self.plot_box(snp_name, probe_name, hgnc_name,
                                      eqtl_type, eqtl_data, index2, row, count,
                                      allele_map, eqtl_interaction_outdir,
                                      self.sex_color_map, self.extension)
                    else:
                        self.plot_inter(snp_name, probe_name, hgnc_name,
                                        eqtl_type, eqtl_data, index2, row,
                                        count, allele_map, self.group_color_map,
                                        eqtl_interaction_outdir, self.extension)
                    count += 1

    @staticmethod
    def create_color_map(colormap):
        major_small = list(Color(colormap["major"]).range_to(Color("#FFFFFF"), 12))[2]
        center_small = list(Color(colormap["center"]).range_to(Color("#FFFFFF"), 12))[2]

        palette = list(Color(colormap["major"]).range_to(Color(major_small), 50)) + \
                  list(Color(major_small).range_to(Color(colormap["center"]), 50)) + \
                  list(Color(colormap["center"]).range_to(Color(center_small), 50)) + \
                  list(Color(center_small).range_to(Color(colormap["minor"]), 51))
        colors = [str(x).upper() for x in palette]
        values = [x / 100 for x in list(range(201))]
        group_color_map = {0.0: colormap["major"], 1.0: colormap["center"], 2.0: colormap["minor"]}
        value_color_map = {}
        for val, col in zip(values, colors):
            value_color_map[val] = col
        return group_color_map, value_color_map

    @staticmethod
    def plot_inter(snp_name, probe_name, hgnc_name, eqtl_type, df, cov_name,
                   zscore, count, allele_map, group_color_map, outdir, extension):
        # calculate axis limits.
        ymin_value = df["expression"].min()
        ymin = ymin_value - abs(ymin_value * 0.2)
        ymax_value = df["expression"].max()
        ymax = ymax_value + abs(ymax_value * 0.6)

        xmin_value = df[cov_name].min()
        xmin = xmin_value - max(abs(xmin_value * 0.1), 0.05)
        xmax_value = df[cov_name].max()
        xmax = xmax_value + max(abs(xmax_value * 0.1), 0.05)

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        label_pos = {0.0: 0.94, 1.0: 0.90, 2.0: 0.86}
        for i, genotype in enumerate([1.0, 0.0, 2.0]):
            subset = df.loc[df["round_geno"] == genotype, :].copy()
            color = group_color_map[genotype]
            allele = allele_map[genotype]

            coef_str = "NA"
            p_str = "NA"
            if len(subset.index) > 1:
                # Regression.
                coef, p = stats.spearmanr(subset["expression"],
                                          subset[cov_name])
                coef_str = "{:.2f}".format(coef)
                # p_str = p_value_to_symbol(p)
                p_str = "p = {:.2e}".format(p)

                # Plot.
                sns.regplot(x=cov_name, y="expression", data=subset,
                            scatter_kws={'facecolors': subset['value_hue'],
#                                         'edgecolors': subset['group_hue'],
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": color, "alpha": 0.75},
                            ax=ax
                            )

            # Add the text.
            ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
            ax.annotate(
                '{}: r = {} [{}]'.format(allele, coef_str, p_str),
                #xy=(0.03, 0.94 - ((i / 100) * 4)),
                xy=(0.03, label_pos[genotype]),
                xycoords=ax.transAxes,
                color=color,
                alpha=0.75,
                fontsize=12,
                fontweight='bold')

        ax.text(0.5, 1.06,
                '{} {}-eQTL Interaction with {} '.format(hgnc_name,
                                                         eqtl_type,
                                                         cov_name),
                fontsize=18, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                'SNPName: {}  ProbeName: {}  '
                'Z-score: {:.2f}'.format(snp_name, probe_name, zscore),
                fontsize=12, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel('{} ({}) expression'.format(probe_name, hgnc_name),
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(cov_name,
                      fontsize=14,
                      fontweight='bold')

        # ax.axvline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)
        # ax.axhline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}_inter_eqtl_{}_{}_{}_{}.{}".format(
                                     count,
                                     snp_name,
                                     probe_name,
                                     hgnc_name,
                                     cov_name,
                                     extension)))
        plt.close()

    @staticmethod
    def plot_box(snp_name, probe_name, hgnc_name, eqtl_type, df, cov_name,
                 zscore, count, allele_map, outdir, sex_colormap, extension):
        palette = None
        if cov_name == "SEX":
            df[cov_name] = df[cov_name].map({0: "Male", 1: "Female", -1: np.nan})
            df.dropna(inplace=True)
            palette = sex_colormap

        # Prepare the figure.
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        sns.swarmplot(x="group", y="expression", hue=cov_name, data=df,
                      dodge=True, color=".20", edgecolor="#808080", size=2,
                      ax=ax)
        sns.boxplot(x="group", y="expression", hue=cov_name, dodge=True,
                    zorder=-1, boxprops=dict(alpha=.8), showfliers=False,
                    data=df, palette=palette, ax=ax)

        # plt.setp(ax.artists, edgecolor='k', facecolor='w')
        # plt.setp(ax.lines, color='k')

        # Set the other aesthetics.
        ax.set_xticks(range(3))
        ax.set_xticklabels([allele_map[0.0], allele_map[1.0], allele_map[2.0]])

        ax.text(0.5, 1.06,
                '{} {}-eQTL Interaction with {} '.format(hgnc_name,
                                                         eqtl_type,
                                                         cov_name),
                fontsize=18, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                'SNPName: {}  ProbeName: {}  '
                'Z-score: {:.2f}'.format(snp_name, probe_name, zscore),
                fontsize=12, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel('{} ({}) expression'.format(probe_name, hgnc_name),
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel('',
                      fontsize=14,
                      fontweight='bold')
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=10)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2],
                  loc='center left', bbox_to_anchor=(1, 0.5))

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}_inter_eqtl_{}_{}_{}_{}.{}".format(
                                     count,
                                     snp_name,
                                     probe_name,
                                     hgnc_name,
                                     cov_name,
                                     extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL matrix shape: {}".format(self.eqtl_df.shape))
        print("  > Genotype matrix shape: {}".format(self.geno_df.shape))
        print("  > Alleles matrix shape: {}".format(self.alleles_df.shape))
        print("  > Expression matrix shape: {}".format(self.expr_df.shape))
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
