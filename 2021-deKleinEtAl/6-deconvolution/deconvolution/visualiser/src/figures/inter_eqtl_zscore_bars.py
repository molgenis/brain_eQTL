"""
File:         inter_eqtl_zscore_bars.py
Created:      2020/03/16
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
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from colour import Color

# Local application imports.
from general.utilities import prepare_output_dir


class IntereQTLZscoreBars:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_zscore_bars')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

        # Extract the required data.
        print("Loading data")
        self.eqtl_df = dataset.get_eqtl_df()
        self.inter_df = dataset.get_inter_cov_zscore_df()
        self.z_score_cutoff = dataset.get_significance_cutoff()

    def start(self):
        print("Plotting interaction eQTL z-score barplots.")
        self.print_arguments()

        colormap = self.create_color_map(self.z_score_cutoff)

        print("Iterating over eQTLs.")
        for i, (index, row) in enumerate(self.eqtl_df.iterrows()):
            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            eqtl_type = row["CisTrans"]

            print("\tWorking on: {}\t{}\t{} [{}/{} "
                  "{:.2f}%]".format(snp_name, probe_name, hgnc_name,
                                    i + 1,
                                    self.eqtl_df.shape[0],
                                    (100 / self.eqtl_df.shape[0]) * (i + 1)))

            # Check if the SNP has an interaction effect.
            interaction_effect = self.inter_df.iloc[:, i].to_frame()
            # if interaction_effect.columns[0] in ["7:144177389:rs6464583:A_G_G", "9:32886201:rs10971181:T_C_C", "22:16717912:rs2078647:A_G_G", "15:33985176:rs11856529:A_G_G"]:

            interaction_effect.reset_index(inplace=True)
            interaction_effect.columns = ["index", "zscore"]
            interaction_effect = interaction_effect.reindex(
                interaction_effect["zscore"].sort_values().index)
            interaction_effect["color"] = [colormap[round(x, 1)] for x in
                                           interaction_effect["zscore"]]

            eqtl_interaction_outdir = os.path.join(self.outdir,
                                                   "{}_{}_{}_{}".format(index,
                                                                        snp_name,
                                                                        probe_name,
                                                                        hgnc_name))
            if not os.path.exists(eqtl_interaction_outdir):
                os.makedirs(eqtl_interaction_outdir)

            self.plot(index, snp_name, probe_name, hgnc_name, eqtl_type,
                      self.z_score_cutoff, interaction_effect,
                      eqtl_interaction_outdir, self.extension)
            self.plot(index, snp_name, probe_name, hgnc_name, eqtl_type,
                      self.z_score_cutoff, interaction_effect,
                      eqtl_interaction_outdir, self.extension, positive=True)

    def create_color_map(self, signif_cutoff):
        min_value = -8.3
        max_value = 385

        blue_values = [x / 10 for x in
                       range(int(min_value * 10), int(signif_cutoff * -10), 1)]
        blue_colors = [x.rgb for x in list(
            Color("#6282EA").range_to(Color("#B9D0F9"), len(blue_values)))]
        blue_colors = self.check_for_neg_rgv(blue_colors)

        neg_black_values = [x / 10 for x in
                            range(int(signif_cutoff * -10), 0, 1)]
        neg_black_colors = [x.rgb for x in list(
            Color("#000000").range_to(Color("#DCDCDC"), len(neg_black_values)))]
        neg_black_colors = self.check_for_neg_rgv(neg_black_colors)

        pos_black_values = [x / 10 for x in
                            range(1, int(signif_cutoff * 10) + 2, 1)]
        pos_black_colors = [x.rgb for x in list(
            Color("#DCDCDC").range_to(Color("#000000"), len(pos_black_values)))]
        pos_black_colors = self.check_for_neg_rgv(pos_black_colors)

        red_values = [x / 10 for x in
                      range(int(signif_cutoff * 10) + 2, int(max_value * 10),
                            1)]
        red_colors = [x.rgb for x in list(
            Color("#F5C4AC").range_to(Color("#DC5F4B"), len(red_values)))]
        red_colors = self.check_for_neg_rgv(red_colors)

        values = blue_values + neg_black_values + [
            0.0] + pos_black_values + red_values
        colors = blue_colors + neg_black_colors + [Color(
            rgb=(0.863, 0.863, 0.863)).rgb] + pos_black_colors + red_colors
        value_color_map = {x: y for x, y in zip(values, colors)}
        return value_color_map

    @staticmethod
    def check_for_neg_rgv(colors):
        copy = colors.copy()
        for i, (r, g, b) in enumerate(copy):
            if r < 0:
                r = 0
            if g < 0:
                g = 0
            if b < 0:
                b = 0
            colors[i] = (r, g, b)

        return colors

    @staticmethod
    def plot(i, snp_name, probe_name, hgnc_name, eqtl_type, z_score_cutoff,
             df, outdir, extension, positive=False):

        subtitle_str = "+/-"
        file_suffix = "all"
        title_distance = 1.5e-2
        if positive:
            title_distance = 4e-2
            subtitle_str = "+"
            file_suffix = "positive"
            df = df.loc[df["zscore"] > 0, :]

        sns.set(rc={'figure.figsize': (12, (.2 * (len(df.index))))})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        g = sns.barplot(x="zscore", y="index", palette=df["color"], data=df,
                        orient="h")
        g.text(0.5, 1 + (7.5e-5 * len(df.index)) + title_distance,
               '{} {}-eQTL Interactions'.format(hgnc_name, eqtl_type),
               fontsize=22, weight='bold', ha='center', va='bottom',
               transform=ax.transAxes)
        g.text(0.5, 1 + (7.5e-5 * len(df.index)),
               'z-score cutoff = {} {:.2f}'.format(subtitle_str,
                                                   abs(z_score_cutoff)),
               fontsize=20, alpha=0.75, ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('covariate',
                     fontsize=18,
                     fontweight='bold')
        g.set_xlabel('z-score',
                     fontsize=18,
                     fontweight='bold')
        plt.axvline(x=z_score_cutoff, ls='--', c='black')

        if not positive:
            plt.axvline(x=-1 * z_score_cutoff, ls='--', c='black')
        ax.tick_params(labelsize=10)
        ax.set_yticks(range(len(df.index)))
        ax.set_yticklabels(df["index"])
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}_inter_eqtl_bars"
                                 "_{}_{}_{}_{}.{}".format(i,
                                                          snp_name,
                                                          probe_name,
                                                          hgnc_name,
                                                          file_suffix,
                                                          extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL matrix shape: {}".format(self.eqtl_df.shape))
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
