"""
File:         simple_eqtl_effect.py
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
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.
from general.utilities import prepare_output_dir, p_value_to_symbol


class SimpleeQTLEffect:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'simple_eqtl_effect')
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
        colormap = dataset.get_colormap()

        # Create color map.
        self.group_color_map, self.value_color_map = self.create_color_map(colormap)

    def start(self):
        print("Plotting simple eQTL plots.")
        self.print_arguments()

        print("Iterating over eQTLs.")
        for i, (index, row) in enumerate(self.eqtl_df.iterrows()):
            # Extract the usefull information from the row.
            p_value = row["PValue"]
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

            # Determine the minor allele frequency.
            minor_allele_frequency = min(zero_geno_count, two_geno_count) / (
                        zero_geno_count + two_geno_count)

            # Add the color.
            data["round_geno"] = data["genotype"].round(2)
            data["value_hue"] = data["round_geno"].map(self.value_color_map)
            data["group_hue"] = data["group"].map(self.group_color_map)
            data.drop(["round_geno"], axis=1, inplace=True)

            # Plot a simple eQTL effect.
            self.plot(index, p_value, snp_name, probe_name, hgnc_name,
                      eqtl_type,
                      data, minor_allele, minor_allele_frequency, allele_map,
                      self.group_color_map, self.outdir, self.extension)

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
    def plot(i, p_value, snp_name, probe_name, hgnc_name, eqtl_type, df,
             minor_allele, minor_allele_frequency, allele_map, group_color_map,
             outdir, extension):
        """
        """
        # Calculate the correlation.
        coef, p = stats.spearmanr(df["genotype"],
                                  df["expression"])

        # Prepare the figure.
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        # Plot the scatter / box plot.
        sns.regplot(x="genotype", y="expression", data=df,
                    scatter=False,
                    # scatter_kws={'facecolors': df['value_hue'],
                    #              'edgecolors': "#808080"},
                    line_kws={"color": "#000000"},
                    ax=ax
                    )
        sns.boxplot(x="group", y="expression", data=df,
                    palette=group_color_map,
                    showfliers=False,
                    zorder=1,
                    #boxprops=dict(alpha=.3),
                    ax=ax)
        # plt.setp(ax.artists, edgecolor='k', facecolor='w')
        # plt.setp(ax.lines, color='k')

        # Set the other aesthetics.
        ax.set_xticks(range(3))
        ax.set_xticklabels([allele_map[0.0], allele_map[1.0], allele_map[2.0]])
        ax.text(0.5, 1.06,
                # '{} {}-eQTL [{}]'.format(hgnc_name, eqtl_type,
                #                          p_value_to_symbol(p_value)),
                '{} {}-eQTL [p = {:.2e}]'.format(hgnc_name, eqtl_type, p_value),
                fontsize=22, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                # 'r = {:.2f} [{}]    minor allele frequency '
                # '{} = {:.2f}'.format(coef,
                #                      p_value_to_symbol(p),
                #                      minor_allele,
                #                      minor_allele_frequency),
                'r = {:.2f} [p = {:.2e}]    minor allele frequency '
                '{} = {:.2f}'.format(coef,
                                     p,
                                     minor_allele,
                                     minor_allele_frequency),
                fontsize=14, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)
        ax.set_ylabel('{} ({}) expression'.format(probe_name, hgnc_name),
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel('SNP {}'.format(snp_name),
                      fontsize=14,
                      fontweight='bold')

        # Safe the plot.
        fig.savefig(os.path.join(outdir, "{}_{}_{}_{}.{}".format(i,
                                                                 snp_name,
                                                                 probe_name,
                                                                 hgnc_name,
                                                                 extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype matrix shape: {}".format(self.geno_df.shape))
        print("  > Alleles matrix shape: {}".format(self.alleles_df.shape))
        print("  > Expression matrix shape: {}".format(self.expr_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
