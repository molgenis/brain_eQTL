"""
File:         simple_eqtl_effect.py
Created:      2020/03/16
Last Changed: 2020/04/20
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
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'simple_eqtl_effect')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.eqtl_df = dataset.get_eqtl_df()
        self.geno_df = dataset.get_geno_df()
        self.expr_df = dataset.get_expr_df()
        self.alleles_df = dataset.get_alleles_df()

        # Create color map.
        self.group_color_map, self.value_color_map = self.create_color_map()

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
            (alleles, minor_allele) = self.alleles_df.iloc[i, :]

            # Determine the genotype order.
            first_allele = snp_name.split(":")[-1].split("_")[0]
            second_allele = snp_name.split(":")[-1].split("_")[1]

            # Check if the major / minor allele gentoypes are correct.
            counts = data["group"].value_counts()
            minor_genotype = counts.idxmin()

            # Determine the minor allele frequency.
            minor_allele_frequency = ((counts[minor_genotype] * 2) +
                                      counts[1.0]) / (data.shape[0] * 2)

            # Flip the alleles.
            if ((minor_allele == first_allele) and not (
                    minor_genotype == 0.0)) \
                    or ((minor_allele == second_allele) and not (
                    minor_genotype == 2.0)):
                # Flip the genotypes in order to get the genotype labels
                # correct.
                data["genotype"] = 2.0 - data["genotype"]
                data["group"] = 2.0 - data["group"]
            allele_map = {0.0: "{}/{}".format(first_allele, first_allele),
                          1.0: "{}/{}".format(first_allele, second_allele),
                          2.0: "{}/{}".format(second_allele, second_allele)}
            data["alleles"] = data["group"].map(allele_map)

            # Add the color.
            data["round_geno"] = data["genotype"].round(2)
            data["value_hue"] = data["round_geno"].map(self.value_color_map)
            data["group_hue"] = data["group"].map(self.group_color_map)
            data.drop(["round_geno"], axis=1, inplace=True)

            # Plot a simple eQTL effect.
            self.plot(i, p_value, snp_name, probe_name, hgnc_name, eqtl_type,
                      data, minor_allele, minor_allele_frequency,
                      first_allele, second_allele,
                      self.outdir)

    @staticmethod
    def create_color_map():
        """
        """
        palette = list(Color("#ABDCA2").range_to(Color("#CBE9C5"), 50)) + \
                  list(Color("#B1C2E1").range_to(Color("#89A3D1"), 50)) + \
                  list(Color("#89A3D1").range_to(Color("#B1C2E1"), 50)) + \
                  list(Color("#F6BCAD").range_to(Color("#F08C72"), 51))
        colors = [str(x).upper() for x in palette]
        values = [x / 100 for x in list(range(201))]
        group_color_map = {0.0: "#ABDCA2", 1.0: "#89A3D1", 2.0: "#F08C72"}
        value_color_map = {}
        for val, col in zip(values, colors):
            value_color_map[val] = col
        return group_color_map, value_color_map

    @staticmethod
    def plot(i, p_value, snp_name, probe_name, hgnc_name, eqtl_type, df,
             minor_allele, minor_allele_frequency, first_allele,
             second_allele, outdir):
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
                    scatter_kws={'facecolors': df['value_hue'],
                                 'edgecolors': "#808080"},
                    ax=ax
                    )
        sns.boxplot(x="group", y="expression", data=df,
                    showfliers=False,
                    zorder=-1,
                    boxprops=dict(alpha=.3),
                    ax=ax)
        plt.setp(ax.artists, edgecolor='k', facecolor='w')
        plt.setp(ax.lines, color='k')

        # Set the other aesthetics.
        ax.set_xticks(range(3))
        ax.set_xticklabels(["{}/{}".format(first_allele, first_allele),
                            "{}/{}".format(first_allele, second_allele),
                            "{}/{}".format(second_allele, second_allele)])
        ax.text(0.5, 1.06,
                '{} {}-eQTL [{}]'.format(hgnc_name, eqtl_type,
                                         p_value_to_symbol(p_value)),
                fontsize=22, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                'r = {:.2f} [{}]    minor allele frquency '
                '{} = {:.2f}'.format(coef,
                                     p_value_to_symbol(p),
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
        fig.savefig(os.path.join(outdir, "{}_{}_{}_{}.png".format(i,
                                                                  snp_name,
                                                                  probe_name,
                                                                  hgnc_name)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype matrix shape: {}".format(self.geno_df.shape))
        print("  > Alleles matrix shape: {}".format(self.alleles_df.shape))
        print("  > Expression matrix shape: {}".format(self.expr_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")