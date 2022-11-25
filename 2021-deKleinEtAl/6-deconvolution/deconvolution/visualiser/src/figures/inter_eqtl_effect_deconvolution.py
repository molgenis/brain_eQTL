"""
File:         inter_eqtl_effect_deconvolution.py
Created:      2020/03/17
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
from colour import Color
import os

# Third party imports.
from scipy import stats
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir, p_value_to_symbol


class IntereQTLEffectDeconvolution:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_effect_deconvolution')
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
        self.celltypes = dataset.get_celltypes()
        self.cellmap_methods = dataset.get_cellmap_methods()
        self.marker_genes = dataset.get_marker_genes()

        # Create color map.
        self.group_color_map, self.value_color_map = self.create_color_map()

    def start(self):
        print("Plotting interaction eQTL plots for deconvolution methods")
        self.print_arguments()

        # Determine the columns of the deconvolution rows from the
        # covariate matrix.
        deconvolution_indices = []
        for index in self.cov_df.index:
            for decon_prefix in [x[0] for x in self.cellmap_methods]:
                if index.startswith(decon_prefix) or index.startswith(
                        self.marker_genes):
                    deconvolution_indices.append(index)
                    break

        # Get the covariates of the marker genes.
        decon_df = self.cov_df.loc[deconvolution_indices, :]

        print("Iterating over eQTLs.")
        for i, (index, row) in enumerate(self.eqtl_df.iterrows()):
            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]

            print("\tWorking on: {}\t{}\t{}\t[{}/{}: "
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
            data.drop(["round_geno"], axis=1, inplace=True)

            # Get the interaction zscores
            interaction_effect = self.inter_df.iloc[:, i].to_frame()
            interaction_effect = interaction_effect.loc[deconvolution_indices,
                                 :]
            interaction_effect.columns = ["zscore"]

            self.plot(snp_name, probe_name, hgnc_name, data, decon_df,
                      interaction_effect, self.celltypes, index, allele_map,
                      self.group_color_map, self.outdir, self.extension)

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
    def plot(snp_name, probe_name, hgnc_name, data, decon_df,
             zscores, celltypes, count, allele_map, group_color_map,
             outdir, extension):
        """
        """
        # Calculate number of rows / columns.
        ncols = int(len(decon_df.index) / len(celltypes))
        nrows = int(len(celltypes))

        sns.set(rc={'figure.figsize': (12 * ncols, 9 * nrows)})
        sns.set_style("ticks")
        fig = plt.figure()
        grid = fig.add_gridspec(ncols=ncols,
                                nrows=nrows)

        # calculate axis limits.
        min = data["expression"].min() * 1.1
        max = data["expression"].max() * 1.5

        # Plot the markers.
        row_index = 0
        col_index = 0
        for celltype in celltypes:
            method_celltype = celltype
            if method_celltype == "Microglia":
                method_celltype = "Macrophage"

            celltype_df = decon_df.loc[decon_df.index.str.contains(celltype) |
                                       decon_df.index.str.contains(
                                           method_celltype), :]

            for cov_name, cov_data in celltype_df.iterrows():
                # Combine the data.
                df = data.copy()
                df = df.merge(cov_data, left_index=True, right_index=True)

                # Creat the subplot.
                ax = fig.add_subplot(grid[row_index, col_index])
                sns.despine(fig=fig, ax=ax)

                # Plot the groups.
                for i, genotype in enumerate([0.0, 1.0, 2.0]):
                    subset = df.loc[df["round_geno"] == genotype, :].copy()
                    color = group_color_map[genotype]
                    allele = allele_map[genotype]

                    coef_str = "NA"
                    p_str = "NA"
                    if len(subset.index) > 1:
                        # Calculate the correlation.
                        coef, p = stats.spearmanr(subset["expression"],
                                                  subset[cov_name])
                        coef_str = "{:.2f}".format(coef)
                        #p_str = p_value_to_symbol(p)
                        p_str = "p = {:.2e}".format(p)

                        # Plot.
                        sns.regplot(x=cov_name, y="expression", data=subset,
                                    scatter_kws={
                                        'facecolors': subset['value_hue'],
#                                        'edgecolors': subset['group_hue'],
                                        'linewidth': 0,
                                        'alpha': 0.75},
                                    line_kws={"color": color, "alpha": 0.75},
                                    ax=ax
                                    )

                    # Add the text.
                    ax.set(ylim=(min, max))
                    ax.annotate(
                        '{}: r = {} [{}]'.format(allele, coef_str, p_str),
                        xy=(0.03, 0.94 - ((i / 100) * 4)),
                        xycoords=ax.transAxes,
                        color=color,
                        alpha=0.75,
                        fontsize=18,
                        fontweight='bold')

                ax.text(0.5, 1.06,
                        cov_name,
                        fontsize=26, weight='bold', ha='center', va='bottom',
                        transform=ax.transAxes)
                ax.text(0.5, 1.02,
                        '[z-score: {:.2f}]'.format(
                            zscores.at[cov_name, "zscore"]),
                        fontsize=20, alpha=0.75, ha='center', va='bottom',
                        transform=ax.transAxes)

                ax.set_ylabel(
                    '{} ({}) expression'.format(probe_name, hgnc_name),
                    fontsize=18,
                    fontweight='bold')
                ax.set_xlabel(cov_name,
                              fontsize=18,
                              fontweight='bold')

                col_index += 1
            row_index += 1
            col_index = 0

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}_inter_eqtl_{}_{}_{}_deconvolution.{}".format(
                                     count,
                                     snp_name,
                                     probe_name,
                                     hgnc_name,
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
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > CellMap Methods: {}".format(self.cellmap_methods))
        print("  > Marker Genes: {}".format(self.marker_genes))
        print("  > Output directory: {}".format(self.outdir))
        print("")
