"""
File:         inter_eqtl_celltype_details.py
Created:      2020/05/20
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
import math
import os

# Third party imports.
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class IntereQTLCelltypeDetails:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_celltype_details')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.eqtl_df = dataset.get_eqtl_df()
        self.geno_df = dataset.get_geno_df()
        self.zscore_df = dataset.get_inter_cov_zscore_df()
        self.tvalue_df = dataset.get_inter_cov_tvalue_df()
        self.cellmap_methods = dataset.get_cellmap_methods()
        self.marker_genes = dataset.get_marker_genes()
        self.z_score_cutoff = dataset.get_significance_cutoff()
        self.colormap = dataset.get_colormap()

    def start(self):
        print("Plotting interaction eQTL radar plots.")
        self.print_arguments()

        methods = self.cellmap_methods
        methods.append((self.marker_genes, ""))

        print("Iterating over eQTLs.")
        for i, (index, row) in enumerate(self.eqtl_df.iterrows()):
            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]

            print("\tWorking on: {}\t{}\t{} [{}/{} "
                  "{:.2f}%]".format(snp_name, probe_name, hgnc_name,
                                    i + 1,
                                    self.eqtl_df.shape[0],
                                    (100 / self.eqtl_df.shape[0]) * (i + 1)))

            # Check if we need to flip the genotypes.
            genotype = self.geno_df.iloc[i, :]
            counts = genotype.value_counts()
            for x in [0.0, 1.0, 2.0]:
                if x not in counts:
                    counts.loc[x] = 0
            zero_geno_count = (counts[0.0] * 2) + counts[1.0]
            two_geno_count = (counts[2.0] * 2) + counts[1.0]
            flip = 1
            if two_geno_count > zero_geno_count:
                flip = -1

            # Prepare output directory.
            eqtl_outdir = os.path.join(self.outdir,
                                       "{}_{}_{}_{}".format(index, snp_name,
                                                            probe_name,
                                                            hgnc_name))
            prepare_output_dir(eqtl_outdir)

            # Iterate over the rows.
            for (prefix, suffix) in methods:
                name = prefix.replace("_", "") + suffix

                tvalues = self.tvalue_df.loc[
                          self.tvalue_df.index.str.startswith(prefix), :].copy()
                tvalues = tvalues.iloc[:, i]
                tvalues = tvalues * flip
                tvalues = tvalues.to_frame()

                zscores = self.zscore_df.loc[
                          self.zscore_df.index.str.startswith(prefix), :].copy()
                zscores = zscores.iloc[:, i].to_frame()

                df = tvalues.merge(zscores, left_index=True, right_index=True)
                df.columns = ["tvalue", "zscore"]
                df.index = ["{}".format(x.replace(prefix, "").replace(suffix, "")) for x in df.index]

                self.plot_forest(snp_name, probe_name, hgnc_name, name, df,
                                 self.z_score_cutoff, eqtl_outdir)

    @staticmethod
    def plot_forest(snp_name, probe_name, hgnc_name, method, data,
                    z_score_cutoff, outdir):
        hue_list = []
        for zscore in data["zscore"]:
            add = "not signif."
            if zscore > z_score_cutoff:
                add = "signif."
            hue_list.append(add)
        data["significant"] = hue_list
        colormap = {'signif.': '#6495ED', 'not signif.': '#808080'}

        abs_max = data["tvalue"].abs().max() * 1.1

        sns.set(rc={'figure.figsize': (8, len(data.index)*0.8)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        g = sns.stripplot(x="tvalue", y=data.index, hue="significant", data=data,
                          size=20,
                          dodge=False,
                          orient="h",
                          palette=colormap,
                          linewidth=1,
                          edgecolor="w",
                          jitter=0,
                          ax=ax)

        plt.title('{} Celltype Interactions'.format(hgnc_name),
                  fontsize=14,
                  weight='bold')
        g.set_ylabel('',
                     fontsize=12,
                     fontweight='bold')
        g.set_xlabel('t-value',
                     fontsize=12,
                     fontweight='bold')
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=12)

        ax.set_xlim(-1*abs_max, abs_max)

        plt.axvline(x=0, ls='-', c='black')

        if data["tvalue"].max() > abs_max:
            plt.axvline(x=abs(z_score_cutoff), ls='--', c='firebrick')
        if data["tvalue"].min() < -1 * abs_max:
            plt.axvline(x=-1 * abs(z_score_cutoff), ls='--', c='firebrick')

        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.tight_layout()
        fig.savefig(os.path.join(outdir, method + "_stripplot.png"))
        plt.close()

    @staticmethod
    def plot_radar(name, df, outdir):
        max_val = math.ceil(df.values.max())
        steps = [round(x, 2) for x in np.linspace(0, max_val, 4)]

        categories = list(df.columns)
        N = len(categories)

        angles = [n / float(N) * 2 * math.pi for n in range(N)]
        angles += angles[:1]

        ax = plt.subplot(111, polar=True)

        ax.set_theta_offset(math.pi / 2)
        ax.set_theta_direction(-1)

        plt.xticks(angles[:-1], categories)
        ax.set_rlabel_position(0)
        plt.yticks(steps[1:], steps[1:], color="grey", size=7)
        plt.ylim(0, max_val)

        colors = ['b', 'r']
        for i, (index, row) in enumerate(df.iterrows()):
            values = row.values.flatten().tolist()
            values += values[:1]
            ax.plot(angles, values, linewidth=1, linestyle='solid', label=index)
            ax.fill(angles, values, colors[i], alpha=0.1)

        plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))

        plt.gcf().canvas.draw()
        labels = []
        half_step = 360.0 / N / 2
        angles = [0.0, 8 * half_step, half_step, 9 * half_step, 2 * half_step]
        print(angles)
        for label, angle in zip(ax.get_xticklabels(), angles):
            print(label, angle)
            x, y = label.get_position()
            lab = ax.text(x, y, label.get_text(),
                          transform=label.get_transform(),
                          ha=label.get_ha(), va=label.get_va())
            lab.set_rotation(angle)
            labels.append(lab)
        ax.set_xticklabels([])

        plt.savefig(os.path.join(outdir, name + "_radarplot.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL matrix shape: {}".format(self.eqtl_df.shape))
        print("  > Genotype matrix shape: {}".format(self.geno_df.shape))
        print("  > Z-score matrix shape: {}".format(self.zscore_df.shape))
        print("  > T-value matrix shape: {}".format(self.tvalue_df.shape))
        print("  > CellMap Methods: {}".format(self.cellmap_methods))
        print("  > Marker Genes: {}".format(self.marker_genes))
        print("  > Z-score cutoff: {}".format(self.z_score_cutoff))
        print("  > Output directory: {}".format(self.outdir))
        print("")
