#!/usr/bin/env python3

"""
File:         plot_cell_fraction_per_cohort.py
Created:      2021/02/06
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
import argparse
import glob
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
__program__ = "Plot Cell Fraction per Cohort"
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
./plot_cell_fraction_per_cohort.py -cf ../../2020-03-12-deconvolution/partial_deconvolution/ALL_TMM_LOG2/IHC_0CPM_LOG2_FILTERED_CC/deconvolution.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/trans/2020-05-26-Cortex-EUR-AFR-noENA-noPCA/gwasupdate/Iteration1-dsZscores -e pdf
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.gte_folder = getattr(arguments, 'gte_folder')
        self.extensions = getattr(arguments, 'extensions')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "cell_fraction_per_cohort")
        self.cohort_info = [
            ("GTE-AFR-AMPAD-MSBB-V2", "AMPAD-MSBB (AFR; n=TMP)"),
            ("GTE-EUR-AMPAD-MAYO-V2", "AMPAD-MAYO (EUR; n=TMP)"),
            ("GTE-EUR-AMPAD-MSBB-V2", "AMPAD-MSBB (EUR; n=TMP)"),
            ("GTE-EUR-AMPAD-ROSMAP-V2", "AMPAD-ROSMAP (EUR; n=TMP)"),
            ("GTE-EUR-BrainGVEX-V2", "BrainGVEX (EUR; n=TMP)"),
            ("GTE-AFR-CMC", "CMC (AFR; n=TMP)"),
            ("GTE-EUR-CMC", "CMC (EUR; n=TMP)"),
            ("GTE-AFR-CMC_HBCC_set1", "CMC_HBCC_set1 (AFR; n=TMP)"),
            ("GTE-AFR-CMC_HBCC_set2", "CMC_HBCC_set2 (AFR; n=TMP)"),
            ("GTE-EUR-CMC_HBCC_set2", "CMC_HBCC_set2 (EUR; n=TMP)"),
            ("GTE-EUR-CMC_HBCC_set3", "CMC_HBCC_set3 (EUR; n=TMP)"),
            ("GTE-EUR-GTEx", "GTEx (EUR; n=TMP)"),
            ("GTE-EUR-GVEX", "GVEX (EUR; n=TMP)"),
            ("GTE-AFR-LIBD_1M", "LIBD_1M (AFR; n=TMP)"),
            ("GTE-EUR-LIBD_1M", "LIBD_1M (EUR; n=TMP)"),
            ("GTE-AFR-LIBD_h650", "LIBD_h650 (AFR; n=TMP)"),
            ("GTE-EUR-LIBD_h650", "LIBD_h650 (EUR; n=TMP)"),
            ("GTE-EUR-NABEC-H550", "NABEC-H550 (EUR; n=TMP)"),
            ("GTE-EUR-NABEC-H610", "NABEC-H610 (EUR; n=TMP)"),
            ("GTE-EUR-UCLA_ASD", "UCLA_ASD (EUR; n=TMP)"),
        ]
        self.cohort_dict = {a: b for a, b in self.cohort_info}

        # Color map.

        self.colormap = {
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fracttions matrix")
        parser.add_argument("-gte",
                            "--gte_folder",
                            type=str,
                            required=True,
                            help="The path to the folder containg 'GTE-*' files.")
        parser.add_argument("-e",
                            "--extensions",
                            type=str,
                            nargs="+",
                            default=["png"],
                            choices=["eps", "pdf", "pgf", "png", "ps", "raw", "rgba", "svg", "svgz"],
                            help="The output file format(s), default: ['png']")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data.")
        cf_df = self.load_file(self.cf_path)

        gte_combined = None
        cohort_names = {}
        for fpath in glob.glob(os.path.join(self.gte_folder, 'GTE-*')):
            gte_filename = os.path.basename(fpath).split(".")[0]
            if gte_filename not in self.cohort_dict.keys():
                print("GTE file '{}' will be skipped.".format(gte_filename))
                continue
            gte_df = self.load_file(fpath, index_col=None, header=None)
            gte_df.columns = ["gene_id", "expr_id"]
            cohort = self.cohort_dict[gte_filename].replace("TMP", str(gte_df.shape[0]))
            gte_df["cohort"] = cohort
            if gte_combined is None:
                gte_combined = gte_df
            else:
                gte_combined = pd.concat([gte_combined, gte_df], axis=0)
            cohort_names[gte_filename] = cohort
        gte_combined.set_index("expr_id", inplace=True)

        print("### Step2 ###")
        print("Merging data.")
        df = cf_df.merge(gte_combined, left_index=True, right_index=True, how="right")
        print(df)
        print(df["cohort"].value_counts())

        print("### Step3 ###")
        print("Setting order.")
        order = []
        for (gte_filename, _) in self.cohort_info:
            order.append(cohort_names[gte_filename])

        print("### Step4 ###")
        print("Adding meta analysis.")
        df_meta = df.copy()
        df_meta["cohort"] = "Meta-Analysis"
        order.append("Meta-Analysis")
        df = pd.concat([df, df_meta], axis=0)
        print(df)
        print(order)

        print("### Step5 ###")
        print("Plotting data.")
        for celltype in cf_df.columns:
            color = self.colormap[celltype]
            self.create_catplot(df, y="cohort", x=celltype, order=order,
                                color=color)

        print("### Step6 ###")
        print("Calculating means")
        means = {}
        for celltype in cf_df.columns:
            data = []
            for cohort in order:
                cohort_data = df.loc[df["cohort"] == cohort, celltype]
                data.append(cohort_data.mean())
            s = pd.Series(data, index=order)
            means[celltype] = s

        print("### Step7 ###")
        print("Calculating t-tests")
        for celltype in cf_df.columns:
            t_test_data = []
            fold_change_data = []
            n_signif_per_cohort = {x: 0 for x in order}
            n_signif = 0
            n_total = 0
            for i, cohort1 in enumerate(order):
                n_signif_cohort_count = 0

                cohort1_data = df.loc[df["cohort"] == cohort1, celltype]
                t_test_row = []
                fold_change_row = []
                for j, cohort2 in enumerate(order):
                    cohort2_data = df.loc[df["cohort"] == cohort2, celltype]
                    t2, p2 = stats.ttest_ind(cohort1_data, cohort2_data)

                    n_total += 1
                    if p2 < 0.05:
                        n_signif_cohort_count += 1
                        n_signif += 1
                        print("\tSign. differnce between {} [{:.2f}] versus {} [{:.2f}]:\t{:.2f} {:.2e}".format(cohort1, means[celltype][cohort1], cohort2, means[celltype][cohort2], t2, p2))
                    t_test_row.append(p2)
                    fold_change_row.append(means[celltype][cohort1] / means[celltype][cohort2])
                n_signif_per_cohort[cohort1] = n_signif_cohort_count
                t_test_data.append(t_test_row)
                fold_change_data.append(fold_change_row)

            print("\tSignificant comparisons per cohort")
            for key, value in n_signif_per_cohort.items():
                print("\t\t{}\t{}".format(key, value))

            print("\t{}/{} tests are significant.".format(n_signif, n_total))
            ct_ttest_df = pd.DataFrame(t_test_data, index=order, columns=order)
            ct_ttest_df.to_csv(os.path.join(self.outdir, "{}_ttest_matrix.txt.gz".format(celltype)),
                               sep="\t", compression="gzip", header=True, index=True)

            ct_foldc_df = pd.DataFrame(fold_change_data, index=order, columns=order)
            print(ct_foldc_df)
            ct_foldc_df.to_csv(os.path.join(self.outdir, "{}_fold_change_matrix.txt.gz".format(celltype)),
                               sep="\t", compression="gzip", header=True, index=True)


    def create_catplot(self, df, x="x", y="y", order=None, color="#000000"):
        sns.set(rc={'figure.figsize': (12, 19)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(1, 2,
                                       gridspec_kw={"width_ratios": [0.1, 0.9]})
        ax1.axis('off')
        sns.despine(fig=fig, ax=ax2)
        sns.violinplot(x=x, y=y, data=df, order=order, color=color, ax=ax2)
        ax2.set_ylabel("")

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "x{}_y{}_violinplot.{}".format(x, y, extension)))
        plt.close()

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None, low_memory=True):
        if path.endswith(".pkl"):
            df = pd.read_pickle(path)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def print_arguments(self):
        print("Arguments:")
        print("  > Cell fraction file: {}".format(self.cf_path))
        print("  > Gene to Expression linking folder: {}".format(self.gte_folder))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
