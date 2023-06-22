#!/usr/bin/env python3

"""
File:         plot_decon_eqtl_pvalues.py
Created:      2022/02/15
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
import math
import glob
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

"""
Syntax:
./plot_decon_eqtl_pvalues.py \
    -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron \
    -o 2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron \
    -e png pdf
    
./plot_decon_eqtl_pvalues.py \
    -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected \
    -o 2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected \
    -e png pdf
"""

# Metadata
__program__ = "Plot Decon-eQTL Nominal P-values"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.indir = getattr(arguments, 'indir')
        self.outfile = getattr(arguments, 'outfile')
        self.extensions = getattr(arguments, 'extensions')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot_decon_eqtl_pvalues')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
            "Astrocyte": "#D55E00",
            "EndothelialCell": "#CC79A7",
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Microglia": "#E69F00",
            "Oligodendrocyte": "#009E73",
            "OtherNeuron": "#2690ce"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-id",
                            "--indir",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of the input path.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the output file")
        parser.add_argument("-e",
                            "--extensions",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading nominal p-value data")
        decon_df = self.load_file(os.path.join(self.indir, "deconvolutionResults.txt.gz"), header=0, index_col=0)
        nominal_pvalues_df = decon_df.loc[:, [x for x in decon_df.columns if x.endswith("_pvalue")]]
        nominal_pvalues_df.columns = [x.replace("_pvalue", "") for x in nominal_pvalues_df.columns]
        print(nominal_pvalues_df)
        cell_types = list(nominal_pvalues_df.columns)
        print("\tcell types: {}".format(", ".join(cell_types)))

        print("Plotting nominal p-value data")
        self.plot(df=nominal_pvalues_df,
                  cell_types=cell_types,
                  xlabel="nominal p-value",
                  appendix="_nominal_pvalues",
                  title="Decon-eQTL nominal p-values")

        del decon_df, nominal_pvalues_df

        print("Loading permutation p-value data")
        perm_pvalues_m_list = []
        perm_pvalues_inpaths = glob.glob(os.path.join(self.indir, "permutation_pvalues_*"))
        perm_pvalues_inpaths.sort(key=self.natural_keys)
        for perm_pvalues_inpath in perm_pvalues_inpaths:
            perm_pvalues_m_list.append(self.load_matrix(perm_pvalues_inpath))
        perm_pvalues_m = np.dstack(perm_pvalues_m_list)
        del perm_pvalues_m_list
        print("\tShape: {}".format(perm_pvalues_m.shape))

        total_n = 0
        total_p_et_one = 0
        for i, ct in enumerate(cell_types):
            data = perm_pvalues_m[:, i, :].flatten()
            n_values = np.size(data)
            p_bt_ninefive = np.sum(data > 0.95)
            print("{}: {:,} / {:,} [{:.0f}%] permuted p-values are > 0.95".format(ct, p_bt_ninefive, n_values, (100 / n_values) * p_bt_ninefive))

            p_et_one = np.sum(data == 1)
            print("{}: {:,} / {:,} [{:.0f}%] permuted p-values are == 1".format(ct, p_et_one, n_values, (100 / n_values) * p_et_one))


            print("{}: median = {:.2f}".format(ct, np.median(data)))

            total_n += n_values
            total_p_et_one += p_et_one
        print("{:,} / {:,} [{:.0f}%] permuted p-values are == 1".format(total_p_et_one, total_n, (100 / total_n) * total_p_et_one))
        print(np.median(perm_pvalues_m.flatten()))

        # print("Plotting permuted p-value data")
        # self.plot(df=perm_pvalues_m,
        #           cell_types=cell_types,
        #           xlabel="permuted p-value",
        #           appendix="_permuted_pvalues",
        #           title="Decon-eQTL permuted p-values")

        log10_perm_pvalues_m = -1 * np.log10(perm_pvalues_m)
        self.plot(df=log10_perm_pvalues_m,
                  cell_types=cell_types,
                  xlabel="permuted -1log10 p-value",
                  appendix="_permuted_pvalues_log10",
                  title="Decon-eQTL -log10 permuted p-values")

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def natural_keys(text):
        return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', text)]

    @staticmethod
    def load_matrix(inpath):
        m = np.load(inpath)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      m.shape))
        return m

    def plot(self, df, cell_types, xlabel="", title="", appendix=""):
        nplots = len(cell_types)
        ncols = math.ceil(np.sqrt(nplots))
        nrows = math.ceil(nplots / ncols)
        exes_plots = (ncols * nrows) - nplots

        sns.set(rc={'figure.figsize': (ncols * 8, nrows * 6)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='none',
                                 sharey='all')

        data = None
        row_index = 0
        col_index = 0
        for i in range(ncols * nrows):
            print(i)
            if nrows == 1 and ncols == 1:
                ax = axes
            elif nrows == 1 and ncols > 1:
                ax = axes[col_index]
            elif nrows > 1 and ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            tmp_xlabel = ""
            if (row_index == (nrows - 1)) or (row_index == (nrows - 2) and col_index >= (ncols - exes_plots)):
                tmp_xlabel = xlabel
            tmp_ylabel = ""
            if col_index == 0:
                tmp_ylabel = "proportion"

            if i < len(cell_types):
                ct = cell_types[i]

                if len(df.shape) == 2:
                    data = df[ct]
                if len(df.shape) == 3:
                    data = df[:, i, :].flatten()

                sns.despine(fig=fig, ax=ax)

                sns.histplot(data=data,
                             stat="probability",
                             ax=ax,
                             kde=True,
                             kde_kws={"cut": 0},
                             color=self.palette[ct])

                ax.set_title(ct,
                             fontsize=24,
                             fontweight='bold')
                ax.set_ylabel(tmp_ylabel,
                              fontsize=20,
                              fontweight='bold')
                ax.set_xlabel(tmp_xlabel,
                              fontsize=20,
                              fontweight='bold')

                ax.annotate(
                    'N = {:,}'.format(np.size(data)),
                    xy=(0.03, 0.9),
                    xycoords=ax.transAxes,
                    color=self.palette[ct],
                    fontsize=20,
                    fontweight='bold'
                )
            else:
                # ax.set_xlabel(tmp_xlabel,
                #               fontsize=20,
                #               fontweight='bold')
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        # Add the main title.
        fig.suptitle(title,
                     fontsize=30,
                     color="#000000",
                     weight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}{}.{}".format(self.outfile, appendix, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Decon input directory: {}".format(self.indir))
        print("  > Output filename: {}".format(self.outfile))
        print("  > Extensions: {}".format(", ".join(self.extensions)))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
