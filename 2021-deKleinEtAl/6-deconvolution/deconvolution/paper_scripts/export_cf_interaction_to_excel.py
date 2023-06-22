#!/usr/bin/env python3

"""
File:         export_cf_interaction_to_excel.py
Created:      2022/05/16
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
import argparse
import glob
import os

# Third party imports.
import os.path

import pandas as pd

# Local application imports.


"""
Syntax:
./export_cf_interaction_to_excel.py -h
    
./export_cf_interaction_to_excel.py \
    -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/ \
    -n 2022-03-31-CortexEUR-and-AFR-subset-trans-NPCs-NegativeToZero-DatasetAndRAMCorrected

"""

# Metadata
__program__ = "Export Cell Fraction interaction to Excel"
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
        self.name = getattr(arguments, 'name')

        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), "export_cf_interaction_to_excel")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-i",
                            "--indir",
                            type=str,
                            required=True,
                            help="The data input directory")
        parser.add_argument("-n",
                            "--name",
                            type=str,
                            required=True,
                            help="The filename.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        with pd.ExcelWriter(os.path.join(self.outdir, '{}.xlsx'.format(self.name))) as writer:
            summary_stats = []
            for subset, n_pcs in (["noENA", 0], ["noENA", 100], ["noENA-noAMPAD", 0], ["noENA-noAMPAD", 80]):
                fpath = os.path.join(self.indir, self.name.replace("subset", subset).replace("NPCs", "{}PCs".format(n_pcs)), "merged_decon_results.txt.gz")
                print(fpath)

                df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)

                for col in [col for col in df.columns if "BH-FDR" in col]:
                    summary_stats.append([subset.replace("no", "no ").replace("-", " "), n_pcs, col.replace(" BH-FDR", ""), df.shape[0], (df[col] < 0.05).sum()])

                df.to_excel(writer,
                            sheet_name="{} {}PCs".format(subset.replace("no", "no ").replace("-", " "), n_pcs),
                            na_rep="NA",
                            index=False,
                            header=True)

            df = pd.DataFrame(summary_stats,
                              columns=["trans-eQTL dataset", "#PCs removed", "cell type", "#eQTLs tested", "#ieQTLs (FDR <0.05)"])
            df.to_excel(writer,
                        sheet_name="Summary stats",
                        na_rep="NA",
                        index=False,
                        header=True)

        print("Saved '{}.xlsx'.".format(self.name))

    def print_arguments(self):
        print("Arguments:")
        print("  > Indir: {}".format(self.indir))
        print("  > Name: {}".format(self.name))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
