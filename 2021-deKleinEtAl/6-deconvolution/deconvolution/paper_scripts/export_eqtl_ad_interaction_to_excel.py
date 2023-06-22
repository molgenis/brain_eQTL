#!/usr/bin/env python3

"""
File:         export_eqtl_ad_interaction_to_excel.py
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
./export_eqtl_ad_interaction_to_excel.py -h
    
./export_eqtl_ad_interaction_to_excel.py \
    -i /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/interaction_mapper/ \
    -n 2022-03-31-CortexEUR-and-AFR-noENA-trans-NPCs-NegativeToZero-DatasetAndRAMCorrected

"""

# Metadata
__program__ = "Export eQTL ~ AD interaction to Excel"
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

        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), "export_eqtl_ad_interaction_to_excel")
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

        appendix_dict = {"": "default",
                         "-DatasetCorrected": "dataset corrected",
                         "-AMPAD": "only AMP-AD samples"}
        with pd.ExcelWriter(os.path.join(self.outdir, '{}.xlsx'.format(self.name))) as writer:
            summary_stats = []
            for n_pcs in [0, 100]:
                df_list = []
                for appendix in ["", "-DatasetCorrected", "-AMPAD"]:
                    fpath = os.path.join(self.indir, self.name.replace("NPCs", "{}PCs".format(n_pcs)) + appendix, "Alzheimerdisease_InteractionResults.txt.gz")
                    print(n_pcs, appendix, fpath)

                    df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
                    df.index = df["SNPName"] + "_" + df["ProbeName"]
                    if df_list:
                        df = df[["FDR"]].copy()
                        df.columns = ["{} FDR".format(appendix_dict[appendix])]
                    else:
                        df = df[["SNPName", "ProbeName", "N", "FDR"]].copy()
                        df.columns = ["SNPName", "ProbeName", "N", "{} FDR".format(appendix_dict[appendix])]
                    df.dropna(inplace=True)
                    df_list.append(df)

                    summary_stats.append(["no ENA", n_pcs, appendix_dict[appendix], df.shape[0], (df["{} FDR".format(appendix_dict[appendix])] < 0.05).sum()])

                df = pd.concat(df_list, axis=1)
                df.to_excel(writer,
                            sheet_name="no ENA {}PCs".format(n_pcs),
                            na_rep="NA",
                            index=False,
                            header=True)

                df = pd.DataFrame(summary_stats, columns=["trans-eQTL dataset", "#PCs removed", "sub-analysis", "#eQTLs tested", "#ieQTLs (FDR <0.05)"])
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
