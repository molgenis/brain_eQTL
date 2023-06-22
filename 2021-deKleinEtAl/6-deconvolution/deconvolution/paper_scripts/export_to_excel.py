#!/usr/bin/env python3

"""
File:         export_to_excel.py
Created:      2022/04/06
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

# Third party imports.
import os.path

import pandas as pd

# Local application imports.


"""
Syntax:
./export_to_excel.py -h

./export_to_excel.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -n no_ENA_0_PCs no_ENA_100_PCs no_ENA_no_AMPAD_0_PCs no_ENA_no_AMPAD_80_PCs \
    -o 2022-04-06-CortexEUR-and-AFR-noENA-trans-eQTLs-DeconQTL
    
################################################################################
    
./export_to_excel.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/paper_scripts/sn_replication/trans/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/single_nucleus_replication.txt.gz /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/paper_scripts/sn_replication/trans/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/single_nucleus_replication.txt.gz /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/paper_scripts/sn_replication/trans/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/single_nucleus_replication.txt.gz /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/paper_scripts/sn_replication/trans/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/single_nucleus_replication.txt.gz \
    -n no_ENA_0_PCs no_ENA_100_PCs no_ENA_no_AMPAD_0_PCs no_ENA_no_AMPAD_80_PCs \
    -o 2022-04-06-CortexEUR-and-AFR-noENA-trans-eQTLs-SingleNucleusReplication
    
################################################################################

"""

# Metadata
__program__ = "Export to Excel"
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
        self.data = getattr(arguments, 'data')
        self.name = getattr(arguments, 'name')
        self.outputfile = getattr(arguments, 'outputfile')

        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), "export_to_excel")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-d",
                            "--data",
                            nargs="*",
                            type=str,
                            required=True,
                            help="The data input path(s)")
        parser.add_argument("-n",
                            "--name",
                            nargs="*",
                            type=str,
                            required=True,
                            help="The sheet names.")
        parser.add_argument("-o",
                            "--outputfile",
                            type=str,
                            required=True,
                            help="The output filename.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        with pd.ExcelWriter(os.path.join(self.outdir, '{}.xlsx'.format(self.outputfile))) as writer:
            for sheet_name, inpath in zip(self.name, self.data):
                print("{}: {}".format(sheet_name.replace("_", " "), inpath))

                df = pd.read_csv(inpath,
                                 sep="\t",
                                 header=0,
                                 index_col=None)

                df.to_excel(writer,
                            sheet_name=sheet_name.replace("_", " "),
                            na_rep="NA",
                            index=False,
                            header=True)

                print("")

        print("Saved '{}.xlsx'.".format(self.outputfile))

    def print_arguments(self):
        print("Arguments:")
        print("  > Data: {}".format(self.data))
        print("  > Name: {}".format(self.name))
        print("  > Output filename: {}".format(self.outputfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
