#!/usr/bin/env python3

"""
File:         preprocess_matrix.py
Created:      2020/11/03
Last Changed: 2020/12/01
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
import glob
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Preprocess matrix"
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
        self.infolder = "/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/2020-10-22-MVochteloo-Copy/"
        self.sample_trans_path = "/groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/meta/ROSMAP_IDkey.csv"
        self.cis_eqtls_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/Iteration1/eQTLProbesFDR0.05-ProbeLevel.txt.gz"
        self.trans_eqtls_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/trans/2020-05-26-Cortex-EUR-AFR-noENA-noPCA/Iteration1/eQTLs-crossMappingEQTLsRemoved-FDR0.05.txt.gz"
        self.outdir = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq"

    def start(self):
        sample_trans_df = pd.read_csv(self.sample_trans_path, sep=",")[["projid", "wgs_id"]]
        sample_trans_df.dropna(axis=0, inplace=True)
        sample_trans_dict = dict(zip(sample_trans_df.iloc[:, 0], sample_trans_df.iloc[:, 1]))

        samples = None
        for inpath in glob.glob(os.path.join(self.infolder, "*.tsv")):
            filename = os.path.basename(inpath.replace(".tsv", ".txt"))
            split_filename = filename.split("_")
            outpath = os.path.join(self.outdir, "_".join([split_filename[0].upper(), split_filename[1]]))
            df = pd.read_csv(inpath, sep="\t", header=0, index_col=0)

            df.columns = [int(x) for x in df.columns]

            col_mask, _ = self.translate_and_filter(df.columns, sample_trans_dict)
            df = df.loc[:, col_mask]
            df.index.name = "-"

            if samples is not None and list(df.columns) != samples:
                print("Samples do not match over all matrices!")
            samples = list(df.columns)

            print(outpath, df.shape)
            df.to_csv(outpath, header=True, index=True, sep="\t")

        data = []
        for sample in samples:
            wgs_id = None
            if sample in sample_trans_dict:
                wgs_id = sample_trans_dict[sample]

            data.append([wgs_id, sample])
        df = pd.DataFrame(data)
        print(df, df.shape)
        df.to_csv(os.path.join(self.outdir, "ROSMAP-scRNAseq-genometoexpressioncoupling.txt"), header=False, index=False, sep="\t")

        snps = set()
        for eqtl_path, eqtl_type in zip([self.cis_eqtls_path, self.trans_eqtls_path], ["cis", "trans"]):
            outdir = os.path.join(self.outdir, "{}_100Perm".format(eqtl_type))
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            eqtl_df = pd.read_csv(eqtl_path, sep="\t")[["SNPName", "ProbeName"]]
            print(outdir, eqtl_df.shape)
            eqtl_df.to_csv(os.path.join(outdir, "ROSMAP-scRNAseq-snpProbe-{}.txt".format(eqtl_type)), header=False, index=False, sep="\t")

            snps.update(eqtl_df["SNPName"])

        pd.Series(list(snps)).to_csv(os.path.join(self.outdir, "ROSMAP-scRNAseq-snps.txt"), header=False, index=False, sep="\t")

    @staticmethod
    def translate_and_filter(values, trans):
        mask = []
        names = []
        for val in values:
            if val in trans:
                mask.append(True)
                names.append(trans[val])
            else:
                mask.append(False)

        return mask, names



if __name__ == '__main__':
    m = main()
    m.start()
