#!/usr/bin/env python3

"""
File:         translate_sample_ids.py
Created:      2021/11/25
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
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Translate Sample IDs"
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
./translate_sample_ids.py -h
"""


class main():
    def __init__(self):
        self.indir =  "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE"
        # self.translate_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_TPM_expression_matrix/MetaBrain_NoENA/data/translate.txt.gz"
        self.translate_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/test_scripts/PsychENCODE/PsychENCODE-linkFile.txt.gz"

        self.outdir = os.path.join(self.indir, 'translated')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        translate_df = self.load_file(self.translate_path, header=0, index_col=None)
        print(translate_df)
        # trans_dict = dict(zip(translate_df.loc[:, "Original ID"], translate_df.loc[:, "Revised ID"]))
        trans_dict = dict(zip(translate_df.loc[:, "psychencode_id"], translate_df.loc[:, "rnaseq_id"]))

        for filename in ["cell_fractions_raw.txt.gz", "cell_fractions_normalized.txt.gz", "DER-02_PEC_Gene_expression_matrix_TPM.txt.gz",  "DER-02_PEC_Gene_expression_matrix_TPM_MGFiltered.txt.gz", "DER-02_PEC_Gene_expression_matrix_TPM_MGFiltered_LOG2.txt.gz"]:
            print(filename)
            df = self.load_file(os.path.join(self.indir, filename), header=0, index_col=0)
            print(df)
            overlap = set(translate_df["psychencode_id"]).intersection(set(df.columns))
            print(len(overlap))
            df = df.loc[:, overlap]
            df.columns = [trans_dict[sample] if sample in trans_dict else sample for sample in df.columns]
            self.save_file(df=df, outpath=os.path.join(self.outdir, filename))

    @staticmethod
    def load_file(inpath, header=None, index_col=None, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        if inpath.endswith(".pkl"):
            df = pd.read_pickle(inpath)
        else:
            df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                             low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))


if __name__ == '__main__':
    m = main()
    m.start()
