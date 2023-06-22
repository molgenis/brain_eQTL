"""
File:         utilities.py
Created:      2020/10/08
Last Changed: 2020/10/15
Author(s):    M.Vochteloo

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
import pandas as pd
import os

# Third party imports.

# Local application imports.


def check_file_exists(file_path):
    return os.path.exists(file_path) and os.path.isfile(file_path)


def load_dataframe(inpath, header, index_col, sep="\t", low_memory=True,
                   nrows=None, skiprows=None, logger=None):
    df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                     low_memory=low_memory, nrows=nrows, skiprows=skiprows)
    if logger is None:
        print("\tLoaded dataframe: {} "
                    "with shape: {}".format(os.path.basename(inpath),
                                            df.shape))
    else:
        logger.info("\tLoaded dataframe: {} "
                    "with shape: {}".format(os.path.basename(inpath),
                                            df.shape))
    return df


def save_dataframe(df, outpath, header, index, sep="\t", logger=None):
    compression = 'infer'
    if outpath.endswith('.gz'):
        compression = 'gzip'

    df.to_csv(outpath, sep=sep, index=index, header=header,
              compression=compression)
    if logger is None:
        print("\tSaved dataframe: {} "
                    "with shape: {}".format(os.path.basename(outpath),
                                            df.shape))
    else:
        logger.info("\tSaved dataframe: {} "
                    "with shape: {}".format(os.path.basename(outpath),
                                            df.shape))


def construct_dict_from_df(df, key, value):
    return dict(zip(df.loc[:, key], df.loc[:, value]))


def prepare_output_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

