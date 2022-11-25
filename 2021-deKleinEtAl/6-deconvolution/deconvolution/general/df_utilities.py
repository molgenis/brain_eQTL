"""
File:         df_utilities.py
Created:      2020/03/19
Last Changed: 2020/04/23
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

# Third party imports.
import pandas as pd

# Local application imports.
from .utilities import get_basename


def load_dataframe(inpath, header, index_col, sep="\t", low_memory=True,
                   nrows=None, skiprows=None):
    """
    Method for reading a comma-separated values (csv) file into a pandas
    DataFrame.

    :param inpath: str, the file to be read.
    :param header: int, row number(s) to use as the column names, and the
                   start of the data.
    :param index_col: int, column(s) to use as the row labels of the DataFrame,
                      either given as string name or column index.
    :param sep: str, delimiter to use.
    :param low_memory: boolean, internally process the file in chunks,
                       resulting in lower memory use while parsing, but
                       possibly mixed type inference.
    :param nrows: int, number of rows of file to read.
    :param skiprows: list, the index of rows to skip.
    :return df: DataFrame, the pandas dataframe.
    """
    df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                     low_memory=low_memory, nrows=nrows, skiprows=skiprows)
    print("\tLoaded dataframe: {} with shape: {}".format(get_basename(inpath),
                                                         df.shape))
    return df


def save_dataframe(df, outpath, header, index, sep="\t"):
    """
    Method for writing an dataframe to a comma-separated values (csv) file.

    :param df: DataFrame, the pandas dataframe.
    :param outpath: str, the filepath for the dataframe.
    :param header: boolean, write out the column names.
    :param index: boolean, write row names (index).
    :param sep: str, field delimiter for the output file.
    """
    compression = 'infer'
    if outpath.endswith('.gz'):
        compression = 'gzip'

    df.to_csv(outpath, sep=sep, index=index, header=header,
              compression=compression)
    print("\tSaved dataframe: {} with shape: {}".format(get_basename(outpath),
                                                        df.shape))
