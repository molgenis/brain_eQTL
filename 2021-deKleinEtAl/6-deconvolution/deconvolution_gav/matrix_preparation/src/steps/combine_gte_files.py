"""
File:         combine_gte_files.py
Created:      2020/10/08
Last Changed: 2020/10/12
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
import glob
import os

# Third party imports.
import pandas as pd

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe


class CombineGTEFiles:
    def __init__(self, settings, log, force, outdir):
        self.inpath = os.path.join(settings["input_directory"],
                                   settings["filename_regex"])
        self.log = log
        self.force = force

        # Prepare an output directory.
        self.outdir = os.path.join(outdir, 'combine_gte_files')
        prepare_output_dir(self.outdir)
        self.outpath = os.path.join(self.outdir, "GTE_combined.txt.gz")

        # Declare variables.
        self.gte = None
        self.sample_dict = None
        self.reverse_sample_dict = None
        self.sample_order = None

    def start(self):
        self.log.info("Starting combining GTE files.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.outpath) and not self.force:
            self.log.info("Skipping step, loading result.")
            self.gte = load_dataframe(inpath=self.outpath, header=None,
                                      index_col=None, logger=self.log)
        else:
            # Load each GTE file.
            self.log.info("Loading GTE files.")
            self.gte = self.combine_files()
            self.save()

        # Construct sample translate dict.
        self.sample_dict = self.create_sample_dict()
        self.reverse_sample_dict = self.create_reverse_sample_dict()
        self.sample_order = self.set_sample_order()

    def combine_files(self):
        combined = None
        for i, infile in enumerate(glob.glob(self.inpath)):
            df = load_dataframe(inpath=infile, header=None, index_col=None,
                                logger=self.log)
            if combined is None:
                combined = df
            else:
                combined = pd.concat([combined, df], axis=0, ignore_index=True)

        # Remove duplicate entries.
        combined.drop_duplicates(inplace=True)

        return combined

    def save(self):
        save_dataframe(df=self.gte, outpath=self.outpath,
                       index=False, header=False, logger=self.log)

    def clear_variables(self):
        self.inpath = None
        self.force = None

    def create_sample_dict(self):
        sample_dict = {}

        for _, (name1, name2) in self.gte.iterrows():
            name1 = str(name1)
            name2 = str(name2)

            if name1 not in sample_dict.keys():
                sample_dict[name1] = name2

        return sample_dict

    def create_reverse_sample_dict(self):
        if self.sample_dict is None:
            return None
        return {v: k for k, v in self.sample_dict.items()}

    def get_outpath(self):
        return self.outpath

    def get_gte(self):
        return self.gte

    def get_sample_dict(self, reverse=False):
        if reverse:
            return self.reverse_sample_dict
        else:
            return self.sample_dict

    def set_sample_order(self):
        return list(self.gte.iloc[:, 1])

    def get_sample_order(self):
        return self.sample_order

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Input files: {}".format(self.inpath))
        self.log.info("  > Output path: {}".format(self.outpath))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
