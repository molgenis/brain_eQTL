"""
File:         combine_eqtlprobes.py
Created:      2020/10/08
Last Changed: 2021/07/08
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
import os

# Third party imports.
import pandas as pd

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe


class CombineEQTLProbes:
    def __init__(self, settings, log, force, outdir):
        self.indir = settings["input_directory"]
        self.iter_dirname = settings["iteration_dirname"]
        self.in_filename = settings["in_filename"]
        self.n_iterations = settings["iterations"]
        self.log = log
        self.force = force

        # Prepare an output directory.
        self.outdir = os.path.join(outdir, 'combine_eqtlprobes')
        prepare_output_dir(self.outdir)
        self.outpath = os.path.join(self.outdir, "eQTLprobes_combined.txt.gz")

        # Declare variables.
        self.eqtl_df = None

    def start(self):
        self.log.info("Starting combining eQTL probe files.")
        self.print_arguments()

        # Check if output file exist.
        if not check_file_exists(self.outpath) or self.force:
            # Load each GTE file.
            self.log.info("Loading eQTLprobes files.")
            self.eqtl_df = self.combine_files()
            self.save()
        else:
            self.log.info("Skipping step.")

    def combine_files(self):
        combined = None
        for i in range(1, self.n_iterations+1):
            infile = os.path.join(self.indir, self.iter_dirname + str(i),
                                  self.in_filename)
            df = load_dataframe(inpath=infile, header=0, index_col=False,
                                logger=self.log)
            df["Iteration"] = i
            if combined is None:
                combined = df
            else:
                combined = pd.concat([combined, df], axis=0, ignore_index=True)

        # Remove duplicate entries.
        combined.drop_duplicates(inplace=True)

        return combined

    def save(self):
        save_dataframe(df=self.eqtl_df, outpath=self.outpath,
                       index=False, header=True, logger=self.log)
        save_dataframe(df=self.eqtl_df.loc[:, ["ProbeName", "SNPName"]], outpath=os.path.join(self.outdir, "snp_gene_list.txt"),
                       index=False, header=True, logger=self.log)

    def clear_variables(self):
        self.indir = None
        self.iter_dirname = None
        self.n_iterations = None
        self.in_filename = None
        self.force = None

    def get_outpath(self):
        return self.outpath

    def get_eqtl_file(self):
        return self.outpath

    def get_eqtl_df(self):
        return self.eqtl_df

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Input directory: {}".format(self.indir))
        self.log.info("  > Iteration directory: {}".format(self.iter_dirname))
        self.log.info("  > N. Iterations: {}".format(self.n_iterations))
        self.log.info("  > Input filename: {}".format(self.in_filename))
        self.log.info("  > Output path: {}".format(self.outpath))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
