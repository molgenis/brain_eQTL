"""
File:         combine_gte_files.py
Created:      2020/10/08
Last Changed: 2021/09/22
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
        self.exclude_files = settings["exclude_files"]
        self.exclude_samples_path = settings["exclude_samples"]
        self.log = log
        self.force = force

        # Prepare an output directory.
        self.outdir = os.path.join(outdir, 'combine_gte_files')
        prepare_output_dir(self.outdir)
        self.outpath = os.path.join(self.outdir, "GTE_combined.txt.gz")

        # Declare variables.
        self.gte_df = None
        self.sample_dict = None
        self.reverse_sample_dict = None
        self.sample_order = None
        self.dataset_to_samples_dict = None

    def start(self):
        self.log.info("Starting combining GTE files.")
        self.print_arguments()

        # Check if GTE output file exist.
        if check_file_exists(self.outpath) and not self.force:
            self.log.info("Skipping step, loading result.")
            self.gte_df = load_dataframe(inpath=self.outpath, header=None, index_col=None, logger=self.log)
        else:
            # Load each GTE file.
            self.log.info("Loading GTE files.")
            self.gte_df = self.combine_files()
            self.save()

        # Construct sample translate dict.
        self.sample_dict = self.create_sample_dict()
        self.sample_order = list(self.gte_df.iloc[:, 1])
        self.dataset_to_samples_dict = self.set_dataset_to_samples_dict()

    def combine_files(self):
        combined = None
        for i, gte_inpath in enumerate(glob.glob(self.inpath)):
            gte_file = os.path.basename(gte_inpath).replace(".txt", "")
            if gte_file in self.exclude_files:
                continue
            df = load_dataframe(inpath=gte_inpath, header=None, index_col=None,
                                logger=self.log)
            df["dataset"] = gte_file
            if combined is None:
                combined = df
            else:
                combined = pd.concat([combined, df], axis=0, ignore_index=True)

        # Remove duplicate entries.
        combined.drop_duplicates(inplace=True)

        # Remove samples.
        if self.exclude_samples_path is not None:
            sample_exclude_df = load_dataframe(inpath=self.exclude_samples_path,
                                               header=None, index_col=None,
                                               logger=self.log)
            pre_shape = combined.shape[0]
            combined = combined.loc[~combined.iloc[:, 1].isin(sample_exclude_df.iloc[:, 0].tolist()), :]
            self.log.warn("\tRemoving '{}' samples".format(pre_shape - combined.shape[0]))

        return combined

    def set_dataset_to_samples_dict(self):
        datasets = list(self.gte_df.iloc[:, 2].unique())

        dataset_to_samples_dict = {}
        for dataset in datasets:
            dataset_to_samples_dict[dataset] = list(self.gte_df.loc[self.gte_df.iloc[:, 2] == dataset, self.gte_df.columns[1]].values)

        return dataset_to_samples_dict

    def create_sample_dict(self):
        sample_dict = {}

        for _, (name1, name2, _) in self.gte_df.iterrows():
            name1 = str(name1)
            name2 = str(name2)

            if name1 not in sample_dict.keys():
                sample_dict[name1] = name2

        return sample_dict

    def save(self):
        save_dataframe(df=self.gte_df, outpath=self.outpath,
                       index=False, header=False, logger=self.log)

        sample_dataset_df = self.gte_df.iloc[:, [1, 2]]
        sample_dataset_df.columns = ["sample", "dataset"]
        save_dataframe(df=sample_dataset_df, outpath=os.path.join(self.outdir, "SampleToDataset.txt.gz"),
                       index=False, header=True, logger=self.log)

    def clear_variables(self):
        self.inpath = None
        self.force = None

    def get_gte_file(self):
        return self.outpath

    def get_gte_df(self):
        return self.gte_df

    def get_sample_order(self):
        return self.sample_order

    def get_sample_dict(self):
        return self.sample_dict

    def get_dataset_to_samples_dict(self):
        return self.dataset_to_samples_dict

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Input files: {}".format(self.inpath))
        self.log.info("  > Exclude files: {}".format(", ".join(self.exclude_files)))
        self.log.info("  > Exclude samples: {}".format(self.exclude_samples_path))
        self.log.info("  > Output path: {}".format(self.outpath))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
