"""
File:         main.py
Created:      2020/03/19
Last Changed: 2020/03/20
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
from pathlib import Path
import pickle
import glob
import os
import re

# Third party imports.
from scipy import stats
import numpy as np
import pandas as pd

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists, get_dirname
from general.local_settings import LocalSettings
from general.utilities import get_leaf_dir, get_basename
from general.df_utilities import load_dataframe, save_dataframe


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, groups, force):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param groups: list, the names of groups to analyse.
        :param force: boolean, whether or not to force to redo each step.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Safe arguments.)
        self.eqtl_inpath = settings.get_setting("eqtl_datafile")
        self.cov_inpath = settings.get_setting("cov_datafile")
        self.data_indir = settings.get_setting("data_dir")
        self.g_data_indir = settings.get_setting("groups_data_dir")
        self.g_inter_indir = settings.get_setting("inter_groups_dir")
        self.inter_regex = settings.get_setting("interaction_regex")
        self.group_ids = self.filter_groups(groups)
        self.celltypes = settings.get_setting("celltypes")
        self.force = force

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir,
                                   settings.get_setting("output_dir"))
        prepare_output_dir(self.outdir)

        # Prepare filenames.
        filenames = settings.get_setting("filenames")
        self.obj_filename = filenames["object"]
        self.eqtl_filename = filenames["eqtl"]
        self.geno_filename = filenames["genotype"]
        self.alleles_filename = filenames["alleles"]
        self.expr_filename = filenames["expression"]
        self.cov_filename = filenames["covariates"]
        self.inter_filename = filenames["interaction"]
        self.markers_filename = filenames["markers"]

    def filter_groups(self, groups):
        # Find the groups that are in both group input directories.
        g_in_data_indir = glob.glob(os.path.join(self.g_data_indir, 'group_*'))
        g_data_ids = set([get_leaf_dir(x) for x in g_in_data_indir])

        g_in_inter_indir = glob.glob(
            os.path.join(self.g_inter_indir, 'group_*'))
        g_inter_ids = set([get_leaf_dir(x) for x in g_in_inter_indir])

        ovelapping_ids = g_data_ids.intersection(g_inter_ids)

        # Subset the selected ids.
        group_ids = []
        if "all" in groups:
            group_ids = list(ovelapping_ids)
        else:
            for group_id in ovelapping_ids:
                if group_id in groups:
                    group_ids.append(group_id)

        return group_ids

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting combining groups.")
        self.print_arguments()

        # Combine the indices of each group and combine the interaction
        # matrix if need be.
        inter_outpath = os.path.join(self.outdir, self.inter_filename)
        snp_mask, sample_mask, inter_df = self.combine_groups(inter_outpath)

        print("\nSubsetting data with masks:")
        print("\tSNP mask:\tlength: {}\tlowest index: {}"
              "\thighest index: {}".format(len(snp_mask),
                                           min(snp_mask),
                                           max(snp_mask)))
        print("\tSample mask:\tlength: {}\tlowest index: {}"
              "\thighest index: {}".format(len(sample_mask),
                                           min(sample_mask),
                                           max(sample_mask)))
        print("")

        # Load the eQTL file if either the marker df or the eqtl df needs to be
        # created.
        markers_outpath = os.path.join(self.outdir, self.markers_filename)
        eqtl_outpath = os.path.join(self.outdir, self.eqtl_filename)
        if not check_file_exists(eqtl_outpath) or \
                not check_file_exists(markers_outpath) \
                or self.force:
            print("Loading eQTL file.")
            eqtl_df = load_dataframe(inpath=self.eqtl_inpath,
                                     header=0,
                                     index_col=None)
            eqtl_df = eqtl_df.iloc[snp_mask, :]

            print("Preparing marker matrix.")
            if not check_file_exists(markers_outpath) or self.force:
                self.create_marker_df(inter_df, eqtl_df, markers_outpath)
            else:
                print("\tSkipping step.")

            print("Preparing eQTL matrix.")
            if not check_file_exists(eqtl_outpath) or self.force:
                save_dataframe(outpath=eqtl_outpath, df=eqtl_df,
                               index=False, header=True)
            else:
                print("\tSkipping step.")
            del eqtl_df

        del inter_df

        print("\nPreparing genotype matrix.")
        geno_outpath = os.path.join(self.outdir, self.geno_filename)
        if not check_file_exists(geno_outpath) or self.force:
            geno_df = load_dataframe(inpath=os.path.join(self.data_indir,
                                                         self.geno_filename),
                                     header=0,
                                     index_col=0)
            geno_df = geno_df.iloc[snp_mask, sample_mask]
            save_dataframe(outpath=geno_outpath, df=geno_df,
                           index=True, header=True)
            del geno_df
        else:
            print("\tSkipping step.")

        print("\nPreparing alleles matrix.")
        alleles_outpath = os.path.join(self.outdir, self.alleles_filename)
        if not check_file_exists(alleles_outpath) or self.force:
            alleles_df = load_dataframe(inpath=os.path.join(self.data_indir,
                                                            self.alleles_filename),
                                        header=0,
                                        index_col=0)
            alleles_df = alleles_df.iloc[snp_mask, :]
            save_dataframe(outpath=alleles_outpath, df=alleles_df,
                           index=True, header=True)
            del alleles_df
        else:
            print("\tSkipping step.")

        print("\nPreparing expression matrix.")
        expr_outpath = os.path.join(self.outdir, self.expr_filename)
        if not check_file_exists(expr_outpath) or self.force:
            expr_df = load_dataframe(inpath=os.path.join(self.data_indir,
                                                         self.expr_filename),
                                     header=0,
                                     index_col=0)
            expr_df = expr_df.iloc[snp_mask, sample_mask]
            save_dataframe(outpath=expr_outpath, df=expr_df,
                           index=True, header=True)
            del expr_df
        else:
            print("\tSkipping step.")

        print("\nPreparing covariate matrix.")
        cov_outpath = os.path.join(self.outdir, self.cov_filename)
        if not check_file_exists(cov_outpath) or self.force:
            cov_df = load_dataframe(inpath=self.cov_inpath,
                                    header=0,
                                    index_col=0)
            cov_df = cov_df.iloc[:, sample_mask].copy()
            save_dataframe(outpath=cov_outpath, df=cov_df,
                           index=True, header=True)
            del cov_df
        else:
            print("\tSkipping step.")

    def combine_groups(self, inter_outpath):
        print("Combining groups.")
        snp_mask = np.array([], dtype=np.int16)
        sample_mask = np.array([], dtype=np.int16)
        inter_df = None
        for i, group_id in enumerate(self.group_ids):
            print("  Working on: {:10s} [{}/{} "
                  "{:.2f}%]".format(group_id, i + 1, len(self.group_ids),
                                    (100 / len(self.group_ids)) * (i + 1)))

            # Define the directory names.
            data_indir = os.path.join(self.g_data_indir, group_id)
            inter_indir = os.path.join(self.g_inter_indir, group_id, 'output')

            # Load the group object.
            with open(os.path.join(data_indir, self.obj_filename), "rb") as f:
                group_object = pickle.load(f)

            # Safe the indices.
            snp_mask = np.append(snp_mask, group_object.get_snp_indices())
            sample_mask = np.append(sample_mask, group_object.get_sample_indices())

            if not check_file_exists(inter_outpath) or self.force:
                # Search for the interaction filename.
                inter_inpath = None
                for path in glob.glob(os.path.join(inter_indir, "*")):
                    if re.match(self.inter_regex, get_basename(path)):
                        inter_inpath = path
                        break
                if inter_inpath is None:
                    print("Interaction matrix not found.")
                    exit()

                # Load the interaction file.
                group_inter_df = load_dataframe(inpath=inter_inpath,
                                                header=0,
                                                index_col=0)

                # Merge them.
                if inter_df is None:
                    inter_df = group_inter_df
                else:
                    inter_df = inter_df.merge(group_inter_df,
                                              left_index=True,
                                              right_index=True)

        print("Preparing interaction matrix.")
        if not check_file_exists(inter_outpath) or self.force:
            # Sort the matrix according to the indices.
            inter_df = inter_df.T
            inter_df["index"] = snp_mask
            inter_df.sort_values(by=['index'], inplace=True)
            inter_df.drop(["index"], axis=1, inplace=True)
            inter_df = inter_df.T

            save_dataframe(df=inter_df, outpath=inter_outpath,
                           index=True, header=True)
        else:
            inter_df = load_dataframe(inpath=inter_outpath,
                                      header=0,
                                      index_col=0)

        # Prepare the masks.
        snp_mask = sorted(list(set(snp_mask)))
        sample_mask = sorted(list(set(sample_mask)))

        return snp_mask, sample_mask, inter_df

    def create_marker_df(self, inter_df, eqtl_df, outpath):
        inter_df = inter_df.T
        eqtl_df = eqtl_df[["SNPName", "ProbeName", "HGNCName"]]

        # Calculate the z-score cutoff.
        z_score_cutoff = stats.norm.ppf(
            0.05 / (inter_df.shape[0] * inter_df.shape[1]))
        gini_cutoff = 0.75

        # Subset on the marker genes.
        marker_cols = []
        for colname in inter_df.columns:
            if ("_" in colname) and (colname.split("_")[1] in self.celltypes):
                marker_cols.append(colname)

        marker_df = inter_df.loc[:, marker_cols]
        del inter_df

        # Create a gini dataframe grouped by celltype.
        gini_df = marker_df.copy()
        gini_df = gini_df.abs()
        zscore_mask = list(gini_df.max(axis=1) >= abs(z_score_cutoff))
        gini_df.columns = [x.split("_")[1] for x in gini_df.columns]
        gini_df = gini_df.T.groupby(gini_df.columns).sum().T

        # Calculate the gini impurity.
        gini_values = gini_df.div(gini_df.sum(axis=1), axis=0).pow(2)
        marker_df["gini_impurity"] = 1 - gini_values.sum(axis=1)
        marker_df["eqtl_celltype"] = gini_values.idxmax(axis=1)
        del gini_df

        # Subset the marker df on gini impurity.
        gini_mask = list(marker_df["gini_impurity"] <= gini_cutoff)
        marker_df = marker_df.loc[zscore_mask and gini_mask, :]
        marker_df.index.name = "-"
        marker_df.reset_index(inplace=True)

        # Subset the eQTL dataframe.
        eqtl_df = eqtl_df.loc[zscore_mask and gini_mask, :]
        eqtl_df.reset_index(drop=True, inplace=True)

        # Merge them together.
        merged_df = pd.concat([marker_df, eqtl_df], axis=1)
        merged_df = merged_df.sort_values(by=['eqtl_celltype', 'gini_impurity'])

        # Save the dataframe.
        save_dataframe(df=merged_df,
                       outpath=outpath,
                       header=True,
                       index=False)

        # Save celltype eqtl's HGNC names.
        print("Writing celltype mediated eQTL files.")
        for celltype in marker_df['eqtl_celltype'].unique():
            subset = merged_df.loc[merged_df['eqtl_celltype'] == celltype, :]
            print("\tCelltype: {:20s} {} genes".format(celltype, len(subset.index)))
            if len(subset.index) > 0:
                genes = ', '.join(subset['HGNCName'].to_list())
                outfile = open(os.path.join(get_dirname(outpath),
                                            '{}.txt'.format(celltype)), "w")
                outfile.write(genes)
                outfile.close()

        return eqtl_df

    def print_arguments(self):
        print("Arguments:")
        print("  > EQTL input path: {}".format(self.eqtl_inpath))
        print("  > Covariate input path: {}".format(self.cov_inpath))
        print("  > Data input directory: {}".format(self.data_indir))
        print("  > Group data input directory: {}".format(self.g_data_indir))
        print("  > Group interaction input directory: {}".format(self.g_inter_indir))
        print("  > Ineraction input file regex: {}".format(self.inter_regex))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Group ids: {}".format(self.group_ids))
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > Force: {}".format(self.force))
        print("")