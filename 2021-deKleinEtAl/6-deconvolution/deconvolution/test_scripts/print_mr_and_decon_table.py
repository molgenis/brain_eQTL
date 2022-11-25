#!/usr/bin/env python3

"""
File:         print_mr_and_decon_table.py
Created:      2021/02/16
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

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Print MR and Decon Table"
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
        self.table_path = "../data/SupplementaryTable8-MR_findings_passing_suggestive_threshold_gaps_filled.xlsx"
        self.cell_types = ["Astrocyte_FDR", "EndothelialCell_FDR", "Macrophage_FDR", "Neuron_FDR", "Oligodendrocyte_FDR"]
        self.traits = {
            'ADHD': ['ADHD'],
            'AD': ['Alzheimer’s disease'],
            'ALS': ['Amyotrophic Lateral Sclerosis'],
            'autism': ['Autism Spectrum Disorder'],
            'bipolar disorder': ['bipolar disorder'],
            'depression': ['Depression (broad)'],
            'epilepsy': ['epilepsy, all documented cases',
                         'focal epilepsy, all documented cases',
                         'focal epilepsy, documented hippocampal sclerosis',
                         'focal epilepsy, documented lesion negative',
                         'focal epilepsy, documented lesion other than hippocampal sclerosis',
                         'generalized epilepsy with tonic-clonic seizures',
                         'generalized epilepsy, all documented cases',
                         'juvenile absence epilepsy',
                         'juvenile myoclonic epilepsy '],
            'frontotemporal dementia': ['Frontotemporal Dementia'],
            'years of schooling and fluid intelligence': ['Intelligence',
                                                          'Years of schooling'],
            'MS': ['multiple sclerosis'],
            'PD': ['Parkinson’s disease'],
            'schizophrenia': ['schizophrenia'],
            'brain volume': ['Caudate volume',
                             'Hippocampus volume',
                             'Intracranial volume',
                             'Pallidum volume',
                             'Putamen volume',
                             'Thalamus volume']}
        self.outcomes = []
        for outcome in self.traits.values():
            self.outcomes.extend(outcome)
        if len(self.outcomes) != len(set(self.outcomes)):
            print("DUPLICATES IN TRAITS DICT")
            exit()

    def start(self):
        self.print_arguments()

        print("")
        print("### Step2 ###")
        print("Loading MR results.")
        sheets_dict = self.load_file(self.table_path, header=0, index_col=0)

        for sheet_name, df in sheets_dict.items():
            for outcome in df["outcome"].unique():
                if outcome not in self.outcomes:
                    print("Missing outcome '{}' in traits dict".format(outcome))
                    exit()

            combined_result_df = None
            stats_data = []
            stats_index = []
            for trait, outcome_list in self.traits.items():
                trait_df = df.loc[df["outcome"].isin(outcome_list), :].copy()
                trait_df.index = trait_df["SNP"] + "_" + trait_df["gene"]
                n_total = trait_df.shape[0]
                if len(trait_df.index) != n_total:
                    print("Unique SNP - gene assumptions doesnt hold.")
                    exit()

                trait_df = trait_df.loc[~trait_df[self.cell_types].isnull().any(axis=1)]
                n_tested = trait_df.shape[0]

                trait_df = trait_df.loc[trait_df[self.cell_types].min(axis=1) < 0.05, :]
                n_hits = trait_df.shape[0]

                decon_df = trait_df[self.cell_types].copy()
                decon_df.columns = [x.split("_")[0] for x in decon_df.columns]

                mr_df = trait_df[["pass bonferroni correction", "coloc.relaxed", "coloc.strict"]].copy()
                mr_df["MR"] = mr_df["pass bonferroni correction"] == 1
                mr_df["coloc"] = (mr_df["coloc.relaxed"] + mr_df["coloc.relaxed"]) >= 1
                mr_df.drop(["pass bonferroni correction", "coloc.relaxed", "coloc.strict"], axis=1, inplace=True)

                del trait_df

                result_df = decon_df[decon_df < 0.05].count().to_frame()
                result_df.columns = [trait]

                stats_data.append([n_total, n_tested, n_hits])
                stats_index.append(trait)

                if n_hits > 0:
                    if combined_result_df is None:
                        combined_result_df = result_df
                    else:
                        combined_result_df = pd.merge(combined_result_df, result_df, left_index=True, right_index=True)

                out_string = self.construct_out_string(trait, n_tested, n_total, n_hits, decon_df, mr_df)
                print("Trait: {}".format(trait))
                print(out_string)
                print("")

            combined_result_df["total"] = combined_result_df.sum(axis=1)
            print(combined_result_df.T)
            stats_df = pd.DataFrame(stats_data, columns=["n total", "n tested", "n hits"], index=stats_index)
            stats_df.loc["total", :] = stats_df.sum(axis=0)
            print(stats_df)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  skiprows=0, sheet_name=None):
        if path.endswith(".xlsx"):
            df = pd.read_excel(path, header=header, index_col=index_col,
                         nrows=nrows, skiprows=skiprows, sheet_name=sheet_name)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, skiprows=skiprows)
        return df

    def construct_out_string(self, trait, n_tested, n_total, n_hits, decon_df, mr_df):
        if n_hits == 0:
            return "We also evaluated {} of the {} {} associated eQTL instruments for an interaction effect with cell type proportions, none of which were significant.".format(n_tested, n_total, trait)
        else:

            part1 = "We also evaluated {} of the {} {} associated eQTL instruments for an interaction effect with cell type proportions.".format(n_tested, n_total, trait)

            #####################################################################

            suffix = ""
            if n_hits > 1:
                suffix = "s"

            decon_hits = {}
            decon_order = [x.lower() + "s" for x in list(decon_df.columns)]
            for index, row in decon_df.iterrows():
                gene = index.split("_")[1]
                signif_cell_types = " and ".join([x.lower() + "s" for x in list(row[row < 0.05].index.values)])
                if signif_cell_types in decon_hits:
                    genes = decon_hits[signif_cell_types]
                    genes.append(gene)
                    decon_hits[signif_cell_types] = genes
                else:
                    decon_hits[signif_cell_types] = [gene]
                if "and" in signif_cell_types and signif_cell_types not in decon_order:
                    decon_order.append(signif_cell_types)

            part2_parts = []
            for cell_type in decon_order:
                if cell_type in decon_hits:
                    genes = decon_hits[cell_type]
                    if len(genes) > 0:
                        if cell_type == "endothelialcells":
                            cell_type = "endothelial cells"
                        part2_parts.append("{} with {} ({})".format(len(genes), cell_type, self.special_join(genes)))

            part2 = "In total, {} eQTL instrument{} where significant ieQTLs, ".format(n_hits, suffix) + self.special_join(part2_parts) + "."

            #####################################################################

            part3 = "None of these ieQTLs, however, overlapped with significant MR and colocalization signals."
            if n_hits == 1:
                part3 = "This ieQTL did, however, not overlap with significant MR and colocalization signals."

            mr_order = ["both MR and colocalization", "MR", "colocalization"]
            mr_hits = {x: [] for x in mr_order}
            mr_has_hit = False
            for index, row in mr_df.iterrows():
                gene = index.split("_")[1]
                key = None
                if row["MR"] and row["coloc"]:
                    key = "both MR and colocalization"
                elif row["MR"] and not row["coloc"]:
                    key = "MR"
                elif not row["MR"] and row["coloc"]:
                    key = "colocalization"

                if key is not None:
                    mr_has_hit = True
                    genes = mr_hits[key]
                    genes.append(gene)
                    mr_hits[key] = genes

            if mr_has_hit:
                part3_parts = []
                for key in mr_order:
                    genes = mr_hits[key]
                    if len(genes) > 0:
                        text = ""
                        prefix = ""
                        if len(genes) > 1:
                            prefix = "s"
                        if key != "both MR and colocalization":
                            text = "only "
                        part3_parts.append("{} eQTL instrument{} passed the significance threshold {}for {} ({})".format(len(genes), prefix, text, key, self.special_join(genes)))
                part3 = "Of these, " + self.special_join(part3_parts) + "."

            #####################################################################

            return part1 + " " + part2 + " " + part3

    def print_arguments(self):
        print("Arguments:")
        print("  > Table path: {}".format(self.table_path))
        print("")

    @staticmethod
    def special_join(iterable):
        if len(iterable) == 0:
            return ""
        elif len(iterable) == 1:
            return iterable[0]
        else:
            return ", ".join(iterable[:-1]) + " and " + iterable[-1]


if __name__ == '__main__':
    m = main()
    m.start()
