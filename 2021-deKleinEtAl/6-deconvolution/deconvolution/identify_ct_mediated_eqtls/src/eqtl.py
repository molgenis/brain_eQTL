"""
File:         eqtl.py
Created:      2020/06/09
Last Changed: 2020/06/11
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

# Third party imports.
import scipy.stats as st
import numpy as np


# Local application imports.


class Eqtl:
    def __init__(self, index, snp_name, probe_name, hgnc_name, iteration,
                 eqtl_zscore, alleles, gwas_ids, traits, signif_cutoff,
                 maf_cutoff, selections, genotype, expression, covariates,
                 inter_zscores, inter_tvalues):
        # General information.
        self.index = index
        self.snp_name = snp_name
        self.probe_name = probe_name
        self.hgnc_name = hgnc_name
        self.iteration = iteration
        self.eqtl_zscore = round(eqtl_zscore, 2)
        self.alleles = alleles
        self.gwas_ids = gwas_ids
        self.traits = traits
        self.signif_cutoff = signif_cutoff
        self.maf_cutoff = maf_cutoff
        self.selections = {key.split("_")[-1].lower(): value for key, value in
                           selections.items()}

        # The interaction data.
        df = self.combine_data(genotype, expression)
        self.covariates = covariates

        # Significance of the interaction.
        inter_zscores = {index.split("_")[-1].lower(): round(value, 2) for
                         index, value in inter_zscores.iteritems()}
        inter_tvalues = {index.split("_")[-1].lower(): round(value, 2) for
                         index, value in inter_tvalues.iteritems()}

        # Variables.
        self.minor_allele = self.alleles[-1]
        self.major_allele = self.alleles[0]
        self.genotype_flip = False

        # Calculations.
        zero_geno_count, two_geno_count, self.maf = self.calculate_maf(genotype)
        self.genotype_flip, self.allele_info = self.set_allele_info(
            zero_geno_count, two_geno_count, alleles)
        if self.genotype_flip:
            df = self.flip_genotype(df)
        self.df = df
        interaction_info = self.assess_interactions(inter_zscores)

        # Combine the celltype information.
        self.cov_info = {}
        for key, value in interaction_info.items():
            value["zscore"] = inter_zscores[key]
            value["tvalue"] = inter_tvalues[key]
            self.cov_info[key] = value

    @staticmethod
    def calculate_maf(s):
        counts = s.round(0).value_counts()
        for x in [0.0, 1.0, 2.0]:
            if x not in counts:
                counts.loc[x] = 0

        zero_geno_count = (counts[0.0] * 2) + counts[1.0]
        two_geno_count = (counts[2.0] * 2) + counts[1.0]
        maf = min(zero_geno_count, two_geno_count) / (
                    zero_geno_count + two_geno_count)
        return zero_geno_count, two_geno_count, maf

    @staticmethod
    def set_allele_info(zero_geno_count, two_geno_count, alleles):
        info = None
        flip = None
        if two_geno_count > zero_geno_count:
            flip = True
            info = {"MINOR": (alleles[0], zero_geno_count),
                    "MAJOR": (alleles[-1], two_geno_count)}
        elif zero_geno_count > two_geno_count:
            flip = False
            info = {"MINOR": (alleles[-1], two_geno_count),
                    "MAJOR": (alleles[0], zero_geno_count)}

        return flip, info

    @staticmethod
    def combine_data(genotype, expression):
        genotype_df = genotype.to_frame()
        expression_df = expression.to_frame()
        data = genotype_df.merge(expression_df,
                                 left_index=True,
                                 right_index=True)
        data.columns = ["genotype", "expression"]
        data["geno_group"] = data["genotype"].round(0)
        del genotype_df, expression_df

        data.replace(-1, np.nan, inplace=True)
        data.dropna(inplace=True)
        return data

    @staticmethod
    def flip_genotype(df):
        df["geno_group"] = 2.0 - df["geno_group"]
        df["genotype"] = 2.0 - df["genotype"]
        return df

    def assess_interactions(self, inter_zscores):
        interaction_info = {}

        cov_df = self.covariates.copy()
        for index, row in cov_df.iterrows():
            celltype = index.split("_")[-1].lower()

            info = {"include": False,
                    "spread": None,
                    "low_end_std": None,
                    "high_end_std": None,
                    "interaction": None,
                    "major_yhat": None,
                    "minor_yhat": None,
                    "direction": None}

            if inter_zscores[celltype] > self.signif_cutoff:
                info["include"] = True

                work_df = self.df.copy()
                work_df[celltype] = row.copy()
                min_fraction = work_df[celltype].min()
                max_fraction = work_df[celltype].max()

                info["spread"] = abs(max_fraction - min_fraction)

                low_end = {}
                high_end = {}
                groups_present = []
                for group in [0.0, 1.0, 2.0]:
                    subset = work_df.loc[work_df["geno_group"] == group, :].copy()

                    if len(subset.index) > 0:
                        slope, intercept, _, _, _ = st.linregress(
                            subset[celltype], subset["expression"])
                        low_end[group] = intercept + subset[celltype].min() * slope
                        high_end[group] = intercept + subset[celltype].max() * slope
                        groups_present.append(group)
                low_end_std = np.std(list(low_end.values()))
                high_end_std = np.std(list(high_end.values()))

                info["low_end_std"] = low_end_std
                info["high_end_std"] = high_end_std

                if low_end_std < high_end_std:
                    info["interaction"] = "positive"
                if low_end_std > high_end_std:
                    info["interaction"] = "negative"

                max_geno = max(groups_present)
                min_geno = min(groups_present)
                info["major_yhat"] = high_end[min_geno]
                info["minor_yhat"] = high_end[max_geno]

                if high_end[max_geno] > high_end[min_geno]:
                    info["direction"] = "up"
                elif high_end[max_geno] < high_end[min_geno]:
                    info["direction"] = "down"

            interaction_info[celltype] = info

        return interaction_info

    def get_data(self):
        data = []
        for key, value in self.cov_info.items():
            selection = self.selections[key]
            if self.maf > self.maf_cutoff and value[
                "interaction"] in selection and value["direction"] is not None:
                data.append(
                    [self.index, self.snp_name, self.probe_name, self.hgnc_name,
                     self.iteration, self.df.shape[0], self.maf,
                     self.eqtl_zscore, value["tvalue"], key,
                     value["interaction"], value["direction"], self.gwas_ids,
                     self.traits])
        return data

    def print_info(self):
        print("eQTL info:")
        print("  > Index: {}".format(self.index))
        print("  > SNP name: {}".format(self.snp_name))
        print("  > Probe name: {}".format(self.probe_name))
        print("  > HGNC name: {}".format(self.hgnc_name))
        print("  > Iteration: {}".format(self.iteration))
        print("  > eQTL z-score: {}".format(self.eqtl_zscore))
        print("  > Alleles: {}".format(self.alleles))
        print("  > GWAS IDs: {}".format(self.gwas_ids))
        print("  > Traits: {}".format(self.traits))
        print("  > DataFrame shape: {}".format(self.df.shape))
        print("  > Genotype flip: {}".format(self.genotype_flip))
        print("  > Allele info: {}".format(self.allele_info))
        print("  > Significance cut-off: {}".format(self.signif_cutoff))
        print("  > Cell type info:")
        for key, value in self.cov_info.items():
            if value["include"]:
                print("\t{:15s} = spread: {:.2f}\tlow-end STD: {:.2f}\t"
                      "high-end STD: {:.2f}\tinteraction: {}\t"
                      "major yhat: {:.2f}\tminor yhat: {:.2f}\t"
                      "direction: {}\tz-score: {:.2f}\t"
                      "t-value: {:.2f}\t".format(key,
                                                 value["spread"],
                                                 value["low_end_std"],
                                                 value["high_end_std"],
                                                 value["interaction"],
                                                 value["major_yhat"],
                                                 value["minor_yhat"],
                                                 value["direction"],
                                                 value["zscore"],
                                                 value["tvalue"]))
        print("")
