"""
File:         saver.py
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
import os

# Third party imports.

# Local application imports.
from general.utilities import prepare_output_dir


class Saver:
    def __init__(self, df, outdir, signif_cutoff, max_url_len, top):
        self.df = self.set_df(df)
        self.outdir = outdir
        self.signif_cutoff = signif_cutoff
        self.max_url_len = max_url_len
        self.top = top

    @staticmethod
    def set_df(df):
        df["abs_inter"] = df["Inter"].abs()
        df.sort_values(by="abs_inter", ascending=False, inplace=True)
        df.drop(["abs_inter"], axis=1, inplace=True)
        return df

    def save_all(self, exclude=None):
        df = self.df.copy()
        if exclude is not None:
            df = df.loc[~df["Covariate"].isin(exclude), :]
        fpath = os.path.join(self.outdir, "all.txt")
        self.save(df, fpath, self.max_url_len, self.signif_cutoff)

    def save_per_iter(self, exclude=None):
        df = self.df.copy()
        if exclude is not None:
            df = df.loc[~df["Covariate"].isin(exclude), :]

        for iter in df["Iteration"].unique():
            iter_df = self.df.loc[self.df["Iteration"] == iter, :].copy()
            iter_df.drop(["Iteration"], axis=1, inplace=True)

            fpath = os.path.join(self.outdir,
                                 "all_iteration{}.txt".format(iter))
            self.save(iter_df, fpath, self.max_url_len, self.signif_cutoff)

    def save_per_group(self):
        indices_of_interest = []

        for interaction in self.df["Interaction"].unique():
            inter_df = self.df.loc[self.df["Interaction"] == interaction, :].copy()
            inter_df.drop(["Interaction"], axis=1, inplace=True)
            if len(inter_df.index) <= 0:
                return indices_of_interest

            out_dir = os.path.join(self.outdir, '{}_interaction'.format(interaction))
            prepare_output_dir(out_dir)

            for covariate in inter_df["Covariate"].unique():
                cov_df = inter_df.loc[inter_df["Covariate"] == covariate, :].copy()
                cov_df.drop(["Covariate"], axis=1, inplace=True)
                if len(cov_df.index) <= 0:
                    continue

                fpath = os.path.join(out_dir,
                                     "{}_{}.txt".format(interaction, covariate))
                self.save(cov_df, fpath, self.max_url_len, self.signif_cutoff)

                for direction in ["up", "down"]:
                    dir_df = cov_df.loc[cov_df["Direction"] == direction, :].copy()
                    dir_df.drop(["Direction"], axis=1, inplace=True)
                    if len(dir_df.index) <= 0:
                        continue

                    fpath = os.path.join(out_dir, "{}_{}_{}.txt".format(interaction,
                                                                        covariate,
                                                                        direction))
                    self.save(dir_df, fpath, self.max_url_len, self.signif_cutoff)

                    indices_of_interest.extend(dir_df["Index"][:self.top])

        indices_of_interest = list(set(indices_of_interest))
        indices_of_interest.sort()
        return indices_of_interest

    @staticmethod
    def save(df, outfile, max_url_len, z_score_cutoff):
        print("\tWriting output file: {}\tlen: {}".format(os.path.basename(outfile), len(df.index)))
        with open(outfile, 'w') as f:
            f.write('\t'.join([str(x) for x in df.columns]))
            f.write("\n")

            url_string = ""
            url_genes = []
            for i, row in df.iterrows():
                f.write('\t'.join([str(x) for x in row]))
                f.write("\n")

                if (len(url_string) + len(row["ProbeName"])) < max_url_len:
                    if row["HGNCName"] not in url_genes:
                        url_string += row["ProbeName"]
                        url_genes.append(row["HGNCName"])

            f.write("\n")

            f.write("Z-score cutoff: {}\n".format(z_score_cutoff))
            f.write("N genes: {}\n".format(len(df.index)))
            f.write("\n")

            unique_snps = set(df["SNPName"])
            f.write("Number of unique SNPs: {}\n".format(len(unique_snps)))
            f.write("Unique SNPs:\n{}\n".format(', '.join(unique_snps)))
            f.write("\n")

            unique_genes = set(df["HGNCName"])
            f.write("Number of unique genes: {}\n".format(len(unique_genes)))
            f.write("Unique genes:\n{}\n".format(', '.join(unique_genes)))
            f.write("\n")

            f.write("Number of URL genes: {}\n".format(len(url_genes)))
            f.write("URL genes:\n{}\n".format(', '.join(url_genes)))
            f.write("\n")

        f.close()

