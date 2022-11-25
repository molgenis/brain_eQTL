#!/usr/bin/env python3

"""
File:         cell_fraction_gene_correlations.py
Created:      2020/09/08
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
from pathlib import Path
import argparse
import gzip
import math
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Cell Fractions Gene Correlations"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
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
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.ge_path = getattr(arguments, 'gene_expression')
        self.nrows = getattr(arguments, 'nrows')
        self.gene_info_path = getattr(arguments, 'gene_info')
        self.gene_filter_path = getattr(arguments, 'gene_filter')
        self.gf_id = getattr(arguments, 'gene_filter_id')
        self.min_corr = getattr(arguments, 'min_corr')
        self.compare_path = getattr(arguments, 'compare')
        self.extension = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "gene_cellcount%_corr")
        self.colormap = {
            "Neuron": "#b38d84",
            "Oligodendrocyte": "#5d9166",
            "EndothelialCell": "#f2a7a7",
            "Microglia": "#e8c06f",
            "Macrophage": "#e8c06f",
            "Astrocyte": "#9b7bb8",
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.infile_basename = self.get_basename(self.ge_path)
        self.coef_outpath = os.path.join(self.outdir, '{}_coefficients.txt.gz'.format(self.infile_basename))
        self.pvalue_outpath = os.path.join(self.outdir, '{}_pvalues.txt.gz'.format(self.infile_basename))

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fracttions matrix")
        parser.add_argument("-ge",
                            "--gene_expression",
                            type=str,
                            required=True,
                            help="The path to the gene expression matrix")
        parser.add_argument("-n",
                            "--nrows",
                            type=int,
                            default=None,
                            required=False,
                            help="The number of genes to analyze. "
                                 "Default: None.")
        parser.add_argument("-gi",
                            "--gene_info",
                            type=str,
                            required=True,
                            help="The path to the gene annotation matrix.")
        parser.add_argument("-f",
                            "--gene_filter",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the gene filter matrix. "
                                 "Default: None")
        parser.add_argument("-gf_id",
                            "--gene_filter_id",
                            type=str,
                            required=False,
                            default=None,
                            help="The column on which to filter.")
        parser.add_argument("-min",
                            "--min_corr",
                            type=float,
                            default=0.5,
                            required=False,
                            help="The minimal correlation a gene must have"
                                 "with cell count % before including"
                                 "it in the summary. Default: 0.5.")
        parser.add_argument("-c",
                            "--compare",
                            type=str,
                            required=False,
                            help="THe path to a comparison file. (tmp)")
        parser.add_argument("-e",
                            "--extension",
                            type=str,
                            choices=["png", "pdf"],
                            default="png",
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    @staticmethod
    def get_basename(path):
        basename = os.path.basename(path)
        file_extensions = ["txt", "gz", "zip"]
        while True:
            found = False
            for extension in file_extensions:
                if basename.endswith('.{}'.format(extension)):
                    basename = basename.replace('.{}'.format(extension), "")
                    found = True
            if not found:
                break

        return basename

    def start(self):
        if os.path.exists(self.coef_outpath) and os.path.exists(self.pvalue_outpath):
            coef_df, pvalue_df = self.load_existing()
        else:
            coef_df, pvalue_df = self.correlate()

        gene_df = self.find_top_correlating_genes(coef_df.copy(), pvalue_df.copy())

        trans_dict = dict(zip(gene_df.loc[:, "ArrayAddress"], gene_df.loc[:, "Symbol"]))
        coef_df.index = [trans_dict[x] for x in coef_df.index]
        if self.compare_path == "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/41586_2019_1195_MOESM8_ESM.csv":
            self.compare_with_article(coef_df)
        elif self.compare_path == "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt":
            self.compare_with_cellmap_profile(coef_df)

    def load_existing(self):
        print("Loading coefficients matrix.")
        coef_df = pd.read_csv(self.coef_outpath, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.coef_outpath),
                                      coef_df.shape))

        print("Loading p-values matrix.")
        pvalue_df = pd.read_csv(self.pvalue_outpath, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.pvalue_outpath),
                                      pvalue_df.shape))

        return coef_df, pvalue_df

    def correlate(self):
        print("Loading cell fractions matrix.")
        cf_df = pd.read_csv(self.cf_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.cf_path),
                                      cf_df.shape))
        celltypes = cf_df.columns

        print("Performing correlations.")
        coefs_data = []
        p_vals_data = []
        indices = []
        columns = None
        with gzip.open(self.ge_path, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % 250 == 0):
                    print("\t{}/{} genes done".format(i, self.nrows))
                if self.nrows is not None and i > self.nrows:
                    break

                splitted_line = np.array(line.decode().strip('\n').split('\t'))
                index = splitted_line[0]
                data = splitted_line[1:]
                if i == 0:
                    columns = data
                    continue

                df = pd.Series(data, index=columns, name=index).astype(np.float64).to_frame()
                corr_df = cf_df.merge(df, left_index=True, right_index=True)

                coefs = []
                p_vals = []
                for ct in celltypes:
                    coef, p = stats.spearmanr(corr_df[index], corr_df[ct])
                    coefs.append(coef)
                    p_vals.append(p)

                coefs_data.append(coefs)
                p_vals_data.append(p_vals)
                indices.append(index)

        f.close()

        # Construct data frame.
        coef_df = pd.DataFrame(coefs_data, index=indices, columns=celltypes)
        pvalue_df = pd.DataFrame(p_vals_data, index=indices, columns=celltypes)

        # Save.
        coef_df.to_csv(self.coef_outpath, sep="\t", index=True, header=True,
                       compression="gzip")
        pvalue_df.to_csv(self.pvalue_outpath, sep="\t", index=True, header=True,
                         compression="gzip")

        return coef_df, pvalue_df

    def find_top_correlating_genes(self, coef_df, pvalue_df):
        print("Combining data.")
        coef_df.reset_index(inplace=True)
        coef_df_m = coef_df.melt(id_vars="index")
        coef_df_m.columns = ["gene", "celltype", "coefficient"]

        pvalue_df.reset_index(inplace=True)
        pvalue_df_m = pvalue_df.melt(id_vars="index")
        pvalue_df_m.columns = ["gene", "celltype", "pvalue"]

        corr_df = coef_df_m.merge(pvalue_df_m, on=["gene", "celltype"])
        corr_df["abs_coefficient"] = corr_df["coefficient"].abs()

        print("\tLoading gene information data.")
        gene_df = pd.read_csv(self.gene_info_path, sep="\t", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.gene_info_path),
                                      gene_df.shape))

        df = pd.merge(corr_df, gene_df, left_on="gene", right_on="ArrayAddress")

        print("Saving complete data frame.")
        df.to_csv(os.path.join(self.outdir, '{}_AllCT_correlations.txt.gz'.format(self.infile_basename)),
                  sep="\t", index=True, header=True, compression="gzip")

        if self.gene_filter_path is not None and self.gf_id is not None:
            print("\tLoading gene filter matrix.")
            filter_df = pd.read_csv(self.gene_filter_path, sep="\t", header=0)
            print("\tLoaded dataframe: {} "
                  "with shape: {}".format(os.path.basename(self.gene_filter_path),
                                          filter_df.shape))

            print("Pre-filter shape: {}".format(df.shape))
            df = df.loc[df["gene"].isin(filter_df[self.gf_id]), :]
            print("Post-filter shape: {}".format(df.shape))

            print("Saving filtered data frame.")
            df.to_csv(os.path.join(self.outdir,
                                   '{}_FilteredCT_correlations.txt.gz'.format(
                                       self.infile_basename)),
                      sep="\t", index=True, header=True, compression="gzip")

        filtered_df = df.loc[df["abs_coefficient"] >= self.min_corr, :].copy()

        print("Visualising chromosoom - cell count % correlation")
        self.chr_pos_barplot(filtered_df)

        print("Splitting best genes per cell type.")
        for ct in filtered_df["celltype"].unique():
            subset = filtered_df.loc[filtered_df["celltype"] == ct, :].copy()
            subset.sort_values(by="abs_coefficient", ascending=False, inplace=True)
            print(subset)
            subset.to_csv(os.path.join(self.outdir, '{}_{}.txt.gz'.format(self.infile_basename, ct)),
                          sep="\t", index=True, header=True, compression="gzip")

        return gene_df

    def chr_pos_barplot(self, df):
        print(df)
        order = [str(x) for x in range(22)] + ["x", "y"]

        sns.set_style("ticks")
        g = sns.catplot(x="abs_coefficient", y="Chr",
                        hue="celltype", col="celltype",
                        palette=self.colormap,
                        data=df, kind="point", join=False,
                        dodge=True, order=order,
                        height=4, aspect=.7)

        g.savefig(os.path.join(self.outdir, "chr_barplot.{}".format(self.extension)))
        plt.close()

    def compare_with_article(self, coef_df):
        print("Loading comparison file.")
        art_df = pd.read_csv(self.compare_path, sep=",", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.compare_path),
                                      art_df.shape))

        merg_df = art_df.merge(coef_df, left_on="gene.name", right_index=True)
        merg_df.sort_index(inplace=True)
        merg_df = merg_df.set_index("gene.name")

        merg_df.to_csv(
            os.path.join(self.outdir, self.get_basename(self.compare_path) + ".txt.gz"),
            sep="\t", index=True, header=True, compression="gzip")

        ct_dict = {"Ex": "Neuron", "In": "Neuron", "Ast": "Astrocyte",
                   "Oli": "Oligodendrocyte", "Opc": "Oligodendrocyte",
                   "Mic": "Macrophage"}
        overlap_dict = {}
        for index, row in merg_df.iterrows():
            celltype = None
            for key, value in ct_dict.items():
                if row["subpopulation"].startswith(key):
                    celltype = value
                    break

            if celltype is None:
                continue

            if celltype in overlap_dict.keys():
                (genes, overlap) = overlap_dict[celltype]
                genes.append(index)

                if abs(row[celltype]) > self.min_corr:
                    overlap.add(index)
                overlap_dict[celltype] = (genes, overlap)
            else:
                overlap_dict[celltype] = ([index], set(index))

        print("CellType\tAll\tUnique\t>{}coeff".format(self.min_corr))
        for key, (genes, overlap) in overlap_dict.items():
            print(key, len(genes), len(set(genes)), len(overlap))

        self.plot_coefficients_per_group(merg_df, "subpopulation", coef_df.columns)

    def compare_with_cellmap_profile(self, coef_df):
        print("Loading CellMap reference profile.")
        sign_df = pd.read_csv(self.compare_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.compare_path),
                                      sign_df.shape))

        sign_df = sign_df.subtract(sign_df.mean(axis=1), axis=0).divide(sign_df.std(axis=1), axis=0).idxmax(axis=1).to_frame()
        sign_df.columns = ["SignatureGene"]
        sign_df.index.name = "GeneSymbol"

        merg_df = sign_df.merge(coef_df, left_index=True, right_index=True)

        merg_df.to_csv(os.path.join(self.outdir, self.get_basename(self.compare_path) + ".txt.gz"),
                       sep="\t", index=True, header=True, compression="gzip")

        self.plot_coefficients_per_group(merg_df, "SignatureGene", coef_df.columns)

    def plot_coefficients_per_group(self, df, col_id, value_vars):

        df_m = df.melt(id_vars=col_id, value_vars=value_vars)
        df_m["abs_coef"] = df_m["value"].abs()

        n_groups = len(df[col_id].unique())

        sns.set(style="ticks")
        order = list(df_m["variable"].unique())
        order.sort()
        g = sns.catplot(x="variable", y="abs_coef", col=col_id,
                        col_wrap=math.ceil(np.sqrt(n_groups)), data=df_m, kind="box",
                        palette=self.colormap, order=order)
        [plt.setp(ax.texts, text="") for ax in g.axes.flat]
        g.set_titles(row_template='{row_name}', col_template='{col_name}')
        for axes in g.axes.flat:
            axes.set_xticklabels(axes.get_xticklabels(), rotation=65,
                                 horizontalalignment='right')

        g.savefig(os.path.join(self.outdir, self.get_basename(self.compare_path) + "_catplot.{}".format(self.extension)))
        plt.close()

if __name__ == '__main__':
    m = main()
    m.start()
