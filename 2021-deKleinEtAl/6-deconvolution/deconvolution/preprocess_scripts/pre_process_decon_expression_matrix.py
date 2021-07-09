#!/usr/bin/env python3

"""
File:         pre_process_decon_expression_matrix.py
Created:      2021/07/06
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
from __future__ import print_function
from pathlib import Path
import argparse
import time
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Pre-process Decon-eQTL Expression Matrix"
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

"""
Syntax: 
./pre_process_decon_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR -p GTE-EUR- -of CortexEUR

./pre_process_decon_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR -p GTE-EUR- -e ENA GVEX -of CortexEUR_noENA_noGVEX

./pre_process_decon_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR -p GTE-AFR- -e ENA -of CortexAFR_noENA
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.tcov_path = getattr(arguments, 'technical_covariates')
        self.gte_path = getattr(arguments, 'gene_to_expression')
        self.gte_prefix = getattr(arguments, 'gte_prefix')
        self.exclude = getattr(arguments, 'exclude')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        if outdir is None:
            outdir = str(Path(__file__).parent.parent)
        self.plot_outdir = os.path.join(outdir, 'pre_process_decon_expression_matrix', outfolder, 'plot')
        self.file_outdir = os.path.join(outdir, 'pre_process_decon_expression_matrix', outfolder, 'data')
        for outdir in [self.plot_outdir, self.file_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        self.file_cohort_dict = {
            "AMPAD-MAYO-V2": "MAYO",
            "CMC_HBCC_set2": "CMC HBCC",
            "GTEx": "GTEx",
            "AMPAD-ROSMAP-V2": "ROSMAP",
            "BrainGVEX-V2": "Brain GVEx",
            "TargetALS": "Target ALS",
            "AMPAD-MSBB-V2": "MSBB",
            "NABEC-H610": "NABEC",
            "LIBD_1M": "LIBD",
            "ENA": "ENA",
            "LIBD_h650": "LIBD",
            "GVEX": "GVEX",
            "NABEC-H550": "NABEC",
            "CMC_HBCC_set3": "CMC HBCC",
            "UCLA_ASD": "UCLA ASD",
            "CMC": "CMC",
            "CMC_HBCC_set1": "CMC HBCC"
        }

        self.palette = {
            "MAYO": "#9c9fa0",
            "CMC HBCC": "#0877b4",
            "GTEx": "#0fa67d",
            "ROSMAP": "#6950a1",
            "Brain GVEx": "#48b2e5",
            "Target ALS": "#d5c77a",
            "MSBB": "#5cc5bf",
            "NABEC": "#6d743a",
            "LIBD": "#e49d26",
            "ENA": "#d46727",
            "GVEX": "#000000",
            "UCLA ASD": "#f36d2a",
            "CMC": "#eae453",
            "NA": "#808080"
        }

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
                            help="show program's version number and exit.")
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the data matrix.")
        parser.add_argument("-t",
                            "--technical_covariates",
                            type=str,
                            required=True,
                            help="The path to the technical covariates matrix.")
        parser.add_argument("-gte",
                            "--gene_to_expression",
                            type=str,
                            required=True,
                            help="The path to the gene-expression link files.")
        parser.add_argument("-p",
                            "--gte_prefix",
                            type=str,
                            required=True,
                            help="The gene-expression link file prefix.")
        parser.add_argument("-e",
                            "--exclude",
                            nargs="*",
                            type=str,
                            default=[],
                            help="The gene-expression link files to exclude.")
        parser.add_argument("-od",
                            "--outdir",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of the output path.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # Construct the output filename.
        filename = os.path.basename(self.data_path).replace(".gz", "").replace(".txt", "")

        # Loading samples.
        print("Loading samples.")
        gte_combined_df = None
        cohort_to_sample = {}
        for infile in glob.glob(os.path.join(self.gte_path, "{}*.txt".format(self.gte_prefix))):
            file = os.path.basename(infile).replace(".txt", "").replace(self.gte_prefix, "")
            if file in self.exclude:
                continue
            gte_df = self.load_file(infile, header=None, index_col=None)
            gte_df["file"] = file
            if gte_combined_df is None:
                gte_combined_df = gte_df
            else:
                gte_combined_df = pd.concat([gte_combined_df, gte_df], axis=0, ignore_index=True)

            cohort_to_sample[file] = gte_df.iloc[:, 1]
        gte_combined_df["cohort"] = gte_combined_df.iloc[:, 2].map(self.file_cohort_dict)
        sample_to_cohort = dict(zip(gte_combined_df.iloc[:, 1], gte_combined_df.iloc[:, 3]))
        samples = gte_combined_df.iloc[:, 1].values.tolist()
        print("\tN samples: {}".format(len(samples)))

        # Create cohort matrix.
        cohort_sample_counts = list(zip(*np.unique(gte_combined_df["file"], return_counts=True)))
        cohort_sample_counts.sort(key=lambda x: -x[1])
        cohorts = [csc[0] for csc in cohort_sample_counts]
        print("\tCohorts: {} [N = {}]".format(", ".join(cohorts), len(cohorts)))

        cohort_df = pd.DataFrame(0, index=samples, columns=cohorts)
        for cohort in cohorts:
            cohort_df.loc[cohort_to_sample[cohort], cohort] = 1
        cohort_df.index.name = "-"

        # Load data.
        print("Loading data.")
        df = self.load_file(self.data_path, header=0, index_col=0)
        print(df)

        print("Step 1: sample selection.")
        print("\tUsing {}/{} samples.".format(len(samples), df.shape[1]))
        df = df.loc[:, samples]

        print("Step 2: remove probes with zero variance.")
        mask = df.std(axis=1) != 0
        print("\tUsing {}/{} probes.".format(np.sum(mask), np.size(mask)))
        df = df.loc[mask, :]

        print("Step 3: log2 transform.")
        min_value = df.min(axis=1).min()
        if min_value <= 0:
            df = np.log2(df - min_value + 1)
        else:
            df = np.log2(df + 1)

        print("Step 4: PCA analysis.")
        self.pca(df=df,
                 sample_to_cohort=sample_to_cohort,
                 plot_appendix="_1_CovariatesRemovedOLS")

        print("Step 4: Construct technical covariate matrix.")
        tcov_df = self.load_file(self.tcov_path, header=0, index_col=0)
        tcov_df = tcov_df.loc[samples, :]
        tcov_df = tcov_df.loc[:, tcov_df.std(axis=0) != 0]

        # replace cohort variables.
        drop_cols = []
        for col in tcov_df.columns:
            if set(tcov_df[col].unique()) == {0, 1}:
                drop_cols.append(col)
        tcov_df.drop(drop_cols, axis=1, inplace=True)
        tcov_df = tcov_df.merge(cohort_df.iloc[:, 1:], left_index=True, right_index=True)

        # filter on VIF.
        pre_columns = tcov_df.columns.tolist()
        tcov_df = self.remove_multicollinearity(tcov_df)
        print("Dropped:\t{}".format(", ".join([x for x in pre_columns if x not in tcov_df.columns])))
        #tcov_df = tcov_df.drop(['PCT_CODING_BASES', 'PCT_MRNA_BASES', 'PF_HQ_ALIGNED_READS'], axis=1)

        # add intercept.
        tcov_df.insert(0, "INTERCEPT", 1)
        print(tcov_df)
        print(tcov_df.columns)

        print("Step 5: remove covariates OLS.")
        corrected_m = np.empty_like(df, dtype=np.float64)
        last_print_time = None
        n_tests = df.shape[0]
        for i in range(n_tests):
            now_time = int(time.time())
            if last_print_time is None or (now_time - last_print_time) >= 10 or (i + 1) == n_tests:
                last_print_time = now_time
                print("\t{}/{} ieQTLs analysed [{:.2f}%]".format((i + 1), n_tests, (100 / n_tests) * (i + 1)))

            ols = OLS(df.iloc[i, :], tcov_df)
            results = ols.fit()
            # print(results.summary())
            corrected_m[i, :] = results.resid

        corrected_df = pd.DataFrame(corrected_m, index=df.index, columns=df.columns)
        del df

        print("Step 6: PCA analysis.")
        self.pca(df=corrected_df,
                 sample_to_cohort=sample_to_cohort,
                 plot_appendix="_2_Log2Transformed_CovariatesRemovedOLS")

        print("Step 7: exp added.")
        corrected_df = np.power(2, corrected_df)
        print("\tSaving file.")
        self.save_file(df=corrected_df, outpath=os.path.join(self.file_outdir, "{}.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt.gz".format(filename)))

        print("Step 8: PCA analysis.")
        self.pca(df=corrected_df,
                 sample_to_cohort=sample_to_cohort,
                 plot_appendix="_3_Log2Transformed_CovariatesRemovedOLS_ExpAdded")

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    @staticmethod
    def reverse_dict(dict):
        out_dict = {}
        seen_keys = set()
        for key, value in dict.items():
            if key in seen_keys:
                print("Key {} has muiltiple values.".format(key))
            seen_keys.add(key)

            if value in out_dict.keys():
                keys = out_dict[value]
                keys.append(key)
                out_dict[value] = keys
            else:
                out_dict[value] = [key]

        return out_dict

    def remove_multicollinearity(self, df, threshold=0.9999):
        indices = np.arange(df.shape[1])
        max_vif = np.inf
        while len(indices) > 1 and max_vif > threshold:
            vif = np.array([self.calc_ols_rsquared(df.iloc[:, indices], ix) for ix in range(len(indices))])
            max_vif = max(vif)

            if max_vif > threshold:
                max_index = np.where(vif == max_vif)[0][0]
                indices = np.delete(indices, max_index)

        return df.iloc[:, indices]

    @staticmethod
    def calc_ols_rsquared(df, idx):
        return OLS(df.iloc[:, idx], df.loc[:, np.arange(df.shape[1]) != idx]).fit().rsquared

    def pca(self, df, sample_to_cohort, plot_appendix=""):
        # samples should be on the columns and genes on the rows.
        zscores = (df - df.mean(axis=0)) / df.std(axis=0)
        pca = PCA(n_components=2)
        pca.fit(zscores)
        expl_variance = {"PC{}".format(i+1): pca.explained_variance_ratio_[i] * 100 for i in range(2)}
        components_df = pd.DataFrame(pca.components_)
        components_df.index = ["Comp{}".format(i + 1) for i, _ in enumerate(components_df.index)]
        components_df.columns = df.columns

        print("Plotting PCA")
        plot_df = components_df.T
        plot_df["cohort"] = plot_df.index.map(sample_to_cohort)
        plot_df["cohort"] = plot_df["cohort"].fillna('NA')
        self.plot(df=plot_df, x="Comp1", y="Comp2", hue="cohort", palette=self.palette,
                  xlabel="PC1 [{:.2f}%]".format(expl_variance["PC1"]),
                  ylabel="PC2 [{:.2f}%]".format(expl_variance["PC2"]),
                  title="PCA - eigenvectors",
                  filename="eigenvectors_plot{}".format(plot_appendix))

    def plot(self, df, x="x", y="y", hue=None, palette=None, xlabel=None,
             ylabel=None, title="", filename="PCA_plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, gridspec_kw={"width_ratios": [0.9, 0.1]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        s=100,
                        linewidth=0,
                        legend=None,
                        palette=palette,
                        ax=ax1)

        ax1.set_title(title,
                      fontsize=20,
                      fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_xlabel(xlabel,
                       fontsize=14,
                       fontweight='bold')

        if palette is not None:
            handles = []
            for label, color in palette.items():
                if label in df[hue].values.tolist():
                    handles.append(mpatches.Patch(color=color, label=label))
            ax2.legend(handles=handles, loc="center")

        #fig.savefig(os.path.join(self.plot_outdir, "{}.pdf".format(filename)))
        fig.savefig(os.path.join(self.plot_outdir, "{}.png".format(filename)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Data: {}".format(self.data_path))
        print("  > Technical covariates: {}".format(self.tcov_path))
        print("  > GtE path: {}".format(self.gte_path))
        print("  >   GtE prefix: {}".format(self.gte_prefix))
        print("  >   Exclude: {}".format(self.exclude))
        print("  > Plot output directory {}".format(self.plot_outdir))
        print("  > File output directory {}".format(self.file_outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
