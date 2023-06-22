#!/usr/bin/env python3

"""
File:         sn_replication.py
Created:      2022/02/10
Last Changed: 2022/05/26
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
import argparse
import os

# Third party imports.
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from adjustText import adjust_text

# Local application imports.

"""
Syntax:
./sn_replication.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron/merged_decon_results.txt.gz \
    -e png pdf
    
./sn_replication.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -s 2022-03-03-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected \
    -t ieQTL replication in ROSMAP snRNA-seq \
    -et cis \
    -e png pdf
    
### trans

./sn_replication.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -s 2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected \
    -t CortexEUR-and-AFR-noENA-noAMPAD-trans-0PCs \
    -et trans \
    -e png pdf
    
./sn_replication.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -s 2022-03-31-CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs-NegativeToZero-DatasetAndRAMCorrected \
    -t CortexEUR-and-AFR-noENA-noAMPAD-trans-80PCs \
    -et trans \
    -e png pdf
    
./sn_replication.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -s 2022-03-31-CortexEUR-and-AFR-noENA-trans-0PCs-NegativeToZero-DatasetAndRAMCorrected \
    -t CortexEUR-and-AFR-noENA-trans-0PCs \
    -et trans \
    -e png pdf
    
./sn_replication.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected/merged_decon_results.txt.gz \
    -s 2022-03-31-CortexEUR-and-AFR-noENA-trans-100PCs-NegativeToZero-DatasetAndRAMCorrected \
    -t CortexEUR-and-AFR-noENA-trans-100PCs \
    -et trans \
    -e png pdf
"""

# Metadata
__program__ = "Decon-eQTL Single-Nucleus replication"
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
        arguments = self.create_argument_parser()
        self.discovery_path = getattr(arguments, 'discovery_path')
        self.subdir = getattr(arguments, 'subdir')
        self.title = " ".join(getattr(arguments, 'title'))
        self.eqtl_type = getattr(arguments, 'eqtl_type')
        self.extensions = getattr(arguments, 'extension')

        # Define the replication data.
        self.sn_eqtl_infolder = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2021-12-23-ROSMAP-scRNAseq/{}_100Perm_ieQTLs/".format(self.eqtl_type)
        self.sn_celltypes = [
            ("AST", "Astrocyte"),
            ("END", "EndothelialCell"),
            ("EX", "Excitatory"),
            ("IN", "Inhibitory"),
            ("MIC", "Microglia"),
            ("OLI", "Oligodendrocyte"),
            ("OPC", "OPC"),
            ("PER", "Pericytes")
        ]
        self.sn_eqtl_filename = "eQTLsFDR-ProbeLevel.txt.gz"
        if self.eqtl_type == "trans":
            self.sn_eqtl_filename = "eQTLsFDR.txt.gz"

        # Set variables.
        self.outdir = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), 'sn_replication', self.eqtl_type)
        if self.subdir:
            self.outdir = os.path.join(self.outdir, self.subdir)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
            "Astrocyte": "#D55E00",
            "EndothelialCell": "#CC79A7",
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Microglia": "#E69F00",
            "Oligodendrocyte": "#009E73",
            "OtherNeuron": "#2690ce"
        }

        self.shared_xlim = None
        self.shared_ylim = None

        matplotlib.rcParams['pdf.fonttype'] = 42

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
        parser.add_argument("-d",
                            "--discovery_path",
                            type=str,
                            required=True,
                            help="The path to the discovery deconvolution "
                                 "results matrix")
        parser.add_argument("-s",
                            "--subdir",
                            type=str,
                            default=None,
                            help="The output subdirectory. Default: None.")
        parser.add_argument("-t",
                            "--title",
                            nargs="*",
                            type=str,
                            default="ieQTL replication in ROSMAP snRNA-seq",
                            help="The plot title. Default: 'ieQTL replication "
                                 "in ROSMAP snRNA-seq'.")
        parser.add_argument("-et",
                            "--eqtl_type",
                            type=str,
                            choices=["cis", "trans"],
                            default="cis",
                            help="The type of eQTLs to map. "
                                 "Default: 'cis'.")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading discovery data")
        discovery_df = self.load_file(self.discovery_path, header=0, index_col=None)
        discovery_df.index = discovery_df["Gene"] + "_" + discovery_df["SNP"]
        print(discovery_df)

        # Select columns
        columns_of_interest = ["Gene", "Gene symbol", "SNP", "Alleles", "Allele assessed", "N", "HW pval", "Minor allele", "MAF", "Overall z-score"]
        colnames = columns_of_interest.copy()
        cell_types = [x.split(" ")[0] for x in discovery_df.columns if "pvalue" in x]
        exclude_in_excel = ["N", "HW pval", "Minor allele", "MAF", "Overall z-score"]
        for ct in cell_types:
            columns_of_interest.append("{} pvalue".format(ct))
            colnames.append("Bulk {} pvalue".format(ct))

            columns_of_interest.append("{} BH-FDR".format(ct))
            colnames.append("Bulk {} FDR".format(ct))

            columns_of_interest.append("{} interaction beta".format(ct))
            colnames.append("Bulk {} interaction beta".format(ct))
        discovery_df = discovery_df.loc[:, columns_of_interest].copy()
        discovery_df.columns = colnames
        print(discovery_df)

        # Create a bulk aa translate dict.
        bulk_aa_dict = dict(zip(discovery_df.index, discovery_df["Allele assessed"]))

        print("Loading replication data")
        replication_df_list = []
        for (abbreviation, full_name) in self.sn_celltypes:
            replication_ct_df = self.load_file(os.path.join(self.sn_eqtl_infolder, abbreviation, self.sn_eqtl_filename),
                                               header=0,
                                               index_col=None)

            # Split the meta beta column.
            beta_values = []
            se_values = []
            for _, row in replication_ct_df.iterrows():
                value = row["Meta-Beta (SE)"]
                beta_values.append(float(value.split(" (")[0]))
                se_values.append(float(value.split(" (")[1].replace(")", "")))
            replication_ct_df["Meta-Beta"] = beta_values
            replication_ct_df["Meta-SE"] = se_values

            # Select the columns of interest.
            replication_ct_df = replication_ct_df[["ProbeName",
                                                   "SNPName",
                                                   "AlleleAssessed",
                                                   "PValue",
                                                   "FDR",
                                                   "Meta-Beta",
                                                   "Meta-SE"]]
            replication_ct_df.columns = ["ProbeName",
                                         "SNPName",
                                         "replication_aa",
                                         "SN {} eQTL pvalue".format(full_name),
                                         "SN {} eQTL FDR".format(full_name),
                                         "SN {} eQTL beta".format(full_name),
                                         "SN {} eQTL se".format(full_name)]
            replication_ct_df.index = replication_ct_df["ProbeName"] + "_" + replication_ct_df["SNPName"]

            # Flip the z-scores.
            replication_ct_df["discovery_aa"] = replication_ct_df.index.map(bulk_aa_dict)
            replication_ct_df.dropna(inplace=True)

            replication_ct_df["flip"] = replication_ct_df["replication_aa"] != replication_ct_df["discovery_aa"]
            replication_ct_df.loc[:, "SN {} eQTL beta".format(full_name)] = replication_ct_df["SN {} eQTL beta".format(full_name)] * replication_ct_df["flip"].map({True: -1, False: 1})

            # save.
            replication_df_list.append(replication_ct_df.loc[:, ["SN {} eQTL pvalue".format(full_name),
                                                                 "SN {} eQTL FDR".format(full_name),
                                                                 "SN {} eQTL beta".format(full_name),
                                                                 "SN {} eQTL se".format(full_name)]].copy())
            exclude_in_excel.append("SN {} eQTL se".format(full_name))

        replication_df = pd.concat(replication_df_list, axis=1)

        print("Merging tables")
        df = discovery_df.merge(replication_df, left_index=True, right_index=True, how="left")
        print(df)

        print("Saving output")
        self.save_file(df=df,
                       outpath=os.path.join(self.outdir, "single_nucleus_replication.txt.gz"),
                       index=False)
        self.save_file(df=df.loc[:, [col for col in df.columns if col not in exclude_in_excel]],
                       outpath=os.path.join(self.outdir, "single_nucleus_replication.xlsx"),
                       index=False,
                       sheet_name="ROSMAP Single-Nucleus")

        # df = self.load_file(os.path.join(self.outdir, "single_nucleus_replication.txt.gz"),
        #                     header=0,
        #                     index_col=False)

        print("Visualizing")
        discovery_ct = set([x.split(" ")[1] for x in df.columns if "Bulk" in x and "FDR" in x])
        replication_ct = set([x.split(" ")[1] for x in df.columns if "SN" in x and "FDR" in x])
        overlap_ct = list(discovery_ct.intersection(replication_ct))
        overlap_ct.sort()

        replication_stats_df = self.plot(df=df,
                                         cell_types=overlap_ct,
                                         title=self.title)
        self.save_file(df=replication_stats_df,
                       outpath=os.path.join(self.outdir,
                                            "replication_stats.txt.gz"))

        # replication_stats_df = self.load_file(os.path.join(self.outdir,  "replication_stats.txt.gz"),
        #                                       header=0,
        #                                       index_col=0)

        print("Replication stats")
        for label in replication_stats_df["label"].unique():
            print("\t{}".format(label))
            stats_df = replication_stats_df.loc[replication_stats_df["label"] == label, :]
            stats_df_mean = stats_df[["variable", "value"]].groupby("variable").mean()
            for index, row in stats_df_mean.iterrows():
                print("\t  {}: {:.2f}".format(index, row["value"]))

            stats_df_sum = stats_df[["variable", "value"]].groupby("variable").sum()
            print("\t  Overall concordance: {:,}/{:,} [{:.2f}%]".format(stats_df_sum.loc["N concordant", "value"],
                                                                        stats_df_sum.loc["N", "value"],
                                                                        (100 / stats_df_sum.loc["N", "value"]) * stats_df_sum.loc["N concordant", "value"]))
            print("")

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
    def save_file(df, outpath, header=True, index=True, sep="\t", na_rep="NA",
                  sheet_name="Sheet1"):
        if outpath.endswith('xlsx'):
            df.to_excel(outpath,
                        sheet_name=sheet_name,
                        na_rep=na_rep,
                        header=header,
                        index=index)
        else:
            compression = 'infer'
            if outpath.endswith('.gz'):
                compression = 'gzip'

            df.to_csv(outpath, sep=sep, index=index, header=header,
                      compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def plot(self, df, cell_types, title=""):
        nrows = 3
        ncols = len(cell_types)

        self.shared_ylim = {i: (0, 1) for i in range(nrows)}
        self.shared_xlim = {i: (0, 1) for i in range(ncols)}

        replication_stats = []

        sns.set(rc={'figure.figsize': (ncols * 8, nrows * 6)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='col',
                                 sharey='row')

        for col_index, ct in enumerate(cell_types):
            print("\tWorking on '{}'".format(ct))

            # Select the required columns.
            plot_df = df.loc[:, ["Gene symbol",
                                 "N",
                                 "MAF",
                                 "Bulk {} pvalue".format(ct),
                                 "Bulk {} FDR".format(ct),
                                 "Bulk {} interaction beta".format(ct),
                                 "SN {} eQTL pvalue".format(ct),
                                 "SN {} eQTL FDR".format(ct),
                                 "SN {} eQTL beta".format(ct),
                                 "SN {} eQTL se".format(ct)]].copy()
            plot_df.columns = ["Gene symbol",
                               "bulk N",
                               "bulk MAF",
                               "bulk pvalue",
                               "bulk FDR",
                               "bulk interaction beta",
                               "SN pvalue",
                               "SN FDR",
                               "SN beta",
                               "SN se"]
            plot_df.dropna(inplace=True)
            plot_df.sort_values(by="bulk pvalue", inplace=True)

            # Calculate the discovery standard error.
            self.pvalue_to_zscore(df=plot_df,
                                  beta_col="bulk interaction beta",
                                  p_col="bulk pvalue",
                                  prefix="bulk ")
            self.zscore_to_beta(df=plot_df,
                                z_col="bulk z-score",
                                maf_col="bulk MAF",
                                n_col="bulk N",
                                prefix="bulk zscore-to-")

            # Convert the interaction beta to log scale.
            plot_df["log bulk interaction beta"] = self.log_modulus_beta(plot_df["bulk interaction beta"])
            plot_df["log SN eqtl beta"] = self.log_modulus_beta(plot_df["SN beta"])

            # Add the hue for the significant replicating ieQTLs.
            plot_df["hue"] = "#808080"
            plot_df.loc[(plot_df["bulk FDR"] <= 0.05) & (plot_df["SN FDR"] <= 0.05), "hue"] = self.palette[ct]
            print(plot_df)

            include_ylabel = False
            if col_index == 0:
                include_ylabel = True

            if col_index == 0:
                for row_index, panel in enumerate(["A", "B", "C"]):
                    axes[row_index, col_index].annotate(
                        panel,
                        xy=(-0.3, 0.9),
                        xycoords=axes[row_index, col_index].transAxes,
                        color="#000000",
                        fontsize=40
                    )

            print("\tPlotting row 1.")
            xlim, ylim, stats1 = self.scatterplot(
                df=plot_df,
                fig=fig,
                ax=axes[0, col_index],
                x="log bulk interaction beta",
                y="log SN eqtl beta",
                xlabel="",
                ylabel="ROSMAP SN log eqtl beta",
                title=ct,
                color=self.palette[ct],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 0, col_index)

            print("\tPlotting row 2.")
            xlim, ylim, stats2 = self.scatterplot(
                df=plot_df.loc[plot_df["bulk FDR"] <= 0.05, :],
                fig=fig,
                ax=axes[1, col_index],
                x="log bulk interaction beta",
                y="log SN eqtl beta",
                xlabel="",
                ylabel="ROSMAP SN log eqtl beta",
                title="",
                color=self.palette[ct],
                include_ylabel=include_ylabel,
                pi1_column="SN pvalue",
                rb_columns=[("bulk zscore-to-beta", "bulk zscore-to-se"), ("SN beta", "SN se")]
            )
            self.update_limits(xlim, ylim, 1, col_index)

            print("\tPlotting row 3.")
            xlim, ylim, stats3 = self.scatterplot(
                df=plot_df.loc[(plot_df["bulk FDR"] <= 0.05) & (plot_df["SN FDR"] <= 0.05), :],
                fig=fig,
                ax=axes[2, col_index],
                x="log bulk interaction beta",
                y="log SN eqtl beta",
                label="Gene symbol",
                xlabel="MetaBrain log interaction beta",
                ylabel="ROSMAP SN log eqtl beta",
                title="",
                color=self.palette[ct],
                include_ylabel=include_ylabel
            )
            self.update_limits(xlim, ylim, 2, col_index)
            print("")

            for stats, label in zip([stats1, stats2, stats3], ["all", "discovery significant", "both significant"]):
                stats_m = stats.melt()
                stats_m["label"] = label
                stats_m["cell type"] = ct
                replication_stats.append(stats_m)

        for (m, n), ax in np.ndenumerate(axes):
            (xmin, xmax) = self.shared_xlim[n]
            (ymin, ymax) = self.shared_ylim[m]

            xmargin = (xmax - xmin) * 0.05
            ymargin = (ymax - ymin) * 0.05

            ax.set_xlim(xmin - xmargin - 1, xmax + xmargin)
            ax.set_ylim(ymin - ymargin, ymax + ymargin)

        # Add the main title.
        fig.suptitle(title,
                     fontsize=40,
                     color="#000000",
                     weight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "sn_replication_plot.{}".format(extension)))
        plt.close()

        # Construct the replication stats data frame.
        replication_stats_df = pd.concat(replication_stats, axis=0)
        replication_stats_df.dropna(inplace=True)

        return replication_stats_df

    @staticmethod
    def pvalue_to_zscore(df, beta_col, p_col, prefix=""):
        p_values = df[p_col].to_numpy()
        zscores = stats.norm.ppf(p_values / 2)
        mask = np.ones_like(p_values)
        mask[df[beta_col] > 0] = -1
        df["{}z-score".format(prefix)] = zscores * mask
        df.loc[df[p_col] == 1, "{}z-score".format(prefix)] = 0
        df.loc[df[p_col] == 0, "{}z-score".format(prefix)] = -40.

    @staticmethod
    def zscore_to_beta(df, z_col, maf_col, n_col, prefix=""):
        chi = df[z_col] * df[z_col]
        a = 2 * df[maf_col] * (1 - df[maf_col]) * (df[n_col] + chi)
        df["{}beta".format(prefix)] = df[z_col] / a ** (1/2)
        df["{}se".format(prefix)] = 1 / a ** (1/2)

    @staticmethod
    def log_modulus_beta(series):
        s = series.copy()
        data = []
        for index, beta in s.T.iteritems():
            data.append(np.log(abs(beta) + 1) * np.sign(beta))
        new_df = pd.Series(data, index=s.index)

        return new_df

    def scatterplot(self, df, fig, ax, x="x", y="y", facecolors=None,
                    label=None, max_labels=15, xlabel="", ylabel="", title="",
                    color="#000000", ci=95, include_ylabel=True,
                    pi1_column=None, rb_columns=None):
        sns.despine(fig=fig, ax=ax)

        if not include_ylabel:
            ylabel = ""

        if facecolors is None:
            facecolors = "#808080"
        else:
            facecolors = df[facecolors]

        n = df.shape[0]
        concordance = np.nan
        n_concordant = np.nan
        coef = np.nan
        pi1 = np.nan
        rb = np.nan

        if n > 0:
            lower_quadrant = df.loc[(df[x] < 0) & (df[y] < 0), :]
            upper_quadrant = df.loc[(df[x] > 0) & (df[y] > 0), :]
            n_concordant = lower_quadrant.shape[0] + upper_quadrant.shape[0]
            concordance = (100 / n) * n_concordant

            if n > 1:
                coef, p = stats.pearsonr(df[x], df[y])

                if pi1_column is not None:
                    pi1 = self.calculate_p1(p=df[pi1_column])

                if rb_columns is not None and n > 2:
                    rb_est = self.calculate_rb(
                        b1=df[rb_columns[0][0]],
                        se1=df[rb_columns[0][1]],
                        b2=df[rb_columns[1][0]],
                        se2=df[rb_columns[1][1]],
                    )
                    rb = rb_est[0]

            sns.regplot(x=x, y=y, data=df, ci=ci,
                        scatter_kws={'facecolors': facecolors,
                                     'edgecolors': "#808080"},
                        line_kws={"color": color},
                        ax=ax
                        )

            if label is not None:
                texts = []
                for i, (_, point) in enumerate(df.iterrows()):
                    if i > max_labels:
                        continue
                    texts.append(ax.text(point[x],
                                         point[y],
                                         str(point[label]),
                                         color=color))
                adjust_text(texts,
                            ax=ax,
                            only_move={'points': 'x',
                                       'text': 'xy',
                                       'objects': 'x'},
                            autoalign='x',
                            expand_text=(1., 1.),
                            expand_points=(1., 1.),
                            arrowprops=dict(arrowstyle='-', color='#808080'))

        ax.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        ax.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)

        y_pos = 0.9
        if n > 0:
            ax.annotate(
                'N = {:,}'.format(n),
                xy=(0.03, 0.9),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )
            y_pos -= 0.05

        if not np.isnan(coef):
            ax.annotate(
                'r = {:.2f}'.format(coef),
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )
            y_pos -= 0.05

        if not np.isnan(concordance):
            ax.annotate(
                'concordance = {:.0f}%'.format(concordance),
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )
            y_pos -= 0.05

        if not np.isnan(pi1):
            ax.annotate(
                '\u03C01 = {:.2f}'.format(pi1),
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )
            y_pos -= 0.05

        if not np.isnan(rb):
            ax.annotate(
                'Rb = {:.2f}'.format(rb),
                xy=(0.03, y_pos),
                xycoords=ax.transAxes,
                color=color,
                fontsize=14,
                fontweight='bold'
            )

        ax.set_title(title,
                     fontsize=22,
                     color=color,
                     weight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        stats_df = pd.DataFrame([[n, n_concordant, concordance, coef, pi1, rb]],
                                columns=["N", "N concordant", "concordance", "pearsonr", "pi1", "Rb"],
                                index=[0])

        return (df[x].min(), df[x].max()), (df[y].min(), df[y].max()), stats_df

    def update_limits(self, xlim, ylim, row, col):
        row_ylim = self.shared_ylim[row]
        if ylim[0] < row_ylim[0]:
            row_ylim = (ylim[0], row_ylim[1])
        if ylim[1] > row_ylim[1]:
            row_ylim = (row_ylim[0], ylim[1])
        self.shared_ylim[row] = row_ylim

        col_xlim = self.shared_xlim[col]
        if xlim[0] < col_xlim[0]:
            col_xlim = (xlim[0], col_xlim[1])
        if xlim[1] > col_xlim[1]:
            col_xlim = (col_xlim[0], xlim[1])
        self.shared_xlim[col] = col_xlim

    @staticmethod
    def calculate_p1(p):
        robjects.r("source('qvalue_truncp.R')")
        p = robjects.FloatVector(p)
        qvalue_truncp = robjects.globalenv['qvalue_truncp']
        pi0 = qvalue_truncp(p)[0]
        return 1 - np.array(pi0)

    @staticmethod
    def calculate_rb(b1, se1, b2, se2, theta=0):
        robjects.r("source('Rb.R')")
        b1 = robjects.FloatVector(b1)
        se1 = robjects.FloatVector(se1)
        b2 = robjects.FloatVector(b2)
        se2 = robjects.FloatVector(se2)
        calcu_cor_true = robjects.globalenv['calcu_cor_true']
        rb = calcu_cor_true(b1, se1, b2, se2, theta)
        return np.array(rb)[0]

    def print_arguments(self):
        print("Arguments:")
        print("  > Discovery:")
        print("    > File path: {}".format(self.discovery_path))
        print("  > Replication:")
        print("    > Input folder: {}".format(self.sn_eqtl_infolder))
        print("    > Cell types: {}".format(", ".join([x[0] for x in self.sn_celltypes])))
        print("    > Input filename: {}".format(self.sn_eqtl_filename))
        print("  > Title: {}".format(self.title))
        print("  > eQTL type: {}".format(self.eqtl_type))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Extensions: {}".format(self.extensions))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
