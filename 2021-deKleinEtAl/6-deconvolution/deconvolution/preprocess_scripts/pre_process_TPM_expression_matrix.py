#!/usr/bin/env python3

"""
File:         pre_process_TPM_expression_matrix.py
Created:      2021/11/23
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
import glob
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Local application imports.

# Metadata
__program__ = "Pre-process TPM Expression Matrix"
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
# Calculon
./pre_process_TPM_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-01-31-raw-count-tables/ -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/gencode.v32.primary_assembly.annotation-genelengths.txt.gz -ms /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_all.txt.gz -ps /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/PsychENCODE_samples.txt.gz -o MetaBrain

./pre_process_TPM_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-01-31-raw-count-tables/ -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/gencode.v32.primary_assembly.annotation-genelengths.txt.gz -ms /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain_STD_cortex.txt.gz -ps /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/PsychENCODE_samples.txt.gz -o MetaBrain_Cortex

# Datasets
TargetALS = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/TargetALS/pipelines/no_patch_chromosomes/results/TargetALS.geneCounts.2019-12-13.txt
Braineac = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/ucl-upload-biogen/pipelines/no_patch_chromosomes/results/Braineac.geneCounts.txt.gz
GTEx = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/GTEx/pipelines/no_patch_chromosomes/results/GTEx.geneCounts.txt 
NABEC = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/NABEC/pipelines/no_patch_chromosomes/results/NABEC.geneCounts.txt.gz 
ENA = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/ENACounts/ENA.geneCounts.txt.gz
BrainGVEx = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/pipelines/no_patch_chromosomes/BrainGVEx/results/BrainGVEx.geneCounts.2019-12-16.txt 
BipSeq = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/pipelines/no_patch_chromosomes/BipSeq/results/BipSeq.geneCounts.2019-12-17.txt 
UCLA_ASD = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/pipelines/no_patch_chromosomes/UCLA_ASD/results/UCLA_ASD.geneCounts.2019-12-15.txt 
CMC_HBCC = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/pipelines/no_patch_chromosomes/CMC_HBCC/results/CMC_HBCC.geneCounts.2020-01-02.txt
CMC = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/CMC/pipelines/no_patch_chromosomes/results/CMC.geneCounts.txt.gz 
MSBB = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/AMP_AD/pipelines/no_patch_chromosomes/MSBB/results/MSBB.geneCounts.2020-01-02.txt.gz
ROSMAP = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/AMP_AD/pipelines/no_patch_chromosomes/ROSMAP/results/ROSMAP.geneCounts.2020-01-02.txt.gz  
Mayo_CBE = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/AMP_AD/pipelines/no_patch_chromosomes/MayoCBE/results/MayoCBE.geneCounts.2020-01-06.txt.gz 
Mayo_TCX = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/AMP_AD/pipelines/no_patch_chromosomes/MayoTCX/results/MayoTCX.geneCounts.2020-01-06.txt.gz 
Brainseq = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/Brainseq/pipelines/no_patch_chromosomes/results/Brainseq.geneCounts.txt
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.gene_info_path = getattr(arguments, 'gene_info')
        self.mb_samples_path = getattr(arguments, 'metabrain_samples')
        self.pe_samples_path = getattr(arguments, 'psychencode_samples')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        outdir = str(Path(__file__).parent.parent)
        self.plot_outdir = os.path.join(outdir, 'pre_process_TPM_expression_matrix', outfolder, 'plot')
        self.file_outdir = os.path.join(outdir, 'pre_process_TPM_expression_matrix', outfolder, 'data')
        for outdir in [self.plot_outdir, self.file_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        self.palette = {
            "TargetALS": "#d5c77a",
            "Braineac": "#808080",
            "GTEx": "#0fa67d",
            "NABEC": "#6d743a",
            "ENA": "#d46727",
            "BrainGVEx": "#48b2e5",
            "BipSeq": "#000000",
            "UCLA_ASD": "#f36d2a",
            "CMC_HBCC": "#0877b4",
            "CMC": "#eae453",
            "MSBB": "#5cc5bf",
            "ROSMAP": "#6950a1",
            "Mayo_CBE": "#9c9fa0",
            "Mayo_TCX": "#9c9fa0",
            "Brainseq": "#e49d26"
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
        parser.add_argument("-g",
                            "--gene_info",
                            type=str,
                            required=True,
                            help="The path to the gene info matrix.")
        parser.add_argument("-ms",
                            "--metabrain_samples",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the MetaBrain samples.")
        parser.add_argument("-ps",
                            "--psychencode_samples",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the PsychENCODE samples.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        mb_samples = None
        if self.mb_samples_path is not None:
            print("Loading samples data.")
            mb_samples_df = self.load_file(self.mb_samples_path, header=0, index_col=None)
            mb_samples = set(mb_samples_df["rnaseq_id"].values.tolist())
            del mb_samples_df

        pe_samples = None
        if self.pe_samples_path is not None:
            print("Loading samples data.")
            pe_samples_df = self.load_file(self.pe_samples_path, header=0, index_col=None)
            pe_samples = set(pe_samples_df["PsychENCODE_id"].values.tolist())
            del pe_samples_df

        # gte_trans_df = self.load_file("/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/E-MTAB-5214.sdrf.txt", sep="\t", header=0, index_col=None)
        # gte_trans_df = gte_trans_df.loc[gte_trans_df["Comment[original body site annotation]"].isin(["Brain - Anterior cingulate cortex (BA24)",
        #                                                                                              "Brain - Caudate (basal ganglia)",
        #                                                                                              "Brain - Cerebellar Hemisphere",
        #                                                                                              "Brain - Cerebellum",
        #                                                                                              "Brain - Cortex",
        #                                                                                              "Brain - Frontal Cortex (BA9)",
        #                                                                                              "Brain - Hippocampus",
        #                                                                                              "Brain - Hypothalamus",
        #                                                                                              "Brain - Nucleus accumbens (basal ganglia)",
        #                                                                                              "Brain - Putamen (basal ganglia)"]), :]
        # gte_trans_df = gte_trans_df.loc[:, ["Scan Name", "Characteristics[individual]"]]
        # gte_trans_df.drop_duplicates(inplace=True)
        # gte_trans_df.index = gte_trans_df["Scan Name"]

        print("Loading data.")
        counts_df_list = []
        stats_df_list = []
        translate_data = []
        found_mb_samples = set()
        found_pe_samples = set()
        special_samples_trans_dict = {"AN11864_ba41-42-22": "AN11864_ba41.42.22",
                                      "UMB1376_ba41-42-22": "UMB1376_ba41.42.22"}
        for filepath in glob.glob(os.path.join(self.data_path, "*.txt.gz")):
            df = self.load_file(filepath, header=0, index_col=0)
            df.dropna(how="all", inplace=True)
            old_names = list(df.columns)
            dataset = None

            psychencode_ids = [None] * df.shape[1]

            filename = os.path.basename(filepath)
            if "TargetALS" in filename:
                df.columns = [re.sub("-", "_", colname) for colname in df.columns]
                df.columns = [re.sub("\\.", "_", colname) for colname in df.columns]
                df.columns = [re.search(".*(HRA_[0-9]+)", colname).group(1) for colname in df.columns]
                dataset = "TargetALS"
            elif "Braineac" in filename:
                df.columns = [re.search(".*(A653.*)", colname).group(1) for colname in df.columns]
                dataset = "Braineac"
            elif "GTEx" in filename:
                # psychencode_ids = [re.search("(.*)_.*", colname).group(1) for colname in df.columns]
                # gte_trans_df = gte_trans_df.loc[[colname for colname in psychencode_ids if colname in gte_trans_df.index], :]
                # gte_trans_dict = dict(zip(gte_trans_df["Scan Name"], gte_trans_df["Characteristics[individual]"]))
                # psychencode_ids = [gte_trans_dict[colname] if colname in gte_trans_dict else None for colname in psychencode_ids]
                # print(len([x for x in psychencode_ids if x is not None]))
                # print(len(set(psychencode_ids)))
                # exit()

                df.columns = [re.search("(.*)_.*", colname).group(1) for colname in df.columns]
                dataset = "GTEx"
            elif "NABEC" in filename:
                df.columns = [re.search(".*_(.*)", colname).group(1) for colname in df.columns]
                dataset = "NABEC"
            elif "ENA" in filename:
                df.columns = [re.search(".*_(.*)", colname).group(1) for colname in df.columns]
                dataset = "ENA"
            elif "BrainGVEx" in filename:
                df.columns = [re.sub("_", "-", colname) for colname in df.columns]
                df.columns = [re.sub("\\.", "-", colname) for colname in df.columns]
                df.columns = [re.sub("^X", "", colname) for colname in df.columns]

                psychencode_ids = [re.search("^([0-9]+-[0-9]+)-[0-9]+-[0-9]+", colname).group(1) for colname in df.columns]

                df.columns = [re.search("^([0-9]+-[0-9]+)-[0-9]+-[0-9]+", colname).group(1) for colname in df.columns]
                dataset = "BrainGVEx"
            elif "BipSeq" in filename:
                psychencode_ids = [re.search("(Br[0-9]+)_R[0-9]+", colname).group(1) for colname in df.columns]

                df.columns = [re.search("Br[0-9]+_(R[0-9]+)", colname).group(1) for colname in df.columns]
                dataset = "BipSeq"
            elif "UCLA_ASD" in filename:
                psychencode_ids = [re.search("([aA-zZ]+[0-9]+)_.+", colname).group(1) for colname in df.columns]

                df.columns = [re.search("[aA-zZ]+[0-9]+_(.+)", colname).group(1) for colname in df.columns]
                df.columns = [special_samples_trans_dict[colname] if colname in special_samples_trans_dict else colname for colname in df.columns]
                dataset = "UCLA_ASD"
            elif "CMC_HBCC" in filename:
                df.columns = [re.search("individualID.*_specimenID.(.*)", colname).group(1) for colname in df.columns]
                dataset = "CMC_HBCC"
            elif "CMC" in filename:
                psychencode_ids = [re.search("(^CMC_[aA-zZ]+_[0-9]+)_.*", colname).group(1) for colname in df.columns]

                df.columns = [re.search("^CMC_[aA-zZ]+_[0-9]+_(.*)", colname).group(1) for colname in df.columns]
                dataset = "CMC"
            elif "MSBB" in filename:
                df.columns = [re.sub(".accepted_hits.sort.coordReadsPerGene.out.tab", "", colname) for colname in df.columns]
                df.columns = [re.search("AMPAD_MSSM_[0-9]+_(.*)", colname).group(1) for colname in df.columns]
                dataset = "MSBB"
            elif "ROSMAP" in filename:
                df.columns = [re.sub("^X", "", colname) for colname in df.columns]
                df.columns = [re.sub("ReadsPerGene.out.tab", "", colname) for colname in df.columns]
                df.columns = [re.search(".*_.*_(.*_.*)", colname).group(1) for colname in df.columns]
                dataset = "ROSMAP"
            elif "MayoCBE" in filename:
                df.columns = [re.sub("^X", "", colname) for colname in df.columns]
                df.columns = [re.search("[0-9]+_CER", colname).group(0) for colname in df.columns]
                dataset = "Mayo_CBE"
            elif "MayoTCX" in filename:
                df.columns = [re.sub("^X", "", colname) for colname in df.columns]
                df.columns = [re.search("[0-9]+_TCX", colname).group(0) for colname in df.columns]
                dataset = "Mayo_TCX"
            elif "Brainseq" in filename:
                psychencode_ids = [re.search("(Br[0-9]+)_R[0-9]+", colname).group(1) for colname in df.columns]

                # these did not get adjusted in the other talbes, so keep the same
                dataset = "Brainseq"
            else:
                print("Unexpected input file.")
                exit()

            new_names = list(df.columns)
            translate_data.extend(list(zip([dataset] * df.shape[1], old_names, new_names, psychencode_ids)))

            if mb_samples is not None:
                print("\t  MetaBrain sample overlap: {} / {}".format(len(mb_samples.intersection(set(new_names))), df.shape[1]))
                found_mb_samples.update(set(new_names))
            if pe_samples is not None and dataset in ["GTEx", "BrainGVEx", "BipSeq", "UCLA_ASD", "CMC", "Brainseq"]:
                print("\t  PsychENCODE sample overlap: {} / {}".format(len(pe_samples.intersection(set(psychencode_ids))), df.shape[1]))
                found_pe_samples.update(set(psychencode_ids))

            if "ENA" in filename:
                # Ends with
                # __no_feature
                # __ambiguous
                # __too_low_aQual
                # __not_aligned
                # __alignment_not_unique
                stats_df = df.iloc[(df.shape[0] - 5):, :]
                stats_df.index = ["N_noFeature", "N_ambiguous", "N_too_low_aQual", "N_unmapped", "N_multimapping"]

                counts_df_list.append(df.iloc[:(df.shape[0] - 5), :])
                stats_df_list.append(stats_df)
            else:
                # Starts with
                # N_unmapped
                # N_multimapping
                # N_noFeature
                # N_ambiguous
                stats_df_list.append(df.iloc[:4, :])
                counts_df_list.append(df.iloc[4:, :])
                # counts_df_list.append(df)

        if mb_samples is not None:
            mb_missing = [sample for sample in mb_samples if sample not in found_mb_samples]
            print("\t  Missing MetaBrain samples [N={}]: {}".format(len(mb_missing), ", ".join(mb_missing)))
        if pe_samples is not None:
            pe_missing = [sample for sample in pe_samples if sample not in found_pe_samples]
            print("\t  Missing PsychENCODE samples [N={}]: {}".format(len(pe_missing), ", ".join(pe_missing)))

        counts_df = pd.concat(counts_df_list, axis=1)
        counts_df.fillna(0, inplace=True)
        print(counts_df)
        stats_df = pd.concat(stats_df_list, axis=1).T
        print(stats_df)
        translate_df = pd.DataFrame(translate_data, columns=["Dataset", "Original ID", "rnaseq_id", "psychencode_id"])
        print(translate_df)

        print("Subset samples.")
        overlap_samples = [sample for sample in counts_df.columns if (sample in mb_samples or sample in pe_samples)]
        counts_df = counts_df.loc[:, overlap_samples]
        stats_df = stats_df.loc[overlap_samples, :]

        print("Saving data")
        self.save_file(df=counts_df, outpath=os.path.join(self.file_outdir, "geneCounts.txt.gz"))
        self.save_file(df=stats_df, outpath=os.path.join(self.file_outdir, "geneCounts_stats.txt.gz"))
        self.save_file(df=translate_df, outpath=os.path.join(self.file_outdir, "translate.txt.gz"), index=None)
        del stats_df, translate_df

        print("Loading data")
        # counts_df = self.load_file(os.path.join(self.file_outdir, "geneCounts.txt.gz"), header=0, index_col=0)
        # translate_df = self.load_file(os.path.join(self.file_outdir, "translate.txt.gz"), header=0, index_col=None)
        # std = dict(zip(translate_df.loc[:, "Revised ID"], translate_df.loc[:, "Dataset"]))

        print("Step 1: remove probes with zero variance.")
        mask = counts_df.std(axis=1) != 0
        print("\tUsing {}/{} probes.".format(np.sum(mask), np.size(mask)))
        counts_df = counts_df.loc[mask, :]

        print("Step 2: remove samples with zero counts / variance.")
        mask = (counts_df.std(axis=0) != 0) & (counts_df.sum(axis=0) != 0)
        print("\tUsing {}/{} samples.".format(np.sum(mask), np.size(mask)))
        counts_df = counts_df.loc[:, mask]

        # print("Step 3: PCA analysis.")
        # self.pca(df=counts_df,
        #          sample_to_dataset=std,
        #          plot_appendix="_geneCounts")

        print("Step 4: Loading gene length")
        gene_info_df = self.load_file(self.gene_info_path, header=0, index_col=0)
        print(gene_info_df)
        missing_genes = [gene for gene in counts_df.index if gene not in gene_info_df.index]
        print(len(missing_genes))
        print(missing_genes)

        # Subset overlap.
        overlap = set(counts_df.index).intersection(set(gene_info_df.index))
        counts_df = counts_df.loc[overlap, :]

        for length_col in ["TotalGeneLength", "MergedExonLength"]:
            tmp_counts_df = counts_df.copy()

            print("Step 5: Calculating TPM values")
            # https://btep.ccr.cancer.gov/question/faq/what-is-the-difference-between-rpkm-fpkm-and-tpm/
            # Divide the read counts by the length of each gene in kilobases. This
            # gives you reads per kilobase (RPK).
            kilo_bases_s = gene_info_df.loc[overlap, length_col] / 1e3
            rpk_df = tmp_counts_df.divide(kilo_bases_s, axis=0)

            # Count up all the RPK values in a sample and divide this number by
            # 1,000,000. This is your “per million” scaling factor.
            pm_scaling_factor = rpk_df.sum(axis=0) / 1e6

            # Divide the RPK values by the “per million” scaling factor.
            # This gives you TPM.
            tpm_df = rpk_df.divide(pm_scaling_factor, axis=1)

            print("\tSaving data")
            self.save_file(df=tpm_df, outpath=os.path.join(self.file_outdir, "geneCounts.TPM.{}.txt.gz".format(length_col)))
            # print("Step 7: PCA analysis.")
            # self.pca(df=tpm_df,
            #          sample_to_dataset=std,
            #          plot_appendix="_geneCounts_TPM")

            print("Step 8: log2 transform.")
            min_value = tpm_df.min(axis=1).min()
            if min_value <= 0:
                tpm_df = np.log2(tpm_df - min_value + 1)
            else:
                tpm_df = np.log2(tpm_df + 1)

            print("\tSaving data")
            self.save_file(df=tpm_df, outpath=os.path.join(self.file_outdir, "geneCounts.TPM.{}.Log2Transformed.txt.gz".format(length_col)))

            # print("Step 9: PCA analysis.")
            # self.pca(df=tpm_df,
            #          sample_to_dataset=std,
            #          plot_appendix="_geneCounts_TPM_Log2Transformed")

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

    def pca(self, df, sample_to_dataset, plot_appendix=""):
        # samples should be on the columns and genes on the rows.
        zscores = (df - df.mean(axis=0)) / df.std(axis=0)
        pca = PCA(n_components=2)
        pca.fit(zscores)
        expl_variance = {"PC{}".format(i+1): pca.explained_variance_ratio_[i] * 100 for i in range(2)}
        components_df = pd.DataFrame(pca.components_)
        components_df.index = ["Comp{}".format(i + 1) for i, _ in enumerate(components_df.index)]
        components_df.columns = df.columns

        print("\tPlotting PCA")
        plot_df = components_df.T
        plot_df["cohort"] = plot_df.index.map(sample_to_dataset)
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
        print("  > Data path: {}".format(self.data_path))
        print("  > Gene info path: {}".format(self.data_path))
        print("  > MetaBrain samples path: {}".format(self.mb_samples_path))
        print("  > PsychENCODE samples path: {}".format(self.pe_samples_path))
        print("  > Plot output directory {}".format(self.plot_outdir))
        print("  > File output directory {}".format(self.file_outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
