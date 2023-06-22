#!/usr/bin/env python3

"""
File:         plot_df_per_column.py
Created:      2021/10/26
Last Changed: 2021/10/28
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
import os

# Third party imports.
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Plot DataFrame per Column"
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
./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2021-08-27-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2021-10-19-Cortex-EUR-mergedMDS-RNAIDs.txt -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/OLD/ContainsDuplicateSamples/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/MetaBrain-allchr-mds-dupsremoved-VariantFilter.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/MetaBrain-allchr-mds-noENA-dupsremoved.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/MetaBrain-allchr-mds-noENA-dupsremoved-outlierremoved.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/MetaBrain-allchr-mds-noENA-dupsremoved-outlierremoved-VariantFilter.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -o MetaBrain-allchr-mds-noENA-dupsremoved-outlierremoved-VariantFilter -e png pdf

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/MetaBrain-allchr-mds-dupsremoved-VariantFilter.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/combine_gte_files/SampleToDataset.txt.gz -o MetaBrain-allchr-mds-dupsremoved-VariantFilter -e png pdf

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/BIOS-allchr-mds-NoGONL-VariantFilter_genotypeID.txt.gz -std /groups/umcg-bios/tmp01/projects/PICALO/2021-10-28-DataPreprocessing/BIOS_GTESubset/1000G_harmonized_BIOS_GenotypeToDataset.txt.gz -o BIOS-allchr-mds-NoGONL-VariantFilter -e png pdf

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/BIOS-allchr-mds-BIOS-GTESubset-VariantFilter_genotypeID.txt.gz -std /groups/umcg-bios/tmp01/projects/PICALO/2021-10-28-DataPreprocessing/BIOS_GTESubset/1000G_harmonized_BIOS_GenotypeToDataset.txt.gz -o BIOS-allchr-mds-BIOS-GTESubset-VariantFilter -e png pdf

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/BIOS-allchr-mds-BIOS-GTESubset-VariantSubsetFilter_genotypeID.txt.gz -std /groups/umcg-bios/tmp01/projects/PICALO/2021-10-28-DataPreprocessing/BIOS_GTESubset/1000G_harmonized_BIOS_GenotypeToDataset.txt.gz -o BIOS-allchr-mds-BIOS-GTESubset-VariantSubsetFilter -e png pdf

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/BIOS-allchr-mds-BIOS-GTESubset-noOutlier-VariantSubsetFilter_genotypeID.txt.gz -std /groups/umcg-bios/tmp01/projects/PICALO/2021-10-28-DataPreprocessing/BIOS_GTESubset/1000G_harmonized_BIOS_GenotypeToDataset.txt.gz -o BIOS-allchr-mds-BIOS-GTESubset-noOutlier-VariantSubsetFilter -e png pdf

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/BIOS-allchr-mds-BIOS-GTESubset-noRNAPhenoNA-VariantSubsetFilter_genotypeID.txt.gz -std /groups/umcg-bios/tmp01/projects/PICALO/2021-10-28-DataPreprocessing/BIOS_GTESubset/1000G_harmonized_BIOS_GenotypeToDataset.txt.gz -o BIOS-allchr-mds-BIOS-GTESubset-noRNAPhenoNA-VariantSubsetFilter -e png pdf

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/preprocess_mds_file/BIOS-allchr-mds-BIOS-GTESubset-noRNAPhenoNA-noOutliers-VariantSubsetFilter_genotypeID.txt.gz -std /groups/umcg-bios/tmp01/projects/PICALO/2021-10-28-DataPreprocessing/BIOS_GTESubset/1000G_harmonized_BIOS_GenotypeToDataset.txt.gz -o BIOS-allchr-mds-BIOS-GTESubset-noRNAPhenoNA-noOutliers-VariantSubsetFilter -e png pdf

./plot_df_per_column.py -d /groups/umcg-bios/tmp01/projects/PICALO/preprocess_scripts/preprocess_mds_file/BIOS-allchr-mds-BIOS-GTESubset-noRNAPhenoNA-VariantSubsetFilter_genotypeID.txt.gz -std /groups/umcg-bios/tmp01/projects/PICALO/data/BIOS_STD.txt.gz -o test -e png pdf

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/preprocess_scripts/pre_process_metabrain_expression_matrix/MetaBrain-CortexEUR-cis/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.PCAOverSamplesEigenvectors.txt.gz -std /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/preprocess_scripts/pre_process_metabrain_expression_matrix/MetaBrain-CortexEUR-cis/data/SampleToDataset.txt.gz -transpose -n 4 -o MetaBrain-CortexEUR-cis -e png pdf 

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/preprocess_scripts/pre_process_metabrain_expression_matrix/MetaBrain-CortexEUR-cis/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.PCAOverSamplesEigenvectors.txt.gz -std /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/preprocess_scripts/pre_process_metabrain_expression_matrix/MetaBrain-CortexEUR-cis/data/SampleToDataset.txt.gz -transpose -n 4 -o MetaBrain-CortexEUR-cis-CovariatesRemovedOLS -e png pdf 

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/preprocess_scripts/pre_process_metabrain_expression_matrix/MetaBrain-CortexEUR-cis-NoENA-NoMDSOutlier/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.PCAOverSamplesEigenvectors.txt.gz -std /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/preprocess_scripts/pre_process_metabrain_expression_matrix/MetaBrain-CortexEUR-cis-NoENA-NoMDSOutlier/data/SampleToDataset.txt.gz -transpose -n 4 -o MetaBrain-CortexEUR-cis-NoENA-NoMDSOutlier -e png pdf 

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/preprocess_scripts/pre_process_metabrain_expression_matrix/MetaBrain-CortexEUR-cis-NoENA-NoMDSOutlier/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.PCAOverSamplesEigenvectors.txt.gz -std /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/preprocess_scripts/pre_process_metabrain_expression_matrix/MetaBrain-CortexEUR-cis-NoENA-NoMDSOutlier/data/SampleToDataset.txt.gz -transpose -n 4 -o MetaBrain-CortexEUR-cis-NoENA-NoMDSOutlierCovariatesRemovedOLS -e png pdf 

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/output/MetaBrain-CortexEUR-cis-Uncorrected-NoENA-NoMDSOutlier-MAF5/components.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis-Uncorrected-NoENA-NoMDSOutlier/combine_gte_files/SampleToDataset.txt.gz -transpose -o MetaBrain-CortexEUR-cis-NoENA-NoMDSOutlier-MAF5-PICs -e png pdf 

./plot_df_per_column.py -d /groups/umcg-biogen/tmp01/output/2020-11-10-PICALO/output/MetaBrain-CortexEUR-cis-Uncorrected-NoENA-NoMDSOutlier-MAF10/components.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis-Uncorrected-NoENA-NoMDSOutlier/combine_gte_files/SampleToDataset.txt.gz -transpose -o MetaBrain-CortexEUR-cis-NoENA-NoMDSOutlier-MAF10-PICs -e png pdf 
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.n_columns = getattr(arguments, 'n_columns')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.transpose = getattr(arguments, 'transpose')
        self.extensions = getattr(arguments, 'extension')
        self.output_filename = getattr(arguments, 'output')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
           "GTE-EUR-AMPAD-MAYO-V2": "#9c9fa0",
           "GTE-EUR-CMC_HBCC_set2": "#0877b4",
           "GTE-EUR-GTEx": "#0fa67d",
           "GTE-EUR-AMPAD-ROSMAP-V2": "#6950a1",
           "GTE-EUR-BrainGVEX-V2": "#48b2e5",
           "GTE-EUR-TargetALS": "#d5c77a",
           "GTE-EUR-AMPAD-MSBB-V2": "#5cc5bf",
           "GTE-EUR-NABEC-H610": "#6d743a",
           "GTE-EUR-LIBD_1M": "#e49d26",
           "GTE-EUR-ENA": "#d46727",
           "GTE-EUR-LIBD_h650": "#e49d26",
           "GTE-EUR-GVEX": "#000000",
           "GTE-EUR-NABEC-H550": "#6d743a",
           "GTE-EUR-CMC_HBCC_set3": "#0877b4",
           "GTE-EUR-UCLA_ASD": "#f36d2a",
           "GTE-EUR-CMC": "#eae453",
           "GTE-EUR-CMC_HBCC_set1": "#eae453",
            "LL": "#9c9fa0",
            "PAN": "#0877b4",
            "LLS_OmniExpr": "#0fa67d",
            "NTR_AFFY": "#6950a1",
            "NTR_GONL": "#48b2e5",
            "CODAM": "#6d743a",
            "RS": "#d46727",
            "LLS_660Q": "#000000",
            "GONL": "#eae453"
        }

        self.dataset_to_cohort = {
            "GTE-EUR-AMPAD-MAYO-V2": "MAYO",
            "GTE-EUR-CMC_HBCC_set2": "CMC HBCC",
            "GTE-EUR-GTEx": "GTEx",
            "GTE-EUR-AMPAD-ROSMAP-V2": "ROSMAP",
            "GTE-EUR-BrainGVEX-V2": "Brain GVEx",
            "GTE-EUR-TargetALS": "Target ALS",
            "GTE-EUR-AMPAD-MSBB-V2": "MSBB",
            "GTE-EUR-NABEC-H610": "NABEC",
            "GTE-EUR-LIBD_1M": "LIBD",
            "GTE-EUR-ENA": "ENA",
            "GTE-EUR-LIBD_h650": "LIBD",
            "GTE-EUR-GVEX": "GVEX",
            "GTE-EUR-NABEC-H550": "NABEC",
            "GTE-EUR-CMC_HBCC_set3": "CMC HBCC",
            "GTE-EUR-UCLA_ASD": "UCLA ASD",
            "GTE-EUR-CMC": "CMC",
            "GTE-EUR-CMC_HBCC_set1": "CMC HBCC"
        }

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

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
                            help="The path to the input data.")
        parser.add_argument("-n",
                            "--n_columns",
                            type=int,
                            default=None,
                            help="The number of columns to plot.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-transpose",
                            action='store_true',
                            help="Combine the created files with force."
                                 " Default: False.")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")
        parser.add_argument("-o",
                            "--output",
                            type=str,
                            default="PlotPerColumn_ColorByCohort",
                            help="The name of the output file.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")
        df = self.load_file(self.data_path, header=0, index_col=0)
        if self.transpose:
            df = df.T
        if self.n_columns is not None:
            df = df.iloc[:, :self.n_columns]
        columns = list(df.columns)

        print("Loading sample to dataset")
        std_df = self.load_file(self.std_path, header=0, index_col=None)
        std_dict = dict(zip(std_df.iloc[:, 0], std_df.iloc[:, 1]))

        print("\tAdding color.")
        overlap = set(df.index.values).intersection(set(std_df.iloc[:, 0].values))
        if len(overlap) != df.shape[0]:
            print("Error, some samples do not have a dataset.")
            exit()
        df["dataset"] = df.index.map(std_dict)

        print("\tPlotting")
        self.plot(df=df,
                  columns=columns,
                  hue="dataset",
                  palette=self.palette,
                  name=self.output_filename)

        print("\tAdding z-score color")
        for name in columns:
            df["{} z-score".format(name)] = (df[name] - df[name].mean()) / df[name].std()
        print(df)
        df["outlier"] = "False"
        sd = 3
        df.loc[(df["{} z-score".format(columns[0])].abs() > sd) | (df["{} z-score".format(columns[1])].abs() > sd) | (df["{} z-score".format(columns[2])].abs() > sd) | (df["{} z-score".format(columns[3])].abs() > sd), "outlier"] = "True"
        print(df)
        outlier_df = df.loc[df["outlier"] == "True", :].copy()
        print(outlier_df)
        print(df.loc[df["outlier"] == "True", "dataset"].value_counts())

        outlier_df.reset_index(inplace=True, drop=False)
        outlier_df = outlier_df.loc[:, ["-"]]
        print(outlier_df)
        outlier_df.to_csv(os.path.join(self.outdir, self.output_filename + "_outliers.txt.gz"), compression="gzip", sep="\t", header=False, index=False)
        exit()

        df["outlier"] = "False"
        df.loc[['1_LLDeep_0096', '1_LLDeep_0987', '1_LLDeep_0978', '1_LLDeep_0306', '1_LLDeep_1192', '1_LLDeep_0173', '1_LLDeep_0665', '1_LLDeep_1206', '1_LLDeep_1331', '1_LLDeep_0200', '1_LLDeep_1142'], "outlier"] = "True"

        print("\tPlotting")
        self.plot(df=df,
                  columns=columns,
                  hue="outlier",
                  palette={"True": "#b22222", "False": "#000000"},
                  name=self.output_filename + "_Outlier")


        #
        #
        # pheno_df = self.load_file("../preprocess_scripts/preprocess_ampad_phenotype_file/ampad_phenotype_df.txt.gz", sep="\t", header=0, index_col=1)
        # print(pheno_df)
        # pheno_df = pheno_df.loc[:, [x for x in pheno_df.columns if x not in ["genotype_id", "dataset"]]]
        # pheno_df.fillna("-", inplace=True)
        #
        # df = df.merge(pheno_df, left_index=True, right_index=True)
        #
        # colors = ["#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#0072B2",
        #           "#D55E00", "#009E73", "#E8E8E8", "#D3D3D3", "#A9A9A9",
        #           "#808080"]
        # for hue in pheno_df.columns:
        #     hue_groups = list(df[hue].unique())
        #
        #     if "-" in hue_groups and len(hue_groups) <= 2:
        #         continue
        #     if len(hue_groups) > len(colors):
        #         continue
        #
        #     palette = {"-": "#000000"}
        #     for i, hue_group in enumerate(hue_groups):
        #         if hue_group == "-":
        #             i -= 1
        #             continue
        #
        #         palette[hue_group] = colors[i % len(colors)]
        #     print(palette)
        #
        #     self.plot(df=df,
        #               columns=columns,
        #               hue=hue,
        #               palette=palette,
        #               name=os.path.basename(self.data_path).split(".")[0] + "_" + hue,
        #               title=hue)


    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def plot(self, df, columns, hue, palette, name, title=""):
        ncols = len(columns)
        nrows = len(columns)

        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='col',
                                 sharey='row',
                                 figsize=(10 * ncols, 10 * nrows))
        sns.set(color_codes=True)
        sns.set_style("ticks")

        for i, y_col in enumerate(columns):
            for j, x_col in enumerate(columns):
                print(i, j)
                ax = axes[i, j]
                if i == 0 and j == (ncols - 1):
                    ax.set_axis_off()
                    if hue is not None and palette is not None:
                        groups_present = df[hue].unique()
                        handles = []
                        added_handles = []
                        for key, value in palette.items():
                            if key in groups_present:
                                label = key
                                if key in self.dataset_to_cohort:
                                    label = self.dataset_to_cohort[key]
                                if value + label not in added_handles:
                                    handles.append(mpatches.Patch(color=value, label=label))
                                    added_handles.append(value + label)
                        ax.legend(handles=handles, loc=4, fontsize=25)

                elif i < j:
                    ax.set_axis_off()
                    continue
                elif i == j:
                    ax.set_axis_off()

                    ax.annotate(y_col,
                                xy=(0.5, 0.5),
                                ha='center',
                                xycoords=ax.transAxes,
                                color="#000000",
                                fontsize=40,
                                fontweight='bold')
                else:
                    sns.despine(fig=fig, ax=ax)

                    sns.scatterplot(x=x_col,
                                    y=y_col,
                                    hue=hue,
                                    data=df,
                                    s=100,
                                    palette=palette,
                                    linewidth=0,
                                    legend=False,
                                    ax=ax)

                    ax.set_ylabel("",
                                  fontsize=20,
                                  fontweight='bold')
                    ax.set_xlabel("",
                                  fontsize=20,
                                  fontweight='bold')

        fig.suptitle(title,
                     fontsize=40,
                     fontweight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(name, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Data path: {}".format(self.data_path))
        print("  > N-columns: {}".format(self.n_columns))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > Transpose: {}".format(self.transpose))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output filename: {}".format(self.output_filename))
        print("  > Output directory {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
