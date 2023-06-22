#!/usr/bin/env python3

"""
File:         gtex_cellmap_expression.py
Created:      2020/06/16
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
import gzip
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "GTEx CellMap Expression"
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
        self.sample_rna_link_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-25-samplelinks/all/links2-ExpressionSamplesLinked.txt"
        self.rna_tissue_link_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/archive/2020-03-01.brain.phenotypes.txt"
        self.expression_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"
        #self.expression_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-02-06-step4-remove-residual-covariates/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz"
        #self.expression_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/output-cortex/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz"
        self.profile_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt"
        self.translate_path = "/groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz"
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'gtex_tissue')
        self.outpath = os.path.join(self.outdir, "SamplesZTransformed.txt.gz")

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        print("Load the sample-rna link file")
        samplerna_df = pd.read_csv(self.sample_rna_link_path, sep="\t")
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.sample_rna_link_path),
                                      samplerna_df.shape))
        samplerna_df = samplerna_df.loc[samplerna_df["MetaCohort"] == "GTEx", :]
        samplerna_dict = dict(zip(samplerna_df.loc[:, "RnaID"], samplerna_df.loc[:, "GenotypeID"]))

        print("Load the rna-tissue link file")
        rnatissue_df = pd.read_csv(self.rna_tissue_link_path, sep="\t")
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.rna_tissue_link_path),
                                      rnatissue_df.shape))
        rnatissue_df = rnatissue_df.loc[rnatissue_df["MetaCohort"] == "GTEx", :]
        rnabroad_dict = dict(zip(rnatissue_df.loc[:, "rnaseq_id"], rnatissue_df.loc[:, "BroadBrainRegion"]))
        rnaspecific_dict = dict(zip(rnatissue_df.loc[:, "rnaseq_id"], rnatissue_df.loc[:, "SpecificBrainRegion"]))

        print("Loading translate matrix")
        trans_df = pd.read_csv(self.translate_path, sep="\t", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.translate_path),
                                      trans_df.shape))
        trans_dict = dict(zip(trans_df.loc[:, "ArrayAddress"], trans_df.loc[:, "Symbol"]))
        del trans_df

        print("Load the profile")
        profile_df = pd.read_csv(self.profile_path, sep="\t", header=0,
                                 index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.profile_path),
                                      profile_df.shape))
        profile_genes = list(profile_df.index)

        # Start reading the file.
        print("Filtering expression file.")
        selection = []
        data = []
        columns = []
        indices = []
        keys = []
        with gzip.open(self.expression_path, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % 1000 == 0):
                    print("\tfound {}/{} genes".format(len(data), len(profile_genes)))

                splitted_line = np.array(line.decode().split('\t'))
                if i == 0:
                    for j, col in enumerate(splitted_line):
                        if col in samplerna_dict.keys() and col in rnabroad_dict.keys() and col in rnaspecific_dict.keys():
                            key = (samplerna_dict[col], rnabroad_dict[col], rnaspecific_dict[col])
                            if key not in keys:
                                selection.append(j)
                                columns.append(key)

                gene = splitted_line[0]
                if gene in trans_dict.keys() and trans_dict[gene] in profile_genes and trans_dict[gene] not in indices:
                    indices.append(trans_dict[gene])
                    data.append([float(x) for x in splitted_line[selection]])

                # if len(data) > 25:
                #     break
        f.close()
        print("Found {}/{} genes for {} samples.".format(len(data), len(profile_genes), len(selection)))

        df = pd.DataFrame(data, index=indices, columns=pd.MultiIndex.from_tuples(columns, names=['sample', 'broad_region', 'specific_region']))
        print(df)

        # Start writing.
        if len(df.index) > 0:
            df.to_csv(self.outpath, sep="\t", compression="gzip", header=True, index=True)

if __name__ == '__main__':
    m = main()
    m.start()
