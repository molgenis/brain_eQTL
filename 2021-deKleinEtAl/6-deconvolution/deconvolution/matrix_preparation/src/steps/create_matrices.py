"""
File:         create_matrices.py
Created:      2020/10/08
Last Changed: 2021/12/07
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
import gzip
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe, construct_dict_from_df


class CreateMatrices:
    def __init__(self, settings, log, sample_dict, sample_order, eqtl_file,
                 eqtl_df, force, outdir):
        self.geno_file = settings["genotype_datafile"]
        self.expr_file = settings["expression_datafile"]
        self.sign_file = settings["signature_datafile"]
        gene_translate_settings = settings["gene_translate"]
        self.gene_info_file = gene_translate_settings["datafile"]
        self.ensg_id = gene_translate_settings["ENSG_column"]
        self.hgnc_id = gene_translate_settings["HGNC_column"]
        self.decon_expr_file = settings["deconvolution_expr_datafile"]
        self.log = log
        self.sample_dict = sample_dict
        self.sample_order = sample_order
        self.eqtl_file = eqtl_file
        self.eqtl_df = eqtl_df
        self.force = force
        self.print_interval = 500

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_matrices')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.geno_outpath = os.path.join(self.outdir, "genotype_table.txt.gz")
        self.alleles_outpath = os.path.join(self.outdir, "genotype_alleles.txt.gz")
        self.expr_outpath = os.path.join(self.outdir, "expression_table.txt.gz")
        self.sign_expr_outpath = os.path.join(self.outdir, "signature_expr_table.txt.gz")

        # Create empty variable.
        self.alleles_df = None
        self.geno_df = None
        self.sign_df = None
        self.gene_info_df = None
        self.expr_df = None
        self.sign_expr_df = None

    def start(self):
        self.log.info("Starting creating matrices.")
        self.print_arguments()

        if self.eqtl_df is None:
            self.eqtl_df = load_dataframe(self.eqtl_file,
                                          header=0,
                                          index_col=None,
                                          logger=self.log)

        self.log.info("Parsing genotype input data.")
        if not check_file_exists(self.geno_outpath) or not check_file_exists(self.alleles_outpath) or self.force:
            alleles_df, geno_df = self.parse_genotype_file()

            self.log.info("Reorder, Filter, and save.")
            self.alleles_df = alleles_df.loc[self.eqtl_df.loc[:, "SNPName"], :]
            save_dataframe(df=self.alleles_df, outpath=self.alleles_outpath,
                           index=True, header=True, logger=self.log)

            self.geno_df = geno_df.loc[self.eqtl_df.loc[:, "SNPName"], self.sample_order]
            save_dataframe(df=self.geno_df, outpath=self.geno_outpath,
                           index=True, header=True, logger=self.log)
        else:
            self.log.info("\tSkipping step.")

        self.log.info("Parsing expression input data.")
        if not check_file_exists(self.expr_outpath) or not check_file_exists(self.sign_expr_outpath) or self.force:
            self.log.info("Loading signature matrix.")
            self.sign_df = load_dataframe(inpath=self.sign_file,
                                          header=0,
                                          index_col=0,
                                          logger=self.log)
            signature_genes = set(self.sign_df.index.to_list())

            self.log.info("Loading gene traslate dict.")
            self.gene_info_df = load_dataframe(inpath=self.gene_info_file,
                                               header=0,
                                               index_col=None,
                                               logger=self.log)
            gene_trans_dict = construct_dict_from_df(self.gene_info_df, self.ensg_id, self.hgnc_id)

            if not check_file_exists(self.expr_outpath) or self.force:
                self.log.info("Parsing expression data.")
                self.expr_df, self.sign_expr_df = self.parse_expression_file(self.expr_file,
                                                                             signature_genes,
                                                                             gene_trans_dict,
                                                                             include_decon=self.decon_expr_file is None)

            if (not check_file_exists(self.sign_expr_outpath) or self.force) and (check_file_exists(self.decon_expr_file)):
                self.log.info("Parsing deconvolution expression data.")
                self.log.warning("Using different expresion file for deconvolution.")
                _, self.sign_expr_df = self.parse_expression_file(self.decon_expr_file,
                                                                  signature_genes,
                                                                  gene_trans_dict,
                                                                  include_expr=False,
                                                                  remove_ens_version=True)

            self.log.info("Reorder, Filter, and save.")
            if self.expr_df is not None:
                self.expr_df = self.expr_df.loc[self.eqtl_df.loc[:, "ProbeName"], self.sample_order]
                save_dataframe(df=self.expr_df, outpath=self.expr_outpath,
                               index=True, header=True, logger=self.log)
            if self.sign_expr_df is not None:
                self.sign_expr_df = self.sign_expr_df.loc[:, self.sample_order]
                save_dataframe(df=self.sign_expr_df, outpath=self.sign_expr_outpath,
                               index=True, header=True, logger=self.log)
        else:
            self.log.info("\tSkipping step.")

    def parse_genotype_file(self):
        interest = set(self.eqtl_df.loc[:, "SNPName"].to_list())

        alleles_data_collection = []
        genotype_data_collection = []
        indices = []
        alleles_columns = []
        genotype_columns = []
        with gzip.open(self.geno_file, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % self.print_interval == 0):
                    self.log.info("\tprocessed {} lines\tfound {}/{} genotype lines.".format(i, len(indices), len(interest)))

                splitted_line = np.array(line.decode().strip('\n').split('\t'))
                index = splitted_line[0]
                alleles_data = splitted_line[1:3]
                genotype_data = splitted_line[3:]
                if i == 0:
                    alleles_columns = alleles_data
                    genotype_columns = [self.sample_dict[x] if x in self.sample_dict else x for x in genotype_data]
                else:
                    if index in interest and index not in indices:
                        indices.append(index)
                        alleles_data_collection.append([x for x in alleles_data])
                        genotype_data_collection.append([float(x) for x in genotype_data])

        f.close()
        self.log.info("\tprocessed all lines\tfound {}/{} genotype lines.".format(
            len(indices), len(interest)))

        alleles_df = pd.DataFrame(alleles_data_collection,
                                  index=indices,
                                  columns=alleles_columns)

        genotype_df = pd.DataFrame(genotype_data_collection,
                                   index=indices,
                                   columns=genotype_columns)

        # Add missing data.
        missing = []
        for snp_name in interest:
            if snp_name not in alleles_df.index:
                missing.append(snp_name)
                alleles_df.loc[snp_name, :] = np.nan
            if snp_name not in genotype_df.index:
                genotype_df.loc[snp_name, :] = np.nan
        self.log.warning("\tMissing SNP's [{}]: {}".format(len(missing), ", ".join(missing)))

        return alleles_df, genotype_df

    def parse_expression_file(self, filepath, signature_genes, gene_trans_dict,
                              include_expr=True, include_decon=True, remove_ens_version=False):
        if not include_expr and not include_decon:
            return None, None

        expression_interest = set(self.eqtl_df.loc[:, "ProbeName"].to_list())

        expression_data_collection = []
        expression_indices = []
        sign_expr_data_collection = []
        sign_expr_indices = []
        with gzip.open(filepath, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % self.print_interval == 0):
                    process_str = "\tprocessed {} lines".format(i)
                    if include_expr:
                          process_str += "\tfound {}/{} expression lines".format(len(expression_indices), len(expression_interest))
                    if include_decon:
                        process_str += "\tfound {}/{} signature genes".format(len(sign_expr_indices), len(signature_genes))
                    self.log.info(process_str)

                splitted_line = np.array(line.decode().strip('\n').split('\t'))
                index = splitted_line[0]
                if remove_ens_version:
                    index = index.split(".")[0]
                data = splitted_line[1:]
                if i == 0:
                    columns = [self.sample_dict[x] if x in self.sample_dict else x for x in data]
                else:
                    # hgnc_symbol = None
                    # if index in gene_trans_dict.keys():
                    #     hgnc_symbol = gene_trans_dict[index]

                    if include_expr:
                        if index in expression_interest and index not in expression_indices:
                            expression_indices.append(index)
                            expression_data_collection.append([float(x) for x in data])
                    if include_decon:
                        # if hgnc_symbol in signature_genes and hgnc_symbol not in sign_expr_indices:
                        #     sign_expr_indices.append(hgnc_symbol)
                        #     sign_expr_data_collection.append([float(x) for x in data])
                        if index in signature_genes and index not in sign_expr_indices:
                            sign_expr_indices.append(index)
                            sign_expr_data_collection.append([float(x) for x in data])
        f.close()
        process_str = "\tprocessed all lines"
        if include_expr:
            process_str += "\tfound {}/{} expression lines".format(
                len(expression_indices), len(expression_interest))
        if include_decon:
            process_str += "\tfound {}/{} signature genes".format(
                len(sign_expr_indices), len(signature_genes))
        self.log.info(process_str)

        expression_df = pd.DataFrame(expression_data_collection,
                                     index=expression_indices,
                                     columns=columns)
        sign_expr_df = pd.DataFrame(sign_expr_data_collection,
                                    index=sign_expr_indices,
                                    columns=columns)

        # Add missing data.
        if include_expr:
            missing = []
            for ensg in expression_interest:
                if ensg not in expression_df.index:
                    missing.append(ensg)
                    expression_df.loc[ensg, :] = np.nan
            self.log.warning("\tExpression missing ENSG ID's [{}]: {}".format(len(missing), ", ".join(missing)))
        if include_decon:
            missing = []
            for hgnc in signature_genes:
                if hgnc not in sign_expr_df.index:
                    missing.append(hgnc)
                    sign_expr_df.loc[hgnc, :] = np.nan
            self.log.warning("\tSignature expression missing HGNC symbols [{}]: {}".format(len(missing), ", ".join(missing)))

        return expression_df, sign_expr_df

    def clear_variables(self):
        self.geno_file = None
        self.expr_file = None
        self.gene_info_file = None
        self.ensg_id = None
        self.hgnc_id = None
        self.force = None

    def get_alleles_df(self):
        return self.alleles_df

    def get_geno_df(self):
        return self.geno_df

    def get_sign_file(self):
        return self.sign_file

    def get_sign_df(self):
        return self.sign_df

    def get_sign_expr_outpath(self):
        return self.sign_expr_outpath

    def get_gene_info_df(self):
        return self.gene_info_df

    def get_expr_file(self):
        return self.expr_outpath

    def get_expr_df(self):
        return self.expr_df

    def get_sign_expr_file(self):
        return self.sign_expr_outpath

    def get_sign_expr_df(self):
        return self.sign_expr_df

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Genotype input file: {}".format(self.geno_file))
        self.log.info("  > Expression input file: {}".format(self.expr_file))
        self.log.info("  > Signature input file: {}".format(self.sign_file))
        self.log.info("  > Gene info input file: {}".format(self.gene_info_file))
        self.log.info("  > Gene info - ENSG column: {}".format(self.ensg_id))
        self.log.info("  > Gene info - HGNC column: {}".format(self.hgnc_id))
        self.log.info("  > Deconvolution expr input file: {}".format(self.decon_expr_file))
        if self.eqtl_df is not None:
            self.log.info("  > eQTL input shape: {}".format(self.eqtl_df.shape))
        else:
            self.log.info("  > eQTL input file: {}".format(self.eqtl_file))
        self.log.info("  > Output directory: {}".format(self.outdir))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
