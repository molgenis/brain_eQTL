"""
File:         main.py
Created:      2020/06/08
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
from __future__ import print_function
from pathlib import Path
import os

# Third party imports.
import pandas as pd

# Local application imports.
from general.utilities import prepare_output_dir
from general.local_settings import LocalSettings
from general.objects.dataset import Dataset
from .eqtl import Eqtl
from .plotter import Plotter
from .saver import Saver


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, name, settings_file, alpha, extensions, interest):
        """
        Initializer of the class.

        :param name: string, the name of the base input/ouput directory.
        :param settings_file: string, the name of the settings file.
        :param alpha: float, the significance cut-off.
        :param extensions: str, the output figure file type extension.
        :param interest: list, the HGNC names to print.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        self.settings = LocalSettings(current_dir, settings_file)
        self.covs = self.settings.get_setting("covariates_to_include")
        self.covs_excl_from_overview = [x.lower() for x in self.settings.get_setting("covariates_excl_from_overview")]
        self.max_url_len = self.settings.get_setting("max_url_length")
        self.maf_cutoff = self.settings.get_setting("maf_cutoff")
        self.include_top_n = self.settings.get_setting("include_top_n")

        # Load the variables.
        self.name = name
        self.alpha = alpha
        self.extensions = extensions
        self.interest = interest

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir, name)
        prepare_output_dir(self.outdir)

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting visualiser.")
        self.print_arguments()

        # Create the dataset object.
        ds = Dataset(name=self.name,
                     settings=self.settings,
                     alpha=self.alpha)

        # Loading data.
        print("Loading dataframes")
        cov_df = ds.get_cov_df()
        inter_zscore_df = ds.get_inter_cov_zscore_df()
        inter_tvalue_df = ds.get_inter_cov_inter_tvalue_df()
        eqtl_df = ds.get_eqtl_df()
        geno_df = ds.get_geno_df()
        alleles_df = ds.get_alleles_df()
        expr_df = ds.get_expr_df()
        signif_cutoff = ds.get_significance_cutoff()

        # Subset the method of interest.
        select_cov_df = cov_df.loc[self.covs.keys(), :].copy()
        select_zscore_df = inter_zscore_df.loc[self.covs.keys(), :].copy().T
        select_tvalue_df = inter_tvalue_df.loc[self.covs.keys(), :].copy().T
        del cov_df, inter_zscore_df, inter_tvalue_df

        if not select_cov_df.index.equals(select_zscore_df.columns) or not select_cov_df.index.equals(select_tvalue_df.columns):
            print("Columns do not match.")
            exit()

        print("Analysing interactions")
        data = []
        for i, (index, row) in enumerate(eqtl_df.iterrows()):
            if (i % 250 == 0) or (i == (eqtl_df.shape[0] - 1)):
                print("\tprocessing {}/{} "
                      "[{:.2f}%]".format(i,
                                         (eqtl_df.shape[0] - 1),
                                         (100 / (eqtl_df.shape[0] - 1)) * i))

            # Get the data.
            genotype = geno_df.iloc[i, :].copy()
            expression = expr_df.iloc[i, :].copy()
            (alleles, _) = alleles_df.iloc[i, :].copy()
            zscores = select_zscore_df.iloc[i, :].copy()
            tvalues = select_tvalue_df.iloc[i, :].copy()

            iteration = None
            if "Iteration" in row.index:
                iteration = row["Iteration"]
            gwas_ids = None
            if "GWASIDS" in row.index:
                gwas_ids = row["GWASIDS"]
            traits = None
            if "Trait" in row.index:
                traits = row["Trait"]

            if max(zscores) > signif_cutoff:
                eqtl = Eqtl(index=i,
                            snp_name=row["SNPName"],
                            probe_name=row["ProbeName"],
                            hgnc_name=row["HGNCName"],
                            iteration=iteration,
                            eqtl_zscore=row["OverallZScore"],
                            gwas_ids=gwas_ids,
                            traits=traits,
                            alleles=alleles,
                            signif_cutoff=signif_cutoff,
                            maf_cutoff=self.maf_cutoff,
                            selections=self.covs,
                            genotype=genotype,
                            expression=expression,
                            covariates=select_cov_df.copy(),
                            inter_zscores=zscores,
                            inter_tvalues=tvalues)
                if self.interest is not None and row["HGNCName"] in self.interest:
                    eqtl.print_info()
                data.extend(eqtl.get_data())
                del eqtl

        # Create the complete dataframe.
        data_df = pd.DataFrame(data, columns=["Index", "SNPName", "ProbeName",
                                              "HGNCName", "Iteration", "N",
                                              "MAF", "eQTL", "Inter",
                                              "Covariate", "Interaction",
                                              "Direction", "GWASIDs", "Traits"])

        # Plot the data.
        print("Creating plots")
        for extension in self.extensions:
            plotter = Plotter(data_df, self.outdir, extension=extension)
            plotter.plot_upsetplot(column="Covariate", id_col="Index", exclude=self.covs_excl_from_overview)
            plotter.plot_upsetplot(column="Iteration", id_col="ID")
            plotter.plot_pie(total=eqtl_df.shape[0], part=len(data_df["Index"].unique()))

        # Save data files.
        print("Saving results")
        saver = Saver(data_df, self.outdir, signif_cutoff, self.max_url_len, self.include_top_n)
        saver.save_all(exclude=self.covs_excl_from_overview)
        saver.save_per_iter(exclude=self.covs_excl_from_overview)
        indices_of_interest = saver.save_per_group()
        print("eQTL indices of interest: {}".format(' '.join([str(x) for x in indices_of_interest])))

    def print_arguments(self):
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("  > Alpha: {}".format(self.alpha))
        print("  > Extensions: {}".format(self.extensions))
        print("  > Interest: {}".format(self.interest))
        print("  > Covariates to include: {}".format(self.covs))
        print("  > Covariates to exclude from overview: {}".format(self.covs_excl_from_overview))
        print("  > Max URL length: {}".format(self.max_url_len))
        print("  > MAF cutoff: {}".format(self.maf_cutoff))
        print("  > Include top n: {}".format(self.include_top_n))
        print("")
