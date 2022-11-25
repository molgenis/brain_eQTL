"""
File:         main.py
Created:      2020/03/13
Last Changed: 2020/06/08
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

# Local application imports.
from general.utilities import prepare_output_dir
from general.local_settings import LocalSettings
from general.objects.dataset import Dataset
from .figures.covariate_clustermap import CovariateClustermap
from .figures.covariate_comparison import CovariateComparison
from .figures.deconvolution_covariate_comparison import DeconvolutionCovariateComparison
from .figures.covariates_explained_by_others import CovariatesExplainedByOthers
from .figures.deconvolution_zscore_comparison import DeconvolutionZscoreComparison
from .figures.simple_eqtl_effect import SimpleeQTLEffect
from .figures.inter_zscore_bars import InterZscoreBars
from .figures.inter_pvalue_boxplot import InterPvalueBoxplot
from .figures.inter_zscore_dist import InterZscoreDist
from .figures.inter_clustermap import InterClusterMap
from .figures.inter_eqtl_zscore_bars import IntereQTLZscoreBars
from .figures.inter_eqtl_effect import IntereQTLEffect
from .figures.inter_eqtl_effect_deconvolution import IntereQTLEffectDeconvolution
from .figures.inter_eqtl_effect_celltype import IntereQTLEffectCelltype
from .figures.inter_eqtl_celltype_details import IntereQTLCelltypeDetails


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, name, settings_file, alpha, plots, top, interest,
                 extension, validate):
        """
        Initializer of the class.

        :param name: string, the name of the base input/ouput directory.
        :param settings_file: string, the name of the settings file.
        :param alpha: float, the significance cut-off.
        :param plots: list, the names of the plots to create.
        :param top: int, the number of top eQTLs to plot.
        :param interest: list, the indices of equals to plot.
        :param extension: str, the output figure file type extension.
        :param validate: boolean, whether or not to validate the input.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        self.settings = LocalSettings(current_dir, settings_file)

        # Load the variables.
        self.name = name
        self.alpha = alpha
        self.plots = plots
        self.top = top
        self.interest = interest
        self.extension = extension
        self.validate = validate

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
                     alpha=self.alpha,
                     nrows=self.top,
                     interest=self.interest)
        if self.validate:
            ds.load_all()

        if ('covariate_clustermap' in self.plots) or \
                ('all' in self.plots):
            print("\n### Covariate Clustermap ###\n")
            cocl = CovariateClustermap(dataset=ds, outdir=self.outdir,
                                       extension=self.extension)
            cocl.start()
            del cocl

        if ('covariate_comparison' in self.plots) or \
                ('all' in self.plots):
            print("\n### Covariate Comparison ###\n")
            coco = CovariateComparison(dataset=ds, outdir=self.outdir,
                                       extension=self.extension)
            coco.start()
            del coco

        if ('covariates_explained_by_others' in self.plots) or \
                ('all' in self.plots):
            print("\n### Covariates Explained By Others ###\n")
            cebo = CovariatesExplainedByOthers(dataset=ds,
                                               outdir=self.outdir,
                                               extension=self.extension)
            cebo.start()
            del cebo

        if ('deconvolution_covariate_comparison' in self.plots) or \
                ('all' in self.plots):
            print("\n### Deconvolution Covariate Comparison ###\n")
            dcc = DeconvolutionCovariateComparison(dataset=ds,
                                                   outdir=self.outdir,
                                                   extension=self.extension)
            dcc.start()
            del dcc

        if ('deconvolution_zscore_comparison' in self.plots) or \
                ('all' in self.plots):
            print("\n### DECONVOLUTION Z-SCORE COMPARISON ###\n")
            dzc = DeconvolutionZscoreComparison(dataset=ds,
                                                outdir=self.outdir,
                                                extension=self.extension)
            dzc.start()
            del dzc

        if ('simple_eqtl_effect' in self.plots) or ('all' in self.plots):
            print("\n### SIMPLE EQTL EFFECT ###\n")
            sef = SimpleeQTLEffect(dataset=ds,
                                   outdir=self.outdir,
                                   extension=self.extension)
            sef.start()
            del sef

        if ('inter_clustermap' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION CLUSTERMAP ###\n")
            icp = InterClusterMap(dataset=ds,
                                  outdir=self.outdir,
                                  extension=self.extension)
            icp.start()
            del icp

        if ('inter_zscore_bars' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION Z-SCORE BARPLOT ###\n")
            izb = InterZscoreBars(dataset=ds,
                                  outdir=self.outdir,
                                  extension=self.extension)
            izb.start()
            del izb

        if ('inter_pvalue_boxplot' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION P-VALUE BOXPLOT ###\n")
            ipb = InterPvalueBoxplot(dataset=ds,
                                     outdir=self.outdir,
                                     extension=self.extension)
            ipb.start()
            del ipb

        if ('inter_zscore_dist' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION Z-SCORE DISTRIBUTION PLOT ###\n")
            izd = InterZscoreDist(dataset=ds,
                                  outdir=self.outdir,
                                  extension=self.extension)
            izd.start()
            del izd

        if ('inter_eqtl_zscore_bars' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION EQTL Z-SCORE BARS ###\n")
            iezb = IntereQTLZscoreBars(dataset=ds,
                                       outdir=self.outdir,
                                       extension=self.extension)
            iezb.start()
            del iezb

        if ('inter_eqtl_effect' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION EQTL EFFECT ###\n")
            iee = IntereQTLEffect(dataset=ds,
                                  outdir=self.outdir,
                                  extension=self.extension)
            iee.start()
            del iee

        if ('inter_eqtl_effect_deconvolution' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION EQTL EFFECT DECONVOLUTION ###\n")
            ieed = IntereQTLEffectDeconvolution(dataset=ds,
                                                outdir=self.outdir,
                                                extension=self.extension)
            ieed.start()
            del ieed

        if ('inter_eqtl_effect_celltype' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION EQTL EFFECT CELLTYPE ###\n")
            ieec = IntereQTLEffectCelltype(dataset=ds,
                                           outdir=self.outdir,
                                           extension=self.extension)
            ieec.start()
            del ieec

        if ('inter_eqtl_celltype_details' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION EQTL EFFECT CELLTYPE DETAILS ###\n")
            iecd = IntereQTLCelltypeDetails(dataset=ds,
                                            outdir=self.outdir,
                                            extension=self.extension)
            iecd.start()
            del iecd

    def print_arguments(self):
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("  > Alpha: {}".format(self.alpha))
        print("  > Plots: {}".format(self.plots))
        print("  > Top: {}".format(self.top))
        print("  > Interest: {}".format(self.interest))
        print("  > Validate: {}".format(self.validate))
        print("")
