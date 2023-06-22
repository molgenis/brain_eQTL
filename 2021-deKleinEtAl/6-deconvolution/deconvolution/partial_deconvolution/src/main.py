"""
File:         main.py
Created:      2020/06/29
Last Changed: 2021/08/04
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
from .data_loader import DataLoader
from .data_filter import DataFilter
from .data_preprocessor import DataPreprocessor
from .perform_deconvolution import PerformDeconvolution
from .data_comparitor import DataComparitor
from .visualiser import Visualiser


class Main:
    def __init__(self, settings, outdir, outsubdir, visualise, plot_ids):
        self.settings = settings
        self.visualise = visualise
        self.plot_ids = plot_ids

        # Prepare output directories.
        current_dir = str(Path(__file__).parent.parent)
        outdir_path = os.path.join(current_dir, outdir)
        outsubdir_path = outdir_path
        if outsubdir is not None:
            outsubdir_path = os.path.join(outdir_path, outsubdir)
        for path in [outdir_path, outsubdir_path]:
            if not os.path.exists(path):
                os.makedirs(path)

        self.settings.set_outdir_path(outdir_path)
        self.settings.set_outsubdir_path(outsubdir_path)

    def start(self):
        # Load the data.
        print("### Loading")
        dl = DataLoader(settings=self.settings)
        dl.work()
        dl.print_info()
        self.settings.set_real_info_per_celltype(dl.get_info_per_celltype())

        # Save.
        self.settings.save_data_settings()

        # Filter the samples.
        print("### Filtering")
        filtered_expression = dl.get_expression()
        datasets = None
        if self.settings.sample_to_dataset_path is not None:
            df = DataFilter(settings=self.settings,
                            raw_expression=filtered_expression)
            df.work()
            df.print_info()
            self.settings.set_filter1_shape_diff(df.get_shape_diff())
            self.settings.set_reference_dataset(df.get_reference_dataset())

            filtered_expression = df.get_filtered_expression()
            datasets = df.get_datasets()

        # Preprocessing.
        print("### Preprocessing")
        dp = DataPreprocessor(settings=self.settings,
                              raw_signature=dl.get_signature(),
                              raw_expression=filtered_expression,
                              datasets=datasets)
        dp.work()
        dp.print_info()

        self.settings.set_filter2_shape_diff(dp.get_shape_diff())
        self.settings.set_sign_shift(dp.get_sign_shift())
        self.settings.set_expr_shift(dp.get_expr_shift())
        self.settings.set_n_samples(dp.get_n_samples())
        self.settings.set_n_genes(dp.get_n_genes())
        self.settings.set_n_ng_per_ct(dp.get_n_mg_per_ct())

        # Partial deconvolution.
        print("### Deconvoluting")
        pf = PerformDeconvolution(settings=self.settings,
                                  signature=dp.get_signature(),
                                  expression=dp.get_expression())
        pf.work()
        pf.print_info()
        self.settings.set_avg_residuals(pf.get_avg_rss())
        self.settings.set_avg_recon_accuracy(pf.get_avg_recon_accuracy())
        self.settings.set_pred_info_per_celltype(pf.get_info_per_celltype())

        # Comparison.
        print("### Comparing")
        dc = DataComparitor(settings=self.settings,
                            deconvolution=pf.get_deconvolution(),
                            ground_truth=dl.get_ground_truth())
        dc.work()
        dc.print_info()
        self.settings.set_comparison_n_samples(dc.get_n_samples())
        self.settings.set_comparison_rss(dc.get_rss())

        # Save.
        self.settings.save_all_settings()

        # Visualising profile.
        if self.visualise:
            print("### Visualising")
            v = Visualiser(settings=self.settings,
                           signature=dp.get_signature(),
                           expression=dp.get_expression(),
                           deconvolution=pf.get_deconvolution(),
                           ground_truth=dl.get_ground_truth(),
                           comparison=dc.get_comparison())
            v.plot_profile_clustermap()
            v.plot_profile_correlations()
            v.plot_profile_stripplot()
            v.plot_profile_boxplot()
            v.plot_deconvolution_clustermap()
            v.plot_deconvolution_correlations()
            v.plot_deconvolution_per_sample()
            #v.plot_deconvolution_distribution()
            v.plot_deconvolution_boxplot()
            #v.plot_ground_truth_distribution()
            v.plot_ground_truth_boxplot()
            v.plot_prediction_comparison()
            v.plot_recon_accuracy_boxplot(s=pf.get_recon_accuracy())

            if self.plot_ids is not None:
                for plot_id in self.plot_ids:
                    print("Plotting {}".format(plot_id))
                    v.plot_violin_comparison(plot_id, df.get_sample_to_dataset_dict())