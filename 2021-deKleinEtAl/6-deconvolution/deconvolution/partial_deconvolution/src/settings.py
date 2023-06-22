"""
File:         settings.py
Created:      2020/06/29
Last Changed: 2021/11/22
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
from datetime import datetime
import json
import os

# Third party imports.

# Local application imports.


class Settings:
    def __init__(self, data_path, signature_path, translate_path,
                 ground_truth_path, sample_to_dataset_path, sample_filter_path,
                 dataset_filter, min_expr, dataset_correction, normalize, zscore,
                 log2, decon_method, extension):
        self.data_path = data_path
        self.signature_path = signature_path
        self.translate_path = translate_path
        self.ground_truth_path = ground_truth_path
        self.sample_to_dataset_path = sample_to_dataset_path
        self.sample_filter_path = sample_filter_path
        if dataset_filter is not None:
            dataset_filter = [x.replace("_", " ").upper() for x in dataset_filter]
        self.dataset_filter = dataset_filter
        self.min_expr = min_expr
        self.dataset_correction = dataset_correction
        self.normalize = normalize
        self.zscore = zscore
        self.log2 = log2
        self.decon_method = decon_method
        self.extension = extension

        self.outdir_path = None
        self.outsubdir_path = None
        self.real_info_per_celltype = None
        self.filter1_shape_diff = None
        self.filter2_shape_diff = None
        self.reference_dataset = None
        self.sign_shift = None
        self.expr_shift = None
        self.n_samples = None
        self.n_genes = None
        self.n_ng_per_ct = None
        self.avg_residuals = None
        self.avg_recon_accuracy = None
        self.pred_info_per_celltype = None
        self.comparison_n_samples = None
        self.comparison_rss = None

    def get_data_path(self):
        return self.data_path

    def get_signature_path(self):
        return self.signature_path

    def get_translate_path(self):
        return self.translate_path

    def get_ground_truth_path(self):
        return self.ground_truth_path

    def get_ground_truth_type(self):
        return os.path.basename(self.ground_truth_path).split(".")[0].replace("_counts", "")

    def get_sample_to_dataset_path(self):
        return self.sample_to_dataset_path

    def get_sample_filter_path(self):
        return self.sample_filter_path

    def get_dataset_filter(self):
        return self.dataset_filter

    def get_min_expr(self):
        return self.min_expr

    def get_dataset_correction(self):
        return self.dataset_correction

    def get_normalize(self):
        return self.normalize

    def get_zscore(self):
        return self.zscore

    def get_log2(self):
        return self.log2

    def get_decon_method(self):
        return self.decon_method

    def get_extension(self):
        return self.extension

    def get_outdir_path(self):
        return self.outdir_path

    def get_outsubdir_path(self):
        return self.outsubdir_path

    def set_outdir_path(self, outdir_path):
        self.outdir_path = outdir_path

    def set_outsubdir_path(self, outsubdir_path):
        self.outsubdir_path = outsubdir_path

    def set_real_info_per_celltype(self, real_info_per_celltype):
        self.real_info_per_celltype = real_info_per_celltype

    def set_filter1_shape_diff(self, filter1_shape_diff):
        self.filter1_shape_diff = filter1_shape_diff

    def set_filter2_shape_diff(self, filter2_shape_diff):
        self.filter2_shape_diff = filter2_shape_diff

    def set_reference_dataset(self, reference_dataset):
        self.reference_dataset = reference_dataset

    def set_sign_shift(self, sign_shift):
        self.sign_shift = sign_shift

    def set_expr_shift(self, expr_shift):
        self.expr_shift = expr_shift

    def set_n_samples(self, n_samples):
        self.n_samples = n_samples

    def set_n_genes(self, n_genes):
        self.n_genes = n_genes

    def set_n_ng_per_ct(self, n_ng_per_ct):
        self.n_ng_per_ct = n_ng_per_ct

    def set_avg_residuals(self, avg_residuals):
        self.avg_residuals = avg_residuals

    def set_avg_recon_accuracy(self, avg_recon_accuracy):
        self.avg_recon_accuracy = avg_recon_accuracy

    def set_pred_info_per_celltype(self, pred_info_per_celltype):
        self.pred_info_per_celltype = pred_info_per_celltype

    def set_comparison_n_samples(self, comparison_n_samples):
        self.comparison_n_samples = comparison_n_samples

    def set_comparison_rss(self, comparison_rss):
        self.comparison_rss = comparison_rss

    def get_title(self):
        return "partial deconvolution using {}".format(self.decon_method)

    def get_subtitle(self):
        avg_resid_str = ""
        if self.avg_residuals is not None:
            avg_resid_str = "{:.2f}".format(self.avg_residuals)

        avg_recon_accuracy_str = ""
        if self.avg_recon_accuracy is not None:
            avg_recon_accuracy_str = "{:.0f}".format(self.avg_recon_accuracy * 100)

        rss_str = ""
        if self.comparison_rss is not None:
            rss_str = "{:.2f}".format(self.comparison_rss)

        return "min.expr.: {}, dataset corr.: {}, norm.: {}, z-score: {}, " \
               "log2: {},\navg.residuals: {}, avg.recon.accuracy.: {}%" \
               " N: {}, RSS: {}".format(self.min_expr,
                                        self.dataset_correction,
                                        self.normalize,
                                        self.zscore,
                                        self.log2,
                                        avg_resid_str,
                                        avg_recon_accuracy_str,
                                        self.comparison_n_samples,
                                        rss_str)

    def save_data_settings(self):
        data = {"date time": datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                "data_path": self.data_path,
                "signature_path": self.signature_path,
                "translate_path": self.translate_path
                }

        self.save_settings(self.outdir_path, 'data_settings', data)

    def save_all_settings(self):
        data = {"date time": datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                "data_path": self.data_path,
                "signature_path": self.signature_path,
                "translate_path": self.translate_path,
                "ground_truth_path": self.ground_truth_path,
                "outdir_path": self.outdir_path,
                "sample_to_dataset_path": self.sample_to_dataset_path,
                "sample_filter_path": self.sample_filter_path,
                "dataset_filter": self.dataset_filter,
                "filter1_shape_diff": self.filter1_shape_diff,
                "dataset_correction": self.dataset_correction,
                "reference_dataset": self.reference_dataset,
                "outsubdir_path": self.outsubdir_path,
                "min_expr": self.min_expr,
                "normalize": self.normalize,
                "zscore": self.zscore,
                "log2": self.log2,
                "decon_method": self.decon_method,
                "extension": self.extension,
                "real_info_per_celltype": self.real_info_per_celltype,
                "filter2_shape_diff": self.filter2_shape_diff,
                "sign_shift": self.sign_shift,
                "expr_shift": self.expr_shift,
                "n_samples": self.n_samples,
                "n_genes": self.n_genes,
                "n_ng_per_ct": self.n_ng_per_ct,
                "avg_residuals": self.avg_residuals,
                "avg_recon_accuracy": self.avg_recon_accuracy,
                "pred_info_per_celltype": self.pred_info_per_celltype,
                "comparison_n_samples": self.comparison_n_samples,
                "comparison_rss": self.comparison_rss}
        self.save_settings(self.outsubdir_path, 'all_settings', data)

    @staticmethod
    def save_settings(path, name, data):
        if path is None:
            print("No output path is set.")
            return

        outpath = os.path.join(path, '{}.json'.format(name))
        with open(outpath, 'w') as f:
            json.dump(data, f)
        f.close()
