"""
File:         perform_deconvolution.py
Created:      2020/06/29
Last Changed: 2021/09/30
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
import os

# Third party imports.
from scipy.optimize import nnls
import pandas as pd
import numpy as np

# Local application imports.


class PerformDeconvolution:
    def __init__(self, settings, signature, expression):
        self.decon_method = settings.get_decon_method()
        self.outdir = settings.get_outsubdir_path()
        self.signature = signature
        self.expression = expression

        self.deconvolution_raw = None
        self.deconvolution = None
        self.rss = None
        self.recon_accuracy = None

    def work(self):
        print("Performing deconvolution using '{}'".format(self.decon_method))
        decon_function = None
        if self.decon_method == "NNLS":
            decon_function = self.nnls
        else:
            print("Unexpected deconvolution method.")
            exit()

        decon_data = []
        rss_data = []
        recon_accuracy_data = []
        for index, sample in self.expression.T.iterrows():
            # Model.
            proportions, rnorm = decon_function(self.signature, sample)

            # Calculate reconstruction accuracy.
            recon_accuracy = self.calc_reconstruction_accuracy(y=sample,
                                                               X=self.signature,
                                                               betas=proportions)

            # Save
            decon_data.append(proportions)
            rss_data.append(rnorm * rnorm)
            recon_accuracy_data.append(recon_accuracy)

        deconvolution_raw = pd.DataFrame(decon_data,
                                         index=self.expression.columns,
                                         columns=self.signature.columns)
        rss = pd.Series(rss_data, index=self.expression.columns)
        recon_accuracy = pd.Series(recon_accuracy_data, index=self.expression.columns)

        # Make weights sum up to one.
        deconvolution = self.perform_sum_to_one(deconvolution_raw)

        # Save.
        self.deconvolution_raw = deconvolution_raw
        self.deconvolution = deconvolution
        self.rss = rss
        self.recon_accuracy = recon_accuracy
        self.save()

    @staticmethod
    def nnls(A, b):
        return nnls(A, b)

    @staticmethod
    def calc_reconstruction_accuracy(y, X, betas):
        y_hat = np.dot(X, betas)
        residuals = y - y_hat
        residuals_norm = np.linalg.norm(residuals)
        y_hat_norm = np.linalg.norm(y)
        return 1 - (residuals_norm * residuals_norm) / (y_hat_norm * y_hat_norm)

    @staticmethod
    def perform_sum_to_one(X):
        print("Making weights sum-to-one")
        return X.divide(X.sum(axis=1), axis=0)

    def get_deconvolution(self):
        return self.deconvolution

    def get_rss(self):
        return self.rss

    def get_recon_accuracy(self):
        return self.recon_accuracy

    def get_avg_rss(self):
        return self.rss.mean()

    def get_avg_recon_accuracy(self):
        return self.recon_accuracy.mean()

    def get_info_per_celltype(self):
        means = self.deconvolution.mean(axis=0).to_dict()
        stds = self.deconvolution.std(axis=0).to_dict()

        info = {}
        for key in means.keys():
            info[key] = (means[key], stds[key])
        return info

    def save(self):
        self.deconvolution_raw.to_csv(os.path.join(self.outdir, 'deconvolution_raw.txt.gz'),
                                      compression="gzip",
                                      sep="\t",
                                      header=True,
                                      index=True)
        print("\tsaved dataframe: deconvolution_raw.txt.gz "
              "with shape: {}".format(self.deconvolution.shape))

        self.deconvolution.to_csv(os.path.join(self.outdir, 'deconvolution.txt.gz'),
                                  compression="gzip",
                                  sep="\t",
                                  header=True,
                                  index=True)
        print("\tsaved dataframe: deconvolution.txt.gz "
              "with shape: {}".format(self.deconvolution.shape))

    def print_info(self):
        print("Average RSS: {:.2f}".format(self.get_avg_rss()))
        print("Average reconstruction accuracy: {:.2f}%".format(self.get_avg_recon_accuracy()))
        print("Average weights per celltype:")
        for celltype, (mean, std) in self.get_info_per_celltype().items():
            print("\t{:20s}: mean = {:.2f} std = {:.2f}".format(celltype, mean, std))
        print("")