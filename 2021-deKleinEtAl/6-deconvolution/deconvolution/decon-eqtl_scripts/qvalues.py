"""
File:         qvalues.py
Created:      2020/07/21
Last Changed: 2020/07/22
Author:       M.Vochteloo

This is a custom written Python implementation of R Bioconductor’s qvalue
package Version 2.22.0. Qvalues is written under the GNU LESSER GENERAL PUBLIC
LICENSE. A copy of the GNU LESSER GENERAL PUBLIC
LICENSE can be found here <https://www.gnu.org/licenses/lgpl-3.0.txt>. This
reimplementation is also licensed under the GNU LESSER GENERAL PUBLIC
LICENSE.

The q-values package in R on which from which this implementation is an
derative is written by John D. Storey, Andrew J. Bass, Alan Dabney and
David Robinson.

Reference:
John D. Storey, Andrew J. Bass, Alan Dabney and David Robinson (2020).
qvalue: Q-value estimation for false discovery rate control. R package
version 2.22.0. http://github.com/jdstorey/qvalue
"""

# Standard imports.

# Third party imports.
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from scipy.stats import norm, gaussian_kde

# Local application imports.


def qvalue(p, fdr_level=None, pfdr=False, lfdr_out=True, pi0=None, **kwargs):
    call = ["fdr_level={}".format(fdr_level) if fdr_level is not None else "",
            "pfdr={}".format(pfdr) if pfdr else "",
            "lfdr_out={}".format(lfdr_out) if not lfdr_out else "",
            "pi0={}".format(pi0) if pi0 is not None else "",
            ", ".join(["{}={}".format(key, value) for key, value in kwargs.items()])
            ]
    call = [x for x in call if x != ""]

    retval = {
        "call": "qvalue(p={}{})".format(p, (", " if len(call) > 0 else "") + ", ".join(call)),
        "pi0": pi0,
        "qvalues": None,
        "pvalues": np.copy(p),
        "lfdr": None,
        "fdr_level": fdr_level,
        "significant": None,
        "pi0_lambda": None,
        "lambda": None,
        "pi0_smooth": None
    }

    qvals_out = np.copy(p)
    lfdr_out_a = np.copy(p)
    rm_na = ~np.isnan(p)
    p = p[rm_na]
    if np.min(p) < 0 or np.max(p) > 1:
        print("p-values not in valid range [0, 1].")
        return retval
    if fdr_level is not None and (fdr_level <= 0 or fdr_level > 1):
        print("'fdr.level' must be in (0, 1].")
        return retval
    if pi0 is not None and (pi0 <= 0 or pi0 > 1):
        print("pi0 is not (0,1]")
        return retval

    if pi0 is None:
        pi0s = pi0est(p=p, **kwargs)
        pi0 = pi0s["pi0"]
        retval.update(pi0s)

    m = np.size(p)
    i = np.arange(m, 0, -1)
    o = np.argsort(p)[::-1]
    ro = np.argsort(o)

    if pfdr:
        rank = i * (1 - np.power((1 - p[o]), m))
    else:
        rank = i
    qvals = pi0 * np.min((np.ones_like(p), pd.Series(p[o] * m / rank).cummin().to_numpy()), axis=0)[ro]
    qvals_out[rm_na] = qvals

    if lfdr_out:
        lfdr_out_a[rm_na] = lfdr(p=p, pi0=pi0, **kwargs)
    else:
        lfdr_out_a = None

    retval.update({"qvalues": qvals, "lfdr": lfdr_out_a})
    if fdr_level is not None:
        retval["significant"] = qvals <= fdr_level

    return retval


def lfdr(p, pi0=None, trunc=True, monotone=True, transf="probit", adj=1.5,
         eps=1e-08, **kwargs):
    if transf not in ["probit", "logit"]:
        print("ERROR: transf invalid")
        return None
    lfdr_out = np.empty_like(p)
    lfdr_out[:] = np.nan
    rm_na = ~np.isnan(p)
    p = p[rm_na]

    if np.min(p) < 0 or np.max(p) > 1:
        print("P-values not in valid range [0,1].")
        return lfdr_out

    if pi0 is None:
        pi0 = pi0est(p=p, **kwargs)["pi0"]

    if transf == "probit":
        p[p < eps] = eps
        p[p > 1 - eps] = 1 - eps
        x = norm.ppf(p)
        mydx, mydy = density(x, adjust=adj)
        y = smooth_spline(x_train=mydx, y_train=mydy, x_test=x)
        lfdr_values = pi0 * norm.pdf(x) / y
    elif transf == "logit":
        x = np.log((p + eps) / (1 - p + eps))
        mydx, mydy = density(x, adjust=adj)
        y = smooth_spline(x_train=mydx, y_train=mydy, x_test=x)
        dx = np.exp(x) / np.power((1 + np.exp(x)), 2)
        lfdr_values = (pi0 * dx) / y
    else:
        print("ERROR: transf must be one of \"probit\" or \"logit\".")
        return lfdr_out

    if trunc:
        lfdr_values[lfdr_values > 1] = 1

    if monotone:
        o = np.argsort(p)
        ro = np.argsort(o)
        lfdr_values = pd.Series(lfdr_values[o]).cummax().to_numpy()[ro]

    lfdr_out[rm_na] = lfdr_values

    return lfdr_out


def pi0est(p, lambda_values=np.arange(0.05, 1, 0.05), pio_method="smoother",
           smooth_df=3, smooth_log_pi0=False, **kwargs):
    retval = {"pi0": None,
              "pi0_lambda": None,
              "lambda": None,
              "pi0_smooth": None}
    if pio_method not in ["smoother", "bootstrap"]:
        print("ERROR: pio method invalid")
        return retval
    if not type(lambda_values) is np.ndarray:
        lambda_values = np.array([lambda_values])

    rm_na = ~np.isnan(p)
    p = p[rm_na]
    m = np.size(p)
    lambda_values = np.sort(lambda_values)
    ll = np.size(lambda_values)

    retval["lambda"] = lambda_values

    if np.min(p) < 0 or np.max(p) > 1:
        print("ERROR: p-values not in valid range [0, 1].")
        return retval
    if 1 < ll < 4:
        print("ERROR: length(lambda)= {}. If length of lambda greater than"
              " 1, you need at least 4 values.".format(ll))
        return retval
    if np.min(lambda_values) < 0 or np.max(lambda_values) >= 1:
        print("ERROR: Lambda must be within [0, 1).")
        return retval

    if ll == 1:
        pi0 = np.mean(p >= lambda_values[0]) / (1 - lambda_values[0])
        pi0_lambda = pi0
        pi0 = min(pi0, 1)
        pi0_smooth = None
    else:
        bins = np.digitize(p, lambda_values)
        bincounts = dict(zip(*np.unique(bins, return_counts=True)))
        frequencies = np.array([bincounts[x] if x in bincounts else 0 for x in np.arange(ll, 0, -1)]) # for some reason bin 0 is excluded
        pi0 = np.cumsum(frequencies) / (m * (1 - lambda_values[::-1]))
        pi0 = pi0[::-1]
        pi0_lambda = pi0
        if pio_method == "smoother":
            if smooth_log_pi0:
                pi0 = np.log(pi0)
            pi0_smooth = smooth_spline(x_train=lambda_values, y_train=pi0, df=smooth_df)
            pi0 = min(pi0_smooth[-1], 1)
            retval["pi0_smooth"] = pi0_smooth
        elif pio_method == "bootstrap":
            minpi0 = np.quantile(pi0, q=0.1)
            W = np.array([np.sum(p >= l) for l in lambda_values], dtype=np.float64)
            mse = (W / (np.power(m, 2) * np.power((1 - lambda_values), 2))) * (1 - W / m) + np.power((pi0 - minpi0), 2)
            pi0 = min(pi0[mse == np.min(mse)][0], 1)
            pi0_smooth = None
        else:
            print("ERROR: pi0.method must be one of \"smoother\" or \"bootstrap\".")
            return retval

    if pi0 <= 0:
        print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")

    retval.update({"pi0": pi0, "pi0_lambda": pi0_lambda, "pi0_smooth": pi0_smooth})

    return retval


def smooth_spline(x_train, y_train, x_test=None, **kwargs):
    r_x_train = robjects.FloatVector(x_train)
    r_y_train = robjects.FloatVector(y_train)
    r_x_test = r_x_train
    if x_test is not None:
        r_x_test = robjects.FloatVector(x_test)
    spline_xy = robjects.r['smooth.spline'](r_x_train, r_y_train, **kwargs)
    return np.array(robjects.r['predict'](spline_xy, r_x_test).rx2('y'))


def density(x, adjust=1., n=512, cut=3):
    """
    https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/density.R
    """
    if len(x) < 2:
        print("need at least 2 data points")
        return None, None

    # Calculate bw.nrd0
    # somehow the std is slightly different in R
    hi = np.std(x)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    lo = min(hi, iqr/1.34)
    if not ((lo == hi) or (lo == abs(x[0])) or (lo == 1)):
        lo = 1
    bw_nrd0_val = 0.9 * lo * np.power(np.size(x), -0.2)

    bw = adjust * bw_nrd0_val
    if bw <= 0:
        print("'bw' is not positive.")
        return None, None

    start = np.min(x) - bw * cut
    stop = np.max(x) + bw * cut
    x_eval = np.linspace(start=start, stop=stop, num=n)
    kernel = gaussian_kde(x, bw_method=bw / hi)
    y_kernel = kernel(x_eval)

    return x_eval, y_kernel


def emp_pvals(stat, stat0, pool=True):
    m = np.size(stat)
    n = np.shape(stat0)[1]

    if pool:
        if stat0.ndim > 2:
            print("stat0 must be a 2D matrix or a 1D vector.")
            return None
        elif stat0.ndim == 2:
            stat0 = stat0.flatten(order='F')

        m0 = np.size(stat0)
        v = np.hstack((np.ones(m, dtype=bool), np.zeros(m0, dtype=bool)))
        v = v[np.argsort(np.hstack((stat, stat0)))[::-1]]
        u = np.arange(1, np.size(v) + 1)
        w = np.arange(1, m + 1)
        p = (u[v == True] - w) / m0
        # R uses method = 'average' but this doesn't work if there are ties.
        p = p[pd.Series(-stat).rank(method="max").to_numpy(dtype=np.int) - 1]
        p[p < 1 / m0] = 1 / m0
    else:
        if stat0.ndim == 1:
            print("stat0 must be a matrix.")
            return None

        if stat0.shape[0] != m:
            print("Number of rows of stat0 must equal length of stat.")
            return None

        if n == m:
            stat0 = np.transpose(stat0)

        stat0 = (stat0 - stat[:, np.newaxis]) >= 0
        p = np.mean(stat0, axis=1)
        p[p < 1 / stat0.shape[1]] = 1 / stat0.shape[1]

    return p


def summary(object, cuts=(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1)):
    # Call
    print("")
    print("Call: {}".format(object["call"]))
    print("")
    # Proportion of nulls
    print("pi0:\t{}".format(object["pi0"]))
    print("")
    # Number of significant values for p-value, q-value and local FDR
    print("Cumulative number of significant calls:")
    print("")
    rm_na = ~np.isnan(object["pvalues"])
    pvalues = object["pvalues"][rm_na]
    qvalues = object["qvalues"][rm_na]
    lfdr_values = object["lfdr"][rm_na]
    counts = pd.DataFrame([[np.sum(x < c) for c in cuts] for x in [pvalues, qvalues, lfdr_values]],
                          index=["p-value", "q-value", "local FDR"],
                          columns=["<{}".format(x) for x in cuts])
    print(counts)
    print("")