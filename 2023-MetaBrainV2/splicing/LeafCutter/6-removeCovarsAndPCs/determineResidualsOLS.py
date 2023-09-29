import sys
import gzip
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from numpy.linalg import LinAlgError

import cProfile
import pstats
profiler = cProfile.Profile()

if len(sys.argv) < 4:
    print("Usage: splicefile covfile pcafile [number_of_columns_to_use_in_pcafile] outfile_residuals.txt.gz")
    sys.exit(0)

spliceFile = sys.argv[1]
covFile = sys.argv[2]
pcaFile = sys.argv[3]
nPCs = int(sys.argv[4])
outfileResiduals = sys.argv[5]

print(spliceFile)
print(covFile)
print(pcaFile)
print("nr of PCs:", sys.argv[4], flush=True)

# Load cov matrix
dfSplice = pd.read_csv(spliceFile, sep='\t', index_col=0)
dfSplice = dfSplice.transpose()
dfSplice.index.name = 'Sample'
print("splice size:", dfSplice.shape[0])

dfCov = pd.read_csv(covFile, sep='\t', index_col=0)
print("cov size:", dfCov.shape[0])
dfCov.index.name = 'Sample'
dfCov = dfCov.reindex(dfSplice.index)

# Load pca matrix, samples on rows, nPCs on cols
useCols = [i for i in range(nPCs + 1)]
dfPCA = pd.read_csv(pcaFile, sep='\t', index_col=0, usecols=useCols)
print("pca size:", dfPCA.shape[0])
dfPCA.index.name = 'Sample'
dfPCA = dfPCA.reindex(dfSplice.index)

dfCov = pd.concat([dfCov, dfPCA], axis=1)
print("merged size", dfCov.shape[0])
print("nr of cov after merge:", dfCov.shape[1])
print(dfCov)


sampleMap = {}
sampleList = []
sampleCount = 0
for sample in dfSplice.index:
    sampleMap[sample] = sampleCount
    sampleList.append(sample)
    sampleCount = sampleCount + 1

# iterate splice events
fCount = 0
lCount = 0
eCount = 0
wCount = 0


def checkColinearity(covariates):
    # check for variance
    toDrop = []
    for column in covariates:
        val = covariates[column].var()
        if val == 0:
            toDrop.append(column)
    tmp = covariates.drop(toDrop, axis=1)
    # check for correlation
    cor = tmp.corr()
    cols = tmp.columns
    toDrop = []
    for i in range(0, len(cols)):
        for j in range(i+1, len(cols)):
            colB = cols[j]
            val = cor.iat[i, j]
            if val > 0.9999:
                toDrop.append(colB)
    tmp = tmp.drop(toDrop, axis=1)

    removed = 1
    i = 0
    while removed > 0:
        removed = 0
        for covar in tmp:
            y_cov = tmp[covar]
            x_cov = tmp.drop(covar, axis=1)
            res = sm.OLS(y_cov.astype(float), x_cov.astype(float)).fit()
            rSquared = res.rsquared
            if rSquared > 0.999:
                # drop covar
                tmp = tmp.drop(covar, axis=1)
                removed = 1
                break
        i = i + 1
    return tmp


outfhResiduals = gzip.open(outfileResiduals, 'wt')
outfhResiduals.write("-\t" + "\t".join(sampleList) + "\n")

dfCov = checkColinearity(dfCov)
print(dfCov.shape)

for event in dfSplice:
    psi = dfSplice[[event]]
    nonna = psi.dropna(how='any')
    if nonna.empty:
        continue
    shape = nonna.shape
    outputResiduals = [np.nan] * len(psi)

    logitY = round(nonna[event])
    logitY = np.array(logitY)
    covars = dfCov.reindex(nonna.index)
    covars = checkColinearity(covars)
    X = covars.to_numpy()
    X = sm.add_constant(X)

    try:
        ols = sm.OLS(logitY, X)
        olsResults = ols.fit(disp=0)
        residuals = olsResults.resid
        idx = nonna.index
        for i in range(0, len(residuals)):
            actidx = sampleMap.get(idx[i])
            outputResiduals[actidx] = residuals[i]
        strout = "\t".join(str(v) for v in outputResiduals)
        outln = event + "\t"+strout+"\n"
        outfhResiduals.write(outln)
        wCount = wCount + 1
    except PerfectSeparationError:
        fCount = fCount + 1
    except LinAlgError:
        lCount = lCount + 1

    eCount = eCount + 1
    if eCount % 50 == 0:
        print("{} lines processed, {} failed, {} linalg errors, {} written".format(eCount, fCount, lCount, wCount))

outfhResiduals.close()
print("Done")
