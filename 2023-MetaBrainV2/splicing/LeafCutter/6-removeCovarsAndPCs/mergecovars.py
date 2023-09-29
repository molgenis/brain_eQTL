import sys
import gzip
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from numpy.linalg import LinAlgError

if len(sys.argv) < 4:
    print("Usage: covfile pcafile [number_of_columns_to_use_in_pcafile] outfile.txt.gz")
    sys.exit(0)

covFile = sys.argv[1]
pcaFile = sys.argv[2]
nPCs = int(sys.argv[3])
outfile = sys.argv[4]

print(covFile)
print(pcaFile)
print("nr of PCs:", sys.argv[4])

dfCov = pd.read_csv(covFile, sep='\t', index_col=0)
print("cov size:", dfCov.shape[0])
dfCov.index.name = 'Sample'

# Load pca matrix, samples on rows, nPCs on cols
useCols = [i for i in range(nPCs + 1)]
dfPCA = pd.read_csv(pcaFile, sep='\t', index_col=0, usecols=useCols)
print("pca size:", dfPCA.shape[0])
dfPCA.index.name = 'Sample'
dfPCA = dfPCA.reindex(dfCov.index)

dfCov = pd.concat([dfCov, dfPCA], axis=1)
print("merged size", dfCov.shape[0])
print("nr of cov after merge:", dfCov.shape[1])
print(dfCov)

dfCov.to_csv(outfile, sep='\t')
