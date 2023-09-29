import glob
import os
import sys

import pandas as pd
from scipy.stats import stats

if len(sys.argv) != 5:
    print("Usage: pca-input-dir z-score-cut-off psi-file-dir output-dir")
    sys.exit(0)

pcaDir = sys.argv[1]
zScoreCutOff = int(sys.argv[2])
psiDir = sys.argv[3]
outputDir = sys.argv[4]

samplesToRemove = set()

pcaFiles = glob.glob(pcaDir + "/*PCs.txt")
for pcaFile in pcaFiles:
    print("parsing ", pcaFile)
    pcaDf = pd.read_csv(pcaFile, sep='\t', usecols=[0, 1, 2], index_col=0)

    zScores = pcaDf.apply(stats.zscore)
    zScoresToExclude = zScores[(abs(zScores.PC1) > zScoreCutOff) | (abs(zScores.PC2) > zScoreCutOff)]
    samplesToExcludeEvent = zScoresToExclude.index
    print(len(samplesToExcludeEvent), "outlier samples")
    samplesToRemove.update(list(samplesToExcludeEvent))

print(len(samplesToRemove), "unique outlier samples")

psiFiles = glob.glob(psiDir + "/*.gz")
for psiFile in psiFiles:
    print("Removing outliers from", psiFile)
    psiDf = pd.read_csv(psiFile, sep="\t", index_col=0)
    psiDfFiltered = psiDf.drop(samplesToRemove, axis=1)

    fileName = os.path.basename(psiFile)
    fileName = os.path.splitext(fileName)[0]
    fileName = os.path.splitext(fileName)[0]
    psiDfFiltered.to_csv(outputDir + "/" + fileName + "_filtered.txt.gz", sep='\t', na_rep="nan")

