import glob
import os
import sys
import gzip
from collections import defaultdict

import pandas as pd

if len(sys.argv) < 3:
    print("Usage: input-dir annotation-file filter-file output-dir")
    sys.exit(0)

inputDir = sys.argv[1]
annotationFilePath = sys.argv[2]
filterFilePath = sys.argv[3]
outputDir = sys.argv[4]

sampleMap = {}
fh = gzip.open('cramToRNAseqId.txt.gz', 'rt')
for line in fh:
    elems = line.strip().split("\t")
    sampleMap[elems[1]] = elems[0]
fh.close()

selectedSamples = set()
filterFile = open(filterFilePath, 'r')
for line in filterFile:
    selectedSamples.add(line.strip())

samplesPerDataset = defaultdict(list)
annotationFile = open(annotationFilePath, 'r')
for line in annotationFile:
    elems = line.strip().split("\t")
    dataset = elems[2]
    sampleId = elems[1]
    if sampleId in selectedSamples:
        rnaSeqId = sampleMap[sampleId]
        samplesPerDataset[dataset].append(rnaSeqId)
annotationFile.close()

countFiles = glob.glob(inputDir + "/*counts.gz.phen*")
for countFile in countFiles:
    print("Reading", countFile)
    countDf = pd.read_csv(countFile, sep="\t", index_col=3)

    for dataset in samplesPerDataset:
        samples = samplesPerDataset[dataset]
        countDataset = countDf.loc[:, countDf.columns.isin(samples)]
        print(countDataset.shape, "Samples")
        filename = os.path.basename(countFile)
        filename = dataset + "_" + filename
        print("Writing", filename)
        countDataset.to_csv(outputDir + "/" + filename, sep='\t', na_rep="nan")

