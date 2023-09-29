import sys
import gzip

import pandas as pd

if len(sys.argv) < 2:
    print("Usage: input-tsv annotation-file output-tsv")
    sys.exit(0)

inputTsvPath = sys.argv[1]
annotationFilePath = sys.argv[2]
outputTsvPath = sys.argv[3]

sampleMap = {}
sampleMap["-"] = "-" # for index column 
fh = gzip.open(annotationFilePath, 'rt')
for line in fh:
    elems = line.strip().split("\t")
    sampleMap[elems[0]] = elems[1]
fh.close()

countDf = pd.read_csv(inputTsvPath, sep="\t")
countDf.columns = [sampleMap[x] for x in countDf.columns]
countDf.to_csv(outputTsvPath, sep="\t", index=False, compression='gzip')
