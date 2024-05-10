import gzip
import sys


# sys.path.append("../../../library/")
from pathlib import Path

path = str(Path(__file__).parent.parent.parent.parent.absolute().__str__() + "/library/")
print("library path: " + path)
sys.path.insert(0, path)

from features.splicefeature2 import SpliceFeature

if len(sys.argv) < 3:
    print("Usage: splicefile output [bpwiggle: default=10]")
    sys.exit()

splicefile = sys.argv[1]
outf = sys.argv[2]
wiggle = 10
if len(sys.argv) > 3:
    wiggle = int(sys.argv[3])


def getfh(file):
    if file.endswith(".gz"):
        return gzip.open(file, 'rt')
    return open(file)

# load junctions and clusters
junctions = []
clusters = {}
print("Loading splice junctions from: " + splicefile)
fh = getfh(splicefile)
fh.readline()
lctr = 0
for line in fh:
    elems = line.strip().split("\t", 2)
    if len(elems) == 1:
        elems = line.strip().split(" ", 2)
    id = elems[0]
    # print(id)
    junction = SpliceFeature.parseLeafCutter(id)
    if junction.chr.isAutosomal():
        junctions.append(junction)
        clusterid = junction.clusterId
        junctionsInCluster = clusters.get(clusterid)
        if junctionsInCluster is None:
            junctionsInCluster = []
        junctionsInCluster.append(junction)
        clusters[clusterid] = junctionsInCluster
    lctr += 1
    if lctr % 50000 == 0:
        print("{} lines parsed, {} loaded, {} clusters".format(lctr, len(junctions), len(clusters)), end='\r')
        break
print("{} lines parsed, {} loaded, {} clusters".format(lctr, len(junctions), len(clusters)), end='\n')
fh.close()

for cluster in clusters.keys():
    junctions = clusters.get(cluster)
    overlap = True
    junctionmap = {}
    while overlap:
        overlap = False
        tmpjunctions = set()
        for i in range(0,len(junctions)):
            junctioni = junctions[i]
            starti = junctioni.start
            stopi = junctioni.stop

            mergedfeature = SpliceFeature()
            mergedfeature.chr = junctioni.chr
            mergedfeature.start = starti
            mergedfeature.stop = stopi
            tmpjunctions.add(mergedfeature)
            for j in range(i+1,len(junctions)):
                junctionj = junctions[j]
                startj = junctionj.start
                stopj = junctionj.stop
                
                wigglefeature = SpliceFeature()
                wigglefeature.chr = junctionj.chr
                wigglefeature.start = junctionj.start - wiggle
                wigglefeature.stop = junctionj.stop + wiggle
                if mergedfeature.overlaps(wigglefeature):
                    # merge features
                    if startj < starti:
                        mergedfeature.start = startj
                    if stopj > stopi:
                        mergedfeature.stop = stopj
                    overlap = True
                    
                else:
                    tmpjunctions.add(junctionj)
        junctions = tmpjunctions
            
            