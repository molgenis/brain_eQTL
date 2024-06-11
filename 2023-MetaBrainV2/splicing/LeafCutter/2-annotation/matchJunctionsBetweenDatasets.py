import gzip
import sys
from pathlib import Path

path = str(Path(__file__).parent.parent.parent.parent.absolute().__str__() + "/library/")
print("library path: " + path)
sys.path.insert(0, path)

from features.splicefeature2 import SpliceFeature
from intervaltree import Interval, IntervalTree

if len(sys.argv) < 4:
    print("Usage: junctions1file junctions2file outputfile")
    sys.exit()

infile1 = sys.argv[1]
infile2 = sys.argv[2]
outfile = sys.argv[3]

def getfh(file,mode):
    if file.endswith(".gz"):
        if mode == 'r':
            return gzip.open(file,'rt')
        elif mode == 'w':
            return gzip.open(file,'wt')
    else:
        if mode == 'r':
            return open(file,'r')
        elif mode == 'w':
            return open(file,'w')
    return None

def getJunctions(file):
    chrs = set()
    junctions = []
    print("Reading junctions from: "+file)
    fh = getfh(file,'r')
    for line in fh:
        junctionId = line.strip().split("\t",2)[0]
        idelems = junctionId.split(":")
        if len(idelems) > 1:
            junction = SpliceFeature.parseLeafCutter(junctionId)
            if not junction.chr.isAutosomal():
                chrs.add(idelems[0])
            else:
                junctions.append(junction)
    if len(chrs) > 0:
        print("Some junctions on these chromosomes were not parsed: ")
        for chr in chrs:
            print(chr)
    print(f"{len(junctions)} junctions loaded")
    return junctions

def toIntervalTree(junctions):
    print(f"Parsing {len(junctions)} junctions into intervaltree")
    intervalTreeByChr = {}
    chrs = set()
    for junction in junctions:
        tree = intervalTreeByChr.get(junction.chr)
        if tree is None:
            tree = IntervalTree()
        tree[junction.start:junction.stop] = junction
        chrs.add(junction.chr)
        intervalTreeByChr[junction.chr] = tree
    
    print("Loaded junctions per chromosome:")
    for chr in chrs:
        arr = intervalTreeByChr.get(chr)
        print(f"{chr}\t{len(arr)}")
    return intervalTreeByChr

def getOverlappingJunctions(feature, wiggle, tree):
    if tree is None:
        return None
    start = feature.start - wiggle
    if start < 0:
        start = 0
    stop = feature.stop + wiggle    
    intervals = tree[start:stop]
    junctions = []
    for interval in intervals:
        junctions.append(interval.data)
    return junctions



junctions1 = getJunctions(infile1)
junctions2 = getJunctions(infile2)
treesByChr = toIntervalTree(junctions2)

wiggle = 25 # 10 basepairs wiggle

nrWithOverlap = 0
nrWithtoutOverlap = 0
nrWithSingleOverlap = 0
nrWithMultipleOverlap = 0
nrTotal = 0
exactmatch = 0
fho = getfh(outfile,'w')
fho.write("Junction1\tBestMatch\tBestMatchSumDistance\tBestMatchStartDistance\tBestMatchStopDistance\tNrOptions\tOtherOptions(sumDistance)\n")
for junction1 in junctions1:
    chr = junction1.chr
    tree = treesByChr.get(chr)
    bestMatch = None
    bestMatchdSum = 100000000000

    if tree is not None: 
        allOptions = []
        junctionToDist = {}
        overlap = getOverlappingJunctions(junction1, wiggle, tree)
        if len(overlap) == 0:
            nrWithtoutOverlap += 1
        else:
            nrWithOverlap += 1
            if len(overlap) == 1:
                nrWithSingleOverlap += 1
                junction2 = overlap[0]
                bestMatch = junction2
                dSum = 0
                dSum += abs(junction1.start - junction2.start)
                dSum += abs(junction1.stop  - junction2.stop)
                bestMatchdSum = dSum
                allOptions.append(bestMatch)
            elif len(overlap) > 1:
                nrWithMultipleOverlap += 1
                # multiple overlaps; need to find best match!
    
                for junction2 in overlap:
                    dSum = 0
                    dSum += abs(junction1.start - junction2.start)
                    dSum += abs(junction1.stop  - junction2.stop)
                    if dSum < bestMatchdSum:
                        bestMatchdSum = dSum
                        bestMatch = junction2
                    allOptions.append(junction2)
                    junctionToDist[junction2] = dSum

        otherOptions = []
        otherOptionsStr = ""
        for option in allOptions:
            if option != bestMatch:
                otherOptions.append(option)
                if len(otherOptionsStr) == 0:
                    distance = junctionToDist.get(option)
                    otherOptionsStr = option.name+f" ({distance})"
                else:
                    otherOptionsStr += "; "+option.name+f" ({distance})"
        if len(otherOptions) == 0:
            otherOptionsStr = "-"
        if bestMatch is not None:
            # check whether best match is actual match
            if junction1.start == bestMatch.start and junction1.stop == bestMatch.stop:
                exactmatch += 1
            dStart = abs(junction1.start-bestMatch.start)
            dStop = abs(junction1.stop-bestMatch.stop)
            fho.write(f"{junction1.name}\t{bestMatch.name}\t{bestMatchdSum}\t{dStart}\t{dStop}\t{len(otherOptions)}\t{otherOptionsStr}\n")
        else:
            fho.write(f"{junction1.name}\t-\t-\t-\t-\t{len(otherOptions)}\t{otherOptionsStr}\n")
        nrTotal += 1
print(f"{nrWithOverlap} have overlap, {nrWithSingleOverlap} with single junction, {nrWithMultipleOverlap} with multiple junction, {nrWithtoutOverlap} have no overlap, {nrTotal} total")
print(f"{exactmatch} with exact match")
fho.close()