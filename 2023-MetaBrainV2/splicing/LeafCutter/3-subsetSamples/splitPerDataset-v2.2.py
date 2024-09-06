import gzip
import sys
import os
import numpy
import argparse
import math
from enum import Enum
from pprint import pprint
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

# argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-c","--counts", help="Leafcutter cluster count file", required=True)
parser.add_argument("-s","--samplerename", help="Replace e.g. BAM ids to something more friendly")
parser.add_argument("-d","--datasetfile", help="Sample to dataset definition file", required=False)
parser.add_argument("-i","--sampleinclude", help="Sample include list (after --samplerename, when applied)")
parser.add_argument("-o","--out", help="Output prefix",required=True)
parser.add_argument("--minNrDatasets", help="Minimal number of datasets per junction",default=2)
parser.add_argument("--minAvgPSI", help="Minimal PSI value",default=0.01)
parser.add_argument("--minStDevPSI", help="Minimal PSI stdev value",default=0.005)
# parser.add_argument("--averageImpute", help="Replace NaN values with average calculated over non-NaN samples",action='store_true')
# parser.add_argument("--averageImputeMissingDatasets", help="Replace NaN values with average calculated over non-NaN samples, also for datasets that have not data at all, or those that didn't pass QC (like the default in LeafCutter)",action='store_true')
parser.add_argument("--averageImputePerDataset", help="Replace NaN values with average calculated over non-NaN samples per dataset",action='store_true')
parser.add_argument("--minObersvationsPerDataset",help="Minimal number of non-NaN samples per dataset (or a proportion thereof)", default=30)
# parser.add_argument("--minProportionOfObservationsPerDataset",help="Minimal number of non-NaN samples per dataset", default=0.4)
parser.add_argument("--minNrOfReads",help="Minimal number counts in the junction denominator per sample (values below this number will be replaced by NaN)",default=10)
# parser.add_argument("--writeExtraFiles",help="Write additional output matrices (raw unfiltered PSI values etc)",action='store_true')
parser.add_argument("--removeNonStandardChr",help="Remove non-autosomal and non-X or non-Y splice events",action='store_true')
parser.add_argument("--removeNonAutosomal",help="Remove non-autosomal splice events",action='store_true')
parser.add_argument("--usePseudoCount",help="Use pseudocount in PSI calculation",action='store_true')
parser.add_argument("--centerAndScale",help="Center by mean, scale by stdev",action='store_true')

args = vars(parser.parse_args())

psifile = args["counts"]
cramfile  = args["samplerename"]
etdfile = args["datasetfile"]  # expression to dataset
samplelimitfile = args["sampleinclude"]  # expression to dataset
outprefix  = args["out"]
minnrdatasets = int(args["minNrDatasets"])

minObersvationsPerDataset = float(args["minObersvationsPerDataset"])
minAvgPsi = float(args["minAvgPSI"])
minStDevPsi = float(args["minStDevPSI"])
maxAvgPsi = 1-(float(args["minAvgPSI"]))
minnrdatasets = int(args["minNrDatasets"])
minNrOfReads = int(args["minNrOfReads"])
threads = 10
debug = False

# imputeAverage = args["averageImpute"]
# imputeAverageMissingDatasets = args["averageImputeMissingDatasets"]
imputeAveragePerDataset = args["averageImputePerDataset"]
# writeExtraFiles = args["writeExtraFiles"]
removeNonAutosomal = args["removeNonAutosomal"]
removeNonStandardChr = args["removeNonStandardChr"]
usePseudoCount = args["usePseudoCount"]
centerAndScale = args["centerAndScale"]

# if imputeAverage and imputeAveragePerDataset:
#     print("Error: cannot average impute over all samples and per dataset at the same time!")
#     sys.exit(-1)

pprint(args)

# main functions
def calcPsi(num, denom):
    # returns list: readct, psi unfiltered, psi filtered
    if denom < 1 or denom < minNrOfReads:
        return numpy.nan
    psi = 0
    if usePseudoCount:
        psi = (num+0.5)/(denom+0.5)
    else:
        psi = (num)/(denom)
    return psi

def meanAndVar(vals):
    sum = 0
    ct = 0
    for v in vals:
        if not numpy.isnan(v):
            sum += v
            ct += 1
    if ct > 1:
        sum /= ct
        varV = var(vals,sum)
        return [sum, varV]
    else:
        return [0,0]

def var(x, meanx):
    var = 0
    nonnan = 0
    for v in x:
        if not numpy.isnan(v):
            var +=  (v - meanx) * (v-meanx)
            nonnan += 1
    if nonnan > 1:
            var /= (nonnan-1)
    else:
            var = 0
    return var

def f(v):
    return str(round(v,3))

def scale(v):
    mean, var = meanAndVar(v)
    stdev = math.sqrt(var)
    return scaleMeanStDev(v,mean,stdev)

def scaleMeanStDev(v, mean, stdev):
    for i in range(0,len(v)):
        if not numpy.isnan(v[i]):
            v[i] -= mean
            v[i] /= stdev
    return v

def loadMap(cramfile,keycol,valcol):
    cramMap = {}
    fh = gzip.open(cramfile,'rt')
    for line in fh:
        elems = line.strip().split("\t")
        cramMap[elems[keycol]] = elems[valcol]
    fh.close()

    # custom mappings
    cramMap["NEUYV496XLP_CGND_HRA_00084-2"] = "HRA_00084"
    cramMap["AN11864_AN11864_ba41-42-22"] = "AN11864_ba41.42.22"

    print("{} sample mappings read from {}".format( len(cramMap), cramfile ) )
    return cramMap

def loadSet(file):
    sampleLimit = set()
    fh = open(samplelimitfile,'r')
    for line in fh:
        id = line.strip()
        sampleLimit.add(id)
    fh.close()
    return sampleLimit

def loadETD(etdfile):
    samplesPerDataset = {}
    sampleToDataset = {}
    selectedSamples = set()
    fh = open(etdfile,'r')
    for line in fh:
        elems = line.strip().split()
        ds = elems[2]
        sa = elems[1]
        if sampleLimit is None or sa in sampleLimit:
            data = samplesPerDataset.get(ds)
            if data is None:
                data = set()
            data.add(sa)
            samplesPerDataset[ds] = data
            selectedSamples.add(sa)
            sampleToDataset[sa] = ds
    fh.close()
    return samplesPerDataset, sampleToDataset, selectedSamples

def fill(arr, val):
    for i in range(0,len(arr)):
        arr[i] = val

def processDataset(d, colsPerDataset, minNrObservationsPerDataset, elems):
    cols = colsPerDataset[d]
    dataset = availableDatasets[d]
    dsMinNrObs = minNrObservationsPerDataset[d]
    psivals = [numpy.nan] * len(cols)
    fill(psivals, numpy.nan)

    passqc = False
    mean = 0
    variance = 0
    stdev = 0
    nrNotMissing = 0
    
    
    for i in range(0,len(cols)):
        col = cols[i]
        num, denom = elems[col].split("/")
        num = float(num)
        denom = float(denom)
        psi = calcPsi(num, denom)
        psivals[i] = psi
        if not numpy.isnan(psi):
            nrNotMissing += 1
            mean += psi
    if nrNotMissing > 0:
        mean /= nrNotMissing
    # Minimal number of non-NaN samples per dataset (if specified < 1, will be interpreted as proportion of total)
    # leafcutter: skip if more than 40% missing, or skip if less than 60% not missing
    # 
    if nrNotMissing < dsMinNrObs:
        fill(psivals, numpy.nan)
    else:
        if imputeAveragePerDataset and nrNotMissing < len(cols):
            for i in range(0,len(cols)):
                if numpy.isnan(psivals[i]):
                    psivals[i] = mean
        
        mean, variance = meanAndVar(psivals)
        #stdev = math.sqrt(variance)
        stdev = numpy.std(psivals)
        
        #if mean < minAvgPsi or mean > maxAvgPsi or stdev < minStDevPsi:
        if stdev < minStDevPsi:
            # wipe values if not passing thresholds
            fill(psivals, numpy.nan)
        else:
            passqc = True
            if centerAndScale:
                psivals = scaleMeanStDev(psivals,mean,stdev)
    if debug:
        print(f"{dataset}\t{len(psivals)}\t{nrNotMissing}\t{mean}\t{stdev}\t{passqc}")
    psivalstr = "\t".join( str(x) for x in psivals )
    return d, psivals, psivalstr, passqc

def processLine(line, colsPerDataset, minNrObservationsPerDataset):
    
    elems = line.strip().split(" ")
    ctr = 0
    rowid = elems[0] # always append first column
    spliceId = rowid.split(":")
    chr = spliceId[0].replace("chr","").lower()
    sexchr = False
    autosomal = False
    if chr == "x" or chr == "y":
        sexchr = True
    else:
        try:
            chr = int(chr)
            if chr > 0 and chr < 23:
                autosomal = True
        except:
            pass

    if removeNonAutosomal and not autosomal:
        return [False, None]
    if removeNonStandardChr and not (autosomal or sexchr):
        return [False, None]

    outln = rowid
    nrDatasetsPassingQC = 0
    
    for d in range(0, len(availableDatasets)):
        
        d2, psivals, psivalstr, passqc = processDataset(d, colsPerDataset, minNrObservationsPerDataset, elems)
        if passqc:
            nrDatasetsPassingQC+=1
        outln += "\t"+psivalstr    

    if nrDatasetsPassingQC < minnrdatasets:
        return [False, None]
    return [True, outln]

cramMap = None
if cramfile is not None:
    cramMap = loadMap(cramfile,0,1)

sampleLimit = None
if samplelimitfile is not None:
    sampleLimit = loadSet(samplelimitfile)
    print("{} sample ids read from {}".format(len(sampleLimit), samplelimitfile))

samplesPerDataset = None
sampleToDataset = None
selectedSamples = None

# link samples to datasets
if etdfile is not None:
    samplesPerDataset, sampleToDataset, selectedSamples = loadETD(etdfile)
    print("{} samples selected from {}. {} datasets".format( len(selectedSamples), etdfile, len(samplesPerDataset.keys()) ) )

# parse header
fh = gzip.open(psifile,'rt')
header = fh.readline().strip().split(" ")

datasetCols = {}
datasetSamples = {}

for i in range(1,len(header)):
    # replace ID with actual RNA-ID
    cram = header[i]
    sample = cram
    if cramMap is not None:
        sample = cramMap.get(cram)
        if sample is None:
            print("Warning: "+cram+" was not found in sample conversion file")
            sample = cram
            ok = True
    else:
        ok =  True

    dataset = sampleToDataset.get(sample)
    if dataset is not None:
        # do something
        colsForDs = datasetCols.get(dataset)
        if colsForDs is None:
            colsForDs = []
        colsForDs.append(i)
        datasetCols[dataset] = colsForDs
        samplesForDs = datasetSamples.get(dataset)
        if samplesForDs is None:
            samplesForDs = []
        samplesForDs.append(sample)
        datasetSamples[dataset] = samplesForDs


availableDatasets = list(datasetCols.keys())
availableDatasets.sort()


print("Number of included samples per included dataset:")
colsPerDataset = []
minNrObservationsPerDataset = []

outHeader = header[0]
totalSamples = 0
for dataset in availableDatasets:
    cols = datasetCols.get(dataset)
    print(f"{dataset}\t{len(cols)}")
    totalSamples += len(cols)
    samples = datasetSamples.get(dataset)
    outHeader +="\t"+"\t".join(samples)
    # outHeader.append(header[cols])
    colsPerDataset.append(cols)
    dsPresentThreshold = minObersvationsPerDataset
    if minObersvationsPerDataset < 1:
        dsPresentThreshold = minObersvationsPerDataset * len(cols)
    minNrObservationsPerDataset.append(dsPresentThreshold)
outHeader = outHeader+"\n"
print(f"Total samples: {totalSamples}")

fhoPsiFiltered = gzip.open(outprefix+"-PSI-filtered.txt.gz",'wt',3)
fhoPsiFiltered.write(outHeader)
lctr = 0
written = 0

pool = ProcessPoolExecutor(max_workers=threads)
nrJobs = 1000
futures = [None] * nrJobs
jobctr = 0
for line in fh:


    futures[jobctr] = pool.submit(processLine,line, colsPerDataset, minNrObservationsPerDataset)
    jobctr += 1
    if jobctr == nrJobs:
        for future in futures:
            write, outln = future.result()
            if write:
                fhoPsiFiltered.write(outln+"\n")
                written += 1
        jobctr = 0
        # print(availableDatasets[d]+"\t"+str(mean)+"\t"+str(stdev)+"\t"+str(len(cols)) +"\t"+str(nrNotMissing) +"\t"+str(nrDatasetsPassingQC) +"\t"+str(dsMinNrObs)) 
    
    lctr += 1
    if lctr % 1000 == 0:
        perc = (written / lctr)*100
        print(f"{lctr} lines parsed, {written} written, {f(perc)} %",end='\r')
if jobctr > 0:
    for i in range(0, jobctr):
        write, outln = futures[i].result()
        if write:
            fhoPsiFiltered.write(outln+"\n")
            written += 1

perc = (written / lctr)*100
print(f"{lctr} lines parsed, {written} written, {f(perc)} %",end='\n')
fhoPsiFiltered.close()
pool.shutdown()