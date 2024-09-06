import gzip
import sys
import os

def getrfh(fname):
    if fname[-3:] == ".gz":
        return gzip.open(fname, 'rt')
    else:
        return open(fname,'r'   )

def getwfh(fname):
    if fname[-3:] == ".gz":
        return gzip.open(fname, 'wt', 3)
    else:
        return open(fname,'w')   

if len(sys.argv) < 4:
    print("Usage: infile.txt.gz sampleToDs.txt outdir")
    sys.exit()

infile = sys.argv[1]
samplefile = sys.argv[2]
outdir = sys.argv[3]


sampleToDataset = {}

fh = getrfh(samplefile)
for line in fh:
    elems = line.strip().split("\t")
    if len(elems) < 3:
        print("Expected three columns: genotype rnaseq dataset")
        sys.exit()
    sample = elems[1]
    dataset = elems[2]
    sampleToDataset[sample] = dataset
fh.close()

fh = getrfh(infile)
header = fh.readline().split()
print(f"{infile} has {len(header)} header size")
datasetSamples = {}
datasetSampleIdx = {}

overlap = 0
for sampleIdx in range(1,len(header)):
    sample = header[sampleIdx]
    if sample in sampleToDataset.keys():
        dataset = sampleToDataset.get(sample)
        colIds = datasetSampleIdx.get(dataset)
        samples = datasetSamples.get(dataset)
        if colIds is None:
            colIds = []
            samples = []
        colIds.append(sampleIdx)
        samples.append(sample)
        datasetSampleIdx[dataset] = colIds
        datasetSamples[dataset] = samples
        overlap+=1
fh.close()

# open dataset outputs
if overlap == 0:
    print("Could not match samples to datasets")
    sys.exit()
else:
    print(f"{overlap} samples could be connected to a dataset")

# convert to arrays
datasets = []
for dataset in datasetSampleIdx.keys():
    datasets.append(dataset)
datasets.sort()

datasetoutputs = []
datasetindexes = []
for i in range(0,len(datasets)):
    dataset = datasets[i]
    samplesForDs = datasetSamples.get(dataset)
    if samplesForDs is None:
        print(f"No samples for dataset {dataset} but defined?")
        sys.exit()
    print(f"{i} --> {dataset} --> {len(samplesForDs)}")

    dsoutdir = outdir+"/"+dataset+"/"
    os.makedirs(dsoutdir,exist_ok=True)
    dsoutfile = dsoutdir+"/"+dataset+".counts.gz"
    ofh = getwfh(dsoutfile)
    
    
    dsheader = header[0]+" " + " ".join(samplesForDs)+"\n"
    ofh.write(dsheader)
    datasetoutputs.append( ofh )
    datasetindexes.append( datasetSampleIdx.get(dataset) )

fh = getrfh(infile)
lctr = 0
header = fh.readline().split()
for line in fh:
    elems = line.strip().split()
    id = elems[0]
    for i in range(0,len(datasets)):    
        cols = datasetindexes[i]
        ofh = datasetoutputs[i]
        ofh.write(id)
        for col in cols:
            ofh.write(" "+elems[col])
        ofh.write("\n")
    lctr += 1
    if lctr % 1000 == 0:
        print(f"{lctr} lines parsed",end='\r' )

fh.close()
print(f"{lctr} lines parsed",end='\n')

for i in range(0,len(datasets)):
    datasetoutputs[i].close()