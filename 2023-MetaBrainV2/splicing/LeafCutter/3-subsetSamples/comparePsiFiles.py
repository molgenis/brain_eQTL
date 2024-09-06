import gzip
import sys
import numpy
import random

if len(sys.argv) < 4:
    print("Usage: file1 file2 outdir")
    sys.exit()


f1 = sys.argv[1]
f2 = sys.argv[2]
outdir = sys.argv[3]


def getIDs(file):
    fh = gzip.open(file,'rt')
    junc1 = set()
    lctr = 0
    for line in fh:
        id = line.split("\t",2)[0]
        junc1.add(id)
        lctr += 1
        if lctr % 1000 == 0:
            print(f"{lctr} lctr",end='\r')
    print(f"{lctr} lctr",end='\n')
    fh.close()
    return junc1

def compareJuncs(junc1,junc2,output):
    overlap = []
    fho = open(output,'w')
    for j in junc1:
        if j in junc2:
            overlap.append(j)
            exampleOverlap = j
        else:
            fho.write(j+"\n")
    fho.close()
    print(f"{len(overlap)} overlap")
    return overlap

def getData(file, overlap):
    query = set()
    query.update(overlap)
    fh = gzip.open(file,'rt')
    print("Parsing: "+file)
    samples = {}
    header = fh.readline().strip().split("\t")
    for i in range(1,len(header)):
        samples[header[i]] = i - 1
    data = {}
    lctr = 0
    for line in fh:
        elems = line.strip().split("\t",2)
        id = elems[0]
        if id in query:
            elems = line.strip().split("\t")
            vals = []
            for i in range(1,len(elems)):
                try:
                    vals.append(float(elems[i]))
                except:
                    vals.append(numpy.nan)
            data[id] = vals
        lctr += 1
        if lctr % 1000 == 0:
            print(f"{lctr} lctr",end='\r')
    print(f"{lctr} lctr",end='\n')
    fh.close()
    print(f"Loaded: {len(data)} from: {file}")
    return samples, data

junc1 = getIDs(f1)
junc2 = getIDs(f2)

print(f"{len(junc1)} ids in 1")
print(f"{len(junc2)} ids in 2")
compareJuncs(junc1,junc2,outdir+"/junc1-vs-junc2-missing.txt")
overlap = compareJuncs(junc2,junc1,outdir+"/junc2-vs-junc1-missing.txt")

# randomly select 10 percent
nrJuncsToSample = int(len(overlap) * 0.1)
sampledJuncs = random.sample(overlap, nrJuncsToSample)
samples1, data1 = getData(f1, sampledJuncs)
print(f"{len(samples1)} samples, {len(data1)} juncs ")
samples2, data2 = getData(f2, sampledJuncs)
print(f"{len(samples2)} samples, {len(data2)} juncs ")

overlapsamples = 0
for sample in samples1:
    id2 = samples2.get(sample)
    if id2 is not None:
        overlapsamples+=1
print(f"{overlapsamples} overlapping samples")

fho = open(outdir+"/comparison-juncs.txt",'w')
fho.write("Junc\tEqual\tDiff\tNrComp\n")
wobble = 1e-6
lctr = 0
for junction in data1.keys():
    vals1 = data1.get(junction)
    vals2 = data2.get(junction)
    equal = 0
    diff = 0
    comp = 0
    for sample in samples1.keys():
        id1 = samples1.get(sample)
        id2 = samples2.get(sample)
        if id1 is not None and id2 is not None:
            comp += 1
            if id2 > len(vals2) -1 or id1 > len(vals2) -1:
                print(f"Junction has some error in f1: {junction} expected {id1} found {len(vals1)}")
                print(f"Junction has some error in f1: {junction} expected {id2} found {len(vals2)}")
                sys.exit()
            v1 = vals1[id1]
            v2 = vals2[id2]
            
            if numpy.isnan(v1) and numpy.isnan(v2):
                equal+=1
            elif (numpy.isnan(v1) and not numpy.isnan(v2)) or (numpy.isnan(v2) and not numpy.isnan(v1)):
                diff +=1 
            else:
                delta = abs(v1 - v2)
                if delta < wobble:
                    equal+=1
                else:
                    diff +=1
    fho.write(f"{junction}\t{equal}\t{diff}\t{comp}\n")
    lctr += 1
    if lctr % 1000 == 0:
        print(f"{lctr} juncs compared",end='\r')
print(f"{lctr} juncs compared",end='\n')
fho.close()