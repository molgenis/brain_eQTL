import sys
import numpy as np
import gzip
import os
import math


if len(sys.argv) < 5:
	print("Indir should have gzipped leafcutter clustering output files.")
	print("Usage: indir datasetlist.txt templateDATASET-CHR outputprefix")
	sys.exit()

indir = sys.argv[1]
datasetlist = sys.argv[2]
template = sys.argv[3]
outputprefix = sys.argv[4]


print("Merge splicing files over datasets.")
print(f"indir: {indir}")
print(f"datasetlist: {datasetlist}")
print(f"template: {template}")
print(f"outputprefix: {outputprefix}")


datasets = []
fh = open(datasetlist,'r')
for line in fh:
	dataset = line.strip()
	if len(dataset) > 0:
		datasets.append(dataset)
fh.close()
print(f"Nr datasets: {len(datasets)}")

def getfh(file):
    if file.endswith(".gz"):
        return gzip.open(file, 'rt')
    return open(file)



for chr in range(1,23):
	print(f"{chr}")
	alldata = {}
	allids = set()
	allsamples = []
	nsamples = {}
	
	for dataset in datasets:
		datasettemplate = template.replace("CHR",str(chr))
		datasettemplate = datasettemplate.replace("DATASET",dataset)
		
		file = indir+"/"+datasettemplate
		
		if not os.path.exists(file):
			print("Could not find file (or non-gzipped version of): "+file+" for dataset "+dataset)
			sys.exit(-1)

		print(file)
		fh = getfh(file)
		header = fh.readline().strip().split()
		dsnsamples = 0

		for i in range(1,len(header)):
			sample = header[i]
			allsamples.append(sample)
			dsnsamples += 1
		data = {}
		lctr = 0
		for line in fh:
			elems = line.strip().split()
			id = elems[0]
			allids.add(id)
			vals = elems[1:len(elems)]
			# if ztransform:
			# 	ztransformvals(vals)
			data[id] = vals
			lctr += 1
			if lctr % 1000 == 0:
				print(f"{lctr} lines",end='\r')
		print(f"{lctr} lines",end='\n')
		alldata[dataset] = data
		nsamples[dataset] = dsnsamples
		print(f"{len(alldata[dataset])} lines -- {dsnsamples} samples")
		fh.close()
	print("Merge...")
	print(f"nrSamples: {len(allsamples)}")
	print(f"splice variants: {len(allids)}")

	
	tmpoutfile = f"{outputprefix}-chr{chr}.txt.gz"

	fho = gzip.open(tmpoutfile,'wt')
	header = "id"
	sampleset = "\t".join(allsamples)
	header = header + "\t" + sampleset
	fho.write(header+"\n")
	idctr = 0
	
    # sort IDs
	allIdsList = []
	for id in allids:
		allIdsList.append(id)
	allIdsList.sort()
	allids = allIdsList
	
	for id in allids:
		outp = [np.nan] * len(allsamples)
		sctr = 0
		for dataset in datasets:
			n = nsamples.get(dataset)
			data = alldata.get(dataset)
			vals = data.get(id)
			if vals is None:
				# replace outp
				sctr += n
			else:
				# replace outp
				for val in vals:
					outp[sctr] = val
					sctr += 1
		outln = id +"\t" +"\t".join(str(x) for x in outp) +"\n"
		fho.write(outln)
		idctr += 1
		if idctr % 1000 == 0:
			print(f"{idctr} lines written",end='\r')
	print(f"{idctr} lines written",end='\n')
	fho.close()

tmpoutfile = f"{outputprefix}.txt.gz"

print(f"concatenating into: {tmpoutfile}")
fho = gzip.open(tmpoutfile,'wt')

for chr in range(1,23):
	infile = f"{outputprefix}-chr{chr}.txt.gz"

	print(f"Reading: {infile}" )
	fh = getfh(infile)
	header = fh.readline()
	if chr == 1:
		fho.write(header)
	for line in fh:
		fho.write(line)
fho.close()
