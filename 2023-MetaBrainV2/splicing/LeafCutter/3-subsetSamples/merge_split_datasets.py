import sys
import numpy as np
import gzip
import os
import math


if len(sys.argv) < 5:
	print("Indir should have gzipped leafcutter clustering output files.")
	print("Usage: indir datasetlist.txt mergeQQ[true/false] outputprefix")
	sys.exit()



indir = sys.argv[1]
datasetlist = sys.argv[2]
mergeqq = sys.argv[3]
if mergeqq.lower() == "true" or mergeqq.lower() == "1":
	mergeqq = True
else:
	mergeqq = False

# ztransform = sys.argv[4]
# if ztransform.lower() == "true" or ztransform.lower() == "1":
# 	ztransform = True
# else:
# 	ztransform = False


outputprefix = sys.argv[4]


print("Merge splicing files over datasets.")
print(f"indir: {indir}")
print(f"datasetlist: {datasetlist}")
print(f"mergeqq: {mergeqq}")
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

def var(vals, mean):
	variance = 0.0
	for v in vals:
		variance += (v - mean) * (v-mean)
	variance /= (len(vals)-1)
	return variance

def ztransformvals(vals):
	tmp = []
	mean = 0
	for v in vals:
		f = float(v)
		mean += f
		tmp.append(f)
	mean /= len(vals)
	variance = var(vals, mean)
	stdev = math.sqrt(variance)
	for i in range(0,len(vals)):
		v = (tmp[i] - mean)/stdev
		vals[i] = str(v)
	return vals

	

for chr in range(1,23):
	print(f"{chr}")
	alldata = {}
	allids = set()
	allsamples = []
	nsamples = {}
	
	for dataset in datasets:
		file = indir+"/"+dataset+"/"+dataset+".phen_chr"+str(chr)
		if not os.path.exists(file):
			file = indir+"/"+dataset+"/"+dataset+".phen_chr"+str(chr)+".gz"
		
		if mergeqq:
			file = indir+'/'+dataset+"/"+dataset+".qqnorm_chr"+str(chr)
			if not os.path.exists(file):
				file = indir+'/'+dataset+"/"+dataset+".qqnorm_chr"+str(chr)+".gz"
		
		if not os.path.exists(file):
			print("Could not find file (or non-gzipped version of): "+file+" for dataset "+dataset)
			sys.exit(-1)

		print(file)
		fh = getfh(file)
		header = fh.readline().strip().split()
		dsnsamples = 0

		for i in range(4,len(header)):
			sample = header[i]
			allsamples.append(sample)
			dsnsamples += 1
		data = {}
		lctr = 0
		for line in fh:
			elems = line.strip().split()
			id = elems[3]
			allids.add(id)
			vals = elems[4:len(elems)]
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

	
	tmpoutfile = f"{outputprefix}-noqqnorm-chr{chr}.txt.gz"
	if mergeqq:
		tmpoutfile = f"{outputprefix}-qqnorm-chr{chr}.txt.gz"
	fho = gzip.open(tmpoutfile,'wt')
	header = "id"
	sampleset = "\t".join(allsamples)
	header = header + "\t" + sampleset
	fho.write(header+"\n")
	idctr = 0
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

tmpoutfile = f"{outputprefix}-noqqnorm.txt.gz"
if mergeqq:
	tmpoutfile = f"{outputprefix}-qqnorm.txt.gz"
print(f"concatenating into: {tmpoutfile}")
fho = gzip.open(tmpoutfile,'wt')

for chr in range(1,23):
	infile = f"{outputprefix}-noqqnorm-chr{chr}.txt.gz"
	if mergeqq:
		infile = f"{outputprefix}-qqnorm-chr{chr}.txt.gz"
	print(f"Reading: {infile}" )
	fh = getfh(infile)
	header = fh.readline()
	if chr == 1:
		fho.write(header)
	for line in fh:
		fho.write(line)
fho.close()
