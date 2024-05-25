import sys
import numpy as np
import gzip
import os
import math

import sys
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from numpy.linalg import LinAlgError


matplotlib.use('Agg')



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
	variance = var(tmp, mean)
	stdev = math.sqrt(variance)
	for i in range(0,len(vals)):
		v = (tmp[i] - mean)/stdev
		vals[i] = str(v)
	return vals

# concatenate and z-transform
#for dataset in datasets:
def concatenateAndZtransformDataset(indir, tmpoutfile, dataset):
	print("Concatenating and Z-transforming dataset: "+indir+"/"+dataset)
	fho = gzip.open(tmpoutfile,'wt')
	print(f"mergeqq: {mergeqq}")
	for chr in range(1,23):
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
		if chr == 1:
			header = fh.readline().strip().split()
			header = "\t".join(header[4:len(header)])
			fho.write("SpliceEvent\t"+header+"\n")
		else:
			fh.readline()
		lctr = 0
		for line in fh:
			elems = line.strip().split()
			id = elems[3]
			vals = elems[4:len(elems)]
			vals = ztransformvals(vals)
			fho.write(id+"\t"+"\t".join(vals)+"\n")
			lctr += 1
			if lctr % 1000 == 0:
				print(f"{lctr} lines parsed",end='\r')
		print(f"{lctr} lines parsed",end='\n')
	fho.close()
	return tmpoutfile

# 
def pca(file, nrcomponents, outprefix):
	print("PCA")
	print("parsing "+file)
	df = pd.read_csv(file, sep='\t',index_col=0)
	print("Read matrix: {} x {}".format(df.shape[0], df.shape[1]))
	print(df)
	#df.set_index('-')
	print("Correlating..")
	cormat = df.corr()
	cormat = cormat.fillna(0)

	print(cormat)

	print("Performing decomposition.")
	pca = PCA(n_components=nrcomponents)
	pca.fit(cormat)
	components = pca.components_
	explVar = pca.explained_variance_ratio_
	#			print(explVar)
	pcadf = pd.DataFrame(components)
	pcadf.columns = cormat.columns
	#			pcadf.index.name = 'Component'
	pcadf = pcadf.transpose()
	colnames = []
	for comp in range(1,nrcomponents+1):
		colnames.append("PC"+str(comp))
	pcadf.columns = colnames
	#			print(pcadf)

	#pcadfmelt = pcadf.melt()
	#			print(pcadfmelt)
	pcadf.to_csv(outprefix+"_PCs.txt", sep='\t')
	return outprefix+"_PCs.txt"
	# fig, ax = plt.subplots()
	# sns.scatterplot(data=pcadf, x="PC1", y="PC2")
	# fig.savefig(outprefix+"_PC1and2.png")
	# plt.close(fig)

	# sns.scatterplot(data=pcadf, x="PC1", y="PC4")
	# fig.savefig(outprefix+"_PC1and4.png")
	# plt.close(fig)

# def checkColinearity(covariates):
#     # check for variance
#     toDrop = []
#     for column in covariates:
#         val = covariates[column].var()
#         if val == 0:
#             toDrop.append(column)
#     tmp = covariates.drop(toDrop, axis=1)
#     # check for correlation
#     cor = tmp.corr()
#     cols = tmp.columns
#     toDrop = []
#     for i in range(0, len(cols)):
#         for j in range(i+1, len(cols)):
#             colB = cols[j]
#             val = cor.iat[i, j]
#             if val > 0.9999:
#                 toDrop.append(colB)
#     tmp = tmp.drop(toDrop, axis=1)

#     removed = 1
#     i = 0
#     while removed > 0:
#         removed = 0
#         for covar in tmp:
#             y_cov = tmp[covar]
#             x_cov = tmp.drop(covar, axis=1)
#             res = sm.OLS(y_cov.astype(float), x_cov.astype(float)).fit()
#             rSquared = res.rsquared
#             if rSquared > 0.999:
#                 # drop covar
#                 tmp = tmp.drop(covar, axis=1)
#                 removed = 1
#                 break
#         i = i + 1
#     return tmp

def regressPCs(spliceFile, pcaFile, nPCs, outfileResiduals):
	print("Regressing PCs...: "+spliceFile)
	dfSplice = pd.read_csv(spliceFile, sep='\t', index_col=0)
	dfSplice = dfSplice.transpose() # ??
	dfSplice.index.name = 'Sample'
	print("splice size:", dfSplice.shape[0])

	# Load pca matrix, samples on rows, nPCs on cols
	useCols = [i for i in range(nPCs + 1)]
	dfPCA = pd.read_csv(pcaFile, sep='\t', index_col=0, usecols=useCols)
	print("pca size:", dfPCA.shape[0])
	dfPCA.index.name = 'Sample'
	dfPCA = dfPCA.reindex(dfSplice.index)
	dfCov = dfPCA

	sampleMap = {}
	sampleList = []
	sampleCount = 0
	for sample in dfSplice.index:
		sampleMap[sample] = sampleCount
		sampleList.append(sample)
		sampleCount = sampleCount + 1

	# iterate splice events
	fCount = 0
	lCount = 0
	eCount = 0
	wCount = 0

	outfhResiduals = gzip.open(outfileResiduals, 'wt')
	outfhResiduals.write("-\t" + "\t".join(sampleList) + "\n")
	for event in dfSplice:
		psi = dfSplice[[event]]
		nonna = psi.dropna(how='any')
		if nonna.empty:
			continue
		shape = nonna.shape
		outputResiduals = [np.nan] * len(psi)

		logitY = nonna[event]
		logitY = np.array(logitY)
		covars = dfCov.reindex(nonna.index)
		# covars = checkColinearity(covars)
		X = covars.to_numpy()
		X = sm.add_constant(X)

		try:
			ols = sm.OLS(logitY, X)
			olsResults = ols.fit(disp=0)
			residuals = olsResults.resid
			idx = nonna.index
			for i in range(0, len(residuals)):
				actidx = sampleMap.get(idx[i])
				outputResiduals[actidx] = residuals[i]
			strout = "\t".join(str(v) for v in outputResiduals)
			outln = event + "\t"+strout+"\n"
			outfhResiduals.write(outln)
			wCount = wCount + 1
		except PerfectSeparationError:
			fCount = fCount + 1
		except LinAlgError:
			lCount = lCount + 1

		eCount = eCount + 1
		if eCount % 50 == 0:
			print("{} lines processed, {} failed, {} linalg errors, {} written".format(eCount, fCount, lCount, wCount))

	outfhResiduals.close()
	print("Done")

def mergeAllDatasetsPerChromosome():
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

def concatenateAllMergedChromosomes():
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

def splitByChr(inputfile, outputprefix):
	chrwriters = []
	fh = gzip.open(inputfile,'rt')
	header = fh.readline()
	for chr in range(1,23):
		outputfile = outputprefix.replace("CHR",str(chr))
		handle = gzip.open(outputfile,'wt')
		handle.write(header)
		chrwriters.append(handle)
	
	for line in fh:
		elems = line.split("\t",2)
		id = elems[1]
		chr = int(id.split(":")) # assume ID is split by 
		chrwriters[chr].write(line)

	for chr in range(1,23):
		chrwriters[chr].close()


datasetstmp = []
datasetstmp.append(datasets[0])
datasets = datasetstmp

for dataset in datasets:
	# concatenate and z-transform
	ztransformedoutput = outputprefix+"/"+dataset+".phen.ztransform.gz"
	concatenateAndZtransformDataset(indir, ztransformedoutput, dataset)
	# determine 15 PCs (in stead of PEER factors)
	pcaoutput = outputprefix+"/"+dataset+".phen.pca"
	pcaoutput = pca(ztransformedoutput, 15, pcaoutput)
	# correct for components
	residualoutput = outputprefix+"/"+dataset+".phen.pca.15PCsRemoved.txt.gz"
	regressPCs(ztransformedoutput,pcaoutput,15,residualoutput)
	# split back into one file per chromosome
	splitoutputprefix = outputprefix+"/"+dataset+".phen.pca.15PCsRemoved-chrCHR.txt.gz"
	splitByChr(residualoutput,splitoutputprefix)

mergeAllDatasetsPerChromosome()
concatenateAllMergedChromosomes()
