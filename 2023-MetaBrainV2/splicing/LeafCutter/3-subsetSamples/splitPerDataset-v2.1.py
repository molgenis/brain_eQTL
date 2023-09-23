import gzip
import sys
import os
import numpy
import argparse
from enum import Enum
from pprint import pprint

#if len(sys.argv) < 5:
#	print("Usage: dataset_perind.counts.gz
#			cramToRNASeq.txt.gz
#			datasetFile.txt
#			sampleSelection.txt
#			outprefix
#			[minnrdatasets=2]
#			[lowQual(true|false)]")
#	sys.exit(-1)


# argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-c","--counts", help="Leafcutter cluster count file", required=True)
parser.add_argument("-s","--samplerename", help="Replace e.g. BAM ids to something more friendly")
parser.add_argument("-d","--datasetfile", help="Sample to dataset definition file", required=True)
parser.add_argument("-i","--sampleinclude", help="Sample include list (after --samplerename, when applied)")
parser.add_argument("-o","--out", help="Output prefix",required=True)
parser.add_argument("--minNrDatasets", help="Minimal number of datasets per junction",default=2)
parser.add_argument("--minAvgPSI", help="Minimal PSI value",default=0.01)
parser.add_argument("--averageImpute", help="Replace NaN values with average calculated over non-NaN samples",action='store_true')
parser.add_argument("--averageImputeMissingDatasets", help="Replace NaN values with average calculated over non-NaN samples, also for datasets that have not data at all, or those that didn't pass QC (like the default in LeafCutter)",action='store_true')
parser.add_argument("--averageImputePerDataset", help="Replace NaN values with average calculated over non-NaN samples per dataset",action='store_true')
parser.add_argument("--minNonNaNPerDataset",help="Minimal number of non-NaN samples per dataset",default=30)
parser.add_argument("--minEventCt",help="Minimal number counts in the junction numerator per sample (values below this number will be replaced by NaN)",default=10)


args = vars(parser.parse_args())

psifile = args["counts"]
cramfile  = args["samplerename"]
etdfile = args["datasetfile"]  # expression to dataset
samplelimitfile = args["sampleinclude"]  # expression to dataset
outprefix  = args["out"]
minnrdatasets = int(args["minNrDatasets"])

minNonNanPerDataset = int(args["minNonNaNPerDataset"])
minAvgPsi = float(args["minAvgPSI"])
maxAvgPsi = 1-(float(args["minAvgPSI"]))
minnrdatasets = int(args["minNrDatasets"])
minEventCt = int(args["minEventCt"])

imputeAverage = args["averageImpute"]
imputeAverageMissingDatasets = args["averageImputeMissingDatasets"]
imputeAveragePerDataset = args["averageImputePerDataset"]

if imputeAverage and imputeAveragePerDataset:
	print("Error: cannot average impute over all samples and per dataset at the same time!")
	sys.exit(-1)

pprint(args)

#if len(sys.argv) > 7:
#	lowq = sys.argv[7]
#	lowq = lowq.lower()
#	if lowq == "true":
#		minAvgPsi = 0
#		maxAvgPsi = 1
		
#		print("!!! Low quality mode !!!")
#	elif low != "false":
#		print("expect true or false for command line argument 7")
#		sys.exit(0)




# main functions


def calcPsi(num, denom):
	# returns list: readct, psi unfiltered, psi filtered
	try:
		num = float(num)
		denom = float(denom)
		if denom < 1:
			return [denom, numpy.nan, numpy.nan]
		elif denom < minEventCt:
			return [denom, num/denom, numpy.nan]	
		elif num < 1:
			return [denom, 0, 0]
		else:
			return [denom, num/denom, num/denom]
	except:
		return [0, numpy.nan, numpy.nan]

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
	
cramMap = None
if cramfile is not None:
	cramMap = {}
	fh = gzip.open(cramfile,'rt')
	for line in fh:
		elems = line.strip().split("\t")
		cramMap[elems[0]] = elems[1]
	fh.close()

	# custom mappings
	cramMap["NEUYV496XLP_CGND_HRA_00084-2"] = "HRA_00084"
	cramMap["AN11864_AN11864_ba41-42-22"] = "AN11864_ba41.42.22"
 
	print("{} sample mappings read from {}".format( len(cramMap), cramfile ) )

sampleLimit = None
if samplelimitfile is not None:
	sampleLimit = set()
	fh = open(samplelimitfile,'r')
	for line in fh:
		id = line.strip()
		sampleLimit.add(id)
	fh.close()

	print("{} sample ids read from {}".format(len(sampleLimit), samplelimitfile))

samplesPerDataset = None
sampleToDataset = None
selectedSamples = None

if etdfile is not None:
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
	print("{} samples selected from {}. {} datasets".format( len(selectedSamples), etdfile, len(samplesPerDataset.keys()) ) )

fh = gzip.open(psifile,'rt')
# get filename 
header = fh.readline().strip().split(" ")
outHeader = []
includedColumns = []

# parse header
ok   			  = True
colCtr  		  = 0	# column counter
relativeColCtr 		  = 0	# relative column counter
includedSamplesPerDataset = {}
relativeColIdsForDataset  = {}
matchingsamples 	  = set()

for col in header:
	if colCtr == 0:
		outHeader.append(col)
	else:
		# replace ID with actual RNA-ID
		sample = col
		if cramMap is not None:
			sample = cramMap.get(col)
			if sample is None:
				print("Warning: "+col+" was not found in sample conversion file")
				sample = col
				ok = True
		else:
			ok =  True

		# determine if column should be in dataset
		if selectedSamples is None or sample in selectedSamples:
			includedColumns.append(colCtr)
			outHeader.append(sample)

			# some dataset bookkeeping			
			# ds = samplesPerDataset.get(sample)
			ds = None
			if sampleToDataset is not None:
				ds = sampleToDataset.get(sample)
				if ds is None:
			#	print("Sample "+sample+" has no dataset annotation")
					ds = "None"
			#	ok = False
				sampleSet = includedSamplesPerDataset.get(ds)
				if sampleSet is None:
					sampleSet = set()
				sampleSet.add(sample)
				includedSamplesPerDataset[ds] = sampleSet

				relativeColIds = relativeColIdsForDataset.get(ds)
				if relativeColIds is None:
					relativeColIds = []
				relativeColIds.append(relativeColCtr)
				relativeColIdsForDataset[ds] = relativeColIds
				matchingsamples.add(sample)
				relativeColCtr += 1
	colCtr+=1

print()
if sampleLimit is not None:
	for sample in sampleLimit:
		if sample not in matchingsamples:
			print("Warning: selected sample not found: "+sample) 

print("{} samples finally included after parsing header of {}".format(len(includedColumns), psifile))
print()
print("Number of included samples per included dataset:")
for dataset in includedSamplesPerDataset.keys():
	print(dataset+"\t"+str( len( includedSamplesPerDataset.get(dataset) ) ) )

if not ok:
	print("Fix warnings above")
	sys.exit(-1)

if len(includedColumns) == 0:
	print("No samples were selected")
	sys.exit(-1)

includedDatasets = includedSamplesPerDataset.keys()

if len(includedDatasets) == 0:
	print("No datasets were defined")
	sys.exit(-1)

includedDatasetsArr = []
for ds in includedDatasets:
	includedDatasetsArr.append(ds)
includedDatasetsArr.sort()
includedDatasets = includedDatasetsArr

# write header
fhoPsiFiltered = gzip.open(outprefix+"-PSI-filtered.txt.gz",'wt')
fhoPsiUnfiltered = gzip.open(outprefix+"-PSI-unfiltered.txt.gz",'wt')
fhoCtsFiltered = gzip.open(outprefix+"-CTS-filtered.txt.gz",'wt')
fhoCtsUnfiltered = gzip.open(outprefix+"-CTS-unfiltered.txt.gz",'wt')
fhoLog = gzip.open(outprefix+"-filterlog.txt.gz",'wt')

outHeader = "\t".join(outHeader)+"\n"

fhoPsiFiltered.write(outHeader)
fhoPsiUnfiltered.write(outHeader)
fhoCtsFiltered.write(outHeader)
fhoCtsUnfiltered.write(outHeader)

logheader = "Id\tWritten\tNrDatasetsOK\tNonNan"
logheader += "\tmeanPsiUnfiltered\tvarPsiUnfiltered"
logheader += "\tmeanCtsUnfiltered\tvarCtsUnfiltered"
logheader += "\tmeanPsiFiltered\tvarPsiFiltered"
for dataset in includedDatasets:
	logheader += "\tMean-"+dataset
	logheader += "\tVar-"+dataset
	logheader += "\tN-"+dataset
	logheader += "\tLowReadCt-"+dataset
	logheader += "\tNrNonNan-"+dataset
	logheader += "\tPassQC-"+dataset
fhoLog.write(logheader+"\n")

print()
print("Processing..")
print()

# iterate input file
lineCtr    = 0
written = 0

debug = 0
debugQueryId = "chr21:16181687:16231056:clu_184_?"
debugQueryId = "NeverGonnaGiveYouUp"

for line in fh:
	elems = line.strip().split(" ")
	ctr = 0
	rowid = elems[0] # always append first column
	if rowid == debugQueryId:
		debug = 1
		print()
	valuesPsiFiltered = []
	valuesPsiUnfiltered = []
	valuesCtsStr = []
	valuesCtsUnfiltered = []
	valuesPerDataset = {}
	overallNonNan = 0

	# only parse included columns 
	for i in includedColumns:
		elem = elems[i]
		valuesCtsStr.append(elem)
		num, denom = elem.split("/")
		readCt, psiUnfilter, psiFilter = calcPsi(num,denom)
		if not numpy.isnan(psiFilter):
			overallNonNan += 1
		valuesPsiFiltered.append( psiFilter ) # PSI with minimum read ct filter applied
		valuesCtsUnfiltered.append( readCt )  # Cts without minimum read ct filter applied
		valuesPsiUnfiltered.append( psiUnfilter ) # PSI without minimum read ct filter applied

	meanPsiUnfiltered, varPsiUnfiltered = meanAndVar(valuesPsiUnfiltered)
	meanCtsUnfiltered, varCtsUnfiltered = meanAndVar(valuesCtsUnfiltered)
		
	# write unfiltered results
	fhoCtsUnfiltered.write(rowid+"\t"+"\t".join(valuesCtsStr) +"\n") # write unfiltered count input for this subset of samples
	fhoPsiUnfiltered.write(rowid+"\t"+"\t".join( [str(x) for x in valuesPsiUnfiltered] ) +"\n") # 

	# apply dataset specific filters
	okdatasets = 0
	logln = ""
	dsOkList = []
	totalNonNan = 0
	for dataset in includedDatasets:
		dsNrNonNan = 0
		dsLowReadCt = 0
		dsMeanPsi   = 0
		dsMeanCt = 0
		dsOK   = False
		dsPsiVals = []
		dsVariancePsi = 0

		relativeColsForDs = relativeColIdsForDataset.get(dataset)
		if relativeColsForDs is None:
			print("ERROR: no column selection for "+dataset)
			sys.exit(-1)

		for i in relativeColsForDs:
			v = valuesPsiFiltered[i]
			c = valuesCtsUnfiltered[i]
			if c < minEventCt:
				dsLowReadCt += 1
			dsPsiVals.append(v)
			if not numpy.isnan(v):
				dsMeanPsi += v
				dsNrNonNan  += 1
				totalNonNan += 1
				dsMeanCt  += c
			else:
				valuesCtsStr[i] = "0/0"
		
		if dsNrNonNan >= minNonNanPerDataset:
			# dsMeanPsi = 0
			if dsNrNonNan > 0:
				dsMeanPsi /= dsNrNonNan
				dsVariancePsi = var(dsPsiVals, dsMeanPsi)
			else:
				dsVariancePsi = 0
			if dsMeanPsi >= minAvgPsi and dsMeanPsi <= maxAvgPsi and dsVariancePsi > 0:
				# values for dataset pass QC
				okdatasets += 1
				dsOK = True
		else:
			dsMeanPsi = 0

		if not dsOK:
			# wipe values for dataset to NaN, but not when we're imputing missing values in bad datasets
			if not imputeAverageMissingDatasets:
				for i in relativeColsForDs:
					valuesPsiFiltered[i] = numpy.nan
					valuesCtsStr[i] = "0/0"
		elif imputeAveragePerDataset: # this requires the dataset to pass QC
			for i in relativeColsForDs:
				if numpy.isnan(valuesPsiFiltered[i]):
					valuesPsiFiltered[i] = dsMeanPsi
		dsOkList.append(dsOK)
		logln += "\t"+f(dsMeanPsi)
		logln += "\t"+f(dsVariancePsi)
		logln += "\t"+str(len(dsPsiVals))
		logln += "\t"+str(dsLowReadCt)
		logln += "\t"+str(dsNrNonNan)
		logln += "\t"+str(dsOK)
		if debug == 1:
			print("{}\t{}\t{}\t{}\t{}\t{}".format(dataset, dsNrNonNan,dsLowReadCt,dsMeanPsi,dsVariancePsi,dsMeanCt))
	#	print(dataset+"\t"+str(nonNan)+"\t"+str(dsok))
		
	meanPsiFiltered, varPsiFiltered = meanAndVar(valuesPsiFiltered)
	#print(rowid+"\t"+str(meanPsiFiltered)+"\t"+str(totalNonNan))

	#print("Preimpute: ")
	#print(rowid+"\t"+"\t".join([str(x) for x in valuesPsiFiltered]) + "\n")
	#print("------")
	# replace missing values with overall average, only for datasets that passed QC
	if imputeAverage:
		for d in range(0,len(includedDatasets)):
			if dsOkList[d] or imputeAverageMissingDatasets:
				dataset = includedDatasets[d]
				relativeColsForDs = relativeColIdsForDataset.get(dataset)
				for i in relativeColsForDs:
					if numpy.isnan(valuesPsiFiltered[i]):
						valuesPsiFiltered[i] = meanPsiFiltered
					#	print(rowid+"\t"+dataset+"\t"+str(valuesPsiFiltered[i]))


	# write if result passes filters
	eventWritten = False
	if okdatasets >= minnrdatasets:
		fhoPsiFiltered.write(rowid+"\t"+"\t".join([str(x) for x in valuesPsiFiltered]) + "\n")
	#	print("Postimpute")
	#	print("----------")
	#	print(rowid+"\t"+"\t".join([str(x) for x in valuesPsiFiltered]) + "\n")
		fhoCtsFiltered.write(rowid+"\t"+"\t".join([str(x) for x in valuesCtsStr]) + "\n")
		written += 1
		eventWritten = True

	#sys.exit(-1)
	loglnStart =  rowid  +"\t"+str(eventWritten)   +"\t"+ str(okdatasets) +"\t"+ str(overallNonNan)
	loglnStart += "\t"+f(meanPsiUnfiltered)  +"\t"+ f(varPsiUnfiltered)			
	loglnStart += "\t"+f(meanCtsUnfiltered)  +"\t"+ f(varCtsUnfiltered)			
	loglnStart += "\t"+f(meanPsiFiltered)    +"\t"+ f(varPsiFiltered)			
	fhoLog.write(loglnStart+"\t"+logln+"\n")

	lineCtr += 1

	if lineCtr % 100 == 0:
		perc = (written/lineCtr) * 100
		print("{} lines parsed, {} written - {}%.".format(lineCtr, written, f(perc)), end='\r', flush=True)
		fhoLog.flush()		
	if debug == 1:
		exit()
print("{} lines parsed, {} written - {}%.".format(lineCtr, written, f(perc)), end='\n')


fhoLog.close()		
fhoPsiFiltered.close()
fhoPsiUnfiltered.close()
fhoCtsUnfiltered.close()
fhoCtsFiltered.close()
