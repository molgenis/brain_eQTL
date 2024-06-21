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
parser.add_argument("-d","--datasetfile", help="Sample to dataset definition file", required=False)
parser.add_argument("-i","--sampleinclude", help="Sample include list (after --samplerename, when applied)")
parser.add_argument("-o","--out", help="Output prefix",required=True)
parser.add_argument("--minNrDatasets", help="Minimal number of datasets per junction",default=2)
parser.add_argument("--minAvgPSI", help="Minimal PSI value",default=0.01)
parser.add_argument("--minStDevPSI", help="Minimal PSI stdev value",default=0.005)
parser.add_argument("--averageImpute", help="Replace NaN values with average calculated over non-NaN samples",action='store_true')
parser.add_argument("--averageImputeMissingDatasets", help="Replace NaN values with average calculated over non-NaN samples, also for datasets that have not data at all, or those that didn't pass QC (like the default in LeafCutter)",action='store_true')
parser.add_argument("--averageImputePerDataset", help="Replace NaN values with average calculated over non-NaN samples per dataset",action='store_true')
parser.add_argument("--minNonNaNPerDataset",help="Minimal number of non-NaN samples per dataset (if specified < 1, will be interpreted as proportion of total)",default=30)
parser.add_argument("--minEventCt",help="Minimal number counts in the junction numerator per sample (values below this number will be replaced by NaN)",default=10)
parser.add_argument("--writeExtraFiles",help="Write additional output matrices (raw unfiltered PSI values etc)",action='store_true')
parser.add_argument("--removeNonStandardChr",help="Remove non-autosomal and non-X or non-Y splice events",action='store_true')
parser.add_argument("--removeNonAutosomal",help="Remove non-autosomal splice events",action='store_true')
parser.add_argument("--usePseudoCount",help="Use pseudocount in PSI calculation",action='store_true')

args = vars(parser.parse_args())

psifile = args["counts"]
cramfile  = args["samplerename"]
etdfile = args["datasetfile"]  # expression to dataset
samplelimitfile = args["sampleinclude"]  # expression to dataset
outprefix  = args["out"]
minnrdatasets = int(args["minNrDatasets"])

minNonNanPerDataset = float(args["minNonNaNPerDataset"])
minAvgPsi = float(args["minAvgPSI"])
minStDevPsi = float(args["minStDevPSI"])
maxAvgPsi = 1-(float(args["minAvgPSI"]))
minnrdatasets = int(args["minNrDatasets"])
minEventCt = int(args["minEventCt"])

imputeAverage = args["averageImpute"]
imputeAverageMissingDatasets = args["averageImputeMissingDatasets"]
imputeAveragePerDataset = args["averageImputePerDataset"]
writeExtraFiles = args["writeExtraFiles"]
removeNonAutosomal = args["removeNonAutosomal"]
removeNonStandardChr = args["removeNonStandardChr"]
usePseudoCount = args["usePseudoCount"]

if imputeAverage and imputeAveragePerDataset:
	print("Error: cannot average impute over all samples and per dataset at the same time!")
	sys.exit(-1)

pprint(args)

# main functions
def calcPsi(num, denom):
	# returns list: readct, psi unfiltered, psi filtered
	try:
		num = float(num)
		denom = float(denom)
		psi = 0
		if denom < 1:
			return [denom, numpy.nan, numpy.nan]
		elif denom < minEventCt:
			if usePseudoCount:
				psi = (num+0.5)/(denom+0.5)
			else:
				psi = num/denom
			return [denom, psi, numpy.nan]	
		elif num < 1:
			return [denom, 0, 0]
		else:
			if usePseudoCount:
				psi = (num+0.5)/(denom+0.5)
			else:
				psi = num/denom
			return [denom, psi, psi]
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

# link samples to datasets
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

# parse the header of the leafcutter count file
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
			ds = "None"
			if sampleToDataset is not None:
				ds = sampleToDataset.get(sample)
				if ds is None:
					print("Warning: sample "+sample+" has no dataset annotation")
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
fhoPsiFiltered = gzip.open(outprefix+"-PSI-filtered.txt.gz",'wt',3)
fhoPsiUnfiltered = None
fhoCtsFiltered = None
fhoCtsUnfiltered = None
fhoLog = None

outHeader = "\t".join(outHeader)+"\n"
fhoPsiFiltered.write(outHeader)

if writeExtraFiles:
	fhoPsiUnfiltered = gzip.open(outprefix+"-PSI-unfiltered.txt.gz",'wt',5)
	fhoCtsFiltered = gzip.open(outprefix+"-CTS-filtered.txt.gz",'wt',5)
	fhoCtsUnfiltered = gzip.open(outprefix+"-CTS-unfiltered.txt.gz",'wt',5)
	fhoLog = gzip.open(outprefix+"-filterlog.txt.gz",'wt',5)

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
lineCtr = 0
written = 0

debug = 0
debugQueryId = "chr21:16181687:16231056:clu_184_?"
debugQueryId = "NeverGonnaGiveYouUp"

valuesPsiFiltered = [0] * len(includedColumns)
valuesPsiUnfiltered = [0] * len(includedColumns)
valuesCtsStr = [0] * len(includedColumns)
valuesCtsUnfiltered = [0] * len(includedColumns)

for line in fh:
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
	skip = False
	
	if removeNonAutosomal and not autosomal:
		skip = True
	if removeNonStandardChr and not (autosomal or sexchr):
		skip = True

	if rowid == debugQueryId:
		debug = 1
		print()
	
	if not skip:

		valuesPerDataset = {}
		overallNonNan = 0

		# only parse included columns 
		colctr = 0
		for i in includedColumns:
			elem = elems[i]
			valuesCtsStr[colctr] = elem
			num, denom = elem.split("/")
			readCt, psiUnfilter, psiFilter = calcPsi(num,denom)
			if not numpy.isnan(psiFilter):
				overallNonNan += 1
			valuesPsiFiltered[colctr] = psiFilter  # PSI with minimum read ct filter applied
			valuesCtsUnfiltered[colctr] = readCt   # Cts without minimum read ct filter applied
			valuesPsiUnfiltered[colctr] = psiUnfilter # PSI without minimum read ct filter applied
			colctr += 1
			
		# write unfiltered results
		if writeExtraFiles:
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
			dsNonNanThreshold = minNonNanPerDataset
			if minNonNanPerDataset < 1:
				dsNonNanThreshold = minNonNanPerDataset * len(dsPsiVals)
			
			if dsNrNonNan >= dsNonNanThreshold:
				# dsMeanPsi = 0
				if dsNrNonNan > 1:
					dsMeanPsi /= dsNrNonNan
					dsVariancePsi = var(dsPsiVals, dsMeanPsi)
				else:
					dsVariancePsi = 0

				# ala leafcutter: impute then calculate stdev
				if imputeAveragePerDataset:
					dsPsiVals = []
					for i in relativeColsForDs:
						if numpy.isnan(valuesPsiFiltered[i]):
							valuesPsiFiltered[i] = dsMeanPsi
						dsPsiVals.append(valuesPsiFiltered[i])
				if dsMeanPsi >= minAvgPsi and dsMeanPsi <= maxAvgPsi and dsVariancePsi > minStDevPsi:
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
			# elif imputeAveragePerDataset: # this requires the dataset to pass QC
			# 	for i in relativeColsForDs:
			# 		if numpy.isnan(valuesPsiFiltered[i]):
			# 			valuesPsiFiltered[i] = dsMeanPsi
			dsOkList.append(dsOK)
			if writeExtraFiles:
				logln += "\t"+f(dsMeanPsi)
				logln += "\t"+f(dsVariancePsi)
				logln += "\t"+str(len(dsPsiVals))
				logln += "\t"+str(dsLowReadCt)
				logln += "\t"+str(dsNrNonNan)
				logln += "\t"+str(dsOK)
			if debug == 1:
				print("{}\t{}\t{}\t{}\t{}\t{}".format(dataset, dsNrNonNan,dsLowReadCt,dsMeanPsi,dsVariancePsi,dsMeanCt))
		#	print(dataset+"\t"+str(nonNan)+"\t"+str(dsok))
			
		
		#print(rowid+"\t"+str(meanPsiFiltered)+"\t"+str(totalNonNan))

		#print("Preimpute: ")
		#print(rowid+"\t"+"\t".join([str(x) for x in valuesPsiFiltered]) + "\n")
		#print("------")
		# replace missing values with overall average, only for datasets that passed QC
		if imputeAverage:
			meanPsiFiltered, varPsiFiltered = meanAndVar(valuesPsiFiltered)
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
			if writeExtraFiles:
				fhoCtsFiltered.write(rowid+"\t"+"\t".join([str(x) for x in valuesCtsStr]) + "\n")
			written += 1
			eventWritten = True

		#sys.exit(-1)
		if writeExtraFiles:
			meanPsiUnfiltered, varPsiUnfiltered = meanAndVar(valuesPsiUnfiltered)
			meanCtsUnfiltered, varCtsUnfiltered = meanAndVar(valuesCtsUnfiltered)
			meanPsiFiltered, varPsiFiltered = meanAndVar(valuesPsiFiltered)
			loglnStart =  rowid  +"\t"+str(eventWritten)   +"\t"+ str(okdatasets) +"\t"+ str(overallNonNan)
			loglnStart += "\t"+f(meanPsiUnfiltered)  +"\t"+ f(varPsiUnfiltered)			
			loglnStart += "\t"+f(meanCtsUnfiltered)  +"\t"+ f(varCtsUnfiltered)			
			loglnStart += "\t"+f(meanPsiFiltered)    +"\t"+ f(varPsiFiltered)			
			fhoLog.write(loglnStart+"\t"+logln+"\n")

	lineCtr += 1

	if lineCtr % 100 == 0:
		perc = (written/lineCtr) * 100
		print("{} lines parsed, {} written - {}%.".format(lineCtr, written, f(perc)), end='\r', flush=True)
		if writeExtraFiles:
			fhoLog.flush()		
	if debug == 1:
		exit()
perc = (written/lineCtr) * 100
print("{} lines parsed, {} written - {}%.".format(lineCtr, written, f(perc)), end='\n')



fhoPsiFiltered.close()

if writeExtraFiles:
	fhoLog.close()		
	fhoPsiUnfiltered.close()
	fhoCtsUnfiltered.close()
	fhoCtsFiltered.close()
