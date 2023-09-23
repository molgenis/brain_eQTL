import sys
import os


tissues = ["Basalganglia-EUR", "Cerebellum-EUR","Cortex-AFR", "Cortex-EAS", "Cortex-EUR","Hippocampus-EUR","Spinalcord-EUR"]
tissues = ["Cortex-EUR", "Cortex-AFR"]
minReadCt = [1,10]
meanImp = [ True, False ]
dsMeanImp = [ True, False ]
avgpsi = [0,0.01]
minNonNaNPerDataset = [1,10,30]

tmpdir = "tmp02"
basedir = "/groups/umcg-biogen/"+tmpdir+"/output/2021-FreezeThree/2021-02-18-splicing/2023-02-11-all-samples-leafcutter-fractionfix/"

outdir = basedir+"/3-subsetSamples/outputLowQual-run2/"


def writeJob(cmd,tissue, settingsStr):
	cmd += "-o "+outdir+"/"+tissue+"-"+settingsStr
	cmdStr = cmd.replace(" -"," \\\n\t-")
						
	# header
	header = "#!/bin/bash\n"
	header += "#SBATCH --ntasks=1\n"
	header += "#SBATCH --time=5:59:00\n"
	header += "#SBATCH --mem=2g\n"
	header += "#SBATCH --cpus-per-task=1\n"
	header += "#SBATCH -o "+outdir+"/"+tissue+"-"+settingsStr+".log\n"
	header += "#SBATCH -e "+outdir+"/"+tissue+"-"+settingsStr+".err\n"
	header += "#SBATCH -J "+tissue+"-"+settingsStr+"\n\n"
				
	cmdStr = header + cmdStr
						

	print(cmdStr)
	fho = open(outdir+"/"+tissue+"-"+settingsStr+".sh",'w')
	fho.write(cmdStr)
	fho.close()
#						os.system(cmd)
for tissue in tissues:
	minds = 2
	if tissue == "Cortex-EAS":
		minds = 1
	for ct in minReadCt:
		for imp in meanImp:
			if imp:
				# emulate standard LeafCutter approach
				# do not bother about dataset specific thresholds, such as minNonNaNPerDataset, avgpsi
				# build command
				minds = 1
				cmd = "ml Python/3.10.4-GCCcore-11.3.0\n\n"
				cmd += "python "+basedir+"/3-subsetSamples/splitPerDataset-v2.1.py "
				cmd += "-c "+basedir+"/1-createClusters/output/meta-brain_perind.counts.gz "
				cmd += "-s "+basedir+"/3-subsetSamples/cramToRNAseqId.txt.gz "
				cmd += "-d /groups/umcg-biogen/"+tmpdir+"/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-"+tissue+"-SampleToDataset.txt "
				cmd += "-i /groups/umcg-biogen/"+tmpdir+"/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-"+tissue+"-SelectedRNASeqSamples.txt "
				cmd += "--minNrDatasets "+str(minds) +" "
				cmd += "--minEventCt "+str(ct) +" "
				
				

				# imputation part
				cmd +="--averageImpute "
				cmd +="--averageImputeMissingDatasets " # this is how the default leafcutter pipeline operates... this will impute values for datasets that don't pass the QC thresholds, which were are otherwise replace with NaN
				settingsStr = "minct"+str(ct)
				settingsStr=settingsStr+"-meanImp"
				writeJob(cmd,tissue,settingsStr)
			else:
				for dsimp in dsMeanImp:
					for minNonNan in minNonNaNPerDataset:
						for avg in avgpsi:
							# build command
							cmd = "ml Python/3.10.4-GCCcore-11.3.0\n\n"
							cmd += "python "+basedir+"/3-subsetSamples/splitPerDataset-v2.1.py "
							cmd += "-c "+basedir+"/1-createClusters/output/meta-brain_perind.counts.gz "
							cmd += "-s "+basedir+"/3-subsetSamples/cramToRNAseqId.txt.gz "
							cmd += "-d /groups/umcg-biogen/"+tmpdir+"/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-"+tissue+"-SampleToDataset.txt "
							cmd += "-i /groups/umcg-biogen/"+tmpdir+"/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-"+tissue+"-SelectedRNASeqSamples.txt "
						
							cmd += "--minNrDatasets "+str(minds) +" "
							cmd += "--minAvgPSI "+str(avg) +" "

							cmd += "--minEventCt "+str(ct) +" "
							cmd += "--minNonNaNPerDataset "+str(minNonNan) + " "

							settingsStr = "minct"+str(ct)+"-avgpsi"+str(avg)+"-minNonNan"+str(minNonNan)
							if dsimp:
								cmd +="--averageImputePerDataset "
								settingsStr+="-dsMeanImp"
							writeJob(cmd,tissue,settingsStr)
							
# for ct in minReadCt:
# 	for imp in meanImp:
# 		for dsimp in dsMeanImp:
# 			if not (imp and dsimp):
				
# 					for minNonNan in minNonNaNPerDataset:
# 						for avg in avgpsi:
# 							minds = 2
# 							if tissue == "Cortex-EAS":
# 								minds = 1

# 							# build command
							
# 							cmd = "ml Python/3.10.4-GCCcore-11.3.0\n\n"
# 							cmd += "python "+basedir+"/3-subsetSamples/splitPerDataset-v2.1.py "
# 							cmd += "-c "+basedir+"/1-createClusters/output/meta-brain_perind.counts.gz "
# 							cmd += "-s "+basedir+"/3-subsetSamples/cramToRNAseqId.txt.gz "
# 							cmd += "-d /groups/umcg-biogen/"+tmpdir+"/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-"+tissue+"-SampleToDataset.txt "
# 							cmd += "-i /groups/umcg-biogen/"+tmpdir+"/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-"+tissue+"-SelectedRNASeqSamples.txt "
						
# 							cmd += "--minNrDatasets "+str(minds) +" "
# 							cmd += "--minAvgPSI "+str(avg) +" "

# 							cmd += "--minEventCt "+str(ct) +" "
# 							cmd += "--minNonNaNPerDataset "+str(minNonNan) + " "

# 							settingsStr = "minct"+str(ct)+"-avgpsi"+str(avg)+"-minNonNan"+str(minNonNan)

# 							if dsimp:
# 								cmd +="--averageImputePerDataset "
# 								settingsStr+="-dsMeanImp"
# 								writeJob(cmd,tissue,settingsStr)
# 							elif imp:
# 								cmd +="--averageImpute "
# 								cmd +="--averageImputeMissingDatasets " # this is how the default leafcutter pipeline operates... this will impute values for datasets that don't pass the QC thresholds, which were are otherwise replace with NaN
# 								settingsTmp=settingsStr+"-meanImp"
# 								writeJob(cmd,tissue,settingsTmp)
# 								# also average impute missing datasets
# 								#settingsTmp=settingsStr+"-meanImpBadDs"
# 								#cmd +="--averageImputeMissingDatasets "
# 								#writeJob(cmd,tissue,settingsTmp)
# 							else:
# 								writeJob(cmd,tissue,settingsStr)
											
						

# #						sys.exit(-1)
