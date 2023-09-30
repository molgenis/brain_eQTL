import sys
import glob
import os

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def pwarn(msg):
	print(f"{bcolors.WARNING}"+msg+f"{bcolors.ENDC}")
def pfail(msg):
	print(f"{bcolors.FAIL}"+msg+f"{bcolors.ENDC}")

if len(sys.argv) < 4:
	print("Usage: batchdir outdir jobdir")
	sys.exit(0)
batchdir = sys.argv[1]
outdir = sys.argv[2]
jobdir = sys.argv[3]

print(batchdir)

files = glob.glob(batchdir+'/*.txt')
# files = files.sort()
files.sort()
print(str(len(files))+" files found")
for file in files:
	basename = os.path.basename(file)
	outfile = outdir + basename.replace(".txt","") + "-TopEffects.txt"
	checkfile = outdir + basename.replace(".txt","") + "-TopEffects.finished"
	if not os.path.exists(checkfile):
		jobname = jobdir+basename.replace(".txt",".sh")
		if not os.path.exists(outfile):
			pwarn(jobname+"\t did not run")
		else:
			print(outfile+"\t exists")
			pfail(jobname+"\t did not finish")
		# resubmit

		cmd = "sbatch " + jobname
		print(cmd)
#		os.system(cmd)
	else:
		# open the file, check whether all genes are present
		gset = set()
		fh = open(file,'rt')
		for line in fh:
			gset.add(line.strip())
		fh.close()

		gset2 = set()
		fh = open(outfile,'rt')
		fh.readline()
		for line in fh:
			gene = line.split("\t")[0]
			gset2.add(gene)
		fh.close()
		
		# check if all genes were tested
#		if len(gset) != len(gset2):
#			print(basename+"\t did not test all genes: {} in batch {} in output".format( len(gset), len(gset2) ))
#	print()	

