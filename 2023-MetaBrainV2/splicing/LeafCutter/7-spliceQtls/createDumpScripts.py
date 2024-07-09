import gzip
import sys
import os
import glob
import argparse

debug = True

def writeJob(exp, gte, gt, template, jobfile, outprefix, logprefix, chr, annotation, groups, jobname, genelist):
	print("Writing job: "+jobfile)
	fh = open(template,'r')
	lines = fh.readlines()
	fh.close()
	fho = open(jobfile,'w')
	for line in lines:
		line = line.replace("JOBNAME",jobname)
		line = line.replace("GENOTYPE",gt)
		line = line.replace("GTE",gte)
		line = line.replace("EXPRESSION",exp)
		line = line.replace("CHROM",str(chr))
		line = line.replace("OUTPREFIX",outprefix)
		line = line.replace("LOGPREFIX",logprefix)
		line = line.replace("ANNOTATION", annotation)
		line = line.replace("BATCHFILE",genelist)
		if groups is not None:
			line = line.replace("GROUPS", groups)
		fho.write(line)
	fho.close()
	
def checkDir(path, remove):
    if os.path.exists(path):
        # delete contents
        if remove:
            files = glob.glob(path+"*")
            for file in files:
                if debug:
                    print("Removing: "+file)
                os.remove(file)
    else:
        if debug:
            print("Creating dir: "+path)
        os.mkdir(path)

def getfh(file, mode):
    if file.endswith(".gz"):
        if mode == 'r':
            return gzip.open(file,'rt')
        elif mode == 'w':
            return gzip.open(file,'wt', 3)
    else:
        return open(file,mode)

def filterGeneList(genelist, groupsfile, outfile, signeffects):
    print("Filtering on significant effects from: "+signeffects)
    fh = getfh(signeffects,'r')
    header = fh.readline().split("\t")
    grouped = False
    if header[0] == "Group":
         grouped = True
    print(f"")
    significant = set()
    for line in fh:
        elems = line.strip().split("\t")
        significant.add(elems[0])
    fh.close()
    print("")

    towrite = significant
    if grouped:
        allowed = set()
        fh = getfh(genelist,'r')
        for line in fh:
            allowed.add(line.strip())
        fh.close()

        fh = getfh(groupsfile,'r')
        tmp = set()
        for line in fh:
            elems = line.strip().split() 
            id = elems[0]
            if id in allowed:
                grp = elems[1]
                if grp in significant:
                    tmp.add(id)
        fh.close()
        towrite = tmp
    fh = open(outfile,'w')
    for id in towrite:
        fh.write(id+"\n")
    fh.close()

    return outfile

parser = argparse.ArgumentParser()
parser.add_argument("--batchnameprefix", dest="batchnameprefix",
	help="Batch name prefix", required=True)

parser.add_argument("--exp", dest="expfile",
	help="Expression file", required=True)

parser.add_argument("--gte", dest="gte",
	help="GTE file", required=True)

parser.add_argument("--vcf", dest="genotype",
	help="VCF genotypes", required=True)

parser.add_argument("--genelist", dest="genelist",
	help="List of genes to include", required=True)

parser.add_argument("--annotation", dest="annotation",
	help="Annotation file", required=True)

parser.add_argument("--groups", dest="groupsfile",
	help="List defining groups of genes/features/etc")

parser.add_argument("--template", dest="template",
	help="Job template", required=True)

parser.add_argument("--nrgenes", dest="nrgenes",
	help="Nr Genes per batch", default=200)

parser.add_argument("--outdir", dest="out",
	help="Output directory", required=True)

parser.add_argument("--significanteffects", dest="significanteffects",
	help="Limit dump to significant genes or groups of genes listed in QTL file", required=False)

parser.add_argument("--debug", dest="debug",
	help="Debug?",action='store_true')

args = vars(parser.parse_args())


debug = args["debug"]

# if len(sys.argv) < 9:
#     print("Usage: createbatches.py batchnameprefix expfile.txt.gz gte.txt genotype.vcf.gz genelist.txt.gz annotation.txt.gz " +
#           "[group_annotation.txt] template.sh nrmaxgenesperbatch outdir")
#     actr = 0
#     for arg in sys.argv:
#         print(str(actr)+"\t"+arg)
#         actr += 1
#     sys.exit(0)

batchnameprefix = args["batchnameprefix"]
expfile = args["expfile"]
gte = args["gte"]
genotype = args["genotype"]
genelist = args["genelist"]
annotation = args["annotation"]
groupsfile = args["groupsfile"]
template = args["template"]
nrgenes = int(args["nrgenes"])
out = args["out"]
signeffects = args["significanteffects"]

abspath = os.path.abspath(out)
checkDir(abspath+"/output/", False)
checkDir(abspath+"/jobs/", False)
checkDir(abspath+"/logs/", False)

actualgenelist = genelist
if signeffects is not None:
    checkDir(abspath+"/genelist/", False)
    genelist = filterGeneList(genelist, groupsfile, abspath+"/genelist/genelist.txt", signeffects)

for chr in range(1,23):
    batchname = "chr"+str(chr)+"-dump"
    jobfile = abspath+"/jobs/"+batchname+".sh"
    outprefix = abspath+"/output/"+batchname
    logprefix = abspath+"/logs/"+batchname
    writeJob(expfile, gte, genotype, template, jobfile, outprefix, logprefix, chr, annotation, groupsfile, batchname, actualgenelist)