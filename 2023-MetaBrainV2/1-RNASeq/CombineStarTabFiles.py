import glob
import gzip
import sys

if len(sys.argv) < 3:
    print("Usage: stardir outfile.txt.gz")
    sys.exit(-1)

parseSJs = False
# if len(sys.argv) == 5:
#     parseSJs = True
#     print("Not implemented.")
#     sys.exit(-1)

def getfh(file):
    if file.endswith(".gz"):
        return gzip.open(file,'rt')
    else:
        return open(file,'r')

def getfho(file):
    if file.endswith(".gz"):
        return gzip.open(file,'wt')
    else:
        return open(file,'w')

indir = sys.argv[1]
outfile = sys.argv[2]

query = "ReadsPerGene.out.tab"
if parseSJs:
    query="SJ.out.tab"
files = glob.glob(indir+"/*"+query)
if len(files) == 0:
    files = glob.glob(indir+"/*"+query+".gz")

nrfiles = len(files)
if files == 0:
    print("No tab files in directory: "+indir)
    sys.exit(0)

files.sort()

print(f"{nrfiles} files found.")
samples = []
genes = set()
data = {}

for file in files:
    print("Parsing: "+file)
    sample = file.split("/")[-1] # get file name
    sample = sample.replace(".gz","")
    sample = sample.replace("."+query,"")
    samples.append(sample)
    sampledata = {}
    fh = getfh(file)
    for line in fh:
        elems = line.strip().split("\t")
        gene = elems[0]
        if gene.startswith("ENSG"):
            ct = elems[1]
            genes.add(gene)
            oldct = sampledata.get(gene)
            if oldct is not None:
                print("Duplicate gene: "+gene+" in sample: "+sample)
            else:
                sampledata[gene] = ct
    fh.close()
    data[sample] = sampledata

nrsamples = len(samples)
nrgenes = len(genes)

genearr = []
for gene in genes:
    genearr.append(gene)
genes = genearr
genearr.sort()

print(f"{nrsamples} samples found.")
print(f"{nrgenes} genes found.")

print(f"Writing {outfile}")

fho = getfho(outfile)
header = "-\t"+"\t".join(samples)+"\n"
fho.write(header)
ctr = 0
for gene in genes:
    
    outct = []
    for sample in samples:
        # do stuff
        sampledata = data.get(sample)
        if sampledata is None:
            print("Sample found, but not loaded: "+sample)
            sys.exit(-1)
        else:
            genect = sampledata.get(gene)
            if genect is None:
                outct.append("nan")
            else:
                outct.append(genect)
    outln = gene+"\t"+"\t".join(outct)+"\n"
    fho.write(outln)
    ctr += 1
    if ctr % 1000 == 0:
        print(f"{ctr} lines written",end='\r')
print(f"{ctr} lines written",end='\n')
fho.close()

