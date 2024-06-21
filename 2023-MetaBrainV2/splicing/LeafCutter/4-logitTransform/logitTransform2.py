import sys
import gzip
import numpy
import math

if len(sys.argv) < 3:
    print("Usage: psi-file.txt.gz outfile-logit.txt.gz")
    sys.exit(0)

psiFile = sys.argv[1]
outfileLogit = sys.argv[2]

def logit(v):
    if numpy.isnan(v):
        return numpy.nan
    if v < 1e-16:
        v = 1e-16
    if v > (1 - 1e-16):
        v = 1 - 1e-16
    v = math.log(v / (1-v))
    return v

fh = gzip.open(psiFile,'rt')
fho = gzip.open(outfileLogit,'wt',5)
header = fh.readline().split("\t")
header[0] = "-"
fho.write("\t".join(header))
lctr = 0
for line in fh:
	elems = line.strip().split("\t")
	outln = elems[0]
	ectr = 0
	for elem in elems:
		if ectr > 0:
			try:
				v = float(elem)
				v = logit(v)
				outln +="\t"+str(v)
			except:
				# print("error parsing: "+elem)
				outln+="\t"+str(numpy.nan)
		ectr+=1
	fho.write(outln+"\n")
	lctr += 1
	if lctr % 1000 == 0:
		print("{} lines parsed ".format(lctr), end='\r', flush=True)
print("{} lines parsed - done".format(lctr), end='\n')

fho.close()
fh.close()

