import gzip
import sys

inf=sys.argv[1]
outf=sys.argv[2]

fh=gzip.open(inf,'rt')
fho=gzip.open(outf,'wt')

fho.write(fh.readline())
for line in fh:
	elems = line.strip().split("\t")
	ds = elems[11].split(";") # datasets
	sz = elems[13].split(";") # samplesizes
	ctr = 0
	nrds = 0
	samplesize = 0
	while ctr < len(ds):
		if ds[ctr] != "-":
			samplesize = samplesize + int(sz[ctr])
			nrds = nrds + 1
		ctr = ctr + 1
	elems[11] = str(nrds)
	elems[12] = "-"
	elems[13] = str(samplesize)
	elems[14] = "-"
	elems[15] = "-"
	elems[17] = "-"
	elems[19] = "-"
	elems[20] = "-"

	fho.write("\t".join(elems)+"\n")
fho.close()
fh.close()
