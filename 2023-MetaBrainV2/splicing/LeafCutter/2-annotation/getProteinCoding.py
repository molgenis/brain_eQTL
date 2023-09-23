import gzip
import sys

if len(sys.argv) < 4:
	print("Usage: annotfile.txt.gz outfile.txt.gz proteincodinglist.txt.gz")

#proteincodingfile = "/groups/umcg-biogen/tmp01/annotation/GeneReference/GencodeV32/gencode.v32.primary_assembly.annotation.collapsedGenes.proteincoding.txt.gz"
proteincodingfile = sys.argv[3]
proteincodingset = set()

fh = gzip.open(proteincodingfile,'rt')
for line in fh:
	proteincodingset.add(line.strip())
fh.close()

annotfile = sys.argv[1]
outfile = sys.argv[2]

print("Getting protein coding splice events")
print("Input: "+annotfile)
print("Output: "+outfile)

fho = gzip.open(outfile,'wt')
fh = gzip.open(annotfile,'rt')
fh.readline()
for line in fh:
	elems = line.strip().split("\t")
	id = elems[1]
	gene = elems[6]
	if gene in proteincodingset:
		fho.write(id+"\n")
fho.close()

