import gzip
import sys

tissue = sys.argv[1]
file = sys.argv[2]
annot = sys.argv[3]

fh = gzip.open(annot, 'rt')
fh.readline()
idToGene = {}
for line in fh:
	elems = line.strip().split("\t")
	id = elems[1]
	gene = elems[6]
	idToGene[id] = gene
fh.close()

fh = gzip.open(file,'rt')
elems = fh.readline().strip().split("\t")
ids = set()
genes = set()
for line in fh: 
	# form:  chr21:9068629:9069445:clu_99_?
	id = line.strip().split("\t",2)[0]
	ids.add(id)
	gene = idToGene.get(id)
	if gene is not None:
		genes.add(gene)
fh.close()

print("{}\t{}\t{}".format(tissue, len(ids), len(genes)))
