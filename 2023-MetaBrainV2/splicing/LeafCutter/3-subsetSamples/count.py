import gzip
import sys

if len(sys.argv) < 4:
	print("Usage: tissue annot file")
	sys.exit()

tissue = sys.argv[1]
annot = sys.argv[2]
file = sys.argv[3]

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
clusters = set()
genes = set()
for line in fh: 
	# form:  chr21:9068629:9069445:clu_99_?
	elems = line.strip().split("\t",2)
	id = elems[0]
	cluster = id.split(":")[-1]
	ids.add(id)
	clusters.add(cluster)
	gene = idToGene.get(id)
	if gene is not None:
		genes.add(gene)
fh.close()

print(f"{tissue}\t{len(ids)} junctions, {len(clusters)} clusters, {len(genes)} genes")
