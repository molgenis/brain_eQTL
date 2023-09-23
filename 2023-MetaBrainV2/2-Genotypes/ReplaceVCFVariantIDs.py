import gzip
import sys

if len(sys.argv) < 4:
	print("Usage: reffile infile outfile")

reffile = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

class variant:

	def __init__(self, allele1, allele2, id):
		self.allele1 = allele1
		self.allele2 = allele2
		self.id = id


def countSharedAlleles(al1, al2):
	# count number of shared alleles
	shared = 0
	i = 0
	while i < 2:
		j = 0
		while j < 2:
				if al1[i] == al2[j]:
						shared = shared + 1
				j = j + 1
		i = i + 1
	return shared

def complement(allele):
	if allele == "A":
		return "T"
	if allele == "T":
		return "A"
	if allele == "C":
		return "G"
	if allele == "G":
		return "C"

def complementBoth(alleles):
	return [complement(alleles[0]),complement(alleles[1])]

def flipalleles(alleles1, assessed1, alleles2, assessed2):
	al1 = alleles1
	al2 = alleles2

	shared = countSharedAlleles(al1,al2)
	if shared < 2:
		# try complement
		al2 = complementBoth(al2)
		assessed2 = complement(assessed2)
		shared = countSharedAlleles(al1,al2)

	if shared == 2:
		if assessed1 == assessed2:
				return 0
		else:
				return 1
	else:
		return -1



positions = set()

fh = gzip.open(infile,'rt')
print("Parsing " +infile)
ctr = 0
for line in fh:
	if not line.startswith('#'):
		elems = line.split("\t",3)
		chr = elems[0]
		pos = elems[1]
		#alleles = [elems[3], elems[4]]
		positions.add(chr+":"+pos)
		ctr += 1
		if ctr % 10000 == 0:
			print("{} lines parsed".format(ctr),end='\r')
print("{} lines parsed".format(ctr),end='\n')
inlines = ctr
fh.close()

print("{} positions in {}".format(len(positions), infile))


posToId = {}
ctr = 0
matches = 0
fh = gzip.open(reffile,'rt')
print("Parsing: "+reffile)
for line in fh:
	if not line.startswith("#"):
		#CHROM  POS     ID      REF     ALT
		elems = line.split("\t", 5)
		chr = elems[0]
		pos = elems[1]
		chrpos = chr+":"+pos
		if chrpos in positions:
			allele1 = elems[3]
			allele2 = elems[4]             # what about multi allelics?
		#     id = chr+":"+pos+":"+elems[2]+":"+allele1+"_"+allele2
			id = elems[2]
			var = variant(allele1, allele2, id)
			idsAtPos = posToId.get(chrpos)
			if idsAtPos is None:
				idsAtPos = []
			idsAtPos.append(var)
			posToId[chrpos] = idsAtPos
			matches += 1
	ctr +=1
	if ctr % 10000 == 0:
		print("{} lines parsed, {} matches sofar.".format(ctr, matches),end='\r')
print("{} lines parsed, {} matches total.".format(ctr, matches),end='\n')
fh.close()

fh = gzip.open(infile,'rt')
fho = gzip.open(outfile,'wt')
print("Replacing IDs in "+infile+". Writing to: "+outfile)
ctr += 0
for line in fh:
	if line.startswith("##"):
		fho.write(line)
	elif line.startswith("#"):
		fho.write("##ID matched with: "+reffile+"\n")
		fho.write(line)
	else:
		# try to match variants
		elems = line.split("\t",9)
		chr = elems[0]
		pos = elems[1]
		chrpos = chr+":"+pos
		idsAtPos = posToId.get(chrpos)
		if idsAtPos is None:
			fho.write(line)
		else:
			# try to match alleles per variant
			a1 = elems[3]
			a2 = elems[4]
			elems = line.split("\t")
			match = False
			for var in idsAtPos:
				if match == False:
					refa1 = var.allele1
					refa2 = var.allele2
					nrshared = countSharedAlleles([a1,a2],[refa1,refa2])
					if nrshared == 2:
						elems[3] = var.id
						match = True
		fho.write("\t".join(elems))
	ctr+=1
	if ctr % 10000 == 0:
		print("{} lines parsed, out of {}.".format(ctr, inlines),end='\r')
print("{} lines parsed, out of {}.".format(ctr, inlines),end='\n')
fho.close()
fh.close()
