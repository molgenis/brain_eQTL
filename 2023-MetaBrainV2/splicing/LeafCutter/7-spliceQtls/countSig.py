import sys

if len(sys.argv) < 3:
	print("Usage: infile.txt outfile.txt")
	sys.exit(0)
file = sys.argv[1]
ofile = sys.argv[2]

fh = open(file,'r')
fho = open(ofile,'w')
header = fh.readline()
fho.write(header)
header = header.strip().split("\t")
ctr = 0
qvcol = -1
for h in header:
	if h.lower() == "qval":
		qvcol = ctr
	ctr = ctr + 1
if qvcol < 0:
	print("no qval col found")
else:
	sig = 0
	tot = 0
	for line in fh:
		elems = line.strip().split("\t")
		qv = float(elems[qvcol])
		if qv < 0.05:
			sig = sig + 1
			fho.write(line)
		tot = tot + 1

#	print("Significant:\t{}\tTotal:\t{}".format(sig, tot))
	print("{}".format(sig))
		
fho.close()
