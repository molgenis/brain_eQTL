import gzip
import sys

if len(sys.argv) < 4:
	print("Usage: countfile.gz junctionincludelist.txt outfile.gz")
	sys.exit()

protcod = set()
fh = open(sys.argv[2],'r')
q = 0
for line in fh:
	protcod.add(line.strip())
	if q < 3:
		print(line.strip())
	q+=1
fh.close()

print(f"{len(protcod)} protein coding clusters")

countfile= sys.argv[1]
outfile  = sys.argv[3]

print("Parsing: "+countfile)
print("Writing: "+outfile)
fh  = gzip.open(countfile,'rt')
fho = gzip.open(outfile,'wt')

fho.write(fh.readline())
lctr = 0
written = 0
for line in fh:
	elems = line.split(" ",2)
	id = elems[0]
	if id in protcod:
		fho.write(line)
		written += 1
	lctr += 1
	if lctr % 1000 == 0:
		print(f"{lctr} lines parsed, {written} written",end='\r')
print(f"{lctr} lines parsed, {written} written",end='\n')
fho.close()
fh.close()
