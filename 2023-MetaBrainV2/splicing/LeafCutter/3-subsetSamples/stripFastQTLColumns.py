import gzip
import sys


if len(sys.argv) < 3:
        print("Usage: input.txt.gz output.txt.gz")
        sys.exit()

input = sys.argv[1]
output = sys.argv[2]

fh = gzip.open(input,'rt')
fho = gzip.open(output, 'wt')

x = 0
for line in fh:
        elems = line.split("\t")
        outln = "\t".join(elems[3:len(elems)])
        x+= 1
        if x % 1000 == 0:
                print(f"{x}",end='\r')
        fho.write(outln)
fho.close()