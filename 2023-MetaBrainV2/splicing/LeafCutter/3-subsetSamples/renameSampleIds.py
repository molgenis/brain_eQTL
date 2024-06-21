import gzip
import sys

if len(sys.argv) < 4:
        print("Usage: input.txt.gz output.txt.gz linkfile.txt.gz")
        sys.exit()

input = sys.argv[1]
output = sys.argv[2]
linkfile = sys.argv[3]

fh = gzip.open(linkfile,'rt')
map = {}
for line in fh:
        elems = line.strip().split("\t")
        map[elems[0]] = elems[1]
fh.close()

fh = gzip.open(input,'rt')
fho = gzip.open(output,'wt',4)
header = fh.readline().strip().split(" ")
match = 0
for i in range(1,len(header)):
        c = header[i]
        r = map.get(c)
        if r is not None:
                header[i] = r
                match += 1
print(f"{match} matched")
fho.write(" ".join(header)+"\n")
lctr = 0
for line in fh:
        fho.write(line)
        lctr += 1
        if lctr % 1000 == 0:
                print(f"{lctr} lines",end='\r')
fho.close()