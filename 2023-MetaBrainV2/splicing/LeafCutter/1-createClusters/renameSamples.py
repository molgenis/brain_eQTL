import gzip
import sys


def getrfh(fname):
    if fname[-3:] == ".gz":
        return gzip.open(fname, 'rt')
    else:
        return open(fname,'r'   )

def getwfh(fname):
    if fname[-3:] == ".gz":
        return gzip.open(fname, 'wt', 5)
    else:
        return open(fname,'w')        

if len(sys.argv) < 3:
    print("Usage: infile.txt.gz cramToRNAId.txt outfile.txt.gz")
    sys.exit()

infile = sys.argv[1]
ctorfile = sys.argv[2]
outfile = sys.argv[3]

ctor = {}
fh = getrfh(ctorfile)
for line in fh:
    elems = line.strip().split("\t")
    r = elems[1]
    c = elems[0]
    ctor[c] = r
fh.close()

print(f"{len(ctor)} cram to RNA links loaded")

fh = getrfh(infile)
sep = "\t"
headerln = fh.readline().strip()
header = headerln.split(sep)
if len(header) == 1:
    sep = " "
header = headerln.split(sep)
if len(header) == 1:
    print("could not determine line separator")
    sys.exit()

headerout = ""
found = 0
mismatch = 0
for sample in header:
    other = ctor.get(sample)
    if other is None:
        other = sample
        print("Could not match "+sample)
        mismatch += 1
    else:
        found+= 1
    if len(headerout) == 0:
        headerout += other  
    else:
        headerout += sep+other

if found == 0:
    print("Could not match any cram id to rna id")
    fh.close()
    sys.exit()
else:
    print(f"{found} matches between cram id and rna id")

# if mismatch > 0:
#     print("Could not match all cram ids. Shutting down")
#     sys.exit()

fho = getwfh(outfile)
fho.write(headerout+"\n")
lctr = 0
for line in fh:
    fho.write(line)
    lctr += 1
    if lctr % 1000 == 0:
        print(f"{lctr} lines parsed",end='\r')
print(f"{lctr} lines parsed. Done.",end='\n')      
fho.close()
fh.close()