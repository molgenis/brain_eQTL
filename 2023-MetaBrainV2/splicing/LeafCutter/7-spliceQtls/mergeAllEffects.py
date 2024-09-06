import gzip
import glob
import sys

if len(sys.argv) < 3:
	print("Usage: indir chr")
	sys.exit()

dir = sys.argv[1]
chr = sys.argv[2]

files = glob.glob(f"{dir}/chr{chr}-batch-*-AllEffects.txt.gz")
print(f"{len(files)} batch output files in {dir}")
outf = f"{dir}/chr{chr}-AllEffects.txt.gz"

print(f"Writing to: {outf}")

fho = gzip.open(outf,'wt',3)

files.sort()

fctr = 1
for file in files:
	print(f"Merging: {file} - {fctr}/{len(files)}")
	fh = gzip.open(file,'rt')
	if fctr == 1:
		fho.write(fh.readline())
	else:
		fh.readline()
	for line in fh:
		fho.write(line)
	fh.close()
	fctr += 1
fho.close()
print("Done")
print()
