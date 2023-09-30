import sys
import gzip
import glob

if len(sys.argv) < 3:
	print("Usage: /dir/with/topeffects/ outfile.txt")
	sys.exit()
dir = sys.argv[1]
out = sys.argv[2]

print(dir)
print(out)

# qvalue = importr("qvalue")

files = glob.glob(dir+"/*-TopEffects.txt")
print("{} top effect files found".format(len(files)) +" in "+dir)
if len(files) > 0:
	header = ""
	fctr = 0
	lines = []
	pvals = []
	for file in files:
		fh = open(file,'r')
		if fctr == 0:
			header = fh.readline().strip()
		else:
			fh.readline()
		for line in fh:
			elems = line.strip().split("\t")
			if len(elems) <= 11:
				print(line)
				continue
			pval = float(elems[11])
			lines.append( [line.strip(), pval] )
		fh.close()
		fctr = fctr + 1
	
	# calculate Qvalues
#	pvalsr = robjects.FloatVector(pvals)
#	qobj = robjects.r['qvalue'](pvalsr)
# 	qvals = qobj.rx2('qvalues')
	# append to lists

	# sort by qvalue
	lines = sorted(lines, key=lambda p: p[1]) 
	fho = open(out, 'w+')
	fho.write(header+"\n")
	for line in lines:
#		fho.write(line[0]+"\t"+str(line[2])+"\n")
		fho.write(line[0]+"\n")
	fho.close()
