import sys
import os

tissues = tissues = ["Basalganglia-EUR", "Cerebellum-EUR", "Cortex-EUR", "Cortex-AFR", "Cortex-EAS", "Hippocampus-EUR", "Spinalcord-EUR"]

numpcs = 30
stepsize = 5

print("Tissue\tPCs\tGenes\tSignificantGenes\tEvents\tSignificantEvents")

for tissue in tissues:
	for pcs in range(0, numpcs+stepsize, stepsize):
		file = "./output/PCACurve/"+tissue+"-"+str(pcs)+"PCs/output/merged-withqval.txt"

		genes = set()
		events = set()
		siggenes = set()
		sigevents = set()
		if os.path.exists(file):
			fh = open(file,'r')
			fh.readline() # skip header
			for line in fh:
				elems = line.strip().split("\t")
				symbol = elems[5]
				event = elems[1]
				qval = float(elems[-1])
				genes.add(symbol)
				events.add(event)
				if qval < 0.05:
					siggenes.add(symbol)
					sigevents.add(event)
			fh.close()
		print("{}\t{}\t{}\t{}\t{}\t{}".format( tissue, pcs, len(genes), len(siggenes), len(events), len(sigevents) ) )
