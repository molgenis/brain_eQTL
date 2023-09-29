import os
import sys
import gzip

# sys.path.append("../../../library/")
from pathlib import Path


path = str(Path(__file__).parent.parent.parent.parent.absolute().__str__()+"/library/")
print("library path: "+path)
sys.path.insert(0,path)


from parsers.GTFAnnotation import GTFAnnotation
from features.splicefeature import SpliceFeature


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'



if len(sys.argv) < 4:
    print("Usage: gtf[.gz] splicefile outfile")
    sys.exit()
gtffile = sys.argv[1]
splicefile = sys.argv[2]
outfile = sys.argv[3]


# gtffile = "/groups/umcg-biogen/tmp02/annotation/Gencode/v32-b38/gencode.v32.annotation.gtf.gz"
# splicefile = ""
# outfile = ""

def getfh(file):
    if file.endswith(".gz"):
        return gzip.open(file,'rt')
    return open(file)

# load junctions and clusters
junctions = []
clusters = {}
print("Loading splice junctions from: "+splicefile )
fh = getfh(splicefile)
fh.readline()
lctr = 0
for line in fh:
    elems = line.strip().split("\t",2)
    if len(elems) == 1:
        elems = line.strip().split(" ",2)
    id = elems[0]
    # print(id)
    junction = SpliceFeature.parseLeafCutter(id)
    junctions.append(junction)
    clusterid = junction.clusterId
    junctionsInCluster = clusters.get(clusterid)
    if junctionsInCluster is None:
        junctionsInCluster = []
    junctionsInCluster.append(junction)
    clusters[clusterid] = junctionsInCluster
    lctr += 1
    if lctr % 50000 == 0:
        print("{} lines parsed, {} loaded, {} clusters".format(lctr,len(junctions),len(clusters)),end='\r')
        # break
print("{} lines parsed, {} loaded, {} clusters".format(lctr,len(junctions),len(clusters)),end='\n')
fh.close()


annotation = GTFAnnotation(gtffile)
genesByChr = annotation.getGenesByChromosome()

# annotate genes within each cluster
fho = gzip.open(outfile+"-nearestgene-annotation.txt.gz",'wt')
fho2 = gzip.open(outfile+"-nearestgene-clusters.txt.gz",'wt')
fho.write("Platform\tEnsembl\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tStrand\tTypeOfGene\n")
jctr = 0
for junction in junctions:
    genes = genesByChr.get(junction.chr)
    nearestGene = None
    nearestGeneDist = 1e10
    for gene in genes:
        dist = gene.absoluteMinimalDistance(junction)
        if nearestGene is None:
            nearestGene = gene
            nearestGeneDist = dist
            # print("new nearest gene: {}\t{}\t{}".format(junction.name,gene.name,dist))
        else:
            if dist < nearestGeneDist:
                nearestGeneDist = dist
                nearestGene = gene
                # print("new nearest gene: {}\t{}\t{}".format(junction.name,gene.name,dist))
    gene = nearestGene
    fho.write("leafcutter\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(junction.name, gene.symbol, junction.chr.getNumber(), junction.start, junction.stop, gene.strand.toStr(), gene.name, gene.type))
    fho2.write(junction.name+"\t"+junction.clusterId+"\n")
    jctr += 1
    if jctr % 1000 == 0:
        print("{} junctions written".format(jctr),end='\r')
        # break
print("{} junctions written".format(jctr),end='\n')

fho2.close()
fho.close()
sys.exit()

# more finegrained annotation
for cluster in clusters.keys():
    junctions = clusters.get(cluster)
    overlappingGenes = set()
    overlappingGeneCtr = {}

    genesPerJunction = {}
    for junction in junctions:
        chr = junction.chr
        genes = genesByChr.get(chr)
        if genes is None:
            break
        for gene in genes:
            if gene.overlaps(junction):
                overlappingGenes.add(gene)
                gctr = overlappingGeneCtr.get(gene)
                if gctr is None:
                    gctr = 0
                gctr += 1
                overlappingGeneCtr[gene] = gctr
                genesForJunction = genesPerJunction.get(junction)
                if genesForJunction is None:
                    genesForJunction = []
                genesForJunction.append(gene)
                genesPerJunction[junction] = genesForJunction

    if len(overlappingGenes) > 1:
        # do some magic to assign the most probable one
        print(bcolors.FAIL+"cluster {} overlaps {} genes and has {} members".format(cluster, len(overlappingGenes), len(junctions)))
        for junction in junctions:
            genesForJunction = genesPerJunction.get(junction)
            genestr = "\tNA"
            genesWithOverlappingExons = set()
            if genesForJunction is not None:
                genestr = ""
                nearestGene = None
                nearestGeneDist = None
                for gene in genesForJunction:
                    overlappingExons = False
                    bpoverlap = gene.bpOverlap(junction)


                    for transcript in gene.transcripts:
                        for exon in transcript.exons:
                            if exon.overlaps(junction):
                                genesWithOverlappingExons.add(gene)
                                overlappingExons=True
                                
                    if overlappingExons:
                        genestr += "\t{}({}-{}); {} - ov: {}".format(gene.symbol, gene.start, gene.stop,gene.strand,bpoverlap) + "-"+bcolors.OKCYAN+"TRUE"+bcolors.ENDC
                    else:
                        genestr += "\t{}({}-{}); {} - ov: {}".format(gene.symbol, gene.start, gene.stop,gene.strand,bpoverlap)+ "-"+bcolors.FAIL+"FALSE"+bcolors.ENDC
                    if nearestGene is None:
                        nearestGene = gene
                        nearestGeneDist = gene.absoluteMinimalDistance(junction)
                    else:
                        # compare the distance to this gene with the distance to the other gene
                        mindist = gene.absoluteMinimalDistance(junction)
                        if mindist < nearestGeneDist:
                            nearestGeneDist = mindist
                            nearestGene = gene

                        
            print(bcolors.OKBLUE+junction.name+ bcolors.ENDC+genestr+"\tNearest:"+bcolors.OKGREEN+nearestGene.symbol+bcolors.ENDC)
        print()
        cctr += 1
        if cctr % 100 == 0:
            sys.exit()
        pass
    elif len(overlappingGenes) == 1:
        print(bcolors.OKGREEN+"cluster {} overlaps 1 gene and has {} members".format(cluster, len(junctions))+ bcolors.ENDC) 
        fho.write("Platform\tEnsembl\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tStrand\n")
    else:
        print(bcolors.OKGREEN+"cluster {} overlaps no genes and has {} members".format(cluster, len(junctions))+ bcolors.ENDC)
        fho.write("Platform\tEnsembl\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tStrand\n")
