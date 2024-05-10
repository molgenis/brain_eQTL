import os
import sys
import gzip

# sys.path.append("../../../library/")
from pathlib import Path

path = str(Path(__file__).parent.parent.parent.parent.absolute().__str__() + "/library/")
print("library path: " + path)
sys.path.insert(0, path)

from parsers.GTFAnnotation import GTFAnnotation
from features.splicefeature2 import SpliceFeature

maxdist = 10000


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
        return gzip.open(file, 'rt')
    return open(file)


# load junctions and clusters
junctions = []
clusters = {}
print("Loading splice junctions from: " + splicefile)
fh = getfh(splicefile)
fh.readline()
lctr = 0
for line in fh:
    elems = line.strip().split("\t", 2)
    if len(elems) == 1:
        elems = line.strip().split(" ", 2)
    id = elems[0]
    # print(id)
    junction = SpliceFeature.parseLeafCutter(id)
    if junction.chr.isAutosomal():
        junctions.append(junction)
        clusterid = junction.clusterId
        junctionsInCluster = clusters.get(clusterid)
        if junctionsInCluster is None:
            junctionsInCluster = []
        junctionsInCluster.append(junction)
        clusters[clusterid] = junctionsInCluster
    lctr += 1
    if lctr % 50000 == 0:
        print("{} lines parsed, {} loaded, {} clusters".format(lctr, len(junctions), len(clusters)), end='\r')
        # break
print("{} lines parsed, {} loaded, {} clusters".format(lctr, len(junctions), len(clusters)), end='\n')
fh.close()

annotation = GTFAnnotation(gtffile)
# genesByChr = annotation.getGenesByChromosome()

# annotate genes within each cluster
fho = gzip.open(outfile + "-nearestgene-annotation.txt.gz", 'wt')
fho4 = gzip.open(outfile + "-overlappingGenesAndExons.txt.gz", 'wt')
fho2 = gzip.open(outfile + "-nearestgene-clusters.txt.gz", 'wt')
fho3 = gzip.open(outfile + "-nearestgene-junctions.txt.gz", 'wt')
fho.write("Platform\tFeatureId\tFeatureChr\tFeatureChrStart\tFeatureChrEnd\tMappedGeneId\tMappedGeneSymbol\tMappedGeneCoordinates\tMappedGeneStrand\tMappedGeneType\tMappedGeneDistance\tFeatureOverlapsMappedGene\n")
fho4.write("Platform\tFeatureId\tFeatureChr\tFeatureChrStart\tGene\tGeneSymbol\tGeneCoords\tGeneStrand\tGeneType\tGeneDistance\tTranscriptName\tTranscriptCoordinates\tTranscriptDistance\tExonName\tExonRank\tExonCoords\tExonDistance\n")
jctr = 0

wiggle = 1000000 # look for genes within 1mb of the junction
outputStrFull = ["-"]*18
outputStrNearestGene = [""]*12
nrJunctionsOverlapGenes = 0
nrJunctionsOverlapExons = 0
nrJunctionsOverlapTranscripts = 0
for junction in junctions:
    #genes = genesByChr.get(junction.chr)
    genes = annotation.getOverlappingGenes(junction, wiggle)
    nearestGene = None
    nearestGeneDist = 1e10
    overlapsExon = False
    overlapsGene = False
    overlapsTranscript = False
    if genes is not None:
        for gene in genes:
            dist = gene.absoluteMinimalDistance(junction)
            if gene.overlaps(junction):
                geneCoords = gene.coordToStr()
                
                outputStrFull[0] = "LeafCutter"
                outputStrFull[1] = junction.name
                outputStrFull[2] = str(junction.chr.getNumber())
                outputStrFull[3] = str(junction.start)
                outputStrFull[4] = str(junction.stop)
                outputStrFull[5] = gene.name
                outputStrFull[6] = gene.symbol
                outputStrFull[7] = geneCoords
                outputStrFull[8] = gene.strand.toStr()
                outputStrFull[9] = gene.type
                outputStrFull[10] = str(dist)
                transcripts = gene.transcripts
                overlapsGene = True
                fho4.write("\t".join(outputStrFull)+"\n")
                for transcript in transcripts:
                    toverlap = transcript.overlaps(junction)
                    if toverlap:
                        tdist = transcript.absoluteMinimalDistance(junction)
                        tcoords = transcript.coordToStr()

                        outputStrFull[11] = transcript.name
                        outputStrFull[12] = tcoords
                        outputStrFull[13] = str(tdist)
                        outputStrFull[14] = "-"
                        outputStrFull[15] = "-"
                        outputStrFull[16] = "-"
                        outputStrFull[17] = "-"
                        fho4.write("\t".join(outputStrFull)+"\n")
                        exons = transcript.exons
                        overlapsTranscript = True
                        for exon in exons:
                            eoverlap = exon.overlaps(junction)
                            if eoverlap:
                                overlapsExon = True
                                ecoords = exon.coordToStr()
                                edist = exon.absoluteMinimalDistance(junction)
                                rank = transcript.getExonRank(exon)
                                if rank is None:
                                    rank = "-"
                                outputStrFull[14] = exon.name
                                outputStrFull[15] = str(rank)
                                outputStrFull[16] = ecoords
                                outputStrFull[17] = str(edist)
                                fho4.write("\t".join(outputStrFull)+"\n")

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
    if nearestGene is not None:
        nearestGeneCoords = nearestGene.coordToStr()
        
        
        outputStrNearestGene[0] = "LeafCutterNearestGene"
        outputStrNearestGene[1] = junction.name
        outputStrNearestGene[2] = str(junction.chr.getNumber())
        outputStrNearestGene[3] = str(junction.start)
        outputStrNearestGene[4] = str(junction.stop)
        outputStrNearestGene[5] = nearestGene.name
        outputStrNearestGene[6] = gene.symbol
        outputStrNearestGene[7] = nearestGeneCoords
        outputStrNearestGene[8] = nearestGene.strand.toStr()
        outputStrNearestGene[9] = nearestGene.type
        outputStrNearestGene[10] = str(nearestGeneDist)
        outputStrNearestGene[11] = str(nearestGene.overlaps(junction))
        fho.write("\t".join(outputStrNearestGene)+"\n")
   
    else:
        print(bcolors.FAIL + "No gene for " + junction.name + " - " + str(junction.chr.getNumber()) + bcolors.ENDC)
        outputStrNearestGene[0] = "LeafCutterNearestGene"
        outputStrNearestGene[1] = "-"
        outputStrNearestGene[2] = str(junction.chr.getNumber())
        outputStrNearestGene[3] = str(junction.start)
        outputStrNearestGene[4] = str(junction.stop)
        outputStrNearestGene[5] = "-"
        outputStrNearestGene[6] = "-"
        outputStrNearestGene[7] = "-"
        outputStrNearestGene[8] = "-"
        outputStrNearestGene[9] = "-"
        outputStrNearestGene[10] = "False"
        fho.write("\t".join(outputStrNearestGene)+"\n")
   
    fho2.write(junction.name + "\t" + junction.clusterId + "\n")
    fho3.write(junction.name + "\n")
    jctr += 1

    if overlapsGene:
        nrJunctionsOverlapGenes+=1
    if overlapsTranscript:
        nrJunctionsOverlapTranscripts+=1
    if overlapsExon:
        nrJunctionsOverlapExons += 1
    if jctr % 1000 == 0:
        print(f"{jctr} / {len(junctions)} junctions written, overlap: {nrJunctionsOverlapGenes} gene {nrJunctionsOverlapTranscripts} transcript {nrJunctionsOverlapExons} exon", end='\r')
        # break
print(f"{jctr} / {len(junctions)} junctions written, overlap: {nrJunctionsOverlapGenes} gene {nrJunctionsOverlapTranscripts} transcript {nrJunctionsOverlapExons} exon", end='\n')
fho4.close()
fho3.close()
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
        print(bcolors.FAIL + "cluster {} overlaps {} genes and has {} members".format(cluster, len(overlappingGenes),
                                                                                      len(junctions)))
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
                                overlappingExons = True

                    if overlappingExons:
                        genestr += "\t{}({}-{}); {} - ov: {}".format(gene.symbol, gene.start, gene.stop, gene.strand,
                                                                     bpoverlap) + "-" + bcolors.OKCYAN + "TRUE" + bcolors.ENDC
                    else:
                        genestr += "\t{}({}-{}); {} - ov: {}".format(gene.symbol, gene.start, gene.stop, gene.strand,
                                                                     bpoverlap) + "-" + bcolors.FAIL + "FALSE" + bcolors.ENDC
                    if nearestGene is None:
                        nearestGene = gene
                        nearestGeneDist = gene.absoluteMinimalDistance(junction)
                    else:
                        # compare the distance to this gene with the distance to the other gene
                        mindist = gene.absoluteMinimalDistance(junction)
                        if mindist < nearestGeneDist:
                            nearestGeneDist = mindist
                            nearestGene = gene

            print(
                bcolors.OKBLUE + junction.name + bcolors.ENDC + genestr + "\tNearest:" + bcolors.OKGREEN + nearestGene.symbol + bcolors.ENDC)
        print()
        cctr += 1
        if cctr % 100 == 0:
            sys.exit()
        pass
    elif len(overlappingGenes) == 1:
        print(bcolors.OKGREEN + "cluster {} overlaps 1 gene and has {} members".format(cluster,
                                                                                       len(junctions)) + bcolors.ENDC)
        fho.write("Platform\tEnsembl\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tStrand\n")
    else:
        print(bcolors.OKGREEN + "cluster {} overlaps no genes and has {} members".format(cluster,
                                                                                         len(junctions)) + bcolors.ENDC)
        fho.write("Platform\tEnsembl\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tStrand\n")
