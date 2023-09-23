import gzip
import sys
from collections import defaultdict

if len(sys.argv) < 4:
    print("usage: splicefile.txt.gz outfile.txt.gz gtf-file")
    sys.exit(0)

# spliceFilePath = "/Users/harmbrugge/Documents/apps/meta-brain-sqtls/files/merged_fractions_logit_filtered_mds_covariates_residuals_chr1.txt"
# outFilePath = "/Users/harmbrugge/Documents/apps/meta-brain-sqtls/files/annot_test.txt"
# geneAnnotPath = "/Users/harmbrugge/Documents/apps/meta-brain-sqtls/files/annotations/gencode.v32.primary_assembly.annotation_chr1.gtf"
spliceFilePath = sys.argv[1]
outFilePath = sys.argv[2]
geneAnnotPath = sys.argv[3]

print("Creating annotation files")
print("Input: " + spliceFilePath)
print("Output: " + outFilePath)

geneAnnotation = gzip.open(geneAnnotPath, 'rt')
for _ in range(5):
    geneAnnotation.readline()

genes = defaultdict(list)
currentGene = dict()
currentTranscript = dict()
currentExon = dict()

print("Reading annotation file")
for line in geneAnnotation:
    elems = line.strip().split("\t")
    chromosome = elems[0]
    annotationType = elems[2]
    start = int(elems[3])
    stop = int(elems[4])
    strand = elems[6]
    if annotationType == "gene":
        currentGene = dict()
        geneInfo = elems[8].split(";")
        currentGene["ensg"] = geneInfo[0].split("\"")[1]
        currentGene["symbol"] = geneInfo[2].split("\"")[1]
        currentGene["start"] = start
        currentGene["stop"] = stop
        currentGene["strand"] = strand
        currentGene["transcripts"] = []
        genes[chromosome].append(currentGene)
    if annotationType == "transcript":
        currentTranscript = dict()
        currentTranscript["start"] = start
        currentTranscript["stop"] = stop
        transcriptInfo = elems[8].split(";")
        transcriptId = transcriptInfo[1].split("\"")[1]
        currentTranscript["id"] = transcriptId
        currentTranscript["exons"] = []
        currentGene["transcripts"].append(currentTranscript)
    if annotationType == "exon":
        currentExon = dict()
        currentExon["start"] = start
        currentExon["stop"] = stop
        exonInfo = elems[8].split(";")
        exonId = exonInfo[7].split("\"")[1]
        currentExon["id"] = exonId
        currentTranscript["exons"].append(currentExon)
geneAnnotation.close()
print("Done reading annotation file")

spliceFile = gzip.open(spliceFilePath, 'rt')
spliceFile.readline()

outFile = gzip.open(outFilePath, 'wt')
outFile.write("Platform\tEnsembl\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tStrand\n")

outFileClusters = gzip.open(outFilePath[:-7] + "_clusters.txt.gz", 'wt')

lineCount = 0
noGeneCount = 0
noExonMoreGenesCount = 0
perfectCount = 0


def determineGene(junction_genes, junction_begin_exons, junction_end_exons):
    for gene in junction_genes:
        for transcript in gene["transcripts"]:
            for exon_id in junction_begin_exons:
                if exon_id in transcript["exons"]:
                    return gene
            for exon_id in junction_end_exons:
                if exon_id in transcript["exons"]:
                    return gene
    return None


for line in spliceFile:
    elems = line.strip().split()
    spliceEventId = elems[0]

    elems = spliceEventId.split(":")
    chromosome = elems[0]
    if chromosome == "chrY" or chromosome == "chrX":
        continue
    chromosome_nr = chromosome[3:]
    junction_begin = int(elems[1])
    junction_end = int(elems[2])
    cluster_id = elems[3][0:-2]

    junction_begin_exons = set()
    junction_end_exons = set()
    junction_genes = list()

    for gene in genes[chromosome]:
        if junction_begin >= gene["start"] and junction_end <= gene["stop"]:
            junction_genes.append(gene)
            for transcript in gene["transcripts"]:
                if junction_begin >= transcript["start"] and junction_end <= transcript["stop"]:
                    for exon in transcript["exons"]:
                        if junction_begin == exon["start"] or junction_begin == exon["stop"]:
                            junction_begin_exons.add(exon["id"])
                        if junction_end == exon["start"] or junction_end == exon["stop"]:
                            junction_end_exons.add(exon["id"])

    if len(junction_genes) == 1:
        outFile.write("splicing\t" + spliceEventId + "\t" + junction_genes[0]["symbol"] + "\t" + chromosome_nr +
                      "\t" + str(junction_begin) + "\t" + str(junction_end) + "\t" + junction_genes[0]["ensg"] + "\t" +
                      junction_genes[0]["strand"] + "\n")
        outFileClusters.write(spliceEventId + "\t" + cluster_id + "\n")
        if len(junction_begin_exons) == 1 and len(junction_end_exons) == 1:
            perfectCount += 1
        elif not junction_begin_exons and not junction_end_exons:
            a = 0
    elif len(junction_genes) >= 1:
        gene = determineGene(junction_genes, junction_begin_exons, junction_end_exons)
        if gene:
            outFile.write("splicing\t" + spliceEventId + "\t" + gene["symbol"] + "\t" + chromosome_nr +
                          "\t" + str(junction_begin) + "\t" + str(junction_end) + "\t" + gene["ensg"] + "\t" +
                          gene["strand"] + "\n")
            outFileClusters.write(spliceEventId + "\t" + cluster_id + "\n")
    else:
        noGeneCount += 1

    lineCount += 1
    if lineCount % 1000 == 0:
        print(lineCount, "lines processed")
        print(noGeneCount, "junction with no overlapping genes")
        print(noExonMoreGenesCount, "multiple genes no exons")
        print(perfectCount, "perfect count")

outFile.close()
outFileClusters.close()

