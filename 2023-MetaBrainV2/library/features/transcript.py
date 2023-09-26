from enum import Enum
from features.feature import Feature

class Transcript(Feature):
    gene = None
    name = None
    exons = None
    exonRanks = None

    def __init__(self, chr, start, stop, strand, name, gene):
        Feature.__init__(self,chr,start,stop,strand)
        self.name = name
        self.gene = gene

    def addExon(self, exon):
        if self.exons is None:
            self.exons = []
        self.exons.append(exon)

    def setExonRank(self, exon, rank):
        if self.exonRanks is None:
            self.exonRanks = {}
        self.exonRanks[exon] = rank

    def describe(self):
        nrExons = 0
        if self.exons is not None:
            nrExons = len(self.exons)
        return "Transcript: {} - {} ({}) - {}:{}:{} - Strand: {} - Exons: {}".format(self.name, self.gene.name, self.gene.symbol, self.chr.name, self.start, self.stop, self.strand.name,nrExons)
