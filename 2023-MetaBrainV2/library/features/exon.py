from enum import Enum
from features.feature import Feature

class Exon(Feature):
    gene = None
    transcripts = None
    name = None

    def __init__(self, chr, start, stop, strand, name, gene):
        Feature.__init__(self,chr,start,stop,strand)
        self.name = name
        self.gene = gene

    def addTranscript(self, transcript):
        if self.transcripts is None:
            self.transcripts = set()
        self.transcripts.add(transcript)

    def describe(self):
        nrTranscripts = 0
        if self.transcripts is not None:
            nrTranscripts = len(self.transcripts)
        return "Exon: {} - {} ({}) - {}:{}:{} - Strand: {} - Transcripts: {}".format(self.name, self.gene.name, self.gene.symbol, self.chr.name, self.start, self.stop, self.strand.name, nrTranscripts)
