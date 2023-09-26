from enum import Enum
from features.feature import Feature

class Gene(Feature):
    name = None
    symbol = None
    transcripts = None
    type = None

    def __init__(self, chr, start, stop, strand, name, symbol):
        Feature.__init__(self,chr,start,stop,strand)
        self.name = name
        self.symbol = symbol
    
    def addTranscript(self, transcript):
        if self.transcripts is None:
            self.transcripts = []
        self.transcripts.append(transcript)

    def setType(self, type):
        self.type = type

    def describe(self):
        nrTranscripts = 0
        if self.transcripts is not None:
            nrTranscripts = len(self.transcripts)
        return "Gene: {} ({}) - {}:{}:{} - Strand: {} - Transcripts: {} - Type: {}".format(self.name, self.symbol, self.chr.name, self.start, self.stop, self.strand.name, nrTranscripts, self.type)
