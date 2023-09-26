from features.feature import Feature


class SpliceFeature(Feature):
    name = None
    genes = None
    transcripts = None
    exons = None

    def __init__(self, chr, start, stop, strand, name):
        Feature.__init__(self,chr,start,stop,strand)
        self.name = name
    
    def addGene(self, gene):
        if self.genes is None:
            self.genes = set()
        self.genes.add(gene)
    
    def addTranscript(self, transcript):
        if self.transcripts is None:
            self.transcripts = set()
        self.transcripts.add(transcript)

    def addTranscript(self, exon):
        if self.exons is None:
            self.exons = set()
        self.transcripts.add(exon)

    def describe(self):
        return "SpliceSite: {}:{}-{} {}".format(self.chr,self.start,self.stop,self.name)