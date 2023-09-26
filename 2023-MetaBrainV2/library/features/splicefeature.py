from features.feature import Feature
from features.chromosome import Chromosome
from features.strand import Strand

class SpliceFeature(Feature):
    name = None
    genes = None
    transcripts = None
    exons = None
    clusterId = None

    def __init__(self, chr, start, stop, strand, name):
        Feature.__init__(self,chr,start,stop,strand)
        self.name = name

    @classmethod
    def parseLeafCutter(self,str):
        idelems = str.split(":")
        if len(idelems) < 4:
            print("Error: "+str+" is not a leafcutter ID") # should throw an exception
            return None
        chr = Chromosome.parse(idelems[0])
        sta = int(idelems[1])
        sto = int(idelems[2])
        obj = SpliceFeature(chr,sta,sto,Strand.NA,str)
        obj.setClusterId(idelems[3])
        return obj

    def setClusterId(self, clu):
        self.clusterId = clu
        
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