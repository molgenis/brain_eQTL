import gzip

from features.feature import Feature
from features.chromosome import Chromosome
from features.exon import Exon
from features.transcript import Transcript
from features.strand import Strand
from features.gene import Gene



class GTFAnnotation:

    genes = []
    transcripts = []
    exons = []

    idToTranscript = {}
    idToExon = {}
    idToGene = {}
    notparsedtypes = set()

    genesPerChr = {}

    def __init__(self, gtffile):      
        self.parse(gtffile)

    def getfh(self,name):
        if name.endswith(".gz"):
            return gzip.open(name, 'rt')
        else:
            return open(name,'r')

    def parse(self, gtffile):
        print("Reading GTF: "+gtffile)
        fh = self.getfh(gtffile)
        lctr = 0
        for line in fh:
            line = line.strip()
            if not line.startswith("#"):
                self.parseln(line)
            lctr+=1
            if lctr % 100000 == 0:
                print("{} lines parsed, {} genes, {} transcripts, {} exons".format(lctr,len(self.genes),len(self.transcripts), len(self.exons)), end="\r")
        print("{} lines parsed, {} genes, {} transcripts, {} exons - Done".format(lctr,len(self.genes),len(self.transcripts), len(self.exons)), end="\n")
        fh.close()
        for type in self.notparsedtypes:
            print("Unknown type of feature in file: "+type)
        if self.genes is not None:
            self.genesPerChr = {}
            for gene in self.genes:
                arr = self.genesPerChr.get(gene.chr)
                if arr is None:
                    arr = []
                arr.append(gene)
                self.genesPerChr[gene.chr] = arr

    def getOrCreateGene(self, chr,start,stop,strand,annotation):
        gid = annotation.get("gene_id")
        if gid is None:
            return None
        gobj = self.idToGene.get(gid)
        if gobj is None:
            symbol = annotation.get("gene_name")
            
            gobj = Gene(chr,start,stop,strand,gid,symbol)
            gtype = annotation.get("gene_type")
            if gtype is not None:
                gobj.setType(gtype)
            self.genes.append(gobj)
            self.idToGene[gid] = gobj
            # print(gobj.describe())
        return gobj

    def getOrCreateTranscript(self, chr,start,stop,strand,annotation):
        gobj = self.getOrCreateGene(chr,start,stop,strand,annotation)
        if gobj is None:
            return None
        tid = annotation.get("transcript_id")
        if tid is None:
            return None
        tobj = self.idToTranscript.get(tid)
        if tobj is None:
            tobj = Transcript(chr,start,stop,strand,tid,gobj)
            self.idToTranscript[tid] = tobj
            gobj.addTranscript(tobj)
            self.transcripts.append(tobj)
            # print(tobj.describe())
        return tobj

    def getOrCreateExon(self, chr,start,stop,strand,annotation):
        eid = annotation.get("exon_id")
        if eid is None:
            return None
        eobj = self.idToExon.get(eid)
        if eobj is None:
            gobj = self.getOrCreateGene(chr,start,stop,strand,annotation)
            if gobj is None:
                return None
            # print("e1: "+gobj.describe())
            tobj = self.getOrCreateTranscript(chr,start,stop,strand,annotation)
            if tobj is None:
                return None
            # print("e2: "+tobj.describe())
            eobj = Exon(chr,start,stop,strand,eid,gobj)
            if tobj is not None:
                tobj.addExon(eobj)
                exonRank = annotation.get("exon_number")
                if exonRank is not None:
                    try:
                        exonRankI = int(exonRank)
                        tobj.setExonRank(eobj,exonRankI)
                    except:
                        pass
                eobj.addTranscript(tobj)
                self.exons.append(eobj)
            self.idToExon[eid] = eobj
        return eobj
        # print(eobj.describe())
        # print(tobj.describe())
        # print(gobj.describe())

    def getGenesByChromosome(self):
        output = {}
        for gene in self.genes:
            chr = gene.chr
            genesInChr = output.get(chr)
            if genesInChr is None:
                genesInChr = []
            genesInChr.append(gene)
            output[chr] = genesInChr
        return output


    def parseln(self, line):
        elems = line.split("\t")
        chr = Chromosome.parse(elems[0])
        type = elems[2].lower()
        start = int(elems[3])
        stop = int(elems[4])
        strand = Strand.parse(elems[6])
        
        annotation = {}
        if len(elems) > 8:
            annotation = self.toDict(elems[8],"; "," ")
        else:
            print("Weird line:")
            print(line)
            sys.exit()

        if type == "gene":
            # print("Gene")
            gid = annotation.get("gene_id")
            if gid is not None:
                gobj = self.idToGene.get(gid)
                if gobj is None:
                    gobj = self.getOrCreateGene(chr,start,stop,strand,annotation)
                    # print(gobj.describe())
                    # sys.exit()
                else:
                    # gene already exists; update coordinates?
                    print(gid+" already exists")
                    pass
        elif type == "transcript":
            # print("Transcript")
            tid = annotation.get("transcript_id")
            if tid is not None:
                tobj = self.idToTranscript.get(tid)
                if tobj is None:
                    tobj = self.getOrCreateTranscript(chr,start,stop,strand,annotation)
                else:
                    # transcript already exists; print warning?
                    print(tid+" already exists")
                    pass
        elif type == "cds":
            # protein coding region; not implemented yet
            pass
        elif type == "start_codon":
            # protein coding region; not implemented yet
            pass
        elif type == "stop_codon":
            # protein coding region; not implemented yet
            pass
        elif type == "utr":
            # non-translated region; not implemented yet
            pass
        elif type == "selenocyteine":
            # protein related region; not implemented yet
            pass                
        elif type == "exon":
            # print("Exon")
            eid = annotation.get("exon_id")

            # Ensembl exon IDs on para-autosomal reqions do not have the _PAR_Y extension. Fix that here
            gid = annotation.get("gene_id")
            if gid is not None:
                if gid.endswith("_PAR_Y"):
                    paraautosomal = True
                    eid = eid+"_PAR_Y"
                    annotation["exon_id"] = eid

            if eid is not None:
                eobj = self.idToExon.get(eid)
                if eobj is None:
                    eobj = self.getOrCreateExon(chr,start,stop,strand,annotation)
                    # sys.exit()
                else:
                    # check if the positions are the same
                    sameCoordinates = self.validateCoordinates(eid,chr,start,stop,strand,eobj, False)
                    if sameCoordinates:
                        # exon already exist, is also used in another transcript?
                        tobj = self.getOrCreateTranscript(chr,start,stop,strand,annotation)
                        tobj.addExon(eobj)
                        eobj.addTranscript(tobj)
                    else:
                        # 
                        # create new exon object. this seems to happen if the GTF was remapped/lifted over
                        itr = 1
                        fixeid = eid+"_"+str(itr)
                        found = False
                        while self.idToExon.get(fixeid) is not None:
                            fixeobj = self.idToExon.get(fixeid)
                            # apparently, a fixed exon exists; check coordinates
                            sameCoordinates = self.validateCoordinates(fixeid,chr,start,stop,strand,fixeobj, False)
                            if sameCoordinates:
                                found = True
                                eobj = fixeobj
                                break
                            itr += 1
                            fixeid = eid+"_"+str(itr)

                        if found:
                            tobj = self.getOrCreateTranscript(chr,start,stop,strand,annotation)
                            tobj.addExon(eobj)
                            eobj.addTranscript(tobj)
                        else:
                            annotation['exon_id'] = fixeid
                            eobj = self.getOrCreateExon(chr,start,stop,strand,annotation)
                    # pass     
        else:
            # print("Unknown type: "+type)
            self.notparsedtypes.add(type)

        
    def toDict(self, annotation, sep1, sep2):
        a = annotation.split(sep1)
        out = {}
        for elem in a:
            elem = elem.split(sep2)
            key = elem[0]
            value = elem[1]
            value = value.replace('"',"")
            out[key] = value
        return out
    
    def validateCoordinates(self,id,chr,start,stop,strand,obj, verbose):
        ok = True
        if obj.chr != chr:
            if verbose:
                print(id+" - chr found: {}, expected: {}".format(obj.chr,chr))
            ok = False
        if obj.start != start:
            if verbose:
                print(id+" - start found: {}, expected: {} - delta: {}".format(obj.start,start, (obj.start-start) ))
            ok = False
        if obj.stop != stop:
            if verbose:
                print(id+" - stop found: {}, expected: {} - delta: {}".format(obj.stop,stop, (obj.stop-stop)))
            ok = False
        if obj.strand != strand:
            if verbose:
                print(id+" - strand found: {}, expected: {}".format(obj.strand,strand))  
            ok = False
        return ok
    
    # function to check whether loading in the GTF created any inconsistencies, potentially caused by duplicate entries with different coordinates
    def validate(self, gtffile):
        fh = self.getfh(gtffile)
        print("Comparing {} to gene model loaded in memory...".format(gtffile))
        for line in fh:
            if not line.startswith("#"):
                elems = line.split("\t")
                chr = Chromosome.parse(elems[0])
                type = elems[2]
                start = int(elems[3])
                stop = int(elems[4])
                strand = Strand.parse(elems[6])
                annotation = {}
                if len(elems) > 8:
                    annotation = self.toDict(elems[8],"; "," ")
                if type == "gene":
                    gid = annotation.get("gene_id")
                    gobj = self.idToGene.get(gid)
                    if gobj is None:
                        print(gid+" not loaded")
                    else:
                        self.validateCoordinates(gid,chr,start,stop,strand,gobj, True)
                        symbol = annotation.get("gene_name")
                        if symbol is not None and gobj.symbol != symbol:
                            print(id+" - symbol found: {}, expected: {}".format(gobj.stop,stop))                                              
                elif type == "transcript":
                    tid = annotation.get("transcript_id")
                    tobj = self.idToTranscript.get(tid)
                    if tobj is None:
                        print(tid+" not loaded")
                    else:
                        self.validateCoordinates(tid,chr,start,stop,strand,tobj, True)
                        gid = annotation.get("gene_id")
                        gobj = self.idToGene.get(gid)
                        if gobj is None:
                            print(tid+" - "+gid+" not loaded")              
                        else:
                            if gobj.transcripts is None:
                                print(tid+" - "+gid+" not assigned to gene as transcript")
                            else:
                                found = False
                                for tobj2 in gobj.transcripts:
                                    if tobj2 == tobj:
                                        found = True
                                if not found:
                                    print(tid+" - "+gid+" not assigned to gene as transcript")
                            if gobj != tobj.gene:
                                print(tid+" - "+gid+" not linked to correct gene object")
                elif type == "exon":
                    eid = annotation.get("exon_id")
                    # Ensembl exon IDs on para-autosomal reqions do not have the _PAR_Y extension. Fix that here
                    gid = annotation.get("gene_id")
                    paraautosomal = False
                    if gid is not None:
                        if gid.endswith("_PAR_Y"):
                            paraautosomal = True
                    if paraautosomal:
                        eid = eid+"_PAR_Y"
                        annotation['exon_id'] = eid
                    eobj = self.idToExon.get(eid)
                    if eobj is None:
                        print(eid+" - not loaded")
                    else:
                        ok = self.validateCoordinates(eid,chr,start,stop,strand,eobj, False)

                        # could be a partially remapped exon; check whether I fixed this.
                        # create new exon object.
                        if not ok:
                            found = False
                            itr = 1
                            fixeid = eid+"_"+str(itr)
                            while self.idToExon.get(fixeid) is not None:
                                fixeobj = self.idToExon.get(fixeid)
                                if fixeobj is None:
                                    print(fixeid+" was never created as fix.")
                                else:
                                    ok = self.validateCoordinates(eid,chr,start,stop,strand,fixeobj, False)
                                    if ok:
                                        eobj = fixeobj
                                        found = True
                                        eid = fixeid
                                itr += 1
                                fixeid = eid+"_"+str(itr)
                            if not found:
                                print(eid+" is a broken exon without a fix :(")

                        tid = annotation.get("transcript_id")
                        tobj = self.idToTranscript.get(tid)
                        gid = annotation.get("gene_id")
                        gobj = self.idToGene.get(gid)
                        if gobj is None:
                            print(eid+" - "+tid+" not loaded")
                        else:
                            if gobj != eobj.gene:
                                print(eid+" - "+gid+" not linked to correct gene object")
                        if tobj is None:
                            print(eid+" - "+tid+" not loaded")
                        else:
                            # check whether exon is present in transcript
                            if tobj.exons is None:
                                print(eid+" - "+tid+" no exons loaded in transcript")
                            else:
                                found = False
                                for eobj2 in tobj.exons:
                                    if eobj2 == eobj:
                                        found = True
                                if not found:
                                    print(eid+" - "+tid+" exon not linked to transcript object")
                        if eobj.transcripts is None:
                            print(eid+" - "+tid+" no transcripts linked to exon")
                        else:
                            found = False
                            for tobj2 in eobj.transcripts:
                                if tobj2 == tobj:
                                    found = True
                            if not found:
                                print(eid+" - "+tid+" transcript not linked to exon")
        fh.close()