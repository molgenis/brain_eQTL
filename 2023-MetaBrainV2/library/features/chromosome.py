from enum import Enum

class Chromosome(Enum):
    
    Chr1 = 1
    Chr2 = 2
    Chr3 = 3
    Chr4 = 4
    Chr5 = 5
    Chr6 = 6    
    Chr7 = 7  
    Chr8 = 8    
    Chr9 = 9    
    Chr10 = 10    
    Chr11 = 11   
    Chr12 = 12    
    Chr13 = 13   
    Chr14 = 14   
    Chr15 = 15   
    Chr16 = 16   
    Chr17 = 17   
    Chr18 = 18   
    Chr19 = 19    
    Chr20 = 20   
    Chr21 = 21
    Chr22 = 22    
    ChrX = 23    
    ChrY = 24    
    ChrMT = 25    
    ChrXY = 26   
    NA = -1

    @classmethod
    def parse(self, chrstr):
        chrstr = chrstr.lower()
        chrstr = chrstr.replace("chr","")
        if chrstr == "m" or chrstr == "mt":
            return Chromosome.ChrMT
        if chrstr == "x":
            return Chromosome.ChrX
        if chrstr == "xy":
            return Chromosome.ChrXY
        if chrstr == "y":
            return Chromosome.ChrY
        if chrstr == "xy":
            return Chromosome.ChrXY
        try:
            chrint = int(chrstr)
            return Chromosome(chrint)
        except:
            return Chromosome.NA  

    def getNumber(self):
        return self.value  
    
    def getName(self):
        return self.name
    
    def isAutosomal(self):
        return (self.value > 0 and self.value < 23)