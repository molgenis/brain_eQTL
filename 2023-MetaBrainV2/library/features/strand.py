from enum import Enum

class Strand(Enum):        
    POS = 1
    NEG = -1
    NA = "NA"

    @classmethod
    def parse(self, chrstr):
        chrstr = chrstr.lower()
        if chrstr == "+" or chrstr == "1":
            return Strand.POS
        if chrstr == "-" or chrstr == "-1":
            return Strand.NEG
        return Strand.NA
    
    def toStr(str):
        if str == Strand.NA:
            return "NA"
        if str == Strand.POS:
            return "+"
        return "-"    