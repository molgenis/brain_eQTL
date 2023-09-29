class Feature:
    chr = None
    start = None
    stop = None
    strand = None

    def __init__(self, chr, start, stop, strand):
        self.chr = chr
        self.start = start
        self.stop = stop
        self.strand = strand

    def overlaps(self, other):
        if self.chr != other.chr:
            return False
        if self.start > other.stop:
            return False
        if self.stop < other.start:
            return False
        return True

    def overlapsStrand(self, other):
        ov = self.overlap(other)
        if ov:
            return self.strand == other.strand

    def bpOverlap(self, other):
        ov = self.overlaps(other)
        if not ov:
            return 0
        
        # |------| self
        #   |---|  other
        if self.start <= other.start and self.stop >= other.stop:
            return other.stop - other.start
        #   |---|  self
        # |------| other
        if other.start <= self.start and other.stop >= self.stop:
            return self.stop - self.start
        
        # |---| self
        #   |---|  other
        if self.stop <= other.stop and self.start <= other.start:
            return self.stop - other.start
        #      |---| self
        #   |---|    other
        if self.start <= other.stop and self.start >= other.start:
            return other.stop - self.start
        
    def absoluteMinimalDistance(self, feat2):
        mind = abs(self.start - feat2.start)
        d = abs(self.start - feat2.stop)
        if d < mind:
            mind = d
        d = abs(self.stop - feat2.start)
        if d < mind:
            mind = d
        d = abs(self.stop - feat2.stop)
        if d < mind:
            mind = d
        return mind