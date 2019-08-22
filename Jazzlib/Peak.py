class Peak:
    """
        Peaks
    
    """    
   
    def __init__(self, start, end, chromosome, peakpoint, peakid, score, parent=1, fdr=1):

        self.start = start

        self.end = end

        self.chromosome = chromosome

        self.peakpoint = peakpoint

        self.peakid = peakid

        self.score = score

        self.fdr = fdr

        self.parent = parent

