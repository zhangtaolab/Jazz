class Hotspot:
    """
        Hotspot

    """
   
    def __init__(self, start, end, chromosome, hotspotid, peaks=list(), score=0, fdr=1):

        self.start = start

        self.end = end

        self.chromosome = chromosome

        self.hotspotid = hotspotid

        self.score = score

        self.fdr = fdr

        self.peaks = peaks

    def addpeak(self, peak):

        self.peaks.append(peak)
