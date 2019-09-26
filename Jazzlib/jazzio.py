import io
from .Peak import *
from .Hotspot import *
from .Peak import *
from .FRegion import *

def peakbedswriter(samplename, peaks):

    bedfilename =samplename+ '_' + 'peak' + ".bed"

    open_bed = io.FileIO(bedfilename, 'w')

    for peak in peaks:

        #bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.hotspotid]
        bedlist = [str(peak.chromosome), str(peak.start), str(peak.end),str(peak.peakid),str(peak.score)]

        linker = "\t"

        outstring = linker.join(bedlist) + "\n"

        open_bed.write(bytes(outstring, encoding = 'utf-8'))

    open_bed.close()


def peakbedgraphswriter(samplename, peaks):

    bedfilename =samplename+ '_' + 'peak' + ".bedgraph"

    open_bed = io.FileIO(bedfilename, 'w')

    for peak in peaks:

        #bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.hotspotid]
        bedlist = [str(peak.chromosome), str(peak.start), str(peak.end),str(peak.score)]

        linker = "\t"

        outstring = linker.join(bedlist) + "\n"

        open_bed.write(bytes(outstring, encoding = 'utf-8'))

    open_bed.close()


def hotspotsbedswriter(samplename, hotspots):

    bedfilename =samplename+ '_' + 'hotspots' + ".bed"

    open_bed = io.FileIO(bedfilename, 'w')

    for hotspot in hotspots:

        #bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.hotspotid]
        bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end)]

        linker = "\t"

        outstring = linker.join(bedlist) + "\n"

        open_bed.write(bytes(outstring, encoding = 'utf-8'))

    open_bed.close()


def hotpeakbedswriter2(samplename, hotspots):

    bedfilename =samplename+ '_' + 'peaks' + ".bed"

    open_bed = io.FileIO(bedfilename, 'w')

    for hotspot in hotspots:

        for peak in hotspot.peaks:

            bedlist = [str(peak.chromosome), str(peak.start), str(peak.end),str(peak.peakid),str(peak.score)]

            linker = "\t"

            outstring = linker.join(bedlist) + "\n"

            open_bed.write(bytes(outstring, encoding = 'utf-8'))

    open_bed.close()


def jazzgffout(samplename, hotspots, peaks, fregion):

    bedfilename =samplename+ '_' + 'peaks_hotspots' + ".gff3"

    open_bed = io.FileIO(bedfilename, 'w')
    linker = "\t"

    frsite = dict()

    for fr in fregion.filted_region:

        (frchrnow,frstartend) = fr.split(":")

        (frstart,frend) = frstartend.split("-")

        for sitenow in range(int(frstart), int(frend)+1):

            if frchrnow in frsite:

                frsite[frchrnow][sitenow] = 1

            else:

                frsite[frchrnow] = dict()

                frsite[frchrnow][sitenow] = 1

    hotspotsinfr = dict()

    for hotspot in hotspots:



        if hotspot.chromosome in frsite:

            for nowsite in range(hotspot.start, hotspot.end+1):

                if nowsite in frsite[frchrnow]:

                    hotspotanno = "ID="+str(hotspot.hotspotid)+";anno=FREGION"

                    hotspotsinfr[hotspot.hotspotid] = 1

                else:

                    hotspotanno = "ID="+str(hotspot.hotspotid)

        else:

            hotspotanno = "ID="+str(hotspot.hotspotid)

        hotspotsstr = [str(hotspot.chromosome), "JAZZ", "gene", str(hotspot.start), str(hotspot.end),
                        '.', '.', '.',hotspotanno
                      ]

        hotspotstring = linker.join(hotspotsstr) + "\n"

        # open_bed.write(hotspotstring)
        open_bed.write(bytes(hotspotstring, encoding='utf-8'))

    for peak in peaks:

        # peakanno = "Parent="+str(peak.parent)+";"+"ID="+str(peak.peakid)

        if peak.parent in hotspotsinfr:
            peakanno = "Parent="+str(peak.parent)+";"+"ID="+str(peak.peakid)+";anno=FREGION"
        else:
            peakanno = "Parent="+str(peak.parent)+";"+"ID="+str(peak.peakid)

        peakstr = [str(peak.chromosome), "JAZZ", "CDS", str(peak.start), str(peak.end),
                                '.', '.', '.',peakanno]

        peakstring = linker.join(peakstr)+"\n"

        #open_bed.write(peakstring)
        open_bed.write(bytes(peakstring, encoding='utf-8'))


    open_bed.close()
