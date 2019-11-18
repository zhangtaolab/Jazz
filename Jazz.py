

import os
import sys
from optparse import OptionParser
import logging
from Jazzlib.FRegion import *
from Jazzlib.localmax import *
from Jazzlib.normalize_ratio import *
from Jazzlib.countreads import *
from Jazzlib.Peak import *
from Jazzlib.sta import *
from Jazzlib.jazzio import *
from Jazzlib.randombg import *
from Jazzlib.hotspotsscan import *
from Jazzlib.Hotspot import *
from Jazzlib.peaksscan import *

def main():

    opt = opt_check(get_optparser())

    if opt.controlfile == "no":

        nocontrol(opt)

    else:

        withcontrol(opt)


def withcontrol(opt):

    try:

        datafile = opt.datafile

        inputfile = opt.controlfile

        jobtype = opt.jobtype

        count_chr = opt.countchr

        maxinsert = opt.maxinsert

        nthreads = opt.nthreads

        bayesfactorthreshold = opt.threshold

        # bayesfactorthreshold = 10

        samplename = opt.samplename

        fdr = opt.fdr


        chipfregion = FRegion(bamfile=datafile, jobtype=jobtype, countchr=count_chr, nthreads=nthreads, maxinsert=maxinsert)

        inputfregion = FRegion(bamfile=inputfile, jobtype=jobtype, countchr=count_chr, nthreads=nthreads, maxinsert=maxinsert)

        ratio = normalize_ratio_input2(fregegion_input=inputfregion, fregion_chip=chipfregion)

        if opt.genomesize:

            print("###chipfregion.adjreads * chipfregion.readlengthmean/chipfregion.countgenomelength,",
                  chipfregion.adjreads, chipfregion.readlengthmean, opt.genomesize)

            gloablumbda = chipfregion.adjreads * chipfregion.readlengthmean / opt.genomesize

        else:

            print("###inputfregion.adjreads,inputfregion.readlengthmean,inputfregion.countgenomelength",
                  inputfregion.adjreads , inputfregion.readlengthmean,inputfregion.countgenomelength)

            gloablumbda = inputfregion.adjreads * inputfregion.readlengthmean/inputfregion.countgenomelength

        windowscare=100000

        hotspots = hotspotsscan_withcontrol(chipfile=datafile,maxinsert=maxinsert, windowscare=windowscare,
                                            countchr=count_chr, inputgloablumbda=gloablumbda,
                                            bayesfactorthreshold=bayesfactorthreshold, nthreads=nthreads,
                                            chipfregion=chipfregion, jobtype=jobtype, ratio=ratio, inputfile=inputfile,
                                            inputfregion=inputfregion)

        peaks = peakscan_control(datafile=datafile,maxinsert=maxinsert, bayesfactorthreshold=bayesfactorthreshold,
                                 nthreads=nthreads,chipfregion=chipfregion, jobtype=jobtype, hotspots=hotspots,
                                 gloablumbda=gloablumbda,inputfile=inputfile,ratio=ratio,inputfregion=inputfregion)

        hotspotsenrich = hotspotsfilter(hotspots=hotspots, peaks=peaks)

        hotspotsbedswriter(hotspots=hotspotsenrich, samplename=samplename)

        peakbedswriter(samplename=samplename,peaks=peaks)

        jazzgffout(samplename=samplename, hotspots=hotspotsenrich, peaks=peaks, fregion=fregion)

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)


def nocontrol(opt):

    try:

        datafile = opt.datafile

        jobtype = opt.jobtype

        count_chr = opt.countchr

        maxinsert = opt.maxinsert

        print ("maxinsert",maxinsert)

        nthreads = opt.nthreads

        bayesfactorthreshold = opt.threshold

        samplename = opt.samplename

        chipfregion = FRegion(bamfile=datafile, jobtype=jobtype, countchr=count_chr, nthreads=nthreads,
                              maxinsert=maxinsert)

        if opt.genomesize:

            print("###chipfregion.adjreads * chipfregion.readlengthmean/chipfregion.countgenomelength,", chipfregion.adjreads, chipfregion.readlengthmean,opt.genomesize)

            gloablumbda = chipfregion.adjreads * chipfregion.readlengthmean / opt.genomesize

        else:

            print("###chipfregion.adjreads * chipfregion.readlengthmean/chipfregion.countgenomelength,", chipfregion.adjreads, chipfregion.readlengthmean,chipfregion.countgenomelength)

            gloablumbda = chipfregion.adjreads * chipfregion.readlengthmean/chipfregion.countgenomelength

        windowscare=100000

        for fregions in chipfregion.filted_region:

            print (fregions)

        hotspots = hotspotsscan_withoutcontrol(file=datafile, maxinsert=maxinsert, windowscare=windowscare, countchr=count_chr,
                                               bayesfactorthreshold=bayesfactorthreshold, nthreads=nthreads,
                                               fregion=chipfregion, jobtype=jobtype, gloablumbda=gloablumbda)

        peaks = peakscan_without_control(datafile=datafile,maxinsert=maxinsert,
                                         bayesfactorthreshold=bayesfactorthreshold, nthreads=nthreads,
                                         fregion=chipfregion,jobtype=jobtype,
                                         hotspots=hotspots, gloablumbda=gloablumbda)

        hotspotsenrich = hotspotsfilter(hotspots=hotspots, peaks=peaks)

        hotspotsbedswriter(hotspots=hotspotsenrich, samplename=samplename)

        peakbedswriter(samplename=samplename,peaks=peaks)

        jazzgffout(samplename=samplename, hotspots=hotspotsenrich, peaks=peaks, fregion=chipfregion)

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)


def get_optparser():

    usage = """usage: %prog <-d datafile> [-n name] [options]
    Example %prog -i nh_sample1.bam -n sample1
    """

    description = "%prog Non-Histone protein banding site identification"

    jazzopt = OptionParser(version="%prog 0.1 20140521", description=description, usage=usage, add_help_option=False)

    jazzopt.add_option("-h", "--help", action="help", help="show this help message and exit.")

    jazzopt.add_option("-d", "--data", dest="datafile", type="string", help='data file, should be sorted bam format')

    jazzopt.add_option("-c", "--control", dest="controlfile", type="string", help='control(input) file, should be sorted bam format', default="no")

    jazzopt.add_option("-n", "--name", dest="samplename", help="NH sample name default=NH_sample", type="string" , default="DH_sample")

    jazzopt.add_option("-t", "--threshold", dest="threshold", type="float", help="peak threshold, default=6.0", default=6.0)

    jazzopt.add_option("--threads", dest="nthreads", type="int", help="threads number or cpu number, default=4", default=4)

    jazzopt.add_option("-w", "--wig", action="store_true", help="whether out put wiggle file, default=False", default=False)

    jazzopt.add_option("-f","--fdr", dest="fdr", type="float",help="using FDR as threshold", default=0.1)

    jazzopt.add_option("-x", "--excludechr", dest="excludechr", help="Don't count those chromosome, strongly suggest skip mitochondrion and chloroplast, example='-x ChrM,ChrC'")

    jazzopt.add_option("-g", "--gff", action="store_true", help="whether out put gff file, default=False", default=False)

    jazzopt.add_option("-j","--jobtype",dest="jobtype",type="string",help="job type, such as nhpaired or nhsingle")

    jazzopt.add_option("-m","--maxinsert",dest="maxinsert",type="int",help="when you use paired library, please set the maxinsert size",default=130)

    jazzopt.add_option("--pe", dest="pe", action="store_true", help="paired-end reads or single-end reads, default=False (single end)", default=False)

    jazzopt.add_option("--genomesize", dest="genomesize", action="int",
                       help="Set genome size", default=False)

    return jazzopt


def opt_check(jazzopt):

    (opt, args) = jazzopt.parse_args()

    if not opt.datafile:

        logging.error("you need input a bam file, '-d nh_sample1.bam -j nhsingle'")

        jazzopt.print_help()

        sys.exit(1)

    if not os.path.isfile (opt.datafile):

        logging.error("No such file: %s" % opt.datafile)

        sys.exit(1)

    dataindexfile1 = opt.datafile + '.bai'

    dataindexfile2 = opt.datafile + '.csi'

    if not (os.path.isfile(dataindexfile1) or os.path.isfile(dataindexfile2)):

        logging.error("Missing bam index file: %s or %s" % (dataindexfile1, dataindexfile2))

        sys.exit(1)

    if not opt.controlfile == "no":

        if not os.path.isfile (opt.controlfile):

            logging.error("No such file: %s" % opt.controlfile)

            sys.exit(1)

        controlindexfile1 = opt.controlfile + '.bai'

        controlindexfile2 = opt.controlfile + '.csi'

        if not (os.path.isfile(controlindexfile1) or os.path.isfile(controlindexfile2)):

            logging.error("Missing bam index file: %s or %s" % (controlindexfile1, controlindexfile2))

            sys.exit(1)

    else:

        opt.controlfile = "no"

    if not (opt.nthreads > 0):

        logging.error("threads number should >=1")

        jazzopt.print_help()

        sys.exit(1)

    if (opt.jobtype):

        if opt.jobtype == 'nhsingle':

            if (opt.maxinsert < 0):

                logging.error("maxinsert size error")

                jazzopt.print_help()

                sys.exit(1)

        elif opt.jobtype == 'nhpaired':

            if (opt.maxinsert < 0):

                logging.error("maxinsert size error")

                jazzopt.print_help()

                sys.exit(1)

        else:

            logging.error("missing or wrong jobtype")

            jazzopt.print_help()

            sys.exit(1)

    else:

        logging.error("missing or wrong jobtype")

        jazzopt.print_help()

        sys.exit(1)

    opt.countchr = list()

    samfile = pysam.Samfile(opt.datafile)

    sam_ref = samfile.references

    for i in sam_ref:

        opt.countchr.append(i)

    if (opt.excludechr):

        excludchr = opt.excludechr.split(',')

        for chri in excludchr:

            if not chri in sam_ref:

                print (chri,'not in the %s file' % opt.datafile)

                print ("try to selcet exclude Chr from", end =" : ")

                print (sam_ref, sep=",")

                jazzopt.print_help()

                sys.exit(1)

            else:

                j = 0

                for n in opt.countchr:

                    if chri == n:

                        del opt.countchr[j]

                    j = j + 1

    return opt

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)

