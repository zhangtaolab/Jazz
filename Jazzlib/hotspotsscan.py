

from .countreads import *
from .cEM_zip import *
from .FRegion import *
from multiprocessing import Pool
from .Hotspot import *
from .sta import *
from .region import *
from .Peak import *


class KeyboardInterruptError(Exception):

    pass


def hotspotsscan_withoutcontrol(file, maxinsert, windowscare,countchr,gloablumbda,
                                bayesfactorthreshold, nthreads, fregion, jobtype):

    pool = Pool(nthreads)

    try:

        pars = list()

        hotspots = list()

        print ("gloablumbda",gloablumbda , "readlengthmean", fregion.readlengthmean)

        bayesfactorthresholdcount = 2

        i = 2

        while True:

            nowbayesfactor = bayesfactor(gloablumbda, i)

            if nowbayesfactor > bayesfactorthreshold:

                break

            bayesfactorthresholdcount = i

            i = i + 1

        print ("bayesfactorthresholdcount", bayesfactorthresholdcount)

        windowsize = 100000

        for chromosmoe in countchr:

            chr_length = fregion.chrs_length[chromosmoe]

            for scare in range(0, int(chr_length/windowsize)+1):

                nowstart = scare*windowsize + 1 -200

                nowend = (scare+1)*windowsize + 200

                if nowend > chr_length:

                    nowend = chr_length

                if nowstart < 1:

                    nowstart = 1

                nowregion = chromosmoe + ":" + str(nowstart) + "-" + str(nowend)

                par = dict()

                par['region'] = nowregion

                par['maxinsert'] = maxinsert

                par['bamfile'] = file

                par['jobtype'] = jobtype

                par['chrlength'] = chr_length

                par['regionchromosome'] = chromosmoe

                par['regionstart'] = nowstart

                par['regionend'] = nowend

                # par['bayesfactordic'] = bayesfactordic

                par['bayesfactorcount'] = bayesfactorthresholdcount

                par['readlengthmean'] = fregion.readlengthmean

                pars.append(par)

        enrichedinthreads = pool.map(hotspot_withoutcontrol_worker, pars)

        chrenrichedpotin = dict()

        for enrichedinthread in enrichedinthreads:

            nowchr = enrichedinthread['chromosome']

            if nowchr in chrenrichedpotin:

                chrenrichedpotin[nowchr].append(enrichedinthread['list'])

            else:

                chrenrichedpotin[nowchr] = list()

                chrenrichedpotin[nowchr].append(enrichedinthread['list'])

        chrhotpars = list()

        for nowchr in chrenrichedpotin:

            hotpar = dict()

            hotpar['chromosome'] = nowchr

            hotpar['preregion'] = chrenrichedpotin[nowchr]

            hotpar['chr_length'] = fregion.chrs_length[chromosmoe]

            hotpar['fregion'] = fregion

            chrhotpars.append(hotpar)

        hotsptosinthreads = pool.map(hotspots_chromsome_merge,chrhotpars)

        for hotinth in hotsptosinthreads:

            for hotspotnow in hotinth:

                hotspots.append(hotspotnow)

        pool.close()

        return hotspots

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

        print ('got exception in Jazzlib.hotspotsscan.hotspotsscan_withoutcontrol: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:

        pool.join()


def hotspot_withoutcontrol_worker(par):

    try:

        maxinsert = par['maxinsert']

        bamfile = par['bamfile']

        jobtype = par['jobtype']

        chromosome = par['regionchromosome']

        nowstart = par['regionstart']

        nowend = par['regionend']

        bayesfactorcount = par['bayesfactorcount']

        readlengthmean = par['readlengthmean']

        datacount = extenddepthcount(bamfile=bamfile, regionchromosome=chromosome, regionstart=nowstart,
                                     regionend=nowend, maxinsert=maxinsert, jobtype=jobtype,
                                     readlengthmean=readlengthmean)

        enrichedlist = dict()

        enrichedlist['chromosome'] = chromosome

        enrichedlist['list'] = list()

        for site in datacount:

            if datacount[site] >= bayesfactorcount:

                enrichedlist['list'].append(site)

        return enrichedlist

    except Exception as e:

        print ('got exception in Jazzlib.hotspotsscan.hotspot_withoutcontrol_worker: %r, terminating the pool' % (e,))

        print ('pool is terminated')

    except KeyboardInterrupt:

         print ("You cancelled the program!")

         sys.exit(1)


def hotspotsscan_withcontrol(chipfile, maxinsert, windowscare,countchr,inputgloablumbda,
                             bayesfactorthreshold, nthreads, chipfregion, jobtype, ratio,
                             inputfile, inputfregion):

    pool = Pool(nthreads)

    try:

        pars = list()

        hotspots = list()

        print ("gloablumbda",inputgloablumbda , "readlengthmean", inputfregion.readlengthmean)

        bayesfactorthresholdcount = 2

        i = 2

        while True:

            nowbayesfactor = bayesfactor(inputgloablumbda, i)

            if nowbayesfactor > bayesfactorthreshold:

                break

            bayesfactorthresholdcount = i

            i = i + 1

        print ("bayesfactorthresholdcount", bayesfactorthresholdcount)

        windowsize = 100000

        for chromosmoe in countchr:

            chr_length = chipfregion.chrs_length[chromosmoe]

            for scare in range(0, int(chr_length/windowsize)+1):

                nowstart = scare*windowsize + 1 -200

                nowend = (scare+1)*windowsize + 200

                if nowend > chr_length:

                    nowend = chr_length

                if nowstart < 1:

                    nowstart = 1

                nowregion = chromosmoe + ":" + str(nowstart) + "-" + str(nowend)

                par = dict()

                par['region'] = nowregion

                par['maxinsert'] = maxinsert

                par['bamfile'] = chipfile

                par['jobtype'] = jobtype

                par['chrlength'] = chr_length

                par['regionchromosome'] = chromosmoe

                par['regionstart'] = nowstart

                par['regionend'] = nowend

                par['ratio'] = ratio

                # par['bayesfactordic'] = bayesfactordic

                par['bayesfactorcount'] = bayesfactorthresholdcount

                par['readlengthmean'] = chipfregion.readlengthmean

                pars.append(par)

        enrichedinthreads = pool.map(hotspot_control_worker, pars)

        chrenrichedpotin = dict()

        for enrichedinthread in enrichedinthreads:

            nowchr = enrichedinthread['chromosome']

            if nowchr in chrenrichedpotin:

                chrenrichedpotin[nowchr].append(enrichedinthread['list'])

            else:

                chrenrichedpotin[nowchr] = list()

                chrenrichedpotin[nowchr].append(enrichedinthread['list'])

        chrhotpars = list()

        for nowchr in chrenrichedpotin:

            hotpar = dict()

            hotpar['chromosome'] = nowchr

            hotpar['preregion'] = chrenrichedpotin[nowchr]

            hotpar['chr_length'] = chipfregion.chrs_length[chromosmoe]

            hotpar['fregion'] = chipfregion

            chrhotpars.append(hotpar)

        hotsptosinthreads = pool.map(hotspots_chromsome_merge, chrhotpars)

        for hotinth in hotsptosinthreads:

            for hotspotnow in hotinth:

                hotspots.append(hotspotnow)

        pool.close()

        pool.close()

        return hotspots

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

        print ('got exception in Jazzlib.hotspotsscan.hotspotsscan_withcontrol: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:

        pool.join()


def hotspot_control_worker(par):

    try:

        maxinsert = par['maxinsert']

        bamfile = par['bamfile']

        jobtype = par['jobtype']

        chromosome = par['regionchromosome']

        nowstart = par['regionstart']

        nowend = par['regionend']

        bayesfactorcount = par['bayesfactorcount']

        readlengthmean = par['readlengthmean']

        ratio = par['ratio']

        datacount = extenddepthcount(bamfile=bamfile, regionchromosome=chromosome, regionstart=nowstart,
                                     regionend=nowend, maxinsert=maxinsert, jobtype=jobtype,
                                     readlengthmean=readlengthmean)

        enrichedlist = dict()

        enrichedlist['chromosome'] = chromosome

        enrichedlist['list'] = list()

        for site in datacount:

            if datacount[site]*ratio >= bayesfactorcount:

                enrichedlist['list'].append(site)

        return enrichedlist

    except Exception as e:

        print ('got exception in Jazzlib.hotspotsscan.hotspot_withoutcontrol_worker: %r, terminating the pool' % (e,))

        print ('pool is terminated')

    except KeyboardInterrupt:

         print ("You cancelled the program!")

         sys.exit(1)



def hotspotsfilter(hotspots, peaks):

    peaksparent = dict()

    for peak in peaks:

        if peak.parent not in peaksparent:

            peaksparent[peak.parent] = 1

    hotspotreturen = list()

    for hotspot in hotspots:

        if hotspot.hotspotid in peaksparent:

            hotspotreturen.append(hotspot)

    return hotspotreturen


def hotspots_chromsome_merge(par):

    try:

        chromosome = par['chromosome']

        preregion = par['preregion']

        chr_length = par['chr_length']

        fregion = par['fregion']

        hotspotslist = list()

        enrichedpotin = dict()

        for regionpoint in preregion:

            for nowsite in regionpoint:

                if not nowsite in enrichedpotin:

                    enrichedpotin[nowsite] = 1

        chrenrichlist = list(enrichedpotin.keys())

        temphotspots = continueregion(chrenrichlist, 2)

        for hotspotstarend in temphotspots:

            hotspotstart = hotspotstarend['start_site']

            hotspotend = hotspotstarend['end_site']

            if hotspotend-hotspotstart < fregion.readlengthmean/2:

                continue

            hotspotid = str(chromosome) + ":" + str(hotspotstart) +"-"+ str(hotspotend)

            hotspot = Hotspot(start=hotspotstart, end=hotspotend, chromosome=chromosome, hotspotid=hotspotid)

            hotspotslist.append(hotspot)

        return hotspotslist

    except Exception as e:

        print ('got exception in Jazzlib.hotspotsscan.hotspots_chromsome_merge: %r, terminating the pool' % (e,))

        print (par)

        print ('pool is terminated')

    except KeyboardInterrupt:

         print ("You cancelled the program!")

         sys.exit(1)