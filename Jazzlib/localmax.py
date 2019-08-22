

from .kernelsmooth import *
from multiprocessing import Pool
from .kernel import *
from .FRegion import *


class KeyboardInterruptError(Exception):

    pass


def get_all_localmax(bamfile, jobtype, maxinsert, nthreads, fregion, countchr, rndth):

    pool = Pool(nthreads)

    try:

        pars = list()

        windowsize = 100000

        adjreads = fregion.adjreads

        totallength = 0

        onesmoothkernel = smooth_kernel(30)

        kermax = max(onesmoothkernel.values())
        #


        for chromosmoe in countchr:

            chr_length = fregion.chrs_length[chromosmoe]

            totallength = totallength + chr_length

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

                par['bamfile'] = bamfile

                par['jobtype'] = jobtype

                par['chrlength'] = chr_length

                par['regionchromosome'] = chromosmoe

                par['regionstart'] = nowstart

                par['regionend'] = nowend

                par['rndth'] = rndth

                pars.append(par)

        avgcount = adjreads/totallength

        threshhold = int(avgcount + 1) * kermax
        ###test threhhold
        #threshhold = avgcount

        # print ("threshhold:", threshhold)

        filted_region = fregion.filted_region

        filted_site = dict()

        for fr in filted_region:

            chromosome, sesite = fr.split(':')

            startsite, endsite = sesite.split('-')

            startsite = int(startsite)

            endsite = int(endsite)

            if chromosome in filted_site:

                for i in range(startsite,endsite):

                    filted_site[chromosome][i] = 1

            else:

                filted_site[chromosome] = dict()

                for i in range(startsite,endsite):

                    filted_site[chromosome][i] = 1

        localmax = dict()

        localmax_worker_returnres = pool.map(localmax_worker, pars)

        for each_worker_res in localmax_worker_returnres:

            for chromosome in each_worker_res:

                for site in each_worker_res[chromosome]:

                    if chromosome in localmax:

                        if each_worker_res[chromosome][site] > threshhold:

                            if chromosome in filted_site:

                                if site in filted_site[chromosome]:

                                    continue

                            localmax[chromosome][site] = each_worker_res[chromosome][site]

                    else:

                        if each_worker_res[chromosome][site]>threshhold:

                            if chromosome in filted_site:

                                if site in filted_site[chromosome]:

                                    continue

                            localmax[chromosome] = dict()

                            localmax[chromosome][site] = each_worker_res[chromosome][site]

        pool.close()

        # print (localmax)

        return localmax

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

        print ('got exception in Jazzlib.localmax.get_all_localmax: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:
        #     print ('joining pool processes')
        pool.join()
            # print ('join complete')


def localmax_worker(par):

    try:

        nowregion = par['region']

        maxinsert = par['maxinsert']

        bamfile = par['bamfile']

        jobtype = par['jobtype']

        chr_length = par['chrlength']

        regionchromosome = par['regionchromosome']

        regionstart = par['regionstart']

        regionend = par['regionend']

        rndth = par['rndth']

        # smoothedscore = regionsmooth(bamfile=bamfile, maxinsert=maxinsert, jobtype=jobtype, region=nowregion,
        #                              chr_length=chr_length)

        smoothedscore = regionsmooth(bamfile=bamfile, maxinsert=maxinsert, jobtype=jobtype,
                                     regionchromosome=regionchromosome,
                                     regionstart=regionstart, regionend=regionend,
                                     chr_length=chr_length)

        localmax = smoothedlocalmax(smoothedscore, rndth)

        return localmax

    except KeyboardInterrupt:

        raise KeyboardInterruptError()

    except Exception as e:

        print ('got exception in Jazzlib.localmax.localmax_worker: %r,' % (e,))


def smoothedlocalmax(smoothedscore, rndth):

    try:

        maxsites = dict()

        startsite = min(smoothedscore['score'].keys())

        endsite = max(smoothedscore['score'].keys())

        chromosome = smoothedscore['chromosome']

        maxsites[chromosome] = dict()

        for nowsite in range(startsite+2, endsite-2):

            if smoothedscore['score'][nowsite] >=rndth:

                if (smoothedscore['score'][nowsite-2]<smoothedscore['score'][nowsite-1]<=smoothedscore['score'][nowsite]>=smoothedscore['score'][nowsite+1]>smoothedscore['score'][nowsite+2]):

                    maxsites[chromosome][nowsite] = smoothedscore['score'][nowsite]

                # print (nowsite)

        return maxsites

    except KeyboardInterrupt:

        raise KeyboardInterruptError()

    except Exception as e:

        print ('got exception in Jazzlib.localmax.smoothedlocalmax: %r,' % (e,))

