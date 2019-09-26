

from numpy import *
from .countreads import *
from multiprocessing import Pool
from .countreads import *
import timeit
import sys


class KeyboardInterruptError(Exception):

    pass


class FRegion:

    def __init__(self, bamfile, nthreads, maxinsert, jobtype, countchr=[]):

        self.bamfile = bamfile

        self.count_chr = countchr

        self.nthreads = nthreads

        self.maxinsert = maxinsert

        self.jobtype = jobtype

        self.__filte_region()

    def filte_region(self):

        bam_file = self.bamfile

        count_chr = self.count_chr

        nthreads = self.nthreads

        jobtype = self.jobtype

        maxinsert = self.maxinsert

        res = filter_region(bamfile=bam_file, count_chr=count_chr, nthreads=nthreads, maxinsert=maxinsert,
                            jobtype=jobtype)

        filted_region = res['filted_region']

        thresh_hold = res['thresh_hold']

        scare_std = res['region_std']

        scare_mean = res['region_mean']

        chr_total_reads = res['chr_total_reads']

        chrs_length = res['chrs_length']

        chrsfrcount = res['chrfrcount']

        filterreadscount = res['filterreadscount']

        totalreads = res['totalreads']

        chruniqlength = res['chruniqlength']

        readlengthmean = res['readlengthmean']

        adjreads = totalreads - filterreadscount

        countgenomelength = 0

        countgenomeuniqlength = 0

        for chromosome in count_chr:

            countgenomelength = countgenomelength + int(chrs_length[chromosome])

            countgenomeuniqlength = countgenomeuniqlength + int(chruniqlength[chromosome])

        self.countgenomelength = countgenomelength

        self.filted_region = filted_region

        self.thresh_hold = thresh_hold

        self.region_std = scare_std

        self.region_mean = scare_mean

        self.chr_total_reads = chr_total_reads

        self.chrs_length = chrs_length

        self.chrsfcount = chrsfrcount

        self.totalreads = totalreads

        self.filterreadscount = filterreadscount

        self.adjreads = adjreads

        self.chruniqlength = chruniqlength

        self.countgenomeuniqlength = countgenomeuniqlength

        self.readlengthmean = readlengthmean

    __filte_region = filte_region


def filter_region(bamfile, count_chr, nthreads, maxinsert, jobtype):

    pool = Pool(nthreads)

    try:

        samfile = pysam.Samfile(bamfile)

        windowsize = 1000

        totalreads = 0

        refere_ncenumber = samfile.nreferences

        ref_lengths = samfile.lengths

        sam_ref = samfile.references

        chrs_length = dict()

        chr_total_reads = dict()

        pars = list()

        chruniqlength = dict()

        chrreadlengthmean = dict()

        for chromosome in count_chr:

            for i in range(refere_ncenumber):

                if sam_ref[i] == chromosome:

                    chr_length = ref_lengths[i]

                    chrs_length[chromosome] = chr_length

                    chrcount = windowcounter(bamfile=bamfile, regionchromosome=chromosome,
                                             regionstart=1, regionend=int(chr_length),
                                             maxinsert=maxinsert,
                                             jobtype=jobtype)

                    chr_total_reads[chromosome] = chrcount

                    totalreads = chrcount + totalreads

        for chromosome in chrs_length:

            par = dict()

            par['chrmosome'] = chromosome

            par['windowsize'] = windowsize

            par['chr_length'] = chrs_length[chromosome]

            par['bamfile'] = bamfile

            par['maxinsert'] = maxinsert

            par['jobtype'] = jobtype

            pars.append(par)

        windowcountlist = list()

        windowregionlist = list()

        chrswindow = pool.map(chrwindow_counter, pars)

        for nowchrcount in chrswindow:

            nowchromosome = nowchrcount['chromosome']

            nowchromosome = str(nowchromosome)

            nowwindowcount = nowchrcount['windowcount']

            nowuniqcount = nowchrcount['uniqcount']

            nowreadslengthmean = nowchrcount['readlengthmean']

            print(nowchromosome, nowreadslengthmean)

            chrreadlengthmean[nowchromosome] = nowreadslengthmean

            chruniqlength[nowchromosome] = nowuniqcount

            for nowscare in nowwindowcount:

                nowstart = nowscare * windowsize + 1

                nowend = (nowscare+1) * windowsize

                if nowend > chrs_length[nowchromosome]:

                    nowend = chrs_length[nowchromosome]

                nowregion = nowchromosome+":"+str(nowstart)+"-"+str(nowend)

                windowcountlist.append(nowwindowcount[nowscare])

                windowregionlist.append(nowregion)

        scare_mean = mean(windowcountlist)

        scare_std = std(windowcountlist)

        print ("mean:", scare_mean, "std",scare_std)

        thresh_hold = scare_mean + 10 * scare_std

        chrsfrcount = 0

        filterreadscount = 0

        filted_region = list()

        for i in range(0, len(windowcountlist)):

            if windowcountlist[i] >= thresh_hold:

                # print (windowregionlist[i]," reads count ", windowcountlist[i])

                filted_region.append(windowregionlist[i])

                filterreadscount = filterreadscount + windowcountlist[i]

        res = dict()

        res['filted_region'] = filted_region

        res['thresh_hold'] = thresh_hold

        res['region_std'] = scare_std

        res['region_mean'] = scare_mean

        res['chr_total_reads'] = chr_total_reads

        res['chrs_length'] = chrs_length

        res['chrfrcount'] = chrsfrcount

        res['filterreadscount'] = filterreadscount

        res['totalreads'] = totalreads

        res['chruniqlength'] = chruniqlength

        # res['chrreadlengthmean'] = chrreadlengthmean

        totallengmean = 0

        totalchrnumber = 0

        for chromsome in count_chr:

            if chromsome in chrreadlengthmean:

                totallengmean = totallengmean + chrreadlengthmean[chromsome]

                totalchrnumber = totalchrnumber + 1

        readlengthmean = totallengmean/totalchrnumber

        res['readlengthmean'] = readlengthmean

        pool.close()

        return res

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

        print ('got exception in Jazzlib.FRegion.filter_region: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:
        #     print ('joining pool processes')
        pool.join()
            # print ('join complete')


def chrwindow_counter(par):

    try:

        chromosome = par['chrmosome']

        windowsize = par['windowsize']

        chr_length = par['chr_length']

        bamfile = par['bamfile']

        maxinsert = par['maxinsert']

        jobtype = par['jobtype']

        windowcount = windowscarecounter(bamfile=bamfile, regionchromosome=chromosome,
                                         regionstart=1, regionend=chr_length,
                                         windowsize=windowsize, maxinsert=maxinsert, jobtype=jobtype)

        uniqcount = uniqsitecount(bamfile=bamfile, regionchromosome=chromosome,
                                  regionstart=1, regionend=chr_length, maxinsert=maxinsert,
                                  jobtype=jobtype)

        readlengthmean = readslengthmean(bamfile=bamfile, regionchromosome=chromosome,
                                        regionstart=1, regionend=chr_length, maxinsert=maxinsert,
                                        jobtype=jobtype)

        chrwindowcount = dict()

        chrwindowcount['windowcount'] = windowcount

        chrwindowcount['chromosome'] = chromosome

        chrwindowcount['uniqcount'] = uniqcount

        chrwindowcount['readlengthmean'] = readlengthmean

        # for debug
        print("in chrwindow_counter", readlengthmean)

        return chrwindowcount

    except KeyboardInterrupt:

        print ("You cancelled the program!")

        sys.exit(1)
