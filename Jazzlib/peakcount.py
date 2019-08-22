

from .Peak import *
from .sta import *
from .Hotspot import *
from .readscounter import *
import sys
from .region import *
from multiprocessing import Pool
from .bgcount import *

class KeyboardInterruptError(Exception):

    pass


def peakcount_nocontrol(hotspots, pvalue, bamfile, initiallength, jobtype, nthreads, fregion,bpc, maxinsert=10000):

    #bpc = get_bpc(bamfile, hotspots, jobtype, fregion.filted_region)

    pars = list()

    for hotspot in hotspots:

        # print (hotspot.chromosome,hotspot.start,hotspot.end,hotspot.hotspotid)

        whether_filter = False

        for i in range(int(hotspot.start/100)-1,int(hotspot.end/100)+1):

            if hotspot.chromosome in fregion.filted_region:
                if i in fregion.filted_region[hotspot.chromosome]:

                    whether_filter = True


        if not whether_filter:

            par = dict()

            par['pvalue'] = pvalue

            par['hotspot'] = hotspot

            par['bamfile'] = bamfile

            par['initiallength'] = initiallength

            par['chrlength'] = fregion.chrs_length[hotspot.chromosome]

            par['bpc'] = bpc

            par['jobtype'] = jobtype

            par['maxinsert'] = maxinsert

            pars.append(par)

    pool=Pool(nthreads)

    chrpeaks = list()

    try:
        hotspotpeaks = pool.map(hotspotpeak_nocontrol, pars)

        for nowpeaks in hotspotpeaks:

            numberofpeak = len(nowpeaks)

            if numberofpeak == 0:

                pass

                # print ("no peak")

            else:

                for peak in nowpeaks:

                    # print ("peak",peak.parent,peak.chromosome,peak.start,peak.end,peak.peakid)

                    chrpeaks.append(peak)

            # print ()
        pool.close()

        return chrpeaks


    except KeyboardInterrupt:
        pool.terminate()
        print ("You cancelled the program!")
        sys.exit(1)
    except Exception as e:
        print ('got exception: %r, terminating the pool' % (e,))
        pool.terminate()
        print ('pool is terminated')
    finally:
        # print ('joining pool processes')
        pool.join()
        # print ('join complete')



def hotspotpeak_nocontrol(par):

    try:
        pvalue = par['pvalue']

        hotspot = par['hotspot']

        bamfile = par['bamfile']

        initiallength = par['initiallength']

        chrlength=par['chrlength']

        bpc = par['bpc']

        maxinsert = par['maxinsert']

        start = hotspot.start

        end = hotspot.end

        hotspotlength = end - start + 1

        chromosome = hotspot.chromosome

        parentid = hotspot.hotspotid

        hotspottype = hotspot.hotspottype

        jobtype = hotspot.jobtype

        # regionstart = start
        #
        # regionend = end

        regionstart = start - hotspotlength

        regionend = end + hotspotlength

        peaks = list()

        if regionstart < 1:

            regionstart = 1

        if regionend>chrlength:

            regionend = chrlength

        hotspotregio = chromosome + ':' + str(regionstart) + '-' + str(regionend)

        readscount = dict()

        if jobtype == 'nhsingle':

            readscount = nhreadscounter(bamfile = bamfile, region=hotspotregio, paired=False, maxinsert = maxinsert)

        elif jobtype == 'nhpaired':

            readscount = nhreadscounter(bamfile = bamfile, region=hotspotregio, paired=True, maxinsert = maxinsert)

        elif jobtype == 'dh':

            readscount = dhreadscounter(bamfile = bamfile, region=hotspotregio)

        else:

            print ("%s count type error!!!!" % jobtype)

            sys.exit(1)

        totalreads = 0

        for i in range(regionstart,regionend):
            if i in readscount:
                # uniqcount = uniqcount + 1
                totalreads = totalreads + readscount[i]

        lambda_mu = (totalreads+0.0) /(regionend - regionstart +1) * initiallength * 3

        #lambda_mu = (totalreads+0.0) /(regionend - regionstart +1) * initiallength

        avg_mu = bpc * initiallength

        lambda_mu = max(lambda_mu, avg_mu)

        peaks = list()
        
        #if end-start+1 < initiallength:
        if end-start+1 < 1:

            pass

        else:

            region_site_pvalue = dict()

            region_point = dict()

            for window_start in range(start,end-initiallength+1):
                window_end = window_start + initiallength
                readsinwindow = 0
                for j in range(window_start,window_end+1):
                    if j in readscount:
                        readsinwindow = readsinwindow + readscount[j]
                now_pvalue = poissonpvalue(x=readsinwindow, mu =lambda_mu)
                # print (window_start,window_end,readsinwindow,lambda_mu,now_pvalue)
                region_site_pvalue[window_start] = now_pvalue

            #filter pvalue

            uniqpoint = list()

            for now_window_start in region_site_pvalue:
                if region_site_pvalue[now_window_start] < pvalue:
                    # print (now_window_start)
                    for now_site in range(now_window_start,now_window_start+initiallength):
                        region_point[now_site] = 1

            uniqpoint = list(region_point.keys())

            if region_point:

                #peaksregion = continueregion(points=uniqpoint, minlength=initiallength)
                peaksregion = continueregion(points=uniqpoint, minlength=2)

                initid = 1

                for now_region in peaksregion:

                    peaks_reads = 0

                    start_site = now_region['start_site']

                    end_site = now_region['end_site']

                    for i in range(start_site, end_site+1):

                        if i in readscount:

                            peaks_reads = peaks_reads + readscount[i]

                    regionlength = end_site - start_site +1

                    region_lambda_mu = (totalreads+0.0) / (regionend - regionstart +1) * regionlength

                    avg_region_lambda = bpc * regionlength

                    region_lambda_mu = max(region_lambda_mu, avg_region_lambda)

                    totalpvalue = 0

                    minsite = 0

                    peakscore = 0

                    peaksite = 1

                    for i in range(start_site,end_site+1):
                        # if i in region_site_pvalue:
                        #     totalpvalue = totalpvalue + region_site_pvalue[i]
                        # else:
                        #     print ("can't find p in", i)
                        if i in readscount:

                            if readscount[i] > peakscore:

                                peakscore = readscount[i]

                                peaksite = i





                    now_region_pvalue = poissonpvalue(x=peaks_reads, mu = region_lambda_mu)

                    if now_region_pvalue < pvalue:

                        peakid = parentid + '.' + str(initid)

                        initid = initid + 1

                        nowpeak = Peak(chromosome=chromosome, start=start_site, end=end_site,
                                       pvalue=float(now_region_pvalue), peakpoint=peaksite,
                                       parent=parentid, peakid=peakid)

                        peaks.append(nowpeak)

        return peaks

    except KeyboardInterrupt:

        raise KeyboardInterruptError()


def peakcount_control(hotspots, pvalue, datafile, controlfile, initiallength, jobtype, nthreads, datafregion,
                    controlfregion, databpc, controlbpc, maxinsert):
    pars = list()

    for hotspot in hotspots:

        # print (hotspot.chromosome,hotspot.start,hotspot.end,hotspot.hotspotid)

        whether_filter = False

        for i in range(int(hotspot.start/100)-1,int(hotspot.end/100)+1):

            if hotspot.chromosome in controlfregion.filted_region:

                if i in controlfregion.filted_region[hotspot.chromosome]:

                    whether_filter = True


        if not whether_filter:

            par = dict()

            par['pvalue'] = pvalue

            par['hotspot'] = hotspot

            par['datafile'] = datafile

            par['controlfile'] = controlfile

            par['initiallength'] = initiallength

            par['chrlength'] = controlfregion.chrs_length[hotspot.chromosome]

            par['databpc'] = databpc

            par['controlbpc'] = controlbpc

            par['jobtype'] = jobtype

            par['datafregion'] = datafregion

            par['controlfregion'] = controlfregion

            par['maxinsert'] = maxinsert

            pars.append(par)

    pool=Pool(nthreads)

    chrpeaks = list()

    try:
        hotspotpeaks = pool.map(hotspotpeak_control, pars)

        for nowpeaks in hotspotpeaks:

            numberofpeak = len(nowpeaks)

            if numberofpeak == 0:

                pass

                # print ("no peak")

            else:

                for peak in nowpeaks:

                    # print ("peak",peak.parent,peak.chromosome,peak.start,peak.end,peak.peakid)

                    chrpeaks.append(peak)

            # print ()
        pool.close()

        return chrpeaks

    except KeyboardInterrupt:
        pool.terminate()
        print ("You cancelled the program!")
        sys.exit(1)
    except Exception as e:
        print ('got exception: %r, terminating the pool' % (e,))
        pool.terminate()
        print ('pool is terminated')
    finally:
        # print ('joining pool processes')
        pool.join()
        # print ('join complete')


def hotspotpeak_control(par):

    try:

        pvalue = par['pvalue']

        hotspot = par['hotspot']

        datafile = par['datafile']

        controlfile = par['controlfile']

        initiallength = par['initiallength']

        chrlength = par['chrlength']

        databpc = par['databpc']

        controlbpc = par['controlbpc']

        jobtype = par['jobtype']

        datafregion = par['datafregion']

        controlfregion = par['controlfregion']

        maxinsert = par['maxinsert']

        start = hotspot.start

        end = hotspot.end

        hotspotlength = end - start + 1

        chromosome = hotspot.chromosome

        parentid = hotspot.hotspotid

        hotspottype = hotspot.hotspottype

        jobtype = hotspot.jobtype

        regionstart = start - initiallength

        regionend = end + initiallength

        datachrtotalreads = datafregion.chr_total_reads[chromosome]

        controlchrtotalreads = controlfregion.chr_total_reads[chromosome]

        readsratio = datachrtotalreads/controlchrtotalreads

        peaks = list()

        if regionstart < 1:

            regionstart = 1

        if regionend>chrlength:

            regionend = chrlength

        hotspotregio = chromosome + ':' + str(regionstart) + '-' + str(regionend)

        datareadscount = dict()

        controlreadscount = dict()

        if jobtype == 'nhsingle':

            datareadscount = nhreadscounter(bamfile = datafile, region=hotspotregio, paired=False, maxinsert=maxinsert)

            controlreadscount = nhreadscounter(bamfile = controlfile, region=hotspotregio, paired=False, maxinsert=maxinsert)

        elif jobtype == 'nhpaired':

            datareadscount = nhreadscounter(bamfile = datafile, region=hotspotregio, paired=True, maxinsert=maxinsert)

            controlreadscount = nhreadscounter(bamfile =  controlfile, region=hotspotregio, paired=True, maxinsert=maxinsert)

        elif jobtype == 'dh':

            datareadscount = dhreadscounter(bamfile = datafile, region=hotspotregio)

            controlreadscount = dhreadscounter(bamfile = controlfile, region=hotspotregio)

        else:

            print ("%s count type error!!!!" % jobtype)

            sys.exit(1)

        datatotalreads = 0

        controltotalreads = 0

        for i in range(regionstart,regionend):
            if i in datareadscount:
                    # uniqcount = uniqcount + 1
                datatotalreads = datatotalreads + datareadscount[i]
            if i in controlreadscount:
                controltotalreads = controltotalreads + controlreadscount[i]

        datalambda_mu = (datatotalreads+0.0) /(regionend - regionstart +1) * initiallength

        controllambda_mu = (controltotalreads+0.0) /(regionend - regionstart +1) * initiallength * readsratio

        dataavg_mu = databpc * initiallength

        controlavg_mu = controlbpc * initiallength * readsratio

        # print ("datalambda_mu, controllambda_mu, dataavg_mu, controlavg_mu", datalambda_mu, controllambda_mu, dataavg_mu, controlavg_mu)

        lambda_mu = max(datalambda_mu, controllambda_mu, dataavg_mu, controlavg_mu)

        #if end-start+1 < initiallength:
        if end-start+1 < 1:

            pass

        else:
            region_site_pvalue = dict()

            region_point = dict()

            for window_start in range(start,end-initiallength+1):
                window_end = window_start + initiallength
                readsinwindow = 0
                for j in range(window_start,window_end+1):
                    if j in controlreadscount:
                            readsinwindow = readsinwindow + controlreadscount[j]
                now_pvalue = poissonpvalue(x=readsinwindow, mu =lambda_mu)
                    # print (window_start,window_end,readsinwindow,lambda_mu,now_pvalue)
                region_site_pvalue[window_start] = now_pvalue

                #filter pvalue

            uniqpoint = list()

            for now_window_start in region_site_pvalue:
                if region_site_pvalue[now_window_start] < pvalue:
                    # print (now_window_start)
                    for now_site in range(now_window_start,now_window_start+initiallength):
                        region_point[now_site] = 1

            uniqpoint = list(region_point.keys())

            if region_point:

                #peaksregion = continueregion(points=uniqpoint, minlength=initiallength)
                peaksregion = continueregion(points=uniqpoint, minlength=2)

                initid = 1

                for now_region in peaksregion:

                    peaks_reads = 0

                    start_site = now_region['start_site']

                    end_site = now_region['end_site']

                    for i in range(start_site, end_site+1):

                        if i in datareadscount:

                            peaks_reads = peaks_reads + datareadscount[i]

                    regionlength = end_site - start_site +1

                    region_lambda_mu = (datatotalreads+0.0) / (regionend - regionstart +1) * regionlength

                    avg_region_lambda = databpc * regionlength

                    region_lambda_mu = max(region_lambda_mu, avg_region_lambda)

                    totalpvalue = 0

                    minsite = 0

                    peakscore = 0

                    peaksite = 1

                    for i in range(start_site,end_site+1):
                        # if i in region_site_pvalue:
                        #     totalpvalue = totalpvalue + region_site_pvalue[i]
                        # else:
                        #     print ("can't find p in", i)
                        if i in datareadscount:

                            if datareadscount[i] > peakscore:

                                peakscore = datareadscount[i]

                                peaksite = i





                    now_region_pvalue = poissonpvalue(x=peaks_reads, mu = region_lambda_mu)

                    if now_region_pvalue < pvalue:

                        peakid = parentid + '.' + str(initid)

                        initid = initid + 1

                        nowpeak = Peak(chromosome=chromosome, start=start_site, end=end_site,
                                           pvalue=float(now_region_pvalue), peakpoint=peaksite,
                                           parent=parentid, peakid=peakid)

                        print (chromosome, start_site, end_site)

                        peaks.append(nowpeak)

        return peaks

    except KeyboardInterrupt:

        raise KeyboardInterruptError()





def main():
    pass

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)