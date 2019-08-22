

from multiprocessing import Pool
from .FRegion import *
import random as rnd
from .kernel import *
from numpy import *

class KeyboardInterruptError(Exception):

    pass


def randombg2(fregion, nthreads, maxinsert, randomthreshold=2, runtime=1000, randomwindow=10000):

    countgenomelength = fregion.countgenomelength

    adjreads = fregion.adjreads

    bg = adjreads/countgenomelength

    return bg


def randombg(fregion, nthreads, maxinsert, randomthreshold=2, runtime=1000, randomwindow=10000):

    pool = Pool(nthreads)

    try:



        countgenomeuniqlength = fregion.countgenomeuniqlength

        adjreads = fregion.adjreads

        countgenomelength = fregion.countgenomelength

        uniqrate = countgenomeuniqlength/countgenomelength

        if uniqrate <0.5:

            uniqrate = uniqrate * 2

        countreads = int(adjreads/countgenomeuniqlength * randomwindow)+1

        onekernel = smooth_kernel(length=maxinsert)

        kernel_score = list()

        pars = list()

        for i in sorted(onekernel):

            kernel_score.append(onekernel[i])

        for j in range(runtime):

            par = dict()

            par['countreads'] = countreads

            par['kernel_score'] = kernel_score

            par['uniqrate'] = uniqrate

            par['randomwindo'] = randomwindow

            par['randomthreshold'] =randomthreshold

            # print (par)

            pars.append(par)

        randths = pool.map(sim_bg_worker, pars)

        thsum = 0

        for randth in randths:

            thsum = thsum + randth

        random_th = thsum/runtime

        pool.close()

        return random_th

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

        print ('got exception in Jazzlib.randombg.randombg: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:
        # print ('joining pool processes')
        pool.join()


def sim_bg_worker(par):

    try:

        countreads = par['countreads']

        kernel_score = par['kernel_score']

        uniqrate = par['uniqrate']

        randomwindow = par['randomwindo']

        randomthreshold = par['randomthreshold']

        totaluniqsite = int(uniqrate * randomwindow)

        rand_reads_count = list()

        region_site = list(range(0, randomwindow))

        for i in range(0, randomwindow):

            rand_reads_count.append(0)

        sim_uniqsite = rnd.sample(region_site, totaluniqsite)


        for k in range(0, countreads):

            rand_number = int(rnd.uniform(0, totaluniqsite))

            rand_reads = sim_uniqsite[rand_number]

            # print (rand_number, rand_reads)

            rand_reads_count[rand_reads] = rand_reads_count[rand_reads] + 1

        smoothed_result = correlate(array(rand_reads_count), kernel_score, "same")

        # scores = list()

        rand_mean = smoothed_result.mean()

        rand_std = smoothed_result.std()

        # total_sum = smoothed_result.sum()
        # print (rand_mean, rand_std, randomthreshold)

        rand_threshhold = rand_mean + randomthreshold * rand_std

        return rand_threshhold


    except KeyboardInterrupt:

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception as e:

        print ('got exception in Jazzlib.randombg.sim_bg_worker: %r, terminating the pool' % (e,))


if __name__ == "__main__":
    try:

        onekernel = smooth_kernel(length=100)

        kernel_score = list()

        pars = list()

        for i in sorted(onekernel):

            kernel_score.append(onekernel[i])

        par = dict()

        par['countreads'] = 100000

        par['kernel_score'] = kernel_score

        par['uniqrate'] = 0.3

        par['randomwindo'] = int(1e5)

        par['randomthreshold'] = 3

        th = sim_bg_worker(par)

        print (th)

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)