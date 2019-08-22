
from .countreads import *
from .kernel import *
import numpy as np


class KeyboardInterruptError(Exception):

    pass


def regionsmooth(bamfile, jobtype, maxinsert, regionchromosome, regionstart, regionend, chr_length):

    try:

        renewstart = regionstart - maxinsert*2

        renewend = regionend + maxinsert*2

        if renewstart < 1:

            renewstart = 1

        if renewend > chr_length:

            renewend = chr_length

        insertsize_middle_site_count = midsiteinsersizecounter(bamfile=bamfile, regionchromosome=regionchromosome,
                                                               regionstart=renewstart, regionend=renewend,
                                                               jobtype=jobtype, maxinsert=maxinsert)


        renewlength = renewend - renewstart + 1

        smoothed_score = np.repeat(0, renewlength)

        for insertlen in insertsize_middle_site_count:

            # print ("count size", insertlen)

            readcount_nowinsertsize = list()

            kernelnow = smooth_kernel(insertlen)

            kernel_score = list()

            for w in sorted(kernelnow):

                kernel_score.append(kernelnow[w])

            for n in range(renewstart, renewend+1):

                nowscore = 0

                if n in insertsize_middle_site_count[insertlen]:

                    nowscore = insertsize_middle_site_count[insertlen][n]

                readcount_nowinsertsize.append(nowscore)

            nowsmoothed = np.correlate(np.array(readcount_nowinsertsize), kernel_score, "same")

            smoothed_score = nowsmoothed + smoothed_score

        outputscore = dict()

        outputscore['chromosome'] = regionchromosome

        outputscore['score'] = dict()

        # print (smoothed_score[0])

        for j in range(0, renewlength):

            nowsite = j + renewstart

            nowscore = smoothed_score[j]

            if regionstart <= nowsite <= regionend:

                outputscore['score'][nowsite] = nowscore

        return outputscore

    except KeyboardInterrupt:

        raise KeyboardInterruptError()
        sys.exit(0)