

from scipy.special import gammaincc
from scipy import math
import scipy.stats as stats
from decimal import Decimal, localcontext
from .Peak import *
import sys

def bionompvalue(x, n, p):

    bionompvalue = 1 - stats.binom.cdf(x, n, p)

    return bionompvalue


def poissonpvalue(x,mu):

    poissonpvalue = Decimal(1) - Decimal(stats.poisson.cdf(x, mu))

    return poissonpvalue



def fdr(pnow, plist, prank):
    #FDR=length(pvalue)*pvalue/rank(pvalue)

    rankofplist = prank

    lengthofplist = len(plist)

    for i in range(0,lengthofplist):

        if plist[i] == pnow:
            now_rank = rankofplist[i]
            fdr = lengthofplist*pnow/now_rank
            fdr = min(1,fdr)
            break

    return fdr


def bayesfactor(locallambda, peakscore):

    try:

        # bayesfactor = 2 * (math.log((gammaincc(peakscore-1, locallambda)*gamma(peakscore-1)), math.e) - (peakscore-1)*math.log(locallambda, math.e) + locallambda)
        #
        # a = (math.log(gammaincc(peakscore-1, locallambda), math.e) )
        # b = math.lgamma(peakscore-1)
        # c=(peakscore-1)*math.log(locallambda, math.e)
        # print (locallambda,peakscore,a,b,c)
        bayesfactor2 = 2 * (math.log(gammaincc(peakscore-1, locallambda), math.e)+math.lgamma(peakscore-1) - (peakscore-1)*math.log(locallambda, math.e) + locallambda)

        return bayesfactor2

    except Exception as e:

        print ('got exception in Jazzlib.sta.bayesfactor: %r,' % (e,))

        print (locallambda, peakscore)

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)

def fdr_control(chippeaks, inputpeaks, fdr):

    fdrpeakdict = dict()

    chipscore = list()

    inputscore = list()

    overlaptedpeak = dict()

    fdrth = -1

    # print ("check fdr")

    for inputpeak in inputpeaks:

        start = inputpeak.start

        end = inputpeak.end

        inputscore.append(inputpeak.score)

        for chippeak in chippeaks:

            if chippeak.chromosome == inputpeak.chromosome:

                if chippeak.peakpoint == inputpeak.peakpoint:

                    overlaptedpeak[chippeak.peakid] = dict()

                    overlaptedpeak[chippeak.peakid]['inputscore'] = inputpeak.score

                    overlaptedpeak[chippeak.peakid]['chipscore'] = chippeak.score

                    chipscore.append(chippeak.score)

                    # print(chippeak.chromosome, chippeak.peakpoint, chippeak.score, inputpeak.score, chippeak.peakid)

    for i in sorted(chipscore):

        # print("score", i)

        chippeakcount = 0.0

        inputpeakcount = 0.0

        for peakid in overlaptedpeak:

            if i <= overlaptedpeak[peakid]['inputscore']:

                inputpeakcount = inputpeakcount + 1

        for chippeak in chippeaks:

            if i <= chippeak.score:

                chippeakcount = chippeakcount + 1

        nowfdr = inputpeakcount/chippeakcount

        # print (i, chippeakcount, inputpeakcount, nowfdr)

        if chippeakcount == 0:

            break

        for peaknow in chippeaks:

            if peaknow.score > i:

                peaknow.fdr = nowfdr


        # if (inputpeakcount/chippeakcount) < fdr:
        #
        #     fdrth = i
        #
        #     break

    return chippeaks



def fdr_control2(chippeaks, inputpeaks, fdr):

    fdrpeakdict = dict()

    chipscore = list()

    inputscore = list()

    overlaptedpeak = dict()

    fdrth = -1

    # print ("check fdr")

    for inputpeak in inputpeaks:

        start = inputpeak.start

        end = inputpeak.end

        inputscore.append(inputpeak.score)

        for chippeak in chippeaks:

            if chippeak.chromosome == inputpeak.chromosome:

                if inputpeak.start <chippeak.peakpoint < inputpeak.end:

                    overlaptedpeak[chippeak.peakid] = dict()

                    overlaptedpeak[chippeak.peakid]['inputscore'] = inputpeak.score

                    overlaptedpeak[chippeak.peakid]['chipscore'] = chippeak.score

                    chipscore.append(chippeak.score)

                    # print(chippeak.chromosome, chippeak.peakpoint, chippeak.score, inputpeak.score, chippeak.peakid)

    for i in sorted(chipscore):

        # print("score", i)

        chippeakcount = 0.0

        inputpeakcount = 0.0

        for peakid in overlaptedpeak:

            if i <= overlaptedpeak[peakid]['inputscore']:

                inputpeakcount = inputpeakcount + 1

        for chippeak in chippeaks:

            if i <= chippeak.score:

                chippeakcount = chippeakcount + 1

        nowfdr = inputpeakcount/chippeakcount

        if chippeakcount == 0:

            break

        for peaknow in chippeaks:

            if peaknow.score > i:

                peaknow.fdr = nowfdr

    return chippeaks


def fdr_bh(peaks):

    b01s = list()

    peakscores = list()

    for peak in peaks:

        b01 = 1/(math.e**(peak.score/2))

        peakscores.append(peak.score)

        b01s.append(b01)

    sortedb01s = sorted(b01s,reverse=True)

    listlength = len(sortedb01s)

    for peak in peaks:

        b01 = 1/(math.e**(peak.score/2))

        rank = 1

        for i in range(0,listlength):

            if sortedb01s[i] == b01:

                rank = i + 1

                break

        fdr = b01*listlength/rank

        peak.fdr = fdr

    return peaks





if __name__ == "__main__":

    try:

        for i in range(100,2000,100):

            for j in range (2,80):


                bs = bayesfactor(locallambda=i, peakscore=j)
                # if bs == 1500:
                #     bs = 'error'
                print ("locallambda:",i, "peakscore",j,"bs",bs)

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)